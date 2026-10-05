"""
gui.py — la interfaz de `pulsos/`: elegir un .xml, leer pulsos, definir el
intervalo de integracion ARRASTRANDO sobre el propio grafico, ver area y FWHM
de cada pulso, y comparar muchos pulsos superpuestos.

No llama a nada de `analisis/gui.py` (ese arrastra `simulacion`, `calibracion`,
`pipeline`, `registro`, que no tienen nada que ver con esta herramienta); el
patron `Lienzo` se copia porque es plomeria generica de Tk, no logica del
proyecto.

La pieza nueva es el `SpanSelector` de la pestaña "Intervalo": no hay ningun
precedente de `matplotlib.widgets` en el repo. Se recrea cada vez que se
redibuja el eje (porque `Lienzo.dibujar` limpia la Figure entera) y sus
`extents` se restauran al `(pre, post)` vigente, para que el intervalo elegido
sobreviva a cambiar de pulso o de vista.

Correr `python pulsos/gui.py`.
"""

import os
import queue
import sys
import threading
import tkinter as tk
from tkinter import filedialog, ttk

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import pulso
from pulsos import rutas, ancho, lote

ARCHIVO_SUELTO = "(archivo suelto)"

# nombre de cada bandera, para decodificar en la tabla de resultados. Los bits
# de pulso.py y los nuevos de ancho.py no colisionan, asi que se combinan con
# un simple OR de bandera_area | bandera_fwhm antes de decodificar.
_NOMBRES_BANDERA = (
    (pulso.MAL_SATURA, "SATURA"),
    (pulso.MAL_TRUNCADO, "TRUNCADO"),
    (pulso.MAL_SIN_PRE, "SIN_PRE"),
    (pulso.MAL_SIN_PICO, "SIN_PICO"),
    (ancho.MAL_SIN_CRUCE_SUBIDA, "SIN_CRUCE_SUBIDA"),
    (ancho.MAL_SIN_CRUCE_BAJADA, "SIN_CRUCE_BAJADA"),
)


def _decodificar_bandera(bits):
    if bits == 0:
        return ""
    return "+".join(nombre for bit, nombre in _NOMBRES_BANDERA if bits & bit)


def _media_dispersion_sem(x):
    """
    `(media, dispersion, sem)` de un array 1D.

    `dispersion` es el std MUESTRAL (ddof=1: divide por N-1) de los valores
    pulso a pulso — cuanto varia el area/FWHM de un evento a otro. Bajo
    irradiacion eso es fisico y esperable (no todos los eventos depositan la
    misma energia), no un problema de la medicion.

    `sem` es el error estandar DE LA MEDIA, `dispersion / sqrt(N)`: que tan
    bien determinado esta el promedio reportado, no cuanto varian los pulsos
    individuales. Es la cantidad que importa para decir si "la media vale X"
    es un numero confiable.

    Con un solo pulso (N=1) `dispersion` no esta definida (ddof=1 divide por
    0): se devuelve 0.0 para no romper el formateo, ya que con N=1 la media
    tampoco tiene ningun respaldo estadistico.
    """
    x = np.asarray(x, float)
    n = x.size
    media = float(np.mean(x))
    if n < 2:
        return media, 0.0, 0.0
    disp = float(np.std(x, ddof=1))
    return media, disp, disp / np.sqrt(n)


def activar_dpi_awareness():
    """
    Le dice a Windows que esta app maneja su propio DPI.

    Sin esto, en pantallas con escalado Windows renderiza a 96 DPI y estira el
    resultado como bitmap: se ve borroso y ninguna opcion de matplotlib lo
    arregla. Hay que llamarlo ANTES de crear la ventana Tk.
    """
    if sys.platform != "win32":
        return
    import ctypes
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(2)
    except Exception:
        try:
            ctypes.windll.user32.SetProcessDPIAware()
        except Exception:
            pass


# ==========================================================================
# Lienzo: figura de matplotlib metida en un frame de Tk
# ==========================================================================
class Lienzo(ttk.Frame):
    def __init__(self, master, dpi=100):
        super().__init__(master)
        self.fig = Figure(figsize=(9, 5.2), dpi=dpi)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        barra = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        barra.update()
        barra.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.vacio("todavia no hay pulsos cargados")

    def vacio(self, texto):
        self.fig.clear()
        a = self.fig.add_subplot(111)
        a.text(0.5, 0.5, texto, ha="center", va="center", fontsize=11,
              color="0.45", transform=a.transAxes)
        a.set_axis_off()
        self.canvas.draw_idle()

    def dibujar(self, fn):
        """`fn(fig)` dibuja sobre la Figure ya limpia. Un error se muestra en
        el propio canvas en vez de tirar abajo la interfaz."""
        self.fig.clear()
        try:
            fn(self.fig)
        except Exception as e:
            self.fig.clear()
            a = self.fig.add_subplot(111)
            a.text(0.5, 0.5, "no se pudo dibujar:\n%s" % e, ha="center",
                  va="center", fontsize=10, color="#b32b2b",
                  transform=a.transAxes)
            a.set_axis_off()
        self.canvas.draw_idle()


# ==========================================================================
# La aplicacion
# ==========================================================================
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Pulsos — FWHM, area e intensidad de pulsos individuales")
        self.geometry("1360x840")
        self.minsize(1080, 680)

        self._configurar_nitidez()

        # --- estado ---
        self.T = self.V = self.relojes = self.scalers = None
        self.canal = None
        self.ruta_actual = None
        self.pre = self.post = None
        self.vent = None            # dict de lote.sugerir_ventana
        self.datos_lote = None      # dict de lote.calcular_lote
        self._span = None

        self.cola = queue.Queue()

        self._construir()
        self.after(80, self._vaciar_cola)
        self._al_elegir_corrida()

    # ------------------------------------------------------------ nitidez --
    def _configurar_nitidez(self):
        dpi = self.winfo_fpixels("1i")
        self.tk.call("tk", "scaling", dpi / 72.0)
        self.dpi_fig = max(int(100 * dpi / 96.0), 100)
        s = ttk.Style()
        if "vista" in s.theme_names():
            s.theme_use("vista")
        s.configure("TLabelframe.Label", foreground="#1f4e79")

    # ----------------------------------------------------------- widgets --
    def _construir(self):
        panel = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        panel.pack(fill=tk.BOTH, expand=True)

        izq = ttk.Frame(panel, padding=8)
        panel.add(izq, weight=0)
        self._construir_panel_datos(izq)

        self.notebook = ttk.Notebook(panel)
        panel.add(self.notebook, weight=1)

        self._construir_tab_intervalo()
        self._construir_tab_resultados()
        self._construir_tab_comparacion()
        self.notebook.bind("<<NotebookTabChanged>>", self._al_cambiar_tab)

    def _construir_panel_datos(self, master):
        caja = ttk.LabelFrame(master, text=" Datos ", padding=8)
        caja.pack(fill=tk.X)

        ttk.Label(caja, text="Corrida:").grid(row=0, column=0, sticky="w")
        claves = sorted(rutas.disponibles()) + [ARCHIVO_SUELTO]
        self.var_corrida = tk.StringVar(value=claves[0] if claves else ARCHIVO_SUELTO)
        combo = ttk.Combobox(caja, textvariable=self.var_corrida, values=claves,
                             state="readonly", width=22)
        combo.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 4))
        combo.bind("<<ComboboxSelected>>", lambda e: self._al_elegir_corrida())

        ttk.Label(caja, text="Archivo (.xml):").grid(row=2, column=0, sticky="w")
        self.var_ruta = tk.StringVar(value="")
        self.entry_ruta = ttk.Entry(caja, textvariable=self.var_ruta, width=24,
                                    state="disabled")
        self.entry_ruta.grid(row=3, column=0, sticky="ew")
        self.btn_elegir = ttk.Button(caja, text="...", width=3, state="disabled",
                                     command=self._elegir_archivo)
        self.btn_elegir.grid(row=3, column=1, sticky="w", padx=(4, 0))

        ttk.Label(caja, text="Pulsos a leer:").grid(row=4, column=0, sticky="w",
                                                    pady=(8, 0))
        self.var_n = tk.IntVar(value=500)
        ttk.Spinbox(caja, from_=50, to=5000, increment=50, textvariable=self.var_n,
                   width=8).grid(row=5, column=0, sticky="w")

        ttk.Label(caja, text="Atenuador (dB):").grid(row=6, column=0, sticky="w",
                                                      pady=(8, 0))
        self.var_atenuador = tk.StringVar(value="0")
        ttk.Entry(caja, textvariable=self.var_atenuador, width=8
                 ).grid(row=7, column=0, sticky="w")
        ttk.Label(caja, text="voltajes ×10 por 20 dB", foreground="gray45",
                 wraplength=190).grid(row=8, column=0, columnspan=2, sticky="w")

        self.btn_cargar = ttk.Button(caja, text="Cargar pulsos",
                                     command=self._cargar_pulsos)
        self.btn_cargar.grid(row=9, column=0, columnspan=2, sticky="ew", pady=(8, 4))

        self.var_estado = tk.StringVar(value="listo")
        ttk.Label(caja, textvariable=self.var_estado, foreground="#1f4e79",
                 wraplength=190).grid(row=10, column=0, columnspan=2, sticky="w")
        self.barra = ttk.Progressbar(caja, length=190, maximum=100.0)
        self.barra.grid(row=11, column=0, columnspan=2, sticky="ew", pady=(2, 0))

        caja.columnconfigure(0, weight=1)

    # -------------------------------------------------------- tab Intervalo --
    def _construir_tab_intervalo(self):
        f = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(f, text="Intervalo")

        ctrl = ttk.Frame(f)
        ctrl.pack(fill=tk.X, pady=(0, 4))

        ttk.Label(ctrl, text="Vista:").pack(side=tk.LEFT)
        self.var_vista = tk.StringVar(value="Pulso individual")
        ttk.Combobox(ctrl, textvariable=self.var_vista, state="readonly", width=16,
                    values=["Pulso individual", "Forma promedio"]
                    ).pack(side=tk.LEFT, padx=(4, 12))
        self._combo_vista_bind = ctrl.winfo_children()[-1]
        self._combo_vista_bind.bind("<<ComboboxSelected>>",
                                    lambda e: self._redibujar_intervalo())

        ttk.Label(ctrl, text="Pulso #:").pack(side=tk.LEFT)
        self.var_idx = tk.IntVar(value=0)
        self.spin_idx = ttk.Spinbox(ctrl, from_=0, to=0, textvariable=self.var_idx,
                                    width=6, command=self._redibujar_intervalo)
        self.spin_idx.pack(side=tk.LEFT, padx=(4, 12))

        ttk.Button(ctrl, text="Sugerir automaticamente",
                  command=self._sugerir_ventana).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(ctrl, text="Calcular area y FWHM para todos",
                  command=self._calcular_lote).pack(side=tk.LEFT)

        self.var_ancho = tk.StringVar(value="sin intervalo definido todavia")
        ttk.Label(f, textvariable=self.var_ancho, foreground="#1f4e79"
                 ).pack(anchor="w", pady=(0, 4))

        self.lienzo_intervalo = Lienzo(f, dpi=self.dpi_fig)
        self.lienzo_intervalo.pack(fill=tk.BOTH, expand=True)

    # ------------------------------------------------------- tab Resultados --
    def _construir_tab_resultados(self):
        f = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(f, text="Resultados")

        self.var_resumen = tk.StringVar(value="todavia no hay nada calculado")
        ttk.Label(f, textvariable=self.var_resumen, foreground="#1f4e79"
                 ).pack(anchor="w", pady=(0, 4))

        # tabla de resumen: N, ventana aplicada, y estadistica de area/FWHM
        # sobre el lote entero (esto es lo que va en la seccion de resultados
        # de un informe, la tabla de abajo es para mirar pulso por pulso).
        caja_resumen = ttk.LabelFrame(f, text=" Resumen del lote ", padding=4)
        caja_resumen.pack(fill=tk.X, pady=(0, 6))
        cols_r = ("campo", "valor")
        self.tabla_resumen = ttk.Treeview(caja_resumen, columns=cols_r,
                                          show="headings", height=7)
        self.tabla_resumen.heading("campo", text="")
        self.tabla_resumen.heading("valor", text="")
        self.tabla_resumen.column("campo", width=170, anchor="w")
        self.tabla_resumen.column("valor", width=320, anchor="w")
        self.tabla_resumen.pack(fill=tk.X)

        marco_tabla = ttk.Frame(f)
        marco_tabla.pack(fill=tk.BOTH, expand=True)
        cols = ("idx", "area", "amplitud", "fwhm", "t_onset", "bandera")
        self.tabla = ttk.Treeview(marco_tabla, columns=cols, show="headings", height=20)
        titulos = {"idx": "#", "area": "area [mV·ns]", "amplitud": "amplitud [mV]",
                  "fwhm": "FWHM [ns]", "t_onset": "t_onset [ns]", "bandera": "bandera"}
        anchos = {"idx": 50, "area": 110, "amplitud": 100, "fwhm": 90,
                 "t_onset": 90, "bandera": 200}
        for c in cols:
            self.tabla.heading(c, text=titulos[c])
            self.tabla.column(c, width=anchos[c], anchor="center")
        self.tabla.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        scroll = ttk.Scrollbar(marco_tabla, orient=tk.VERTICAL, command=self.tabla.yview)
        scroll.pack(side=tk.LEFT, fill=tk.Y)
        self.tabla.configure(yscrollcommand=scroll.set)
        self.tabla.bind("<<TreeviewSelect>>", self._al_elegir_fila)

    # ------------------------------------------------------ tab Comparacion --
    def _construir_tab_comparacion(self):
        f = ttk.Frame(self.notebook, padding=6)
        self.notebook.add(f, text="Comparacion")

        ctrl = ttk.Frame(f)
        ctrl.pack(fill=tk.X, pady=(0, 4))

        ttk.Label(ctrl, text="Superponer:").pack(side=tk.LEFT)
        self.var_n_comp = tk.IntVar(value=60)
        ttk.Spinbox(ctrl, from_=5, to=200, increment=5, textvariable=self.var_n_comp,
                   width=6).pack(side=tk.LEFT, padx=(4, 12))

        self.var_alinear = tk.BooleanVar(value=True)
        ttk.Checkbutton(ctrl, text="alinear por onset", variable=self.var_alinear
                       ).pack(side=tk.LEFT, padx=(0, 8))
        self.var_normalizar = tk.BooleanVar(value=False)
        ttk.Checkbutton(ctrl, text="normalizar por amplitud",
                       variable=self.var_normalizar).pack(side=tk.LEFT, padx=(0, 8))

        ttk.Button(ctrl, text="Actualizar", command=self._redibujar_comparacion
                  ).pack(side=tk.LEFT)

        self.lienzo_comparacion = Lienzo(f, dpi=self.dpi_fig)
        self.lienzo_comparacion.pack(fill=tk.BOTH, expand=True)

    # ------------------------------------------------------------ archivo --
    def _al_elegir_corrida(self):
        clave = self.var_corrida.get()
        if clave == ARCHIVO_SUELTO:
            self.entry_ruta.configure(state="normal")
            self.btn_elegir.configure(state="normal")
            self.var_atenuador.set("0")
        else:
            c = rutas.corrida(clave)
            self.var_ruta.set(rutas.relativa(c.muestra))
            self.entry_ruta.configure(state="disabled")
            self.btn_elegir.configure(state="disabled")
            # el catalogo (analisis/rutas.py) trae el atenuador SOLO para las
            # corridas donde se anoto (co60 = 20 dB); si no esta registrado
            # se deja en 0 y es responsabilidad de quien mide completarlo.
            atenuador = getattr(c, "atenuador", None)
            self.var_atenuador.set("%g" % atenuador if atenuador else "0")

    def _elegir_archivo(self):
        # "Todos" va PRIMERO: los binarios del DRS4 de la sesion del 1-sep-2026
        # no tienen extension en el nombre (`Co60`, `Cs137`, `B1700V`), asi que
        # con el filtro "*.xml" la carpeta se ve vacia. `analisis.drs` decide el
        # formato por los primeros bytes, no por el nombre.
        ruta = filedialog.askopenfilename(
            parent=self, title="Elegir el archivo de la muestra (XML o binario)",
            filetypes=[("Datos del DRS4 (XML o binario)", "*.*"),
                       ("XML del DRSOSC", "*.xml")])
        if ruta:
            self.var_ruta.set(ruta)

    def _ruta_elegida(self):
        clave = self.var_corrida.get()
        if clave == ARCHIVO_SUELTO:
            ruta = self.var_ruta.get().strip()
            return ruta if ruta else None
        return rutas.corrida(clave).muestra

    # -------------------------------------------------------------- carga --
    def _cargar_pulsos(self):
        ruta = self._ruta_elegida()
        if not ruta or not os.path.exists(ruta):
            self.var_estado.set("elegi un .xml valido primero")
            return
        n = int(self.var_n.get())
        self.btn_cargar.configure(state="disabled")
        self.var_estado.set("leyendo %s..." % os.path.basename(ruta))
        self.barra.configure(mode="determinate", value=0)

        def callback(frac, n_ev):
            self.cola.put(dict(tipo="avance_carga", frac=frac, n=n_ev))

        def trabajo():
            try:
                T, V, relojes, scalers, canal = lote.cargar_pulsos(
                    ruta, n=n, callback=callback)
                self.cola.put(dict(tipo="cargado", ruta=ruta, T=T, V=V,
                                   relojes=relojes, scalers=scalers, canal=canal))
            except Exception as e:
                self.cola.put(dict(tipo="error", etapa="carga", error=str(e)))

        threading.Thread(target=trabajo, daemon=True).start()

    def _calcular_lote(self):
        if self.T is None:
            self.var_estado.set("carga pulsos primero")
            return
        if self.pre is None or self.post is None:
            self.var_estado.set("define el intervalo en la pestana Intervalo")
            return
        try:
            atenuador_db = float(self.var_atenuador.get().replace(",", ".") or 0.0)
        except ValueError:
            self.var_estado.set("atenuador invalido, uso 0 dB")
            atenuador_db = 0.0
        self.btn_cargar.configure(state="disabled")
        self.var_estado.set("calculando area y FWHM de %d pulsos..." % len(self.T))
        self.barra.configure(mode="determinate", value=0)
        pre, post = self.pre, self.post

        def callback(frac, n_ev):
            self.cola.put(dict(tipo="avance_calculo", frac=frac, n=n_ev))

        def trabajo():
            try:
                d = lote.calcular_lote(self.T, self.V, pre, post, callback=callback,
                                       atenuador_db=atenuador_db)
                self.cola.put(dict(tipo="calculado", datos=d))
            except Exception as e:
                self.cola.put(dict(tipo="error", etapa="calculo", error=str(e)))

        threading.Thread(target=trabajo, daemon=True).start()

    # -------------------------------------------------------------- cola --
    def _vaciar_cola(self):
        try:
            while True:
                m = self.cola.get_nowait()
                self._procesar_mensaje(m)
        except queue.Empty:
            pass
        self.after(80, self._vaciar_cola)

    def _procesar_mensaje(self, m):
        tipo = m["tipo"]
        if tipo == "avance_carga":
            self.barra["value"] = 100.0 * m["frac"]
            self.var_estado.set("leyendo... %d eventos vistos" % m["n"])
        elif tipo == "avance_calculo":
            self.barra["value"] = 100.0 * m["frac"]
            self.var_estado.set("calculando... %d/%d" % (m["n"], len(self.T)))
        elif tipo == "cargado":
            self.btn_cargar.configure(state="normal")
            self.T, self.V = m["T"], m["V"]
            self.relojes, self.scalers = m["relojes"], m["scalers"]
            self.canal, self.ruta_actual = m["canal"], m["ruta"]
            self.datos_lote = None
            self.spin_idx.configure(to=len(self.T) - 1)
            self.var_idx.set(0)
            self.var_estado.set("%d pulsos leidos (canal %s)"
                               % (len(self.T), self.canal))
            self.barra["value"] = 100.0
            self._sugerir_ventana()
        elif tipo == "calculado":
            self.btn_cargar.configure(state="normal")
            self.datos_lote = m["datos"]
            self.var_estado.set("listo: %d pulsos calculados" % len(self.T))
            self.barra["value"] = 100.0
            self._llenar_tabla()
            self._redibujar_comparacion()
        elif tipo == "error":
            self.btn_cargar.configure(state="normal")
            self.var_estado.set("error en %s: %s" % (m["etapa"], m["error"]))

    # --------------------------------------------------------- ventana --
    def _sugerir_ventana(self):
        if self.T is None:
            return
        try:
            self.vent = lote.sugerir_ventana(self.T, self.V)
        except Exception as e:
            self.var_estado.set("no pude sugerir ventana: %s" % e)
            return
        self.pre, self.post = self.vent["pre"], self.vent["post"]
        self._redibujar_intervalo()

    def _al_seleccionar_span(self, xmin, xmax):
        self.pre = max(-float(xmin), 0.0)
        self.post = max(float(xmax), 0.0)
        self.var_ancho.set("pre = %.1f ns   post = %.1f ns   ancho = %.1f ns"
                          % (self.pre, self.post, self.pre + self.post))

    # ---------------------------------------------------- redibujar tabs --
    def _al_cambiar_tab(self, evento):
        actual = self.notebook.tab(self.notebook.select(), "text")
        if actual == "Intervalo":
            self._redibujar_intervalo()
        elif actual == "Comparacion":
            self._redibujar_comparacion()

    def _pulso_alineado(self, i):
        """`(tau, d, t_onset)` del pulso `i`: mismo criterio de localizacion
        que `pulso.area_evento`/`ancho.ancho_media_altura` (preparar, pico,
        onset, linea base), reusado aca solo para DIBUJAR."""
        t, v = self.T[i], self.V[i]
        tau, v = pulso.preparar(t, v)
        base0 = float(np.median(v))
        i_p, _, amp0, signo = pulso.detectar_pico(tau, v, base0)
        i_on, t_on = pulso.detectar_onset(tau, v, base0, i_p, signo)
        base, ruido, n_base = pulso.linea_base(tau, v, hasta=t_on - 4.0)
        if n_base < 10:
            base = base0
        i_p, _, amp, signo = pulso.detectar_pico(tau, v, base)
        i_on, t_on = pulso.detectar_onset(tau, v, base, i_p, signo)
        return tau - t_on, signo * (v - base), t_on

    def _redibujar_intervalo(self):
        if self.T is None:
            self.lienzo_intervalo.vacio("todavia no hay pulsos cargados")
            return

        def fn(fig):
            ax = fig.add_subplot(111)
            if self.var_vista.get() == "Forma promedio" and self.vent is not None:
                ax.plot(self.vent["tau_rel"], self.vent["forma"], lw=1.4, color="C0")
                ax.set_title("forma promedio alineada por onset (%d pulsos)"
                            % self.vent["n_usados"])
                ax.set_ylabel("forma normalizada")
            else:
                i = int(np.clip(self.var_idx.get(), 0, len(self.T) - 1))
                d_t, d_v, _ = self._pulso_alineado(i)
                ax.plot(d_t, d_v, lw=1.1, color="C0")
                ax.set_title("pulso #%d" % i)
                ax.set_ylabel("v - base, orientado [mV]")
            ax.axhline(0, color="k", lw=0.6)
            ax.set_xlabel("t - t_onset [ns]")
            ax.grid(alpha=0.3)

            self._span = SpanSelector(
                ax, self._al_seleccionar_span, "horizontal",
                useblit=True, interactive=True, drag_from_anywhere=True,
                props=dict(alpha=0.25, facecolor="C0"))
            if self.pre is not None and self.post is not None:
                self._span.extents = (-self.pre, self.post)
                self.var_ancho.set("pre = %.1f ns   post = %.1f ns   ancho = %.1f ns"
                                  % (self.pre, self.post, self.pre + self.post))

        self.lienzo_intervalo.dibujar(fn)

    def _redibujar_comparacion(self):
        if self.T is None:
            self.lienzo_comparacion.vacio("todavia no hay pulsos cargados")
            return
        n = int(np.clip(self.var_n_comp.get(), 1, len(self.T)))
        alinear = self.var_alinear.get()
        normalizar = self.var_normalizar.get()
        datos = self.datos_lote

        def fn(fig):
            if datos is not None:
                ax1 = fig.add_subplot(1, 3, 1)
                ax2 = fig.add_subplot(1, 3, 2)
                ax3 = fig.add_subplot(1, 3, 3)
            else:
                ax1 = fig.add_subplot(1, 1, 1)
                ax2 = ax3 = None

            for i in range(n):
                if alinear:
                    d_t, d_v, _ = self._pulso_alineado(i)
                else:
                    tau, v = pulso.preparar(self.T[i], self.V[i])
                    d_t, d_v = tau, v - float(np.median(v))
                if normalizar:
                    pico = np.max(np.abs(d_v))
                    if pico > 0:
                        d_v = d_v / pico
                ax1.plot(d_t, d_v, lw=0.5, alpha=0.35)
            if self.pre is not None and self.post is not None and alinear:
                ax1.axvspan(-self.pre, self.post, color="C2", alpha=0.12,
                          label="ventana aplicada")
                ax1.legend(fontsize=8)
            ax1.set_title("%d pulsos superpuestos" % n)
            ax1.set_xlabel("t - t_onset [ns]" if alinear else "t [ns]")
            ax1.set_ylabel("forma normalizada" if normalizar else "v - base [mV]")
            ax1.grid(alpha=0.3)

            if datos is not None:
                bandera = datos["bandera_area"] | datos["bandera_fwhm"]
                ok = bandera == 0
                ax2.hist(datos["fwhm"][ok], bins=40, color="C0")
                ax2.set_title("FWHM (%d pulsos limpios)" % int(ok.sum()))
                ax2.set_xlabel("FWHM [ns]"); ax2.set_ylabel("cuentas")
                ax2.grid(alpha=0.3)

                ax3.hist(datos["area"][ok], bins=40, color="C1")
                ax3.set_title("area")
                ax3.set_xlabel("area [mV·ns]"); ax3.set_ylabel("cuentas")
                ax3.grid(alpha=0.3)

        self.lienzo_comparacion.dibujar(fn)

    # --------------------------------------------------------------- tabla --
    def _llenar_tabla(self):
        self.tabla.delete(*self.tabla.get_children())
        d = self.datos_lote
        n = len(d["area"])
        factor = d.get("factor_atenuacion", 1.0)
        db = d.get("atenuador_db", 0.0)
        sufijo = " (corr. ×%.3g)" % factor if factor != 1.0 else ""
        self.tabla.heading("area", text="area [mV·ns]%s" % sufijo)
        self.tabla.heading("amplitud", text="amplitud [mV]%s" % sufijo)
        for i in range(n):
            bits = int(d["bandera_area"][i] | d["bandera_fwhm"][i])
            self.tabla.insert("", "end", iid=str(i), values=(
                i, "%.1f" % d["area"][i], "%.1f" % d["amplitud"][i],
                "%.2f" % d["fwhm"][i], "%.1f" % d["t_onset"][i],
                _decodificar_bandera(bits)))
        ok = (d["bandera_area"] == 0) & (d["bandera_fwhm"] == 0)
        n_ok = int(ok.sum())
        self.var_resumen.set(
            "%d pulsos · %d sin bandera · ventana pre=%.1f post=%.1f ns"
            % (n, n_ok, d["pre"], d["post"]))

        self.tabla_resumen.delete(*self.tabla_resumen.get_children())
        if n_ok > 0:
            area_ok, fwhm_ok = d["area"][ok], d["fwhm"][ok]
            m_area, disp_area, sem_area = _media_dispersion_sem(area_ok)
            m_fwhm, disp_fwhm, sem_fwhm = _media_dispersion_sem(fwhm_ok)
            filas = [
                ("pulsos leidos", "%d" % n),
                ("sin ninguna bandera", "%d (%.0f%%)" % (n_ok, 100.0 * n_ok / n)),
                ("intervalo de integracion",
                 "[onset − %.1f, onset + %.1f] ns  (ancho %.1f ns, igual para "
                 "todos los pulsos)" % (d["pre"], d["post"], d["pre"] + d["post"])),
                ("atenuador corregido",
                 "%.0f dB → area/amplitud ×%.3g (voltaje real a la salida del "
                 "detector, antes del atenuador)" % (db, factor) if factor != 1.0
                 else "0 dB — valores tal como entran al digitalizador, sin corregir"),
                ("area — mediana", "%.1f mV·ns" % np.median(area_ok)),
                ("area — dispersion (std, pulso a pulso)", "%.1f mV·ns" % disp_area),
                ("area — media ± error estandar", "%.1f ± %.2f mV·ns" % (m_area, sem_area)),
                ("FWHM — mediana", "%.2f ns" % np.median(fwhm_ok)),
                ("FWHM — dispersion (std, pulso a pulso)", "%.2f ns" % disp_fwhm),
                ("FWHM — media ± error estandar", "%.2f ± %.3f ns" % (m_fwhm, sem_fwhm)),
            ]
        else:
            filas = [("pulsos leidos", "%d" % n),
                     ("sin ninguna bandera", "0 — nada limpio para resumir")]
        for campo, valor in filas:
            self.tabla_resumen.insert("", "end", values=(campo, valor))

    def _al_elegir_fila(self, evento):
        sel = self.tabla.selection()
        if not sel:
            return
        i = int(sel[0])
        self.var_idx.set(i)
        self.var_vista.set("Pulso individual")
        self.notebook.select(0)
        self._redibujar_intervalo()


def main():
    activar_dpi_awareness()
    App().mainloop()


if __name__ == "__main__":
    main()
