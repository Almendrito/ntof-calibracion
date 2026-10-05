"""
visor.py — interfaz para VER la calibracion ya medida.

    python analisis/visor.py

No calcula calibraciones: eso lo hace `gui.py`, que corre la cadena y guarda en
el registro. Este visor lee lo que ya esta guardado y lo pone junto, que es lo
que hace falta cuando hay 16 corridas, 4 detectores, 3 voltajes y 2 fuentes y la
pregunta ya no es "cuanto da esta corrida" sino "cierran entre si".

Cuatro pestanas, una por pregunta:

  1. Calibracion (TEST A)   canal contra energia, por punto de operacion. Con
                            dos fuentes se puede liberar la ordenada; el chequeo
                            de pedestal-contra-ganancia mira TODOS los puntos a
                            la vez, porque con uno solo las dos hipotesis son
                            indistinguibles.
  2. Ganancia (TEST B)      D contra alto voltaje, por detector y fuente. Un PMT
                            da D ~ V^b con b entre 6 y 10; si sale fuera de ese
                            rango, entre las corridas cambio algo mas que el
                            voltaje.
  3. Espectros              el espectro neto de cada corrida con los bordes que
                            el ajuste puso, para mirar si el ajuste cae donde
                            uno esperaria.
  4. Archivos               que hay en el catalogo: formato (XML o binario),
                            canal, eventos, tiempo VIVO y -lo que rompio la
                            sesion del 01-sep- si la compuerta de integracion
                            entra en el registro.

Los calculos que tardan (leer pulsos, armar espectros) corren en un hilo aparte
y se comunican por una cola. El hilo trabajador nunca toca widgets.
"""

import os
import queue
import sys
import threading
import traceback

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# gui.py declara el DPI awareness al importarse, ANTES de crear el root de Tk,
# y ademas trae el lienzo de matplotlib y el capturador de informes. Se reusan
# en vez de copiarlos.
from analisis.gui import Lienzo, texto_de              # noqa: E402

import tkinter as tk                                    # noqa: E402
from tkinter import ttk                                 # noqa: E402

from analisis import (area, drs, espectro, fisica, registro, rutas,  # noqa: E402
                      voltaje)

FUENTE_MONO = ("Consolas", 9)
GRIS = "#5b6b7c"
ROJO = "#b32b2b"


def _texto(padre, **kw):
    """Un Text monoespaciado con scroll, que es como se leen estos informes."""
    cont = ttk.Frame(padre)
    t = tk.Text(cont, wrap="none", font=FUENTE_MONO, relief="flat",
                background="#fbfcfd", **kw)
    sy = ttk.Scrollbar(cont, orient="vertical", command=t.yview)
    sx = ttk.Scrollbar(cont, orient="horizontal", command=t.xview)
    t.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
    sy.pack(side=tk.RIGHT, fill=tk.Y)
    sx.pack(side=tk.BOTTOM, fill=tk.X)
    t.pack(fill=tk.BOTH, expand=True)
    return cont, t


def _poner(t, s):
    t.configure(state="normal")
    t.delete("1.0", tk.END)
    t.insert(tk.END, s)
    t.configure(state="disabled")


# ==========================================================================
class Visor(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Calibracion por borde Compton - visor")
        self.geometry("1500x900")
        self.cola = queue.Queue()
        self.cache_espectros = {}

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=(8, 0))
        self._tab_calibracion()
        self._tab_ganancia()
        self._tab_espectros()
        self._tab_archivos()

        f = ttk.Frame(self, padding=(10, 4, 10, 8))
        f.pack(fill=tk.X)
        self.v_estado = tk.StringVar(value="listo")
        ttk.Label(f, textvariable=self.v_estado, foreground=GRIS).pack(side=tk.LEFT)
        self.barra = ttk.Progressbar(f, mode="indeterminate", length=180)
        self.barra.pack(side=tk.RIGHT)

        self.after(80, self._bombear)
        self.refrescar_calibracion()
        self.refrescar_ganancia()

    # ---------------------------------------------------------- infraestructura
    def _trabajar(self, fn, cuando_termine, texto="calculando..."):
        """Corre `fn()` en un hilo y entrega el resultado por la cola."""
        self.v_estado.set(texto)
        self.barra.start(12)

        def correr():
            try:
                self.cola.put(("ok", cuando_termine, fn()))
            except Exception as e:                       # noqa: BLE001
                self.cola.put(("error", cuando_termine,
                               "%s\n\n%s" % (e, traceback.format_exc())))

        threading.Thread(target=correr, daemon=True).start()

    def _bombear(self):
        try:
            while True:
                estado, destino, dato = self.cola.get_nowait()
                self.barra.stop()
                self.v_estado.set("listo" if estado == "ok" else "error")
                destino(dato) if estado == "ok" else self._error(dato)
        except queue.Empty:
            pass
        self.after(80, self._bombear)

    def _error(self, msg):
        for t in (self.txt_esp, self.txt_arch):
            try:
                _poner(t, msg)
            except Exception:
                pass

    # ============================================== 1. calibracion (test A) ===
    def _tab_calibracion(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="1 · Calibracion (test A)")

        arriba = ttk.Frame(tab)
        arriba.pack(fill=tk.X, padx=6, pady=6)
        ttk.Button(arriba, text="Actualizar",
                   command=self.refrescar_calibracion).pack(side=tk.LEFT)
        self.v_escala = tk.StringVar(value="keVee")
        for e in ("keVee", "keV"):
            ttk.Radiobutton(arriba, text=e, value=e, variable=self.v_escala,
                            command=self.refrescar_calibracion).pack(side=tk.LEFT,
                                                                     padx=4)
        ttk.Label(arriba, foreground=GRIS,
                  text="   Test A: DOS FUENTES en el mismo punto de operacion. "
                       "Con una sola, el mapa por el origen hace que los bordes "
                       "caigan exactos sobre la recta y el test no dice nada."
                  ).pack(side=tk.LEFT)

        cols = ("punto", "fuente", "D", "err", "borde", "chi2/ndf", "sigma/E",
                "corrida", "fecha")
        self.tv_cal = ttk.Treeview(tab, columns=cols, show="headings", height=9)
        for c, a in zip(cols, (150, 75, 90, 80, 100, 75, 75, 130, 130)):
            self.tv_cal.heading(c, text=c)
            self.tv_cal.column(c, width=a, anchor="w")
        self.tv_cal.pack(fill=tk.X, padx=6)
        self.tv_cal.bind("<<TreeviewSelect>>", lambda _e: self._dibujar_test_a())

        pw = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        cont, self.txt_cal = _texto(pw, width=66)
        pw.add(cont, weight=0)
        self.lz_cal = Lienzo(pw)
        self.lz_cal.nombre = 'test A'
        pw.add(self.lz_cal, weight=1)

    def refrescar_calibracion(self):
        esc = self.v_escala.get()
        reg = registro.cargar()
        self.grupos_cal = registro.por_punto(reg, escala_pedida=esc,
                                             metodo="global")
        for i in self.tv_cal.get_children():
            self.tv_cal.delete(i)
        for punto, entradas in sorted(self.grupos_cal.items()):
            for k, r in entradas:
                self.tv_cal.insert(
                    "", tk.END, iid="%s||%s" % (punto, k),
                    values=(punto, r["fuente"], "%.4f" % r["D"],
                            "%.4f" % r["err_D"],
                            "%.1f" % (r["D"] * r["E_ref"]),
                            "%.2f" % r["chi2_ndf"],
                            "%.4f" % r.get("sigma_rel", float("nan")),
                            r.get("corrida", ""), r.get("guardado", "")))
        _poner(self.txt_cal, texto_de(registro.informe, escala_pedida=esc,
                                      metodo="global"))
        self._dibujar_test_a()

    def _dibujar_test_a(self):
        sel = self.tv_cal.selection()
        punto = sel[0].split("||")[0] if sel else None
        grupos = getattr(self, "grupos_cal", {})
        if not grupos:
            self.lz_cal.vacio("el registro esta vacio: corre `gui.py` o "
                              "`pipeline.py` primero")
            return
        if punto is None or punto not in grupos:
            punto = sorted(grupos)[0]
        entradas = grupos[punto]
        lin = registro.linealidad(entradas)

        def dibujar(fig):
            a = fig.add_subplot(111)
            E = np.array([r["E_ref"] for _, r in entradas], float)
            ch = np.array([r["D"] * r["E_ref"] for _, r in entradas], float)
            er = np.array([r["err_D"] * r["E_ref"] for _, r in entradas], float)
            fu = [r["fuente"] for _, r in entradas]
            a.errorbar(E, ch, yerr=np.maximum(er, 1e-9), fmt="o", ms=7,
                       color="tab:blue", capsize=3, zorder=3)
            for x, y, f in zip(E, ch, fu):
                a.annotate("  %s" % f, (x, y), fontsize=9, color=GRIS)
            x = np.linspace(0, max(E.max() * 1.15, 1.0), 50)
            if lin is not None:
                a.plot(x, lin["D_origen"] * x, "-", color="tab:green", lw=1.4,
                       label="por el origen: D = %.4f" % lin["D_origen"])
                a.plot(x, lin["ordenada"] + lin["D_libre"] * x, "--",
                       color="tab:red", lw=1.4,
                       label="ordenada libre: D = %.4f, c0 = %+.1f"
                             % (lin["D_libre"], lin["ordenada"]))
                a.legend(fontsize=9)
            else:
                a.plot(x, (ch[0] / E[0]) * x, "-", color="tab:green", lw=1.2,
                       label="por el origen (una sola fuente: sin test)")
                a.legend(fontsize=9)
            a.set_xlabel("energia del borde Compton [%s]" % self.v_escala.get())
            a.set_ylabel("canal [mV*ns]")
            a.set_title("Test A - %s" % punto)
            a.grid(alpha=0.3)
            a.set_xlim(left=0)
            a.set_ylim(bottom=0)

        self.lz_cal.dibujar(dibujar)

    # ================================================= 2. ganancia (test B) ===
    def _tab_ganancia(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="2 · Ganancia D(V) (test B)")
        arriba = ttk.Frame(tab)
        arriba.pack(fill=tk.X, padx=6, pady=6)
        ttk.Button(arriba, text="Actualizar",
                   command=self.refrescar_ganancia).pack(side=tk.LEFT)
        ttk.Label(arriba, text="  serie:").pack(side=tk.LEFT)
        self.cb_serie = ttk.Combobox(arriba, state="readonly", width=34)
        self.cb_serie.pack(side=tk.LEFT, padx=4)
        self.cb_serie.bind("<<ComboboxSelected>>",
                           lambda _e: self._dibujar_test_b())
        ttk.Label(arriba, foreground=GRIS,
                  text="   Test B: la MISMA fuente a varios HV. Un PMT da "
                       "D ~ V^b con b entre 6 y 10; la respuesta esperada NO "
                       "es una recta.").pack(side=tk.LEFT)

        pw = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        cont, self.txt_gan = _texto(pw, width=66)
        pw.add(cont, weight=0)
        self.lz_gan = Lienzo(pw)
        self.lz_gan.nombre = 'test B'
        pw.add(self.lz_gan, weight=1)

    def refrescar_ganancia(self):
        esc = self.v_escala.get() if hasattr(self, "v_escala") else "keVee"
        self.series = voltaje.puntos_del_registro(escala_pedida=esc,
                                                  metodo="global")
        nombres = ["%s / %s" % (d, f) for (d, f) in sorted(self.series)]
        self.cb_serie["values"] = nombres
        if nombres and self.cb_serie.get() not in nombres:
            self.cb_serie.set(nombres[0])
        if not nombres:
            self.cb_serie.set("")
        self._dibujar_test_b()

    def _dibujar_test_b(self):
        nombre = self.cb_serie.get()
        series = getattr(self, "series", {})
        if not nombre or not series:
            self.lz_gan.vacio("hacen falta al menos 2 voltajes de la misma "
                              "fuente y detector en el registro")
            _poner(self.txt_gan,
                   "El test B necesita la MISMA fuente medida a varios "
                   "altos voltajes, con el mismo detector.")
            return
        det, fue = [x.strip() for x in nombre.split("/")]
        d = series[(det, fue)]
        V, D, eD = d["V"], d["D"], d["err"]
        if V.size < 2:
            self.lz_gan.vacio("%s: un solo voltaje" % nombre)
            _poner(self.txt_gan, "Con un solo voltaje no hay ley de ganancia.")
            return
        res = voltaje.ajustar(V, D, eD)
        txt = texto_de(voltaje.informe, res, d["etiquetas"])
        avisos = voltaje.avisos_de_cadena([c for c in d["claves"] if c])
        if avisos:
            partes = ["", "  AVISOS DE CADENA DE SENAL"]
            partes += ["     %s" % a for a in avisos]
            txt += os.linesep.join(partes) + os.linesep
        _poner(self.txt_gan, txt)

        def dibujar(fig):
            a1 = fig.add_subplot(211)
            a2 = fig.add_subplot(212, sharex=a1)
            x = np.linspace(V.min() * 0.99, V.max() * 1.01, 200)
            a1.errorbar(V, D, yerr=eD, fmt="o", ms=7, color="tab:blue",
                        capsize=3, zorder=3, label="medido")
            a1.plot(x, voltaje.evaluar(res, x), "-", color="tab:red", lw=1.5,
                    label="potencia  b = %.2f +- %.2f   (chi2 %.1f / %d gdl)"
                          % (res["b"], res["err_b"], res["chi2_potencia"],
                             res["gdl"]))
            a1.plot(x, res["ordenada"] + res["pendiente"] * x, "--",
                    color="tab:green", lw=1.2,
                    label="recta  (chi2 %.1f)" % res["chi2_recta"])
            a1.set_ylabel("D [mV*ns/%s]" % self.v_escala.get())
            a1.set_title("Test B - %s%s" % (nombre,
                         "" if res["b_plausible"] else
                         "   [b FUERA del rango 6-10 de un PMT]"))
            a1.legend(fontsize=9)
            a1.grid(alpha=0.3)

            a2.plot(V, 100 * res["residuos_potencia"] / D, "o-",
                    color="tab:red", ms=5, lw=1, label="potencia")
            a2.plot(V, 100 * res["residuos_recta"] / D, "s--",
                    color="tab:green", ms=5, lw=1, label="recta")
            a2.axhline(0, color="0.4", lw=0.8)
            a2.set_xlabel("alto voltaje [V]")
            a2.set_ylabel("residuo [%]")
            a2.legend(fontsize=9)
            a2.grid(alpha=0.3)
            fig.tight_layout()

        self.lz_gan.dibujar(dibujar)


    # ======================================================= 3. espectros ====
    def _tab_espectros(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="3 · Espectros")
        arriba = ttk.Frame(tab)
        arriba.pack(fill=tk.X, padx=6, pady=6)
        ttk.Label(arriba, text="corrida:").pack(side=tk.LEFT)
        self.cb_corrida = ttk.Combobox(arriba, state="readonly", width=22,
                                       values=sorted(rutas.disponibles()))
        self.cb_corrida.pack(side=tk.LEFT, padx=4)
        if self.cb_corrida["values"]:
            self.cb_corrida.set(self.cb_corrida["values"][0])
        self.cb_corrida.bind("<<ComboboxSelected>>",
                             lambda _e: self.mostrar_espectro())
        ttk.Button(arriba, text="Mostrar",
                   command=self.mostrar_espectro).pack(side=tk.LEFT, padx=4)
        self.v_log = tk.BooleanVar(value=False)
        ttk.Checkbutton(arriba, text="log", variable=self.v_log,
                        command=self.mostrar_espectro).pack(side=tk.LEFT, padx=8)
        ttk.Label(arriba, foreground=GRIS,
                  text="   si la corrida no tiene el espectro cacheado hay que "
                       "parsearla (~20 s el binario, ~90 s el XML grande)"
                  ).pack(side=tk.LEFT)

        pw = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        cont, self.txt_esp = _texto(pw, width=52)
        pw.add(cont, weight=0)
        self.lz_esp = Lienzo(pw)
        self.lz_esp.nombre = 'espectros'
        pw.add(self.lz_esp, weight=1)

    def mostrar_espectro(self):
        clave = self.cb_corrida.get()
        if not clave:
            return
        if clave in self.cache_espectros:
            self._pintar_espectro((clave, self.cache_espectros[clave]))
            return

        def calcular():
            try:
                e = espectro.cargar("espectro_%s_global" % clave)
            except Exception:                            # noqa: BLE001
                c = rutas.corrida(clave)
                v = area.calibrar_ventana(c.muestra, verbose=False)
                d_m = area.areas_de_archivo(c.muestra, pre=v["pre"],
                                            post=v["post"], verbose=False)
                d_f = (area.areas_de_archivo(c.fondo, pre=v["pre"],
                                             post=v["post"], verbose=False)
                       if c.fondo and os.path.exists(c.fondo) else None)
                e = espectro.espectro_neto(d_m, d_f, bins=600, verbose=False)
                e["ventana"] = v
            return clave, e

        self._trabajar(calcular, self._pintar_espectro,
                       "armando el espectro de %s..." % clave)

    def _pintar_espectro(self, dato):
        clave, e = dato
        self.cache_espectros[clave] = e
        c = rutas.corrida(clave)
        r = self._entrada_registro(clave)

        lineas = ["corrida : %s" % clave, "        : %s" % c.descripcion,
                  "punto   : %s" % c.punto(),
                  "cadena  : %s" % c.cadena(),
                  "formato : %s" % drs.formato(c.muestra),
                  "",
                  "fondo   : %s" % ("restado" if e.get("hay_fondo")
                                    else "NO HAY (espectro crudo)"),
                  "pulsos  : %s" % format(int(e.get("n_muestra", 0)), ",d"),
                  "vivo    : %.1f%% (%.0f s de %.0f de pared)"
                  % (100 * e.get("frac_viva_muestra", float("nan")),
                     e.get("t_vivo_muestra", 0), e.get("t_pared_muestra", 0)),
                  "binning : %d bins hasta %.0f mV*ns"
                  % (e.get("bins", 0), e.get("rango_max", 0))]
        if r:
            lineas += ["", "CALIBRACION GUARDADA",
                       "  D        = %.4f +- %.4f mV*ns/%s"
                       % (r["D"], r["err_D"], "keVee" if "keVee" in
                          str(r.get("rama", "")) or r.get("escala") == "keVee"
                          else "keV"),
                       "  sigma/E  = %.4f  -> FWHM %.1f%%"
                       % (r.get("sigma_rel", float("nan")),
                          100 * r.get("fwhm_rel", float("nan"))),
                       "  chi2/ndf = %.3f" % r["chi2_ndf"],
                       "  ventana  = [%.0f, %.0f] mV*ns"
                       % tuple(r.get("ventana", (0, 0)))]
            for E, ch in zip(r.get("bordes_E", []), r.get("bordes_ch", [])):
                lineas.append("  borde %7.1f keV -> %8.1f mV*ns" % (E, ch))
        else:
            lineas += ["", "sin calibracion en el registro para esta corrida"]
        _poner(self.txt_esp, "\n".join(lineas))

        def dibujar(fig):
            a = fig.add_subplot(111)
            a.step(e["canales"], e["neto"], where="mid", lw=0.8, color="0.35",
                   label="neto")
            if r:
                lo, hi = r.get("ventana", (None, None))
                if lo is not None:
                    a.axvspan(lo, hi, color="tab:orange", alpha=0.12,
                              label="ventana de ajuste")
                for E, ch in zip(r.get("bordes_E", []), r.get("bordes_ch", [])):
                    a.axvline(ch, color="tab:red", lw=1.3)
                    a.annotate("%.0f keV" % E, (ch, 0.92), xycoords=("data",
                                                                     "axes fraction"),
                               rotation=90, fontsize=8, color="tab:red",
                               ha="right", va="top")
            if self.v_log.get():
                a.set_yscale("log")
                a.set_ylim(bottom=0.5)
            a.set_xlabel("area del pulso [mV*ns]")
            a.set_ylabel("cuentas netas")
            a.set_title("%s   (%s)" % (clave, c.punto()))
            a.grid(alpha=0.3)
            a.legend(fontsize=9)

        self.lz_esp.dibujar(dibujar)

    def _entrada_registro(self, clave):
        for _, r in registro.cargar().items():
            if r.get("corrida") == clave and r.get("metodo") == "global":
                return r
        return None

    # ======================================================== 4. archivos ====
    def _tab_archivos(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="4 · Archivos")
        arriba = ttk.Frame(tab)
        arriba.pack(fill=tk.X, padx=6, pady=6)
        ttk.Button(arriba, text="Revisar el catalogo",
                   command=self.revisar_archivos).pack(side=tk.LEFT)
        ttk.Label(arriba, foreground=GRIS,
                  text="   lee la cabecera y 600 pulsos de cada corrida "
                       "(~10 s en total). La columna 'entran' es la fraccion "
                       "de eventos cuya compuerta CABE en el registro."
                  ).pack(side=tk.LEFT)

        cols = ("corrida", "formato", "canal", "eventos", "post", "pedia",
                "carga", "entran", "muestreo", "archivo")
        self.tv_arch = ttk.Treeview(tab, columns=cols, show="headings", height=17)
        for c, a in zip(cols, (130, 90, 60, 80, 65, 65, 70, 70, 85, 330)):
            self.tv_arch.heading(c, text=c)
            self.tv_arch.column(c, width=a, anchor="w")
        self.tv_arch.pack(fill=tk.BOTH, expand=True, padx=6)
        self.tv_arch.tag_configure("recorte", foreground=ROJO)

        cont, self.txt_arch = _texto(tab, height=9)
        cont.pack(fill=tk.X, padx=6, pady=6)
        _poner(self.txt_arch,
               "Apreta 'Revisar el catalogo'.\n\n"
               "Lo que hay que mirar es la columna 'entran': si la compuerta de\n"
               "integracion no cabe en el registro, los eventos salen truncados\n"
               "y pierden una fraccion de carga que depende de CUANDO llego el\n"
               "pulso, que es justo el error que la compuerta anclada al onset\n"
               "existe para evitar. En la sesion del 01-sep-2026 el trigger\n"
               "quedo en t ~ 123 ns de un registro de 200 ns y el pulso necesita\n"
               "83 ns para converger: `area.calibrar_ventana` recorta el post y\n"
               "lo avisa.")

    def revisar_archivos(self):
        def calcular():
            filas = []
            for clave in sorted(rutas.disponibles()):
                c = rutas.corrida(clave)
                try:
                    info = drs.metadatos(c.muestra)
                    v = area.calibrar_ventana(c.muestra, canal=info["canal"],
                                              verbose=False)
                    filas.append(dict(
                        clave=clave, formato=info.get("formato", "?"),
                        canal=info["canal"],
                        eventos=info.get("n_eventos"),
                        post=v["post"], pedia=v["post_pedido"],
                        carga=v["fraccion"],
                        entran=v["frac_eventos_enteros"],
                        muestreo=info["muestreo_GSs"],
                        archivo=rutas.relativa(c.muestra),
                        recortado=v["recortado"]))
                except Exception as e:                   # noqa: BLE001
                    filas.append(dict(clave=clave, error=str(e)[:120]))
            return filas

        self._trabajar(calcular, self._pintar_archivos,
                       "leyendo las cabeceras del catalogo...")

    def _pintar_archivos(self, filas):
        for i in self.tv_arch.get_children():
            self.tv_arch.delete(i)
        recortadas = []
        for f in filas:
            if "error" in f:
                self.tv_arch.insert("", tk.END, values=(f["clave"], "ERROR", "",
                                                        "", "", "", "", "", "",
                                                        f["error"]))
                continue
            if f["recortado"]:
                recortadas.append(f)
            self.tv_arch.insert(
                "", tk.END, tags=("recorte",) if f["recortado"] else (),
                values=(f["clave"], f["formato"], f["canal"],
                        format(f["eventos"], ",d") if f["eventos"] else "",
                        "%.1f" % f["post"], "%.1f" % f["pedia"],
                        "%.1f%%" % (100 * f["carga"]),
                        "%.1f%%" % (100 * f["entran"]),
                        "%.2f GS/s" % f["muestreo"], f["archivo"]))
        ok = [f for f in filas if "error" not in f]
        s = ["%d corridas accesibles: %d en XML, %d en binario."
             % (len(ok), sum(1 for f in ok if f["formato"] == "XML del DRSOSC"),
                sum(1 for f in ok if str(f["formato"]).startswith("binario"))),
             ""]
        if recortadas:
            s += ["%d corridas necesitan RECORTAR la compuerta (en rojo):"
                  % len(recortadas)]
            for f in recortadas:
                s.append("   %-14s el pulso pedia %.1f ns, se uso %.1f -> "
                         "recoge %.1f%% de la carga y entra el %.1f%% de los "
                         "eventos" % (f["clave"], f["pedia"], f["post"],
                                      100 * f["carga"], 100 * f["entran"]))
            s += ["",
                  "Para COMPARAR dos corridas hay que usar la MISMA compuerta:",
                  "el `post` que mide cada archivo depende de su propio pulso, y",
                  "dos corridas con post distinto recogen fracciones distintas de",
                  "la carga, lo que entra derecho en D. Se fuerza con la opcion",
                  "`post_fijo` de `pipeline.calibrar_completo`."]
        else:
            s += ["Ninguna corrida necesita recortar la compuerta."]
        _poner(self.txt_arch, "\n".join(s))


def main():
    app = Visor()
    app.mainloop()


if __name__ == "__main__":
    main()
