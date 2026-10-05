"""
gui.py — el paper, paso a paso, calculando en vivo.

    python replica_paper/gui.py

Para que sirve: ver COMO se llega a cada numero de Kong et al. (IEEE TNS 2025),
ecuacion por ecuacion, con el razonamiento de cada paso al lado y los controles
para mover lo que el paper deja sin especificar.

Catorce pestanas, en el orden del propio paper:

    -  Mapa                     el diagrama de la Fig. 2, y donde vive cada ec.
    0  Datos                    donde ARRANCA el diagrama de la Fig. 2
    1  Ec. 1                    el borde Compton, impresa contra corregida
    2  Ec. 2 y 10               Klein-Nishina, y la Fig. 1 del paper
    3  Ec. 8-9                  lineshape -> canal + convolucion gaussiana
    4  Ec. 15                   NCC: seleccion automatica del intervalo
    5  Seccion II-A             recocido simulado con SUS constantes
    6  Fig. 2                   la cadena completa sobre una corrida
    7  Convencionales           los dos metodos contra los que compara
    8  Fig. 5                   robustez frente al intervalo
    9  Seccion II-C             las "diez mediciones independientes"
   10  Auditoria                los numeros publicados, contra si mismos
   11  Seccion III              la recta canal-energia
   12  Cierre                   verdad conocida, de punta a punta

Tres reglas de esta interfaz
----------------------------
1. NO CALCULA NADA POR SU CUENTA. Todo lo que muestra sale de llamar a los
   modulos de `replica_paper/`, que son los que replican el paper. Si un numero
   aparece en pantalla es porque una de esas funciones lo devolvio recien.

2. NO HAY CACHE DE RESULTADOS. Cada "Recalcular" rehace el calculo desde el
   espectro. Es mas lento y es a proposito: asi es imposible que quede en
   pantalla un numero de otros datos o de otros parametros. Lo unico que se
   reusa es el `.areas_gate.npz` que `analisis/area.py` guarda al lado de cada
   archivo de datos, que es el PARSEO del archivo (con clave = tamano + fecha),
   no un resultado del analisis.

3. LO QUE NO ESTA EN EL PAPER SE DICE. Cada pestana separa la cita literal del
   paper de las decisiones que hubo que tomar porque el paper no las escribe.
   Esas van marcadas [NO esta en el paper] y son, casi siempre, los controles.

Que se verifica y como: `python replica_paper/verificar_gui.py` corre las
catorce pestanas sin `mainloop` y despues compara sus numeros contra
`informe.txt`, que genero `correr.py`. Es lo que sostiene la regla 1.

El calculo corre en un hilo aparte y se comunica por una cola. El hilo
trabajador NUNCA toca widgets ni matplotlib: manda datos y el hilo de Tk dibuja.
"""

import contextlib
import io
import os
import queue
import sys
import threading
import time
import traceback
import warnings

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"


# --------------------------------------------------------------------------
# Nitidez en Windows: ANTES de crear el root de Tk
# --------------------------------------------------------------------------
def activar_dpi_awareness():
    """
    Declara el proceso consciente del DPI y devuelve True si se pudo.

    Sin esto Windows renderiza la ventana a 96 DPI y la ESTIRA como bitmap
    cuando la pantalla tiene escalado, y todo se ve borroso. Ninguna opcion de
    matplotlib lo arregla: hay que declararlo antes de que exista la ventana.
    Despues hay que escalar a mano, porque si no todo queda chico.
    """
    if sys.platform != "win32":
        return None
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(2)      # PER_MONITOR_AWARE
        return True
    except Exception:
        try:
            import ctypes
            ctypes.windll.user32.SetProcessDPIAware()
            return True
        except Exception:
            return False


activar_dpi_awareness()

import tkinter as tk                                        # noqa: E402
from tkinter import ttk                                     # noqa: E402

import matplotlib                                           # noqa: E402
matplotlib.use("TkAgg")
from matplotlib.figure import Figure                        # noqa: E402
from matplotlib.backends.backend_tkagg import (             # noqa: E402
    FigureCanvasTkAgg, NavigationToolbar2Tk)

from replica_paper import auditoria                         # noqa: E402
from replica_paper import convencionales as conv            # noqa: E402
from replica_paper import correr                            # noqa: E402
from replica_paper import datos as dat                      # noqa: E402
from replica_paper import ecuaciones as ec                  # noqa: E402
from replica_paper import metodo                            # noqa: E402
from replica_paper import ncc as mod_ncc                    # noqa: E402
from replica_paper import resultados as R                   # noqa: E402
from replica_paper import rutas                             # noqa: E402
from replica_paper import sa as mod_sa                      # noqa: E402

FONDO_TXT = "#f7f7f7"
COL = {"NCC-SA": "tab:red", "Gauss FWHM": "tab:green", "Diferencia": "tab:blue"}
LSS = {"NCC-SA": "-", "Gauss FWHM": "--", "Diferencia": ":"}


def texto_de(fn, *a, **kw):
    """Captura lo que una funcion de informe IMPRIME, sin duplicar su formato."""
    b = io.StringIO()
    with contextlib.redirect_stdout(b):
        fn(*a, **kw)
    return b.getvalue()


def formato_tiempo(segundos):
    s = max(0.0, float(segundos))
    if s < 90:
        return "%.0f s" % s
    if s < 3600:
        return "%d min %02d s" % (int(s // 60), int(s % 60))
    return "%d h %02d min" % (int(s // 3600), int((s % 3600) // 60))


def etiqueta_corrida(clave):
    """Texto del combo: clave, descripcion y si tiene fondo propio."""
    c = dat.corrida(clave)
    return "%s  -  %s%s" % (clave, c.descripcion,
                            "" if c.fondo else "   [SIN FONDO]")


def clave_de_etiqueta(txt):
    return txt.split("  -  ")[0].strip()


# ==========================================================================
#  Estado: que corrida esta elegida y con que binning
# ==========================================================================
class Estado:
    """
    Lo que las pestanas comparten. Deliberadamente NO guarda espectros.

    Guardar el espectro seria una cache, y la regla de esta interfaz es que
    cada calculo se rehaga entero. Lo que se comparte es COMO pedirlo: la clave
    de la corrida y el binning. `espectro()` lo vuelve a armar cada vez.
    """

    def __init__(self):
        hay = dat.con_fondo() or dat.disponibles()
        self.clave = hay[0] if hay else None
        cfg = dat.BINNING.get(self.clave, dict(bins=400, rango_max=None))
        self.bins = int(cfg["bins"])
        self.rango_max = cfg["rango_max"]

    def fijar(self, clave):
        """Cambia de corrida y trae el binning que `datos.BINNING` fija para ella."""
        self.clave = clave
        cfg = dat.BINNING.get(clave, dict(bins=400, rango_max=None))
        self.bins = int(cfg["bins"])
        self.rango_max = cfg["rango_max"]

    def espectro(self):
        """(canal, cuentas) recien armado. Sin cache: ver el docstring del modulo."""
        if self.clave is None:
            raise RuntimeError("no hay ninguna corrida accesible en el catalogo")
        return dat.espectro(self.clave, bins=self.bins, rango_max=self.rango_max)

    def descripcion(self):
        if self.clave is None:
            return "sin datos"
        return dat.corrida(self.clave).descripcion


# ==========================================================================
#  Panel de texto con secciones
# ==========================================================================
class Panel(ttk.Frame):
    """Texto con titulos de seccion, scroll y una fuente legible."""

    def __init__(self, master, ancho=52, mono=False, fondo=FONDO_TXT, alto=None):
        super().__init__(master)
        fam = "Consolas" if mono else "Segoe UI"
        kw = dict(width=ancho, wrap="none" if mono else "word", relief="flat",
                  font=(fam, 9), background=fondo, padx=6, pady=6,
                  spacing1=1, spacing3=2)
        if alto:
            kw["height"] = alto
        self.txt = tk.Text(self, **kw)
        sy = ttk.Scrollbar(self, orient="vertical", command=self.txt.yview)
        self.txt.configure(yscrollcommand=sy.set)
        if mono:
            sx = ttk.Scrollbar(self, orient="horizontal", command=self.txt.xview)
            self.txt.configure(xscrollcommand=sx.set)
            sx.pack(side="bottom", fill="x")
        sy.pack(side="right", fill="y")
        self.txt.pack(side="left", fill="both", expand=True)
        self.txt.tag_configure("h", font=(fam, 9, "bold"), foreground="#1f4e79",
                               spacing1=9, spacing3=3)
        self.txt.tag_configure("cita", font=("Georgia", 9, "italic"),
                               foreground="#33383d", lmargin1=14, lmargin2=14)
        self.txt.tag_configure("marca", font=(fam, 9, "bold"),
                               foreground="#a4501a")
        self.txt.tag_configure("mal", foreground="#b32b2b")
        self.txt.tag_configure("mono", font=("Consolas", 9))
        self.txt.configure(state="disabled")

    def poner(self, bloques):
        """bloques: lista de (tag, texto); tag en {None,h,cita,marca,mal,mono}."""
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        for tag, s in bloques:
            self.txt.insert("end", s.rstrip("\n") + "\n",
                            () if tag is None else (tag,))
        self.txt.configure(state="disabled")
        self.txt.yview_moveto(0.0)

    def texto(self, s):
        self.poner([(None, s)])


# ==========================================================================
#  Pestana: controles a la izquierda, figura al centro, el porque a la derecha
# ==========================================================================
class Pestana(ttk.Frame):
    titulo = "pestana"
    auto = True                 # calcular sola al abrirla por primera vez
    necesita_datos = False      # si hace falta cargar el espectro de una corrida
    ancho_res = 46              # ancho del panel de resultado (izquierda)
    ancho_porque = 54           # ancho del panel del porque (derecha)
    centro = "figura"           # "figura" o "ambos" (figura + texto ancho abajo)

    # --- contenido del panel derecho, declarativo -------------------------
    paper = []                  # [(seccion, cita literal), ...]
    codigo = ""                 # que funcion de replica_paper/ hace esto
    abiertas = []               # [(que no dice el paper, que se hizo), ...]
    mirar = ""                  # que hay que mirar en el grafico

    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.calculada = False
        self._corriendo = False
        self._ultimo = None
        self._rango = (0.0, 1.0)

        # ---------------- columna izquierda ----------------
        izq = ttk.Frame(self)
        izq.pack(side="left", fill="y", padx=(6, 3), pady=6)

        self.ctrl = ttk.LabelFrame(izq, text="Controles")
        self.ctrl.pack(side="top", fill="x")
        self.construir_controles(self.ctrl)
        b = ttk.Frame(self.ctrl)
        b.pack(side="top", fill="x", pady=(6, 4))
        self.boton = ttk.Button(b, text="Recalcular", command=self.lanzar)
        self.boton.pack(side="right", padx=6)

        marco_res = ttk.LabelFrame(izq, text="Resultado (recalculado recien)")
        marco_res.pack(side="top", fill="both", expand=True, pady=(6, 0))
        self.p_res = Panel(marco_res, ancho=self.ancho_res, mono=True,
                           fondo="#ffffff")
        self.p_res.pack(fill="both", expand=True)

        # ---------------- columna derecha ----------------
        der = ttk.LabelFrame(self, text="El porque de este paso")
        der.pack(side="right", fill="y", padx=(3, 6), pady=6)
        self.p_por = Panel(der, ancho=self.ancho_porque)
        self.p_por.pack(fill="both", expand=True)
        self.p_por.poner(self.bloques_porque())

        # ---------------- centro ----------------
        centro = ttk.Frame(self)
        centro.pack(side="left", fill="both", expand=True, pady=6)
        if self.centro == "ambos":
            pw = ttk.PanedWindow(centro, orient="vertical")
            pw.pack(fill="both", expand=True)
            arriba = ttk.Frame(pw)
            abajo = ttk.Frame(pw)
            pw.add(arriba, weight=3)
            pw.add(abajo, weight=4)
            self.p_sal = Panel(abajo, ancho=84, mono=True, fondo="#ffffff")
            self.p_sal.pack(fill="both", expand=True)
        else:
            arriba = centro
            self.p_sal = None
        self.fig = Figure(figsize=(8.2, 5.4), dpi=app.dpi_fig)
        self.canvas = FigureCanvasTkAgg(self.fig, master=arriba)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        NavigationToolbar2Tk(self.canvas, arriba).update()
        self.mensaje("todavia no se calculo nada")

    # ------------------------------------------------------------------
    #  panel derecho
    # ------------------------------------------------------------------
    def bloques_porque(self):
        """Arma el panel de la derecha con los atributos declarativos."""
        b = []
        if self.paper:
            b.append(("h", "QUE DICE EL PAPER"))
            for seccion, cita in self.paper:
                b.append((None, seccion))
                b.append(("cita", '"' + cita.strip() + '"'))
        if self.codigo:
            b.append(("h", "QUE CORRE ESTA PESTANA"))
            b.append(("mono", self.codigo.strip()))
        if self.abiertas:
            b.append(("h", "LO QUE EL PAPER NO DICE"))
            b.append(("marca", "[NO esta en el paper]  son los controles de la "
                               "izquierda: moverlos es ver de que depende el "
                               "resultado."))
            for que, como in self.abiertas:
                b.append((None, "• " + que))
                b.append((None, "   " + como))
        if self.mirar:
            b.append(("h", "QUE MIRAR"))
            b.append((None, self.mirar.strip()))
        return b

    # ------------------------------------------------------------------
    #  helpers de controles
    # ------------------------------------------------------------------
    def _fila(self, frame, etiqueta):
        f = ttk.Frame(frame)
        f.pack(side="top", fill="x", padx=6, pady=2)
        ttk.Label(f, text=etiqueta, width=15, anchor="w").pack(side="left")
        return f

    def entrada(self, frame, etiqueta, valor, ancho=10):
        v = tk.StringVar(value=str(valor))
        f = self._fila(frame, etiqueta)
        ttk.Entry(f, textvariable=v, width=ancho).pack(side="left")
        return v

    def combo(self, frame, etiqueta, valores, valor=None, ancho=18, cmd=None):
        v = tk.StringVar(value=valor if valor is not None else valores[0])
        f = self._fila(frame, etiqueta)
        c = ttk.Combobox(f, textvariable=v, values=list(valores), width=ancho,
                         state="readonly")
        c.pack(side="left")
        if cmd is not None:
            c.bind("<<ComboboxSelected>>", cmd)
        return v

    def check(self, frame, etiqueta, valor=True):
        v = tk.BooleanVar(value=bool(valor))
        ttk.Checkbutton(frame, text=etiqueta, variable=v).pack(
            side="top", anchor="w", padx=8, pady=1)
        return v

    def nota(self, frame, texto):
        ttk.Label(frame, text=texto, wraplength=int(230 * self.app.escala),
                  foreground="#5a6570", justify="left").pack(
                      side="top", anchor="w", padx=8, pady=(4, 2))

    def separador(self, frame):
        ttk.Separator(frame, orient="horizontal").pack(
            side="top", fill="x", padx=6, pady=5)

    # ------------------------------------------------------------------
    #  salidas
    # ------------------------------------------------------------------
    def res(self, texto):
        """Escribe en el panel de resultado (izquierda). Hilo principal."""
        self.p_res.texto(texto)

    def salida(self, texto):
        """Escribe en el panel ancho de abajo (solo si centro == 'ambos')."""
        if self.p_sal is not None:
            self.p_sal.texto(texto)

    def mensaje(self, texto, color="0.45"):
        self.fig.clear()
        a = self.fig.add_subplot(111)
        a.text(0.5, 0.5, texto, ha="center", va="center", fontsize=11,
               color=color, wrap=True, transform=a.transAxes)
        a.set_axis_off()
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    #  hilo de calculo
    # ------------------------------------------------------------------
    def leer_controles(self):
        """
        Pasa los controles a valores de Python. Corre en el HILO PRINCIPAL.

        Existe para que `trabajo()` no toque ni una variable de Tk: lee de
        `self.cfg`, que se arma aca justo antes de largar el hilo. Leer un
        StringVar desde otro hilo es exactamente la clase de carrera que
        aparece una vez cada veinte corridas y no se puede reproducir.
        """
        return {}

    def trabajo(self):
        """Se ejecuta EN OTRO HILO. No tocar widgets aca. Devuelve un resultado."""
        return None

    def pintar(self, res):
        """Se ejecuta en el hilo principal con lo que devolvio trabajo()."""

    def progreso(self, texto, frac=None):
        """
        Aviso a la barra de estado DESDE el hilo de calculo.

        Tkinter no es thread-safe: se manda por la misma cola que el resultado
        y el hilo principal lo desvia. `frac=None` deja la barra indeterminada,
        que para lo que se le pide -ver que el programa no se colgo- alcanza.
        """
        if frac is not None:
            a, b = self._rango
            frac = a + max(0.0, min(float(frac), 1.0)) * (b - a)
        self.app.cola.put((self, ("progreso", texto, frac)))

    def fase(self, a, b):
        """Declara que lo que viene ocupa el tramo [a, b] de la barra."""
        self._rango = (float(a), float(b))

    def espectro(self):
        """El espectro de la corrida elegida, recien armado, con aviso."""
        self.progreso("leyendo areas y armando el espectro...")
        e = self.app.estado.espectro()
        self.progreso("espectro listo: %d pulsos" % e["n_pulsos"])
        return e

    def lanzar(self):
        if self._corriendo:
            return
        try:
            self.cfg = self.leer_controles()
        except Exception as e:                               # noqa: BLE001
            self.res("Los controles no se pudieron leer:\n\n%s" % e)
            self.mensaje("Hay un control con un valor invalido\n"
                         "(ver el panel izquierdo)", color="#b32b2b")
            return
        self._corriendo = True
        self.calculada = True
        self._rango = (0.0, 1.0)
        self.boton.configure(state="disabled")
        self.app.estado_txt.set("calculando: %s ..." % self.titulo)
        self.app.progreso_iniciar(self)
        t0 = time.time()

        def worker():
            try:
                r = ("ok", self.trabajo(), time.time() - t0)
            except Exception:
                r = ("error", traceback.format_exc(), time.time() - t0)
            self.app.cola.put((self, r))

        threading.Thread(target=worker, daemon=True).start()

    def recibir(self, res):
        self._corriendo = False
        self.boton.configure(state="normal")
        self.app.progreso_terminar(self)
        estado, carga = res[0], res[1]
        dt = res[2] if len(res) > 2 else 0.0
        if estado == "error":
            self._ultimo = None
            self.mensaje("Error al calcular\n(el detalle esta en el panel "
                         "de la izquierda)", color="#b32b2b")
            self.res(carga)
            self.app.estado_txt.set("error en %s" % self.titulo)
            return
        self._ultimo = carga
        self.fig.clear()
        try:
            self.pintar(carga)
            # tight_layout avisa cuando hay twinx o colorbar; el aviso es
            # cosmetico y ensuciaria la consola en cada recalculo.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                self.fig.tight_layout()
        except Exception:
            self.res(traceback.format_exc())
            self.mensaje("Error al dibujar\n(ver el panel izquierdo)",
                         color="#b32b2b")
        self.canvas.draw()
        self.app.estado_txt.set("%s: listo (%s)" % (self.titulo,
                                                    formato_tiempo(dt)))

    def repintar(self):
        """Redibuja con el ultimo resultado (cambios de nitidez o de tamano)."""
        if self._ultimo is None or self._corriendo:
            return
        self.fig.clear()
        try:
            self.pintar(self._ultimo)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                self.fig.tight_layout()
        except Exception:
            pass
        self.canvas.draw_idle()

    def al_cambiar_corrida(self):
        """Gancho: la corrida cambio. Se llama en el hilo principal."""
        if not self.necesita_datos:
            return
        self.calculada = False
        self._ultimo = None
        self.mensaje("Cambio la corrida.\nApreta «Recalcular».")
        self.res("")


# ==========================================================================
#  0 - Datos: donde ARRANCA el diagrama de la Fig. 2
# ==========================================================================
class TabDatos(Pestana):
    titulo = "0 · Datos"
    necesita_datos = True
    ancho_res = 50

    paper = [("Fig. 2, primer bloque despues de Start:",
              "Energy Spectrum (Channel, Counts)")]
    codigo = ("replica_paper/datos.py -> espectro()\n"
              "  presta de analisis/: area.calibrar_ventana,\n"
              "  area.areas_de_archivo, espectro.espectro_neto")
    abiertas = [
        ("como se llega del archivo del DRS4 al histograma",
         "El paper no lo trata: su diagrama EMPIEZA en (canal, cuentas). Ese "
         "tramo se toma prestado de analisis/ (area con compuerta anclada al "
         "onset, tiempo vivo por el scaler, resta de fondo con binning comun) "
         "en vez de reescribirlo peor, y queda aislado en datos.py."),
        ("el binning",
         "datos.BINNING lo fija por corrida y NO lo deja en automatico: unos "
         "pocos eventos de area enorme (pile-up, cosmicos) estiran el rango "
         "automatico y dejan toda la fisica en los primeros bins. El numero de "
         "bins ademas importa para la Fig. 5: los intervalos mas angostos del "
         "paper cubren el 6% del canal del borde."),
    ]
    mirar = """
El eje horizontal NO es un canal de ADC: es el AREA del pulso en mV*ns. Todo lo
demas del paper funciona igual, porque el metodo solo pide que el eje sea
proporcional a la energia depositada.

Panel de arriba: muestra y fondo YA normalizados por tiempo vivo (los factores
k_m y k_f del resultado). Panel de abajo: el neto, que es lo que se ajusta.

Los negativos del neto se conservan a proposito. Son fluctuaciones estadisticas
legitimas y recortarlas sesga el ajuste hacia arriba.
"""

    def construir_controles(self, f):
        claves = dat.disponibles()
        self.v_cor = self.combo(
            f, "corrida", [etiqueta_corrida(k) for k in claves],
            etiqueta_corrida(self.app.estado.clave) if claves else "",
            ancho=34, cmd=self.al_elegir)
        self.v_bins = self.entrada(f, "bins", self.app.estado.bins)
        self.v_rango = self.entrada(
            f, "rango_max [mV*ns]",
            "auto" if self.app.estado.rango_max is None
            else "%g" % self.app.estado.rango_max)
        self.nota(f, "Sin fondo propio no hay espectro NETO: el catalogo lo "
                     "marca [SIN FONDO] y espectro_neto devuelve el crudo. "
                     "El Co-60 tiene DOS bordes solapados, asi que sus "
                     "ajustes tardan ~10 veces mas que los del Cs-137.")

    def al_elegir(self, _e=None):
        """Al cambiar de corrida, traer el binning que datos.BINNING le fija."""
        self.app.estado.fijar(clave_de_etiqueta(self.v_cor.get()))
        self.v_bins.set(str(self.app.estado.bins))
        self.v_rango.set("auto" if self.app.estado.rango_max is None
                         else "%g" % self.app.estado.rango_max)

    def leer_controles(self):
        est = self.app.estado
        anterior = est.clave
        est.clave = clave_de_etiqueta(self.v_cor.get())
        est.bins = int(self.v_bins.get())
        r = self.v_rango.get().strip().lower()
        est.rango_max = None if r in ("auto", "none", "") else float(r)
        if est.clave != anterior:
            self.app.propagar_corrida(excepto=self)
        self.app.actualizar_titulo()
        return {}

    def trabajo(self):
        return self.espectro()

    def pintar(self, e):
        a1 = self.fig.add_subplot(211)
        a2 = self.fig.add_subplot(212, sharex=a1)
        ch = e["canales"]
        a1.step(ch, e["muestra"], where="mid", lw=0.9, color="tab:blue",
                label="muestra / tiempo vivo")
        if e["hay_fondo"]:
            a1.step(ch, e["fondo"], where="mid", lw=0.9, color="tab:orange",
                    label="fondo / tiempo vivo")
        a1.set_yscale("log")
        a1.set_ylabel("cuentas por segundo vivo")
        a1.legend(fontsize=8)
        a1.grid(alpha=0.3)
        a1.set_title("%s  (%s)" % (e["clave"], e["punto"]), fontsize=10)

        a2.fill_between(ch, e["cuentas"] - e["error"], e["cuentas"] + e["error"],
                        step="mid", color="0.8", label="+-1 sigma")
        a2.step(ch, e["cuentas"], where="mid", lw=0.9, color="k", label="neto")
        a2.axhline(0, color="0.5", lw=0.7)
        a2.set_xlabel("canal = area del pulso [mV*ns]")
        a2.set_ylabel("cuentas netas")
        a2.legend(fontsize=8)
        a2.grid(alpha=0.3)

        c = e["corrida"]
        L = ["corrida        %s" % e["clave"],
             "descripcion    %s" % c.descripcion,
             "punto de oper. %s" % e["punto"],
             "fuente         %s" % e["fuente"],
             "",
             "lineas gamma y su borde (Ec. 1 corregida):"]
        for eg in e["lineas"]:
            L.append("   %9.3f keV  ->  E_max = %7.1f keV"
                     % (eg, float(ec.e_max(eg))))
        L += ["",
              "pulsos          %d" % e["n_pulsos"],
              "bins            %d  (%.2f mV*ns por bin)"
              % (e["bins"], e["canales"][1] - e["canales"][0]),
              "rango_max       %.0f mV*ns" % e["rango_max"],
              "fondo propio    %s" % ("si" if e["hay_fondo"] else "NO"),
              "k_m (muestra)   %.4f" % e["k_m"],
              "k_f (fondo)     %.4f" % e["k_f"],
              "",
              "El neto es  muestra*k_m - fondo*k_f, con k = 1/tiempo vivo.",
              "El fondo entra con k_f y no con el tiempo de pared: en el",
              "Co-60 la placa estuvo viva el 40% del tiempo.",
              "",
              "cuentas netas totales   %.0f" % np.sum(e["cuentas"]),
              "bins con neto negativo  %d de %d"
              % (int(np.sum(e["cuentas"] < 0)), e["cuentas"].size)]
        self.res("\n".join(L))


# ==========================================================================
#  1 - Ec. 1: el borde Compton
# ==========================================================================
class TabEc1(Pestana):
    titulo = "1 · Ec. 1 · Borde"
    ancho_res = 52

    paper = [("Ec. 1, tal como se imprime:",
              "E_max = E_gamma / (1 + 2 E_gamma / m_e c^2)")]
    codigo = ("replica_paper/ecuaciones.py -> e_max_impresa(), e_max()\n"
              "replica_paper/auditoria.py  -> ec1()")
    abiertas = [
        ("cuanto vale m_e c^2",
         "El paper la llama 'electron rest mass energy' y no da el valor. 511 "
         "keV es el que reproduce la Fig. 1 y los bordes rotulados en la "
         "Fig. 3."),
    ]
    mirar = """
LA ECUACION IMPRESA ESTA MAL. Lo que escribe es la energia del FOTON
retrodispersado, no la del electron de retroceso. La correcta lleva 2*alpha en
el numerador, con alpha = E_gamma / m_e c^2:

    E_max = 2 alpha E_gamma / (1 + 2 alpha)

No hace falta creerle a nadie: el propio paper rotula los bordes de su Fig. 3
en 477, 639, 963, 1061 y 1117 keV, y esos numeros salen de la CORREGIDA. La
impresa daria 184, 196, 210, 213 y 214 keV.

O sea que el paper usa la formula corregida y imprime la otra. La tabla del
resultado tiene las dos columnas al lado del rotulo publicado.

En el grafico: la curva roja es la impresa, la negra la corregida, y los puntos
son los cinco bordes que el paper rotula en su Fig. 3.
"""

    def construir_controles(self, f):
        self.v_eg = self.entrada(f, "E_gamma [keV]", "661.657")
        self.v_me = self.entrada(f, "m_e c^2 [keV]", "511.0")
        self.v_max = self.entrada(f, "eje hasta [keV]", "2000")
        self.nota(f, "Los puntos del grafico son los bordes que el paper "
                     "ROTULA en su Fig. 3 (auditoria.FIG3), no un ajuste.")

    def leer_controles(self):
        return dict(eg=float(self.v_eg.get()), me=float(self.v_me.get()),
                    xmax=float(self.v_max.get()))

    def trabajo(self):
        c = self.cfg
        x = np.linspace(1.0, c["xmax"], 1200)
        return dict(x=x, imp=ec.e_max_impresa(x, c["me"]),
                    cor=ec.e_max(x, c["me"]), cfg=c,
                    tabla=texto_de(auditoria.ec1))

    def pintar(self, r):
        c = r["cfg"]
        a = self.fig.add_subplot(111)
        a.plot(r["x"], r["cor"], "k-", lw=1.6, label="Ec. 1 corregida  "
               "$2\\alpha E_\\gamma/(1+2\\alpha)$")
        a.plot(r["x"], r["imp"], "-", color="#b32b2b", lw=1.4,
               label="Ec. 1 impresa  $E_\\gamma/(1+2\\alpha)$")
        a.plot(r["x"], r["x"], ":", color="0.6", lw=1.0,
               label="$E_{max}=E_\\gamma$ (imposible)")
        rot = {661.657: 477, 834.848: 639, 1173.228: 963,
               1274.537: 1061, 1332.492: 1117}
        a.plot(list(rot.keys()), list(rot.values()), "o", color="tab:blue",
               ms=7, zorder=5, label="rotulos de la Fig. 3 del paper")
        eg = c["eg"]
        a.axvline(eg, color="tab:green", lw=1.0, ls="--")
        a.plot([eg], [float(ec.e_max(eg, c["me"]))], "*", color="tab:green",
               ms=14, zorder=6, label="E_gamma elegida")
        a.set_xlabel("energia de la linea gamma  $E_\\gamma$ [keV]")
        a.set_ylabel("energia del borde Compton  $E_{max}$ [keV]")
        a.grid(alpha=0.3)
        a.legend(fontsize=8, loc="upper left")

        alpha = eg / c["me"]
        L = ["E_gamma        %.3f keV" % eg,
             "alpha = Eg/mec2  %.5f" % alpha,
             "",
             "Ec. 1 impresa    %8.2f keV   <- foton retrodispersado",
             "Ec. 1 corregida  %8.2f keV   <- electron (la buena)",
             ]
        L[3] = L[3] % float(ec.e_max_impresa(eg, c["me"]))
        L[4] = L[4] % float(ec.e_max(eg, c["me"]))
        L += ["cociente         %8.2f" % (float(ec.e_max(eg, c["me"]))
                                          / max(float(ec.e_max_impresa(eg, c["me"])), 1e-9)),
              "",
              "La suma de las dos da E_gamma, que es lo que",
              "tiene que pasar: el foton retrodispersado se",
              "lleva lo que el electron no.",
              "   %.2f + %.2f = %.2f"
              % (float(ec.e_max(eg, c["me"])),
                 float(ec.e_max_impresa(eg, c["me"])), eg),
              "", "-" * 46, "",
              r["tabla"]]
        self.res("\n".join(L))


# ==========================================================================
#  2 - Ec. 2 y 10: Klein-Nishina
# ==========================================================================
class TabEc10(Pestana):
    titulo = "2 · Ec. 2 y 10 · Klein-Nishina"
    ancho_res = 52

    # Valores del eje de la Fig. 1 del paper, leidos de la propia figura.
    # Son los que ecuaciones.py cita en su docstring y README.md tabula.
    FIG1 = {500.0: (1.01, 1.68), 662.0: (0.58, 1.12), 1000.0: (0.255, 0.645),
            2000.0: (0.064, 0.283), 5000.0: (0.010, 0.11)}

    paper = [("Ec. 2: Klein-Nishina diferencial en angulo solido.", "d sigma / d Omega"),
             ("Ec. 10: la misma, pasada a energia del electron por el jacobiano.",
              "d sigma/dE = A pi r0^2/(alpha^2 m0c^2) * "
              "[ 2 + (E/(hv-E))^2 ( 1/alpha^2 - 2(hv-E)/(alpha E) + (hv-E)/hv ) ]"),
             ("Sobre r0 el paper dice:",
              "approximately equal to 2.818e-15 m")]
    codigo = ("replica_paper/ecuaciones.py -> dsigma_domega(), dsigma_de()")
    abiertas = [
        ("en que unidades esta el eje vertical de la Fig. 1",
         "El paper no lo rotula. Con r0 en cm y m0c^2 en keV la Ec. 10 da "
         "cm^2/keV y reproduce el eje: 1.02e-27 en E->0 para E0 = 500 keV. "
         "Es la comprobacion de que la Ec. 10 esta bien impresa."),
    ]
    mirar = """
ESTA ECUACION SI ESTA BIEN IMPRESA, y la tabla del resultado lo demuestra:
compara lo que la Ec. 10 calcula contra los valores del eje de la Fig. 1 del
paper, en E->0 y en el pico.

Tres cosas de la forma, que son las que despues hace o no hace cada metodo:

  · el PLATEAU: casi plano a baja energia;
  · la ACUMULACION antes del borde: la curva SUBE hacia E_max. Por esto el
    50% de altura y el minimo de la derivada se sesgan (pestana 7);
  · el CORTE en E_max: cae a cero de golpe. Eso es el borde Compton, y la
    convolucion de la Ec. 8 es lo unico que lo suaviza.

El chequeo Ec. 2 -> Ec. 10 integra las dos por su propia variable (angulo
solido y energia del electron) y compara. Tienen que dar la misma seccion
eficaz total, porque una es la otra cambiada de variable. Si el cociente da 1,
el jacobiano del paper esta bien.

A no es un parametro fisico: absorbe el prefactor. Para el NCC y para el ajuste
solo importa la FORMA.
"""

    def construir_controles(self, f):
        self.v_e0 = self.entrada(f, "E0 [keV]", "500, 662, 1000, 2000, 5000",
                                 ancho=22)
        self.v_n = self.entrada(f, "puntos", "40001")
        self.v_xmax = self.entrada(f, "eje hasta [keV]", "5000")
        self.v_chk = self.check(f, "verificar Ec. 2 -> Ec. 10", True)
        self.nota(f, "El chequeo integra las dos ecuaciones y compara la "
                     "seccion eficaz total. Tarda ~1 s.")

    def leer_controles(self):
        e0 = [float(s) for s in self.v_e0.get().replace(";", ",").split(",")
              if s.strip()]
        return dict(e0=e0, n=int(self.v_n.get()), xmax=float(self.v_xmax.get()),
                    chk=bool(self.v_chk.get()))

    def trabajo(self):
        c = self.cfg
        E = np.linspace(0.0, c["xmax"], c["n"])
        curvas, tabla = [], []
        for e0 in c["e0"]:
            y = ec.dsigma_de(E, e0)
            y0 = float(ec.dsigma_de(np.array([1e-9]), e0)[0])
            curvas.append((e0, y))
            tabla.append((e0, y0, float(y.max()), float(ec.e_max(e0)),
                          self.FIG1.get(e0)))
        chk = []
        if c["chk"]:
            self.progreso("integrando la Ec. 2 y la Ec. 10...")
            for e0 in c["e0"]:
                th = np.linspace(0.0, np.pi, 200001)
                s_om = 2 * np.pi * np.trapezoid(
                    ec.dsigma_domega(th, e0) * np.sin(th), th)
                Em = float(ec.e_max(e0))
                Ee = np.linspace(1e-6, Em, 400001)
                s_e = float(np.trapezoid(ec.dsigma_de(Ee, e0), Ee))
                chk.append((e0, float(s_om), s_e))
        return dict(E=E, curvas=curvas, tabla=tabla, chk=chk, cfg=c)

    def pintar(self, r):
        a = self.fig.add_subplot(111)
        for e0, y in r["curvas"]:
            a.plot(r["E"], y, lw=1.4, label="$E_0$ = %g keV" % e0)
            a.axvline(float(ec.e_max(e0)), color="0.75", lw=0.7, ls=":")
        a.set_xlim(0, r["cfg"]["xmax"])
        a.set_xlabel("energia del electron de retroceso  E [keV]")
        a.set_ylabel(r"d$\sigma$/dE  [cm$^2$/keV]")
        a.set_title("Fig. 1 del paper, recalculada con la Ec. 10", fontsize=10)
        a.grid(alpha=0.35)
        a.legend(fontsize=8)

        L = ["Ec. 10 contra el eje de la Fig. 1 del paper",
             "(los valores del paper estan en unidades de 1e-27)", "",
             " E0 [keV]   E->0      pico      E_max    paper (E->0/pico)"]
        for e0, y0, pk, em, pap in r["tabla"]:
            s = "  %7.0f  %8.4f  %8.4f  %8.1f" % (e0, y0 * 1e27, pk * 1e27, em)
            s += ("    %.3f / %.3f" % pap) if pap else "    -"
            L.append(s)
        if r["chk"]:
            L += ["", "-" * 46, "",
                  "Ec. 2 -> Ec. 10 (el jacobiano del paper)",
                  "integral de cada una, en cm^2:", "",
                  " E0 [keV]   por angulo    por energia   cociente"]
            for e0, so, se in r["chk"]:
                L.append("  %7.0f  %.6e  %.6e  %.6f"
                         % (e0, so, se, se / so))
            L += ["", "Cociente 1 = la Ec. 10 es la Ec. 2 cambiada",
                  "de variable. El paper cita el jacobiano y no lo",
                  "muestra; esto lo verifica con numeros."]
        self.res("\n".join(L))


# ==========================================================================
#  3 - Ec. 8 y 9: lineshape -> canal + convolucion
# ==========================================================================
class TabForward(Pestana):
    titulo = "3 · Ec. 8-9 · Forward model"
    ancho_res = 50

    paper = [("Ec. 8: la lineshape llevada a canal y convolucionada.",
              "E_eff = E_e * B / E_max"),
             ("Ec. 9: el nucleo del ensanchamiento.",
              "G(E) = 1/(sigma sqrt(2 pi)) exp(-E^2 / (2 sigma^2))"),
             ("Ec. 3 y Tabla I: de donde sale que sea gaussiano.",
              "sigma_total = sqrt(sigma_LY^2 + sigma_gain^2 + sigma_elec^2)")]
    codigo = ("replica_paper/ecuaciones.py -> modelo(), kernel_gaussiano()")
    abiertas = [
        ("hasta donde se evalua la convolucion",
         "La Ec. 8 es una integral sobre TODO E_e, pero un ajuste se hace "
         "dentro de una ventana. Convolucionar solo adentro es un error "
         "numerico: np.convolve rellena con ceros y pierde los eventos de "
         "abajo que el ensanchamiento deberia meter hacia adentro. `pad` "
         "extiende la grilla al canal 0 y 4 sigma por arriba, convoluciona y "
         "recorta. Destildalo para ver el efecto."),
        ("como se normaliza el nucleo",
         "A suma 1, no a integral 1, porque se usa sobre un histograma: asi la "
         "convolucion conserva las CUENTAS y A queda desacoplado de sigma. "
         "Subir sigma ensancha y baja el pico sin cambiar el total."),
    ]
    mirar = """
EL MAPA CANAL-ENERGIA DE LA EC. 8 PASA POR EL ORIGEN. E_eff = E_e * B / E_max
no tiene ordenada: area cero corresponde a energia cero. Esto vuelve mas tarde
en dos lugares y conviene tenerlo presente:

  · en la pestana 4, porque fija donde cae el borde dentro de la ventana;
  · en la pestana 6, porque D = B / E_max queda determinado por B, y la Ec. 11
    ademas deja D libre. El modelo del paper se sobre-parametriza respecto de
    su propia Ec. 8.

En el grafico: la linea gris es la lineshape SIN ensanchar (el escalon de la
Ec. 10 llevado a canal) y la negra es despues de la Ec. 9. La linea punteada
marca B, o sea donde esta el borde de verdad. Fijate que el maximo de la curva
convolucionada NO esta en B, y que el 50% de la altura tampoco: eso es todo el
problema que el paper resuelve.

Con `pad` destildado, la curva roja muestra lo que pasa si la convolucion se
evalua solo dentro de la ventana: el modelo pierde cuentas cerca del borde
izquierdo, y un ajuste lo compensa achicando sigma.
"""

    def construir_controles(self, f):
        self.v_eg = self.entrada(f, "E_gamma [keV]", "661.657")
        self.v_B = self.entrada(f, "B [canal]", "900")
        self.v_sr = self.entrada(f, "sigma/B", "0.12")
        self.v_lo = self.entrada(f, "ventana desde", "400")
        self.v_hi = self.entrada(f, "ventana hasta", "1400")
        self.v_bins = self.entrada(f, "bins", "500")
        self.v_pad = self.check(f, "pad (extender la grilla)", True)

    def leer_controles(self):
        return dict(eg=float(self.v_eg.get()), B=float(self.v_B.get()),
                    sr=float(self.v_sr.get()), lo=float(self.v_lo.get()),
                    hi=float(self.v_hi.get()), bins=int(self.v_bins.get()),
                    pad=bool(self.v_pad.get()))

    def trabajo(self):
        c = self.cfg
        sig = c["sr"] * c["B"]
        ch = np.linspace(c["lo"], c["hi"], c["bins"])
        paso = ch[1] - ch[0]
        emax = float(ec.e_max(c["eg"]))
        crudo = ec.dsigma_de(ch * emax / c["B"], c["eg"], A=1.0)
        con_pad = ec.modelo(ch, c["eg"], 1.0, c["B"], sig, pad=True)
        sin_pad = ec.modelo(ch, c["eg"], 1.0, c["B"], sig, pad=False)
        k = ec.kernel_gaussiano(sig, paso, largo_senal=ch.size)
        # conservacion: la misma A con tres sigma distintos
        cons = []
        for f_ in (0.5, 1.0, 2.0):
            y = ec.modelo(ch, c["eg"], 1.0, c["B"], f_ * sig, pad=True)
            cons.append((f_ * sig, float(np.sum(y) * paso)))
        return dict(ch=ch, crudo=crudo, con=con_pad, sin=sin_pad, k=k,
                    paso=paso, sig=sig, emax=emax, cons=cons, cfg=c)

    def pintar(self, r):
        c = r["cfg"]
        a1 = self.fig.add_subplot(211)
        a2 = self.fig.add_subplot(212)
        esc = r["crudo"].max() if r["crudo"].max() > 0 else 1.0
        a1.plot(r["ch"], r["crudo"] / esc, color="0.65", lw=1.2,
                label="Ec. 10 en canal (sin ensanchar)")
        a1.plot(r["ch"], r["con"] / esc, "k-", lw=1.8,
                label="Ec. 8 con pad (lo correcto)")
        if not c["pad"]:
            a1.plot(r["ch"], r["sin"] / esc, "-", color="#b32b2b", lw=1.4,
                    label="Ec. 8 SIN pad (artefacto numerico)")
        a1.axvline(c["B"], color="tab:blue", ls="--", lw=1.2,
                   label="B = canal del borde")
        a1.set_ylabel("normalizado al plateau")
        a1.set_xlabel("canal")
        a1.grid(alpha=0.3)
        a1.legend(fontsize=8)
        a1.set_title("Ec. 8:  sigma/B = %.3f  ->  sigma = %.1f canales"
                     % (c["sr"], r["sig"]), fontsize=10)

        x = (np.arange(r["k"].size) - r["k"].size // 2) * r["paso"]
        a2.plot(x, r["k"], "o-", ms=2.5, color="tab:purple", lw=1.1)
        a2.set_xlabel("canal relativo")
        a2.set_ylabel("Ec. 9  G(E)")
        a2.grid(alpha=0.3)
        a2.set_title("el nucleo, muestreado y normalizado a suma 1 "
                     "(suma = %.6f, %d puntos)" % (r["k"].sum(), r["k"].size),
                     fontsize=9)

        i_pico = int(np.argmax(r["con"]))
        mitad = 0.5 * r["con"].max()
        der = r["ch"][i_pico:]
        yd = r["con"][i_pico:]
        cruce = float(np.interp(-mitad, -yd, der)) if yd.size > 2 else float("nan")
        L = ["B (verdadero)          %8.1f canales" % c["B"],
             "sigma                  %8.1f canales" % r["sig"],
             "E_max de la linea      %8.1f keV" % r["emax"],
             "D = B / E_max          %8.4f canal/keV" % (c["B"] / r["emax"]),
             "",
             "maximo del modelo en   %8.1f  (%+.1f%% de B)"
             % (r["ch"][i_pico], 100 * (r["ch"][i_pico] / c["B"] - 1)),
             "50%% de altura en       %8.1f  (%+.1f%% de B)"
             % (cruce, 100 * (cruce / c["B"] - 1)),
             "",
             "Ninguno de los dos cae en B. Ese corrimiento es el",
             "sesgo que la pestana 7 mide sobre los dos metodos",
             "convencionales.",
             "",
             "-" * 44, "",
             "conservacion de cuentas (A fija, sigma variable):"]
        for s, tot in r["cons"]:
            L.append("   sigma = %7.1f   integral = %.6e" % (s, tot))
        L += ["", "El nucleo suma 1, asi que cambiar sigma no cambia",
              "el area: A y sigma quedan desacoplados.",
              "", "-" * 44, "",
              "efecto del pad, dentro de la ventana:"]
        d = r["con"] - r["sin"]
        rel = 100 * float(np.sum(np.abs(d))) / max(float(np.sum(r["con"])), 1e-30)
        L += ["   cuentas que el pad recupera: %.2f%% del total" % rel,
              "   maxima diferencia local:     %.2f%% del pico"
              % (100 * float(np.max(np.abs(d))) / max(float(r["con"].max()), 1e-30)),
              "",
              "El pad no agrega nada al modelo: corrige la",
              "evaluacion. np.convolve rellena con ceros fuera de",
              "la ventana, y esas cuentas si existen."]
        self.res("\n".join(L))


def tres_metodos(canales, cuentas, lineas, intervalo=None, sigma_rel=None,
                 ligar_origen=True, **extra):
    """
    Los tres metodos sobre el mismo espectro y la MISMA ventana.

    Es `correr.tres_metodos` con las ambiguedades del paper abiertas: esta
    interfaz necesita poder cambiar el objetivo del SA (Ec. 14 contra el error
    NCC), la version de la Ec. 12, el pulido local y las ligaduras, que es
    justamente de lo que se trata cada pestana. Dos cosas NO cambian, para que
    los numeros de la interfaz sean los mismos que los del informe:

      * los dos metodos convencionales corren en la ventana que eligio el NCC,
        no en una propia (es lo unico que hace justa la comparacion);
      * con dos lineas gamma, `ligar_origen` arranca en True, o sea D =
        B1/E_max,1 como impone la Ec. 8. Es el default de `correr.py`. La
        pestana 6 es la que lo deja destildar, porque ahi es el punto.
    """
    kw = dict(extra)
    if len(lineas) == 2:
        kw.setdefault("ligar_origen", ligar_origen)
    r = metodo.ajustar(canales, cuentas, lineas, intervalo=intervalo,
                       sigma_rel=sigma_rel, **kw)
    iv = r["intervalo"]
    try:
        g = conv.metodo_gauss(canales, cuentas, iv)["borde"]
    except Exception:                                        # noqa: BLE001
        g = float("nan")
    try:
        d = conv.metodo_diferencia(canales, cuentas, iv)["borde"]
    except Exception:                                        # noqa: BLE001
        d = float("nan")
    return {"NCC-SA": correr._borde_de(r), "Gauss FWHM": g, "Diferencia": d}, r


def _tabla_barrido(barrido, elegido):
    L = ["  sigma_rel de S2    rho        B (NCC)"]
    for sr, rho, B in barrido:
        L.append("      %.4f       %.5f   %8.1f%s"
                 % (sr, rho, B, "   <- elegido" if sr == elegido else ""))
    return L


# ==========================================================================
#  4 - Ec. 15: NCC, seleccion automatica del intervalo
# ==========================================================================
class TabNCC(Pestana):
    titulo = "4 · Ec. 15 · NCC"
    necesita_datos = True
    auto = True                 # ~3 s: se puede calcular sola al abrirla
    ancho_res = 50

    paper = [("Seccion II-B, literal:",
              "we first generate 100 reference points (S2) by uniformly "
              "sampling the predicted lineshape from the Compton plateau to "
              "the edge minimum. For experimental data, each non-zero point "
              "serves as a potential edge start position, with subsequent "
              "points forming test intervals that are interpolated to 100 "
              "points (S1). The interval showing maximum Pearson correlation "
              "(rho) between S1 and S2 is identified as the optimal Compton "
              "edge region."),
             ("Ec. 15:", "rho de Pearson entre S1 y S2")]
    codigo = ("replica_paper/ncc.py -> plantilla(), buscar_intervalo(),\n"
              "                        barrer_sigma_rel()\n"
              "replica_paper/ecuaciones.py -> pearson()  (Ec. 15)")
    abiertas = [
        ("con que sigma se construye S2",
         "Una Klein-Nishina SIN ensanchar no tiene descenso: se corta de "
         "golpe, asi que 'to the edge minimum' no existiria. S2 tiene que "
         "estar convolucionada con ALGO, y sigma es justo lo que el SA ajusta "
         "despues. Con 'barrer' gana el de mayor rho, que es el unico criterio "
         "que el paper da."),
        ("donde empieza el plateau y donde termina el minimo del borde",
         "ini_rel y fin_rel, en unidades de E_max. NO es cosmetico: FIJA EL "
         "MAPA. Si S2 abarca [ini*E_max, fin*E_max] y gana el intervalo "
         "[lo, hi], el borde cae en B = lo + (1-ini)/(fin-ini)*(hi-lo). Mover "
         "el span de la plantilla mueve B en proporcion."),
        ("el largo minimo de un intervalo candidato",
         "min_bins. Pearson ignora escala y offset, asi que sobre un punado de "
         "bins de cola vacia encuentra correlacion alta con cualquier forma."),
        ("cuantos pares (inicio, fin) se prueban",
         "El barrido es O(N^2) en bins. n_max fija cuantos inicios y cuantos "
         "finales se muestrean; rho varia suave con los extremos, asi que "
         "muestrear no cambia el ganador."),
    ]
    mirar = """
QUE HACE PEARSON ACA: compara FORMA. Es invariante a offset y a escala, asi que
no le importa cuantas cuentas hay ni donde esta el cero; le importa el perfil.
Por eso puede elegir la ventana ANTES de ajustar nada.

Panel de arriba: el espectro, la ventana ganadora sombreada y B, el canal del
borde que se deduce del span de la plantilla.

Panel de abajo a la izquierda: S1 y S2 estandarizados, o sea lo que Pearson ve
de verdad. Si las dos curvas se pisan, rho es alto.

Panel de abajo a la derecha: con 'barrer', rho y B contra el sigma de la
plantilla; sin 'barrer', el mapa de rho sobre todos los pares (inicio, fin)
probados.

UNA TRAMPA QUE CONVIENE VER EN VIVO: fin_rel depende de sigma_rel (fin = 1 +
4.5 sigma_rel si se deja en auto), asi que al barrer sigma_rel se estan
comparando rho medidos sobre CONJUNTOS DE VENTANAS DISTINTOS. rho tiende a
crecer con sigma_rel simplemente porque la plantilla se vuelve mas suave.
Quedarse con el de mayor rho es lo unico que el paper ofrece, pero es un
criterio confundido: por eso se muestra el barrido entero y no solo el ganador.
Fijando fin_rel a mano se puede ver cuanto de la subida era eso.
"""

    def construir_controles(self, f):
        self.v_barrer = self.check(f, "barrer sigma_rel (gana el de mayor rho)",
                                   True)
        self.v_sr = self.entrada(f, "sigma_rel fijo", "0.16")
        self.v_ini = self.entrada(f, "ini_rel", "%.2f" % mod_ncc.INI_REL)
        self.v_fin = self.entrada(f, "fin_rel", "auto")
        self.v_min = self.entrada(f, "min_bins", str(mod_ncc.MIN_BINS))
        self.v_n = self.entrada(f, "n_puntos (S1,S2)", str(mod_ncc.N_PUNTOS))
        self.v_nmax = self.entrada(f, "n_max del barrido", "220")
        self.nota(f, "fin_rel 'auto' = 1 + 4.5*sigma_rel. Con dos lineas la "
                     "plantilla se arma con la de MENOR energia: en la Fig. 2 "
                     "el NCC va ANTES de fijar Theta, o sea antes de que "
                     "exista D.")

    def leer_controles(self):
        fin = self.v_fin.get().strip().lower()
        return dict(barrer=bool(self.v_barrer.get()),
                    sr=float(self.v_sr.get()), ini=float(self.v_ini.get()),
                    fin=None if fin in ("auto", "none", "") else float(fin),
                    minb=int(self.v_min.get()), n=int(self.v_n.get()),
                    nmax=int(self.v_nmax.get()))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        eg = e["lineas"][0]
        kw = dict(ini_rel=c["ini"], fin_rel=c["fin"], min_bins=c["minb"],
                  n=c["n"], n_max=c["nmax"])
        if c["barrer"]:
            self.progreso("barriendo el sigma de la plantilla S2...")
            cand = mod_ncc.barrer_sigma_rel(e["canales"], e["cuentas"], eg,
                                            metodo.SIGMAS_S2, **kw)
            r = max(cand, key=lambda x: x["rho"])
            barrido = [(x["sigma_rel"], x["rho"], x["B"]) for x in cand]
        else:
            self.progreso("barriendo intervalos (Seccion II-B)...")
            r = mod_ncc.buscar_intervalo(e["canales"], e["cuentas"], eg,
                                         sigma_rel=c["sr"], **kw)
            barrido = [(r["sigma_rel"], r["rho"], r["B"])]
        return dict(e=e, r=r, barrido=barrido, eg=eg, cfg=c)

    def pintar(self, d):
        e, r, c = d["e"], d["r"], d["cfg"]
        gs = self.fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0])
        a0 = self.fig.add_subplot(gs[0, :])
        a1 = self.fig.add_subplot(gs[1, 0])
        a2 = self.fig.add_subplot(gs[1, 1])

        a0.step(e["canales"], e["cuentas"], where="mid", lw=0.9, color="0.45",
                label="neto")
        a0.axvspan(r["lo"], r["hi"], color="tab:orange", alpha=0.16,
                   label="intervalo ganador")
        a0.axvline(r["B"], color="tab:red", lw=1.5,
                   label="B = %.1f (deducido del span de S2)" % r["B"])
        # la plantilla, puesta en canal con el mismo mapa
        xs = r["lo"] + (r["x_rel"] - r["ini_rel"]) / (r["fin_rel"] - r["ini_rel"]) \
            * (r["hi"] - r["lo"])
        m = (e["canales"] >= r["lo"]) & (e["canales"] <= r["hi"])
        esc = np.max(e["cuentas"][m]) / max(np.max(r["s2"]), 1e-30) if m.any() else 1.0
        a0.plot(xs, r["s2"] * esc, color="tab:blue", lw=1.6,
                label="S2 (Ec. 8-10) escalada")
        a0.set_xlabel("canal [mV*ns]")
        a0.set_ylabel("cuentas netas")
        a0.grid(alpha=0.3)
        a0.legend(fontsize=8)
        a0.set_title("%s   rho = %.5f   sigma_rel de S2 = %.3f"
                     % (e["clave"], r["rho"], r["sigma_rel"]), fontsize=10)

        tramo = e["cuentas"][r["i"]:r["j"] + 1]
        u = np.linspace(0, 1, r["n_puntos"])
        s1 = np.interp(u, np.linspace(0, 1, tramo.size), tramo)
        z = lambda v: (v - v.mean()) / (v.std() if v.std() > 0 else 1.0)
        a1.plot(u, z(s1), lw=1.2, color="0.3", label="S1 (datos)")
        a1.plot(u, z(r["s2"]), lw=1.4, color="tab:blue", label="S2 (Ec. 8-10)")
        a1.set_xlabel("punto reinterpolado")
        a1.set_ylabel("estandarizado")
        a1.set_title("Ec. 15: lo que Pearson compara", fontsize=9)
        a1.grid(alpha=0.3)
        a1.legend(fontsize=8)

        if c["barrer"] and len(d["barrido"]) > 1:
            sr = [x[0] for x in d["barrido"]]
            rho = [x[1] for x in d["barrido"]]
            Bs = [x[2] for x in d["barrido"]]
            a2.plot(sr, rho, "o-", color="tab:purple", lw=1.2, ms=4)
            a2.set_xlabel("sigma_rel de la plantilla S2")
            a2.set_ylabel("rho", color="tab:purple")
            a2.tick_params(axis="y", labelcolor="tab:purple")
            a2b = a2.twinx()
            a2b.plot(sr, Bs, "s--", color="tab:red", lw=1.0, ms=3.5)
            a2b.set_ylabel("B [mV*ns]", color="tab:red")
            a2b.tick_params(axis="y", labelcolor="tab:red")
            a2.axvline(r["sigma_rel"], color="0.6", ls=":", lw=1.0)
            a2.set_title("rho y B contra el sigma de S2", fontsize=9)
            a2.grid(alpha=0.3)
        else:
            mp = r["mapa"]
            sc = a2.scatter(mp[:, 0], mp[:, 1] - mp[:, 0], c=mp[:, 2], s=5,
                            cmap="viridis")
            a2.plot([r["lo"]], [r["hi"] - r["lo"]], "*", color="red", ms=14)
            self.fig.colorbar(sc, ax=a2, label="rho")
            a2.set_xlabel("inicio del intervalo")
            a2.set_ylabel("ancho")
            a2.set_title("%d pares (inicio, fin) probados" % mp.shape[0],
                         fontsize=9)

        emax = float(ec.e_max(d["eg"]))
        L = ["fuente        %s" % e["fuente"],
             "plantilla     E_gamma = %.3f keV  ->  E_max = %.1f keV"
             % (d["eg"], emax),
             "",
             "INTERVALO GANADOR",
             "  lo, hi        %.1f, %.1f mV*ns" % (r["lo"], r["hi"]),
             "  ancho         %d bins" % (r["j"] - r["i"] + 1),
             "  rho           %.6f" % r["rho"],
             "  sigma_rel S2  %.4f" % r["sigma_rel"],
             "  ini_rel       %.3f" % r["ini_rel"],
             "  fin_rel       %.3f" % r["fin_rel"],
             "  pasos         inicio %d, fin %d"
             % (r["paso_inicio"], r["paso_fin"]),
             "",
             "DE AHI SALE EL BORDE",
             "  frac del borde dentro de la ventana:",
             "     (1 - ini)/(fin - ini) = %.4f"
             % ((1 - r["ini_rel"]) / (r["fin_rel"] - r["ini_rel"])),
             "  B = lo + frac*(hi - lo) = %.1f mV*ns" % r["B"],
             "  sigma semilla = sigma_rel*B = %.1f" % r["sigma"],
             "  D = B / E_max = %.4f mV*ns/keV" % (r["B"] / emax),
             "",
             "Ese D es solo la SEMILLA. El numero del ajuste sale",
             "de la pestana 6.",
             ""]
        if len(d["barrido"]) > 1:
            L += ["-" * 44, "",
                  "barrido del sigma de la plantilla:"]
            L += _tabla_barrido(d["barrido"], r["sigma_rel"])
            L += ["",
                  "OJO: fin_rel depende de sigma_rel, asi que estos",
                  "rho no se midieron sobre el mismo conjunto de",
                  "ventanas. No son estrictamente comparables."]
        self.res("\n".join(L))


# ==========================================================================
#  5 - Seccion II-A: recocido simulado
# ==========================================================================
class TabSA(Pestana):
    titulo = "5 · II-A · Recocido"
    necesita_datos = True
    auto = False
    ancho_res = 52

    paper = [("Seccion II-A, literal:",
              "The SA process initiates with parameter vector x0 at "
              "T_init = 2000 and cooling coefficient alpha = 0.95. During each "
              "temperature stage T_k, the algorithm executes 100 Markov chain "
              "iterations by generating new states x_new = x_current + dx "
              "where dx_i ~ U(-0.1 d_i, 0.1 d_i) with d_i defining the "
              "parameter search range. State transitions follow the Metropolis "
              "criterion: new configurations are accepted with probability "
              "P = exp(-dE/T_k) where dE = E(x_new) - E(x_current) and E(x) "
              "represents the NCC-based matching error objective function. "
              "After completing Markov chains at current temperature, the "
              "system cools geometrically via T_k+1 = alpha*T_k "
              "(alpha in [0.85, 0.99]). This thermodynamic cycle iterates "
              "until meeting dual termination criteria: T_k < 1e-6 or energy "
              "stabilization criteria (|dE_mean| < 0.001 over 10 consecutive "
              "cycles)."),
             ("Ec. 14, en cambio, define la loss asi:",
              "loss = sum_i [S1(E_i) - S_model(E_i; p)]^2")]
    codigo = ("replica_paper/sa.py -> recocido(), multiarranque(),\n"
              "                       costo_ec14(), costo_ncc()")
    abiertas = [
        ("cual es el objetivo: la Ec. 14 o el error NCC",
         "El paper dice las dos cosas y son objetivos distintos. Estan los dos "
         "implementados y se elige con el control. Con el costo NCC, Pearson "
         "ignora la escala, asi que las amplitudes quedan literalmente "
         "INDETERMINADAS."),
        ("que significa T_init = 2000",
         "Metropolis compara T con dE, o sea T vive en las UNIDADES DEL COSTO, "
         "y el costo no esta fijado (punto anterior). La Ec. 14 es una suma de "
         "cuadrados sobre cuentas, asi que vale muchos ordenes de magnitud mas "
         "que el error NCC, que esta acotado en [0, 2]. El costo final de cada "
         "caso esta en la tabla del resultado: comparalo con T = 2000. 'auto' "
         "calibra T midiendo la dispersion del costo alrededor de x0."),
        ("de donde sale x0",
         "El paper no lo dice. Aca sale del NCC (pestana 4) y, con dos bordes, "
         "de un barrido grueso de (B1, D, sigma1) con las amplitudes exactas "
         "por nnls."),
        ("que es 'd_i, the parameter search range'",
         "Las cotas. Se derivan de la semilla del NCC: metodo.COTAS_REL."),
        ("que es dE_mean",
         "El paper no lo define. Aca es la diferencia entre el costo MEDIO de "
         "un ciclo de temperatura y el del ciclo anterior, que es la lectura "
         "literal de 'energy stabilization over 10 consecutive cycles'."),
    ]
    mirar = """
Panel de arriba: la corrida del recocido. La linea azul es el costo MEDIO de
cada ciclo de temperatura (lo que el criterio de corte mira), la roja el mejor
costo visto, y la escala de la derecha es T bajando geometricamente.

Panel de abajo: el ajuste que quedo, sobre la ventana del NCC.

Con 'comparar los 4 casos' se corre la misma cadena con los dos objetivos y con
T = 2000 y T = auto, y se mide la ACEPTACION. Ahi se ve de un vistazo por que
T_init = 2000 no es una constante transferible: es un numero con unidades, y el
paper no fija las unidades.

El SA es ESTOCASTICO: una sola corrida no es un resultado. Por eso corre desde
varios arranques (el paper reporta una) y el panel de resultado dice cuantos de
ellos llegaron al mismo optimo. Si esa fraccion es baja, el numero depende de
la suerte.
"""

    def construir_controles(self, f):
        self.v_obj = self.combo(f, "objetivo", ["Ec. 14 (cuadrados)",
                                                "error NCC (1 - rho)"],
                                ancho=20)
        self.v_t = self.entrada(f, "T_init", "2000")
        self.v_alfa = self.entrada(f, "alfa (enfriado)", str(mod_sa.ALFA))
        self.v_it = self.entrada(f, "iters por T", str(mod_sa.ITERS))
        self.v_paso = self.entrada(f, "paso_rel", str(mod_sa.PASO_REL))
        self.v_sem = self.entrada(f, "arranques", "3")
        self.v_sr = self.entrada(f, "sigma_rel de S2", "auto")
        self.v_pulir = self.check(f, "pulido local (el lazo de la Fig. 2)", True)
        self.v_cmp = self.check(f, "comparar los 4 casos (tarda 4x)", False)
        self.nota(f, "T_init acepta 'auto'. alfa entre 0.85 y 0.99 segun el "
                     "paper; el default 0.95 es el suyo.")

    def leer_controles(self):
        t = self.v_t.get().strip().lower()
        sr = self.v_sr.get().strip().lower()
        return dict(obj="ec14" if self.v_obj.get().startswith("Ec") else "ncc",
                    t=("auto" if t == "auto" else float(t)),
                    alfa=float(self.v_alfa.get()), it=int(self.v_it.get()),
                    paso=float(self.v_paso.get()), sem=int(self.v_sem.get()),
                    sr=(None if sr in ("auto", "", "none") else float(sr)),
                    pulir=bool(self.v_pulir.get()),
                    cmp=bool(self.v_cmp.get()))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        sa_kw = dict(historia=True, t_init=c["t"], alfa=c["alfa"],
                     iters=c["it"], paso_rel=c["paso"])
        # ligar_origen=True con dos lineas, igual que correr.py y que el resto
        # de la interfaz: la variante de D libre (Ec. 11) vive en la pestana 6,
        # que es donde se la compara con esta.
        extra = dict(objetivo=c["obj"], pulir=c["pulir"], semillas_sa=c["sem"],
                     sa_kw=sa_kw)
        if len(e["lineas"]) == 2:
            extra["ligar_origen"] = True
        self.progreso("NCC + recocido (%s)..." % c["obj"])
        r = metodo.ajustar(e["canales"], e["cuentas"], e["lineas"],
                           sigma_rel=c["sr"], **extra)
        comp = []
        if c["cmp"]:
            # La ventana y el sigma de la plantilla se FIJAN a los del ajuste de
            # arriba: si cambiaran, se estaria comparando tambien la eleccion de
            # ventana y no el objetivo ni la temperatura.
            for k, (obj, t0) in enumerate([("ec14", 2000.0), ("ec14", "auto"),
                                           ("ncc", 2000.0), ("ncc", "auto")]):
                self.progreso("comparando: objetivo %s, T = %s" % (obj, t0),
                              (k + 1) / 5.0)
                kw = dict(sa_kw)
                kw["t_init"] = t0
                otro = dict(extra)
                otro.update(objetivo=obj, pulir=False, semillas_sa=1, sa_kw=kw)
                q = metodo.ajustar(e["canales"], e["cuentas"], e["lineas"],
                                   intervalo=r["intervalo"],
                                   sigma_rel=r["sigma_rel_s2"], **otro)
                comp.append((obj, t0, q["sa"]))
        return dict(e=e, r=r, comp=comp, cfg=c)

    def pintar(self, d):
        e, r, c = d["e"], d["r"], d["cfg"]
        h = r["sa"].get("historia")
        a1 = self.fig.add_subplot(211)
        if h is not None and len(h):
            k = np.arange(1, h.shape[0] + 1)
            a1.plot(k, h[:, 1], color="tab:blue", lw=1.2,
                    label="costo medio del ciclo")
            a1.plot(k, h[:, 2], color="tab:red", lw=1.2, label="mejor costo")
            # El costo NCC (1 - rho) puede llegar a cero: la escala log lo
            # tiraria abajo. La Ec. 14 vive en 1e5-1e9 y sin log no se ve nada.
            if np.all(h[:, 1:] > 0):
                a1.set_yscale("log")
            a1.set_xlabel("ciclo de temperatura")
            a1.set_ylabel("costo (%s)" % ("Ec. 14" if c["obj"] == "ec14"
                                          else "1 - rho"))
            a1b = a1.twinx()
            a1b.plot(k, h[:, 0], color="0.55", ls="--", lw=1.0)
            if np.all(h[:, 0] > 0):
                a1b.set_yscale("log")
            a1b.set_ylabel("T", color="0.45")
            a1.grid(alpha=0.3)
            a1.legend(fontsize=8, loc="upper right")
            a1.set_title("recocido: %d ciclos, aceptacion %.3f, corte por %s"
                         % (r["sa"]["ciclos"], r["sa"]["aceptacion"],
                            r["sa"]["corte"]), fontsize=9)
        else:
            a1.text(0.5, 0.5, "sin historia", ha="center", transform=a1.transAxes)
            a1.set_axis_off()

        a2 = self.fig.add_subplot(212)
        a2.step(r["canales"], r["datos"], where="mid", lw=0.9, color="0.45",
                label="datos en la ventana")
        a2.plot(r["canales"], r["modelo"], color="tab:red", lw=1.6,
                label="modelo del SA")
        B = correr._borde_de(r)
        a2.axvline(B, color="k", ls="--", lw=1.1, label="B = %.1f" % B)
        a2.set_xlabel("canal [mV*ns]")
        a2.set_ylabel("cuentas netas")
        a2.grid(alpha=0.3)
        a2.legend(fontsize=8)

        L = ["objetivo      %s" % ("Ec. 14 (minimos cuadrados)"
                                   if c["obj"] == "ec14" else "1 - rho (NCC)"),
             "T_init        %s" % c["t"],
             "alfa          %.4f" % c["alfa"],
             "iters por T   %d" % c["it"],
             "paso_rel      %.3f" % c["paso"],
             "",
             "CORRIDA DEL SA",
             "  ciclos          %d" % r["sa"]["ciclos"],
             "  evaluaciones    %d" % r["sa"]["evaluaciones"],
             "  aceptacion      %.4f" % r["sa"]["aceptacion"],
             "  corte por       %s" % r["sa"]["corte"],
             "  T final         %.3e" % r["sa"]["T_final"],
             "",
             "ARRANQUES (el SA es estocastico)",
             "  arranques       %d" % r["dispersion"]["n_arranques"],
             "  en el optimo    %d de %d" % (r["dispersion"]["n_en_el_optimo"],
                                             r["dispersion"]["n_arranques"]),
             "  pulido local    %s" % ("si" if c["pulir"] else "NO"),
             "",
             "RESULTADO"]
        if r["tipo"] == "un_borde":
            L += ["  B               %.1f mV*ns" % r["B"],
                  "  sigma           %.1f  (sigma/B = %.4f)"
                  % (r["sigma"], r["sigma"] / r["B"]),
                  "  D = B/E_max     %.4f mV*ns/keV" % r["D"]]
        else:
            L += ["  B1              %.1f mV*ns" % r["B1"],
                  "  B2              %.1f mV*ns" % r["B2"],
                  "  D               %.4f mV*ns/keV" % r["D"],
                  "  sigma1, sigma2  %.1f, %.1f" % (r["sigma1"], r["sigma2"]),
                  "  A2/A1           %.3f" % (r["A2"] / max(r["A1"], 1e-30))]
        L += ["  chi2/ndf        %.3f" % r["chi2_ndf"], ""]
        if d["comp"]:
            L += ["-" * 48, "",
                  "LOS 4 CASOS (misma ventana, 1 arranque, sin pulir)", "",
                  " objetivo   T_init    aceptacion  ciclos    costo final",
                  ]
            for obj, t0, q in d["comp"]:
                L.append("  %-8s  %-8s   %8.4f   %5d   %.4e"
                         % ("Ec. 14" if obj == "ec14" else "NCC",
                            t0 if isinstance(t0, str) else "%.0f" % t0,
                            q["aceptacion"], q["ciclos"], q["costo"]))
            L += ["",
                  "Leer asi: compara la columna del COSTO contra",
                  "T = 2000. Metropolis acepta una subida dE con",
                  "probabilidad exp(-dE/T), asi que lo unico que",
                  "importa es el tamano relativo de los dos. La",
                  "Ec. 14 es una suma de cuadrados sobre cuentas y el",
                  "error NCC esta acotado en [0, 2]: la misma T no",
                  "puede significar lo mismo en los dos.",
                  "",
                  "La constante impresa solo tiene sentido en las",
                  "unidades del costo que se elija, y el paper no",
                  "elige. Lo mismo pasa con |dE_mean| < 0.001."]
        self.res("\n".join(L))


# ==========================================================================
#  6 - Fig. 2: la cadena completa
# ==========================================================================
class TabCadena(Pestana):
    titulo = "6 · Fig. 2 · La cadena"
    necesita_datos = True
    auto = False
    ancho_res = 54

    paper = [("Fig. 2, el diagrama de flujo completo:",
              "Start -> Energy Spectrum (Channel, Counts) -> Implement NCC "
              "Automatic Interval Selection -> Calculate Pearson Correlation "
              "Coefficient -> Maximum Data Interval of rho? -> Set Theta = "
              "{A1, A2, B1, sigma1, sigma2, D} -> Use SA to Solve Local "
              "Optimal Problems -> Optimal Parameter Solution? -> Fit "
              "Energy-Band Structure Relationship -> Stop"),
             ("Ec. 12, sobre la separacion de los dos bordes:",
              "CE_ch,2 = B1 + D*dE, con dE = E2 - E1 energy difference "
              "between gamma lines")]
    codigo = ("replica_paper/metodo.py -> ajustar(), ajustar_un_borde(),\n"
              "                           ajustar_dos_bordes()\n"
              "replica_paper/correr.py -> diagnostico_dos_bordes()")
    abiertas = [
        ("el lazo 'Optimal Parameter Solution? No -> modify'",
         "El diagrama lo dibuja y no dice que es. Aca es un pulido local "
         "(least_squares) desde CADA salida del SA, quedandose con la de menor "
         "loss. Sin ese lazo el resultado depende del ultimo salto aceptado."),
        ("que ventana usan los metodos convencionales",
         "El paper no da regla. Aca corren en la MISMA que eligio el NCC, que "
         "es lo unico que hace justa la comparacion de la Fig. 4."),
    ]
    mirar = """
LA EC. 12 IMPRESA ESTA MAL, y esta pestana lo deja probar. Dice que dE es la
diferencia entre las LINEAS GAMMA (159.3 keV en el Co-60). Pero D convierte
energia del ELECTRON a canal, asi que lo que se separa en canal son los BORDES:
154.7 keV. Cambiando el control se ve cuanto mueve a D.

Donde la errata se ve limpia no es aca sino en los datos publicados del paper
(pestana 10): sus dos bordes del Co-60 caen en los canales 469 y 543, o sea 74
canales, y la formula impresa predice 76.0 mientras la corregida predice 73.8.

EL HALLAZGO PROPIO DE CORRER ESTO SOBRE NUESTROS DATOS. La Ec. 11 pone D en
Theta, o sea LIBRE. Pero la Ec. 8 del mismo paper ya define el mapa por el
ORIGEN, y con ese mapa D = B1/E_max,1 queda determinado por B1. Dejar D libre
ademas sobre-parametriza el modelo respecto de su propia Ec. 8.

Tildar 'diagnostico del solapamiento' corre las dos variantes y las compara. Lo
que decide si el modelo de dos bordes sirve es la separacion de los bordes en
unidades de sigma: el cierre con verdad conocida de metodo.py dice que a 2
sigma D sale sesgada ~-14% con +-15% de dispersion, y que recien a 4 sigma es
utilizable.

En el panel de residuos: (dato - modelo) / sqrt(modelo), que es lo mismo que
suma el chi2/ndf que se reporta.
"""

    def construir_controles(self, f):
        self.v_sr = self.entrada(f, "sigma_rel de S2", "auto")
        self.v_iv = self.entrada(f, "intervalo", "auto (NCC)", ancho=14)
        self.v_obj = self.combo(f, "objetivo", ["Ec. 14 (cuadrados)",
                                                "error NCC (1 - rho)"], ancho=20)
        self.v_sem = self.entrada(f, "arranques del SA", "auto")
        self.v_pulir = self.check(f, "pulido local", True)
        self.separador(f)
        ttk.Label(f, text="solo con DOS lineas gamma:").pack(
            side="top", anchor="w", padx=8)
        self.v_de = self.combo(f, "Ec. 12, dE de", ["bordes (corregida)",
                                                    "gammas (impresa)"], ancho=20)
        self.v_lig = self.check(f, "D = B1/E_max,1 (lo que impone la Ec. 8)",
                                True)
        self.v_diag = self.check(f, "diagnostico del solapamiento", True)
        self.nota(f, "intervalo: 'auto (NCC)' o dos numeros separados por coma. "
                     "Con una sola linea la Ec. 12 y las ligaduras no aplican.")

    def leer_controles(self):
        s = self.v_iv.get().strip().lower()
        iv = None
        if s and not s.startswith("auto"):
            p = [float(x) for x in s.replace(";", ",").split(",")]
            iv = (p[0], p[1])
        sr = self.v_sr.get().strip().lower()
        sem = self.v_sem.get().strip().lower()
        return dict(sr=(None if sr in ("auto", "", "none") else float(sr)),
                    iv=iv,
                    obj="ec14" if self.v_obj.get().startswith("Ec") else "ncc",
                    sem=(None if sem in ("auto", "", "none") else int(sem)),
                    pulir=bool(self.v_pulir.get()),
                    de="bordes" if self.v_de.get().startswith("bordes")
                    else "gammas",
                    lig=bool(self.v_lig.get()), diag=bool(self.v_diag.get()))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        # 'auto' = el default de metodo.py, que es 3 arranques con un borde y 8
        # con dos. Con dos el SA no llega desde un solo arranque (ver el demo
        # de metodo.py); con uno, 8 solo cuesta tiempo.
        extra = dict(objetivo=c["obj"], pulir=c["pulir"])
        if c["sem"] is not None:
            extra["semillas_sa"] = c["sem"]
        if len(e["lineas"]) == 2:
            extra.update(delta_e=c["de"], ligar_origen=c["lig"])
        self.progreso("NCC -> SA -> pulido (Fig. 2)...")
        bordes, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"],
                                 intervalo=c["iv"], sigma_rel=c["sr"], **extra)
        diag = None
        if len(e["lineas"]) == 2 and c["diag"]:
            self.progreso("diagnostico del solapamiento (3 ajustes mas)...")
            diag = correr.diagnostico_dos_bordes(e, r, verbose=False)
        return dict(e=e, r=r, bordes=bordes, diag=diag, cfg=c)

    def pintar(self, d):
        e, r, bordes = d["e"], d["r"], d["bordes"]
        gs = self.fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
        a1 = self.fig.add_subplot(gs[0])
        a2 = self.fig.add_subplot(gs[1], sharex=a1)

        a1.step(e["canales"], e["cuentas"], where="mid", lw=0.8, color="0.55",
                label="neto")
        a1.axvspan(r["intervalo"][0], r["intervalo"][1], color="tab:orange",
                   alpha=0.12, label="ventana del NCC")
        a1.plot(r["canales"], r["modelo"], color="tab:red", lw=1.8,
                label="modelo (Ec. 8-13) del SA")
        if r["tipo"] == "dos_bordes":
            c1, c2, _ = ec.componentes_dos_bordes(
                r["canales"], r["e_gamma"][0], r["e_gamma"][1], r["A1"],
                r["A2"], r["B1"], r["sigma1"], r["sigma2"], r["D"],
                delta_e=r["delta_e"])
            a1.plot(r["canales"], c1, "--", color="tab:purple", lw=1.1,
                    label="componente %.1f keV" % r["E_borde"][0])
            a1.plot(r["canales"], c2, "--", color="tab:brown", lw=1.1,
                    label="componente %.1f keV" % r["E_borde"][1])
            a1.axvline(r["B2"], color="tab:brown", lw=1.0, ls=":")
        for m in correr.METODOS:
            if np.isfinite(bordes[m]):
                a1.axvline(bordes[m], color=COL[m], ls=LSS[m], lw=1.4,
                           label="%s: %.1f" % (m, bordes[m]))
        a1.set_ylabel("cuentas netas")
        a1.grid(alpha=0.3)
        a1.legend(fontsize=7.5, ncol=2)
        a1.tick_params(labelbottom=False)
        a1.set_title("%s  (%s)   chi2/ndf = %.3f"
                     % (e["clave"], e["punto"], r["chi2_ndf"]), fontsize=10)

        var = np.maximum(np.abs(r["modelo"]), 1.0)
        rr = (r["datos"] - r["modelo"]) / np.sqrt(var)
        a2.axhline(0, color="0.5", lw=0.8)
        a2.axhspan(-1, 1, color="0.85")
        a2.step(r["canales"], rr, where="mid", lw=0.8, color="k")
        a2.set_xlabel("canal [mV*ns]")
        a2.set_ylabel("residuo\nnormalizado")
        a2.grid(alpha=0.3)

        emax = [float(ec.e_max(x)) for x in e["lineas"]]
        L = ["corrida     %s  (%s)" % (e["clave"], e["punto"]),
             "fuente      %s  ->  bordes %s keV"
             % (e["fuente"], ", ".join("%.1f" % x for x in emax)),
             "",
             "1) NCC (Seccion II-B)",
             "   intervalo    [%.0f, %.0f] mV*ns"
             % (r["intervalo"][0], r["intervalo"][1]),
             "   rho          %.5f" % r["rho"],
             "   sigma_rel S2 %.3f" % r["sigma_rel_s2"],
             "",
             "2) SA (Seccion II-A) + pulido"]
        if r["tipo"] == "un_borde":
            L += ["   B            %.1f mV*ns" % r["B"],
                  "   sigma        %.1f  (sigma/B = %.4f -> %.1f%% FWHM)"
                  % (r["sigma"], r["sigma"] / r["B"],
                     100 * 2.3548 * r["sigma"] / r["B"]),
                  "   D = B/E_max  %.4f mV*ns/keV" % r["D"]]
        else:
            sep = (r["B2"] - r["B1"]) / max(r["sigma1"], 1e-9)
            L += ["   B1           %.1f mV*ns" % r["B1"],
                  "   B2           %.1f mV*ns" % r["B2"],
                  "   D            %.4f mV*ns/keV" % r["D"],
                  "   sigma1,2     %.1f, %.1f" % (r["sigma1"], r["sigma2"]),
                  "   A2/A1        %.3f" % (r["A2"] / max(r["A1"], 1e-30)),
                  "   separacion   %.2f sigma1" % sep,
                  "   Ec. 12 con   dE de %s = %.1f keV"
                  % (r["delta_e"], r["dE"]),
                  "   D ligada     %s" % ("si (Ec. 8)" if r["ligar_origen"]
                                          else "NO (Ec. 11)")]
        L += ["   chi2/ndf     %.3f" % r["chi2_ndf"],
              "   arranques    %d, en el optimo %d"
              % (r["dispersion"]["n_arranques"],
                 r["dispersion"]["n_en_el_optimo"]),
              "",
              "3) los tres metodos, MISMA ventana",
              "   %-12s %10s   contra NCC-SA" % ("metodo", "borde")]
        b0 = bordes["NCC-SA"]
        for m in correr.METODOS:
            v = bordes[m]
            dif = ("      -" if m == "NCC-SA"
                   else ("   %+6.1f%%" % (100 * (v / b0 - 1))
                         if np.isfinite(v) and b0 else "     nan"))
            L.append("   %-12s %10.1f%s" % (m, v, dif))
        L += ["", "-" * 50, ""]
        L += _tabla_barrido(r["barrido_sigma_rel"], r["sigma_rel_s2"])

        if d["diag"] is not None:
            g = d["diag"]
            L += ["", "-" * 50, "",
                  "DIAGNOSTICO DEL SOLAPAMIENTO (Ec. 11-13)", "",
                  " variante           B1       B2       D      A2/A1  chi2/ndf"]
            for k, nombre in ((False, "D libre (Ec. 11)"),
                              (True, "D = B1/Emax (Ec.8)"),
                              ("gammas", "Ec. 12 impresa")):
                f = g[k]
                L.append(" %-18s %7.1f  %7.1f  %.4f  %6.2f   %.2f"
                         % (nombre, f["B1"], f["B2"], f["D"],
                            f["A2"] / max(f["A1"], 1e-30), f["chi2_ndf"]))
            dif = 100 * abs(g[False]["D"] / g[True]["D"] - 1)
            sep = (g[True]["B2"] - g[True]["B1"]) / max(g[True]["sigma1"], 1e-9)
            L += ["",
                  " D cambia %.0f%% segun se ligue o no, con chi2" % dif,
                  " casi igual. Y el cociente de amplitudes toma",
                  " valores sin sentido fisico: el Co-60 emite las",
                  " dos lineas con intensidades practicamente",
                  " iguales (la simulacion da 0.95).",
                  "",
                  " Los dos bordes quedan a %.2f sigma1." % sep,
                  " metodo.py, con verdad conocida, dice que a 2",
                  " sigma D sale sesgada -14% +- 15% y que recien",
                  " a 4 sigma es utilizable."]
        self.res("\n".join(L))


# ==========================================================================
#  7 - los dos metodos convencionales
# ==========================================================================
class TabConvencionales(Pestana):
    titulo = "7 · Convencionales"
    necesita_datos = True
    auto = False
    ancho_res = 54

    paper = [("Metodo 1, 'Gaussian FWHM Method':",
              "The first method employs Gaussian function fitting to the "
              "spectral descent slope, defining the channel position at 50% of "
              "the maximum amplitude as the Compton edge"),
             ("Metodo 2, 'Difference Method':",
              "The second method calculates the first derivative of the "
              "smoothed measured energy spectrum, identifying the channel "
              "corresponding to the derivative minimum as the Compton edge"),
             ("Y el diagnostico que el paper hace de ellos:",
              "the assumed ideal Gaussian distribution significantly deviates "
              "from the observed asymmetric broadening, leading to systematic "
              "errors in edge localization")]
    codigo = ("replica_paper/convencionales.py -> metodo_gauss(),\n"
              "        metodo_diferencia(), region_flanco(),\n"
              "        sensibilidad_ventana()")
    abiertas = [
        ("con que se suaviza y con que ventana",
         "Savitzky-Golay, que es lo estandar y lo que la Fig. 5 del paper "
         "muestra como 'Smoothed Data'. La ventana es un control, y la tabla "
         "de sensibilidad mide cuanto la mueve."),
        ("que region usan los dos metodos",
         "El paper no da regla. Por defecto se les da la MISMA que eligio el "
         "NCC: si cada uno eligiera la suya se estaria comparando la eleccion "
         "de ventana y no el metodo. 'region_flanco' es la alternativa "
         "automatica que no depende del NCC."),
        ("que es 'the maximum amplitude'",
         "Aca es la amplitud de la gaussiana AJUSTADA, que es la unica "
         "cantidad bien definida despues de un ajuste gaussiano. El borde "
         "queda en mu + sigma*sqrt(2 ln(1/fraccion))."),
        ("la fraccion",
         "El texto dice 50% y el resumen dice '50% 70%'. Queda como control."),
    ]
    mirar = """
EL SESGO DE ESTOS DOS METODOS NO ES ESTADISTICO, ES DEL METODO. Con 'sesgo con
verdad conocida' se genera un borde en una posicion CONOCIDA y se mide cuanto
se corren. La ultima columna de esa tabla es el mismo metodo sobre el espectro
SIN RUIDO: si el sesgo fuera estadistico, ahi seria cero, y no lo es.

La causa es comun a los dos: la lineshape de Klein-Nishina no es un escalon,
tiene una acumulacion creciente ANTES del borde (se ve en la pestana 2), y ni
un escalon gaussiano ni el minimo de la derivada la contienen. Esa es la
justificacion cuantitativa de ajustar la lineshape entera, que es lo que hace
el paper.

Panel de arriba: la gaussiana ajustada al flanco y el punto al 50% de su
amplitud. Panel del medio: el espectro suavizado y su derivada, con el minimo.
Panel de abajo (si se pidio): el sesgo contra la resolucion.

La curva de la diferencia CON ruido salta de un signo al otro mientras la
misma curva sin ruido sube suave. No es contradiccion: la derivada del borde es
un pozo cuyo ancho escala con la resolucion, y la ventana de suavizado aca es
fija (el paper no dice otra cosa). Cuando el pozo se ensancha mas que la
ventana aparecen minimos locales casi iguales y el ruido elige cual gana. La
tabla de sensibilidad a la ventana, en el resultado, mide lo mismo.

DOS METODOS ESTABLES PUEDEN SER LOS DOS INEXACTOS. En el Co-60 de este proyecto
el metodo gaussiano tiene MUY poca dispersion entre bloques y esta corrido; el
NCC-SA dispersa mucho mas. Poca dispersion no es exactitud, y esa distincion es
justamente la que la Fig. 4 del paper mezcla en un solo grafico.
"""

    SIGMAS_SINT = (0.03, 0.06, 0.09, 0.12, 0.18)

    def construir_controles(self, f):
        self.v_reg = self.combo(f, "region", ["la del NCC-SA",
                                              "region_flanco (sin NCC)",
                                              "manual"], ancho=20)
        self.v_man = self.entrada(f, "manual lo,hi", "600, 1400", ancho=14)
        self.v_fr = self.entrada(f, "fraccion", "0.5")
        self.v_ven = self.entrada(f, "ventana Sav-Golay", str(conv.VENTANA))
        self.v_ord = self.entrada(f, "orden", str(conv.ORDEN))
        self.v_sint = self.check(f, "sesgo con verdad conocida", True)
        self.nota(f, "El sesgo sintetico usa un borde en el canal 600 con "
                     "E_gamma = 661.657 keV y ruido de Poisson, igual que el "
                     "demo de convencionales.py.")

    def leer_controles(self):
        p = [float(x) for x in self.v_man.get().replace(";", ",").split(",")]
        return dict(reg=self.v_reg.get(), man=(p[0], p[1]),
                    fr=float(self.v_fr.get()), ven=int(self.v_ven.get()),
                    orden=int(self.v_ord.get()), sint=bool(self.v_sint.get()))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        r = None
        if c["reg"].startswith("la del"):
            self.progreso("corriendo el NCC-SA para tomarle la ventana...")
            _, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"])
            region = r["intervalo"]
        elif c["reg"].startswith("region_flanco"):
            region = conv.region_flanco(e["canales"], e["cuentas"],
                                        ventana=c["ven"])
        else:
            region = c["man"]
        g = conv.metodo_gauss(e["canales"], e["cuentas"], region, c["fr"])
        df = conv.metodo_diferencia(e["canales"], e["cuentas"], region,
                                    ventana=c["ven"], orden=c["orden"])
        sens = conv.sensibilidad_ventana(e["canales"], e["cuentas"], region)
        sint = []
        if c["sint"]:
            self.progreso("midiendo el sesgo con verdad conocida...")
            rng = np.random.default_rng(5)
            eg, B_real = 661.657, 600.0
            ch = np.linspace(2.5, 1400.0, 280)
            for k, sr in enumerate(self.SIGMAS_SINT):
                limpio = ec.modelo(ch, eg, 1.0, B_real, sr * B_real)
                limpio = limpio / limpio.max() * 800.0
                dd = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
                reg = conv.region_flanco(ch, dd)
                reg0 = conv.region_flanco(ch, limpio)
                sint.append((sr,
                             conv.metodo_gauss(ch, dd, reg, 0.5)["borde"],
                             conv.metodo_gauss(ch, dd, reg, 0.7)["borde"],
                             conv.metodo_diferencia(ch, dd, reg)["borde"],
                             conv.metodo_gauss(ch, limpio, reg0, 0.5)["borde"],
                             conv.metodo_diferencia(ch, limpio, reg0)["borde"],
                             B_real))
                self.progreso("sesgo sintetico %d/%d"
                              % (k + 1, len(self.SIGMAS_SINT)),
                              (k + 1) / len(self.SIGMAS_SINT))
        return dict(e=e, r=r, g=g, df=df, region=region, sens=sens, sint=sint,
                    cfg=c)

    def pintar(self, d):
        e, g, df, c = d["e"], d["g"], d["df"], d["cfg"]
        n = 3 if d["sint"] else 2
        a1 = self.fig.add_subplot(n, 1, 1)
        a2 = self.fig.add_subplot(n, 1, 2)

        lo, hi = d["region"]
        m = (e["canales"] >= lo) & (e["canales"] <= hi)
        a1.step(e["canales"], e["cuentas"], where="mid", lw=0.8, color="0.6",
                label="neto")
        x = e["canales"][m]
        a1.plot(x, g["A"] * np.exp(-(x - g["mu"]) ** 2 / (2 * g["sigma"] ** 2)),
                color="tab:green", lw=1.7, label="gaussiana al flanco")
        a1.axhline(c["fr"] * g["A"], color="tab:green", ls=":", lw=1.0,
                   label="%.0f%% de la amplitud" % (100 * c["fr"]))
        a1.axvline(g["borde"], color="tab:green", ls="--", lw=1.4,
                   label="Gauss FWHM: %.1f" % g["borde"])
        if d["r"] is not None:
            a1.axvline(correr._borde_de(d["r"]), color="tab:red", lw=1.4,
                       label="NCC-SA: %.1f" % correr._borde_de(d["r"]))
        a1.axvspan(lo, hi, color="tab:orange", alpha=0.10, label="region")
        a1.set_ylabel("cuentas netas")
        a1.grid(alpha=0.3)
        a1.legend(fontsize=7.5, ncol=2)
        a1.set_title("metodo 1: gaussiana al flanco", fontsize=10)

        a2.plot(df["canales"], df["suave"], color="0.35", lw=1.2,
                label="suavizado (Sav-Golay %d)" % df["ventana"])
        a2.set_ylabel("cuentas")
        a2.grid(alpha=0.3)
        a2b = a2.twinx()
        a2b.plot(df["canales"], df["derivada"], color="tab:blue", lw=1.1,
                 label="derivada")
        a2b.axvline(df["borde"], color="tab:blue", ls=":", lw=1.6)
        a2b.set_ylabel("d(cuentas)/d(canal)", color="tab:blue")
        a2b.tick_params(axis="y", labelcolor="tab:blue")
        a2.set_xlabel("canal [mV*ns]")
        a2.set_title("metodo 2: minimo de la derivada -> %.1f mV*ns"
                     % df["borde"], fontsize=10)
        a2.legend(fontsize=8, loc="upper left")

        if d["sint"]:
            a3 = self.fig.add_subplot(n, 1, 3)
            sr = [t[0] for t in d["sint"]]
            p = lambda v, B: 100 * (v / B - 1)
            a3.plot(sr, [p(t[1], t[6]) for t in d["sint"]], "o-",
                    color="tab:green", label="gauss 50%")
            a3.plot(sr, [p(t[3], t[6]) for t in d["sint"]], "s-",
                    color="tab:blue", label="diferencia")
            a3.plot(sr, [p(t[4], t[6]) for t in d["sint"]], "o--",
                    color="tab:green", alpha=0.5, label="gauss 50%, sin ruido")
            a3.plot(sr, [p(t[5], t[6]) for t in d["sint"]], "s--",
                    color="tab:blue", alpha=0.5, label="diferencia, sin ruido")
            a3.axhline(0, color="k", lw=0.8)
            a3.set_xlabel("sigma/B del borde sintetico (resolucion)")
            a3.set_ylabel("sesgo [%]")
            a3.grid(alpha=0.3)
            a3.legend(fontsize=7.5, ncol=2)
            a3.set_title("el sesgo es del METODO: sin ruido tampoco es cero",
                         fontsize=10)

        L = ["corrida    %s  (%s)" % (e["clave"], e["punto"]),
             "region     [%.0f, %.0f] mV*ns  (%s)" % (lo, hi, c["reg"]),
             "",
             "METODO 1 - gaussiana al flanco",
             "   A       %.1f" % g["A"],
             "   mu      %.1f" % g["mu"],
             "   sigma   %.1f" % g["sigma"],
             "   borde = mu + sigma*sqrt(2 ln(1/%.2f)) = %.1f"
             % (c["fr"], g["borde"]),
             "",
             "METODO 2 - minimo de la derivada",
             "   ventana Sav-Golay %d, orden %d" % (c["ven"], c["orden"]),
             "   borde   %.1f  (derivada = %.2f)" % (df["borde"],
                                                     df["derivada_min"]),
             ""]
        if d["r"] is not None:
            b0 = correr._borde_de(d["r"])
            L += ["NCC-SA (misma ventana)   %.1f" % b0,
                  "   Gauss FWHM  %+6.2f%%" % (100 * (g["borde"] / b0 - 1)),
                  "   Diferencia  %+6.2f%%" % (100 * (df["borde"] / b0 - 1)),
                  ""]
        L += ["-" * 50, "",
              "sensibilidad a la ventana de suavizado",
              "(el paper no la especifica):",
              "   ventana    borde"]
        for v, b in d["sens"]:
            L.append("   %5d    %9.1f" % (v, b))
        if d["sint"]:
            L += ["", "-" * 50, "",
                  "SESGO CON VERDAD CONOCIDA  (B real = 600)", "",
                  " sigma/B  gauss50  gauss70   dif    | sin ruido: g50 / dif"]
            for sr, g5, g7, dd, g0, d0, B in d["sint"]:
                p = lambda v: 100 * (v / B - 1)
                L.append("  %.2f   %+6.1f%%  %+6.1f%%  %+6.1f%%  |  %+6.1f%% / %+6.1f%%"
                         % (sr, p(g5), p(g7), p(dd), p(g0), p(d0)))
            L += ["",
                  "La columna de la derecha es el MISMO metodo sobre",
                  "el espectro sin ruido: el sesgo no es estadistico."]
        self.res("\n".join(L))


# ==========================================================================
#  8 - Fig. 5: robustez frente al intervalo de ajuste
# ==========================================================================
class TabRobustez(Pestana):
    titulo = "8 · Fig. 5 · Robustez"
    necesita_datos = True
    auto = False
    ancho_res = 56

    paper = [("Seccion III, sobre la Fig. 5:",
              "the calibration results remained highly consistent"),
             ("Y el epigrafe de la propia Fig. 5 usa BC408, mientras el texto "
              "dice:", "calibrating EJ309 spectra (137Cs source)")]
    codigo = ("replica_paper/correr.py -> INTERVALOS_FIG5, tres_metodos()\n"
              "replica_paper/resultados.py -> robustez_intervalo()\n"
              "replica_paper/auditoria.py -> fig5()  (los suyos)")
    abiertas = [
        ("cuanto es 'highly consistent'",
         "El paper no da numero. El que sale de sus propios seis paneles esta "
         "en el resultado, calculado igual que el nuestro."),
        ("hasta donde se puede achicar un intervalo",
         "min_bins. Los intervalos mas angostos del paper cubren el 6% del "
         "canal del borde y sobre sus espectros son ~100 bins; sobre los "
         "nuestros, ~30. Por debajo de min_bins no hay con que ajustar tres "
         "parametros ni con que derivar, asi que se saltean, y el mismo corte "
         "se aplica a los tres metodos para no favorecer a ninguno."),
    ]
    mirar = """
Los seis intervalos de la Fig. 5 del paper estan en canales absolutos de SU
detector ([1500,1600] ... [800,2200] con el borde en 1535). Aca se los pasa a
unidades del borde y se los aplica al nuestro, que es la unica forma de que la
comparacion signifique algo.

EN LOS DOS CASOS EL INTERVALO MAS ANGOSTO ES EL QUE SE APARTA, y ahi esta casi
toda la dispersion. Por eso el resultado reporta las dos cifras: con todos los
intervalos y sin el mas angosto. El paper reporta cero cifras.

Lo que se ve en el grafico de arriba: cada metodo contra el indice del
intervalo (1 = el mas angosto). Una linea plana es un metodo robusto frente a
la eleccion de ventana. Ojo: robusto no quiere decir exacto. La pestana 7 mide
lo otro.
"""

    def construir_controles(self, f):
        self.v_min = self.entrada(f, "min_bins", str(correr.MIN_BINS_ROBUSTEZ))
        self.v_sr = self.entrada(f, "sigma_rel de S2", "auto")
        self.nota(f, "Corre el NCC-SA una vez para tener el borde de "
                     "referencia, y despues los tres metodos en cada uno de "
                     "los seis intervalos: siete cadenas completas. Medido: "
                     "~2 min en el Cs-137. Las ventanas angostas son las "
                     "lentas, porque el SA corre muchos mas ciclos.")

    def leer_controles(self):
        sr = self.v_sr.get().strip().lower()
        return dict(minb=int(self.v_min.get()),
                    sr=(None if sr in ("auto", "", "none") else float(sr)))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        self.progreso("ajuste de referencia (para tener el borde)...")
        bordes, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"],
                                 sigma_rel=c["sr"])
        B = bordes["NCC-SA"]
        paso = float(e["canales"][1] - e["canales"][0])
        filas, saltados = [], []
        n = len(correr.INTERVALOS_FIG5)
        for k, (a_, b_) in enumerate(correr.INTERVALOS_FIG5):
            iv = (a_ * B, b_ * B)
            nb = int(round((iv[1] - iv[0]) / paso))
            if nb < c["minb"]:
                saltados.append((k, iv, nb))
                continue
            self.progreso("intervalo %d/%d: [%.0f, %.0f]"
                          % (k + 1, n, iv[0], iv[1]), (k + 1) / float(n))
            try:
                bb, _ = tres_metodos(e["canales"], e["cuentas"], e["lineas"],
                                     intervalo=iv, sigma_rel=r["sigma_rel_s2"])
            except Exception:                                # noqa: BLE001
                bb = {m: float("nan") for m in correr.METODOS}
            # k es el indice del panel de la Fig. 5 del paper. Se guarda porque
            # si algun intervalo se saltea por min_bins, numerar 1..N los que
            # quedaron desalinearia la comparacion con sus seis paneles.
            filas.append((k, iv, bb, nb))
        rob = {}
        ivs = [f[1] for f in filas]
        for m in correr.METODOS:
            tabla = {tuple(iv): bb[m] for _, iv, bb, _ in filas}
            rob[m] = R.robustez_intervalo(lambda iv, t=tabla: t[tuple(iv)], ivs)
            v = rob[m]["valores"][1:]
            ok = np.isfinite(v)
            rob[m]["sin_el_mas_angosto"] = (
                float((v[ok].max() - v[ok].min()) / np.mean(v[ok]))
                if ok.sum() > 1 else float("nan"))
        return dict(e=e, r=r, bordes=bordes, filas=filas, saltados=saltados,
                    rob=rob, B=B, minb=c["minb"],
                    texto_paper=texto_de(auditoria.fig5))

    def pintar(self, d):
        e = d["e"]
        a1 = self.fig.add_subplot(211)
        a2 = self.fig.add_subplot(212)
        idx = [f[0] + 1 for f in d["filas"]]
        for m in correr.METODOS:
            v = np.array([b for _, b in d["rob"][m]["filas"]], float)
            a1.plot(idx, v / d["B"], "o-", color=COL[m],
                    lw=1.2, ms=5, label="%s  (disp %.1f%%)"
                    % (m, 100 * d["rob"][m]["dispersion_rel"]))
        vp = np.array([b for _, b in auditoria.FIG5], float)
        a1.plot(range(1, vp.size + 1), vp / vp.mean(), "^--", color="0.45",
                lw=1.0, ms=5,
                label="los 6 paneles del paper (disp %.2f%%)"
                % (100 * (vp.max() - vp.min()) / vp.mean()))
        for j, (k, _, nb) in enumerate(d["saltados"]):
            a1.plot([k + 1], [1.0], "x", color="0.55", ms=10, mew=2,
                    label=("salteado: menos de %d bins" % d["minb"])
                    if j == 0 else None)
        a1.axhline(1.0, color="0.7", lw=0.8)
        a1.set_xlabel("intervalo de la Fig. 5   (1 = el mas angosto)")
        a1.set_ylabel("borde / borde de referencia")
        a1.grid(alpha=0.3)
        a1.legend(fontsize=7.5)
        a1.set_title("%s: cuanto mueve el borde la eleccion de la ventana"
                     % e["clave"], fontsize=10)

        a2.step(e["canales"], e["cuentas"], where="mid", lw=0.8, color="0.55")
        cm = matplotlib.colormaps["viridis"]
        for j, (k, iv, _, nb) in enumerate(d["filas"]):
            y = (0.92 - 0.09 * j) * np.nanmax(e["cuentas"])
            a2.plot(iv, [y, y], lw=3.0,
                    color=cm(j / max(len(d["filas"]) - 1, 1)))
            a2.text(iv[1], y, "  #%d, %d bins" % (k + 1, nb), fontsize=7,
                    va="center")
        a2.axvline(d["B"], color="tab:red", lw=1.3,
                   label="borde de referencia %.1f" % d["B"])
        a2.set_xlabel("canal [mV*ns]")
        a2.set_ylabel("cuentas netas")
        a2.grid(alpha=0.3)
        a2.legend(fontsize=8)

        L = ["borde de referencia (NCC-SA)  %.1f mV*ns" % d["B"],
             "bin del histograma            %.2f mV*ns"
             % (e["canales"][1] - e["canales"][0]),
             "",
             " #  intervalo [mV*ns]     bins   " +
             "".join("%-12s" % m for m in correr.METODOS)]
        for k, iv, bb, nb in d["filas"]:
            L.append(" %d [%7.0f, %7.0f]  %5d   " % (k + 1, iv[0], iv[1], nb) +
                     "".join("%-12.1f" % bb[m] for m in correr.METODOS))
        for k, iv, nb in d["saltados"]:
            L.append(" %d [%7.0f, %7.0f]  %5d   salteado (< min_bins)"
                     % (k + 1, iv[0], iv[1], nb))
        L += ["",
              " %-22s" % "dispersion (todos)" +
              "".join("%-12s" % ("%.2f%%" % (100 * d["rob"][m]["dispersion_rel"]))
                      for m in correr.METODOS),
              " %-22s" % "  sin el mas angosto" +
              "".join("%-12s" % ("%.2f%%" % (100 * d["rob"][m]["sin_el_mas_angosto"]))
                      for m in correr.METODOS),
              "", "-" * 54, "",
              "LOS SEIS PANELES DEL PAPER, con el mismo criterio:",
              d["texto_paper"]]
        self.res("\n".join(L))


# ==========================================================================
#  9 - Seccion II-C: las "diez mediciones independientes"
# ==========================================================================
class TabBloques(Pestana):
    titulo = "9 · II-C · Las 10 mediciones"
    necesita_datos = True
    auto = False
    ancho_res = 56

    paper = [("Seccion II-C:",
              "Ten independent pulse-height spectrum measurements per "
              "scintillator-source combination yielded 80 datasets ... error "
              "bars representing standard deviations across measurements")]
    codigo = ("replica_paper/datos.py -> bloques()\n"
              "replica_paper/correr.py -> analizar() (la parte de bloques)")
    abiertas = [
        ("nosotros tenemos UNA corrida por fuente, no diez",
         "bloques() parte los eventos de la corrida en n tramos consecutivos y "
         "los trata como las diez mediciones. Es lo unico que se puede hacer "
         "con lo medido, y NO es lo mismo."),
    ]
    mirar = """
OJO CON QUE PRUEBA ESTO. La dispersion entre bloques es SOLO la estadistica de
conteo. Las diez mediciones del paper se tomaron por separado, asi que sus
barras incluyen ademas la variacion de montaje, de posicion de la fuente y de
deriva de ganancia. La barra de aca es una COTA INFERIOR de la de ellos.

Y va en la otra direccion tambien: cada bloque tiene la DECIMA parte de las
cuentas, asi que su dispersion estadistica es ~sqrt(10) veces la de una
medicion completa. Por eso el resultado da el SD crudo y el SD/sqrt(n), y
ninguno de los dos es directamente su barra.

QUE MIRAR DE VERDAD: cual de los tres metodos dispersa mas. En el Cs-137 (una
sola linea) el NCC-SA y el gaussiano dispersan parecido y el de la diferencia
seis veces mas. En el Co-60 (dos bordes solapados) el NCC-SA se va al 40%,
que es la senal de que Theta esta sub-determinado a esa resolucion, no de que
el optimizador ande mal.

La ventana y el sigma de la plantilla se FIJAN a los de la corrida completa
para todos los bloques. Si cada bloque eligiera su propia ventana, la
dispersion mediria las dos cosas juntas.
"""

    def construir_controles(self, f):
        self.v_n = self.entrada(f, "bloques", "10")
        self.v_sr = self.entrada(f, "sigma_rel de S2", "auto")
        self.nota(f, "Es el calculo mas caro de la interfaz: son n+1 cadenas "
                     "NCC-SA completas. En el Cs-137 tarda ~1 min y en el "
                     "Co-60 ~8 min.")

    def leer_controles(self):
        sr = self.v_sr.get().strip().lower()
        return dict(n=int(self.v_n.get()),
                    sr=(None if sr in ("auto", "", "none") else float(sr)))

    def trabajo(self):
        c = self.cfg
        e = self.espectro()
        self.progreso("corrida completa (referencia)...")
        bordes, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"],
                                 sigma_rel=c["sr"])
        self.progreso("partiendo la corrida en %d bloques..." % c["n"])
        bl = dat.bloques(e["clave"], n=c["n"], bins=e["bins"],
                         rango_max=e["rango_max"])
        vals = {m: [] for m in correr.METODOS}
        espectros = []
        for k, b in enumerate(bl):
            self.progreso("bloque %d/%d" % (k + 1, c["n"]), (k + 1) / float(c["n"]))
            try:
                bb, _ = tres_metodos(b["canales"], b["cuentas"], b["lineas"],
                                     intervalo=r["intervalo"],
                                     sigma_rel=r["sigma_rel_s2"])
            except Exception:                                # noqa: BLE001
                bb = {m: float("nan") for m in correr.METODOS}
            for m in correr.METODOS:
                vals[m].append(bb[m])
            espectros.append((b["canales"], b["cuentas"]))
        disp = {}
        for m in correr.METODOS:
            v = np.array(vals[m], float)
            ok = np.isfinite(v)
            disp[m] = dict(valores=v,
                           media=float(np.mean(v[ok])) if ok.any() else np.nan,
                           sd=float(np.std(v[ok], ddof=1)) if ok.sum() > 1 else np.nan,
                           n=int(ok.sum()))
        return dict(e=e, r=r, bordes=bordes, disp=disp, espectros=espectros,
                    cfg=c)

    def pintar(self, d):
        e, c = d["e"], d["cfg"]
        a1 = self.fig.add_subplot(211)
        a2 = self.fig.add_subplot(212)
        for m in correr.METODOS:
            v = d["disp"][m]["valores"]
            k = np.arange(1, v.size + 1)
            a1.plot(k, v, "o-", color=COL[m], lw=1.0, ms=5, label=m)
            a1.axhline(d["bordes"][m], color=COL[m], ls=":", lw=1.0)
            s = d["disp"][m]
            if np.isfinite(s["sd"]):
                a1.fill_between([0.5, v.size + 0.5],
                                s["media"] - s["sd"], s["media"] + s["sd"],
                                color=COL[m], alpha=0.10)
        a1.set_xlabel("bloque")
        a1.set_ylabel("borde [mV*ns]")
        a1.set_xlim(0.5, c["n"] + 0.5)
        a1.grid(alpha=0.3)
        a1.legend(fontsize=8)
        a1.set_title("%s: las %d 'mediciones', con la MISMA ventana"
                     % (e["clave"], c["n"]), fontsize=10)

        for k, (ch, y) in enumerate(d["espectros"]):
            a2.step(ch, y, where="mid", lw=0.6, alpha=0.55,
                    color=matplotlib.colormaps["viridis"](
                        k / max(len(d["espectros"]) - 1, 1)))
        a2.axvspan(d["r"]["intervalo"][0], d["r"]["intervalo"][1],
                   color="tab:orange", alpha=0.10, label="ventana comun")
        a2.set_xlabel("canal [mV*ns]")
        a2.set_ylabel("cuentas netas del bloque")
        a2.grid(alpha=0.3)
        a2.legend(fontsize=8)

        npb = e["n_pulsos"] / float(c["n"])
        L = ["corrida        %s  (%s)" % (e["clave"], e["punto"]),
             "bloques        %d, de ~%.0f pulsos cada uno" % (c["n"], npb),
             "ventana comun  [%.0f, %.0f]" % d["r"]["intervalo"],
             "",
             " metodo        completa    media    SD (%)      SD/sqrt(n) (%)"]
        for m in correr.METODOS:
            s = d["disp"][m]
            b0 = abs(d["bordes"][m]) if np.isfinite(d["bordes"][m]) else np.nan
            L.append(" %-12s %9.1f %8.1f  %7.2f (%5.2f%%) %7.2f (%5.2f%%)"
                     % (m, d["bordes"][m], s["media"], s["sd"],
                        100 * s["sd"] / b0, s["sd"] / np.sqrt(c["n"]),
                        100 * s["sd"] / np.sqrt(c["n"]) / b0))
        L += ["",
              "bloques utiles: " + ", ".join(
                  "%s %d/%d" % (m, d["disp"][m]["n"], c["n"])
                  for m in correr.METODOS),
              "",
              "La columna SD/sqrt(%d) escala la dispersion a la" % c["n"],
              "estadistica de una corrida completa, que es lo que",
              "tiene cada una de las diez mediciones del paper.",
              "Sigue sin ser su barra: la de ellos incluye montaje,",
              "posicion de la fuente y deriva de ganancia."]
        self.res("\n".join(L))


# ==========================================================================
#  10 - auditoria de los numeros publicados
# ==========================================================================
class TabAuditoria(Pestana):
    titulo = "10 · Auditoría del paper"
    auto = True
    ancho_res = 40
    centro = "ambos"

    paper = [("Seccion III, sobre la Fig. 3:",
              "the standard deviation associated with the NCC-SA calibration "
              "consistently remains below that of the full energy peak fitting "
              "method, indicating higher precision"),
             ("Resumen y conclusion, sobre la Fig. 4:",
              "reduced mean errors by 80% and 85% (BC408) and 44% and 66% "
              "(EJ309) compared to the difference method and the Gaussian "
              "method")]
    codigo = ("replica_paper/auditoria.py -> fig3(), rotulos(), ec1(),\n"
              "        ec12(), coseno(), fig5(), fig4(), inconsistencias()\n"
              "replica_paper/resultados.py -> ajuste_lineal(), similitudes()")
    abiertas = [
        ("que es 'SD' en el recuadro de la Fig. 3",
         "El paper no lo escribe. Se dedujo reproduciendo la figura desde sus "
         "propios puntos publicados: es el RMS de los residuos con ddof = 0. "
         "La prueba es que el ajuste devuelve las rectas impresas a 4-5 "
         "cifras."),
        ("entre que vectores es la similitud coseno",
         "El paper reporta 99.9993% y 99.9999% y nunca lo dice. Se prueban "
         "cinco lecturas razonables y ninguna reproduce los DOS valores."),
    ]
    mirar = """
ACA NO HAY NINGUNA MEDICION NUESTRA. Todo sale de los numeros que el paper
publica: los puntos rotulados uno por uno en su Fig. 3, los seis paneles de su
Fig. 5 y las frases de su texto. No hubo que digitalizar pixeles: los puntos
estan rotulados, y la prueba de que se leyeron bien es que el ajuste devuelve
las rectas impresas a 4-5 cifras.

LO QUE SALIO:

1. SD y R2 estan INTERCAMBIADOS entre azul y rojo en el recuadro, en los dos
   paneles. Con la asignacion correcta, la frase de la Seccion III dice lo
   contrario de lo que muestran sus propios datos.
2. Siete rotulos cruzados en el panel de NaI (Na22 donde va Cs137, etc.). No
   cambian ningun ajuste, pero atribuyen lineas a la fuente equivocada.
3. La Ec. 1 impresa esta mal y sus propios rotulos siguen la corregida.
4. La Ec. 12 impresa esta mal, y se mide en sus propios canales: 74 medidos
   contra 76.0 de la impresa y 73.8 de la corregida.
5. La similitud coseno no es verificable.
6. La robustez de la Fig. 5 tiene numero: 1.24%, o 0.33% sin el intervalo mas
   angosto.
7. La Fig. 4 se deja SIN auditar a proposito: el unico soporte de sus 80/85% y
   44/66% es el grafico, y leerlo tendria una precision del mismo orden que
   los valores que se quieren leer. Digitalizarlo daria una cifra inventada.
"""

    BLOQUES = [("todo", None), ("Fig. 3: recta, SD y R2", "fig3"),
               ("Fig. 3: rotulos cruzados", "rotulos"),
               ("Ec. 1", "ec1"), ("Ec. 12", "ec12"),
               ("similitud coseno", "coseno"), ("Fig. 5: robustez", "fig5"),
               ("Fig. 4: por que no es auditable", "fig4"),
               ("otras inconsistencias", "inconsistencias")]

    def construir_controles(self, f):
        self.v_b = self.combo(f, "bloque", [t for t, _ in self.BLOQUES],
                              ancho=26)
        self.nota(f, "El grafico reproduce la Fig. 3 desde los puntos que el "
                     "propio paper rotula (auditoria.FIG3). Las rectas de "
                     "abajo son el ajuste de esos puntos.")

    def leer_controles(self):
        d = dict(self.BLOQUES)
        return dict(fn=d[self.v_b.get()])

    def trabajo(self):
        fn = self.cfg["fn"]
        txt = texto_de(auditoria.todo if fn is None else getattr(auditoria, fn))
        ajustes = {}
        for det, d in auditoria.FIG3.items():
            ajustes[det] = dict(
                fotopico=auditoria._linea("fotopico", d["fotopico"]),
                borde=auditoria._linea("borde", d["borde"]),
                impreso=d["impreso"])
        return dict(txt=txt, ajustes=ajustes)

    def pintar(self, d):
        dets = list(d["ajustes"])
        for k, det in enumerate(dets):
            a = self.fig.add_subplot(1, len(dets), k + 1)
            g = d["ajustes"][det]
            for serie, col, nom in (("fotopico", "tab:blue", "fotopico"),
                                    ("borde", "tab:red", "borde Compton")):
                r = g[serie]
                a.plot(r["canal"], r["energia"], "o", color=col, ms=6,
                       label="%s  (SD %.3f)" % (nom, r["sd"]))
                x = np.linspace(0, r["canal"].max() * 1.08, 30)
                a.plot(x, r["m"] * x + r["b"], "-", color=col, lw=1.1)
            imp = g["impreso"]
            a.text(0.03, 0.97,
                   "recuadro del paper:\n SD azul %.6f\n SD rojo %.6f\n\n"
                   "lo que dan sus puntos:\n azul %.6f\n rojo %.6f"
                   % (imp["sd_azul"], imp["sd_rojo"], g["fotopico"]["sd"],
                      g["borde"]["sd"]),
                   transform=a.transAxes, va="top", fontsize=7.5,
                   family="monospace", bbox=dict(fc="w", ec="0.7", alpha=0.9))
            a.set_title("Fig. 3, %s" % det, fontsize=10)
            a.set_xlabel("canal")
            a.set_ylabel("energia [keV]")
            a.grid(alpha=0.3)
            a.legend(fontsize=7.5, loc="lower right")

        self.salida(d["txt"])
        L = ["Los cuatro ajustes, contra lo impreso:", ""]
        for det, g in d["ajustes"].items():
            imp = g["impreso"]
            L += [det,
                  "  fotopico  y = %.4f x %+.3f" % (g["fotopico"]["m"],
                                                    g["fotopico"]["b"]),
                  "            impreso  %.2f x %+.2f" % imp["fotopico"],
                  "  borde     y = %.4f x %+.3f" % (g["borde"]["m"],
                                                    g["borde"]["b"]),
                  "            impreso  %.2f x %+.2f" % imp["borde"],
                  "  SD azul   %.6f  (recuadro %.6f)"
                  % (g["fotopico"]["sd"], imp["sd_azul"]),
                  "  SD rojo   %.6f  (recuadro %.6f)"
                  % (g["borde"]["sd"], imp["sd_rojo"]),
                  ""]
        L += ["El informe completo esta en el panel de abajo."]
        self.res("\n".join(L))


# ==========================================================================
#  11 - Seccion III: la recta canal-energia, y lo que no se puede replicar
# ==========================================================================
class TabSeccionIII(Pestana):
    titulo = "11 · Sec. III · La recta"
    necesita_datos = True
    auto = False
    ancho_res = 44
    centro = "ambos"

    paper = [("Fig. 2, ultimo bloque antes de Stop:",
              "Fit Energy-Band Structure Relationship"),
             ("Seccion III: la recta se reporta con pendiente, ordenada, SD y "
              "R^2, y ademas con", "cosine similarity of the energy-channel "
              "relationship")]
    codigo = ("replica_paper/resultados.py -> ajuste_lineal(),\n"
              "        rectas_por_punto(), similitudes(), errores_de_metodo()\n"
              "replica_paper/correr.py -> puntos_de_calibracion(),\n"
              "        seccion_iii(), fig4(), no_replicable()")
    abiertas = [
        ("con que se agrupan los puntos de una recta",
         "El paper usa un detector y una ganancia por panel y no lo discute. "
         "Aca la recta se ajusta SIEMPRE dentro de un mismo punto de operacion "
         "(detector + alto voltaje): mezclar voltajes no da una recta "
         "canal-energia, da la ley de ganancia del fotodetector, que es otra "
         "medicion."),
        ("a que borde se le asigna el canal cuando la fuente tiene dos lineas",
         "Al de MENOR energia, que es la componente que el ajuste conserva. "
         "Misma regla para los tres metodos."),
    ]
    mirar = """
ESTA ES LA PARTE QUE NO SE PUEDE CERRAR CON LO MEDIDO, y el motivo es de
medicion, no de analisis.

El paper arma cada recta con CUATRO fuentes medidas con el mismo detector y la
misma ganancia: 5 o 6 bordes por recta. Nosotros tenemos una fuente por punto
de operacion (Co-60 a 1550 V, Cs-137 a 1600 V). Con un borde no hay recta; con
dos, la recta pasa exacta y SD = 0, R2 = 1 no miden nada.

Y hay un detalle estructural que conviene no perder: la Ec. 8 mapea canal a
energia POR EL ORIGEN, asi que todos los bordes que salen de UNA MISMA FUENTE
caen exactos sobre una recta por el origen. La recta solo mide algo si los
puntos vienen de fuentes DISTINTAS.

Sobre la similitud coseno, una advertencia que el cierre (pestana 12) vuelve
concreta: entre dos calibraciones cuyas pendientes difieren un 3.0%, la
similitud coseno de las energias reconstruidas da 0.99998812. La metrica es
insensible por construccion. Hay que reportar residuos en keV y la ordenada con
su error.

El panel de abajo tiene el informe completo, incluida la lista de lo que no se
puede replicar con el motivo de cada punto.
"""

    def construir_controles(self, f):
        self.v_cuales = self.combo(f, "corridas", ["las que tengan fondo",
                                                   "solo la elegida"], ancho=22)
        self.v_diag = self.check(f, "diagnostico del solapamiento", True)
        self.nota(f, "Corre la cadena NCC-SA completa en cada corrida. Con las "
                     "dos y el diagnostico tildado tarda ~1 min.")

    def leer_controles(self):
        cuales = (dat.con_fondo() if self.v_cuales.get().startswith("las")
                  else [self.app.estado.clave])
        return dict(claves=cuales, diag=bool(self.v_diag.get()))

    def trabajo(self):
        c = self.cfg
        analisis = []
        n = len(c["claves"])
        for k, clave in enumerate(c["claves"]):
            self.fase(k / float(n), (k + 1) / float(n))
            self.progreso("%s: espectro..." % clave, 0.0)
            e = dat.espectro(clave, bins=dat.BINNING.get(clave, {}).get("bins"),
                             rango_max=dat.BINNING.get(clave, {}).get("rango_max"))
            self.progreso("%s: NCC -> SA..." % clave, 0.35)
            bordes, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"])
            diag = None
            if c["diag"] and len(e["lineas"]) == 2:
                self.progreso("%s: diagnostico del solapamiento..." % clave, 0.6)
                diag = correr.diagnostico_dos_bordes(e, r, verbose=False)
            analisis.append(dict(clave=clave, espectro=e, ajuste=r,
                                 bordes=bordes, dispersion=None, robustez={},
                                 punto=e["punto"], lineas=e["lineas"],
                                 diagnostico=diag))
        self.fase(0.0, 1.0)
        self.progreso("armando la Seccion III...", 0.95)
        puntos = correr.puntos_de_calibracion(analisis)
        rectas = correr.seccion_iii(puntos, verbose=False)
        txt = (texto_de(correr.seccion_iii, puntos)
               + texto_de(correr.fig4, analisis)
               + texto_de(correr.no_replicable, analisis))
        return dict(analisis=analisis, puntos=puntos, rectas=rectas, txt=txt)

    def pintar(self, d):
        a = self.fig.add_subplot(111)
        marcas = {"NCC-SA": "o", "Gauss FWHM": "s", "Diferencia": "^"}
        for p in d["puntos"]:
            if not np.isfinite(p["canal"]):
                continue
            a.plot(p["canal"], p["energia"], marcas[p["metodo"]],
                   color=COL[p["metodo"]], ms=8, mfc="none", mew=1.6)
            a.annotate("%s\n%s" % (p["punto"], p["fuente"]),
                       (p["canal"], p["energia"]), textcoords="offset points",
                       xytext=(8, -4), fontsize=7.5, color="0.35")
        for m in correr.METODOS:
            a.plot([], [], marcas[m], color=COL[m], mfc="none", mew=1.6,
                   label=m)
        for (punto, met), r in sorted(d["rectas"].items()):
            x = np.linspace(0, float(np.max(r["canal"])) * 1.15, 40)
            a.plot(x, r["m"] * x + r["b"], "--", lw=1.0, color=COL[met])
        if not d["rectas"]:
            a.text(0.03, 0.96,
                   "Ningun punto de operacion tiene dos o mas bordes medidos:\n"
                   "no hay recta que ajustar. Es la Fig. 3 del paper, y no se\n"
                   "puede hacer con lo medido. El motivo esta abajo.",
                   transform=a.transAxes, va="top", fontsize=9,
                   bbox=dict(fc="#fff6e5", ec="#c8a45a"))
        a.set_xlabel("canal del borde [mV*ns]")
        a.set_ylabel("energia del borde [keV]")
        a.set_title("Fig. 3 del paper, con nuestras corridas", fontsize=10)
        a.grid(alpha=0.3)
        a.legend(fontsize=8, loc="lower right")

        self.salida(d["txt"])
        L = ["puntos de calibracion obtenidos:", "",
             " punto de operacion  metodo       canal    E [keV]"]
        for p in d["puntos"]:
            L.append(" %-18s  %-11s %8.1f  %7.1f"
                     % (p["punto"], p["metodo"], p["canal"], p["energia"]))
        L += ["", "rectas ajustadas: %d" % len(d["rectas"]),
              "", "El informe completo esta en el panel de abajo."]
        self.res("\n".join(L))


# ==========================================================================
#  12 - cierre con verdad conocida
# ==========================================================================
class TabCierre(Pestana):
    titulo = "12 · Cierre con verdad conocida"
    auto = False
    ancho_res = 44
    centro = "ambos"

    paper = [("Lo que el paper afirma de la Fig. 4, y que este cierre pone a "
              "prueba:",
              "reduced mean errors by 80% and 85% (BC408) and 44% and 66% "
              "(EJ309) compared to the difference method and the Gaussian "
              "method")]
    codigo = "replica_paper/correr.py -> cierre()"
    abiertas = [
        ("D, la resolucion y la semilla del ruido",
         "Son de este cierre, no del paper: definen la verdad contra la que se "
         "mide. Moverlas es ver a partir de que resolucion el metodo deja de "
         "funcionar."),
    ]
    mirar = """
ESTE ES EL UNICO LUGAR DONDE SE PUEDE MEDIR SI EL METODO ACIERTA, porque es el
unico donde se sabe la respuesta.

Genera dos fuentes -Cs-137 y Co-60- al MISMO punto de operacion, con la misma D
conocida, que es exactamente la situacion que hoy no esta medida. Despues corre
la Seccion III entera: la recta, el SD, el R2, la similitud coseno y los
errores relativos.

LO QUE SALE: la afirmacion central del paper se reproduce. El NCC-SA recupera
la pendiente con un error del orden de 0.1%, mientras el metodo de la
diferencia y el gaussiano se van varios por ciento. Y sobre el Co-60, que tiene
los dos bordes solapados, el metodo gaussiano se corre mucho mas que el NCC-SA.

Tambien queda a la vista la insensibilidad de la similitud coseno: entre
calibraciones que difieren varios por ciento en pendiente sigue dando 0.9999.

Ojo con una cosa: las "diez mediciones" de este cierre son diez realizaciones
de Poisson del MISMO espectro, o sea solo estadistica de conteo. No incluyen
montaje ni deriva. Lo mismo vale para la pestana 9.
"""

    def construir_controles(self, f):
        self.v_D = self.entrada(f, "D verdadera", "1.30")
        self.v_sr = self.entrada(f, "sigma/B", "0.07")
        self.v_sem = self.entrada(f, "semilla", "2026")
        self.nota(f, "Tarda ~1 min: son dos fuentes x (1 + 10) cadenas "
                     "NCC-SA, mas los dos metodos convencionales en cada una.")

    def leer_controles(self):
        return dict(D=float(self.v_D.get()), sr=float(self.v_sr.get()),
                    sem=int(self.v_sem.get()))

    def trabajo(self):
        c = self.cfg
        self.progreso("generando las dos fuentes y corriendo la Seccion III...")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            analisis, rectas = correr.cierre(D_real=c["D"], sigma_rel=c["sr"],
                                             semilla=c["sem"])
        return dict(analisis=analisis, rectas=rectas, txt=buf.getvalue(),
                    cfg=c)

    def pintar(self, d):
        c = d["cfg"]
        n = len(d["analisis"])
        for k, a in enumerate(d["analisis"]):
            ax = self.fig.add_subplot(2, n, k + 1)
            e = a["espectro"]
            ax.step(e["canales"], e["cuentas"], where="mid", lw=0.8,
                    color="0.6", label="sintetico + Poisson")
            ax.plot(e["canales"], e["limpio"], color="k", lw=1.0,
                    label="sin ruido")
            ax.axvline(a["canal_verdadero"], color="k", ls="-", lw=1.6,
                       label="verdad: %.1f" % a["canal_verdadero"])
            for m in correr.METODOS:
                if np.isfinite(a["bordes"][m]):
                    ax.axvline(a["bordes"][m], color=COL[m], ls=LSS[m], lw=1.3,
                               label="%s %+.2f%%"
                               % (m, 100 * (a["bordes"][m]
                                            / a["canal_verdadero"] - 1)))
            ax.set_title(a["clave"], fontsize=10)
            ax.set_xlabel("canal")
            ax.set_ylabel("cuentas")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7)

        ax = self.fig.add_subplot(2, 1, 2)
        anchura = 0.26
        x = np.arange(n)
        for i, m in enumerate(correr.METODOS):
            v = [100 * (a["bordes"][m] / a["canal_verdadero"] - 1)
                 for a in d["analisis"]]
            ax.bar(x + (i - 1) * anchura, v, anchura, color=COL[m], label=m)
        ax.axhline(0, color="k", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels([a["clave"] for a in d["analisis"]])
        ax.set_ylabel("error del borde [%]")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)
        ax.set_title("error contra la verdad (D = %.4f, sigma/B = %.3f)"
                     % (c["D"], c["sr"]), fontsize=10)

        self.salida(d["txt"])
        L = ["VERDAD",
             "  D          %.4f mV*ns/keV" % c["D"],
             "  sigma/B    %.3f" % c["sr"],
             "  m = 1/D    %.4f keV/canal" % (1.0 / c["D"]),
             "  ordenada   0 (recta por el origen)",
             "",
             "ERROR DEL BORDE, por fuente:", ""]
        for a in d["analisis"]:
            L.append("  %s  (verdad %.1f)" % (a["clave"], a["canal_verdadero"]))
            for m in correr.METODOS:
                b = a["bordes"][m]
                L.append("     %-12s %8.1f  %+7.2f%%   SD10 %6.1f"
                         % (m, b, 100 * (b / a["canal_verdadero"] - 1),
                            a["dispersion"][m]["sd"]))
            L.append("")
        L += ["PENDIENTE RECUPERADA", "",
              "  metodo        m ajustada   error contra la verdad"]
        for (p, met), r in sorted(d["rectas"].items()):
            L.append("  %-12s %10.4f    %+7.2f%%"
                     % (met, r["m"], 100 * (r["m"] * c["D"] - 1)))
        L += ["", "El informe completo esta en el panel de abajo."]
        self.res("\n".join(L))


# ==========================================================================
#  Mapa: el diagrama de la Fig. 2 y donde vive cada ecuacion
# ==========================================================================
class TabMapa(Pestana):
    titulo = "· Mapa del paper"
    auto = True
    ancho_res = 40
    centro = "ambos"

    # (numero, que es, donde esta implementada, que pestana la muestra)
    ECUACIONES = [
        ("1", "borde Compton E_max(E_gamma)", "ecuaciones.e_max", "1"),
        ("2", "Klein-Nishina en angulo solido", "ecuaciones.dsigma_domega", "2"),
        ("3", "sigma_total, suma en cuadratura", "ecuaciones.sigma_total", "3"),
        ("Tabla I", "los tres terminos de sigma", "(escrita, sin evaluar)", "3"),
        ("8", "lineshape -> canal + convolucion", "ecuaciones.modelo", "3"),
        ("9", "el nucleo gaussiano", "ecuaciones.kernel_gaussiano", "3"),
        ("10", "Klein-Nishina en energia del e-", "ecuaciones.dsigma_de", "2"),
        ("11", "Theta = {A1,A2,B1,s1,s2,D}", "metodo.ajustar_dos_bordes", "6"),
        ("12", "CE_ch,2 = B1 + D*dE", "ecuaciones.delta_e_bordes", "6"),
        ("13", "suma de las dos componentes", "ecuaciones.modelo_dos_bordes", "6"),
        ("14", "loss de minimos cuadrados", "sa.costo_ec14", "5"),
        ("15", "Pearson entre S1 y S2", "ecuaciones.pearson", "4"),
    ]

    ERRATAS = [
        ("Ec. 1", "da la energia del foton retrodispersado, no la del "
                  "electron. Falta el 2*alpha en el numerador.",
         "sus propios rotulos (477, 639, 963, 1061, 1117 keV) siguen la "
         "corregida", "1 y 10"),
        ("Ec. 12", "define dE como la diferencia entre las LINEAS GAMMA; D "
                   "convierte energia del electron, asi que va la diferencia "
                   "de los BORDES.",
         "sus dos bordes del Co-60 estan a 74 canales; la impresa predice 76.0 "
         "y la corregida 73.8", "6 y 10"),
        ("Ec. 14 vs Metropolis",
         "la Ec. 14 define la loss como minimos cuadrados y el parrafo de "
         "Metropolis dice que E(x) es el 'NCC-based matching error'. Son "
         "objetivos distintos y hay que elegir.",
         "la Ec. 14 es una suma de cuadrados sobre cuentas y el error NCC esta "
         "acotado en [0,2], asi que T_init = 2000 no significa lo mismo en los "
         "dos (la pestana 5 lo mide)", "5"),
        ("Fig. 3", "SD y R2 estan intercambiados entre azul y rojo, en los dos "
                   "paneles.",
         "el ajuste de sus propios puntos reproduce las rectas impresas a 4-5 "
         "cifras y devuelve los SD cruzados", "10"),
        ("Fig. 3", "siete rotulos de fuente cruzados en el panel de NaI.",
         "las energias estan bien; los nombres de fuente no", "10"),
    ]

    FLUJO = [("Energy Spectrum\n(Channel, Counts)", "0"),
             ("Implement NCC Automatic\nInterval Selection", "4"),
             ("Calculate Pearson\nCorrelation Coefficient", "4"),
             ("Maximum Data\nInterval of rho?", "4"),
             ("Set Theta =\n{A1,A2,B1,s1,s2,D}", "6"),
             ("Use SA to Solve Local\nOptimal Problems", "5"),
             ("Optimal Parameter\nSolution?", "6"),
             ("Fit Energy-Band\nStructure Relationship", "11")]

    codigo = ("todos los modulos de replica_paper/:\n"
              "  ecuaciones · ncc · sa · convencionales · metodo\n"
              "  datos · resultados · auditoria · correr")
    mirar = """
ESTA PESTANA ES EL INDICE. El diagrama es la Fig. 2 del paper tal como esta
impresa, y al lado de cada bloque dice que pestana lo calcula.

La regla de esta carpeta: aca solo entra lo que el paper dice. Plantilla
Klein-Nishina analitica (Ec. 2 y 10), no Monte Carlo. Nada de modelo extendido,
mezcla de scattering, termino de exceso a bajo canal, chi2 de Pearson iterado,
MCMC ni presupuesto de sistematicos medido: todo eso vive en codigo/ y
analisis/ y es EXTRA respecto del paper.

Lo unico prestado de afuera es el paso archivo -> (canal, cuentas), que es
justamente donde ARRANCA el diagrama. Ese prestamo esta aislado en datos.py.

En el panel de abajo: cada ecuacion del paper con la funcion que la implementa
y la pestana que la muestra, y las erratas encontradas con la evidencia de cada
una. Ninguna errata se afirma sin un numero que la sostenga, y en los cinco
casos el numero sale de los PROPIOS datos publicados del paper.
"""

    def construir_controles(self, f):
        self.nota(f, "Nada que configurar: esta pestana es el indice. Los "
                     "numeros de la derecha de cada bloque son las pestanas "
                     "de arriba.")

    def trabajo(self):
        # Devuelve un dict y no None para que `repintar()` pueda redibujar al
        # cambiar la nitidez: el criterio de ahi es "hay un ultimo resultado".
        return dict(corridas=dat.disponibles(), con_fondo=dat.con_fondo())

    def pintar(self, d):
        # Todo en coordenadas de EJE, no de datos: asi el diagrama se acomoda
        # solo cuando cambia el tamano de la ventana o la nitidez, en vez de
        # dejar los rotulos a media pantalla de las cajas.
        a = self.fig.add_subplot(111)
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        a.set_axis_off()
        T = a.transAxes
        n = len(self.FLUJO)
        y0, y1 = 0.855, 0.115
        ds = (y0 - y1) / (n - 1)
        flecha = dict(arrowstyle="-|>", color="#1f4e79", lw=1.2)

        def caja(y, txt, fc, ec, fs, pad):
            a.text(0.42, y, txt, ha="center", va="center", fontsize=fs,
                   transform=T,
                   bbox=dict(boxstyle="round,pad=%.2f" % pad, fc=fc, ec=ec,
                             lw=1.2))

        for i, (txt, pest) in enumerate(self.FLUJO):
            y = y0 - i * ds
            rombo = txt.strip().endswith("?")
            caja(y, txt, "#fdf3e3" if rombo else "#eaf1f8",
                 "#c8a45a" if rombo else "#1f4e79", 8.5, 0.5 if rombo else 0.45)
            a.text(0.99, y, "pestana %s" % pest, ha="right", va="center",
                   fontsize=8, color="#a4501a", family="monospace",
                   transform=T)
            if i < n - 1:
                a.annotate("", xy=(0.42, y - ds + 0.030),
                           xytext=(0.42, y - 0.030), xycoords=T, textcoords=T,
                           arrowprops=flecha)
        caja(y0 + ds * 0.75, "Start", "#dfe8ef", "0.5", 9, 0.4)
        a.annotate("", xy=(0.42, y0 + 0.028),
                   xytext=(0.42, y0 + ds * 0.75 - 0.022), xycoords=T,
                   textcoords=T, arrowprops=flecha)
        caja(y1 - ds * 0.75, "Stop", "#dfe8ef", "0.5", 9, 0.4)
        a.annotate("", xy=(0.42, y1 - ds * 0.75 + 0.022),
                   xytext=(0.42, y1 - 0.028), xycoords=T, textcoords=T,
                   arrowprops=flecha)
        # El titulo va a la izquierda y no arriba del centro: arriba del
        # centro pisa la caja de Start cuando la figura es baja.
        a.text(0.01, 0.97, "Fig. 2\ndel paper", fontsize=11, va="top",
               ha="left", color="#1f4e79", weight="bold", transform=T)

        L = ["LAS ECUACIONES DEL PAPER", "",
             " ec.      que es                            implementada en"
             "                    pestana"]
        for num, que, donde, pest in self.ECUACIONES:
            L.append(" %-8s %-33s %-34s %s" % (num, que, donde, pest))
        L += ["",
              " Las Ec. 4 a 7 del paper son el desarrollo hacia la Ec. 8 y no "
              "definen nada que",
              " haya que calcular aparte.",
              "",
              "", "ERRATAS Y CONTRADICCIONES ENCONTRADAS", ""]
        for k, (donde, que, evidencia, pest) in enumerate(self.ERRATAS, 1):
            L += [" %d. %s: %s" % (k, donde, que),
                  "    evidencia: %s" % evidencia,
                  "    se ve en la pestana %s" % pest, ""]
        ambiguas = [
            ("con que sigma se construye S2", "ncc.plantilla",
             "se barre una grilla; gana el de mayor rho"),
            ("donde empieza el plateau y termina el borde",
             "ncc.INI_REL / FIN_REL",
             "[%.2f E_max, (1 + 4.5 sigma_rel) E_max]" % mod_ncc.INI_REL),
            ("largo minimo de un intervalo candidato", "ncc.MIN_BINS",
             "%d bins" % mod_ncc.MIN_BINS),
            ("cotas del SA (los d_i)", "metodo.COTAS_REL",
             "derivadas de la semilla del NCC"),
            ("de donde sale x0", "metodo.semilla_por_barrido",
             "barrido grueso de (B1, D, sigma1) con nnls"),
            ("el lazo 'Optimal Parameter Solution?'", "metodo._cerrar",
             "pulido local desde cada salida del SA"),
            ("con que se suaviza (metodo de la diferencia)",
             "convencionales.VENTANA",
             "Savitzky-Golay de %d bins, orden %d" % (conv.VENTANA, conv.ORDEN)),
            ("que region usan los convencionales", "correr.tres_metodos",
             "la misma que elige el NCC"),
            ("que es 'the maximum amplitude'", "convencionales.metodo_gauss",
             "la amplitud de la gaussiana ajustada"),
        ]
        L += ["", "LO QUE EL PAPER NO ESPECIFICA Y HUBO QUE DECIDIR", "",
              " %-44s %-28s %s" % ("ambiguedad", "donde", "que se hizo")]
        for que, donde, como in ambiguas:
            L.append(" %-44s %-28s %s" % (que, donde, como))
        self.salida("\n".join(L))
        self.res("\n".join([
            "Kong, Zhang, Lin, Qiu, Chen, Wang.",
            "The Automatic Calibration Method of the",
            "Compton Edge Based on Normalized",
            "Cross-correlation and Simulated Annealing",
            "Algorithm.",
            "IEEE Trans. Nucl. Sci., 2025.",
            "arXiv:2507.16282 [physics.ins-det]",
            "",
            "El PDF esta en docs/Paper ref.pdf",
            "",
            "-" * 38,
            "",
            "ecuaciones implementadas   %d + la Tabla I"
            % (len(self.ECUACIONES) - 1),
            "  (de las 15 del paper: las",
            "   Ec. 4 a 7 son el desarrollo",
            "   hacia la Ec. 8)",
            "erratas con evidencia      %d" % len(self.ERRATAS),
            "",
            "corridas accesibles ahora  %d" % len(d["corridas"]),
            "con fondo propio           %d" % len(d["con_fondo"]),
            "",
            "El detalle esta en el panel de abajo.",
        ]))


# ==========================================================================
#  App
# ==========================================================================
PESTANAS = [TabMapa, TabDatos, TabEc1, TabEc10, TabForward, TabNCC, TabSA,
            TabCadena, TabConvencionales, TabRobustez, TabBloques,
            TabAuditoria, TabSeccionIII, TabCierre]


class App(tk.Tk):
    def __init__(self):
        activar_dpi_awareness()              # antes de crear la ventana
        super().__init__()
        self.estado = Estado()
        self.cola = queue.Queue()
        self.estado_txt = tk.StringVar(value="listo")
        self.pct = tk.StringVar(value="")
        self._trabajos = set()
        self._t0 = 0.0
        self._modo = "indeterminado"
        self.configurar_nitidez()
        self.actualizar_titulo()

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True)
        self.pestanas = []
        for clase in PESTANAS:
            p = clase(self.nb, self)
            self.nb.add(p, text=clase.titulo)
            self.pestanas.append(p)

        barra = ttk.Frame(self)
        barra.pack(side="bottom", fill="x")
        ttk.Label(barra, textvariable=self.estado_txt, anchor="w",
                  relief="sunken").pack(side="left", fill="x", expand=True)
        self.barra_prog = ttk.Progressbar(barra, mode="determinate",
                                          maximum=100.0,
                                          length=int(180 * self.escala))
        self.barra_prog.pack(side="left", padx=(8, 4))
        ttk.Label(barra, textvariable=self.pct, width=20,
                  anchor="w").pack(side="left")
        ttk.Label(barra, text="nitidez").pack(side="left", padx=(10, 3))
        self.var_nit = tk.StringVar(value="%.2f" % self.escala)
        cb = ttk.Combobox(barra, textvariable=self.var_nit, width=6,
                          state="readonly",
                          values=["1.00", "1.25", "1.50", "1.75", "2.00", "2.50"])
        cb.pack(side="left", padx=(0, 8))
        cb.bind("<<ComboboxSelected>>", self.cambiar_nitidez)

        self.nb.bind("<<NotebookTabChanged>>", self.al_cambiar)
        self.after(100, self.revisar_cola)
        self.after(200, lambda: self.al_cambiar(None))

    # ------------------------------------------------------------------
    def actualizar_titulo(self):
        self.title("Kong et al. 2025 paso a paso  -  %s" %
                   self.estado.descripcion())

    def propagar_corrida(self, excepto=None):
        """
        Marca como pendientes las pestanas que dependen de los datos.

        No se recalculan solas: varias tardan minutos y recalcular las diez de
        golpe seria peor que dejar que se abra la que interesa. Lo que si hay
        que hacer es BORRAR el grafico viejo, para que no quede en pantalla un
        resultado de otra corrida.
        """
        for p in self.pestanas:
            if p is excepto or p._corriendo:
                continue
            p.al_cambiar_corrida()
        self.estado_txt.set("corrida cambiada: cada pestana se recalcula al "
                            "abrirla o con «Recalcular»")

    # ------------------------------------------------------------------
    def progreso_iniciar(self, pest):
        """
        Arranca la barra cuando una pestana lanza su calculo.

        Empieza INDETERMINADA (la barrita que va y viene) porque la mayoria de
        los pasos no sabe de antemano cuanto falta, y para lo que se le pide
        -saber que el programa no se colgo- eso ya alcanza. En cuanto llega el
        primer progreso(..., frac) se pasa sola a determinada, con porcentaje y
        tiempo restante.
        """
        self._trabajos.add(pest)
        self._t0 = time.time()
        self._modo = "indeterminado"
        self.barra_prog.configure(mode="indeterminate")
        self.barra_prog.start(50)
        self.pct.set("trabajando...")

    def progreso_actualizar(self, frac):
        if not self._trabajos or frac is None:
            return
        frac = max(0.0, min(float(frac), 1.0))
        if self._modo != "determinado":
            self.barra_prog.stop()
            self.barra_prog.configure(mode="determinate", maximum=100.0)
            self._modo = "determinado"
        self.barra_prog["value"] = 100.0 * frac
        t = time.time() - self._t0
        # el ETA sale de una regla de tres: antes de tener algo de historia
        # (3% y 2 s) es puro ruido y conviene no mostrarlo.
        if frac > 0.03 and t > 2.0:
            self.pct.set("%.0f%% - faltan ~%s"
                         % (100 * frac, formato_tiempo(t * (1 - frac) / frac)))
        else:
            self.pct.set("%.0f%%" % (100 * frac))

    def progreso_terminar(self, pest):
        self._trabajos.discard(pest)
        if self._trabajos:
            return
        self.barra_prog.stop()
        self.barra_prog.configure(mode="determinate", maximum=100.0)
        self.barra_prog["value"] = 0.0
        self._modo = "indeterminado"
        self.pct.set("")

    # ------------------------------------------------------------------
    def configurar_nitidez(self):
        """
        Ajusta Tk y matplotlib al DPI real de la pantalla.

        Con DPI awareness activo Windows ya no estira la ventana, pero entonces
        todo queda CHICO si no escalamos nosotros. Dos perillas: `tk scaling`
        (las fuentes estan en puntos) y figure.dpi de matplotlib.
        """
        try:
            dpi = float(self.winfo_fpixels("1i"))
        except Exception:
            dpi = 96.0
        self.escala = max(1.0, min(dpi / 96.0, 3.0))
        try:
            self.tk.call("tk", "scaling", dpi / 72.0)
        except Exception:
            pass
        self.dpi_fig = 100.0 * self.escala
        matplotlib.rcParams.update({
            "figure.dpi": self.dpi_fig, "savefig.dpi": 150, "font.size": 9,
            "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
            "lines.antialiased": True, "path.simplify": True,
            "path.simplify_threshold": 0.8, "agg.path.chunksize": 10000,
        })
        ancho = min(int(1560 * self.escala), self.winfo_screenwidth() - 40)
        alto = min(int(900 * self.escala), self.winfo_screenheight() - 80)
        self.geometry("%dx%d" % (ancho, alto))
        self.minsize(1180, 700)

    def cambiar_nitidez(self, _e=None):
        try:
            f = float(self.var_nit.get())
        except ValueError:
            return
        self.escala = f
        self.dpi_fig = 100.0 * f
        matplotlib.rcParams["figure.dpi"] = self.dpi_fig
        for p in self.pestanas:
            p.fig.set_dpi(self.dpi_fig)
            p.repintar()
        self.estado_txt.set("nitidez x%.2f (dpi de figura = %.0f)"
                            % (f, self.dpi_fig))

    # ------------------------------------------------------------------
    def al_cambiar(self, _e):
        p = self.pestanas[self.nb.index(self.nb.select())]
        if p.calculada or p._corriendo:
            return
        if p.auto:
            p.lanzar()
        else:
            self.estado_txt.set("%s: apreta «Recalcular» (este calculo tarda)"
                                % p.titulo)

    def revisar_cola(self):
        try:
            while True:
                pest, res = self.cola.get_nowait()
                if res[0] == "progreso":
                    self.estado_txt.set("%s: %s" % (pest.titulo, res[1]))
                    self.progreso_actualizar(res[2] if len(res) > 2 else None)
                    continue
                pest.recibir(res)
        except queue.Empty:
            pass
        self.after(100, self.revisar_cola)


def main():
    rutas.salida_utf8()
    activar_dpi_awareness()
    App().mainloop()


if __name__ == "__main__":
    main()
