"""
gui.py — interfaz de calibracion con simulaciones Monte Carlo.

    python analisis/gui.py

Que hace: elegir los datos (del catalogo o archivos sueltos), correr la cadena
completa mostrando en que etapa va, y dejar el resultado guardado en el registro
por detector.

Que modelo corre esta interfaz
------------------------------
Dos cosas la separan de los scripts de linea de comandos, y las dos estan
arriba como constantes (`COMPUERTA_DEFECTO` y `EXCESO`, con el porque al lado):

  * la compuerta de integracion arranca en 12 ns despues del onset, y es un
    campo editable de la ventana en vez de un valor por archivo;
  * el modelo es la plantilla convolucionada y NADA MAS. Sin la componente
    exponencial de bajo canal no hay `lambda` ni amplitud `C`, y el recocido
    busca en dos dimensiones.

`analisis/calibracion.py` y `analisis/calibrar_todo.py` siguen con el modelo
completo: si dos numeros no coinciden, es por esto.

Por que la barra de avance esta hecha asi
-----------------------------------------
El parseo del XML es ~40 s y todo el resto de la cadena junto ~3 s, o sea que el
95% del tiempo transcurre en un solo paso que, sin barra, no dice nada. Tres
decisiones que salen de eso:

  * la fraccion del parseo se mide en BYTES leidos, no en eventos: cuantos
    eventos tiene el archivo no se sabe hasta terminarlo (`drs_xml.iter_eventos`);
  * los pesos de la barra global se calculan ANTES de empezar, mirando el tamano
    de cada XML y si sus areas ya estan cacheadas (`pipeline.repartir`);
  * hay una fila por etapa con su tiempo. Una barra sola dice "va por el 30%";
    la lista dice ademas cual etapa esta corriendo y cuanto tardo cada una, que
    es lo que distingue "lento" de "colgado" y, si algo falla, donde fallo.

El calculo corre en un hilo aparte y se comunica por una cola. El hilo
trabajador NUNCA toca widgets ni matplotlib: manda datos y el hilo de Tk dibuja.

Como se muestran los graficos
-----------------------------
En el docstring de `Lienzo`: por que la figura se encuadra a mano antes de cada
dibujo (si no, TkAgg la dibuja un 25% mas grande que la imagen sobre la que la
vuelca y el resto se RECORTA sin decir nada) y como se agranda un panel solo.
`verificaciones/_verif_gui_encuadre.py` comprueba las dos cosas.
"""

import io
import os
import queue
import re
import sys
import textwrap
import threading
import time
import traceback
import contextlib

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# --------------------------------------------------------------------------
# Nitidez en Windows: ANTES de crear el root de Tk
# --------------------------------------------------------------------------
def activar_dpi_awareness():
    """
    Declara el proceso consciente del DPI y devuelve los px/pulgada reales.

    Sin esto Windows renderiza la ventana a 96 DPI y la ESTIRA como bitmap
    cuando la pantalla tiene escalado, y todo se ve borroso. Ninguna opcion de
    matplotlib lo arregla: hay que declararlo antes de que exista la ventana.
    Despues hay que escalar a mano, porque si no todo queda chico.
    """
    if sys.platform != "win32":
        return None
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(2)   # PER_MONITOR_AWARE
        return True
    except Exception:
        try:
            ctypes.windll.user32.SetProcessDPIAware()
            return True
        except Exception:
            return False


activar_dpi_awareness()

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

# Tipografia mas chica que la de fabrica (10 pt). Las pestanas son grillas de
# 2x2 dentro de un panel de ~900x600 px, o sea paneles de 450x300: con los
# tamanos por defecto los numeros de los ejes se tocan entre si y el titulo de
# la fila de abajo se mete en el eje de la de arriba. Va a nivel de MODULO y no
# dentro de la App para que las figuras salgan iguales cuando las dibuja
# `verificaciones/_verif_gui_pestanas.py`, que no crea ventana.
matplotlib.rcParams.update({"font.size": 9.0, "axes.titlesize": 10.0,
                            "axes.labelsize": 9.0, "xtick.labelsize": 8.0,
                            "ytick.labelsize": 8.0, "legend.fontsize": 8.0,
                            "figure.titlesize": 11.0,
                            "axes.titlepad": 4.0, "axes.labelpad": 2.0})


def ajustar(fig):
    """
    Que matplotlib recalcule el espaciado en CADA dibujado, no una sola vez.

    Aca habia `fig.tight_layout()`, que resuelve el espaciado con el tamano que
    la figura tiene EN ESE MOMENTO y no vuelve a intervenir. El lienzo de Tk
    cambia de tamano con la ventana y el redibujado por resize no lo llama de
    nuevo: en una grilla de 2x2 los titulos y los rotulos de los ejes terminan
    pisandose en cuanto la ventana no mide lo que media al calcular.

    El motor "constrained" resuelve lo mismo pero en cada draw, asi que el
    espaciado sigue al tamano real del panel. Ademas tiene en cuenta las barras
    de color y los ejes gemelos (`twinx`), que `tight_layout` maneja peor.

    ⚠ Hay que llamarla ANTES de crear los ejes, y por eso es la primera linea
    de cada figura en vez de la ultima. `fig.colorbar` se fija si la figura
    tiene motor de layout: si no lo tiene, mete la barra en un gridspec
    anidado con `height_ratios = [0, 1, 0]`, y encender el motor despues hace
    que el layout constrained divida por ese cero (ZeroDivisionError al
    dibujar, no al construir, o sea lejos de donde esta la causa).
    """
    fig.set_layout_engine("constrained")
    return fig


from analisis import (rutas, fisica, drs, pulso, simulacion, calibracion,
                      pipeline, voltaje,
                      registro, ncc as mod_ncc)

MANUAL = "(archivos sueltos)"
ESCALAS = {"luz  ->  keVee  (Lvis_keVee)": "Lvis_keVee",
           "energia depositada  ->  keV  (Edep_keV)": "Edep_keV"}

# --------------------------------------------------------------------------
# El modelo que corre esta interfaz
# --------------------------------------------------------------------------
# Compuerta de integracion [ns] DESPUES del onset de cada pulso. Es el valor
# con el que arranca el campo "Compuerta tras el onset"; se edita en la
# ventana antes de cada calibracion, y vaciarlo hace que cada archivo mida la
# suya como antes.
#
# 12 ns es corto a proposito: con pulsos cortos, seguir integrando despues de
# que el pulso termino no agrega senal, agrega la linea base y su ruido. Lo
# que se paga es que la compuerta no recoge la carga entera, asi que D deja de
# ser la ganancia de pulso completo: es la ganancia MEDIDA EN ESTA VENTANA.
# Para comparar dos corridas eso da igual mientras la ventana sea la misma
# (que es justamente para lo que el valor esta fijo y a la vista); para citar
# una D absoluta, no.
COMPUERTA_DEFECTO = 12.0

# Sin la componente exponencial de bajo canal: el modelo es la plantilla
# convolucionada y nada mas, S = A * [P(ch/D) (*) G(sigma)].
#
# Que se gana y que se pierde, para no descubrirlo despues:
#   - se gana un parametro menos que ajustar (se van lambda y la amplitud C),
#     y el recocido busca en dos dimensiones en vez de tres;
#   - se pierde la unica defensa contra el exceso de cuentas a bajo canal, que
#     es real y la plantilla no lo tiene. Con el exceso adentro el ajuste podia
#     arrancar en 0.25*B; sin el hay que dejarlo AFUERA de la ventana, y eso lo
#     decide `inicio_rel` (el campo de al lado). Ojo con bajarlo: sobre el
#     Co-60 el modelo desnudo da D = 1.2230 con inicio_rel = 0.45 y 1.2422 con
#     0.75, o sea que la ventana pasa a ser un sistematico del ~1.5% que se
#     elige a mano.
EXCESO = False


def tope_sobre_el_exceso(e, r):
    """
    Escala vertical del espectro, ignorando el pico de bajo canal.

    Antes el corte era `canales > 300` mV*ns, escrito a mano. Eso es una
    constante de CANAL, o sea que vale para un detector y una ganancia (§7.40,
    la misma trampa de siempre): con el detector D, que da D ~ 0.5 mV*ns/keV,
    el espectro entero vive debajo de 300 y `max()` recibia una secuencia
    vacia. La figura moria con "max() arg is an empty sequence" y lo que se
    veia era el cartel gris de `Lienzo.dibujar`.

    El corte tiene que ser relativo a lo que se quiere excluir: el exceso vive
    muy por debajo del borde, asi que un cuarto del primer borde lo deja afuera
    en cualquier ganancia.
    """
    corte = 0.25 * min(r["bordes_ch"])
    y = e["neto"][e["canales"] > corte]
    return float(np.max(y)) if y.size else float(np.max(e["neto"]))


def texto_de(fn, *a, **kw):
    """Captura lo que una funcion de informe IMPRIME, sin duplicar su formato."""
    b = io.StringIO()
    with contextlib.redirect_stdout(b):
        fn(*a, **kw)
    return b.getvalue()


# ==========================================================================
# Lienzo: figura de matplotlib metida en una pestana
# ==========================================================================
def trazar(texto):
    """
    Una linea en `.cache/gui_traza.log`, con la hora y flush inmediato.

    Es la traza de `App._traza`, sacada a nivel de modulo para que tambien la
    escriban los lienzos, que no conocen a la App. Si el proceso muere de golpe
    no queda excepcion que leer, pero SI queda este archivo, y la ultima linea
    dice hasta donde llego. El flush no es opcional: sin el, lo ultimo -que es
    justo lo que interesa- se pierde en el buffer cuando el proceso se va.
    """
    try:
        os.makedirs(rutas.CACHE, exist_ok=True)
        with open(os.path.join(rutas.CACHE, "gui_traza.log"), "a",
                  encoding="utf-8") as f:
            f.write("%s  %s\n" % (time.strftime("%H:%M:%S"), texto))
            f.flush()
    except OSError:
        pass


class _CanvasEncuadrado(FigureCanvasTkAgg):
    """
    TkAgg que, en el instante de renderizar, pone de acuerdo figura, imagen y
    widget, y deja constancia en la traza si tuvo que hacerlo.

    Por que hace falta ademas de `Lienzo._encuadrar`
    -----------------------------------------------
    El recorte de 7.62 tiene una sola mecanica: Agg renderiza lo que mide la
    FIGURA (pulgadas x dpi) y lo vuelca sobre una PhotoImage que mide lo que
    dijo el ultimo `<Configure>`. La PhotoImage tiene tamano fijo, asi que si
    el render es mas grande Tk lo recorta por la derecha y por abajo, sin
    error y sin aviso. `Lienzo._encuadrar` evita que los dos se desacoplen
    reaccionando a `<Configure>` y a `<Map>`, pero eso depende del ORDEN en
    que Tk entrega los eventos, y hay un reporte (14-sep-2026, con un solo
    monitor) que no se pudo reproducir con dos monitores por ningun camino.

    En vez de perseguir el orden, se verifica en el unico momento que
    importa: justo antes de renderizar. Si los tres tamanos coinciden no
    hace nada; si no, corrige los dos que dependen de la figura -su tamano en
    pulgadas y la PhotoImage- y escribe en la traza los tres numeros, que es
    exactamente lo que le faltaba al reporte para poder diagnosticarlo.
    Y una linea por cada tamano NUEVO al que se dibujo cada lienzo, para que
    la traza diga con que tamano se vio cada figura.
    """

    MINIMO = 40      # px: por debajo el widget todavia no esta acomodado

    def __init__(self, figure, master, nombre=None):
        super().__init__(figure, master=master)
        self._nombre = nombre or (lambda: "?")
        self._registrado = None    # (ancho, alto, dpi) del ultimo render anotado
        self._anotar_id = None
        self.forzados = 0          # cuantas veces hubo que corregir

    def _medidas(self):
        w = self._tkcanvas
        an, al = w.winfo_width(), w.winfo_height()
        render = self.get_width_height(physical=True)
        foto = (self._tkphoto.width(), self._tkphoto.height())
        return (an, al), render, foto

    def _acordar(self):
        try:
            (an, al), render, foto = self._medidas()
        except tk.TclError:                  # la ventana se esta cerrando
            return
        if an <= self.MINIMO or al <= self.MINIMO:
            return
        # 1 px de tolerancia en el render: `int(bbox.width)` puede quedar uno
        # abajo por redondeo de `an / dpi * dpi`, y eso no recorta nada
        render_bien = abs(render[0] - an) <= 1 and abs(render[1] - al) <= 1
        if render_bien and foto == (an, al):
            return
        self.forzados += 1
        trazar("!!! reencuadre forzado en %s: widget %dx%d, render %dx%d, "
               "imagen %dx%d, dpi %.1f, ratio %s"
               % (self._nombre(), an, al, render[0], render[1],
                  foto[0], foto[1], self.figure.dpi, self._device_pixel_ratio))
        if not render_bien:
            dpi = self.figure.dpi
            self.figure.set_size_inches(an / dpi, al / dpi, forward=False)
        if foto != (an, al):
            # lo mismo que hace `FigureCanvasTk.resize` al recibir <Configure>
            self._tkcanvas.delete(self._tkcanvas_image_region)
            self._tkphoto.configure(width=an, height=al)
            self._tkcanvas_image_region = self._tkcanvas.create_image(
                an // 2, al // 2, image=self._tkphoto)

    def draw(self):
        self._acordar()
        super().draw()
        # La anotacion va en reposo y no en cada draw: mientras se arrastra el
        # borde de la ventana hay un render por <Configure>, y lo que interesa
        # es el tamano en que quedo, no los cien intermedios.
        try:
            if self._anotar_id is not None:
                self._tkcanvas.after_cancel(self._anotar_id)
            self._anotar_id = self._tkcanvas.after(500, self._anotar)
        except tk.TclError:
            pass

    def _anotar(self):
        self._anotar_id = None
        try:
            (an, al), render, foto = self._medidas()
        except tk.TclError:
            return
        clave = (an, al, round(self.figure.dpi, 1))
        if clave == self._registrado or an <= self.MINIMO:
            return
        self._registrado = clave
        trazar("render %s: widget %dx%d, render %dx%d, imagen %dx%d, dpi %.1f, "
               "ratio %s" % (self._nombre(), an, al, render[0], render[1],
                             foto[0], foto[1], clave[2], self._device_pixel_ratio))


class Lienzo(ttk.Frame):
    """
    Una figura de matplotlib en una pestana, que se dibuja cuando SE VE.

    Por que no se dibuja siempre
    ----------------------------
    Una pestana que Tk todavia no mapeo mide **1x1 px** (medido: el lienzo de
    una pestana de atras devuelve `winfo_width() == 1` hasta que se la
    selecciona). Pintar ahi no muestra nada a nadie y ademas no es gratis:
    TkAgg vuelca el render en una PhotoImage de Tk, y esa imagen la termina
    dibujando Tk **por su cuenta**, sin ningun frame de Python en el medio.
    Cuando el proceso se cae en esa parte, `faulthandler` deja una pila donde
    el hilo principal esta solamente en `mainloop`: no hay a quien culpar y no
    hay mensaje. Eso es exactamente lo que dejo el cierre del 03-sep-2026.

    Cambiarlo por dibujado perezoso no prueba cual era la causa -60 vueltas
    pintando en lienzos de 1x1 no alcanzaron para reproducir la caida-, pero
    saca de la ventana siete blits por corrida que no le servian a nadie, y
    deja el trabajo pesado en el momento en que la pestana tiene tamano de
    verdad. De paso, terminar una calibracion pasa a costar UN dibujo y no
    siete.

    Por que la figura se ENCUADRA a mano antes de cada dibujo
    ---------------------------------------------------------
    TkAgg escala la figura por el *device pixel ratio* de la pantalla: en
    Windows lo saca de `tk scaling / (96/72)` y, como la App declara
    `tk scaling = dpi/72` para que Tk no salga borroso (7.24), en esta
    pantalla vale **1.25**. Pero lo aplica al recibir `<Map>`, o sea DESPUES
    de que `FigureCanvasTk.resize` ya fijo el tamano en pulgadas con el dpi
    viejo, y la unica manera que tiene de reconciliar las dos cosas es
    pedirle al widget que crezca:

        self._tkcanvas.configure(width=ancho_fisico, height=alto_fisico)

    Adentro de un Notebook dentro de un Panedwindow ese pedido no se cumple
    -el tamano lo manda la ventana-, asi que no llega ningun `<Configure>`
    nuevo y la figura queda mas grande que la imagen sobre la que se vuelca.
    El blit **recorta**, sin error y sin aviso. Medido con la ventana en
    1550x830:

        widget 1029x739 px  ·  figura 8.23x5.91 in a 156.25 dpi = 1286x923 px

    o sea el 20% de la derecha y el 20% de abajo no se dibujaban nunca. Eso
    es lo que se veia como "los graficos salen cortados".

    El pedido de tamano tenia ademas un segundo efecto, que se leia como otro
    problema: el Panedwindow terminaba pidiendo 2490x1194 px, mas que la
    pantalla, y en `pack` los hijos que se empaquetan al final son los que se
    quedan sin lugar. Con la ventana sin maximizar desaparecian la barra de
    progreso de abajo y esta misma barra de herramientas.

    Los dos se arreglan con lo mismo: el lienzo **PIDE** `PEDIDO` (un tamano
    chico y fijo) y se conforma con el que le den, y `_encuadrar` le da a la
    figura el tamano exacto del widget con el dpi que la figura tenga en ese
    momento. Corolario, en `App._configurar_nitidez`: el dpi de la figura va
    en 100 -el nominal-, porque el 1.25 de la pantalla lo pone TkAgg.

    Y como todo eso depende del orden en que Tk entrega `<Configure>` y
    `<Map>`, el canvas es un `_CanvasEncuadrado`: comprueba los tres tamanos
    justo antes de renderizar, corrige si hace falta y lo anota en la traza.
    Es la red de seguridad para el caso que no se pudo reproducir (7.63).

    Agrandar un panel
    -----------------
    Aun sin recortes, cuatro paneles en 1029x739 px son ~500x350 cada uno.
    Doble clic sobre uno lo lleva a ocupar la figura entera; doble clic de
    nuevo -o el selector de la barra de abajo- vuelve a la grilla. No se
    recalcula NADA: se guardan las posiciones, se apaga el motor de layout y
    se reubica el panel elegido. Eso es lo unico que lo hace usable en el
    paisaje del chi2, donde volver a llamar a la funcion de dibujo costaria
    ~6300 evaluaciones del modelo.
    """

    MINIMO = 40           # px: por debajo de esto el lienzo no se pinta
    PEDIDO = (360, 240)   # px que el lienzo PIDE; el que ocupa lo da la ventana
    TODOS = "todos los paneles  ·  doble clic para agrandar uno"

    def __init__(self, master, dpi=100):
        super().__init__(master)
        # 12x7 pulgadas es solo la proporcion de partida: `_encuadrar` le pone
        # el tamano del widget antes de cada dibujado.
        self.fig = Figure(figsize=(12, 7), dpi=dpi)
        self.canvas = _CanvasEncuadrado(self.fig, self,
                                        nombre=lambda: self.nombre or "?")
        self.widget = self.canvas.get_tk_widget()

        # La barra de herramientas y el selector de panel comparten fila para
        # no gastar dos veces el alto. El selector NO va adentro de la barra:
        # `NavigationToolbar2Tk._rescale` recorre sus hijos y avisa por log de
        # cualquiera que no sea un boton o una etiqueta de Tk.
        fila = ttk.Frame(self)
        fila.pack(side=tk.BOTTOM, fill=tk.X)
        self.barra = NavigationToolbar2Tk(self.canvas, fila, pack_toolbar=False)
        self.barra.update()
        self.barra.pack(side=tk.LEFT)
        self.v_panel = tk.StringVar(value=self.TODOS)
        self.cb_panel = ttk.Combobox(fila, textvariable=self.v_panel,
                                     state="readonly", width=42)
        self.cb_panel.pack(side=tk.RIGHT, padx=6)
        self.cb_panel.bind("<<ComboboxSelected>>", self._al_elegir_panel)

        self.widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._pendiente = None       # dibujo que espera a que la pestana se vea
        self.nombre = ""             # para la traza
        self._paneles = []           # los ejes "de verdad" de la figura
        self._gemelos = {}           # panel -> ejes que comparten su lugar (twinx)
        self._barras = {}            # panel -> eje de su barra de color
        self._solo = None            # indice del panel agrandado, o None
        self._guardado = None        # posiciones de la grilla, para volver

        # cuando la pestana se muestra (o la ventana cambia de tamano) el
        # widget recibe <Configure>: ahi es cuando el dibujo pendiente puede
        # hacerse con el tamano real. <Map> ademas es donde TkAgg cambia el
        # dpi, asi que hay que reencuadrar despues de el (add="+" garantiza
        # que nuestro manejador corra segundo).
        self.widget.bind("<Configure>", self._al_configurar, add="+")
        self.widget.bind("<Map>", self._al_configurar, add="+")
        self.canvas.mpl_connect("button_press_event", self._al_clic)
        self._encuadrar()
        self.vacio("todavia no hay nada que mostrar")

    # ---------------------------------------------------------- visible --
    def visible(self):
        w = self.widget
        try:
            return (bool(w.winfo_ismapped()) and w.winfo_width() > self.MINIMO
                    and w.winfo_height() > self.MINIMO)
        except tk.TclError:                  # la ventana se esta cerrando
            return False

    def _encuadrar(self):
        """
        Que la figura mida EXACTAMENTE lo que mide el widget que la muestra.

        Las dos mitades del arreglo que describe el docstring de la clase: el
        lienzo no pide mas lugar del que se le puede dar, y la figura toma el
        tamano del widget con el dpi que tenga EN ESE MOMENTO -que no es el
        que se le paso al construirla, porque TkAgg lo multiplica por el
        device pixel ratio-.
        """
        try:
            an, al = self.widget.winfo_width(), self.widget.winfo_height()
            pedido = (self.widget.winfo_reqwidth(),
                      self.widget.winfo_reqheight())
        except tk.TclError:                  # la ventana se esta cerrando
            return
        if pedido != self.PEDIDO:
            self.widget.configure(width=self.PEDIDO[0], height=self.PEDIDO[1])
        if an <= self.MINIMO or al <= self.MINIMO:
            return
        dpi = self.fig.dpi
        if (abs(self.fig.get_figwidth() - an / dpi) > 1e-6
                or abs(self.fig.get_figheight() - al / dpi) > 1e-6):
            # forward=False: el tamano no hay que devolverselo al widget, es
            # el widget el que lo esta imponiendo.
            self.fig.set_size_inches(an / dpi, al / dpi, forward=False)

    def _pintar(self):
        if self.visible():
            self._encuadrar()
            self.canvas.draw_idle()

    def _al_configurar(self, _e=None):
        self._encuadrar()
        if self._pendiente is not None:
            self.asegurar()
        else:
            self._pintar()

    # ------------------------------------------------------------ dibujo --
    def vacio(self, texto):
        self.fig.clear()
        ajustar(self.fig)
        self._encuadrar()
        a = self.fig.add_subplot(111)
        a.text(0.5, 0.5, texto, ha="center", va="center", fontsize=11,
               color="0.45", transform=a.transAxes)
        a.set_axis_off()
        self._paneles, self._gemelos, self._barras = [], {}, {}
        self._guardado = None
        self.cb_panel["values"] = [self.TODOS]
        self.v_panel.set(self.TODOS)
        self._pintar()

    def dibujar(self, fn):
        self.fig.clear()
        ajustar(self.fig)
        self._encuadrar()
        self._guardado = None
        self._barras = {}
        # Se anota que barra de color pertenece a que panel MIENTRAS se
        # dibuja: despues no hay como saberlo (el eje de una barra solo se
        # distingue por su etiqueta "<colorbar>", que no dice de quien es), y
        # sin eso agrandar un panel dejaria su barra escondida o suelta.
        original = self.fig.colorbar

        def colorbar(mappable, cax=None, ax=None, **kw):
            cb = original(mappable, cax=cax, ax=ax, **kw)
            if cax is None and isinstance(ax, Axes):
                self._barras[ax] = cb.ax
            return cb

        self.fig.colorbar = colorbar
        try:
            fn(self.fig)
        except Exception as e:                                  # noqa: BLE001
            # El cartel evita que una figura rota mate la ventana, pero un
            # cartel que nadie mira es un error perdido: la traza va al log.
            registrar_error(type(e), e, e.__traceback__,
                            "dibujo de la pestana %s" % (self.nombre or "?"))
            self.fig.clear()
            a = self.fig.add_subplot(111)
            a.text(0.5, 0.5, "no se pudo dibujar:\n%s\n\n(la traza quedo en "
                             "%s)" % (e, os.path.basename(RUTA_ERRORES)),
                   ha="center", va="center", fontsize=10, color="#b32b2b",
                   transform=a.transAxes)
            a.set_axis_off()
        finally:
            del self.fig.colorbar          # vuelve el metodo de la clase
        self._inventariar()
        self._aplicar_panel()
        self._pintar()

    # ----------------------------------------------------------- diferido --
    def programar(self, fn, texto="calculando..."):
        """Deja el dibujo listo; se hace ahora si se ve, y si no al abrirla."""
        self._pendiente = fn
        if self.visible():
            self.asegurar()
        else:
            self.vacio(texto)

    def asegurar(self):
        """Dibuja lo pendiente, si lo hay y si la pestana ya tiene tamano."""
        if self._pendiente is None or not self.visible():
            return
        fn, self._pendiente = self._pendiente, None
        self.dibujar(fn)

    def limpiar(self, texto):
        """
        Descarta el dibujo pendiente y deja un cartel.

        Hace falta al empezar una corrida nueva: un pendiente de la corrida
        ANTERIOR que no se llego a mirar se dibujaria despues, con los datos
        viejos, en una ventana que ya muestra otro resultado.
        """
        self._pendiente = None
        self.vacio(texto)

    # ------------------------------------------------- agrandar un panel --
    def _inventariar(self):
        """
        Cuales de los ejes de la figura son PANELES, y cual es de quien.

        Tres clases de eje conviven en estas figuras y solo la primera es un
        panel: los de la grilla, los `twinx` (que comparten el lugar de su
        panel y llevan el segundo eje y) y los de las barras de color. Los
        gemelos se detectan por posicion -un `twinx` nace con el mismo
        rectangulo que su panel- y las barras por su etiqueta.
        """
        self._paneles, self._gemelos = [], {}
        for a in self.fig.axes:
            if a.get_label() == "<colorbar>":
                continue
            caja = tuple(round(v, 4) for v in a.get_position().bounds)
            for p in self._paneles:
                if tuple(round(v, 4) for v in p.get_position().bounds) == caja:
                    self._gemelos.setdefault(p, []).append(a)
                    break
            else:
                self._paneles.append(a)

        nombres = [self.TODOS]
        for k, a in enumerate(self._paneles):
            # los paneles de tabla escriben su titulo ADENTRO de los ejes (ver
            # `tabla`), asi que `get_title` no lo tiene
            t = (getattr(a, "_nombre_panel", "") or a.get_title()
                 or "panel %d" % (k + 1)).splitlines()[0]
            nombres.append("%d · %s"
                           % (k + 1, textwrap.shorten(t, 44, placeholder="…")))
        self.cb_panel["values"] = nombres
        if self._solo is not None and self._solo >= len(self._paneles):
            self._solo = None
        self._sincronizar_selector()

    def _sincronizar_selector(self):
        vals = list(self.cb_panel["values"])
        i = 0 if self._solo is None else self._solo + 1
        self.v_panel.set(vals[i] if i < len(vals) else self.TODOS)

    def _al_elegir_panel(self, _e=None):
        i = self.cb_panel.current()
        self.elegir_panel(None if i <= 0 else i - 1)

    def _al_clic(self, ev):
        """Doble clic: agranda el panel de abajo del cursor, o vuelve."""
        if not getattr(ev, "dblclick", False):
            return
        if str(getattr(self.barra, "mode", "")):     # zoom o desplazamiento
            return
        if self._solo is not None:
            self.elegir_panel(None)
            return
        for k, a in enumerate(self._paneles):
            x0, y0, x1, y1 = a.bbox.extents
            if x0 <= ev.x <= x1 and y0 <= ev.y <= y1:
                self.elegir_panel(k)
                return

    def elegir_panel(self, k):
        """`k` = indice del panel a agrandar, o None para volver a la grilla."""
        self._solo = k
        self._sincronizar_selector()
        self._aplicar_panel()
        self._pintar()

    def _aplicar_panel(self):
        if self._solo is None or not self._paneles:
            self._restaurar()
        else:
            self._agrandar(min(self._solo, len(self._paneles) - 1))

    def _agrandar(self, k):
        fig = self.fig
        if self._guardado is None:
            self._guardado = [(a, a.get_position().bounds, a.get_visible(),
                               a.get_anchor()) for a in fig.axes]
        elegido = self._paneles[k]
        familia = [elegido] + list(self._gemelos.get(elegido, ()))
        cb = self._barras.get(elegido)
        # El motor de layout hay que APAGARLO: si no, en el proximo dibujado
        # vuelve a repartir la figura entre todos los ejes de la grilla y
        # deshace la posicion que se acaba de poner.
        fig.set_layout_engine("none")
        for a in fig.axes:
            a.set_visible(any(a is f for f in familia) or a is cb)
        if not elegido.axison:
            # panel de tabla: el titulo y la nota van adentro de los ejes, asi
            # que no hay que reservarle margen a nada
            caja = [0.02, 0.03, 0.96, 0.94]
        else:
            # margen derecho: 0.03 alcanza para el ultimo numero del eje x,
            # pero un twinx o una barra de color necesitan su propio lugar
            der = 0.09 if (cb is not None or self._gemelos.get(elegido)) else 0.03
            caja = [0.07, 0.10, 1.0 - 0.07 - der, 0.81]
        for a in familia:
            a.set_position(caja)
        if cb is not None:
            # Un eje con aspecto fijo (el `imshow` de la matriz de
            # correlacion) no ocupa la caja que se le da, sino la mayor que
            # entre adentro con su relacion de lados. `apply_aspect` la
            # calcula sin necesidad de dibujar, y es la que hay que seguir: si
            # no, la barra de color queda al otro lado de la figura y con el
            # triple de alto que el mapa que describe.
            ref = caja
            if elegido.get_aspect() != "auto":
                # centrado: con una barra de color al lado, matplotlib le deja
                # el ancla en el este, y en una caja ancha eso deja todo el
                # hueco de un solo lado
                elegido.set_anchor("C")
                elegido.apply_aspect()
                ref = elegido.get_position().bounds
            cb.set_position([min(ref[0] + ref[2] + 0.015, 0.955), ref[1],
                             0.02, ref[3]])

    def _restaurar(self):
        if self._guardado is None:
            return
        for a, caja, visible, ancla in self._guardado:
            a.set_position(caja)
            a.set_visible(visible)
            a.set_anchor(ancla)
        self._guardado = None
        ajustar(self.fig)


def tabla(ax, filas, columnas, anchos, titulo=None, resaltar=(), nota=None):
    """
    Una tabla de matplotlib con el formato que usan todas las pestanas.

    Existe para que las tablas de las pestanas que las usan salgan iguales, y
    porque `ax.table` sin `bbox` y sin apagar el autosize deja las celdas de
    altura arbitraria segun cuanto texto entre.

    `resaltar` son indices de fila (contando desde 0 en `filas`) que se pintan
    de fondo: se usa para la fila del arranque ganador del recocido, que es lo
    unico que hay que encontrar de un vistazo en esa tabla.

    El titulo y la nota se dibujan ADENTRO de los ejes, y la tabla se corre
    para dejarles el hueco. No es capricho: un titulo de eje largo con
    `loc="left"` se sale del ancho de sus ejes y se le monta al titulo del
    panel de al lado, y constrained layout no lo acorta ni lo parte en lineas.
    Metiendo los dos adentro, lo unico que el motor de layout ve es la caja de
    los ejes.

    Ojo con el `bbox` de `ax.table`: es (x0, y0, ANCHO, ALTO), no
    (x0, y0, x1, y1). Pasarle (0, 0.09, 1, 1) -que es lo que parece razonable
    para "dejar 9% abajo"- da una tabla de alto 1 que arranca en 0.09, o sea
    que se sale por ARRIBA de los ejes y tapa el titulo.
    """
    ax.set_axis_off()
    # El titulo va adentro de los ejes (ver mas abajo), asi que `get_title`
    # queda vacio y el selector de paneles de `Lienzo` no tendria como nombrar
    # a este panel. Se lo deja anotado aparte.
    ax._nombre_panel = titulo or ""
    # La nota se parte en lineas cortas, y no es cosmetico: un `ax.text` cuenta
    # como decoracion del eje, asi que constrained layout mide el tight bbox
    # CON el texto adentro y, si el texto es mas ancho que el eje, achica el
    # eje hasta que entre. Con una nota de 100 caracteres en una linea la tabla
    # quedaba de un quinto del panel.
    if nota:
        nota = "\n".join(textwrap.fill(l, 46) for l in nota.splitlines())
    y0, y1 = 0.0, 1.0
    if titulo:
        ax.text(0.0, 1.0, titulo, transform=ax.transAxes, va="top",
                fontsize=10, color="#1f4e79")
        y1 = 0.945
    if nota:
        ax.text(0.0, 0.0, nota, transform=ax.transAxes, va="bottom",
                fontsize=7.5, color="#5b6b7c")
        # 0.048 por linea: a 7.5 pt con interlineado normal cada linea mide
        # ~0.14 pulgadas, y en la ventana chica el eje de una grilla 2x2 mide
        # ~2.9, o sea ~5% por linea. Reservar de menos hace que la nota se
        # meta encima de la ultima fila de la tabla.
        y0 = 0.035 + 0.048 * (nota.count("\n") + 1)
    t = ax.table(cellText=filas, colLabels=columnas, colWidths=anchos,
                 cellLoc="left", colLoc="left", loc="upper left",
                 bbox=[0, y0, 1, y1 - y0])
    t.auto_set_font_size(False)
    # El cuerpo sale del ANCHO REAL del panel. matplotlib no achica el texto
    # de una celda para que entre: lo RECORTA, sin avisar y sin dejar rastro
    # -asi se veian "punto INICIAI" y "A (VARPR"-, y `auto_set_font_size(True)`
    # tampoco sirve porque solo sabe achicar: una vez que encogio no vuelve a
    # crecer aunque la ventana se agrande.
    #
    # Aca decia que medir el ancho no servia "porque la figura se dibuja una
    # vez y despues el lienzo la reescala". Eso valia mientras la figura se
    # construia con el tamano de fabrica (12x7): ahora `Lienzo._encuadrar` le
    # pone el tamano del widget ANTES de llamar a la funcion de dibujo, asi
    # que lo que se mide aca es el tamano de verdad.
    #
    # El 1.28 convierte el ancho que la grilla le asigna al panel en el que
    # constrained layout le va a dar de verdad (medido: 0.405 -> 0.517 de la
    # figura en una grilla de 2x2 con razon 1.35:1). Y 5.8 pulgadas es el
    # ancho al que 8 pt entran sin recortarse, que es el caso conocido: la
    # ventana maximizada.
    ancho = 1.28 * ax.get_position().width * ax.figure.get_figwidth()
    t.set_fontsize(min(8.0, max(5.8, 8.0 * ancho / 5.8)))
    for (fila, col), celda in t.get_celld().items():
        celda.set_edgecolor("#dde3e8")
        if fila == 0:
            celda.set_text_props(weight="bold", color="#1f4e79")
            celda.set_facecolor("#eef3f7")
        elif (fila - 1) in resaltar:
            celda.set_facecolor("#e8f4ec")
            celda.set_text_props(weight="bold")
    return t


# ==========================================================================
# Panel de etapas
# ==========================================================================
class PanelEtapas(ttk.LabelFrame):
    """Una fila por etapa: nombre, barra propia y estado/tiempo."""

    COLORES = {"pendiente": "#8a97a3", "corriendo": "#1f4e79",
               "hecho": "#1f7a3d", "cache": "#1f7a3d",
               "cancelado": "#c4630a", "error": "#b32b2b"}

    def __init__(self, master):
        super().__init__(master, text=" Etapas ", padding=(8, 6))
        self.filas = {}
        for i, (clave, nombre) in enumerate(pipeline.ETAPAS):
            ttk.Label(self, text="%d" % (i + 1), width=2,
                      foreground="#8a97a3").grid(row=i, column=0, sticky="w")
            lab = ttk.Label(self, text=nombre, width=25, anchor="w")
            lab.grid(row=i, column=1, sticky="w", padx=(2, 6))
            bar = ttk.Progressbar(self, length=90, maximum=100.0)
            bar.grid(row=i, column=2, sticky="ew", pady=1)
            est = ttk.Label(self, text="—", width=17, anchor="w",
                            foreground=self.COLORES["pendiente"])
            est.grid(row=i, column=3, sticky="w", padx=(6, 0))
            self.filas[clave] = dict(lab=lab, bar=bar, est=est, indet=False)
        self.columnconfigure(2, weight=1)

    def reiniciar(self):
        for f in self.filas.values():
            self._determinada(f)
            f["bar"]["value"] = 0
            f["est"].config(text="—", foreground=self.COLORES["pendiente"])
            f["lab"].config(foreground="")

    def _determinada(self, f):
        if f["indet"]:
            f["bar"].stop()
            f["bar"].config(mode="determinate")
            f["indet"] = False

    def comenzar(self, clave):
        f = self.filas.get(clave)
        if not f:
            return
        f["est"].config(text="corriendo...", foreground=self.COLORES["corriendo"])
        f["lab"].config(foreground=self.COLORES["corriendo"])

    def avance(self, clave, frac, texto=""):
        f = self.filas.get(clave)
        if not f:
            return
        if frac is None:
            # Sin fraccion conocida, la barra va y viene: no informa cuanto
            # falta, pero si informa que el programa esta vivo, que es la mitad
            # de lo que se le pide.
            if not f["indet"]:
                f["bar"].config(mode="indeterminate")
                f["bar"].start(14)
                f["indet"] = True
        else:
            self._determinada(f)
            f["bar"]["value"] = 100.0 * frac
        if texto:
            f["est"].config(text=texto[:17])

    def terminar(self, clave, segundos, nota="", estado="hecho"):
        f = self.filas.get(clave)
        if not f:
            return
        self._determinada(f)
        f["bar"]["value"] = 100.0 if estado in ("hecho", "cache") else f["bar"]["value"]
        # 17 caracteres: el detalle largo (cuantos pulsos se leyeron) ya sale
        # en el informe y en la barra de abajo; lo que hace falta ver por fila
        # es el TIEMPO y si fue cache.
        txt = "%.1f s%s" % (segundos, "  ·cache" if "cache" in nota else "")
        f["est"].config(text=txt[:17], foreground=self.COLORES[estado])
        f["lab"].config(foreground="")

    def marcar(self, clave, estado, texto):
        f = self.filas.get(clave)
        if not f:
            return
        self._determinada(f)
        f["est"].config(text=texto[:17], foreground=self.COLORES[estado])


# ==========================================================================
# La aplicacion
# ==========================================================================
class App(tk.Tk):
    # Opcion del selector de la pestana Registro que superpone los detectores.
    COMPARAR = "Comparar detectores"
    GANANCIA = "Ganancia: log-log y relativa al A"

    def __init__(self):
        super().__init__()
        self.title("Calibracion por borde Compton  ·  Monte Carlo")
        self.geometry("1500x900")
        self.minsize(1080, 700)
        # Maximizada de entrada: lo que se mira en esta ventana son grillas de
        # 2x2, y en 1500x900 cada panel queda de ~430x250 px, que es donde los
        # rotulos empiezan a estorbar. Con la pantalla completa el mismo
        # contenido respira sin cambiar nada del dibujo.
        try:
            self.state("zoomed")               # Windows y la mayoria de los WM
        except tk.TclError:
            try:
                self.attributes("-zoomed", True)
            except tk.TclError:
                pass

        self.cola = queue.Queue()
        self.cancelar = threading.Event()
        self.hilo = None
        self.salida = None
        self._etapa_actual = None

        self._configurar_nitidez()
        self._construir()
        self._vaciar_id = self.after(80, self._vaciar_cola)
        self._al_elegir_corrida()
        self._avisar_fallo_anterior()
        # Pantalla, ventana y escalado a la traza, y cada vez que la ventana
        # cambie de tamano o de monitor. Con el reporte de 7.63 no habia forma
        # de saber en que condiciones se habia dibujado lo que se veia mal.
        self._cfg_id = None
        self._ultima_ventana = None
        self.bind("<Configure>", self._al_configurar_ventana, add="+")
        self.after_idle(self._trazar_entorno)

    def _avisar_fallo_anterior(self):
        """
        Si la sesion anterior se cerro de golpe, DECIRLO al abrir.

        Un cierre por violacion de acceso no deja excepcion ni mensaje: lo
        unico que queda es el volcado de `faulthandler`, y hasta ahora nadie
        avisaba de que existiera. La consecuencia practica era que el programa
        "se cerraba solo y sin decir nada" teniendo la pila escrita en disco.
        """
        if not FALLO_ANTERIOR:
            return
        self._traza("!!! la sesion anterior termino en un cierre abrupto")
        self.v_aviso.set(
            "La sesion anterior se cerro de golpe (no fue un error de Python: "
            "no hay excepcion que mostrar). La pila quedo en %s y el detalle "
            "esta abajo, en Resultado."
            % rutas.relativa(RUTA_FALLO_ANTERIOR))
        self._escribir(
            "=" * 70 + "\nLA SESION ANTERIOR SE CERRO DE GOLPE\n" + "=" * 70
            + "\nEsto es el volcado de faulthandler, no una excepcion de "
              "Python.\nLa ultima linea de gui_traza.log dice hasta donde "
              "habia llegado.\n\n" + FALLO_ANTERIOR + "\n")

    # ------------------------------------------------------------ nitidez --
    def _configurar_nitidez(self):
        dpi = self.winfo_fpixels("1i")
        self.tk.call("tk", "scaling", dpi / 72.0)
        # El dpi de las figuras va en el NOMINAL. Aca decia
        # `max(int(100 * dpi / 96), 100)` -o sea 125 en esta pantalla- para
        # compensar el escalado de Windows, pero TkAgg ya lo compensa solo:
        # multiplica el dpi de la figura por `tk scaling / (96/72)`, que con
        # la linea de arriba da exactamente ese mismo 1.25. El factor entraba
        # dos veces, la figura terminaba a 156 dpi y todo se dibujaba un 25%
        # mas grande de lo pedido, o sea con un 36% menos de superficie util
        # (y recortado, ver `Lienzo`).
        self.dpi_fig = 100
        s = ttk.Style()
        if "vista" in s.theme_names():
            s.theme_use("vista")
        s.configure("TLabelframe.Label", foreground="#1f4e79")
        s.configure("Titulo.TLabel", foreground="#1f4e79",
                    font=("Segoe UI", 10, "bold"))

    # ----------------------------------------------------------- widgets --
    def _construir(self):
        # La barra de abajo se empaqueta PRIMERO, y no al final. `pack` reparte
        # el espacio en el orden en que se lo piden y al que llega ultimo, si
        # ya no queda, no le toca nada: con la ventana sin maximizar la barra
        # de progreso y el boton de Cancelar no se dibujaban (medido: la
        # columna izquierda pide 907 px de alto y la ventana tenia 830).
        self._barra_inferior()

        raiz = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        raiz.pack(fill=tk.BOTH, expand=True, padx=8, pady=(8, 4))

        izq, dentro = self._columna_izquierda(raiz)
        raiz.add(izq, weight=0)
        self._panel_entradas(dentro)

        # La lista de etapas va ACA, debajo de las entradas, y no arriba de los
        # graficos. Medido en esta pantalla con la ventana maximizada
        # (1920x991): con la lista cruzada arriba el lienzo queda de 11.95 x
        # 5.53 pulgadas, o sea filas de 2.75 pulgadas en una grilla de 2x2, que
        # es donde los rotulos empiezan a chocar. La lista se lleva ~170 px de
        # ALTO y solo necesita ~430 de ancho, y la columna izquierda tiene ese
        # ancho y le sobra alto: moverla le devuelve el 25% de altura al
        # grafico y no le saca nada a nadie.
        self.izq = dentro
        self.etapas = PanelEtapas(dentro)
        self.etapas.pack(fill=tk.X, pady=(10, 0))

        der = ttk.Panedwindow(raiz, orient=tk.VERTICAL)
        raiz.add(der, weight=1)
        self.der = der
        self._panel_resultados(der)

    def _columna_izquierda(self, padre):
        """
        La columna de controles, con desplazamiento vertical si no entra.

        Un `Panedwindow` reparte el ANCHO entre sus paneles, pero el alto se
        lo impone la ventana: si los controles piden mas alto del que hay -en
        este equipo piden 907 px y la ventana sin maximizar mide 830-, lo que
        sobra se corta abajo y no hay forma de llegar a el. Adentro de un
        Canvas si la hay, y encima la barra de desplazamiento solo aparece
        cuando de verdad hace falta, asi que con la ventana maximizada no se
        nota que este.

        Devuelve `(contenedor, marco)`: el contenedor va al Panedwindow, y los
        paneles se cuelgan del marco.
        """
        cont = ttk.Frame(padre, padding=(0, 0, 8, 0))
        lienzo = tk.Canvas(cont, highlightthickness=0, borderwidth=0,
                           width=10, height=10)
        barra = ttk.Scrollbar(cont, orient=tk.VERTICAL, command=lienzo.yview)
        lienzo.configure(yscrollcommand=barra.set)
        lienzo.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        marco = ttk.Frame(lienzo)
        ventana = lienzo.create_window((0, 0), window=marco, anchor="nw")

        pendiente = []

        def revisar():
            # En reposo y no en el acto: el alto del Canvas se lee con
            # `winfo_height` y adentro de un manejador de <Configure> todavia
            # puede ser el de antes, asi que decidir ahi si la barra hace
            # falta la dejaba puesta con la ventana ya grande (medido: hacia
            # falta UN <Configure> mas para que se escondiera).
            del pendiente[:]
            if not lienzo.winfo_exists():
                return
            caja = lienzo.bbox("all")
            lienzo.configure(scrollregion=caja or (0, 0, 0, 0))
            # el ancho lo manda el contenido (un Canvas no propaga el tamano
            # de lo que tiene adentro, hay que pedirselo a mano); el alto, la
            # ventana
            ancho = marco.winfo_reqwidth()
            if lienzo.winfo_reqwidth() != ancho:
                lienzo.configure(width=ancho)
            falta = bool(caja) and caja[3] > lienzo.winfo_height()
            if falta and not barra.winfo_ismapped():
                barra.pack(side=tk.RIGHT, fill=tk.Y)
            elif not falta and barra.winfo_ismapped():
                barra.pack_forget()
                lienzo.yview_moveto(0.0)

        def programar(_e=None):
            if not pendiente:
                pendiente.append(lienzo.after_idle(revisar))

        def al_estirar(e):
            if lienzo.itemcget(ventana, "width") != str(e.width):
                lienzo.itemconfigure(ventana, width=e.width)
            programar()

        marco.bind("<Configure>", programar)
        lienzo.bind("<Configure>", al_estirar)
        self._col_lienzo = lienzo
        # La rueda se ata al TOPLEVEL y el manejador decide por donde esta el
        # cursor, que es lo mismo que hace TkAgg con la suya: atarla al Canvas
        # no alcanza, porque los eventos los recibe el widget de abajo del
        # cursor -un combo, una etiqueta- y no su contenedor.
        self.bind("<MouseWheel>", self._rueda, add="+")
        return cont, marco

    def _rueda(self, ev):
        """La rueda desplaza la columna solo si el cursor esta sobre ella."""
        # Se compara contra el RECTANGULO del Canvas en coordenadas de
        # pantalla y no con `winfo_containing`, que es lo que usa TkAgg para
        # lo mismo: `winfo_containing` devuelve None cuando la ventana no esta
        # al frente, y entonces la rueda deja de andar sin motivo aparente.
        c = self._col_lienzo
        try:
            x, y = c.winfo_rootx(), c.winfo_rooty()
            dentro = (x <= ev.x_root < x + c.winfo_width()
                      and y <= ev.y_root < y + c.winfo_height())
        except tk.TclError:
            return
        if dentro:
            c.yview_scroll(-1 if ev.delta > 0 else 1, "units")

    def _panel_entradas(self, p):
        f = ttk.LabelFrame(p, text=" Datos ", padding=(10, 8))
        f.pack(fill=tk.X)

        ttk.Label(f, text="Corrida").grid(row=0, column=0, sticky="w")
        self.v_corrida = tk.StringVar()
        claves = sorted(rutas.disponibles()) + [MANUAL]
        self.cb_corrida = ttk.Combobox(f, textvariable=self.v_corrida,
                                       values=claves, state="readonly", width=28)
        self.cb_corrida.grid(row=0, column=1, columnspan=2, sticky="ew", pady=2)
        self.cb_corrida.set(claves[0] if claves else MANUAL)
        self.cb_corrida.bind("<<ComboboxSelected>>",
                             lambda e: self._al_elegir_corrida())

        self.v_arch = {}
        # "Datos del DRS4 (*.*)" va PRIMERO a proposito: los binarios de la
        # sesion del 1-sep-2026 no tienen extension en el nombre (`Co60`,
        # `Cs137`, `B1700V`), asi que con el filtro "*.xml" por defecto la
        # carpeta se ve VACIA y parece que no hay datos. El formato lo decide
        # `drs.formato` mirando los primeros 4 bytes, no el nombre.
        TIPOS_DATOS = (("Datos del DRS4 (XML o binario)", "*.*"),
                       ("XML del DRSOSC", "*.xml"))
        filas = (("muestra", "Muestra", TIPOS_DATOS),
                 ("fondo", "Fondo (opcional)", TIPOS_DATOS),
                 ("sim", "Simulacion .root", (("ROOT", "*.root"),)))
        for i, (k, etiqueta, tipos) in enumerate(filas, start=1):
            ttk.Label(f, text=etiqueta).grid(row=i, column=0, sticky="w")
            v = tk.StringVar()
            self.v_arch[k] = v
            e = ttk.Entry(f, textvariable=v, width=24)
            e.grid(row=i, column=1, sticky="ew", pady=2)
            ttk.Button(f, text="...", width=3,
                       command=lambda kk=k, tt=tipos: self._elegir(kk, tt)
                       ).grid(row=i, column=2, padx=(4, 0))

        # Que formato y que canal se detectaron. Es el diagnostico mas barato que
        # existe y evita las dos fallas silenciosas conocidas: leer el canal
        # equivocado (7.36) y no ver los binarios porque el nombre no dice nada.
        self.v_formato = tk.StringVar(value="")
        ttk.Label(f, textvariable=self.v_formato, wraplength=270,
                  foreground="#5b6b7c", justify="left").grid(
            row=4, column=0, columnspan=3, sticky="w", pady=(4, 0))

        ttk.Separator(f, orient="horizontal").grid(row=5, column=0, columnspan=3,
                                                   sticky="ew", pady=8)

        ttk.Label(f, text="Detector").grid(row=6, column=0, sticky="w")
        self.v_detector = tk.StringVar()
        ttk.Entry(f, textvariable=self.v_detector).grid(
            row=6, column=1, columnspan=2, sticky="ew", pady=2)

        # El HV va aparte del nombre del detector: D no es del detector, es del
        # detector EN UN PUNTO DE OPERACION. Sin este campo dos corridas del
        # mismo tubo a voltajes distintos entrarian a la misma recta.
        ttk.Label(f, text="Alto voltaje [V]").grid(row=7, column=0, sticky="w")
        self.v_hv = tk.StringVar()
        ttk.Entry(f, textvariable=self.v_hv).grid(
            row=7, column=1, columnspan=2, sticky="ew", pady=2)

        ttk.Label(f, text="Fuente").grid(row=8, column=0, sticky="w")
        self.v_fuente = tk.StringVar(value="Co-60")
        ttk.Combobox(f, textvariable=self.v_fuente,
                     values=sorted(fisica.FUENTES), state="readonly").grid(
            row=8, column=1, columnspan=2, sticky="ew", pady=2)

        ttk.Label(f, text="Escala").grid(row=9, column=0, sticky="w")
        self.v_escala = tk.StringVar(value=list(ESCALAS)[0])
        ttk.Combobox(f, textvariable=self.v_escala, values=list(ESCALAS),
                     state="readonly").grid(row=9, column=1, columnspan=2,
                                            sticky="ew", pady=2)
        f.columnconfigure(1, weight=1)

        # ---- opciones -----------------------------------------------------
        g = ttk.LabelFrame(p, text=" Ajuste ", padding=(10, 8))
        g.pack(fill=tk.X, pady=(8, 0))
        self.v_metodo = tk.StringVar(value="global")
        self.v_bins = tk.IntVar(value=500)
        self.v_inicio = tk.DoubleVar(value=0.45)
        self.v_arranques = tk.IntVar(value=8)
        # Cuanto se integra DESPUES del onset de cada pulso, en ns. El default
        # es COMPUERTA_DEFECTO y se edita aca, corrida por corrida; vaciar el
        # campo devuelve el comportamiento viejo (cada archivo mide la suya
        # sobre la convergencia de su propia integral acumulada).
        #
        # Que el valor este fijo y a la vista importa para COMPARAR: el `post`
        # que mide cada archivo sale de su propio pulso, y dos corridas que
        # integran ventanas distintas recogen fracciones distintas de la carga,
        # lo que entra derecho en D. Medido en el Cs-137 del detector A, cuyas
        # tres corridas midieron 55, 98 y 124 ns cada una por su cuenta: la ley
        # de ganancia daba b = 10.4 con residuos de +-14%; con compuerta comun,
        # b = 10.1 con +-1.5%.
        self.v_post = tk.StringVar(value="%g" % COMPUERTA_DEFECTO)
        campos = (("Metodo del area", ttk.Combobox, self.v_metodo,
                   dict(values=["global", "convergencia"], state="readonly")),
                  ("Bines del espectro", ttk.Spinbox, self.v_bins,
                   dict(from_=100, to=2000, increment=50)),
                  ("inicio_rel", ttk.Spinbox, self.v_inicio,
                   dict(from_=0.20, to=0.90, increment=0.05, format="%.2f")),
                  ("Arranques del SA", ttk.Spinbox, self.v_arranques,
                   dict(from_=1, to=32, increment=1)),
                  ("Compuerta tras el onset [ns]", ttk.Entry, self.v_post, {}))
        for i, (etiqueta, clase, var, kw) in enumerate(campos):
            ttk.Label(g, text=etiqueta).grid(row=i, column=0, sticky="w")
            clase(g, textvariable=var, width=14, **kw).grid(
                row=i, column=1, sticky="ew", pady=2)
        g.columnconfigure(1, weight=1)

        self.v_forzar = tk.BooleanVar(value=False)
        ttk.Checkbutton(g, text="usar la simulacion aunque no\ncoincida con la fuente",
                        variable=self.v_forzar).grid(row=len(campos), column=0,
                                                     columnspan=2, sticky="w",
                                                     pady=(6, 0))

        h = ttk.Frame(p)
        h.pack(fill=tk.X, pady=(10, 0))
        self.b_calibrar = ttk.Button(h, text="Calibrar", command=self._calibrar)
        self.b_calibrar.pack(fill=tk.X, ipady=4)

        self.v_aviso = tk.StringVar(value="")
        ttk.Label(p, textvariable=self.v_aviso, wraplength=290,
                  foreground="#b32b2b", justify="left").pack(fill=tk.X,
                                                             pady=(8, 0))

    def _panel_resultados(self, p):
        self.nb = ttk.Notebook(p)
        p.add(self.nb, weight=1)

        # --- Resultado: informe + figura ---
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Resultado")
        pw = ttk.Panedwindow(tab, orient=tk.VERTICAL)
        pw.pack(fill=tk.BOTH, expand=True)
        cont = ttk.Frame(pw)
        # 9 lineas y no 14: el informe tiene ~29 y hay que desplazarlo de todas
        # formas, asi que lo que se gana con mas alto es cero y lo que se pierde
        # es figura. El divisor se puede arrastrar si hace falta leerlo entero.
        self.txt = tk.Text(cont, height=9, wrap="none",
                           font=("Consolas", 9), relief="flat",
                           background="#fbfcfd")
        sy = ttk.Scrollbar(cont, orient="vertical", command=self.txt.yview)
        self.txt.configure(yscrollcommand=sy.set)
        sy.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt.pack(fill=tk.BOTH, expand=True)
        pw.add(cont, weight=0)
        self.lz_resultado = Lienzo(pw, self.dpi_fig)
        self.lz_resultado.nombre = "resultado"
        pw.add(self.lz_resultado, weight=1)

        # --- diagnosticos ---
        # Las cuatro numeradas son la cadena en el orden en que corre: que
        # entro (datos y plantilla), como se eligio la ventana (NCC), como se
        # busco el optimo (recocido) y cuanto mejoro respecto del punto de
        # partida. Las otras tres son transversales.
        self.lienzos = {}
        for clave, titulo in (("ventana", "Ventana de integracion"),
                              ("espectro", "1 · Datos y simulacion"),
                              ("ncc", "2 · NCC e intervalo"),
                              ("sa", "3 · Semilla y recocido"),
                              ("comparacion", "4 · Inicial vs final"),
                              ("parametros", "Θ y correlaciones"),
                              ("paisaje", "Paisaje del χ²")):
            lz = Lienzo(self.nb, self.dpi_fig)
            lz.nombre = clave
            self.nb.add(lz, text=titulo)
            self.lienzos[clave] = lz

        # --- registro ---
        self._tab_registro()
        self.nb.bind("<<NotebookTabChanged>>", self._al_cambiar_pestana)

    def _tab_registro(self):
        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Registro")
        arriba = ttk.Frame(tab)
        arriba.pack(fill=tk.X, padx=6, pady=6)
        ttk.Button(arriba, text="Actualizar",
                   command=self._refrescar_registro).pack(side=tk.LEFT)
        ttk.Button(arriba, text="Borrar la seleccionada",
                   command=self._borrar_registro).pack(side=tk.LEFT, padx=6)
        ttk.Label(arriba, text="Ver:").pack(side=tk.LEFT, padx=(12, 4))
        self.v_ver_reg = tk.StringVar(value=self.COMPARAR)
        self.cb_ver_reg = ttk.Combobox(arriba, textvariable=self.v_ver_reg,
                                       state="readonly", width=22,
                                       values=[self.COMPARAR])
        self.cb_ver_reg.pack(side=tk.LEFT)
        self.cb_ver_reg.bind("<<ComboboxSelected>>",
                             lambda _e: self._fig_linealidad())
        ttk.Label(arriba, text="Test A: dos FUENTES al mismo punto de "
                               "operacion.   Test B: la misma FUENTE a varios HV.",
                  foreground="#5b6b7c").pack(side=tk.LEFT, padx=12)

        cols = ("punto", "fuente", "escala", "D", "err", "borde", "chi2", "fecha")
        self.tv = ttk.Treeview(tab, columns=cols, show="headings", height=8)
        anchos = (170, 80, 60, 90, 80, 100, 70, 120)
        for c, a in zip(cols, anchos):
            self.tv.heading(c, text=c)
            self.tv.column(c, width=a, anchor="w")
        self.tv.pack(fill=tk.X, padx=6)

        # Sin el panel de texto del informe: le quitaba ~58 columnas al lienzo
        # y los dos tests quedaban apretados. El informe sigue en
        # `python analisis/registro.py`.
        self.lz_registro = Lienzo(tab, self.dpi_fig)
        self.lz_registro.nombre = "registro"
        self.lz_registro.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

    def _barra_inferior(self):
        f = ttk.Frame(self, padding=(10, 4, 10, 8))
        # side=BOTTOM explicito: esta barra se empaqueta ANTES que el resto
        # (ver `_construir`), asi que sin decirlo se iria arriba de todo.
        f.pack(side=tk.BOTTOM, fill=tk.X)
        self.bar_global = ttk.Progressbar(f, maximum=100.0)
        self.bar_global.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.v_estado = tk.StringVar(value="listo")
        ttk.Label(f, textvariable=self.v_estado, width=52, anchor="w").pack(
            side=tk.LEFT, padx=10)
        self.b_cancelar = ttk.Button(f, text="Cancelar", state="disabled",
                                     command=self._cancelar)
        self.b_cancelar.pack(side=tk.RIGHT)
        # La lista de etapas solo dice algo mientras el calculo corre. Ya no le
        # quita alto al grafico (vive en la columna izquierda), asi que esto es
        # para ventanas chicas, donde la columna izquierda no alcanza para las
        # entradas y la lista a la vez. No se pliega sola al terminar a
        # proposito: los tiempos por etapa son parte del diagnostico, y una
        # ventana que se reacomoda sin que nadie la toque se lee como un error.
        self.v_ver_etapas = tk.BooleanVar(value=True)
        ttk.Checkbutton(f, text="ver etapas", variable=self.v_ver_etapas,
                        command=self._alternar_etapas).pack(side=tk.RIGHT,
                                                            padx=10)

    def _alternar_etapas(self):
        """Muestra u oculta la lista de etapas, sin recalcular nada."""
        if self.v_ver_etapas.get():
            self.etapas.pack(fill=tk.X, pady=(10, 0))
        else:
            self.etapas.pack_forget()

    # ------------------------------------------------------------ entradas --
    def _elegir(self, clave, tipos):
        r = filedialog.askopenfilename(title="Elegir %s" % clave,
                                       filetypes=list(tipos))
        if r:
            self.v_arch[clave].set(r)
            self.v_corrida.set(MANUAL)
            self._refrescar_formato()

    def _refrescar_formato(self):
        """
        Que formato y que canal tiene lo elegido, en una linea.

        Cuesta milisegundos (el XML lee el primer evento, el binario la
        cabecera) y ataja las dos fallas que en este proyecto no dan error sino
        un numero equivocado: el canal mal leido (7.36) y un archivo que no es
        del DRS4. Sin fondo tambien lo dice, porque cambia lo que significa el
        espectro.
        """
        partes = []
        for clave in ("muestra", "fondo"):
            ruta = self.v_arch[clave].get().strip()
            if not ruta:
                partes.append("%s: -" % clave if clave != "fondo"
                              else "fondo: no hay (el neto es la muestra cruda)")
                continue
            if not os.path.exists(ruta):
                partes.append("%s: NO EXISTE" % clave)
                continue
            try:
                fmt = drs.formato(ruta)
                canal, tag = drs.detectar_canal(ruta)
                partes.append("%s: %s, %s/%s" % (clave, fmt, canal,
                                                 tag or "sin scaler"))
            except Exception as e:                             # noqa: BLE001
                partes.append("%s: ilegible (%s)" % (clave, str(e)[:50]))
        self.v_formato.set('\n'.join(partes))

    def _al_elegir_corrida(self):
        k = self.v_corrida.get()
        if k == MANUAL:
            self._refrescar_formato()
            return
        c = rutas.corrida(k)
        self.v_arch["muestra"].set(c.muestra or "")
        self.v_arch["fondo"].set(c.fondo or "")
        self.v_arch["sim"].set(c.sim or "")
        self.v_detector.set(c.detector or k)
        self.v_hv.set("" if c.hv is None else "%g" % c.hv)
        if c.fuente in fisica.FUENTES:
            self.v_fuente.set(c.fuente)
        self.v_aviso.set("")
        self._refrescar_formato()

    def _corrida_actual(self):
        """Arma la Corrida con lo que hay en el panel, venga del catalogo o no."""
        m = self.v_arch["muestra"].get().strip()
        if not m or not os.path.exists(m):
            raise ValueError("falta el archivo de la muestra")
        fondo = self.v_arch["fondo"].get().strip() or None
        sim = self.v_arch["sim"].get().strip() or None
        det = self.v_detector.get().strip()
        if not det:
            raise ValueError("hay que nombrar el detector: es la clave con la "
                             "que se guarda en el registro, y sin ella no se "
                             "puede saber despues que dos calibraciones son "
                             "del mismo aparato.")
        fuente = self.v_fuente.get()
        k = self.v_corrida.get()
        if k == MANUAL:
            k = re.sub(r"[^A-Za-z0-9_.-]", "_", "%s_%s" % (det, fuente))
        hv = self.v_hv.get().strip()
        try:
            hv = float(hv) if hv else None
        except ValueError:
            raise ValueError("el alto voltaje tiene que ser un numero (o vacio "
                             "si no se registro): '%s'" % hv)
        base = rutas.CORRIDAS.get(self.v_corrida.get())
        return rutas.Corrida(k, "%s · %s" % (det, fuente), muestra=m,
                             fondo=fondo, fuente=fuente, detector=det, sim=sim,
                             hv=hv,
                             # el resto de la cadena solo se conoce si la
                             # corrida viene del catalogo; con archivos sueltos
                             # queda sin registrar, y `misma_cadena` avisa
                             atenuador=getattr(base, "atenuador", None),
                             umbral=getattr(base, "umbral", None),
                             muestreo=getattr(base, "muestreo", None))

    def _opciones(self):
        post = self.v_post.get().strip()
        try:
            post = float(post) if post else None
        except ValueError:
            raise ValueError("la compuerta tiene que ser un numero de ns "
                             "(o vacia para que cada archivo mida la suya): "
                             "'%s'" % post)
        if post is not None and post <= 0:
            raise ValueError("la compuerta tiene que ser positiva: '%g'" % post)
        return dict(metodo=self.v_metodo.get(),
                    rama=ESCALAS[self.v_escala.get()],
                    bins=int(self.v_bins.get()),
                    inicio_rel=float(self.v_inicio.get()),
                    n_arranques=int(self.v_arranques.get()),
                    post_fijo=post,
                    # `exceso` va explicito y no en None: None significa "lo que
                    # diga calibracion.EXCESO", que es True, y entonces cambiar
                    # el modelo de la interfaz obligaria a cambiar el default
                    # del modulo y con el todo lo demas que lo usa.
                    exceso=EXCESO,
                    forzar_simulacion=bool(self.v_forzar.get()))

    # -------------------------------------------------------------- correr --
    def _calibrar(self):
        if self.hilo is not None and self.hilo.is_alive():
            return
        try:
            corrida = self._corrida_actual()
            opciones = self._opciones()
        except ValueError as e:
            self.v_aviso.set(str(e))
            return
        self.v_aviso.set("")
        self.etapas.reiniciar()
        self.bar_global["value"] = 0
        self.txt.delete("1.0", tk.END)
        for lz in self.lienzos.values():
            lz.limpiar("calculando...")
        self.lz_resultado.limpiar("calculando...")
        self.salida = None
        self.cancelar.clear()
        self.b_calibrar.config(state="disabled")
        self.b_cancelar.config(state="normal")
        self.v_estado.set("preparando...")
        self.hilo = threading.Thread(target=self._trabajar,
                                     args=(corrida, opciones), daemon=True)
        self.hilo.start()

    def _traza(self, texto):
        """Una linea por etapa en `.cache/gui_traza.log` (ver `trazar`)."""
        trazar(texto)

    # ------------------------------------------------ ventana a la traza --
    def _dpi_monitor(self):
        """DPI del monitor donde esta la ventana AHORA (Windows 10+), o None."""
        if sys.platform != "win32":
            return None
        try:
            import ctypes
            return int(ctypes.windll.user32.GetDpiForWindow(
                int(self.wm_frame(), 16)))
        except Exception:                                   # noqa: BLE001
            return None

    def _trazar_entorno(self):
        """Pantalla, ventana y escalado, una vez que la ventana se acomodo."""
        try:
            partes = ["pantalla %dx%d" % (self.winfo_screenwidth(),
                                          self.winfo_screenheight())]
            if sys.platform == "win32":
                import ctypes
                partes.append("%d monitor(es)"
                              % ctypes.windll.user32.GetSystemMetrics(80))
            partes.append("tk scaling %.3f (%.1f px/in)"
                          % (float(self.tk.call("tk", "scaling")),
                             self.winfo_fpixels("1i")))
        except Exception:                                   # noqa: BLE001
            return
        self._traza("entorno: " + ", ".join(partes))

    def _al_configurar_ventana(self, ev):
        # `bind` sobre la raiz recibe el <Configure> de TODOS los widgets (la
        # raiz esta en las bindtags de cada uno): solo interesa el propio.
        if ev.widget is not self:
            return
        if self._cfg_id is not None:
            self.after_cancel(self._cfg_id)
        self._cfg_id = self.after(600, self._trazar_ventana)

    def _trazar_ventana(self):
        self._cfg_id = None
        try:
            ahora = (self.winfo_width(), self.winfo_height(), self.state(),
                     self.winfo_rootx(), self.winfo_rooty(), self._dpi_monitor())
        except tk.TclError:
            return
        if ahora == self._ultima_ventana:
            return
        self._ultima_ventana = ahora
        self._traza("ventana: %dx%d en (%d, %d), %s, dpi del monitor %s"
                    % (ahora[0], ahora[1], ahora[3], ahora[4], ahora[2],
                       ahora[5] if ahora[5] is not None else "?"))

    def _sink(self, m):
        """Lo que manda el pipeline: va a la cola de la GUI y a la traza."""
        if m.get("tipo") == "etapa":
            self._traza("etapa: %s" % m.get("clave", "?"))
        elif m.get("tipo") == "fin_etapa":
            self._traza("  fin %s (%.1f s) %s" % (m.get("clave", "?"),
                                                  m.get("segundos", 0.0),
                                                  m.get("nota", "")))
        self.cola.put(m)

    def _trabajar(self, corrida, opciones):
        """Corre en el HILO TRABAJADOR. Solo manda mensajes; no toca widgets."""
        try:
            self._traza("--- calibrar %s | muestra %s | opciones %s"
                        % (corrida.clave, os.path.basename(corrida.muestra or "-"),
                           opciones))
            pr = pipeline.Progreso(sink=self._sink, cancelar=self.cancelar)
            pr.pesos = pipeline.repartir(corrida, opciones)
            self.cola.put(dict(tipo="pesos", pesos=pr.pesos))
            salida = pipeline.calibrar_completo(corrida, opciones, pr)
            self._traza("--- terminado, D = %.4f"
                        % salida["resultado"]["D"])
            self.cola.put(dict(tipo="listo", salida=salida, corrida=corrida))
        except pipeline.Cancelado:
            self._traza("--- cancelado")
            self.cola.put(dict(tipo="cancelado"))
        except Exception as e:
            self._traza("--- ERROR: %s: %s" % (type(e).__name__, e))
            self.cola.put(dict(tipo="error", mensaje=str(e),
                               traza=traceback.format_exc()))

    def _cancelar(self):
        self.cancelar.set()
        self.v_estado.set("cancelando...")

    # ------------------------------------------------------------- la cola --
    def _vaciar_cola(self):
        try:
            while True:
                self._mensaje(self.cola.get_nowait())
        except queue.Empty:
            pass
        self._vaciar_id = self.after(80, self._vaciar_cola)

    def destroy(self):
        # El `after` periodico sobrevive a la ventana: si el proceso sigue vivo
        # (los verificadores abren mas de una App), Tk lo dispara contra un
        # comando que ya no existe y lo imprime como error.
        for atr in ("_vaciar_id", "_cfg_id"):
            ident = getattr(self, atr, None)
            if ident is not None:
                try:
                    self.after_cancel(ident)
                except tk.TclError:
                    pass
        super().destroy()

    def _mensaje(self, m):
        t = m["tipo"]
        if t == "pesos":
            self._escribir("pesos de la barra: %s\n" % ", ".join(
                "%s %.0f%%" % (k, 100 * v)
                for k, v in sorted(m["pesos"].items(), key=lambda x: -x[1])))
        elif t == "etapa":
            # se recuerda cual es la etapa viva: los avisos de avance no
            # repiten la clave, y sin esto la barra de una etapa se actualizaria
            # con el avance de otra.
            self._etapa_actual = m["clave"]
            self.etapas.comenzar(m["clave"])
            self.v_estado.set(m["texto"])
        elif t == "avance":
            self.etapas.avance(self._etapa_viva(), m["frac"], m.get("texto", ""))
            if m.get("global_") is not None:
                self.bar_global["value"] = 100.0 * m["global_"]
            eta = m.get("eta")
            self.v_estado.set("%s%s" % (
                m.get("texto", ""),
                "   ·   faltan %s" % _tiempo(eta) if eta else ""))
        elif t == "fin_etapa":
            nota = m.get("nota", "")
            estado = "cache" if "cache" in nota else "hecho"
            self.etapas.terminar(m["clave"], m["segundos"], nota, estado)
            self._etapa_actual = None
        elif t == "log":
            self._escribir(m["texto"] + "\n")
        elif t == "listo":
            self._terminado(m["salida"], m["corrida"])
        elif t == "cancelado":
            self._fin_corrida("cancelado por el usuario")
            if getattr(self, "_etapa_actual", None):
                self.etapas.marcar(self._etapa_actual, "cancelado", "cancelado")
        elif t == "error":
            self._fin_corrida("ERROR: %s" % m["mensaje"].splitlines()[0][:60])
            if getattr(self, "_etapa_actual", None):
                self.etapas.marcar(self._etapa_actual, "error", "falló aquí")
            self._escribir("\n" + "=" * 70 + "\nERROR\n" + "=" * 70 + "\n"
                           + m["mensaje"] + "\n\n" + m["traza"])

    def _etapa_viva(self):
        return getattr(self, "_etapa_actual", None)

    def _escribir(self, s):
        self.txt.insert(tk.END, s)
        self.txt.see(tk.END)

    def _fin_corrida(self, estado):
        self.b_calibrar.config(state="normal")
        self.b_cancelar.config(state="disabled")
        self.v_estado.set(estado)

    # --------------------------------------------------------- el resultado --
    def _terminado(self, salida, corrida):
        self.salida = salida
        r = salida["resultado"]
        self.bar_global["value"] = 100.0
        self._fin_corrida("listo · D = %.4f ± %.4f mV·ns/%s"
                          % (r["D"], r["err_D"], _u(r)))
        self._escribir("\n" + texto_de(calibracion.informe, r))
        registro.guardar(r, corrida, detector=corrida.detector)
        self._refrescar_registro()
        # Las figuras se PROGRAMAN, no se dibujan todas aca. La que esta al
        # frente se hace en el acto; el resto, al abrir su pestana. Ver el
        # docstring de `Lienzo`: las pestanas de atras miden 1x1 px, o sea que
        # dibujarlas no le mostraba nada a nadie y ademas metia siete blits de
        # Tk por corrida en el camino donde el proceso se cayo sin mensaje.
        self._programar(self.lz_resultado, "resultado", self._fig_resultado)
        for clave, metodo in (("ventana", self._fig_ventana),
                              ("espectro", self._fig_espectro),
                              ("ncc", self._fig_ncc),
                              ("sa", self._fig_recocido),
                              ("comparacion", self._fig_comparacion),
                              ("parametros", self._fig_parametros),
                              ("paisaje", self._fig_paisaje)):
            self._programar(self.lienzos[clave], clave, metodo)

    def _programar(self, lz, clave, metodo):
        """Deja la figura lista y anota en la traza cuando se dibuja."""
        salida = self.salida
        lz.nombre = clave

        def dibujo(fig):
            # Una linea ANTES y otra DESPUES: si el proceso se cae dibujando,
            # la ultima linea de la traza dice exactamente cual figura estaba
            # en curso. Sin esto, un cierre despues de calcular parecia -y asi
            # se reporto- un problema del recocido.
            self._traza("dibujando %s" % clave)
            metodo(fig, salida)
            self._traza("  dibujado %s" % clave)

        lz.programar(dibujo)

    def _al_cambiar_pestana(self, _e=None):
        """Dibuja lo que quedo pendiente de la pestana que se acaba de abrir."""
        try:
            titulo = self.nb.tab(self.nb.select(), "text")
        except tk.TclError:
            return
        self._traza("pestana: %s" % titulo)
        # `asegurar` no hace nada en los lienzos que no se ven, asi que
        # llamarlos a todos es mas simple -y mas robusto- que averiguar cual
        # widget corresponde a la solapa seleccionada.
        for lz in list(self.lienzos.values()) + [self.lz_resultado,
                                                 self.lz_registro]:
            lz.asegurar()
        if titulo == "Registro":
            self._refrescar_registro()

    # ------------------------------------------------------------- figuras --
    def _fig_resultado(self, fig, s):
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        r, e, p = s["resultado"], s["espectro"], s["problema"]
        u = _u(r)
        ax = fig.add_subplot(121)
        ax.step(e["canales"], e["neto"], where="mid", lw=1, color="0.55",
                label="espectro neto")
        ax.plot(p.canales, p.prediccion(calibracion.parametros(r)), lw=1.6,
                color="C1", label="modelo ajustado")
        tope = tope_sobre_el_exceso(e, r)
        for k, (E, ch, ech) in enumerate(zip(r["bordes_E"], r["bordes_ch"],
                                             r["err_bordes_ch"])):
            ax.axvline(ch, color="C3", ls="--", lw=1.2)
            ax.axvspan(ch - ech, ch + ech, color="C3", alpha=.25)
            # escalonadas: con dos bordes a ~190 mV*ns las etiquetas se pisan
            ax.annotate("%.0f mV·ns\n(%.0f %s)" % (ch, E, u), xy=(ch, 0),
                        xytext=(ch + 22, (0.88 - 0.24 * k) * tope),
                        fontsize=8, color="C3")
        ax.set_xlim(0, 1.6 * max(r["bordes_ch"]))
        ax.set_title("donde caen los bordes", fontsize=10)
        ax.set_xlabel("area [mV·ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

        ax = fig.add_subplot(122)
        En = e["canales"] / r["D"]
        ax.step(En, e["neto"] * r["D"], where="mid", lw=1, color="0.55",
                label="medido, calibrado")
        ax.plot(p.canales / r["D"],
                p.prediccion(calibracion.parametros(r)) * r["D"], lw=1.6,
                color="C1", label="modelo")
        for E in r["bordes_E"]:
            ax.axvline(E, color="C3", ls="--", lw=1.2)
        y0 = 0.35 * tope_sobre_el_exceso(e, r) * r["D"]
        ax.errorbar([r["E_ref"]], [y0], xerr=[[r["sigma_E"]], [r["sigma_E"]]],
                    fmt="o", color="C0", capsize=4,
                    label="±1σ = %.0f %s" % (r["sigma_E"], u))
        ax.set_xlim(0, 1.6 * max(r["bordes_E"]))
        ax.set_title("EL ESPECTRO CALIBRADO", fontsize=10)
        ax.set_xlabel("energia [%s]" % u); ax.set_ylabel("cuentas por %s" % u)
        ax.legend(fontsize=8); ax.grid(alpha=.3)

    def _fig_ventana(self, fig, s):
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        v = s["ventana"]
        tau, forma = v["tau_rel"], v["forma"]
        _, frac, info = pulso.curva_convergencia(tau, forma)
        ax = fig.add_subplot(121)
        ax.plot(tau, forma, lw=1.4, color="C0")
        ax.axhline(0, color="k", lw=.6)
        ax.axvline(-v["pre"], color="C2", ls="--", label="t_ini = %.1f ns" % -v["pre"])
        ax.axvline(v["post"], color="C3", ls="--", label="t_fin = %.1f ns" % v["post"])
        ax.set_title("pulso promedio alineado (%d pulsos)" % v["n_usados"],
                     fontsize=10)
        ax.set_xlabel("t − onset [ns]"); ax.set_ylabel("forma normalizada")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

        ax = fig.add_subplot(122)
        ax.plot(tau, 100 * frac, lw=1.6, color="C0")
        ax.axhline(100, color="k", ls=":", lw=1)
        ax.axvline(-v["pre"], color="C2", ls="--")
        ax.axvline(v["post"], color="C3", ls="--")
        ax.set_ylim(-8, 118)
        ax.set_title("integral acumulada: recoge %.2f%% de la carga"
                     % (100 * v["fraccion"]), fontsize=10)
        ax.set_xlabel("t − onset [ns]"); ax.set_ylabel("carga acumulada [%]")
        ax.grid(alpha=.3)

    def _fig_espectro(self, fig, s):
        """
        Pagina 1: lo que ENTRA al ajuste, medido y simulado, en una vista.

        Arriba el dato: muestra y fondo con el binning comun (en log, porque el
        fondo esta uno o dos ordenes por debajo) y el neto con sus barras de
        error, con la ventana de ajuste marcada.

        Abajo la plantilla de Geant4: primero en su propio eje de energia, que
        es como sale del `.root`, y despues superpuesta al dato usando la D
        ajustada para llevar el eje medido a energia. Esa superposicion es el
        control mas directo de que la plantilla corresponde a esta fuente: si
        los bordes no caen juntos, no hay ajuste que lo arregle (es la falla
        silenciosa de usar el .root de otra fuente).
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        e, r, p = s["espectro"], s["resultado"], s["problema"]
        cen, pesos = s["plantilla"]
        u = _u(r)

        # ---- (1) muestra y fondo, mismo eje ------------------------------
        ax = fig.add_subplot(221)
        ax.step(e["canales"], e["muestra"], where="mid", lw=1, label="muestra")
        if e["hay_fondo"]:
            ax.step(e["canales"], e["fondo"], where="mid", lw=1,
                    label="fondo · k = %.4f" % e["k_f"])
        else:
            ax.plot([], [], " ", label="sin fondo: el neto es el crudo")
        ax.set_yscale("log")
        ax.set_title("el dato crudo: %s pulsos\n%.0f s vivos de %.0f de pared"
                     % (format(e["n_muestra"], ",d"), e["t_vivo_muestra"],
                        e["t_pared_muestra"]), fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("cuentas")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

        # ---- (2) el neto, con la ventana de ajuste -----------------------
        ax = fig.add_subplot(222)
        ax.errorbar(e["canales"], e["neto"], yerr=e["error"], fmt="none",
                    ecolor="0.8", lw=.5)
        ax.step(e["canales"], e["neto"], where="mid", lw=1.2, color="C2",
                label="neto")
        ax.axvspan(p.lo, p.hi, color="C0", alpha=.12,
                   label="ventana de ajuste [%.0f, %.0f]" % (p.lo, p.hi))
        ax.axhline(0, color="k", lw=.8)
        for ch in r["bordes_ch"]:
            ax.axvline(ch, color="C3", ls="--", lw=1)
        ax.set_xlim(0, min(e["canales"][-1], 1.9 * max(r["bordes_ch"])))
        ax.set_title("espectro neto y donde se ajusta\n%d bines de %.1f mV·ns"
                     % (e["bins"], e["rango_max"] / e["bins"]), fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("cuentas")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

        # ---- (3) la plantilla, en energia --------------------------------
        ax = fig.add_subplot(223)
        paso_E = float(np.median(np.diff(cen)))
        ax.step(cen, pesos / max(pesos.max(), 1.0), where="mid", lw=1.2,
                color="C1", label="plantilla %s" % r["rama"])
        for E in r["bordes_E"]:
            ax.axvline(E, color="C3", ls="--", lw=1)
            ax.annotate("%.0f %s" % (E, u), xy=(E, 1.18), xytext=(-3, 0),
                        textcoords="offset points", rotation=90, fontsize=7.5,
                        color="C3", ha="right", va="top")
        lineas = ", ".join("%.0f keV" % e0 for e0, _ in s.get("lineas_sim", []))
        ax.set_ylim(0, 1.22)
        # tres lineas y no dos: en una sola, los 60 caracteres de la segunda
        # son mas anchos que el panel, y un titulo que no entra no se acorta
        # ni se parte solo, se sale de los ejes y se recorta contra el borde
        # de la figura (se veia "2,079 eventos" en vez de "852,079 eventos")
        ax.set_title("la simulacion tal como sale del .root\n"
                     "%s eventos, bines de %.1f keV\ngammas primarios %s"
                     % (format(int(pesos.sum()), ",d"), paso_E, lineas or "?"),
                     fontsize=9)
        ax.set_xlabel("energia [%s]" % u)
        ax.set_ylabel("cuentas (norm.)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

        # ---- (4) los dos juntos, en energia ------------------------------
        # La plantilla se dibuja con la resolucion ajustada y tambien casi sin
        # borronear (un bin, que es lo mas angosto que se puede dibujar sin que
        # el remuestreo la vuelva un peine de picos). La diferencia entre las
        # dos es lo que hace el detector; la que queda a baja energia entre la
        # borroneada y el dato es el exceso de bajo canal, que la plantilla no
        # tiene y el modelo agrega aparte.
        ax = fig.add_subplot(224)
        En = e["canales"] / r["D"]
        paso_ch = float(np.median(np.diff(e["canales"])))
        pl = mod_ncc.plantilla_en_canal(e["canales"], cen, pesos, r["D"],
                                        r["sigma_rel"], r["E_ref"])
        pl_fina = mod_ncc.plantilla_en_canal(
            e["canales"], cen, pesos, r["D"],
            paso_ch / (r["D"] * r["E_ref"]), r["E_ref"])
        ax.step(En, e["neto"], where="mid", lw=1, color="0.45", label="medido")
        ax.plot(En, r["A"] * pl_fina, lw=1, ls=":", color="C0",
                label="plantilla sin resolucion")
        ax.plot(En, r["A"] * pl, lw=1.6, color="C1",
                label="plantilla ⊗ σ = %.0f %s" % (r["sigma_E"], u))
        ax.axvspan(p.lo / r["D"], p.hi / r["D"], color="C0", alpha=.10)
        for E in r["bordes_E"]:
            ax.axvline(E, color="C3", ls="--", lw=1)
        ax.set_xlim(0, 1.6 * max(r["bordes_E"]))
        dentro = e["neto"][e["canales"] >= p.lo]
        if dentro.size:
            tope = float(np.max(dentro))
            ax.set_ylim(-0.05 * tope, 1.45 * tope)
        ax.set_title("dato y simulacion en el MISMO eje\nD = %.4f mV·ns/%s"
                     % (r["D"], u), fontsize=10)
        ax.set_xlabel("energia [%s]" % u)
        ax.set_ylabel("cuentas")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

    def _fig_ncc(self, fig, s):
        """
        Pagina 2: como se eligio el intervalo, y cuanto vale la correlacion.

        El NCC barre UN parametro -el canal del borde B, que fija la ventana
        entera porque el mapa pasa por el origen- y para cada valor calcula el
        Pearson entre la forma del dato y la de la plantilla. Los cuatro
        paneles son: la curva ρ(B) del barrido, la ventana que gano dibujada
        sobre el espectro, las DOS curvas que Pearson compara de verdad (ya
        estandarizadas: Pearson resta la media y divide por la desviacion, o
        sea que no ve ni offset ni escala) y el diagrama de dispersion de una
        contra la otra, que es la lectura geometrica de ρ: cuanto mas se
        parece a una recta, mas cerca de 1.
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        r, e, n, p = s["resultado"], s["espectro"], s["ncc"], s["problema"]
        cen, pesos = s["plantilla"]

        # las dos curvas del barrido ganador, reconstruidas igual que en
        # `ncc.buscar_intervalo`: plantilla en canal con la D de ese B, dato
        # suavizado, y las dos remuestreadas a los mismos n puntos.
        pl = mod_ncc.plantilla_en_canal(e["canales"], cen, pesos, n["D"],
                                        n["sigma_rel"], n["E_ref"])
        ys = mod_ncc.suavizar(e["neto"], 21)
        m = (e["canales"] >= n["lo"]) & (e["canales"] <= n["hi"])
        u = np.linspace(n["lo"], n["hi"], n["n_pts"])
        s1 = np.interp(u, e["canales"][m], ys[m])
        s2 = np.interp(u, e["canales"][m], pl[m])
        rho_rec = mod_ncc.pearson(s1, s2)

        def z(v):
            d = float(np.std(v))
            return (v - np.mean(v)) / (d if d > 0 else 1.0)

        # ---- (1) el barrido ----------------------------------------------
        ax = fig.add_subplot(221)
        ax.plot(n["Bs"], n["rhos"], lw=1.4, color="C0")
        ax.axvline(n["B"], color="C3", ls="--",
                   label="B = %.0f mV·ns  ·  ρ = %.5f" % (n["B"], n["rho"]))
        ax.plot([n["B"]], [n["rho"]], "o", color="C3", ms=6)
        ax.set_ylim(max(0.0, float(n["rhos"].max()) - 0.35), 1.005)
        ax.set_title("barrido de UN parametro:\nρ contra el canal del borde",
                     fontsize=10)
        ax.set_xlabel("canal del borde B [mV·ns]")
        ax.set_ylabel("ρ (Pearson)")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=.3)

        # cuanto castiga ρ equivocarse en B, en el propio panel: la curva sola
        # no permite leer si el maximo es agudo o una meseta.
        filas = []
        for f in (0.90, 0.95, 1.00, 1.05, 1.10):
            B = f * n["B"]
            filas.append("%+5.0f%%  B=%7.0f  D=%.4f  ρ=%.5f"
                         % (100 * (f - 1), B, B / n["E_ref"],
                            float(np.interp(B, n["Bs"], n["rhos"]))))
        ax.text(0.02, 0.03, "ρ alrededor del maximo\n" + "\n".join(filas),
                transform=ax.transAxes, fontsize=7, family="monospace",
                va="bottom", bbox=dict(fc="white", ec="#dde3e8", alpha=.9))

        # ---- (2) la ventana sobre el espectro ----------------------------
        ax = fig.add_subplot(222)
        ax.step(e["canales"], e["neto"], where="mid", lw=1, color="0.55",
                label="neto")
        esc = np.sum(e["neto"][m]) / max(np.sum(pl[m]), 1e-9)
        ax.plot(e["canales"], pl * esc, lw=1.3, color="C1",
                label="plantilla (σ_rel = %.2f del barrido)" % n["sigma_rel"])
        ax.axvspan(n["lo"], n["hi"], color="C0", alpha=.13,
                   label="ventana NCC [%.0f, %.0f]" % (n["lo"], n["hi"]))
        if p.lo < n["lo"] - 1:
            # el ajuste arranca mas abajo que el NCC a proposito: el exceso de
            # bajo canal esta modelado, asi que conviene ajustar el flanco
            # entero en vez de esquivarlo (ver modelo.armar)
            ax.axvspan(p.lo, n["lo"], color="C2", alpha=.13,
                       label="extension del ajuste (exceso modelado)")
        ax.axvline(n["B"], color="C3", ls="--", lw=1.2, label="B")
        ax.set_xlim(0, min(e["canales"][-1], 1.9 * n["B"]))
        ax.set_title("el intervalo elegido, sin intervencion", fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("cuentas")
        ax.legend(fontsize=7.5)
        ax.grid(alpha=.3)

        # ---- (3) lo que Pearson compara ----------------------------------
        ax = fig.add_subplot(223)
        ax.plot(u, z(s1), lw=1.4, color="C0", label="dato suavizado (S1)")
        ax.plot(u, z(s2), lw=1.4, color="C1", label="plantilla (S2)")
        ax.set_title("las dos formas estandarizadas\n%d puntos, uno por bin"
                     % n["n_pts"], fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("z-score")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

        # ---- (4) ρ, geometricamente --------------------------------------
        ax = fig.add_subplot(224)
        z1, z2 = z(s1), z(s2)
        ax.scatter(z2, z1, s=10, alpha=.6, color="C0")
        lim = [min(z1.min(), z2.min()) - 0.2, max(z1.max(), z2.max()) + 0.2]
        ax.plot(lim, lim, ls="--", lw=1, color="C3", label="ρ = 1 seria esto")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_title("ρ, geometricamente: %.5f\n(recalculado aqui: %.5f)"
                     % (n["rho"], rho_rec), fontsize=10)
        ax.set_xlabel("plantilla, z-score")
        ax.set_ylabel("dato, z-score")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

    # nombre y unidad de cada parametro LIBRE que puede aparecer en `p.libres`
    # (D y sigma_rel siempre estan; lam y c0 son condicionales al modelo)
    _NOMBRE_LIBRE = {"D": "D  (ganancia)", "sigma_rel": "σ_rel  (resolucion)",
                     "lam": "λ_exceso  (bajo canal)",
                     "c0": "c0  (ordenada, test del origen)"}

    def _fig_parametros(self, fig, s):
        """
        El vector Theta completo del ajuste: un parametro por fila, con su
        valor, su error (de la curvatura del chi2, ya escalado por
        sqrt(chi2/ndf)) y cuanto pesa ese error en relativo. Al lado, la matriz
        de CORRELACION entre esos parametros: un error grande en D que viene
        acompanado de una correlacion fuerte con sigma (tipico aca, ver
        calibracion.informe) no es lo mismo que un error grande aislado, y la
        tabla sola no lo muestra.

        `A` y `C_exceso` (las amplitudes) se listan aparte: no son parametros
        que el recocido busque, se PROYECTAN exactas (VARPRO) en cada
        evaluacion, asi que no tienen curvatura ni entran en la matriz.

        Las unidades y las aclaraciones van pegadas al nombre en vez de en una
        columna propia: cinco columnas de texto no entran en medio panel y
        matplotlib no las achica, las RECORTA sin avisar.
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        r, p = s["resultado"], s["problema"]
        u = _u(r)

        libres = r["libres"]
        x = r["x_libres"]
        err = r["err_libres"]
        cov = np.asarray(r["cov_libres"], float)
        unidad = {"D": " [mV·ns/%s]" % u, "sigma_rel": "  (adim.)",
                  "lam": " [mV·ns]", "c0": " [mV·ns]"}

        filas = []
        for nombre, v, e in zip(libres, x, err):
            rel = "%.2f %%" % (100 * abs(e / v)) if v else "-"
            filas.append((self._CORTO.get(nombre, nombre)
                          + unidad.get(nombre, ""),
                          "%.5g" % v, "± %.3g" % e, rel))
        filas.append(("", "", "", ""))
        filas.append(("A [cuentas]", "%.5g" % r["A"], "proyectada", "VARPRO"))
        if r.get("exceso"):
            filas.append(("C_exceso [cuentas]", "%.5g" % r["C_exceso"],
                          "proyectada", "nnls"))
        filas.append(("", "", "", ""))
        # las tres de abajo no son parametros: son los indicadores de si el
        # resultado es de fiar (describe el modelo a los datos, es reproducible
        # el optimo, correlaciona la ventana automatica)
        filas.append(("χ²/ndf  (%d gdl)" % r["ndf"], "%.3f" % r["chi2_ndf"],
                      "", "bondad"))
        filas.append(("dispersion SA (%d arr.)" % r["n_arranques"],
                      "%.3g" % r["dispersion_multiarranque"], "",
                      "%.3f %%" % (100 * r["dispersion_multiarranque"]
                                   / r["D"])))
        filas.append(("ρ NCC (ventana auto)", "%.5f" % r["rho"], "", ""))

        gs = fig.add_gridspec(1, 2, width_ratios=(1.45, 1.0))
        ax = fig.add_subplot(gs[0, 0])
        nombres = ", ".join(
            self._NOMBRE_LIBRE.get(k, k).replace("  (", " = ").rstrip(")")
            for k in libres)
        tabla(ax, filas, ("parametro", "valor", "error", "rel."),
              (0.44, 0.20, 0.22, 0.14), titulo="vector Θ del ajuste",
              nota=nombres + ".")

        ax = fig.add_subplot(gs[0, 1])
        n = len(libres)
        diag = np.sqrt(np.maximum(np.diag(cov), 1e-300))
        corr = cov / np.outer(diag, diag)
        im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
        etiquetas = [self._CORTO.get(k, k) for k in libres]
        ax.set_xticks(range(n)); ax.set_xticklabels(etiquetas, fontsize=8)
        ax.set_yticks(range(n)); ax.set_yticklabels(etiquetas, fontsize=8)
        for i in range(n):
            for j in range(n):
                ax.text(j, i, "%+.2f" % corr[i, j], ha="center", va="center",
                        fontsize=8.5,
                        color="white" if abs(corr[i, j]) > 0.6 else "black")
        ax.set_title("correlacion entre\nparametros libres", fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # --------------------------------------------------- el recocido (3) --
    # Etiqueta corta de cada parametro libre, para tablas y ejes.
    _CORTO = {"D": "D", "sigma_rel": "σ_rel", "lam": "λ_exc", "c0": "c0"}

    def _fig_recocido(self, fig, s):
        """
        Pagina 3: de donde arranca el recocido, por donde pasa y donde termina.

        Cuatro paneles, que son las cuatro preguntas que hay que hacerle a un
        optimizador estocastico:

        (1) DE DONDE. El NCC entrega la semilla (D = B/E_ref y el σ_rel del
            barrido) y con ella se arman las COTAS de la busqueda. Ojo con un
            detalle que la tabla dice explicitamente: los N arranques NO salen
            de la semilla, salen de puntos al azar dentro de las cotas
            (`sa.multiarranque` con x0=None). Eso es a proposito: si todos
            arrancaran del mismo punto, que todos terminen igual no probaria
            nada sobre el minimo global.
        (2) COMO BAJA. El costo de cada arranque ciclo a ciclo. La linea gruesa
            es el ganador; las finas, los otros. Que todas terminen apiladas es
            la evidencia practica de que el minimo es global -y la dispersion
            entre ellas es un numero que se reporta.
        (3) LA MECANICA. Temperatura y tasa de aceptacion del ganador. Al
            principio se acepta casi todo (la cadena pasea), al final casi
            nada (se comporta como un optimizador local). Si la aceptacion
            arranca en ~1 o en ~0, T_init esta mal calibrada.
        (4) POR DONDE. El camino en el plano (D, σ_rel). Esto es lo que
            distingue "bajo en linea recta" de "cruzo de un valle a otro", que
            es la unica razon por la que se usa recocido y no un optimizador
            local.
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        r, det = s["resultado"], s.get("sa") or {}
        corridas = det.get("corridas")
        if not corridas:
            raise ValueError("esta corrida no guardo el detalle del recocido")
        p = s["problema"]
        libres = det.get("libres", r["libres"])
        semilla = np.asarray(det["semilla"], float)
        cotas = np.asarray(det["cotas"], float)
        mejor = int(det["mejor"])
        u = _u(r)
        ndf = float(p.ndf)

        # ---- (1) la semilla, las cotas y los arranques -------------------
        # La columna izquierda va mas ancha que la derecha: de este lado hay
        # tablas de 4 columnas de texto, y con la mitad justa del ancho los
        # rotulos se cortan (matplotlib no los achica, los recorta).
        gs = fig.add_gridspec(2, 2, width_ratios=(1.35, 1.0))
        ax = fig.add_subplot(gs[0, 0])
        filas = []
        for k, nombre in enumerate(libres):
            filas.append((self._CORTO.get(nombre, nombre),
                          "%.5g" % semilla[k],
                          "[%.4g, %.4g]" % (cotas[k, 0], cotas[k, 1]),
                          "%.5g" % r["x_libres"][k]))
        filas.append(("χ²/ndf", "%.3f" % (det["chi2_semilla"] / ndf), "",
                      "%.3f" % r["chi2_ndf"]))
        filas.append(("", "", "", ""))
        for i, m in enumerate(corridas):
            corte = "estable" if m["motivo"].startswith("estab") else "T min"
            filas.append(("#%d" % (i + 1),
                          "D₀ = %.4f" % m["x0"][0],
                          "%d cic · %s" % (m["ciclos"], corte),
                          "D = %.4f" % m["x_pulido"][0]))
        tabla(ax, filas,
              ("parametro", "INICIAL", "cotas · corte", "FINAL"),
              (0.17, 0.26, 0.30, 0.27),
              titulo="la semilla, las cotas y los %d arranques" % len(corridas),
              resaltar=(len(libres) + 2 + mejor,),
              nota="arrancan AL AZAR dentro de las cotas, no en la semilla: "
                   "por eso que todos terminen igual dice algo.")

        # ---- (2) el costo, ciclo a ciclo ---------------------------------
        ax = fig.add_subplot(gs[0, 1])
        for i, m in enumerate(corridas):
            h = m["historia"]
            if i == mejor:
                continue
            ax.semilogy(h[:, 2] / ndf, lw=.8, color="0.72")
        h = corridas[mejor]["historia"]
        # Los limites se ponen a mano: `costo` devuelve 1e30 cuando un
        # parametro se sale de lo fisico (sigma o lambda negativos), y un solo
        # valor asi en escala log aplasta todo lo demas contra el piso.
        vistos = [c["historia"][:, 2] for c in corridas] + [h[:, 1]]
        vistos = np.concatenate(vistos) / ndf
        vistos = vistos[np.isfinite(vistos) & (vistos < 1e10)]
        if vistos.size:
            ax.set_ylim(0.92 * r["chi2_ndf"], 1.25 * float(vistos.max()))
        ax.semilogy(h[:, 1] / ndf, lw=.9, color="C0", alpha=.7,
                    label="ganador: costo actual")
        ax.semilogy(h[:, 2] / ndf, lw=1.8, color="C3",
                    label="ganador: mejor hasta aqui")
        ax.axhline(det["chi2_semilla"] / ndf, color="C1", ls=":", lw=1.2,
                   label="semilla del NCC (%.3f)" % (det["chi2_semilla"] / ndf))
        ax.axhline(r["chi2_ndf"], color="k", ls="--", lw=1,
                   label="final, ya pulido (%.3f)" % r["chi2_ndf"])
        ax.set_title("como baja el costo: %d arranques" % len(corridas),
                     fontsize=10)
        ax.set_xlabel("ciclo de temperatura")
        ax.set_ylabel("χ²/ndf")
        ax.legend(fontsize=7.5)
        ax.grid(alpha=.3, which="both")

        # ---- (3) temperatura y aceptacion --------------------------------
        ax = fig.add_subplot(gs[1, 0])
        ax.plot(h[:, 3], lw=1.4, color="C0", label="fraccion aceptada")
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("ciclo de temperatura")
        ax.set_ylabel("fraccion aceptada", color="C0")
        a2 = ax.twinx()
        a2.semilogy(h[:, 0], lw=1.3, color="C1")
        a2.set_ylabel("temperatura (unidades de χ²)", color="C1")
        ax.set_title("la mecanica del enfriamiento\nT_init = %.3g, calibrada"
                     % corridas[mejor]["T_init"], fontsize=10)
        ax.legend(fontsize=8, loc="center right")
        ax.grid(alpha=.3)

        # ---- (4) el camino en el plano (D, sigma) ------------------------
        ax = fig.add_subplot(gs[1, 1])
        for i, m in enumerate(corridas):
            hx = m["historia_x"]
            if i == mejor or hx.size == 0:
                continue
            ax.plot(hx[:, 0], hx[:, 1], lw=.6, color="0.72", alpha=.55,
                    zorder=1)
            ax.plot([m["x0"][0]], [m["x0"][1]], "s", ms=4, color="0.55",
                    zorder=2)
        hx = corridas[mejor]["historia_x"]
        ax.plot(hx[:, 0], hx[:, 1], lw=.9, color="0.35", zorder=3)
        sc = ax.scatter(hx[:, 0], hx[:, 1], c=np.arange(len(hx)), s=16,
                        cmap="plasma", zorder=4)
        fig.colorbar(sc, ax=ax, label="ciclo de temperatura")
        ax.plot([corridas[mejor]["x0"][0]], [corridas[mejor]["x0"][1]], "s",
                ms=8, color="C0", zorder=5, label="arranque del ganador")
        ax.plot([semilla[0]], [semilla[1]], "*", ms=13, color="C1", zorder=5,
                label="semilla del NCC")
        # el final de CADA arranque: que los 8 puntos rojos caigan uno sobre
        # otro es la version visual de la dispersion del multiarranque
        for m in corridas:
            ax.plot([m["x_pulido"][0]], [m["x_pulido"][1]], "o", ms=3.5,
                    color="#e02020", zorder=5)
        ax.plot([r["D"]], [r["sigma_rel"]], "o", ms=8, color="#e02020",
                zorder=6, label="optimo de los %d (D = %.4f)"
                                % (len(corridas), r["D"]))
        ax.set_xlim(cotas[0, 0], cotas[0, 1])
        ax.set_ylim(cotas[1, 0], cotas[1, 1])
        ax.set_title("por donde paseo la busqueda", fontsize=10)
        ax.set_xlabel("D [mV·ns/%s]" % u)
        ax.set_ylabel("σ_rel")
        ax.legend(fontsize=7.5)
        ax.grid(alpha=.3)

    # ------------------------------------------------ inicial vs final (4) --
    def _fig_comparacion(self, fig, s):
        """
        Pagina 4: el punto de partida contra el punto final, sobre los datos.

        La semilla del NCC no es un numero cualquiera: sale de correlacionar
        formas, y en este proyecto cae a ~1-2% del resultado. Lo que esta
        pagina permite ver es exactamente CUANTO agrego el ajuste sobre esa
        estimacion, y en que parte del espectro lo agrego.

        El panel de abajo a la derecha es un corte del χ² a lo largo del
        SEGMENTO que une los dos vectores. No es el paisaje (eso es la otra
        pestana): es la recta que va de uno al otro, parametrizada por t, con
        t = 0 la semilla y t = 1 el optimo. Que el minimo caiga en t = 1 y que
        la curva siga subiendo despues es la manera mas barata de ver que el
        ajuste no se quedo a medio camino ni se paso de largo.
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        r, e, p, det = s["resultado"], s["espectro"], s["problema"], s.get("sa") or {}
        if not det:
            raise ValueError("esta corrida no guardo el detalle del recocido")
        libres = det.get("libres", r["libres"])
        x0 = np.asarray(det["semilla"], float)
        x1 = np.asarray(r["x_libres"], float)
        u = _u(r)

        # las dos predicciones, dentro de la ventana de ajuste
        y0 = p.prediccion(tuple(x0))
        chi0 = float(np.sum(((p.dato - y0) / p.error) ** 2))
        y1 = p.prediccion(tuple(x1))
        chi1 = float(np.sum(((p.dato - y1) / p.error) ** 2))

        # ---- (1) los dos modelos sobre el dato ---------------------------
        gs = fig.add_gridspec(2, 2, width_ratios=(1.25, 1.0))
        ax = fig.add_subplot(gs[0, 0])
        ax.step(e["canales"], e["neto"], where="mid", lw=.9, color="0.78",
                label="neto (fuera de la ventana)")
        ax.errorbar(p.canales, p.dato, yerr=p.error, fmt="none", ecolor="0.65",
                    lw=.5)
        ax.step(p.canales, p.dato, where="mid", lw=1.1, color="0.35",
                label="dato ajustado")
        ax.plot(p.canales, y0, lw=1.3, ls="--", color="C1",
                label="semilla NCC: D = %.4f, χ²/ndf = %.2f"
                      % (x0[0], chi0 / p.ndf))
        ax.plot(p.canales, y1, lw=1.7, color="C3",
                label="final: D = %.4f, χ²/ndf = %.2f" % (x1[0], chi1 / p.ndf))
        ax.set_xlim(0, min(e["canales"][-1], 1.5 * p.hi))
        # la escala la fija lo que se AJUSTA: con el pico de bajo canal
        # adentro, las dos curvas quedan pegadas y no se distingue ninguna
        tope = float(np.max(p.dato))
        ax.set_ylim(-0.06 * tope, 1.45 * tope)
        ax.set_title("el modelo antes y despues de ajustar", fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("cuentas")
        ax.legend(fontsize=7.5)
        ax.grid(alpha=.3)

        # ---- (2) los residuos de los dos ---------------------------------
        ax = fig.add_subplot(gs[0, 1])
        ax.fill_between(p.canales, -1, 1, color="0.88", zorder=0)
        ax.axhline(0, color="k", lw=.8)
        ax.step(p.canales, (p.dato - y0) / p.error, where="mid", lw=1,
                color="C1", label="semilla (χ² = %.0f)" % chi0)
        ax.step(p.canales, (p.dato - y1) / p.error, where="mid", lw=1.2,
                color="C3", label="final (χ² = %.0f, %d gdl)" % (chi1, p.ndf))
        ax.set_title("residuos normalizados:\ndonde mejoro el ajuste",
                     fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("(dato − modelo)/error")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

        # ---- (3) los dos vectores, parametro por parametro ---------------
        ax = fig.add_subplot(gs[1, 0])
        # La unidad va pegada al nombre y el cambio relativo junto al cambio
        # en sigmas: cuatro columnas entran en el panel de una ventana chica y
        # seis no (matplotlib recorta el texto de la celda sin avisar).
        unidad = {"D": " [mV·ns/%s]" % u, "sigma_rel": "", "lam": " [mV·ns]",
                  "c0": " [mV·ns]"}
        filas = []
        for k, nombre in enumerate(libres):
            cambio = ("%+.2f %%" % (100 * (x1[k] - x0[k]) / x0[k])
                      if x0[k] else "-")
            # cuantos sigma del error final se movio: si un parametro se
            # movio menos que su propia barra, el ajuste no lo cambio de verdad
            err = r["err_libres"][k]
            if err:
                cambio += "  (%+.1fσ)" % ((x1[k] - x0[k]) / err)
            filas.append((self._CORTO.get(nombre, nombre)
                          + unidad.get(nombre, ""),
                          "%.5g" % x0[k], "%.5g" % x1[k], cambio))
        filas.append(("", "", "", ""))
        filas.append(("χ²  (%d gdl)" % p.ndf, "%.1f" % chi0, "%.1f" % chi1,
                      "%+.1f" % (chi1 - chi0)))
        filas.append(("χ²/ndf", "%.3f" % (chi0 / p.ndf),
                      "%.3f" % (chi1 / p.ndf), ""))
        filas.append(("A [cuentas]", "%.4g" % p.chi2_de(tuple(x0))[1],
                      "%.4g" % r["A"], ""))
        tabla(ax, filas, ("parametro", "semilla", "final", "cambio  (en σ)"),
              (0.34, 0.18, 0.18, 0.30),
              titulo="cuanto movio el ajuste a la semilla del NCC")

        # ---- (4) el corte del chi2 entre los dos vectores ----------------
        ax = fig.add_subplot(gs[1, 1])
        ts = np.linspace(-0.35, 1.6, 90)
        chis = np.array([p.costo(tuple(x0 + t * (x1 - x0))) for t in ts]) / p.ndf
        ax.plot(ts, chis, lw=1.5, color="C0")
        ax.plot([0], [chi0 / p.ndf], "s", ms=8, color="C1",
                label="t = 0: semilla del NCC")
        ax.plot([1], [chi1 / p.ndf], "o", ms=8, color="C3",
                label="t = 1: optimo")
        ax.axvline(1, color="C3", ls="--", lw=.9)
        ax.set_title("corte del χ² sobre el segmento\nsemilla → optimo",
                     fontsize=10)
        ax.set_xlabel("t  (0 = semilla, 1 = optimo)")
        ax.set_ylabel("χ²/ndf")
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)

    def _fig_paisaje(self, fig, s):
        """
        El paisaje del χ² en el plano (D, σ_rel), con el camino del recocido.

        Es un CORTE, no el χ² perfilado: los otros parametros libres (λ del
        exceso, y c0 si se libero) quedan congelados en el optimo, que es lo
        que `x_base` le dice a `Problema.mapa`. Pasarselo no es opcional en
        cuanto el modelo tiene mas de dos parametros -y el de por defecto los
        tiene-: sin eso el costo se evalua con λ = None y revienta.
        """
        ajustar(fig)   # antes de crear los ejes: ver `ajustar`
        p, r, det = s["problema"], s["resultado"], s.get("sa") or {}
        x = np.asarray(r["x_libres"], float)
        D, sg = r["D"], r["sigma_rel"]
        Ds = np.linspace(0.6 * D, 1.6 * D, 90)
        Ss = np.linspace(0.02, 0.30, 70)
        Z = p.mapa(Ds, Ss, x_base=x) / p.ndf
        ax = fig.add_subplot(121)
        # La escala de color se corta en 5 veces el minimo. Con el rango
        # completo el chi2 de las paredes (100 veces el del fondo) se lleva
        # todo el contraste y el valle -que es lo que hay que mirar- queda de
        # un solo color, con los contornos de 1, 2 y 3 sigma escondidos
        # adentro.
        z0 = float(np.nanmin(Z))
        im = ax.pcolormesh(Ds, Ss, np.log10(Z), shading="auto", cmap="viridis",
                           vmin=np.log10(0.97 * z0), vmax=np.log10(5.0 * z0))
        ax.contour(Ds, Ss, Z * p.ndf,
                   levels=[r["chi2"] + 2.30, r["chi2"] + 6.17, r["chi2"] + 11.8],
                   colors="w", linewidths=.9)
        # el camino del arranque ganador, encima del paisaje: es donde se ve si
        # el recocido tuvo que cruzar una loma o si bajo derecho
        corridas = det.get("corridas")
        if corridas:
            hx = corridas[int(det["mejor"])]["historia_x"]
            if hx.size:
                ax.plot(hx[:, 0], hx[:, 1], lw=.9, color="w", alpha=.75,
                        zorder=3, label="camino del recocido")
        ax.plot([D], [sg], "o", color="#e02020", ms=7, label="optimo")
        ax.plot([r["semilla_D"]], [det.get("semilla", (0, 0.09))[1]], "*",
                color="C1", ms=12, label="semilla del NCC")
        fig.colorbar(im, ax=ax, label="log10(χ²/ndf)")
        ax.set_title("el paisaje del χ²\n(corte: λ y c0 en el optimo)",
                     fontsize=10)
        ax.set_xlabel("D [mV·ns/%s]" % _u(r))
        ax.set_ylabel("σ_rel")
        # leyenda oscura: el camino se dibuja en blanco para que se vea sobre
        # el mapa, y en una leyenda blanca no se veria nada
        ax.legend(fontsize=8, facecolor="0.22", edgecolor="0.22",
                  labelcolor="w")

        ax = fig.add_subplot(122)
        ax.axhline(0, color="k", lw=.8)
        ax.fill_between(p.canales, -1, 1, color="0.85", zorder=0)
        res = (p.dato - p.prediccion(calibracion.parametros(r))) / p.error
        ax.step(p.canales, res, where="mid", lw=1, color="C0")
        ax.set_title("residuos normalizados  (χ²/ndf = %.2f, %d gdl)"
                     % (r["chi2_ndf"], r["ndf"]), fontsize=10)
        ax.set_xlabel("area [mV·ns]")
        ax.set_ylabel("(dato − modelo)/error")
        ax.grid(alpha=.3)

    # ------------------------------------------------------------ registro --
    def _refrescar_registro(self):
        reg = registro.cargar()
        self.tv.delete(*self.tv.get_children())
        for k, r in sorted(reg.items()):
            self.tv.insert("", tk.END, iid=k, values=(
                registro.punto(r), r["fuente"], _u(r),
                "%.4f" % r["D"], "± %.4f" % r["err_D"],
                "%.1f mV·ns" % (r["D"] * r["E_ref"]),
                "%.2f" % r["chi2_ndf"], r.get("guardado", "")))
        self._fig_linealidad()

    def _borrar_registro(self):
        sel = self.tv.selection()
        if not sel:
            return
        if messagebox.askyesno("Borrar", "¿Borrar %d entrada(s) del registro?"
                               % len(sel)):
            for k in sel:
                registro.borrar(k)
            self._refrescar_registro()

    # Colores de la pestana Registro. Okabe-Ito: se distinguen entre si y con
    # daltonismo, y ninguna categoria repite color dentro de un mismo panel.
    _COL_DET = {"A": "#0072B2", "B": "#E69F00", "C": "#009E73", "D": "#CC79A7"}
    _COL_HV = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"]
    _COL_FUENTE = {"Co-60": "#D55E00", "Cs-137": "#56B4E9"}
    _MARCA_FUENTE = {"Co-60": "o", "Cs-137": "s"}
    _LINEA_HV = ["-", "--", ":", "-."]

    @staticmethod
    def _det_hv(punto):
        """'A @ 1600 V' -> ('A', '1600 V'). Sin HV registrado, el HV es ''."""
        det, _, hv = punto.partition(" @ ")
        return det.strip(), hv.strip()

    def _fig_linealidad(self):
        """
        Los dos tests, uno al lado del otro: canal(E) y D(V).

        Las leyendas NO van adentro de los ejes: con cuatro detectores a tres
        voltajes eran ~20 entradas que tapaban los puntos. Cada cosa se
        codifica con UNA propiedad (color, marcador, tipo de linea), la clave
        va en un cuadro a la derecha y los numeros en tablas debajo de ella.
        El selector de arriba elige un detector o la comparacion.
        """
        grupos = registro.por_punto(escala_pedida="keVee")
        conjuntos = [(p, e) for p, e in grupos.items()
                     if registro.linealidad(e) is not None]
        curvas = {k: v for k, v in voltaje.puntos_del_registro().items()
                  if v["V"].size >= 2}

        dets = sorted({self._det_hv(p)[0] for p, _ in conjuntos}
                      | {d for d, _ in curvas})
        opciones = ([self.COMPARAR]
                    + ([self.GANANCIA] if curvas else [])
                    + ["Detector %s" % d for d in dets])
        self.cb_ver_reg["values"] = opciones
        if self.v_ver_reg.get() not in opciones:
            self.v_ver_reg.set(self.COMPARAR)
        ver = self.v_ver_reg.get()
        if ver == self.GANANCIA:
            self._fig_ganancia(curvas)
            return
        solo = None if ver == self.COMPARAR else ver.replace("Detector ", "")

        if solo is not None:
            conjuntos = [(p, e) for p, e in conjuntos
                         if self._det_hv(p)[0] == solo]
            curvas = {k: v for k, v in curvas.items() if k[0] == solo}

        if not conjuntos and not curvas:
            self.lz_registro.vacio(
                "Todavia no hay ningun test posible%s.\n\n"
                "TEST A (el eje, canal = D*E): necesita DOS fuentes\n"
                "medidas en el MISMO punto de operacion. Los dos\n"
                "bordes del Co-60 no sirven: el mapa por el origen los\n"
                "ata a la misma D (Ec. 12 del paper).\n\n"
                "TEST B (la ganancia, D(V)): necesita la MISMA fuente\n"
                "medida a dos o mas voltajes."
                % ("" if solo is None else " para el detector %s" % solo))
            return

        lins = [(p, registro.linealidad(e)) for p, e in conjuntos]
        ajustes = [(k, v, voltaje.ajustar(v["V"], v["D"], v["err"]))
                   for k, v in sorted(curvas.items())]
        hvs = sorted({self._det_hv(p)[1] for p, _ in lins})
        dets_a = sorted({self._det_hv(p)[0] for p, _ in lins}
                        | {k[0] for k, _, _ in ajustes})
        fuentes = sorted({f for _, l in lins for f in l["fuentes"]}
                         | {k[1] for k, _, _ in ajustes})

        def col_det(d):
            return self._COL_DET.get(d, "0.3")

        def col_hv(hv):
            return self._COL_HV[hvs.index(hv) % len(self._COL_HV)]

        def ls_hv(hv):
            return self._LINEA_HV[hvs.index(hv) % len(self._LINEA_HV)]

        def marca(f):
            return self._MARCA_FUENTE.get(f, "D")

        def relleno(f, c):
            # Cs-137 hueco y Co-60 lleno: se distinguen aunque el marcador
            # quede chico o tapado por la barra de error.
            return "white" if f == "Cs-137" else c

        # ---- la clave (lo que antes era la leyenda) ---------------------
        from matplotlib.lines import Line2D
        gris = "0.35"
        clave = []
        if solo is None:
            clave += [Line2D([], [], color=col_det(d), lw=2.5,
                             label="detector %s" % d) for d in dets_a]
            if lins:
                clave += [Line2D([], [], color=gris, ls=ls_hv(h),
                                 label="%s (test A)" % h) for h in hvs]
        else:
            if lins:
                clave += [Line2D([], [], color=col_hv(h), lw=2.5,
                                 label="%s (test A)" % h) for h in hvs]
                clave += [Line2D([], [], color=gris, ls="-",
                                 label="A: por el origen"),
                          Line2D([], [], color=gris, ls="--",
                                 label="A: ordenada libre")]
            if ajustes:
                clave += [Line2D([], [], color=self._COL_FUENTE.get(f, "0.3"),
                                 lw=2.5, label="%s (test B)" % f)
                          for f in fuentes
                          if any(k[1] == f for k, _, _ in ajustes)]
                clave += [Line2D([], [], color=gris, ls="-",
                                 label="B: potencia"),
                          Line2D([], [], color=gris, ls="--",
                                 label="B: recta")]
        clave += [Line2D([], [], color=gris, ls="none", marker=marca(f),
                         mfc=relleno(f, gris), ms=6, label=f)
                  for f in fuentes]

        # ---- las tablas de numeros --------------------------------------
        # Con un detector, una fila por punto o por serie. Comparando, eso son
        # 16 filas que no entran al lado del grafico: se pivotea a detector x
        # HV (test A) y detector x fuente (test B), 3-4 filas cada una.
        nota_b = None
        if solo is not None:
            filas_a = [[p, "%.4f" % l["D_origen"],
                        "%+.0f ± %.0f" % (l["ordenada"], l["err_ordenada"])]
                       for p, l in lins]
            cols_a = (["punto", "D origen", "ordenada"], [0.40, 0.26, 0.34])
            tit_a = "Test A  [mV·ns/keVee, mV·ns]"
            filas_b = [["%s %s" % (d, f), "%.2f ± %.2f" % (r["b"], r["err_b"]),
                        "si" if r["decide"] else "no"]
                       for (d, f), v, r in ajustes]
            cols_b = (["serie", "b (D ∝ V^b)", "decide"], [0.42, 0.36, 0.22])
        else:
            orden = {}
            for p, l in lins:
                det, hv = self._det_hv(p)
                orden.setdefault(det, {})[hv] = "%+.0f±%.0f" % (
                    l["ordenada"], l["err_ordenada"])
            filas_a = [[d] + [orden[d].get(h, "—") for h in hvs]
                       for d in sorted(orden)]
            cols_a = (["det"] + [h.replace(" V", "") + " V" for h in hvs],
                      [0.16] + [0.84 / max(len(hvs), 1)] * len(hvs))
            tit_a = "Test A: ordenada libre [mV·ns]"
            fb = [f for f in fuentes if any(k[1] == f for k, _, _ in ajustes)]
            tb = {}
            for (d, f), v, r in ajustes:
                tb.setdefault(d, {})[f] = "%.2f±%.2f%s" % (
                    r["b"], r["err_b"], "" if r["decide"] else "*")
            filas_b = [[d] + [tb[d].get(f, "—") for f in fb]
                       for d in sorted(tb)]
            cols_b = (["det"] + fb, [0.16] + [0.84 / max(len(fb), 1)] * len(fb))
            if any(not r["decide"] for _, _, r in ajustes):
                nota_b = "* el span de HV no alcanza para decidir"

        def dibujar(fig):
            ajustar(fig)   # antes de crear los ejes: ver `ajustar`
            paneles = int(bool(lins)) + int(bool(ajustes))
            ncol = 1 if solo is not None else 2
            alturas = [-(-len(clave) // ncol) + 1.5]
            if filas_a:
                alturas.append(len(filas_a) + 2.5)
            if filas_b:
                alturas.append(len(filas_b) + 2.5 + (1 if nota_b else 0))
            gs = fig.add_gridspec(len(alturas), paneles + 1,
                                  width_ratios=[1] * paneles + [0.62],
                                  height_ratios=alturas)
            k = 0
            if lins:
                ax = fig.add_subplot(gs[:, k]); k += 1
                Emax = max(float(np.max(l["E"])) for _, l in lins)
                x = np.linspace(0, 1.15 * Emax, 50)
                for p, l in lins:
                    det, hv = self._det_hv(p)
                    c = col_det(det) if solo is None else col_hv(hv)
                    ls = ls_hv(hv) if solo is None else "-"
                    ax.plot(x, l["D_origen"] * x, color=c, ls=ls, lw=1.3)
                    if solo is not None:
                        ax.plot(x, l["ordenada"] + l["D_libre"] * x,
                                color=c, ls="--", lw=1.0, alpha=.8)
                    for e_, c_, er_, f_ in zip(l["E"], l["canal"], l["err"],
                                              l["fuentes"]):
                        ax.errorbar(e_, c_, yerr=er_, fmt=marca(f_), ms=6,
                                    capsize=3, color=c, mfc=relleno(f_, c))
                ax.plot([0], [0], "+", color="k", ms=8)
                ax.set_xlim(0, 1.15 * Emax)
                ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
                ax.set_title("TEST A · el eje canal-energia")
                ax.set_xlabel("energia del borde [keVee]")
                ax.set_ylabel("canal del borde [mV·ns]")
                ax.grid(alpha=.3)
            if ajustes:
                ax = fig.add_subplot(gs[:, k]); k += 1
                for (d, f), v, r in ajustes:
                    c = col_det(d) if solo is None \
                        else self._COL_FUENTE.get(f, "0.3")
                    Vf = np.linspace(v["V"].min() * .99, v["V"].max() * 1.01,
                                     120)
                    ax.errorbar(v["V"], v["D"], yerr=v["err"], fmt=marca(f),
                                ms=6, capsize=3, color=c, mfc=relleno(f, c))
                    ax.plot(Vf, voltaje.evaluar(r, Vf), color=c, lw=1.3)
                    if solo is not None:
                        ax.plot(Vf, r["ordenada"] + r["pendiente"] * Vf,
                                color=c, ls="--", lw=1.0, alpha=.8)
                ax.set_title("TEST B · ganancia contra voltaje")
                ax.set_xlabel("alto voltaje [V]")
                ax.set_ylabel("D [mV·ns/keVee]")
                ax.grid(alpha=.3)

            # la columna de la derecha: clave arriba, numeros abajo
            fila = 0
            axk = fig.add_subplot(gs[fila, -1]); fila += 1
            axk.set_axis_off()
            axk._nombre_panel = "clave"
            axk.legend(handles=clave, loc="upper left", frameon=False,
                       fontsize=8, handlelength=2.4, borderaxespad=0,
                       ncol=ncol, columnspacing=1.0,
                       title=("clave" if solo is None
                              else "detector %s" % solo),
                       title_fontsize=9, alignment="left")
            if filas_a:
                axa = fig.add_subplot(gs[fila, -1]); fila += 1
                tabla(axa, filas_a, cols_a[0], cols_a[1], titulo=tit_a)
            if filas_b:
                axb = fig.add_subplot(gs[fila, -1]); fila += 1
                tabla(axb, filas_b, cols_b[0], cols_b[1],
                      titulo="Test B: exponente b (D ∝ V^b)", nota=nota_b)
        self.lz_registro.dibujar(dibujar)

    def _fig_ganancia(self, curvas, ref="A"):
        """
        El test B en escala log-log, y la ganancia de cada detector relativa
        a la del detector `ref`.

        Log-log y no semilog: la ley es de POTENCIA, D = D0 * (V/V0)^b, asi
        que log D = log D0 + b * log(V/V0) es una recta de pendiente b. Una
        exponencial (D ∝ e^(kV)) seria recta en semilog; en un span de 10% de
        HV las dos se parecen, pero la que se ajusta es la potencia.

        El cociente D_det / D_ref se toma solo a voltajes medidos en los dos
        detectores con la MISMA fuente (no hace falta extrapolar), y la curva
        es el cociente de las dos potencias ajustadas: (V/V0)^(b_det - b_ref)
        en log-log, o sea que si el cociente cambia con V es porque los b son
        distintos.
        """
        from matplotlib.lines import Line2D
        from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter

        ajustes = {k: (v, voltaje.ajustar(v["V"], v["D"], v["err"]))
                   for k, v in sorted(curvas.items())}
        dets = sorted({d for d, _ in ajustes})
        fuentes = sorted({f for _, f in ajustes})
        ls_fuente = {"Co-60": "-", "Cs-137": "--"}

        def col(d):
            return self._COL_DET.get(d, "0.3")

        def marca(f):
            return self._MARCA_FUENTE.get(f, "D")

        def relleno(f, c):
            return "white" if f == "Cs-137" else c

        # cocientes contra el detector de referencia, fuente por fuente
        cocientes = {}
        for (d, f), (v, r) in ajustes.items():
            if d == ref or (ref, f) not in ajustes:
                continue
            vr, rr = ajustes[(ref, f)]
            comunes = [x for x in v["V"] if x in set(vr["V"])]
            if not comunes:
                continue
            i = [list(v["V"]).index(x) for x in comunes]
            j = [list(vr["V"]).index(x) for x in comunes]
            q = v["D"][i] / vr["D"][j]
            eq = q * np.hypot(v["err"][i] / v["D"][i], vr["err"][j] / vr["D"][j])
            lo, hi = min(comunes), max(comunes)
            Vf = np.linspace(lo, hi, 60)
            cocientes[(d, f)] = dict(V=np.array(comunes), q=q, err=eq, Vf=Vf,
                                     qf=voltaje.evaluar(r, Vf)
                                     / voltaje.evaluar(rr, Vf),
                                     db=r["b"] - rr["b"],
                                     err_db=np.hypot(r["err_b"], rr["err_b"]))

        clave = [Line2D([], [], color=col(d), lw=2.5, label="detector %s" % d)
                 for d in dets]
        clave += [Line2D([], [], color="0.35", ls=ls_fuente.get(f, "-"),
                         marker=marca(f), mfc=relleno(f, "0.35"), ms=6,
                         label=f) for f in fuentes]

        dets_q = sorted({d for d, _ in cocientes})
        filas = []
        for d in dets_q:
            fila = [d]
            for f in fuentes:
                c = cocientes.get((d, f))
                fila.append("—" if c is None else "%.3f→%.3f" % (c["q"][0], c["q"][-1]))
            filas.append(fila)
        filas_db = []
        for d in dets_q:
            fila = [d]
            for f in fuentes:
                c = cocientes.get((d, f))
                fila.append("—" if c is None
                            else "%+.2f±%.2f" % (c["db"], c["err_db"]))
            filas_db.append(fila)
        anchos = [0.16] + [0.84 / max(len(fuentes), 1)] * len(fuentes)

        def ejes_hv(ax, Vs):
            ticks = sorted({float(x) for x in Vs})
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.set_xlabel("alto voltaje [V]  (escala log)")

        def dibujar(fig):
            ajustar(fig)   # antes de crear los ejes: ver `ajustar`
            alturas = [-(-len(clave) // 2) + 1.5]
            if filas:
                alturas += [len(filas) + 2.5, len(filas_db) + 2.5]
            gs = fig.add_gridspec(len(alturas), 3, width_ratios=[1, 1, 0.62],
                                  height_ratios=alturas)

            # ---- log-log: la potencia es una recta de pendiente b
            ax = fig.add_subplot(gs[:, 0])
            todos_V = []
            for (d, f), (v, r) in ajustes.items():
                c = col(d)
                Vf = np.linspace(v["V"].min() * .99, v["V"].max() * 1.01, 60)
                ax.errorbar(v["V"], v["D"], yerr=v["err"], fmt=marca(f), ms=6,
                            capsize=3, color=c, mfc=relleno(f, c))
                ax.plot(Vf, voltaje.evaluar(r, Vf), color=c, lw=1.2,
                        ls=ls_fuente.get(f, "-"))
                todos_V += list(v["V"])
            ax.set_xscale("log"); ax.set_yscale("log")
            ejes_hv(ax, todos_V)
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_formatter(ScalarFormatter())
            ax.tick_params(axis="y", which="minor", labelsize=7)
            ax.set_ylabel("D [mV·ns/keVee]  (escala log)")
            ax.set_title("TEST B en log-log: la pendiente es b")
            ax.grid(alpha=.3, which="both")

            # ---- cociente contra el detector de referencia
            ax = fig.add_subplot(gs[:, 1])
            if cocientes:
                ax.axhline(1.0, color=col(ref), lw=1.5)
                for (d, f), c in sorted(cocientes.items()):
                    cc = col(d)
                    ax.errorbar(c["V"], c["q"], yerr=c["err"], fmt=marca(f),
                                ms=6, capsize=3, color=cc, mfc=relleno(f, cc))
                    ax.plot(c["Vf"], c["qf"], color=cc, lw=1.2,
                            ls=ls_fuente.get(f, "-"))
                ax.set_xscale("log")
                ejes_hv(ax, [x for c in cocientes.values() for x in c["V"]])
            else:
                ax.text(0.5, 0.5, "sin voltajes en comun con el detector %s"
                        % ref, ha="center", va="center", color="0.45",
                        transform=ax.transAxes)
            ax.set_ylabel("D_det / D_%s  (misma fuente, mismo HV)" % ref)
            ax.set_title("ganancia relativa al detector %s" % ref)
            ax.grid(alpha=.3)

            # ---- clave y tablas
            axk = fig.add_subplot(gs[0, 2])
            axk.set_axis_off()
            axk._nombre_panel = "clave"
            axk.legend(handles=clave, loc="upper left", frameon=False,
                       fontsize=8, handlelength=2.6, borderaxespad=0, ncol=2,
                       columnspacing=1.0, title="clave (linea = potencia)",
                       title_fontsize=9, alignment="left")
            if filas:
                ax1 = fig.add_subplot(gs[1, 2])
                tabla(ax1, filas, ["det"] + fuentes, anchos,
                      titulo="D_det / D_%s: menor HV → mayor HV" % ref)
                ax2 = fig.add_subplot(gs[2, 2])
                tabla(ax2, filas_db, ["det"] + fuentes, anchos,
                      titulo="b_det − b_%s" % ref)
        self.lz_registro.dibujar(dibujar)


def _u(r):
    return "keVee" if r["rama"].startswith("Lvis") else "keV"


def _tiempo(s):
    if s is None:
        return ""
    s = int(round(s))
    return "%d s" % s if s < 60 else "%d min %02d s" % (s // 60, s % 60)


# --------------------------------------------------------------------------
# Que un error se VEA
# --------------------------------------------------------------------------
# Tkinter atrapa las excepciones que ocurren dentro de un callback (un boton, un
# `after`, el vaciado de la cola), las imprime en stderr y sigue andando. Si la
# GUI se lanzo SIN consola —doble clic, pythonw, un acceso directo— ese texto no
# va a ninguna parte: la ventana se queda a medio dibujar o el boton deja de
# responder, y desde afuera eso se ve igual que un cierre abrupto. O sea que la
# falla mas dificil de diagnosticar de este programa es una que el programa ya
# sabe, pero no cuenta.
#
# `report_callback_exception` es el gancho que Tk deja justo para esto. Aca se
# escribe la traza en `.cache/gui_error.log` y ademas se muestra, para que
# SIEMPRE haya algo que copiar y pegar.
RUTA_ERRORES = os.path.join(rutas.CACHE, "gui_error.log")
RUTA_FALLO = os.path.join(rutas.CACHE, "gui_fallo.log")
RUTA_FALLO_ANTERIOR = os.path.join(rutas.CACHE, "gui_fallo_anterior.log")

# Volcado del cierre abrupto ANTERIOR, si lo hubo. Lo rescata
# `activar_faulthandler` antes de truncar el archivo, y la ventana lo muestra
# al abrirse: un cierre sin mensaje deja de ser sin mensaje.
FALLO_ANTERIOR = None


def _rescatar_fallo():
    """
    Guarda el volcado de la sesion anterior antes de pisarlo, y lo devuelve.

    `faulthandler` necesita un archivo abierto en modo escritura, o sea que
    arrancar la interfaz BORRABA la unica evidencia que dejaba el cierre
    anterior. Si el usuario relanzaba el programa -que es lo primero que
    cualquiera hace- la pila se perdia y no quedaba nada que mirar.

    La cabecera son 2 lineas: mas que eso significa que hubo volcado.
    """
    try:
        if not os.path.exists(RUTA_FALLO):
            return None
        with open(RUTA_FALLO, encoding="utf-8", errors="replace") as f:
            texto = f.read()
        if len(texto.splitlines()) <= 2:
            return None                       # cierre limpio, nada que contar
        os.replace(RUTA_FALLO, RUTA_FALLO_ANTERIOR)
        return texto
    except OSError:
        return None


def activar_faulthandler():
    """
    Deja escrita la pila si el proceso MUERE sin excepcion.

    Hay fallas que no son excepciones de Python y que por lo tanto ningun
    try/except puede ver: una violacion de acceso en una extension en C, o —lo
    mas probable con estos datos— un EXCEPTION_IN_PAGE_ERROR de Windows cuando
    el disco USB se duerme o se desconecta mientras `drs_bin` esta leyendo el
    archivo por `np.memmap`. En los dos casos la ventana desaparece de golpe y
    no queda ni traza ni mensaje.

    `faulthandler` engancha esas senales a bajo nivel y vuelca la pila de TODOS
    los hilos —que es lo que importa, porque el trabajo corre en uno aparte— al
    archivo que se le pase. Es la unica forma de tener algo que leer despues de
    un cierre abrupto.
    """
    global FALLO_ANTERIOR
    try:
        import faulthandler
        FALLO_ANTERIOR = _rescatar_fallo()
        os.makedirs(rutas.CACHE, exist_ok=True)
        f = open(RUTA_FALLO, "w", encoding="utf-8")
        f.write("Si el programa se cierra solo, la pila queda abajo.\n"
                "Si el archivo queda solo con esta linea, el cierre no fue un "
                "fallo de este tipo.\n")
        f.flush()
        faulthandler.enable(file=f, all_threads=True)
        return f                      # se conserva para que no lo cierre el GC
    except Exception:                                       # noqa: BLE001
        return None


_FALLO = activar_faulthandler()


def registrar_error(exc, valor, traza, contexto="callback de Tk"):
    """Deja la traza en disco y la devuelve como texto, para mostrarla."""
    texto = "".join(traceback.format_exception(exc, valor, traza))
    cabecera = "%s  %s  (%s)" % ("=" * 70, time.strftime("%Y-%m-%d %H:%M:%S"),
                                 contexto)
    try:
        os.makedirs(rutas.CACHE, exist_ok=True)
        with open(RUTA_ERRORES, "a", encoding="utf-8") as f:
            f.write("\n%s\n%s\n" % (cabecera, texto))
    except OSError:
        pass                          # no poder loguear no puede ser el error
    return texto


def _reportar_excepcion(self, exc, valor, traza):
    texto = registrar_error(exc, valor, traza)
    try:
        self._escribir("\n" + "=" * 70 + "\nERROR EN LA INTERFAZ\n"
                       + "=" * 70 + "\n" + texto)
        self._fin_corrida("ERROR: %s" % str(valor)[:60])
    except Exception:                                       # noqa: BLE001
        pass                          # la ventana puede estar a medio construir
    messagebox.showerror(
        "Error en la interfaz",
        "%s\n\nLa traza completa quedo en:\n%s"
        % (str(valor)[:300], RUTA_ERRORES))


App.report_callback_exception = _reportar_excepcion


def main():
    if _FALLO is None:
        activar_faulthandler()
    try:
        app = App()
    except Exception as e:                                  # noqa: BLE001
        # Un fallo ANTES de que exista la ventana no tiene donde mostrarse: no
        # hay pestanas ni messagebox todavia. Queda en el log igual.
        registrar_error(type(e), e, e.__traceback__, "arranque de la ventana")
        print("No se pudo abrir la interfaz. La traza quedo en %s"
              % RUTA_ERRORES)
        raise
    app.mainloop()


if __name__ == "__main__":
    main()
