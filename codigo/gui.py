"""
gui.py - Visor por etapas del analisis de calibracion.

Una pestana por etapa. Cada una tiene:
  - CONTROLES arriba: los parametros que de verdad mueven el resultado.
  - GRAFICO al centro (con la barra de zoom/pan de matplotlib).
  - EXPLICACION a la derecha: que es lo que estas viendo y que mirar.

La idea no es mostrar PNGs ya hechos sino poder RECALCULAR: mover sigma, la
ventana, la semilla del SA, y ver que se rompe. Ahi es donde se ven las mejoras.

Los calculos pesados (leer los XML, NCC, SA, forward folding) corren en un hilo
aparte para que la ventana no se congele; el estado se ve abajo a la izquierda.

La primera pestana ("Archivos") es la que elige QUE se analiza: el .root de la
simulacion y los dos .xml (muestra y fondo). Lo que se elija ahi vale para todas
las demas pestanas y se recuerda entre sesiones en .cache/archivos.json.

Uso:  python gui.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/gui.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.abspath(__file__))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas


import json
import os
import queue
import sys
import threading
import traceback
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from datos_io import spectrum
from etapas.physics import compton_edge_energy, dsigma_dE, Fuente, FUENTES
from etapas.model import forward_model, punto_mitad_altura
from etapas.fit_single import ajustar_borde_unico, _espectro_sintetico
from etapas.ncc import (buscar_intervalo, buscar_intervalo_simulado,
                 plantilla_referencia, plantilla_simulada)
from etapas.clasicos import comparar as comparar_clasicos
from etapas.fit_double import (ajustar_dos_bordes, componentes, _amplitudes_optimas,
                        delta_E, modelo_dos_bordes)
import calibracion as cal

def activar_dpi_awareness():
    """
    Le dice a Windows que esta app maneja su propio DPI.

    Sin esto, en pantallas con escalado (125%, 150%, 200%) Windows renderiza la
    ventana a 96 DPI y despues la ESTIRA como bitmap: por eso se ve borrosa, y
    ninguna configuracion de matplotlib lo arregla. Hay que llamarlo ANTES de
    crear la ventana Tk.
    """
    if sys.platform != "win32":
        return
    import ctypes
    try:                                    # Windows 10 1703+: per-monitor v2
        ctypes.windll.shcore.SetProcessDpiAwareness(2)
    except Exception:
        try:
            ctypes.windll.user32.SetProcessDPIAware()   # fallback, Vista+
        except Exception:
            pass


ARCHIVO_MUESTRA = rutas.MUESTRA
ARCHIVO_FONDO = rutas.FONDO
ARCHIVO_SIM = rutas.SIM
ARBOL_SIM = "ComptonSpectrum"
FUENTE_DEF = "Co-60"


def formato_tiempo(segundos):
    """'45 s' / '2 min 10 s' / '1 h 05 min'. Para el ETA de la barra."""
    seg = max(int(round(segundos)), 0)
    if seg < 60:
        return "%d s" % seg
    if seg < 3600:
        return "%d min %02d s" % (seg // 60, seg % 60)
    return "%d h %02d min" % (seg // 3600, (seg % 3600) // 60)


class SinSimulacion(RuntimeError):
    """
    La pestana necesita la simulacion Geant4 y no hay ninguna cargada.

    Es una condicion ESPERADA desde que el .root es opcional (se puede calibrar
    solo con Klein-Nishina), asi que se muestra como aviso y no como error con
    traceback: ver `Pestana.recibir`.
    """

RUTA_PREFS = os.path.join(rutas.CACHE, "archivos.json")


def cargar_preferencias():
    """
    Ultima seleccion de archivos, si sigue existiendo.

    Se guarda aparte de los caches de datos porque es preferencia de la GUI, no
    resultado: si el archivo se movio o se borro, se vuelve callado al default
    en vez de reventar al arrancar.
    """
    elegidos = dict(muestra=ARCHIVO_MUESTRA, fondo=ARCHIVO_FONDO,
                    sim=ARCHIVO_SIM, arbol=ARBOL_SIM, fuente=FUENTE_DEF,
                    lineas=None)
    try:
        with open(RUTA_PREFS, encoding="utf-8") as f:
            guardado = json.load(f)
    except Exception:
        return elegidos
    for clave in ("muestra", "fondo"):
        ruta = rutas.absoluta(guardado.get(clave))
        if ruta and os.path.exists(ruta):
            elegidos[clave] = ruta
    # la simulacion es OPCIONAL: "" es una eleccion valida (ruta analitica)
    if "sim" in guardado:
        sim = rutas.absoluta(guardado["sim"] or "")
        elegidos["sim"] = sim if (not sim or os.path.exists(sim)) else ARCHIVO_SIM
    if guardado.get("arbol"):
        elegidos["arbol"] = guardado["arbol"]
    if guardado.get("fuente"):
        elegidos["fuente"] = guardado["fuente"]
        elegidos["lineas"] = guardado.get("lineas")
    return elegidos


def guardar_preferencias(datos):
    try:
        os.makedirs(os.path.dirname(RUTA_PREFS), exist_ok=True)
        with open(RUTA_PREFS, "w", encoding="utf-8") as f:
            json.dump({"muestra": rutas.relativa(datos.archivo_muestra),
                       "fondo": rutas.relativa(datos.archivo_fondo),
                       "sim": rutas.relativa(datos.archivo_sim),
                       "arbol": datos.arbol_sim,
                       "fuente": datos.fuente.nombre,
                       "lineas": list(datos.fuente.lineas)}, f, indent=1)
    except Exception:
        pass                                # no poder guardar no es un error


def info_archivo(ruta):
    """(existe, texto corto) para la etiqueta de estado de cada archivo."""
    import time
    if not ruta:
        return False, "sin elegir"
    if not os.path.exists(ruta):
        return False, "NO EXISTE"
    st = os.stat(ruta)
    fecha = time.strftime("%d-%b %H:%M", time.localtime(st.st_mtime))
    txt = f"{st.st_size/1e6:,.0f} MB · {fecha}"
    if ruta.lower().endswith(".root"):
        return True, txt
    # El FORMATO se pregunta por los primeros bytes y no por la extensión: los
    # binarios del DRS4 de la sesión del 1-sep-2026 no tienen extensión en el
    # nombre (se llaman `Co60`, `Cs137`), así que gatear esto por ".xml" los
    # dejaba afuera sin decir nada — se veía como un archivo del que no se sabe
    # ni el canal.
    fmt = spectrum.formato(ruta)
    txt += f"  · {fmt}"
    # el canal sale de mirar el primer evento (o la cabecera, en el binario): es
    # la falla más fácil de cometer (las corridas de junio están en CHN2 y las
    # de agosto en CHN1) y la más difícil de diagnosticar después.
    try:
        canal, tag = spectrum.detectar_canal(ruta)
        txt += f"  · {canal}/{tag or 'sin scaler'}"
    except Exception:
        txt += ("  · ⚠ sin <CHNn> legible" if fmt == "xml"
                else "  · ⚠ cabecera binaria ilegible")
    txt += ("  · con caché" if os.path.exists(ruta + ".areas_gate.npz")
            else "  · SIN caché (parseo lento la 1ª vez)")
    return True, txt


def inspeccionar_root(ruta):
    """
    Lista los TTree de un .root con sus ramas, SIN leer los datos.

    uproot abre solo las cabeceras, asi que esto cuesta milisegundos aunque el
    archivo pese 31 MB. Sirve para dos cosas: elegir el arbol (el default
    'ComptonSpectrum' es el nombre que usa esta simulacion, no una ley) y avisar
    si faltan las ramas que el analisis necesita.
    """
    import uproot
    arboles = {}
    with uproot.open(ruta) as f:
        for clave in f.keys(cycle=False):
            try:
                obj = f[clave]
            except Exception:
                continue
            if hasattr(obj, "num_entries") and hasattr(obj, "keys"):
                arboles[clave] = (list(obj.keys()), int(obj.num_entries))
    return arboles


# ======================================================================
#  cache compartido de datos (para no releer los XML en cada pestana)
# ======================================================================
class Datos:
    def __init__(self, muestra=ARCHIVO_MUESTRA, fondo=ARCHIVO_FONDO,
                 sim=ARCHIVO_SIM, arbol=ARBOL_SIM, fuente=None):
        self._cache = {}
        self._sim = {}
        self._borde = {}
        self.archivo_muestra = muestra
        self.archivo_fondo = fondo
        self.archivo_sim = sim
        self.arbol_sim = arbol
        self.fuente = fuente or Fuente(FUENTE_DEF)
        self.forzar_sim = False      # usar el .root aunque sea de otra fuente
        # Gancho de progreso: `Pestana.lanzar` deja aca su propio `progreso`,
        # asi que el parseo de los XML (el unico paso que tarda ~40 s sin decir
        # nada) alimenta la barra sin que cada pestana tenga que enterarse.
        self.avisar = None
        self._lineas_sim = {}
        self._comp = {}

    # ------------------------------------------------------------------
    def hay_simulacion(self):
        """La simulacion es opcional: sin ella queda la ruta Klein-Nishina."""
        return bool(self.archivo_sim) and os.path.exists(self.archivo_sim)

    def lineas_sim(self):
        """
        Gammas primarios del .root (rama E0_keV), o None si no se pueden leer.

        Cacheado por (archivo, arbol): leer la rama de 850 k entradas cuesta
        ~0.2 s, pero se pregunta en toda pestana que toque la simulacion.
        """
        if not self.hay_simulacion():
            return None
        clave = (self.archivo_sim, self.arbol_sim)
        if clave not in self._lineas_sim:
            from ajuste.sim_template import lineas_primarias
            try:
                self._lineas_sim[clave] = lineas_primarias(self.archivo_sim,
                                                           arbol=self.arbol_sim)
            except Exception:
                self._lineas_sim[clave] = None
        return self._lineas_sim[clave]

    def compat_simulacion(self):
        """
        (compatible, mensaje): si la plantilla Geant4 es de ESTA fuente.

        Es el chequeo que faltaba al hacer la fuente configurable. Elegir
        Cs-137 y dejar puesto el .root de Co-60 no rompe nada visible: el
        forward folding devuelve un D perfectamente formateado y equivocado,
        porque la plantilla tiene sus bordes en 963 y 1118 keV mientras los
        datos tienen uno solo en 477. Cuesta una rama leer E0_keV y evita
        exactamente ese modo de falla.
        """
        from ajuste.sim_template import comparar_con_fuente
        if not self.hay_simulacion():
            return False, "no hay simulacion cargada."
        return comparar_con_fuente(self.lineas_sim(), self.fuente.lineas)

    def exigir_sim_compatible(self):
        """Corta con `SinSimulacion` si la plantilla es de otra fuente."""
        ok, msg = self.compat_simulacion()
        if ok or self.forzar_sim:
            return
        raise SinSimulacion(
            "La simulacion no corresponde a la fuente elegida.\n\n" + msg
            + "\n\nOpciones: cargar el .root de %s, quitar la simulacion (la "
            "ruta analitica calibra igual: pestanas 4 y 6), o tildar «usar la "
            "simulacion igual» en la pestana «Archivos»." % self.fuente.nombre)

    def componentes(self, rama="Lvis_keVee", bins=175):
        """
        (centros, h_1scattering, h_multiple, w_simulacion) de `modelo_ext`.

        Vive aca y no en cada pestana porque el barrido de sistematicos de la
        pestana de Calibracion pide el mismo histograma media docena de veces,
        y releer el .root cada vez es lo unico lento de ese barrido.
        """
        self.exigir_sim_compatible()
        clave = (rama, bins, self.techo_sim())
        if clave not in self._comp:
            from ajuste import modelo_ext as me
            self._comp[clave] = me.cargar_componentes(
                ruta=self.archivo_sim, arbol=self.arbol_sim, rama=rama,
                bins=bins, rango=(0.0, self.techo_sim()))
        return self._comp[clave]

    def fijar_fuente(self, fuente):
        """
        Cambia la fuente (Co-60, Cs-137, ...). Devuelve True si cambio.

        No toca el cache del espectro —la fuente no interviene en como se lee el
        XML— pero SI el de la simulacion, porque el techo del histograma se saca
        de la linea mas energetica.
        """
        viejo = (self.fuente.nombre, self.fuente.lineas)
        if (fuente.nombre, fuente.lineas) == viejo:
            return False
        self.fuente = fuente
        self._sim.clear()
        self._borde.clear()
        self._comp.clear()
        return True

    def borde_estimado(self, rango_max=2500.0, num_canales=500, umbral=0.05,
                       metodo="gate"):
        """
        B del borde de la PRIMERA linea, estimado por NCC con Klein-Nishina.

        Para que: las ventanas de ajuste por defecto (300-1900 mV*ns) valen para
        el SiPM con Co-60 y para nada mas. Con otro detector o otra fuente hay
        que ponerlas en escala, y la unica forma de hacerlo sin tipear numeros a
        mano es preguntarle al propio espectro donde esta el borde. Es la misma
        idea de la Seccion II-B del paper, usada aca como regla de encuadre.

        Se cachea porque lo piden varias pestanas y cuesta ~1 s.
        """
        clave = (rango_max, num_canales, umbral, metodo,
                 self.fuente.nombre, self.fuente.lineas)
        if clave not in self._borde:
            from etapas.ncc import buscar_intervalo
            r = self.espectro(rango_max, num_canales, umbral, metodo)
            can, neto = r["canales"], r["neto"]
            suave = spectrum.suavizar(neto, ventana=21, orden=3)
            tope = float(can[-1])
            self._borde[clave] = buscar_intervalo(
                can, suave, self.fuente.lineas[0], region=(0.08 * tope, tope),
                paso=1, n_puntos="auto", cuentas_crudas=neto)
        return self._borde[clave]

    def ventana_ajuste(self, manual=None, lo_rel=0.25, hi_rel=1.55):
        """
        Ventana de ajuste [mV*ns]. Si `manual` no es None, esa; si no, del NCC.

        Los factores estan calibrados contra la Etapa 9: con B = 1229 mV*ns
        (Co-60, area con compuerta) dan (307, 1905), o sea los 300-1900 que se
        venian tipeando a mano. La diferencia es que ahora se adaptan solos al
        detector y a la fuente.
        """
        if manual is not None:
            return (float(manual[0]), float(manual[1]))
        B = float(self.borde_estimado()["edge_ch"])
        tope = float(self.espectro()["canales"][-1])
        return (lo_rel * B, min(hi_rel * B, tope))

    def techo_sim(self):
        """
        Tope del histograma de la simulacion [keV].

        1400 es el valor con el que se hicieron las Etapas 7-9 (Co-60 llega a
        1332.5), y se conserva para no mover el binning de esos resultados. Solo
        sube si la fuente tiene lineas por encima.
        """
        return max(1400.0, 1.2 * max(self.fuente.lineas))

    # ------------------------------------------------------------------
    def fijar_archivos(self, muestra, fondo, sim, arbol=ARBOL_SIM):
        """
        Cambia los archivos a analizar y tira los caches EN MEMORIA.

        Los caches de areas de spectrum.py viven en disco al lado de cada XML y
        llevan la huella (tamano + mtime + parametros) del archivo, asi que esos
        no hay que tocarlos: cada XML tiene el suyo y se invalidan solos. Los de
        aca si, porque estan indexados por (rango, bins, umbral, metodo) y no
        saben de que archivo salieron.

        Devuelve True si algo cambio de verdad.
        """
        nuevo = (muestra, fondo, sim, arbol)
        viejo = (self.archivo_muestra, self.archivo_fondo, self.archivo_sim,
                 self.arbol_sim)
        self.archivo_muestra, self.archivo_fondo = muestra, fondo
        self.archivo_sim, self.arbol_sim = sim, arbol
        if nuevo == viejo:
            return False
        self._cache.clear()
        self._sim.clear()
        self._borde.clear()
        self._comp.clear()
        self._lineas_sim.clear()
        return True

    def huella(self):
        """
        Identidad corta de los tres archivos (nombre + tamano + fecha).

        Entra en las claves de los caches EN DISCO (por ejemplo las cadenas de
        MCMC de la pestana 8): sin esto, cambiar de datos y reusar un .npz viejo
        daria un posterior que no corresponde a lo que se esta mirando. Incluye
        la FUENTE por la misma razon: el mismo espectro ajustado contra Cs-137 o
        contra Co-60 da D distinto.
        """
        import hashlib
        h = hashlib.md5()
        h.update(repr((self.fuente.nombre, self.fuente.lineas)).encode("utf-8"))
        for ruta in (self.archivo_muestra, self.archivo_fondo, self.archivo_sim):
            try:
                st = os.stat(ruta)
                marca = f"{os.path.basename(ruta)}|{st.st_size}|{int(st.st_mtime)}"
            except OSError:
                marca = f"{ruta}|?"
            h.update(marca.encode("utf-8"))
        return h.hexdigest()[:8]

    def espectro(self, rango_max=2500.0, num_canales=500, umbral=0.05,
                 metodo="gate"):
        """
        Espectro neto. `metodo="gate"` es el area buena (linea base pre-onset +
        compuerta, CLAUDE.md 7.33) y normaliza por tiempo VIVO (7.34);
        `metodo="legacy"` reproduce el analisis viejo, con el area rota, para
        poder ver el antes/despues. Ojo que el eje cambia de escala: con el area
        buena todo pasa entre 0 y ~1900 mV*ns, con la vieja llegaba a 8000.
        """
        clave = (rango_max, num_canales, umbral, metodo)
        if clave not in self._cache:
            avisar, cb = self.avisar, None
            if avisar is not None:
                def cb(frac, texto):
                    avisar("leyendo los datos — " + texto, frac)
            # `archivo_fondo` vacío es una elección válida: las corridas del
            # 1-sep-2026 no tienen fondo propio. `espectro_neto` lo entiende y
            # devuelve neto = muestra con k_f = 0.
            self._cache[clave] = spectrum.espectro_neto(
                self.archivo_muestra, self.archivo_fondo or None,
                umbral_ruido=umbral,
                num_canales=num_canales, rango_max=rango_max,
                recortar_negativos=False, metodo=metodo, callback=cb)
        return self._cache[clave]

    def simulacion(self, rama="Lvis_keVee", filtro_nombre="todo", bins=700):
        """Histograma de la simulacion; `SinSimulacion` si no hay .root."""
        """
        Histograma de la simulacion. `bins` entra en la clave del cache porque
        el SA y el MCMC piden 175 (mismo resultado, 6x mas rapido: ver
        sim_template.cargar_simulacion) y las otras pestanas siguen con 700.
        """
        if not self.hay_simulacion():
            raise SinSimulacion(
                "No hay simulación Geant4 cargada.\n\n"
                "Elegí un .root en la pestaña «Archivos», o quedate con la ruta "
                "analítica (Klein-Nishina): las pestañas 4 · NCC y "
                "6 · Ajuste analítico calibran sin simulación.")
        self.exigir_sim_compatible()
        techo = self.techo_sim()
        clave = (rama, filtro_nombre, bins, techo)
        if clave not in self._sim:
            from ajuste.sim_template import cargar_simulacion
            f = {"todo": None,
                 "1 scattering": lambda d: d["NCompton"] == 1,
                 "multiple": lambda d: d["NCompton"] > 1}[filtro_nombre]
            self._sim[clave] = cargar_simulacion(self.archivo_sim, rama=rama,
                                                 filtro=f, bins=bins,
                                                 arbol=self.arbol_sim,
                                                 rango=(0.0, techo))
        return self._sim[clave]


# ======================================================================
#  pestana base
# ======================================================================
class Pestana(ttk.Frame):
    titulo = "pestana"
    explicacion = ""
    n_ejes = 1
    auto = True             # calcular al abrirla por primera vez

    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.calculada = False
        self._corriendo = False
        self._rango = (0.0, 1.0)     # sub-rango de la barra para el paso actual

        # --- controles ---
        self.ctrl = ttk.LabelFrame(self, text="Controles")
        self.ctrl.pack(side="top", fill="x", padx=6, pady=(6, 3))
        self.construir_controles(self.ctrl)
        ttk.Button(self.ctrl, text="Recalcular",
                   command=self.lanzar).pack(side="right", padx=6, pady=4)

        cuerpo = ttk.Frame(self)
        cuerpo.pack(side="top", fill="both", expand=True, padx=6, pady=3)

        # --- explicacion (derecha) ---
        panel = ttk.LabelFrame(cuerpo, text="Que estas viendo")
        panel.pack(side="right", fill="y", padx=(6, 0))
        self.txt = tk.Text(panel, width=46, wrap="word", relief="flat",
                           font=("Segoe UI", 9), background="#f7f7f7",
                           padx=4, pady=4)
        self.txt.pack(side="top", fill="both", expand=True, padx=6, pady=6)
        self.escribir(self.explicacion)

        # --- grafico ---
        marco = ttk.Frame(cuerpo)
        marco.pack(side="left", fill="both", expand=True)
        self.fig = Figure(figsize=(8.5, 5.6), dpi=app.dpi_fig)
        self.canvas = FigureCanvasTkAgg(self.fig, master=marco)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        NavigationToolbar2Tk(self.canvas, marco).update()

    # ---------- utilidades ----------
    def fuente(self):
        """La fuente elegida en la pestana Archivos (physics.Fuente)."""
        return self.app.datos.fuente

    def exigir_simulacion(self):
        """Corta la tarea con un aviso limpio si no hay .root cargado."""
        if not self.app.datos.hay_simulacion():
            raise SinSimulacion(
                "«%s» necesita la simulación Geant4 y no hay ninguna cargada.\n\n"
                "Elegí un .root en la pestaña «Archivos», o usá la ruta "
                "analítica (Klein-Nishina): las pestañas 4 · NCC y "
                "6 · Ajuste analítico no la necesitan." % self.titulo)

    def al_cambiar_fuente(self, fuente):
        """Gancho: se llama en el hilo principal cuando cambia la fuente."""
        pass

    def titulo_actual(self):
        """Texto de la solapa; puede depender de la fuente (ver TabAnalitico)."""
        return self.titulo

    def escribir(self, texto):
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        self.txt.insert("1.0", texto.strip())
        self.txt.configure(state="disabled")

    def construir_controles(self, frame):
        pass

    def trabajo(self):
        """Se ejecuta EN OTRO HILO. No tocar widgets aca. Devuelve un resultado."""
        return None

    def progreso(self, texto, frac=None):
        """
        Manda un aviso a la barra de estado DESDE el hilo de calculo.

        Tkinter no es thread-safe, asi que no se puede tocar el widget desde el
        hilo: se manda por la misma cola que el resultado y el hilo principal lo
        desvia en `App.revisar_cola`. Hace falta porque el MCMC tarda minutos y
        una barra de estado congelada parece un cuelgue.

        frac : avance de 0 a 1, o None si este paso no sabe cuanto falta (ahi la
               barra queda en modo indeterminado, que igual muestra que el
               programa esta vivo). Se re-escala al sub-rango declarado con
               `fase`, para que un paso que reporta su propio 0-100% (el parseo
               del XML) no haga saltar la barra del resto de la tarea.
        """
        if frac is not None:
            a, b = self._rango
            frac = a + max(0.0, min(float(frac), 1.0)) * (b - a)
        self.app.cola.put((self, ("progreso", texto, frac)))

    def fase(self, a, b):
        """Declara que lo que viene ocupa el tramo [a, b] de la barra."""
        self._rango = (float(a), float(b))

    def pintar(self, res):
        """Se ejecuta en el hilo principal con lo que devolvio trabajo()."""
        pass

    def lanzar(self):
        if self._corriendo:
            return
        self._corriendo = True
        self.calculada = True
        self._rango = (0.0, 1.0)
        self.app.datos.avisar = self.progreso
        self.app.estado.set(f"calculando: {self.titulo} ...")
        self.app.progreso_iniciar(self)

        def worker():
            try:
                res = ("ok", self.trabajo())
            except SinSimulacion as e:
                res = ("aviso", str(e))
            except Exception:
                res = ("error", traceback.format_exc())
            self.app.cola.put((self, res))

        threading.Thread(target=worker, daemon=True).start()

    def recibir(self, res):
        self._corriendo = False
        self.app.progreso_terminar(self)
        estado, carga = res
        if estado == "aviso":
            self._ultimo = None
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, carga.split("\n")[0], ha="center", va="center",
                    color="0.35", fontsize=11, wrap=True)
            ax.set_axis_off()
            self.canvas.draw()
            self.escribir(carga + "\n\n" + self.explicacion.strip())
            self.app.estado.set("falta la simulación")
            return
        if estado == "error":
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Error al calcular\n(ver panel derecho)",
                    ha="center", va="center", color="crimson")
            ax.set_axis_off()
            self.canvas.draw()
            self.escribir("ERROR:\n\n" + carga)
            self.app.estado.set("error")
            return
        self._ultimo = carga            # guardado para repintar sin recalcular
        self.fig.clear()
        self.pintar(carga)
        self.fig.tight_layout()
        self.canvas.draw()
        self.app.estado.set("listo")

    def repintar(self):
        """Redibuja con el ultimo resultado (para cambios de nitidez/tamano)."""
        carga = getattr(self, "_ultimo", None)
        if carga is None or self._corriendo:
            return
        self.fig.clear()
        try:
            self.pintar(carga)
            self.fig.tight_layout()
        except Exception:
            pass
        self.canvas.draw_idle()


# ======================================================================
#  Archivos - que datos entran al analisis
# ======================================================================
class TabArchivos(Pestana):
    titulo = "Archivos"
    explicacion = """
DE DONDE SALEN LOS DATOS.

Acá se elige QUÉ se analiza, y de acá sale la configuración de todo el resto:

• Muestra — la corrida con la fuente radiactiva.
• Fondo   — la corrida sin fuente, que se resta. Se puede dejar VACÍO: las
  corridas del 1-sep-2026 se tomaron sin fondo propio, y prestarles el de otra
  sesión sería peor que no restar nada (el fondo depende del punto de
  operación: a 1700 V el mismo cósmico deja el triple de área). Sin fondo el
  "neto" es la muestra cruda, y la etiqueta lo dice.
• Simulación (.root) — el Geant4 que sirve de plantilla (Etapas 7-9). OPCIONAL,
  y tiene que ser de la MISMA fuente: se verifica leyendo sus gammas primarios.
• Fuente — las líneas gamma de la muestra. Es lo que fija los bordes Compton,
  el ancla E_ref del mapa por el origen, qué modelo se ajusta (dos bordes
  solapados con dos líneas, uno solo con una) y cómo se cierran las ventanas
  del NCC. Hay presets (Cs-137, Co-60, Na-22, …) y "manual" para tipear una o
  dos energías.

Los dos formatos del DRSOSC valen y NO hay que decir cuál es: el XML de junio y
del 26-ago, y el BINARIO de "Save as binary" del 1-sep, que además viene sin
extensión en el nombre (`Co60`, `Cs137`, `B1700V`). Se decide por los primeros
4 bytes del archivo, y la etiqueta de abajo dice qué se detectó. Por eso el
diálogo abre mostrando TODOS los archivos: con el filtro "*.xml" las carpetas
de la sesión nueva se ven vacías.

La selección se guarda en .cache/archivos.json y se recuerda la próxima vez.

QUE MIRAR ANTES DE SEGUIR
• Que las dos corridas sean del MISMO detector y la MISMA ganancia. Si no, la
  resta de fondo no significa nada: el mismo mV·ns sería otra energía.
• El tiempo VIVO de cada una (§7.34). Si el scaler de la placa dice mucho más
  que los eventos guardados, la corrida perdió triggers y normalizar por tiempo
  de pared sobre-resta el fondo.
• El percentil 99.9 del área: RANGO_MAX tiene que dejar aire ARRIBA del borde
  (§7.6). Si el nuevo detector tiene otra ganancia, el 2500 mV·ns de default ya
  no sirve y hay que ajustarlo en cada pestaña.
• Que el .root tenga las ramas que el análisis usa: Lvis_keVee, Edep_keV y
  NCompton. Si falta NCompton no se puede separar scattering simple de múltiple,
  que es lo que mide `w` en la Etapa 9.

OJO
Cambiar los archivos o la fuente invalida lo ya calculado: las demás pestañas se
marcan como pendientes y hay que apretar "Recalcular" en cada una. Los números
de referencia que aparecen en las explicaciones (D = 1.26 mV·ns/keVee, bordes en
1215 y 1410) son del Co-60 con ESTE detector; con otros archivos o con otra
fuente son sólo un ejemplo. El número que vale para lo que tenés cargado lo
calcula la pestaña ★ Calibración.
"""

    # "Todos" va PRIMERO a propósito: los binarios del DRS4 de la sesión del
    # 1-sep-2026 no tienen extensión (`Co60`, `Cs137`, `B1700V`), así que con el
    # filtro "*.xml" por defecto la carpeta se ve VACÍA y parece que no hay
    # datos. El formato se detecta después por los primeros bytes.
    TIPOS_XML = [("Datos del DRS4 (XML o binario)", "*.*"),
                 ("XML del DRSOSC", "*.xml")]
    TIPOS_ROOT = [("ROOT de Geant4", "*.root"), ("Todos", "*.*")]
    RAMAS_ESPERADAS = ("Lvis_keVee", "Edep_keV", "NCompton")

    def construir_controles(self, frame):
        d = self.app.datos
        self.var = {"muestra": tk.StringVar(value=d.archivo_muestra),
                    "fondo": tk.StringVar(value=d.archivo_fondo),
                    "sim": tk.StringVar(value=d.archivo_sim)}
        self.arbol = tk.StringVar(value=d.arbol_sim)
        self.rango = tk.DoubleVar(value=2500.0)
        self.parsear = tk.BooleanVar(value=False)
        self.nombre_fuente = tk.StringVar(value=d.fuente.nombre)
        self.forzar = tk.BooleanVar(value=d.forzar_sim)
        self.eg1 = tk.DoubleVar(value=d.fuente.lineas[0])
        self.eg2 = tk.DoubleVar(value=d.fuente.lineas[1] if d.fuente.doble else 0.0)
        self.etiq = {}

        # grid adentro de un sub-frame: el LabelFrame de la clase base usa pack
        # (le mete el boton "Recalcular" a la derecha) y Tk no deja mezclar los
        # dos administradores en el mismo contenedor.
        caja = ttk.Frame(frame)
        caja.pack(side="left", fill="x", expand=True, padx=4, pady=2)
        filas = [("muestra", "Muestra (con fuente)", self.TIPOS_XML),
                 ("fondo", "Fondo (sin fuente)", self.TIPOS_XML),
                 ("sim", "Simulación Geant4 .root (opcional)", self.TIPOS_ROOT)]
        for i, (clave, etiqueta, tipos) in enumerate(filas):
            ttk.Label(caja, text=etiqueta).grid(row=i, column=0, sticky="w", pady=1)
            ttk.Entry(caja, textvariable=self.var[clave], width=54).grid(
                row=i, column=1, sticky="we", padx=4)
            ttk.Button(caja, text="Examinar…",
                       command=lambda c=clave, t=tipos: self.elegir(c, t)).grid(
                row=i, column=2, padx=2)
            lb = ttk.Label(caja, text="", foreground="#666")
            lb.grid(row=i, column=3, sticky="w", padx=(6, 0))
            self.etiq[clave] = lb
        caja.columnconfigure(1, weight=1)

        pie = ttk.Frame(caja)
        pie.grid(row=len(filas), column=0, columnspan=4, sticky="w", pady=(5, 0))
        ttk.Label(pie, text="TTree").pack(side="left", padx=(0, 3))
        self.cb_arbol = ttk.Combobox(pie, textvariable=self.arbol, width=16,
                                     values=[d.arbol_sim])
        self.cb_arbol.pack(side="left")
        ttk.Label(pie, text="RANGO_MAX (vista previa)").pack(side="left", padx=(12, 3))
        ttk.Entry(pie, textvariable=self.rango, width=7).pack(side="left")
        ttk.Checkbutton(pie, text="parsear los XML aunque no haya caché",
                        variable=self.parsear).pack(side="left", padx=12)

        # --- LA FUENTE ---
        # De acá sale cuántas líneas gamma tiene lo medido, y con eso qué modelo
        # aplica: dos bordes solapados (Ec. 11-13) con Co-60, un solo borde con
        # Cs-137. Nadie más en el código tiene las energías cableadas: todos le
        # preguntan a physics.Fuente.
        pie2 = ttk.Frame(caja)
        pie2.grid(row=len(filas) + 1, column=0, columnspan=4, sticky="w",
                  pady=(6, 0))
        ttk.Label(pie2, text="Fuente").pack(side="left", padx=(0, 3))
        cb = ttk.Combobox(pie2, textvariable=self.nombre_fuente, width=9,
                          state="readonly", values=list(FUENTES) + ["manual"])
        cb.pack(side="left")
        cb.bind("<<ComboboxSelected>>", self.al_elegir_preset)
        for etiq, var in (("E_γ1 [keV]", self.eg1), ("E_γ2", self.eg2)):
            ttk.Label(pie2, text=etiq).pack(side="left", padx=(8, 2))
            e = ttk.Entry(pie2, textvariable=var, width=8)
            e.pack(side="left")
            e.bind("<FocusOut>", lambda _e: self.refrescar_etiquetas())
        ttk.Label(pie2, text="(E_γ2 = 0 → una sola línea)",
                  foreground="#666").pack(side="left", padx=6)
        self.lbl_fuente = ttk.Label(pie2, text="", foreground="#444")
        self.lbl_fuente.pack(side="left", padx=12)

        botones = ttk.Frame(caja)
        botones.grid(row=len(filas) + 2, column=0, columnspan=4, sticky="w",
                     pady=(7, 0))
        ttk.Button(botones, text="Aplicar y recargar",
                   command=self.aplicar).pack(side="left")
        ttk.Button(botones, text="Volver a los de por defecto",
                   command=self.por_defecto).pack(side="left", padx=6)
        ttk.Button(botones, text="Quitar la simulación (solo Klein-Nishina)",
                   command=self.quitar_sim).pack(side="left", padx=6)
        # La plantilla Geant4 ES la fuente: si el .root simula otros gammas, el
        # forward folding devuelve un D equivocado sin quejarse. Por eso se
        # bloquea, y destrabarlo tiene que ser un acto deliberado.
        ttk.Checkbutton(botones, text="usar la simulación igual (aunque sea de "
                                      "otra fuente)",
                        variable=self.forzar).pack(side="left", padx=10)
        self.refrescar_etiquetas()

    # ---------------- la fuente ----------------
    def al_elegir_preset(self, _evt=None):
        nombre = self.nombre_fuente.get()
        if nombre in FUENTES:
            E = FUENTES[nombre]
            self.eg1.set(E[0])
            self.eg2.set(E[1] if len(E) > 1 else 0.0)
        self.refrescar_etiquetas()

    def leer_fuente(self):
        """
        La Fuente que describen los controles.

        Si el nombre es un preset pero las energías fueron editadas a mano, la
        fuente pasa a llamarse "manual": el nombre no puede mentir sobre las
        líneas, porque es lo que se guarda en las preferencias y lo que aparece
        en los títulos.
        """
        lineas = [float(self.eg1.get())]
        try:
            E2 = float(self.eg2.get())
        except (tk.TclError, ValueError):
            E2 = 0.0
        if E2 > 0:
            lineas.append(E2)
        nombre = self.nombre_fuente.get()
        fu = Fuente(nombre if nombre in FUENTES else "manual", lineas)
        if nombre in FUENTES and not fu.es_preset():
            fu = Fuente("manual", lineas)
        return fu

    def quitar_sim(self):
        """Deja el análisis en la ruta analítica: sin .root."""
        self.var["sim"].set("")
        self.refrescar_etiquetas()

    # ---------------- seleccion ----------------
    def elegir(self, clave, tipos):
        actual = self.var[clave].get().strip()
        base = os.path.dirname(rutas.absoluta(actual)) if actual else rutas.DATOS
        ruta = filedialog.askopenfilename(
            parent=self, title=f"Elegir archivo de {clave}",
            initialdir=base if os.path.isdir(base) else rutas.DATOS,
            filetypes=tipos)
        if not ruta:
            return
        # se guarda RELATIVA a la raiz del proyecto si esta adentro: asi la
        # preferencia sigue valiendo si la carpeta se mueve de lugar, y la
        # etiqueta no muestra media pantalla de ruta de OneDrive.
        self.var[clave].set(rutas.relativa(os.path.normpath(ruta)))
        if clave == "sim":
            self.inspeccionar_sim(ruta)
        self.refrescar_etiquetas()

    def inspeccionar_sim(self, ruta):
        """Llena el combo de TTree con lo que el .root realmente tiene dentro."""
        try:
            arboles = inspeccionar_root(ruta)
        except Exception as e:
            messagebox.showwarning("Simulación",
                                   f"No se pudo abrir el .root:\n{e}")
            return
        if not arboles:
            messagebox.showwarning("Simulación", "El .root no tiene ningún TTree.")
            return
        nombres = list(arboles)
        self.cb_arbol.configure(values=nombres)
        if self.arbol.get() not in nombres:
            # el que tenga las ramas del analisis; si ninguna, el primero
            mejor = next((n for n in nombres
                          if "Edep_keV" in arboles[n][0]
                          or "Lvis_keVee" in arboles[n][0]), nombres[0])
            self.arbol.set(mejor)

    def refrescar_etiquetas(self):
        for clave, var in self.var.items():
            ruta = var.get().strip()
            if clave == "sim" and not ruta:      # la simulacion es opcional
                self.etiq[clave].configure(
                    text="sin simulación → solo Klein-Nishina", foreground="#666")
                continue
            if clave == "fondo" and not ruta:    # el fondo tambien lo es
                self.etiq[clave].configure(
                    text="sin fondo → el «neto» es la muestra cruda",
                    foreground="#a06000")
                continue
            ok, txt = info_archivo(ruta)
            self.etiq[clave].configure(text=txt,
                                       foreground="#666" if ok else "crimson")
        try:
            self.lbl_fuente.configure(text=self.leer_fuente().descripcion())
        except Exception as e:
            self.lbl_fuente.configure(text=f"⚠ {e}")

    def por_defecto(self):
        self.var["muestra"].set(ARCHIVO_MUESTRA)
        self.var["fondo"].set(ARCHIVO_FONDO)
        self.var["sim"].set(ARCHIVO_SIM)
        self.arbol.set(ARBOL_SIM)
        self.nombre_fuente.set(FUENTE_DEF)
        self.al_elegir_preset()

    def aplicar(self):
        rutas = {k: v.get().strip().strip('"') for k, v in self.var.items()}
        # La simulacion puede quedar vacia (ruta analitica) y el FONDO tambien:
        # las corridas binarias del 1-sep-2026 no tienen uno propio, y prestarles
        # el de otra sesion mete una forma ajena en el espectro (el fondo depende
        # del punto de operacion). La muestra es lo unico obligatorio.
        faltan = ["muestra"] if (not rutas["muestra"]
                                 or not os.path.exists(rutas["muestra"])) else []
        for opcional in ("fondo", "sim"):
            if rutas[opcional] and not os.path.exists(rutas[opcional]):
                faltan.append(opcional)
        if faltan:
            messagebox.showerror(
                "Archivos", "No existe(n):\n  " + "\n  ".join(
                    f"{k}: {rutas[k] or '(vacío)'}" for k in faltan))
            return
        if rutas["fondo"] and rutas["muestra"] == rutas["fondo"]:
            messagebox.showerror("Archivos",
                                 "Muestra y fondo son el mismo archivo: el neto "
                                 "daría cero.")
            return
        try:
            fu = self.leer_fuente()
        except Exception as e:
            messagebox.showerror("Fuente", str(e))
            return
        cambio = self.app.datos.fijar_archivos(rutas["muestra"], rutas["fondo"],
                                               rutas["sim"],
                                               self.arbol.get().strip() or ARBOL_SIM)
        cambio |= self.app.datos.fijar_fuente(fu)
        forzar = bool(self.forzar.get())
        if forzar != self.app.datos.forzar_sim:
            self.app.datos.forzar_sim = forzar
            cambio = True
        self.nombre_fuente.set(fu.nombre)
        guardar_preferencias(self.app.datos)
        self.app.actualizar_titulo()
        self.app.propagar_fuente()
        if cambio:
            self.app.reiniciar_pestanas(excepto=self)
        self.parsear.set(True)                # el usuario pidio cargar: se carga
        self.refrescar_etiquetas()
        self.lanzar()

    # ---------------- lectura ----------------
    def lanzar(self):
        # Los tk.Var se leen ACA, en el hilo principal. Tkinter no es
        # thread-safe: `trabajo()` corre en otro hilo y ahi un .get() puede
        # reventar ("main thread is not in main loop") o leer basura.
        self._pedido = {"rango": float(self.rango.get()),
                        "parsear": bool(self.parsear.get()),
                        "arbol": self.arbol.get().strip() or ARBOL_SIM}
        super().lanzar()

    def trabajo(self):
        d = self.app.datos
        pedido = getattr(self, "_pedido", {"rango": 2500.0, "parsear": False,
                                           "arbol": d.arbol_sim})
        rep = {"archivos": [], "sim": {"hist": {}, "faltan": []}, "esp": None,
               "areas": {}, "arbol": d.arbol_sim}
        for etiqueta, ruta in (("Muestra", d.archivo_muestra),
                               ("Fondo", d.archivo_fondo),
                               ("Simulación", d.archivo_sim)):
            if etiqueta == "Simulación" and not ruta:
                rep["archivos"].append(
                    (etiqueta, "(ninguna)", True,
                     "ruta analítica: Klein-Nishina, sin forward folding"))
                continue
            rep["archivos"].append((etiqueta, ruta) + info_archivo(ruta))
        for clave, ruta in (("muestra", d.archivo_muestra),
                            ("fondo", d.archivo_fondo)):
            if not ruta:                       # el fondo es opcional
                rep.setdefault("canal", {})[clave] = ("—", "")
                continue
            try:
                rep.setdefault("canal", {})[clave] = spectrum.detectar_canal(ruta)
            except Exception as e:
                rep.setdefault("canal", {})[clave] = ("?", f"{e}")

        rep["fuente"] = d.fuente
        rep["forzar"] = d.forzar_sim
        if d.hay_simulacion():
            self.progreso("gammas primarios de la simulación …", 0.02)
            rep["lineas_sim"] = d.lineas_sim()
            rep["compat"] = d.compat_simulacion()
        if not d.hay_simulacion():
            rep["sim"]["sin_sim"] = True          # ruta analítica, sin .root
        else:
            self.progreso("leyendo la simulación …", 0.05)
            try:
                arboles = inspeccionar_root(d.archivo_sim)
                rep["sim"]["arboles"] = arboles
                ramas = arboles.get(d.arbol_sim, ([], 0))[0]
                rep["sim"]["n_sim"] = arboles.get(d.arbol_sim, ([], 0))[1]
                rep["sim"]["faltan"] = [r for r in self.RAMAS_ESPERADAS
                                        if r not in ramas]
            except Exception as e:
                rep["sim"]["error"] = f"{type(e).__name__}: {e}"
            for rama in ("Lvis_keVee", "Edep_keV"):
                try:
                    rep["sim"]["hist"][rama] = d.simulacion(rama, "todo",
                                                            bins=350)
                except Exception as e:
                    rep["sim"]["hist"][rama] = None
                    rep["sim"].setdefault("errores_rama", []).append(
                        f"{rama}: {type(e).__name__}")

        entradas = [r for r in (d.archivo_muestra, d.archivo_fondo) if r]
        hay_cache = all(os.path.exists(r + ".areas_gate.npz") for r in entradas)
        if hay_cache or pedido["parsear"]:
            # el parseo se lleva casi todo el tiempo: le damos el tramo 10-90%
            # de la barra y adentro `Datos.espectro` reporta su propio avance
            # en bytes leidos (spectrum._iter_eventos).
            self.fase(0.10, 0.90)
            self.progreso("leyendo los datos …" if hay_cache else
                          "parseando los datos (la primera vez tarda ~40 s) …",
                          0.0)
            rep["esp"] = d.espectro(pedido["rango"], 500, 0.05, "gate")
            self.fase(0.0, 1.0)
            self.progreso("percentiles y tiempo vivo …", 0.92)
            for clave, ruta in (("muestra", d.archivo_muestra),
                                ("fondo", d.archivo_fondo)):
                if not ruta:                   # el fondo es opcional
                    continue
                areas, dur = spectrum.leer_areas_y_tiempo(ruta, 0.05,
                                                          metodo="gate",
                                                          verbose=False)
                pct = (np.percentile(areas, [50, 99, 99.9])
                       if len(areas) else np.zeros(3))
                rep["areas"][clave] = {"n": len(areas), "dur": dur, "pct": pct,
                                       "max": float(areas.max()) if len(areas) else 0.0}
        else:
            rep["falta_cache"] = True
        return rep

    # ---------------- dibujo ----------------
    def pintar(self, rep):
        ax1, ax2 = self.fig.subplots(1, 2)

        r = rep["esp"]
        if r is None:
            ax1.text(0.5, 0.5, "Sin caché de áreas.\nTildá «parsear los XML»\n"
                               "y apretá «Recalcular».",
                     ha="center", va="center", color="crimson")
            ax1.set_axis_off()
        else:
            c = r["canales"]
            ax1.step(c, r["muestra"], where="mid", color="tab:blue", alpha=0.45,
                     lw=1, label="muestra")
            ax1.step(c, r["fondo"], where="mid", color="gray", alpha=0.6, lw=1,
                     label="fondo (escalado)")
            ax1.step(c, r["neto"], where="mid", color="crimson", lw=1.3,
                     label="neto")
            ax1.set_xlabel("Área del pulso [mV·ns]")
            ax1.set_ylabel("Cuentas normalizadas")
            ax1.set_title("Vista previa del espectro medido", fontsize=10)
            ax1.legend(fontsize=8)
            ax1.grid(alpha=0.3)

        fu = rep.get("fuente") or self.fuente()
        hist = rep["sim"]["hist"]
        dibujado = False
        for rama, color in (("Lvis_keVee", "tab:blue"), ("Edep_keV", "tab:green")):
            h = hist.get(rama)
            if h is None:
                continue
            cen, pes = h
            ax2.step(cen, pes, where="mid", color=color, lw=1.1, label=rama)
            dibujado = True
        if dibujado:
            titulo = "Simulación (líneas: bordes Compton de %s)" % fu.nombre
        elif rep["sim"].get("sin_sim"):
            # Sin .root, el panel muestra la PLANTILLA ANALITICA: es literalmente
            # lo que se va a ajustar por la ruta Klein-Nishina, asi que sirve
            # igual de referencia visual.
            E = np.linspace(1.0, 1.25 * max(fu.bordes), 2000)
            for Eg, color in zip(fu.lineas, ("tab:blue", "tab:green")):
                ax2.plot(E, dsigma_dE(E, Eg), color=color, lw=1.6,
                         label="Klein-Nishina  E_γ = %.1f keV" % Eg)
            titulo = "Sin simulación: plantilla analítica (%s)" % fu.nombre
            dibujado = True
        else:
            titulo = None
        if dibujado:
            for E in fu.bordes:
                ax2.axvline(E, color="0.4", ls="--", lw=0.9)
            ax2.set_xlabel("Energía [keV]")
            ax2.set_ylabel("Eventos simulados" if hist else "dσ/dE_e (u.a.)")
            ax2.set_title(titulo, fontsize=10)
            ax2.legend(fontsize=8)
            ax2.grid(alpha=0.3)
        else:
            ax2.text(0.5, 0.5, "No se pudo leer la simulación\n(ver panel derecho)",
                     ha="center", va="center", color="crimson")
            ax2.set_axis_off()

        self.escribir(self._informe(rep) + "\n\n" + self.explicacion.strip())

    def _informe(self, rep):
        L = ["ARCHIVOS EN USO\n"]
        for etiqueta, ruta, ok, txt in rep["archivos"]:
            L.append(f"{etiqueta}: {os.path.basename(ruta)}")
            L.append(f"   {txt}")
        fu = rep.get("fuente") or self.fuente()
        L.append("")
        L.append("FUENTE: " + fu.descripcion())
        L.append("   modelo que aplica: "
                 + ("dos bordes solapados (Ec. 11-13)" if fu.doble
                    else "un solo borde (Ec. 8-10)"))
        L.append("   ancla del mapa por el origen: E_ref = %.1f keV" % fu.E_ref)
        if fu.nota():
            L.append("   " + fu.nota())
        for a in fu.avisos():
            L.append("   ⚠ " + a)
        s = rep["sim"]
        L.append("")
        if s.get("sin_sim"):
            L.append("SIMULACIÓN: ninguna cargada → ruta analítica "
                     "(Klein-Nishina). Las pestañas 7, 8 y 10 quedan "
                     "desactivadas; 4 y 6 funcionan igual.")
        elif "error" in s:
            L.append(f"SIMULACIÓN: ERROR — {s['error']}")
        else:
            L.append(f"SIMULACIÓN: TTree «{rep['arbol']}», "
                     f"{s.get('n_sim', 0):,} eventos")
            if s.get("faltan"):
                L.append("   ⚠ faltan ramas: " + ", ".join(s["faltan"]))
            otros = [n for n in s.get("arboles", {}) if n != rep["arbol"]]
            if otros:
                L.append("   otros TTree: " + ", ".join(otros))
            lin = rep.get("lineas_sim")
            if lin:
                L.append("   gammas primarios: "
                         + ", ".join("%.1f keV (%.0f%%)" % (E, 100 * f)
                                     for E, f in lin))
            ok, msg = rep.get("compat", (True, ""))
            if ok:
                L.append("   ✓ " + msg)
            elif rep.get("forzar"):
                L.append("   ⚠ FORZADA: " + msg)
            else:
                L.append("   ⛔ " + msg)
                L.append("   → las pestañas 7, 8 y 10 quedan bloqueadas; 4 y 6 "
                         "calibran igual con Klein-Nishina.")
        for e in s.get("errores_rama", []):
            L.append(f"   ⚠ {e}")

        r = rep["esp"]
        L.append("")
        if r is None:
            L.append("ESPECTRO: no leído (falta el caché de áreas).")
            L.append("Tildá «parsear los XML» y apretá «Recalcular».")
            return "\n".join(L)
        can = rep.get("canal", {})
        if can:
            cm, cf = can.get("muestra", ("?", ""))[0], can.get("fondo", ("?", ""))[0]
            L.append(f"CANAL DEL DRSOSC: muestra {cm}, fondo {cf}"
                     + ("" if (cm == cf or cf == "—")
                        else "   ⚠ ¡son distintos!"))
            L.append("")
        L.append("ESPECTRO NETO (área con compuerta, tiempo vivo)")
        for clave, etiqueta in (("muestra", "Muestra"), ("fondo", "Fondo")):
            a = rep["areas"].get(clave)
            tv = r["t_vivo_muestra"] if clave == "muestra" else r["t_vivo_fondo"]
            tp = r["t_muestra"] if clave == "muestra" else r["t_fondo"]
            if a is None:
                continue
            L.append(f"   {etiqueta}: {a['n']:,} pulsos · "
                     f"t_pared {tp:.0f} s · t_vivo {tv:.0f} s "
                     f"({100*tv/max(tp,1e-9):.0f}%)")
            L.append(f"      área  mediana {a['pct'][0]:.0f} · "
                     f"p99 {a['pct'][1]:.0f} · p99.9 {a['pct'][2]:.0f} · "
                     f"máx {a['max']:.0f} mV·ns")
        L.append(f"   k_muestra = {r['k_m']:.3f}   k_fondo = {r['k_f']:.3f}"
                 f"   (peso con que entra cada corrida)")
        nd = r["n_descartados"]
        L.append(f"   descartados por calidad: {nd[0]:,} de la muestra, "
                 f"{nd[1]:,} del fondo")
        p99 = rep["areas"].get("muestra", {}).get("pct", [0, 0, 0])[1]
        if p99:
            L.append(f"   → RANGO_MAX sugerido ≈ {50*np.ceil(1.6*p99/50):.0f}"
                     f" mV·ns (1.6·p99: deja aire ARRIBA del borde, §7.6, y no "
                     f"se lo come la cola de pile-up, §7.7)")
        return "\n".join(L)


# ======================================================================
#  Etapa 0 - espectro medido
# ======================================================================
class TabEspectro(Pestana):
    titulo = "0 · Espectro"
    explicacion = """
ETAPA 0 — de los XML crudos al espectro neto.

Azul: la muestra de Co-60. Gris: el fondo, escalado al mismo tiempo VIVO.
Rojo: el NETO (muestra - fondo), que es lo que se ajusta.

QUE MIRAR
• El plateau Compton entre ~600 y ~1200 mV·ns y la caída que muere en ~1600:
  ESO es el borde. No hay fotopico porque el plástico es de bajo Z.
• Los dos bordes caen en 1229 y 1427 mV·ns (D = 1.276 mV·ns/keVee).
• El pico enorme abajo de ~150 mV·ns NO es Compton: es el umbral de trigger.
• El fondo aporta poco: con tiempo vivo entra con k = 0.287, no 0.709.

PRUEBA A ROMPERLO — el hallazgo más caro del proyecto
Cambiá "área" de `gate` a `legacy` y subí RANGO_MAX a 8000. Eso reproduce el
análisis viejo, con la línea base promediada sobre los primeros 40 ns cuando el
pulso arranca a los 25 (§7.33): el espectro se estira 3.4x y la caída del borde
pasa de ~300 a ~2000 mV·ns de ancho. Ese ensanchamiento era el 43% de FWHM que
reportábamos como resolución del detector. No lo era.

También: poné RANGO_MAX en 46000 para ver el problema de binning de §7.7.
"""

    def construir_controles(self, frame):
        self.rango = tk.DoubleVar(value=2500.0)
        self.nbins = tk.IntVar(value=500)
        self.umbral = tk.DoubleVar(value=0.05)
        self.log = tk.BooleanVar(value=False)
        self.metodo = tk.StringVar(value="gate")
        for etiq, var, ancho in [("RANGO_MAX [mV·ns]", self.rango, 8),
                                 ("bins", self.nbins, 6),
                                 ("umbral [mV]", self.umbral, 6)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=ancho).pack(side="left")
        ttk.Label(frame, text="área").pack(side="left", padx=(10, 2))
        ttk.Combobox(frame, textvariable=self.metodo, width=8, state="readonly",
                     values=["gate", "legacy"]).pack(side="left")
        ttk.Checkbutton(frame, text="log", variable=self.log,
                        command=self.lanzar).pack(side="left", padx=10)

    def trabajo(self):
        return self.app.datos.espectro(self.rango.get(), self.nbins.get(),
                                       self.umbral.get(), self.metodo.get())

    def pintar(self, r):
        ax = self.fig.add_subplot(111)
        c = r["canales"]
        ax.step(c, r["muestra"], where="mid", color="tab:blue", alpha=0.45, lw=1,
                label="muestra (%s)" % self.fuente().nombre)
        ax.step(c, r["fondo"], where="mid", color="gray", alpha=0.6, lw=1,
                label="fondo (escalado)")
        ax.step(c, r["neto"], where="mid", color="crimson", lw=1.3, label="neto")
        if self.log.get():
            ax.set_yscale("log")
        ax.set_xlabel("Área del pulso [mV·ns]")
        ax.set_ylabel("Cuentas normalizadas")
        ax.set_title(f"Espectro neto ({r.get('metodo','gate')}) — bins de "
                     f"{c[1]-c[0]:.1f} mV·ns, {r['neto'].sum():.0f} cuentas netas, "
                     f"t_vivo {r.get('t_vivo_muestra', 0):.0f} s")
        ax.legend(); ax.grid(alpha=0.3)


# ======================================================================
#  Etapa 1 - Klein-Nishina
# ======================================================================
class TabFisica(Pestana):
    titulo = "1 · Klein-Nishina"
    explicacion = """
ETAPA 1 — la única teoría que entra: dσ/dE_e (Ec. 10).

Es el espectro de electrones de retroceso Compton para cada línea gamma. Un
plateau, una acumulación afilada, y un ESCALÓN en E_max. Sin ajustar nada.

E_max = 2α·E_γ/(1+2α),  α = E_γ/m_ec²

⚠️ La Ec. 1 impresa del paper está MAL: da E_γ/(1+2α), que es la energía del
fotón retrodispersado, no la del electrón. La de acá es la correcta y reproduce
los valores del propio paper (963 y 1118 keV para Co-60).

QUE MIRAR
• Las dos líneas del Co-60 dan bordes en 963.4 y 1118.1 keV: separados solo
  16%. Esa cercanía es la raíz de casi todos los problemas posteriores.
• El escalón es VERTICAL. Todo lo que en los datos reales no es vertical viene
  de la resolución del detector y del scattering múltiple.
"""

    def construir_controles(self, frame):
        fu = self.fuente()
        self.eg1 = tk.DoubleVar(value=fu.lineas[0])
        self.eg2 = tk.DoubleVar(value=fu.lineas[1] if fu.doble else 0.0)
        for etiq, var in [("E_γ 1 [keV]", self.eg1), ("E_γ 2 [keV]", self.eg2)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=8).pack(side="left")
        ttk.Label(frame, text="(0 = sin segunda línea)",
                  foreground="#666").pack(side="left", padx=6)

    def al_cambiar_fuente(self, fuente):
        self.eg1.set(fuente.lineas[0])
        self.eg2.set(fuente.lineas[1] if fuente.doble else 0.0)

    def trabajo(self):
        egs = [float(self.eg1.get())]
        try:
            E2 = float(self.eg2.get())
        except (tk.TclError, ValueError):
            E2 = 0.0
        if E2 > 0:
            egs.append(E2)
        return egs

    def pintar(self, egs):
        ax = self.fig.add_subplot(111)
        E = np.linspace(0, max(egs) * 1.15, 4000)
        for Eg, col in zip(egs, ("tab:blue", "tab:green")):
            Emax = compton_edge_energy(Eg)
            ax.plot(E, dsigma_dE(E, Eg), color=col, lw=1.8,
                    label=f"E_γ = {Eg:.1f} keV  →  borde {Emax:.1f} keV")
            ax.axvline(Emax, color=col, ls="--", alpha=0.6)
        ax.set_xlabel("Energía del electrón de retroceso  E_e  [keV]")
        ax.set_ylabel("dσ/dE_e  (forma, u.a.)")
        if len(egs) > 1:
            sep = 100 * (compton_edge_energy(egs[1]) /
                         compton_edge_energy(egs[0]) - 1)
            ax.set_title("Lineshape teórica — los bordes están separados "
                         f"{sep:.1f}%")
        else:
            ax.set_title("Lineshape teórica — una sola línea: un solo borde, "
                         "sin nada que desenredar")
        ax.legend(); ax.grid(alpha=0.3)


# ======================================================================
#  Etapa 2 - forward model
# ======================================================================
class TabForward(Pestana):
    titulo = "2 · Forward model"
    explicacion = """
ETAPA 2 — la gaussiana derrite el escalón (Ec. 8-9).

La lineshape teórica se lleva a espacio de canal (parámetro B = canal del
borde) y se convoluciona con una gaussiana de ancho σ: eso es lo que el
detector realmente mide.

QUE MIRAR
• Subí σ y mirá cómo el escalón vertical se vuelve una rampa. El kernel está
  NORMALIZADO (suma 1), así que σ ensancha sin cambiar el área: A y σ quedan
  desacoplados.
• La línea roja es el método clásico del 50% del máximo. Fijate que NO cae en
  B: ese es el sesgo que este método corrige.
• Con σ/B ≈ 0.18 (lo que pide tu detector) el borde casi desaparece. Por eso
  hace falta ajustar la forma completa y no buscar "la bajada".

PAD DE LA CONVOLUCIÓN (§7.14)
Destildá "pad" y mirá el borde izquierdo: sin extender la grilla, np.convolve
rellena con ceros y el modelo se cae. Eso hacía que el ajuste subestimara σ
un 18%.
"""

    def construir_controles(self, frame):
        self.sigma = tk.DoubleVar(value=40.0)
        self.edge = tk.DoubleVar(value=600.0)
        self.pad = tk.BooleanVar(value=True)
        ttk.Label(frame, text="B [canal]").pack(side="left", padx=(8, 2))
        ttk.Entry(frame, textvariable=self.edge, width=7).pack(side="left")
        ttk.Label(frame, text="σ").pack(side="left", padx=(10, 2))
        esc = ttk.Scale(frame, from_=0, to=200, variable=self.sigma,
                        orient="horizontal", length=200,
                        command=lambda _=None: self.lanzar())
        esc.pack(side="left")
        self.lbl = ttk.Label(frame, text="40")
        self.lbl.pack(side="left", padx=4)
        ttk.Checkbutton(frame, text="pad", variable=self.pad,
                        command=self.lanzar).pack(side="left", padx=10)

    def trabajo(self):
        return float(self.edge.get()), float(self.sigma.get()), bool(self.pad.get())

    def pintar(self, args):
        B, s, pad = args
        self.lbl.configure(text=f"{s:.0f}")
        ax = self.fig.add_subplot(111)
        can = np.linspace(0, 1000, 1024)
        Eg = self.fuente().lineas[0]
        ax.plot(can, forward_model(can, Eg, 1.0, B, 0.0), color="0.5", lw=1.2,
                label="ideal (σ=0)")
        y = forward_model(can, Eg, 1.0, B, s, pad=pad)
        ax.plot(can, y, "tab:blue", lw=2.2, label=f"σ = {s:.0f} canales")
        ax.axvline(B, color="k", ls="--", label=f"B real = {B:.0f}")
        if s > 0:
            ch50 = punto_mitad_altura(can, y)
            ax.axvline(ch50, color="crimson", ls=":", lw=2,
                       label=f"50% del máx = {ch50:.0f}  (sesgo {ch50-B:+.0f})")
        ax.set_xlabel("Canal"); ax.set_ylabel("S_model (u.a.)")
        ax.set_title(f"Forward model — σ/B = {s/B:.3f}"
                     f"{'' if pad else '   [SIN pad: mirá el borde izquierdo]'}")
        ax.legend(); ax.grid(alpha=0.3)


# ======================================================================
#  Etapa 3 - ajuste de un borde (y su fragilidad)
# ======================================================================
class TabUnBorde(Pestana):
    titulo = "3 · Un borde (frágil)"
    explicacion = """
ETAPA 3 — ajuste local de UN borde, y por qué no alcanza.

Sobre datos SINTÉTICOS con verdad conocida (B=600, σ=40) se ajusta el forward
model con least_squares desde el punto inicial que vos elijas.

QUE MIRAR
• Con un arranque razonable, recupera B casi exacto.
• Movés B₀ a 420 y σ₀ a 120 y el ajuste se va a B≈343: un error de -257 sobre
  600. least_squares es LOCAL, baja siempre, y se queda en el primer mínimo.
  El mínimo global sí está en 600, pero desde ahí no lo ve.
• Esa fragilidad es la motivación de las Etapas 4 (NCC, para sembrarlo bien) y
  5 (SA, para no depender del arranque).

La curva gris de fondo es el paisaje de la función de costo contra B: ahí se
ve el mínimo global y los locales donde se queda pegado.
"""

    def construir_controles(self, frame):
        self.b0 = tk.DoubleVar(value=420.0)
        self.s0 = tk.DoubleVar(value=120.0)
        for etiq, var in [("B₀ inicial", self.b0), ("σ₀ inicial", self.s0)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=7).pack(side="left")

    def trabajo(self):
        Eg = self.fuente().lineas[0]
        can = np.linspace(0, 1000, 1024)
        y = _espectro_sintetico(can, [(Eg, 5000.0, 600.0, 40.0)],
                                semilla=1).astype(float)
        iv = (270.0, 762.0)
        m = (can >= iv[0]) & (can <= iv[1])
        p0 = [2500.0, float(self.b0.get()), float(self.s0.get())]
        fit = ajustar_borde_unico(can, y, Eg, iv, p0=p0)

        # paisaje del costo contra B (con A y sigma reoptimizados groseramente)
        Bs = np.linspace(300, 750, 90)
        costo = []
        for B in Bs:
            mejor = np.inf
            for s in np.linspace(10, 150, 25):
                f = forward_model(can[m], Eg, 1.0, B, s)
                den = float(f @ f)
                A = (float(f @ y[m]) / den) if den > 0 else 0.0
                r = A * f - y[m]
                mejor = min(mejor, float(r @ r))
            costo.append(mejor)
        return can, y, iv, p0, fit, Bs, np.array(costo), Eg

    def pintar(self, res):
        can, y, iv, p0, fit, Bs, costo, Eg = res
        ax1 = self.fig.add_subplot(121)
        ax2 = self.fig.add_subplot(122)

        ax1.step(can, y, where="mid", color="0.75", lw=1, label="datos sintéticos")
        xf = np.linspace(iv[0], iv[1], 500)
        ax1.plot(xf, forward_model(xf, Eg, *p0[:1], p0[1], p0[2]), color="orange",
                 ls=":", lw=1.8, label=f"arranque (B₀={p0[1]:.0f})")
        ax1.plot(xf, forward_model(xf, Eg, fit["A"], fit["edge_ch"], fit["sigma"]),
                 "tab:blue", lw=2.2, label=f"ajuste (B={fit['edge_ch']:.0f})")
        ax1.axvline(600, color="k", ls="--", label="B real = 600")
        ax1.set_xlim(200, 900)
        ax1.set_xlabel("Canal"); ax1.set_ylabel("Cuentas")
        err = fit["edge_ch"] - 600
        ax1.set_title(f"Ajuste local — error {err:+.0f} canales"
                      f"{'  ✗' if abs(err) > 15 else '  ✓'}")
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

        ax2.semilogy(Bs, costo, color="0.4", lw=1.6)
        ax2.axvline(600, color="k", ls="--", label="B real (mínimo global)")
        ax2.axvline(p0[1], color="orange", ls=":", lw=2, label="arranque")
        ax2.axvline(fit["edge_ch"], color="tab:blue", ls=":", lw=2, label="donde terminó")
        ax2.set_xlabel("B [canal]"); ax2.set_ylabel("costo mínimo (log)")
        ax2.set_title("Paisaje de la función de costo")
        ax2.legend(fontsize=8); ax2.grid(alpha=0.3)


# ======================================================================
#  Etapa 4 - NCC
# ======================================================================
class TabNCC(Pestana):
    titulo = "4 · NCC"
    explicacion = """
ETAPA 4 — el intervalo lo elige el método, no vos (Ec. 15).

Se compara la FORMA del espectro medido contra la plantilla teórica, con el
coeficiente de Pearson. Ambas señales se reinterpolan a 100 puntos, así que ρ
es inmune a offset, a amplitud y a escala horizontal: solo compara forma.

QUE MIRAR
• Panel izq.: la zona verde es el intervalo elegido AUTOMÁTICAMENTE.
• Panel der. arriba: S1 (datos) contra S2 (teoría), ambos estandarizados. Eso
  es literalmente lo que Pearson compara.
• Panel der. abajo: ρ contra σ_rel. Tiene un máximo, así que el NCC también
  estima la resolución. En los datos reales el máximo cae en σ_rel ≈ 0.18.

MODO ORIGEN vs LIBRE (§7.10)
"libre" es el barrido 2D del paper. Probalo: sobre datos reales elige ventanas
degeneradas porque el barrido tiene un grado de libertad que el modelo (Ec. 8,
mapa por el origen) no tiene. "origen" barre solo B y es el default.

PLANTILLA (§7.28) [fuera del paper]
Con "Geant4" la referencia deja de ser Klein-Nishina y pasa a ser TU simulación:
ya trae los dos bordes, el scattering múltiple y el quenching. Si la ventana se
va a usar para comparar datos contra simulación, esta es la que corresponde.
Ahí B ya no es "el borde" sino el ancla del mapa por el origen: D₀ = B/963.4.

n_puntos (§7.27)
"auto" = un punto por bin, sin inventar puntos. Mirá el diagnóstico "ρ y B vs
n_puntos": de 100 a 571 puntos B se mueve 0.4% y se estabiliza. Interpolar a
POCOS puntos no mejora nada, sólo INFLA ρ promediando ruido.

OJO: acá inicio_rel = 0.45, el valor del paper. Por eso el D₀ que sale es ~1%
más bajo que el de la pestaña 8: la ventana entra en la zona de bajo canal
contaminada (§7.32) y el NCC no tiene cómo darse cuenta. La pestaña 8 usa 0.62.
"""

    PLANTILLAS = {"Klein-Nishina": None,
                  "Geant4 Lvis": ("Lvis_keVee", "todo"),
                  "Geant4 Edep": ("Edep_keV", "todo")}

    def construir_controles(self, frame):
        self.modo = tk.StringVar(value="origen")
        # ojo: `fuente` (la radiactiva) es de la app; esto elige de donde salen
        # los datos, reales o sinteticos.
        self.origen_datos = tk.StringVar(value="reales")
        self.plantilla = tk.StringVar(value="Klein-Nishina")
        self.npuntos = tk.StringVar(value="100")
        self.diag = tk.StringVar(value="ρ vs σ_rel")
        for etiq, var, vals, ancho in [
                ("datos", self.origen_datos, ["reales", "sintéticos"], 10),
                ("modo", self.modo, ["origen", "libre"], 8),
                ("plantilla", self.plantilla, list(self.PLANTILLAS), 13),
                ("n_puntos", self.npuntos, ["50", "100", "200", "400", "auto"], 6),
                ("diagnóstico", self.diag,
                 ["ρ vs σ_rel", "ρ y B vs n_puntos"], 17)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Combobox(frame, textvariable=var, width=ancho, state="readonly",
                         values=vals).pack(side="left")

    def _buscar(self, can, suave, neto, region, **kw):
        """Un barrido, con la plantilla que este elegida en el combo."""
        fu = self.fuente()
        nombre = self.plantilla.get()
        comun = dict(region=region, modo=self.modo.get(), cuentas_crudas=neto,
                     min_cuentas=200.0)
        comun.update(kw)
        if self.PLANTILLAS[nombre] is None:
            return buscar_intervalo(can, suave, fu.lineas[0], **comun)
        rama, filtro = self.PLANTILLAS[nombre]
        cen, pes = self.app.datos.simulacion(rama, filtro, bins=175)
        return buscar_intervalo_simulado(can, suave, cen, pes,
                                         E_gamma=fu.lineas[0], E_ref=fu.E_ref,
                                         E_ref2=fu.E_ref2, **comun)

    def trabajo(self):
        if self.PLANTILLAS[self.plantilla.get()] is not None:
            self.exigir_simulacion()
        if self.origen_datos.get() == "reales":
            r = self.app.datos.espectro()
            can, neto = r["canales"], r["neto"]
            suave = spectrum.suavizar(neto, ventana=21, orden=3)
            region = (250.0, 2500.0)
            verdad = None
        else:
            can = np.linspace(0, 1000, 1024)
            neto = _espectro_sintetico(can,
                                       [(self.fuente().lineas[0], 5000.0,
                                         600.0, 40.0)],
                                       semilla=1).astype(float)
            suave = neto
            region = None
            verdad = 600.0

        n_txt = self.npuntos.get()
        n_pts = "auto" if n_txt == "auto" else int(n_txt)
        self.progreso("barrido NCC ...", 0.5)
        ncc = self._buscar(can, suave, neto, region, paso=1, n_puntos=n_pts)

        # --- diagnostico de la derecha abajo ---
        if self.diag.get().startswith("ρ vs"):
            grilla = []
            for s in (0.03, 0.06, 0.09, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25):
                rr = self._buscar(can, suave, neto, region, sigma_rel=s, paso=2,
                                  n_puntos=n_pts)
                grilla.append((s, rr["rho"], rr["edge_ch"]))
            eje_diag = "σ_rel"
        else:
            grilla = []
            for n in (25, 50, 100, 200, 400, 800, "auto"):
                rr = self._buscar(can, suave, neto, region, paso=1, n_puntos=n)
                grilla.append((rr["n_puntos"], rr["rho"], rr["edge_ch"]))
            grilla.sort()
            eje_diag = "n_puntos"

        # la plantilla ganadora, para el panel de Pearson
        if self.PLANTILLAS[self.plantilla.get()] is None:
            _, S2 = plantilla_referencia(self.fuente().lineas[0],
                                         sigma_rel=ncc["sigma_rel"],
                                         n_puntos=ncc["n_puntos"],
                                         inicio_rel=ncc["inicio_rel"],
                                         fin_rel=ncc["fin_rel"])
        else:
            rama, filtro = self.PLANTILLAS[self.plantilla.get()]
            cen, pes = self.app.datos.simulacion(rama, filtro, bins=175)
            _, S2 = plantilla_simulada(cen, pes, sigma_rel=ncc["sigma_rel"],
                                       n_puntos=ncc["n_puntos"],
                                       inicio_rel=ncc["inicio_rel"],
                                       fin_rel=ncc["fin_rel"],
                                       E_ref=ncc["E_ref"])
        return can, neto, suave, ncc, np.array(grilla), verdad, S2, eje_diag

    def pintar(self, res):
        can, neto, suave, ncc, grilla, verdad, S2, eje_diag = res
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1.5, 1])
        ax1 = self.fig.add_subplot(gs[:, 0])
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax3 = self.fig.add_subplot(gs[1, 1])

        ax1.step(can, neto, where="mid", color="0.75", lw=1, label="neto")
        ax1.plot(can, suave, color="crimson", lw=1.5, label="suavizado (para NCC)")
        ax1.axvspan(ncc["lo"], ncc["hi"], color="green", alpha=0.13,
                    label=f"intervalo NCC (ρ={ncc['rho']:.4f})")
        ax1.axvline(ncc["edge_ch"], color="green", ls=":", lw=2,
                    label=f"B = {ncc['edge_ch']:.0f}")
        if verdad:
            ax1.axvline(verdad, color="k", ls="--", label=f"B real = {verdad:.0f}")
        # escalar a la zona del borde: si no, el pico de ruido de bajo canal
        # (10x mas alto) aplasta todo lo que interesa
        dentro = (can >= ncc["lo"]) & (can <= ncc["hi"])
        if dentro.any():
            ax1.set_ylim(min(0, suave[dentro].min() * 1.2),
                         suave[dentro].max() * 1.6)
            ax1.set_xlim(max(0, ncc["lo"] - 0.35 * (ncc["hi"] - ncc["lo"])),
                         ncc["hi"] * 1.02)
        ax1.set_xlabel("Canal [mV·ns]"); ax1.set_ylabel("Cuentas")
        ax1.set_title(f"Intervalo automático — {self.plantilla.get()}, modo "
                      f"'{ncc['modo']}'\nσ_rel = {ncc['sigma_rel']:.3f},  "
                      f"n = {ncc['n_puntos']},  D₀ = {ncc['D_seed']:.3f} mV·ns/keV",
                      fontsize=9)
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

        n = ncc["n_puntos"]
        xq = np.linspace(ncc["lo"], ncc["hi"], n)
        S1 = np.interp(xq, can, suave)
        z = lambda v: (v - np.mean(v)) / np.std(v)
        ax2.plot(z(S1), color="0.35", lw=1.4, label="S1 (datos)")
        ax2.plot(z(S2), "tab:green", lw=2, label=f"S2 ({self.plantilla.get()})")
        ax2.set_title(f"Lo que compara Pearson: ρ = {ncc['rho']:.4f}", fontsize=10)
        ax2.set_xlabel(f"punto (1..{n})"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

        if eje_diag == "σ_rel":
            ax3.plot(grilla[:, 0], grilla[:, 1], "o-", color="tab:purple")
            ax3.axvline(ncc["sigma_rel"], color="green", ls=":", lw=2)
            ax3.set_xlabel("σ_rel de la plantilla"); ax3.set_ylabel("ρ máximo")
            ax3.set_title("El NCC también estima la resolución", fontsize=10)
        else:
            # la meseta de n_puntos (§7.27): rho a la izquierda, B a la derecha
            ax3.semilogx(grilla[:, 0], grilla[:, 1], "o-", color="tab:purple")
            ax3.set_xlabel("n_puntos"); ax3.set_ylabel("ρ", color="tab:purple")
            g = ax3.twinx()
            g.semilogx(grilla[:, 0], grilla[:, 2], "s--", color="tab:blue", ms=4)
            g.set_ylabel("B [mV·ns]", color="tab:blue")
            ax3.set_title("Más puntos NO agregan información:\nB se estabiliza y ρ "
                          "sólo baja lo que estaba inflado", fontsize=9)
        ax3.grid(alpha=0.3)


# ======================================================================
#  Etapa 6 - ajuste analitico: uno o dos bordes segun la fuente
# ======================================================================
class TabAnalitico(Pestana):
    titulo = "6 · Ajuste analítico"
    explicacion = """
ETAPA 6 — calibración con la plantilla TEÓRICA (Klein-Nishina), sin simulación.

Lo que ajusta depende de cuántas líneas gamma tenga la fuente elegida en la
pestaña «Archivos»:

CON DOS LÍNEAS (Co-60) — bordes solapados, Ec. 11-13
Θ = {A1, A2, B1, σ1, σ2, D}, con CE_ch,2 = B1 + D·ΔE.
⚠️ ΔE es la diferencia de los BORDES (154.7 keV en Co-60), no la de las líneas
gamma (159.3), como dice el paper: D convierte energía del ELECTRÓN a canal.
• Verde y magenta son las dos componentes. En Co-60 se solapan casi por
  completo: B2−B1 queda por debajo de 2σ1. Ahí está el problema de fondo.
• Rojo punteado: el ajuste de un solo borde. Queda corrido ~6%.

CON UNA LÍNEA (Cs-137) — un solo borde, Ec. 8-10
No hay nada que desenredar: D = B/E_max, directo. Es el caso FÁCIL, y por eso
una segunda fuente vale tanto: da un punto de calibración limpio (§9-bis.3).
• Gris punteado: el método clásico del 50% de altura, para comparar. Está
  sesgado porque el flanco no es simétrico.

LAS LIGADURAS (§7.16, §7.23) — sólo tienen sentido con dos líneas
• "origen": D = B1/E_max,1, o sea el mapa pasa por el origen (área 0 ↔ E 0).
  Sin ella D se va un 8% según la ventana.
• "σ ligadas": σ2 = σ1·√(E2/E1), que es como escala la resolución.
Destildalas y mirá el "offset en E=0" del título: si se va a miles de mV·ns, el
ajuste está usando libertad que la física no le da.

LA VENTANA
"NCC automática" la elige correlacionando contra la plantilla (Etapa 4), así que
se adapta sola a la escala del detector; en "manual" mandan lo/hi.

⚠️ EL SESGO CONOCIDO DE ESTE NÚMERO
La plantilla analítica es scattering SIMPLE. El múltiple mueve cuentas por
encima del borde y el ajuste lo compensa corriendo B hacia arriba: medido en el
Co-60 de este proyecto, +4.3% contra el modelo con simulación (§7.18). Si hay un
.root cargado, el número bueno es el de la pestaña 10.
"""

    def titulo_actual(self):
        return "6 · Dos bordes" if self.fuente().doble else "6 · Un borde"

    def construir_controles(self, frame):
        self.ventana = tk.StringVar(value="NCC automática")
        self.lo = tk.DoubleVar(value=600.0)
        self.hi = tk.DoubleVar(value=1900.0)
        self.origen = tk.BooleanVar(value=True)
        self.sigmas = tk.BooleanVar(value=False)
        ttk.Label(frame, text="ventana").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.ventana, width=15, state="readonly",
                     values=["NCC automática", "manual"]).pack(side="left")
        for etiq, var in [("lo", self.lo), ("hi", self.hi)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=7).pack(side="left")
        ttk.Checkbutton(frame, text="origen ligado", variable=self.origen,
                        command=self.lanzar).pack(side="left", padx=10)
        ttk.Checkbutton(frame, text="σ ligadas", variable=self.sigmas,
                        command=self.lanzar).pack(side="left", padx=4)

    def trabajo(self):
        fu = self.fuente()
        r = self.app.datos.espectro()
        can, neto, err = r["canales"], r["neto"], r["error"]
        self.progreso("NCC: ventana y semilla …", 0.2)
        ncc = self.app.datos.borde_estimado()
        iv = ((float(self.lo.get()), float(self.hi.get()))
              if self.ventana.get() == "manual" else (ncc["lo"], ncc["hi"]))
        m = (can >= iv[0]) & (can <= iv[1])
        B0 = float(ncc["edge_ch"])

        self.progreso("ajuste …", 0.7)
        uno = ajustar_borde_unico(can, neto, fu.lineas[0], iv, p0=ncc["p0"])
        # los tres metodos CLASICOS sobre la misma ventana, para comparar contra
        # lo que hace el paper (Fig. 4). No son el resultado: son la vara.
        self.progreso("métodos clásicos …")
        # la region se ancla al borde del NCC, NO a la ventana de ajuste: esta
        # ultima entra en el exceso de bajo canal y el metodo de la diferencia
        # se va ahi (ver el aviso en clasicos.comparar).
        clas = comparar_clasicos(can, neto, fu.E_ref, B_ref=B0, errores=err)
        salida = {"can": can, "neto": neto, "err": err, "iv": iv, "uno": uno,
                  "fuente": fu, "B_ncc": B0, "rho": float(ncc["rho"]),
                  "clasicos": clas,
                  "mitad": float(clas["50% de altura"]["borde"])}

        if fu.doble:
            # B1_rango sale del NCC, no de una constante: el 4357 que habia aca
            # era de la escala vieja del area (§7.33) y arrastraba el ajuste a
            # D = 1.90 en vez de 1.32.
            dob = ajustar_dos_bordes(can, neto, iv, E_gamma1=fu.lineas[0],
                                     E_gamma2=fu.lineas[1],
                                     ligar_origen=bool(self.origen.get()),
                                     ligar_sigmas=bool(self.sigmas.get()),
                                     B1_rango=(0.60 * B0, 1.25 * B0))
            f1, f2, _ = componentes(can[m], dob["B1"], dob["sigma1"],
                                    dob["sigma2"], dob["D"],
                                    fu.lineas[0], fu.lineas[1])
            rr = (dob["A1"] * f1 + dob["A2"] * f2 - neto[m]) / err[m]
            n_par = 6 - int(self.origen.get()) - int(self.sigmas.get())
            salida.update(dob=dob, chi2=float(rr @ rr),
                          ndf=max(int(m.sum()) - n_par, 1), D=dob["D"])
        else:
            mod = forward_model(can[m], fu.lineas[0], uno["A"], uno["edge_ch"],
                                uno["sigma"])
            rr = (mod - neto[m]) / err[m]
            salida.update(dob=None, chi2=float(rr @ rr),
                          ndf=max(int(m.sum()) - 3, 1),
                          D=uno["edge_ch"] / fu.E_ref)
        return salida

    def pintar(self, s):
        fu, can, neto, iv = s["fuente"], s["can"], s["neto"], s["iv"]
        ax = self.fig.add_subplot(111)
        z = (can >= max(iv[0] - 0.6 * (iv[1] - iv[0]), 0)) & (can <= iv[1] * 1.05)
        xf = np.linspace(iv[0], iv[1], 700)
        ax.step(can[z], neto[z], where="mid", color="0.75", lw=1, label="neto")
        ax.axvspan(iv[0], iv[1], color="green", alpha=0.05)
        uno = s["uno"]

        if s["dob"] is not None:
            dob = s["dob"]
            f1, f2, _ = componentes(xf, dob["B1"], dob["sigma1"], dob["sigma2"],
                                    dob["D"], fu.lineas[0], fu.lineas[1])
            ax.plot(xf, dob["A1"] * f1, "--", color="tab:green", lw=1.6,
                    label=f"f1 ({fu.bordes[0]:.0f} keV) — B1 = {dob['B1']:.0f}")
            ax.plot(xf, dob["A2"] * f2, "--", color="tab:purple", lw=1.6,
                    label=f"f2 ({fu.bordes[1]:.0f} keV) — B2 = {dob['B2']:.0f}")
            ax.plot(xf, dob["A1"] * f1 + dob["A2"] * f2, "tab:blue", lw=2.3,
                    label="S_comp (Ec. 13)")
            ax.plot(xf, forward_model(xf, fu.lineas[0], uno["A"], uno["edge_ch"],
                                      uno["sigma"]),
                    ":", color="crimson", lw=1.8,
                    label=f"un borde (B = {uno['edge_ch']:.0f}, sesgado)")
            titulo = (f"D = {s['D']:.3f} mV·ns/keV   |   "
                      f"χ²/ndf = {s['chi2']/s['ndf']:.3f}   |   "
                      f"offset en E=0 = {dob['offset_ch']:+.0f} mV·ns   |   "
                      f"separación = "
                      f"{(dob['B2']-dob['B1'])/dob['sigma1']:.2f}σ")
        else:
            ax.plot(xf, forward_model(xf, fu.lineas[0], uno["A"], uno["edge_ch"],
                                      uno["sigma"]),
                    "tab:blue", lw=2.3,
                    label=f"un borde — B = {uno['edge_ch']:.0f} mV·ns")
            ax.axvline(uno["edge_ch"], color="tab:blue", ls="--", alpha=0.6)
            for nombre, color, trazo in (
                    ("50% de altura", "0.45", ":"),
                    ("diferencia", "tab:olive", "-."),
                    ("FWHM gaussiano", "tab:brown", (0, (4, 2)))):
                b = s["clasicos"].get(nombre, {}).get("borde", float("nan"))
                if np.isfinite(b):
                    ax.axvline(b, color=color, ls=trazo, lw=1.5,
                               label=f"clásico · {nombre} = {b:.0f}")
            titulo = (f"D = {s['D']:.4f} mV·ns/keV   |   "
                      f"χ²/ndf = {s['chi2']/s['ndf']:.3f}   |   "
                      f"σ/B = {uno['sigma']/uno['edge_ch']:.3f}   |   "
                      f"FWHM = {235.5*uno['sigma']/uno['edge_ch']:.1f}%")
        ax.set_xlabel("Área del pulso [mV·ns]"); ax.set_ylabel("Cuentas")
        ax.set_title(titulo)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        L = [f"FUENTE: {fu.descripcion()}",
             f"ventana [{iv[0]:.0f}, {iv[1]:.0f}] mV·ns  "
             f"({self.ventana.get()}, ρ = {s['rho']:.4f})",
             f"D = {s['D']:.4f} mV·ns/keV   →  {1000/s['D']:.1f} keV por "
             f"1000 mV·ns",
             f"χ²/ndf = {s['chi2']/s['ndf']:.3f}  ({s['ndf']} gdl)"]
        if s["dob"] is None:
            L.append(f"B = {uno['edge_ch']:.1f} mV·ns  ↔  E_max = "
                     f"{fu.E_ref:.1f} keV")
        L.append("")
        L.append("MÉTODOS CLÁSICOS (la vara del paper; región anclada al "
                 "borde del NCC, no la ventana de ajuste)")
        for nombre, v in s["clasicos"].items():
            if not np.isfinite(v.get("D", float("nan"))):
                L.append(f"   {nombre}: no ajustó")
                continue
            extra = ""
            if "fwhm_rel" in v:
                extra = f"   FWHM = {100*v['fwhm_rel']:.1f}%"
            elif "ventana" in v:
                extra = f"   suavizado {v['ventana']} bins"
            L.append(f"   {nombre:<15} D = {v['D']:.4f} "
                     f"({100*(v['D']/s['D'] - 1):+.1f}%){extra}")
        L.append("")
        L.append(self.explicacion.strip())
        self.escribir("\n".join(L))


# ======================================================================
#  Etapa 5 - simulated annealing
# ======================================================================
class TabSA(Pestana):
    titulo = "5 · Recocido (SA)"
    auto = False
    explicacion = """
ETAPA 5 — optimizador GLOBAL (Sección II-A). Tarda ~20-60 s: apretá
"Recalcular".

Acepta a veces movimientos que EMPEORAN el costo, con probabilidad exp(-ΔE/T).
Esa es la puerta para salir de un mínimo local. Al enfriar deja de aceptarlos.

QUE MIRAR
• Panel izq.: la curva gris (costo actual) SUBE a veces. Eso es el algoritmo
  funcionando, no un error. La azul (mejor hasta ahora) solo baja.
• Panel centro: T baja geométricamente y la tasa de aceptación cae con ella.
• Panel der.: dónde terminó, contra el arranque y contra la verdad.

DOS COSAS QUE ENCONTRAMOS (§7.22)
• "T_init = 2000" del paper no es transferible: T se compara contra ΔE, así que
  depende de las UNIDADES del costo. Medido: 2.1e9 con la Ec. 14 y 3.19 con el
  error NCC. Por eso acá T_init se calibra automáticamente.
• Con el objetivo NCC la amplitud es INDETERMINABLE (Pearson la ignora).

Probá el objetivo "NCC" y el arranque malo: ambos llegan a B ≈ 600.
"""

    def construir_controles(self, frame):
        self.objetivo = tk.StringVar(value="mínimos cuadrados (Ec. 14)")
        self.semilla = tk.IntVar(value=1)
        self.b0 = tk.DoubleVar(value=420.0)
        ttk.Label(frame, text="objetivo").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.objetivo, width=26, state="readonly",
                     values=["mínimos cuadrados (Ec. 14)",
                             "NCC / 1-ρ (Metropolis)"]).pack(side="left")
        for etiq, var in [("B₀", self.b0), ("semilla", self.semilla)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(10, 2))
            ttk.Entry(frame, textvariable=var, width=6).pack(side="left")

    def trabajo(self):
        from etapas.sa import (recocido_simulado, costo_minimos_cuadrados, costo_ncc,
                        modelo_un_borde, proyectar_amplitud)
        can = np.linspace(0, 1000, 1024)
        Eg = self.fuente().lineas[0]
        y = _espectro_sintetico(can, [(Eg, 5000.0, 600.0, 40.0)],
                                semilla=1).astype(float)
        lo, hi = 270.0, 762.0
        m = (can >= lo) & (can <= hi)
        xv, yv = can[m], y[m]
        mod = modelo_un_borde(Eg)
        usa_ncc = self.objetivo.get().startswith("NCC")
        E = costo_ncc(xv, yv, mod) if usa_ncc else costo_minimos_cuadrados(xv, yv, mod)
        p0 = [2500.0, float(self.b0.get()), 120.0]
        cotas = ([0.0, lo, 1.0], [5e4, hi, hi - lo])
        r = recocido_simulado(E, p0, cotas, semilla=int(self.semilla.get()))
        x = proyectar_amplitud(xv, yv, mod, r["x"], 0) if usa_ncc else r["x"]
        loc = ajustar_borde_unico(can, y, Eg, (lo, hi), p0=p0)
        return can, y, (lo, hi), r, x, loc, p0, usa_ncc

    def pintar(self, res):
        can, y, iv, r, x, loc, p0, usa_ncc = res
        ax1 = self.fig.add_subplot(131)
        ax2 = self.fig.add_subplot(132)
        ax3 = self.fig.add_subplot(133)
        h = r["hist"]

        ax1.plot(h["actual"], color="0.65", lw=0.9, label="costo actual")
        ax1.plot(h["mejor"], "tab:blue", lw=1.8, label="mejor hasta ahora")
        ax1.set_yscale("log"); ax1.set_xlabel("ciclo de temperatura")
        ax1.set_ylabel("E(x)")
        ax1.set_title("El costo actual SUBE a veces\n(así escapa mínimos locales)",
                      fontsize=10)
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

        ax2.plot(h["T"], color="crimson", lw=1.8)
        ax2.set_yscale("log"); ax2.set_ylabel("T", color="crimson")
        ax2.set_xlabel("ciclo")
        g = ax2.twinx()
        g.plot(h["aceptacion"], color="steelblue", lw=1.2)
        g.set_ylabel("tasa de aceptación", color="steelblue")
        ax2.set_title(f"T_init calibrado = {r['T_inicial']:.3g}", fontsize=10)
        ax2.grid(alpha=0.3)

        xf = np.linspace(iv[0], iv[1], 500)
        ax3.step(can, y, where="mid", color="0.75", lw=1, label="datos")
        Eg = self.fuente().lineas[0]
        ax3.plot(xf, forward_model(xf, Eg, x[0], x[1], x[2]), "tab:blue", lw=2.2,
                 label=f"SA (B = {x[1]:.0f})")
        ax3.plot(xf, forward_model(xf, Eg, loc["A"], loc["edge_ch"], loc["sigma"]),
                 "--", color="crimson", lw=1.8,
                 label=f"local (B = {loc['edge_ch']:.0f})")
        ax3.axvline(600, color="k", ls="--", label="B real = 600")
        ax3.axvline(p0[1], color="orange", ls=":", label=f"arranque {p0[1]:.0f}")
        ax3.set_xlim(200, 900); ax3.set_xlabel("Canal")
        ax3.set_title(f"Desde el mismo arranque malo\n"
                      f"{'objetivo NCC' if usa_ncc else 'objetivo Ec. 14'}, "
                      f"{r['n_evals']} evals", fontsize=10)
        ax3.legend(fontsize=8); ax3.grid(alpha=0.3)


# ======================================================================
#  Etapa 7 - forward folding con Geant4
# ======================================================================
class TabSimulacion(Pestana):
    titulo = "7 · Geant4"
    auto = False
    explicacion = """
ETAPA 7 — forward folding con la simulación. Tarda ~20 s la primera vez.

En vez de la lineshape analítica se usa TU simulación Geant4 como plantilla:
    S(ch) = A · Σ_j h_j · G(ch; D·E_j, D·σ_E(E_j)),  σ_E = α·√E
Solo 3 parámetros, y la plantilla ya contiene el scattering múltiple, el escape
y el quenching de Birks.

QUE MIRAR
• Panel der. arriba: EL ESPECTRO SIMULADO, tal cual sale del Geant4 (852079
  eventos). Der. abajo: el zoom a la zona de los bordes.
• La roja (1 solo scattering) cae a cero de golpe en 963 y 1118: eso es
  Klein-Nishina. La verde (todos los eventos) sigue con cuentas POR ENCIMA del
  borde, hasta 1332 keV: ese relleno es el 25% de scattering múltiple, y es lo
  que el modelo analítico no puede describir.
• La azul (Lvis) tiene el escalón corrido hacia abajo respecto de la verde: son
  los ~4.6% del quenching de Birks.
• Panel izq.: la curva gris es la simulación llevada a canal SIN la resolución
  del detector. Comparala con las de color: todo ese derretimiento es σ.
• Las tres plegadas quedan casi encima de los datos (χ²/ndf ≈ 0.96 las tres)
  pero dan D distintos. Ahí está el error sistemático real.

EL RESULTADO (§7.18)
Klein-Nishina y "1 scattering" coinciden (4.34 y 4.32) porque Klein-Nishina ES
scattering simple. Con la simulación completa D baja a 4.05: el método del
paper sobreestima D ~6% en este detector.
"""

    def construir_controles(self, frame):
        self.ventana = tk.StringVar(value="auto (NCC)")
        self.lo = tk.DoubleVar(value=300.0)
        self.hi = tk.DoubleVar(value=1900.0)
        self.log = tk.BooleanVar(value=False)
        ttk.Label(frame, text="ventana").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.ventana, width=11, state="readonly",
                     values=["auto (NCC)", "manual"]).pack(side="left")
        for etiq, var in [("lo", self.lo), ("hi", self.hi)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=7).pack(side="left")
        ttk.Checkbutton(frame, text="log en el espectro simulado", variable=self.log,
                        command=self.repintar).pack(side="left", padx=12)

    def trabajo(self):
        self.exigir_simulacion()
        from ajuste.sim_template import ajustar_forward_folding, modelo_simulado
        fu = self.fuente()
        r = self.app.datos.espectro()
        can, neto, err = r["canales"], r["neto"], r["error"]
        iv = self.app.datos.ventana_ajuste(
            manual=(float(self.lo.get()), float(self.hi.get()))
            if self.ventana.get() == "manual" else None)
        salida = []
        for etiq, rama, filtro in [("Lvis (luz, con Birks)", "Lvis_keVee", "todo"),
                                   ("Edep (energía cruda)", "Edep_keV", "todo"),
                                   ("Edep, 1 solo scattering", "Edep_keV", "1 scattering")]:
            cen, pes = self.app.datos.simulacion(rama, filtro)
            fit = ajustar_forward_folding(can, neto, cen, pes, iv, errores=err,
                                          E_ref=fu.E_ref)
            fit["centros"], fit["pesos"], fit["etiqueta"] = cen, pes, etiq
            fit["curva"] = modelo_simulado(can, cen, pes, fit["A"], fit["D"],
                                           fit["alpha"], fit["beta"])
            # la MISMA simulacion llevada a canal pero SIN la resolucion del
            # detector: muestra cuanto trabajo hace la gaussiana.
            # Se ensancha exactamente UN BIN (alpha=0, beta=ancho_bin/D, o sea
            # sigma_ch = ancho_bin). Con sigma->0 de verdad sale un peine de
            # picos: los bins de la simulacion (2 keV ~ 8.5 mV*ns) no coinciden
            # con los del histograma (10 mV*ns) y aparece aliasing, que es un
            # artefacto del binning y no fisica.
            ancho_bin = float(can[1] - can[0])
            fit["curva_ideal"] = modelo_simulado(can, cen, pes, fit["A"], fit["D"],
                                                 0.0, ancho_bin / fit["D"])
            salida.append(fit)
        return can, neto, iv, salida, fu

    ESTILOS = [("tab:blue", "-", 2.4), ("tab:green", "--", 2.0),
               ("crimson", ":", 2.0)]

    def pintar(self, res):
        can, neto, iv, fits, fu = res
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1.35, 1])
        ax1 = self.fig.add_subplot(gs[:, 0])
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax3 = self.fig.add_subplot(gs[1, 1])

        # ---- (izq) datos + los tres modelos plegados ----
        # la ventana de dibujo se saca de la de AJUSTE: con 1200-7500 fijo (la
        # escala del area vieja) el panel quedaba vacio.
        z = (can >= 0.35 * iv[0]) & (can <= min(1.25 * iv[1], can[-1]))
        ax1.step(can[z], neto[z], where="mid", color="0.75", lw=1, label="datos (neto)")
        principal = fits[0]
        ax1.plot(can[z], principal["curva_ideal"][z], color="0.35", lw=1.2, alpha=0.9,
                 label="simulación sin la resolución del detector")
        for f, (c, ls, lw) in zip(fits, self.ESTILOS):
            ax1.plot(can[z], f["curva"][z], color=c, ls=ls, lw=lw,
                     label=f"{f['etiqueta']}: D={f['D']:.3f}±{f['err_D']:.3f}")
        ax1.axvspan(iv[0], iv[1], color="green", alpha=0.06)
        ax1.set_xlabel("Área del pulso [mV·ns]"); ax1.set_ylabel("Cuentas")
        ax1.set_title(f"Forward folding — χ²/ndf ≈ {principal['chi2_red']:.2f}",
                      fontsize=11)
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

        # ---- (der. arriba) EL ESPECTRO SIMULADO, en energia ----
        for f, (c, ls, lw) in zip(fits, self.ESTILOS):
            ax2.plot(f["centros"], f["pesos"], color=c, ls=ls, lw=1.5,
                     label=f["etiqueta"])
        for E in fu.bordes:
            ax2.axvline(E, color="k", ls="--", alpha=0.55)
        if self.log.get():
            ax2.set_yscale("log")
            ax2.set_ylim(max(1.0, min(f["pesos"].max() for f in fits) * 1e-4), None)
        else:
            ax2.set_ylim(0, max(f["pesos"].max() for f in fits) * 1.12)
        ax2.set_xlim(0, 1.15 * max(fu.lineas))
        ax2.set_xlabel("Energía depositada [keV]"); ax2.set_ylabel("eventos / bin")
        ax2.set_title("Espectro SIMULADO (Geant4)", fontsize=10)
        ax2.legend(fontsize=7); ax2.grid(alpha=0.3)

        # ---- (der. abajo) zoom a la zona de los bordes ----
        for f, (c, ls, lw) in zip(fits, self.ESTILOS):
            ax3.plot(f["centros"], f["pesos"], color=c, ls=ls, lw=1.6,
                     label=f["etiqueta"])
        zoom = (0.85 * fu.bordes[0], 1.22 * fu.bordes[-1])
        for E in fu.bordes:
            ax3.axvline(E, color="k", ls="--", alpha=0.55)
            ax3.annotate(f"borde {E:.0f}", (E, 0.97),
                         xycoords=("data", "axes fraction"),
                         rotation=90, fontsize=7, ha="right", va="top")
        m = ((fits[0]["centros"] >= zoom[0]) & (fits[0]["centros"] <= zoom[1]))
        tope = max(f["pesos"][m].max() for f in fits) * 1.15
        ax3.set_xlim(*zoom); ax3.set_ylim(0, tope)
        ax3.set_xlabel("Energía depositada [keV]"); ax3.set_ylabel("eventos / bin")
        ax3.set_title("Zoom: el escalón, y lo que lo rellena", fontsize=10)
        ax3.grid(alpha=0.3)


# ======================================================================
#  Etapa 8 - cadena completa: NCC -> SA -> MCMC
# ======================================================================
class TabAjusteGlobal(Pestana):
    titulo = "8 · NCC→SA→MCMC"
    auto = False
    explicacion = """
ETAPA 8 — la cadena completa, sin números tipeados a mano. Tarda ~5 min la
primera vez (después queda en caché). Apretá "Recalcular".

  1. NCC elige la VENTANA correlacionando contra la simulación, y de yapa da la
     semilla D₀ = B/963.4.
  2. SA busca el óptimo GLOBAL de {A, D, α} dentro de esa ventana, sin la grilla
     de arranque que necesita least_squares. Después se pule localmente.
  3. MCMC mapea el posterior entero: mediana, 68%, y la correlación real.

QUE MIRAR
• Panel izq.: las tres plantillas ajustadas y los residuos normalizados abajo.
  Si el modelo es bueno, los residuos son ruido entre ±2 sin estructura.
• Panel der. arriba: la traza del SA. Termina en el mismo punto que
  least_squares (compará la tabla): eso PRUEBA que el óptimo es global, que era
  lo que no sabíamos.
• Panel der. abajo: el plano D-α. Cada nube es una plantilla. Las elipses rojas
  punteadas son la covarianza linealizada del χ²: si coinciden con las nubes,
  las barras de error de la Etapa 7 eran correctas.

LO IMPORTANTE (y es incómodo)
La distancia ENTRE las nubes es ~10 veces el ancho de cada nube. O sea: el
error que domina NO es el estadístico que mide el MCMC (±0.6%), es elegir qué
plantilla es la correcta (±4%). Un posterior angosto no es una calibración
precisa. Por eso el resultado se reporta con las dos barras separadas.

SESGO DE NEYMAN (§7.31)
El χ² pesa con σ del DATO, así que un bin que fluctúa bajo pesa más y arrastra
el ajuste hacia abajo. Con ~25 cuentas/bin eso es -0.5% en D, MÁS que el error
estadístico. Destildá "Pearson" para ver el sesgo en vivo.

inicio_rel Y EL EXCESO A BAJO CANAL (§7.32)
Los datos tienen un EXCESO sobre la simulación abajo de ~2400 mV·ns (1.4x a
1400, 4.8x a 1100) que muere hacia 2600. El NCC no puede verlo: compara FORMA,
y no tiene forma de saber que ahí hay física que la plantilla no contiene. Con
inicio_rel = 0.45 (el default del paper) la ventana arranca en ~1850 y se come
esa zona; el ajuste la absorbe ensanchando α y bajando D un 2.5%. Por eso el
default acá es 0.62. Bajalo a 0.45 y mirá cómo se va D.
"""

    RAMAS = {"Lvis (keVee)": ("Lvis_keVee", "todo", "tab:blue"),
             "Edep": ("Edep_keV", "todo", "tab:green"),
             "Edep 1 scatt.": ("Edep_keV", "1 scattering", "crimson")}

    def construir_controles(self, frame):
        self.usar = {k: tk.BooleanVar(value=True) for k in self.RAMAS}
        self.ventana = tk.StringVar(value="NCC automática")
        self.lo = tk.DoubleVar(value=300.0)
        self.hi = tk.DoubleVar(value=1900.0)
        self.inicio_rel = tk.DoubleVar(value=0.62)
        self.pasos = tk.IntVar(value=12000)
        self.cadenas = tk.IntVar(value=4)
        self.semilla = tk.IntVar(value=1)
        self.pearson = tk.BooleanVar(value=True)
        self.cache = tk.BooleanVar(value=True)

        for k in self.RAMAS:
            ttk.Checkbutton(frame, text=k, variable=self.usar[k]).pack(side="left",
                                                                      padx=(6, 0))
        ttk.Label(frame, text="ventana").pack(side="left", padx=(10, 2))
        ttk.Combobox(frame, textvariable=self.ventana, width=15, state="readonly",
                     values=["NCC automática", "manual"]).pack(side="left")
        for etiq, var, ancho in [("inicio_rel", self.inicio_rel, 5),
                                 ("lo", self.lo, 6), ("hi", self.hi, 6),
                                 ("pasos", self.pasos, 7),
                                 ("cadenas", self.cadenas, 3),
                                 ("semilla", self.semilla, 3)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=ancho).pack(side="left")
        ttk.Checkbutton(frame, text="Pearson", variable=self.pearson).pack(
            side="left", padx=(10, 0))
        ttk.Checkbutton(frame, text="caché", variable=self.cache).pack(side="left")

    # ------------------------------------------------------------------
    def _clave_cache(self, nom, iv):
        # la huella de los archivos va PRIMERO: si se cambian los datos desde la
        # pestana "Archivos", las cadenas viejas no se reusan (serian de otro
        # experimento).
        return (f"{self.app.datos.huella()}_{nom}_{iv[0]:.0f}_{iv[1]:.0f}_"
                f"{self.pasos.get()}_{self.cadenas.get()}_{self.semilla.get()}_"
                f"{int(self.pearson.get())}")

    def trabajo(self):
        self.exigir_simulacion()
        import os
        from ajuste import mcmc
        from datos_io.spectrum import varianza_predicha
        from ajuste.sim_template import ajustar_sa, ajustar_pearson, modelo_simulado

        fu = self.fuente()
        r = self.app.datos.espectro()
        can, neto, err = r["canales"], r["neto"], r["error"]
        suave = spectrum.suavizar(neto, ventana=21, orden=3)
        carpeta = rutas.CACHE
        os.makedirs(carpeta, exist_ok=True)

        elegidas = [k for k in self.RAMAS if self.usar[k].get()] or ["Lvis (keVee)"]
        salida = []
        for i, nombre in enumerate(elegidas):
            rama, filtro, color = self.RAMAS[nombre]
            self.fase(i / len(elegidas), (i + 1) / len(elegidas))
            self.progreso(f"[{i+1}/{len(elegidas)}] {nombre}: cargando "
                          f"simulación …", 0.02)
            cen, pes = self.app.datos.simulacion(rama, filtro, bins=175)

            # --- 1) la ventana ---
            if self.ventana.get() == "manual":
                iv = (float(self.lo.get()), float(self.hi.get()))
                ncc = None
            else:
                self.progreso(f"[{i+1}/{len(elegidas)}] {nombre}: NCC …", 0.08)
                tope = float(can[-1])
                ncc = buscar_intervalo_simulado(can, suave, cen, pes,
                                                E_gamma=fu.lineas[0],
                                                E_ref=fu.E_ref, E_ref2=fu.E_ref2,
                                                region=(0.1 * tope, tope), paso=1,
                                                n_puntos="auto",
                                                inicio_rel=float(self.inicio_rel.get()))
                iv = (ncc["lo"], ncc["hi"])
            m = (can >= iv[0]) & (can <= iv[1])
            var_fn = varianza_predicha(r["fondo"][m], r["k_m"], r["k_f"])

            # --- 2) el optimo global ---
            self.progreso(f"[{i+1}/{len(elegidas)}] {nombre}: recocido "
                          f"simulado …", 0.20)
            sa = ajustar_sa(can, neto, cen, pes, iv, errores=err,
                            D0=(ncc["D_seed"] if ncc else None),
                            alpha0=(ncc["alpha_seed"] if ncc else 6.0),
                            E_ref=fu.E_ref)
            p0 = [sa["A_sa"], sa["D_sa"], sa["alpha_sa"]]
            if self.pearson.get():
                fit = ajustar_pearson(can, neto, cen, pes, iv, var_fn,
                                      errores=err, p0=p0, E_ref=fu.E_ref)
            else:
                fit = sa["fit"]

            # --- 3) el posterior ---
            clave = self._clave_cache(nombre, iv)
            ruta = os.path.join(carpeta, f"mcmc_{clave}.npz")
            muestras = None
            if self.cache.get() and os.path.exists(ruta):
                try:
                    z = np.load(ruta)
                    muestras, R_hat, n_eff, acep = (z["muestras"], z["R_hat"],
                                                    z["n_eff"], z["aceptacion"])
                    self.progreso(f"[{i+1}/{len(elegidas)}] {nombre}: MCMC "
                                  f"en caché", 0.95)
                except Exception:
                    muestras = None
            if muestras is None:
                def cb(k, nk, paso, npas, nombre=nombre, i=i):
                    # el MCMC es el 60% final del trabajo de cada plantilla
                    avance = 0.40 + 0.60 * (k + paso / npas) / nk
                    self.progreso(f"[{i+1}/{len(elegidas)}] {nombre}: MCMC cadena "
                                  f"{k+1}/{nk}, paso {paso}/{npas}", avance)
                post = mcmc.posterior_calibracion(
                    r, cen, pes, iv, fit, n_cadenas=int(self.cadenas.get()),
                    n_pasos=int(self.pasos.get()), semilla=int(self.semilla.get()),
                    callback=cb, usar_pearson=bool(self.pearson.get()))
                muestras = post["muestras"]
                R_hat, n_eff, acep = post["R_hat"], post["n_eff"], post["aceptacion"]
                if self.cache.get():
                    try:
                        np.savez_compressed(ruta, muestras=muestras, R_hat=R_hat,
                                            n_eff=n_eff, aceptacion=acep)
                    except Exception:
                        pass

            D = muestras[:, 1]
            curva = modelo_simulado(can, cen, pes, fit["A"], fit["D"], fit["alpha"],
                                    fit["beta"])
            salida.append({
                "nombre": nombre, "color": color, "iv": iv, "ncc": ncc,
                "sa": sa, "fit": fit, "curva": curva, "muestras": muestras,
                "R_hat": R_hat, "n_eff": n_eff, "aceptacion": acep,
                "D_med": float(np.median(D)),
                "D_p16": float(np.percentile(D, 16)),
                "D_p84": float(np.percentile(D, 84)),
                "al_med": float(np.median(muestras[:, 2])),
                "corr": float(np.corrcoef(muestras[:, 1], muestras[:, 2])[0, 1]),
            })
        return r, salida, fu

    # ------------------------------------------------------------------
    def pintar(self, res):
        from ajuste import mcmc
        r, fits, fu = res
        can, neto, err = r["canales"], r["neto"], r["error"]
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[2.4, 1])
        ax1 = self.fig.add_subplot(gs[0, 0])
        axr = self.fig.add_subplot(gs[1, 0], sharex=ax1)
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax3 = self.fig.add_subplot(gs[1, 1])

        # ---- (izq) datos, modelos y residuos ----
        iv0 = fits[0]["iv"]
        z = (can >= 0.4 * iv0[0]) & (can <= min(1.1 * iv0[1], can[-1]))
        ax1.step(can[z], neto[z], where="mid", color="0.75", lw=1, label="datos (neto)")
        for f in fits:
            ax1.plot(can[z], f["curva"][z], color=f["color"], lw=2,
                     label=f"{f['nombre']}: D = {f['D_med']:.3f}"
                           f" ± {0.5*(f['D_p84']-f['D_p16']):.3f}")
            m = (can >= f["iv"][0]) & (can <= f["iv"][1])
            axr.plot(can[m], (f["curva"][m] - neto[m]) / err[m], ".", ms=2,
                     color=f["color"], alpha=0.7)
        # residuos FUERA de la ventana (extrapolando el modelo): ahi se ve el
        # exceso de datos a bajo canal que la plantilla no contiene (§7.32)
        fuera = z & (can < min(f["iv"][0] for f in fits))
        if fuera.any():
            axr.plot(can[fuera], (fits[0]["curva"][fuera] - neto[fuera]) / err[fuera],
                     ".", ms=2, color="0.55", alpha=0.6)
        fa = fits[0]
        ax1.axvspan(fa["iv"][0], fa["iv"][1], color="green", alpha=0.06,
                    label=f"ventana NCC [{fa['iv'][0]:.0f}, {fa['iv'][1]:.0f}]")
        ax1.set_ylabel("Cuentas"); ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
        ax1.set_title(f"Ventana por NCC, óptimo por SA, error por MCMC — "
                      f"χ²/ndf = {fa['fit']['chi2_red']:.3f}", fontsize=10)
        axr.axhline(0, color="k", lw=0.8)
        for s in (-2, 2):
            axr.axhline(s, color="0.6", lw=0.7, ls=":")
        axr.set_ylim(-8, 4.5); axr.set_ylabel("(modelo − dato)/σ")
        axr.axvspan(fa["iv"][0], fa["iv"][1], color="green", alpha=0.06)
        axr.set_xlabel("Área del pulso [mV·ns]   (gris = fuera de la ventana)")
        axr.grid(alpha=0.3)

        # ---- (der. arriba) la traza del SA + la tabla ----
        h = fa["sa"]["hist"]
        ax2.plot(h["actual"], color="0.65", lw=0.8, label="costo actual")
        ax2.plot(h["mejor"], "tab:blue", lw=1.8, label="mejor hasta ahora")
        ax2.set_yscale("log"); ax2.set_xlabel("ciclo de temperatura")
        ax2.set_ylabel("χ²")
        ax2.set_title(f"SA sobre {fa['nombre']} ({fa['sa']['n_evals']} evals)",
                      fontsize=9)
        ax2.grid(alpha=0.3)
        # aire arriba para que la tabla no tape la curva del "mejor hasta ahora"
        y0, y1 = ax2.get_ylim()
        ax2.set_ylim(y0, y1 * 12)
        ax2.legend(fontsize=7, loc="lower right")
        filas = [f"{'':7s}{'SA':>9s}{'pulido':>9s}{'MCMC':>9s}",
                 f"{'D':7s}{fa['sa']['D_sa']:>9.4f}{fa['fit']['D']:>9.4f}"
                 f"{fa['D_med']:>9.4f}",
                 f"{'α':7s}{fa['sa']['alpha_sa']:>9.3f}"
                 f"{fa['fit']['alpha']:>9.3f}{fa['al_med']:>9.3f}"]
        ax2.text(0.02, 0.97, "\n".join(filas), transform=ax2.transAxes,
                 fontsize=7.5, family="monospace", va="top",
                 bbox=dict(fc="white", ec="0.7", alpha=0.9))

        # ---- (der. abajo) el plano D-alpha ----
        for f in fits:
            D, AL = f["muestras"][:, 1], f["muestras"][:, 2]
            H, cx, cy, niveles = mcmc.contornos_2d(D, AL, bins=50)
            if niveles:
                ax3.contour(cx, cy, H, levels=niveles, colors=f["color"],
                            linewidths=[0.9, 1.9][:len(niveles)])
            ax3.plot(f["D_med"], f["al_med"], "o", color=f["color"], ms=5,
                     label=f"{f['nombre']}  ρ(D,α)={f['corr']:+.2f}")
            cov = f["fit"].get("cov")
            if cov is not None and np.all(np.isfinite(cov)):
                ex, ey = mcmc.elipse_covarianza([f["fit"]["D"], f["fit"]["alpha"]],
                                                cov[1:3, 1:3], n_sigma=1.0)
                ax3.plot(ex, ey, ls="--", lw=1.2, color="crimson", alpha=0.8)
            ax3.plot(f["sa"]["D_sa"], f["sa"]["alpha_sa"], "x", color="k", ms=6)
        ax3.set_xlabel("D [mV·ns/keV]"); ax3.set_ylabel("α [√keV]")
        ax3.set_title("Posterior 68/95 vs covarianza linealizada (rojo)", fontsize=9)
        ax3.legend(fontsize=6.5, loc="upper right"); ax3.grid(alpha=0.3)

        # ---- el resumen numerico, al panel de texto ----
        L = ["RESULTADO DE LA CADENA\n"]
        for f in fits:
            em, ep = f["D_med"] - f["D_p16"], f["D_p84"] - f["D_med"]
            L.append(f"{f['nombre']}")
            L.append(f"   ventana  = [{f['iv'][0]:.0f}, {f['iv'][1]:.0f}] mV·ns"
                     + (f"   (ρ = {f['ncc']['rho']:.4f})" if f["ncc"] else "   (manual)"))
            L.append(f"   D        = {f['D_med']:.4f}  -{em:.4f} +{ep:.4f} mV·ns/keV")
            L.append(f"   α        = {f['al_med']:.3f} √keV   "
                     f"→ {100*2.355*f['al_med']/np.sqrt(fu.E_ref):.1f}% FWHM a "
                     f"{fu.E_ref:.0f} keV")
            L.append(f"   χ²/ndf   = {f['fit']['chi2_red']:.3f}"
                     f"   corr(D,α) = {f['corr']:+.3f}")
            L.append(f"   R̂        = " + " ".join(f"{v:.4f}" for v in f["R_hat"])
                     + f"   n_eff(D) = {f['n_eff'][1]:.0f}"
                     f"   acept = {np.mean(f['aceptacion']):.2f}")
            L.append("")
        if len(fits) > 1:
            Ds = np.array([f["D_med"] for f in fits])
            est = np.mean([0.5 * (f["D_p84"] - f["D_p16"]) for f in fits])
            L.append(f"ESTADÍSTICO (MCMC)  : ±{est:.4f}  ({100*est/Ds.mean():.1f}%)")
            L.append(f"SISTEMÁTICO (plantilla): ±{0.5*(Ds.max()-Ds.min()):.4f}  "
                     f"({100*0.5*(Ds.max()-Ds.min())/Ds.mean():.1f}%)")
            L.append(f"→ el sistemático es {0.5*(Ds.max()-Ds.min())/est:.0f}x el "
                     f"estadístico.")
        L.append("\n" + self.explicacion.strip())
        self.escribir("\n".join(L))


# ======================================================================
#  Etapa 9a - el observable: como se mide el area
# ======================================================================
class TabArea(Pestana):
    titulo = "9 · El área"
    auto = False
    explicacion = """
ETAPA 9a — el error más caro del proyecto estaba en el OBSERVABLE.

Hasta la Etapa 8 el área se calculaba así: línea base = promedio de las primeras
204 muestras (40 ns) del registro, y la integral sobre los 200 ns completos.
Pero el pulso ARRANCA A LOS 25 ns y tiene el pico a los 48. O sea la ventana de
línea base se comía la subida del propio pulso.

QUE MIRAR
• Panel 1 (forma media): con la línea base corta el pre-pulso queda en cero;
  con la de 40 ns queda LEVANTADO. Ese escalón, multiplicado por los 200 ns de
  integración, es lo que se le restaba a cada evento.
• Panel 2 (integral/amplitud): con línea base limpia la distribución es angosta
  (std 0.84 ns). Con la vieja, std 11.3 y un 26% de los eventos con el signo
  INVERTIDO, que el abs() reciclaba como buenos.
• Panel 3 (posición del pulso): el jitter del trigger mueve el pico entre 18 y
  54 ns, así que cada evento perdía una fracción DISTINTA de su carga. Eso es un
  ensanchamiento multiplicativo: era el 43% de FWHM que llamábamos resolución.
• Panel 4 (espectro): el antes y el después. El eje se comprime 3.4x y aparece
  un borde limpio donde antes había una caída de 2000 mV·ns de ancho.

EL ARREGLO (spectrum._area_pulso_gate)
Onset por evento (cruce del 10% retrocediendo desde el pico), línea base =
mediana hasta 4 ns antes del onset, integral en [onset−6, onset+55] ns.
Verificado: el pulso medio devuelve el 100% de su área en 80 ns y área/amplitud
es constante al 1% entre 20 y 220 mV de amplitud.
"""

    def construir_controles(self, frame):
        self.n_ev = tk.IntVar(value=8000)
        ttk.Label(frame, text="eventos a analizar").pack(side="left", padx=(8, 2))
        ttk.Entry(frame, textvariable=self.n_ev, width=8).pack(side="left")
        ttk.Label(frame, text="(~1 s cada 2000)").pack(side="left", padx=8)

    def trabajo(self):
        N = max(int(self.n_ev.get()), 500)
        self.fase(0.0, 0.85)
        self.progreso("leyendo pulsos crudos ...", 0.0)
        P_corta = np.zeros(1024)
        P_larga = np.zeros(1024)
        n_forma = 0
        dt = 0.198
        s_corta, s_larga, t_pico = [], [], []
        for k, (ts, pares, sc) in enumerate(
                spectrum._iter_eventos(self.app.datos.archivo_muestra)):
            if k >= N:
                break
            if k % 250 == 0:
                self.progreso("pulsos crudos: %d de %d" % (k, N), k / N)
            # `pares` es una lista de cadenas "t,v" en el XML y un array
            # (n, 2) en el binario; `pares_a_pulso` unifica los dos.
            t, v = spectrum.pares_a_pulso(pares)
            if t is None or t.size != 1024:
                continue
            # TIEMPO RELATIVO al registro: el DRSOSC escribe t respecto de la
            # celda de trigger, asi que t[0] se mueve entre -78 y +121 ns
            # (spectrum._area_pulso_gate, punto 3). Con t absoluto las ventanas
            # de abajo agarran cualquier cosa: en las corridas del 26-ago-2026 el
            # pico cae fuera de (30, 60) ns y no entraba NINGUN evento.
            tau = t - t[0]
            dt = float(np.median(np.diff(tau)))
            n_base = max(int(round(12.0 / dt)), 5)      # 0-12 ns, antes del pulso
            n_larga = int(len(v) * 0.2)                 # el 20% de la vieja
            b_corta = np.median(v[:n_base])
            b_larga = v[:n_larga].mean()
            vc = v - b_corta
            i = int(np.argmin(vc))
            a = -vc.min()
            t_pico.append(tau[i])
            if a < 20 or not (8 < tau[i] < tau[-1] - 60):
                continue
            s_corta.append(np.trapezoid(vc, x=tau) / a)
            s_larga.append(np.trapezoid(v - b_larga, x=tau) / a)
            P_corta += vc / a
            P_larga += (v - b_larga) / a
            n_forma += 1
        self.fase(0.0, 1.0)
        self.progreso("histogramas ...", 0.88)
        esp_new = self.app.datos.espectro(2500.0, 500, 0.05, "gate")
        esp_old = self.app.datos.espectro(8000.0, 800, 0.05, "legacy")
        return {"P_corta": P_corta / max(n_forma, 1),
                "P_larga": P_larga / max(n_forma, 1),
                "t": np.arange(1024) * dt, "t_base": 1024 * dt * 0.2,
                "s_corta": np.array(s_corta), "s_larga": np.array(s_larga),
                "t_pico": np.array(t_pico), "n": n_forma,
                "esp_new": esp_new, "esp_old": esp_old}

    def pintar(self, r):
        ax = self.fig.subplots(2, 2).ravel()

        a = ax[0]
        a.plot(r["t"], r["P_corta"], color="C0", lw=1.4,
               label="línea base 0–12 ns (limpia)")
        a.plot(r["t"], r["P_larga"], color="C3", lw=1.4,
               label="línea base del 20%% del registro (la vieja, %.0f ns)"
                     % r["t_base"])
        a.axvspan(0, r["t_base"], color="crimson", alpha=0.08)
        a.axhline(0, color="k", lw=0.7)
        a.set_xlim(0, 120)
        a.set_xlabel("t [ns]")
        a.set_ylabel("v / amplitud")
        a.set_title("Forma media normalizada (%d pulsos)" % r["n"])
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

        a = ax[1]
        b = np.linspace(-30, 30, 120)
        a.hist(r["s_corta"], bins=b, histtype="step", color="C0",
               label="limpia: std %.2f" % r["s_corta"].std())
        a.hist(r["s_larga"], bins=b, histtype="step", color="C3",
               label="vieja: std %.2f" % r["s_larga"].std())
        a.axvline(0, color="k", lw=0.8)
        pos = float(np.mean(r["s_larga"] > 0))
        a.set_xlabel("integral / amplitud [ns]   (el pulso es negativo)")
        a.set_ylabel("eventos")
        a.set_title("Área por unidad de amplitud — %.0f%% con el signo dado vuelta"
                    % (100 * pos))
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

        a = ax[2]
        a.hist(r["t_pico"], bins=100, range=(0, 1024 * (r["t"][1] - r["t"][0])),
               color="0.5")
        a.axvspan(0, r["t_base"], color="crimson", alpha=0.15,
                  label="ventana de línea base vieja")
        a.set_xlabel("t del pico, relativo al registro [ns]")
        a.set_ylabel("eventos")
        a.set_title("Jitter del trigger: la pérdida era distinta en cada evento")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

        a = ax[3]
        n, o = r["esp_new"], r["esp_old"]
        a.step(n["canales"], n["neto"], where="mid", color="C0", lw=1.2,
               label="área con compuerta")
        a.step(o["canales"], o["neto"], where="mid", color="C3", lw=1, alpha=0.8,
               label="área vieja")
        a.set_xlim(0, 8000)
        a.set_xlabel("Área del pulso [mV·ns]")
        a.set_ylabel("Cuentas netas")
        a.set_title("El espectro, antes y después")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)


# ======================================================================
#  Etapa 9b - modelo extendido: los sistematicos como parametros
# ======================================================================
# ----------------------------------------------------------------------
#  el ajuste de la Etapa 9, como funcion
# ----------------------------------------------------------------------
def ajuste_extendido(datos, r, iv, rama="Lvis_keVee", w_libre=True, w_fijo=None,
                     exceso=True, forma="exp", progreso=None):
    """
    Corre `modelo_ext` sobre los datos y la FUENTE que estan cargados.

    Estaba adentro de la pestana 10. Salio afuera porque la pestana de
    Calibracion necesita correr el mismo ajuste media docena de veces
    (fijando w, cambiando la forma del exceso, moviendo la ventana) para armar
    el presupuesto de error: ahi es donde se miden los sistematicos, y no tiene
    sentido tener dos copias de la receta que puedan divergir.

    Todo lo que depende de la fuente sale de `datos.fuente`: E_ref para el mapa
    por el origen y para reportar la resolucion, y el techo del histograma de
    la simulacion. Todo lo que depende del DETECTOR sale del borde estimado por
    NCC (`datos.borde_estimado`), asi que no hay ningun canal escrito a mano.

    forma : "exp" o "pot", la forma que se le supone al exceso de bajo canal.
            Es una variable GLOBAL de modelo_ext (la usa `forma_exceso`), asi
            que se cambia y se restaura aca. Ojo: el modelo que queda en el
            resultado hay que evaluarlo con la misma forma con que se ajusto.
    """
    from ajuste import modelo_ext as me
    datos.exigir_sim_compatible()
    fu = datos.fuente
    can, neto, err = r["canales"], r["neto"], r["error"]
    B = float(datos.borde_estimado()["edge_ch"])
    # el borde por NCC pone en escala las ventanas de arranque: con el Co-60 de
    # la Etapa 9 (B = 1229) da (922, 1905) y (983, 1082), o sea los 900-1900 y
    # 1000-1100 que estaban tipeados a mano.
    iv_limpia = (0.75 * B, min(1.55 * B, float(can[-1])))

    def var_de(intervalo):
        m = (can >= intervalo[0]) & (can <= intervalo[1])
        return spectrum.varianza_predicha(r["fondo"][m], r["k_m"], r["k_f"])

    if progreso:
        progreso("cargando la simulacion ...")
    centros, h1, hm, w_nat = datos.componentes(rama, bins=175)

    forma_previa = me.FORMA_EXCESO
    me.FORMA_EXCESO = forma
    try:
        p0 = np.array([1.0, B / fu.E_ref, 3.5, w_nat, 0.0, 0.285 * B, 0.0, 0.0])
        z = (can > 0.80 * B) & (can < 0.88 * B)
        p0[me.I_A] = neto[z].mean() / me.modelo_ext(can, centros, h1, hm,
                                                    p0)[z].mean()
        if progreso:
            progreso("semilla en la zona limpia ...")
        f0 = me.ajustar(can, neto, centros, h1, hm, iv_limpia, p0,
                        ["A", "D", "alpha"], var_fn=var_de(iv_limpia),
                        errores=err, E_ref=fu.E_ref)
        libres = ["A", "D", "alpha"]
        p = f0["p"].copy()
        if w_libre:
            libres.append("w")
        elif w_fijo is not None:
            p[me.I_W] = float(w_fijo)
        if exceso:
            p[me.I_C] = 80.0
            p[me.I_LAM] = 150.0 if forma == "pot" else 0.285 * B
            libres += ["C", "lambda"]
        if progreso:
            progreso("ajuste completo ...")
        f = me.ajustar(can, neto, centros, h1, hm, iv, p, libres,
                       var_fn=var_de(iv), errores=err, E_ref=fu.E_ref)
    finally:
        me.FORMA_EXCESO = forma_previa

    f["semilla"] = f0
    f["w_nat"] = w_nat
    f["forma"] = forma
    f["comp"] = (centros, h1, hm, w_nat)
    f["var_de"] = var_de
    f["rama"] = rama
    return f


class TabModeloExt(Pestana):
    titulo = "10 · Modelo ext."
    auto = False
    explicacion = """
ETAPA 9b — meter los sistemáticos DENTRO del modelo.

  S(ch) = A·[w·h_1sc + (1−w)·h_multi] ⊗ G(σ_E) + C·exp(−ch/λ)

La Etapa 8 reportaba ±2.7% de sistemático comparando tres plantillas, y cortaba
la ventana por arriba del exceso a bajo canal. Acá las dos cosas son parámetros:
w es cuánto scattering SIMPLE piden los datos, y C·exp(−ch/λ) es el exceso.

QUE MIRAR
• Al liberar w y modelar el exceso, la ventana puede bajar de 900 a 300 mV·ns
  (más estadística) sin que D se mueva.
• w sale 0.796 ± 0.070 contra el 0.751 de la simulación: COINCIDEN (0.6σ). La
  geometría simulada reproduce bien cuánto scattering múltiple hay. Fijar w = 1
  (lo que asume la plantilla analítica del paper) empeora el ajuste en
  Δχ² = +13 y sube D un 3.4%.
• El posterior D–w tiene corr = +0.96. Eso es lo importante: la ignorancia sobre
  la plantilla ya está DENTRO del error de D, no colgada aparte.
• Los residuos son planos: el exceso no está tapando estructura del borde.

CONTROL
Destildá "w libre" y fijá w en 1 (solo 1 scattering) o en 0.751 (simulación
completa): eso es lo que el análisis viejo llamaba sistemático de plantilla.
"""

    def construir_controles(self, frame):
        self.rama = tk.StringVar(value="Lvis_keVee")
        self.lo = tk.DoubleVar(value=300.0)
        self.hi = tk.DoubleVar(value=1900.0)
        self.w_libre = tk.BooleanVar(value=True)
        self.w_fijo = tk.DoubleVar(value=0.751)
        self.exceso = tk.BooleanVar(value=True)
        self.mcmc = tk.BooleanVar(value=False)
        self.ventana = tk.StringVar(value="auto (NCC)")
        ttk.Label(frame, text="ventana").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.ventana, width=11, state="readonly",
                     values=["auto (NCC)", "manual"]).pack(side="left")
        ttk.Label(frame, text="rama").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.rama, width=11, state="readonly",
                     values=["Lvis_keVee", "Edep_keV"]).pack(side="left")
        for etiq, var in [("lo", self.lo), ("hi", self.hi)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=7).pack(side="left")
        ttk.Checkbutton(frame, text="w libre", variable=self.w_libre).pack(
            side="left", padx=8)
        ttk.Entry(frame, textvariable=self.w_fijo, width=6).pack(side="left")
        ttk.Checkbutton(frame, text="exceso", variable=self.exceso).pack(
            side="left", padx=8)
        ttk.Checkbutton(frame, text="MCMC", variable=self.mcmc).pack(
            side="left", padx=8)

    def trabajo(self):
        self.exigir_simulacion()
        from ajuste import modelo_ext as me
        fu = self.fuente()
        d = self.app.datos
        r = d.espectro(2500.0, 500, 0.05, "gate")
        can, neto, err = r["canales"], r["neto"], r["error"]
        iv = d.ventana_ajuste(
            manual=(float(self.lo.get()), float(self.hi.get()))
            if self.ventana.get() == "manual" else None)

        f = ajuste_extendido(d, r, iv, rama=self.rama.get(),
                             w_libre=bool(self.w_libre.get()),
                             w_fijo=float(self.w_fijo.get()),
                             exceso=bool(self.exceso.get()),
                             progreso=self.progreso)
        centros, h1, hm, w_nat = f["comp"]
        libres, f0 = f["libres"], f["semilla"]

        post = None
        if self.mcmc.get():
            from ajuste.mcmc import correr_cadenas, resumen
            self.fase(0.30, 1.0)
            self.progreso("MCMC ...", 0.0)
            lp, _, _ = me.log_posterior(can, neto, centros, h1, hm, iv, f["p"],
                                        libres, var_fn=f["var_de"](iv),
                                        errores=err)
            cov = np.asarray(f["cov"], float)
            if cov.ndim != 2 or not np.all(np.isfinite(cov)):
                cov = np.diag((0.02 * np.abs(f["q"]) + 1e-9) ** 2)

            def cb(k, nk, paso, npas):
                self.progreso("MCMC cadena %d/%d, paso %d/%d"
                              % (k + 1, nk, paso, npas),
                              (k + paso / npas) / nk)
            post = resumen(correr_cadenas(lp, f["q"], cov, n_cadenas=4,
                                          n_pasos=6000, dispersion=1.0,
                                          callback=cb),
                           nombres=libres)
        return {"r": r, "f": f, "f0": f0, "post": post, "libres": libres,
                "centros": centros, "h1": h1, "hm": hm, "w_nat": w_nat,
                "iv": iv, "fuente": fu}

    def pintar(self, d):
        from ajuste import modelo_ext as me
        r, f = d["r"], d["f"]
        can, neto, err = r["canales"], r["neto"], r["error"]
        ax = self.fig.subplots(1, 2)

        z = (can >= 0.25 * d["iv"][0]) & (can <= min(1.15 * d["iv"][1], can[-1]))
        a = ax[0]
        a.errorbar(can[z], neto[z], yerr=err[z], fmt=".", ms=2, lw=0.5,
                   color="0.65", label="datos (neto)")
        a.plot(can[z], me.modelo_ext(can[z], d["centros"], d["h1"], d["hm"],
                                     f["p"]), color="C0", lw=1.8,
               label="modelo: D = %.3f" % f["D"])
        p = f["p"].copy()
        if p[me.I_C] > 0:
            a.plot(can[z], p[me.I_C] * np.exp(-can[z] / p[me.I_LAM]), "k--",
                   lw=1, label="exceso empírico")
            q = p.copy()
            q[me.I_C] = 0.0
            a.plot(can[z], me.modelo_ext(can[z], d["centros"], d["h1"], d["hm"],
                                         q), color="C3", ls=":", lw=1.2,
                   label="solo Compton")
        a.axvspan(d["iv"][0], d["iv"][1], color="green", alpha=0.06)
        a.set_xlabel("Área del pulso [mV·ns]")
        a.set_ylabel("Cuentas netas")
        a.set_title("chi2/ndf = %.3f (%d gdl)" % (f["chi2_red"], f["ndf"]))
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

        a = ax[1]
        post = d["post"]
        if post is not None and "w" in d["libres"]:
            from ajuste.mcmc import contornos_2d
            iD, iw = d["libres"].index("D"), d["libres"].index("w")
            x, y = post["muestras"][:, iD], post["muestras"][:, iw]
            H, cx, cy, niv = contornos_2d(x, y)
            a.plot(x[::20], y[::20], ".", ms=1, alpha=0.25, color="0.6")
            if niv:
                a.contour(cx, cy, H, levels=niv, colors="C0")
            a.axhline(d["w_nat"], color="C3", ls="--", lw=1,
                      label="w de la simulación = %.3f" % d["w_nat"])
            a.set_xlabel("D [mV·ns/keV]")
            a.set_ylabel("w (fracción de 1 scattering)")
            a.set_title("Posterior   corr(D,w) = %+.2f" % post["corr"][iD, iw])
            a.legend(fontsize=8)
        else:
            m = (can >= d["iv"][0]) & (can <= d["iv"][1])
            mod = me.modelo_ext(can[m], d["centros"], d["h1"], d["hm"], f["p"])
            a.plot(can[m], (neto[m] - mod) / np.maximum(err[m], 1), ".", ms=3)
            a.axhline(0, color="k", lw=1)
            a.axhspan(-1, 1, color="green", alpha=0.12)
            a.set_ylim(-6, 6)
            a.set_xlabel("Área del pulso [mV·ns]")
            a.set_ylabel("(dato − modelo)/σ")
            a.set_title("Residuos (tildá MCMC para ver el posterior D–w)")
        a.grid(alpha=0.3)

        L = ["D      = %.4f ± %.4f mV·ns/keV   (%.0f keV por 1000 mV·ns)"
             % (f["D"], f["err"].get("D", float("nan")), 1000 / f["D"]),
             "alpha  = %.3f  → FWHM(%.0f) = %.1f%%"
             % (f["alpha"], f["E_ref"], 100 * f["fwhm_rel_ref"]),
             "w      = %.3f%s   (simulación: %.3f)"
             % (f["w"], (" ± %.3f" % f["err"]["w"]) if "w" in f["err"]
                else "  (FIJO)", d["w_nat"]),
             "chi2/ndf = %.3f   (%d gdl)" % (f["chi2_red"], f["ndf"]),
             "semilla 3 par en zona limpia: D = %.4f" % d["f0"]["D"]]
        if post is not None:
            iD = d["libres"].index("D")
            L.append("posterior: D = %.4f −%.4f +%.4f   R_hat = %.4f"
                     % (post["mediana"][iD], post["err_menos"][iD],
                        post["err_mas"][iD], post["R_hat"][iD]))
        L.append("\n" + self.explicacion.strip())
        self.escribir("\n".join(L))


# ======================================================================
#  Resumen / calibracion
# ======================================================================
RUTA_REGISTRO = os.path.join(rutas.CACHE, "calibraciones.json")


def cargar_registro():
    """
    Las calibraciones ya medidas, una por (fuente, archivos).

    Por que un registro y no una constante: una recta de calibracion se hace
    con VARIAS fuentes, y cada una se mide en una corrida distinta, otro dia.
    Sin un lugar donde dejar el resultado de cada una, comparar Cs-137 con
    Co-60 obliga a copiar numeros a mano de una sesion a otra, que es
    exactamente como se cuelan los numeros viejos (CLAUDE.md 7.40).
    """
    try:
        with open(RUTA_REGISTRO, encoding="utf-8") as f:
            reg = json.load(f)
        return reg if isinstance(reg, dict) else {}
    except Exception:
        return {}


def guardar_en_registro(clave, entrada):
    """Agrega o pisa una entrada y devuelve el registro completo."""
    reg = cargar_registro()
    reg[clave] = entrada
    try:
        os.makedirs(os.path.dirname(RUTA_REGISTRO), exist_ok=True)
        with open(RUTA_REGISTRO, "w", encoding="utf-8") as f:
            json.dump(reg, f, indent=1, ensure_ascii=False)
    except Exception:
        pass
    return reg


def recta_por_origen(puntos):
    """
    Ajuste ponderado de area = D*E por el ORIGEN sobre varias fuentes.

    `puntos` = [(E_borde, area_borde, sigma_area), ...] -> (D, sigma_D).

    Ojo con que significa: cada fuente ya entrego SU D con el mapa por el
    origen, asi que su punto cae exacto sobre su propia recta. Lo que este
    ajuste hace es promediarlos, y lo que hay que mirar es la DISPERSION entre
    fuentes: eso, y no el ajuste, es el test de linealidad.
    """
    E = np.array([p[0] for p in puntos], float)
    A = np.array([p[1] for p in puntos], float)
    s = np.maximum(np.array([p[2] for p in puntos], float), 1e-12)
    peso = 1.0 / s ** 2
    den = float(np.sum(peso * E * E))
    if den <= 0:
        return float("nan"), float("nan")
    return float(np.sum(peso * E * A) / den), float(1.0 / np.sqrt(den))


class TabCalibracion(Pestana):
    titulo = "★ Calibración"
    auto = False

    # Sistematico de la DEFINICION DEL AREA, medido en la Etapa 9 barriendo
    # compuerta (45/55/70/90 ns), onset (5/10/20%) y descarte por calidad:
    # D entre 1.2536 y 1.2799 sobre 1.261, o sea +-1.04%. No se re-mide aca
    # porque cada variante obliga a re-parsear los XML (~40 s cada una); se
    # arrastra como fraccion y se dice de donde sale.
    ERR_AREA_REL = 0.0104

    explicacion = """
EL RESULTADO, calculado sobre los archivos y la FUENTE que están cargados.

Ya no es un cartel con los números del Co-60: corre la cadena entera sobre lo
que elegiste en «Archivos» y arma el presupuesto de error midiéndolo, no
citándolo. Lo que cambia con la fuente:

• el borde de referencia E_max (963.4 keV en Co-60, 477.3 en Cs-137) y por lo
  tanto el mapa por el origen D = B/E_max;
• el modelo analítico: dos bordes solapados si hay dos líneas, uno solo si hay
  una;
• la ventana de ajuste, que sale del NCC y se adapta al detector;
• la plantilla Geant4, que TIENE que ser de esta fuente. Si el .root cargado
  simula otros gammas, la pestaña lo dice y usa la ruta analítica en vez de
  devolver un número equivocado con buena cara.

LAS DOS RUTAS
• Klein-Nishina (siempre disponible): scattering SIMPLE. En el Co-60 de este
  proyecto da +4.3% de más que el modelo con simulación, porque el scattering
  múltiple mueve cuentas por encima del borde y el ajuste lo compensa
  corriendo B hacia arriba (§7.18). Ese sesgo se midió acá; en otro detector
  no tiene por qué valer lo mismo.
• Modelo extendido (Etapa 9, necesita el .root): plantilla Geant4 con la mezcla
  de scattering w ajustable + el exceso de bajo canal modelado. Es el número
  bueno cuando está disponible.

DOS ESCALAS, no mezclarlas: el SiPM ve LUZ, no energía. La plantilla Lvis da
keVee (keV electrón equivalente, la que sirve para la rama de neutrones) y la
Edep da keV depositados. Su cociente reproduce el Birks de la simulación, pero
eso NO es un control independiente: Lvis es Edep por el quenching (§7).

QUE MIRAR EN EL GRÁFICO
• Izq.: todas las variantes sobre los MISMOS datos. La banda es el sistemático
  total. Si una variante se va lejos, ahí está el problema, no en la
  estadística.
• Der.: la recta de calibración con los bordes de esta fuente, y los puntos de
  las OTRAS fuentes ya medidas (se guardan en .cache/calibraciones.json). Con
  dos fuentes o más, que los D coincidan es el test de linealidad; con una
  sola, la linealidad se está asumiendo.

EL PRESUPUESTO DE ERROR
  estadística        — del ajuste, ya incluye w y el exceso por correlación
  ventana            — barrido del límite inferior
  forma del exceso   — exponencial contra ley de potencia: describen los datos
                       igual de bien y dan D distinto. Suele ser el que manda.
  definición del área— ±1.04%, medido en la Etapa 9 (§7.35). Es el único que se
                       arrastra en vez de medirse acá: cada variante obliga a
                       re-parsear los XML.
"""

    def construir_controles(self, frame):
        self.rama = tk.StringVar(value="Lvis_keVee")
        self.ventana = tk.StringVar(value="auto (NCC)")
        self.lo = tk.DoubleVar(value=300.0)
        self.hi = tk.DoubleVar(value=1900.0)
        self.sist = tk.BooleanVar(value=True)
        ttk.Label(frame, text="escala principal").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.rama, width=11, state="readonly",
                     values=["Lvis_keVee", "Edep_keV"]).pack(side="left")
        ttk.Label(frame, text="ventana").pack(side="left", padx=(8, 2))
        ttk.Combobox(frame, textvariable=self.ventana, width=11, state="readonly",
                     values=["auto (NCC)", "manual"]).pack(side="left")
        for etiq, var in [("lo", self.lo), ("hi", self.hi)]:
            ttk.Label(frame, text=etiq).pack(side="left", padx=(8, 2))
            ttk.Entry(frame, textvariable=var, width=7).pack(side="left")
        ttk.Checkbutton(frame, text="barrer sistemáticos", variable=self.sist
                        ).pack(side="left", padx=10)
        ttk.Button(frame, text="Exportar…", command=self.exportar).pack(
            side="left", padx=4)
        ttk.Button(frame, text="Vaciar registro", command=self.vaciar).pack(
            side="left", padx=2)

    # ------------------------------------------------------------------
    def exportar(self):
        """Deja el registro completo en un .json fuera de .cache."""
        ruta = filedialog.asksaveasfilename(
            parent=self, title="Exportar las calibraciones medidas",
            defaultextension=".json", initialfile="calibraciones.json",
            filetypes=[("JSON", "*.json"), ("Todos", "*.*")])
        if not ruta:
            return
        try:
            with open(ruta, "w", encoding="utf-8") as f:
                json.dump(cargar_registro(), f, indent=1, ensure_ascii=False)
            messagebox.showinfo("Exportar", "Guardado en:\n" + ruta)
        except Exception as e:
            messagebox.showerror("Exportar", str(e))

    def vaciar(self):
        if not messagebox.askyesno(
                "Registro", "¿Borrar TODAS las calibraciones guardadas?\n\n"
                "Se pierden los puntos de las otras fuentes, que son los que "
                "permiten ver la linealidad."):
            return
        try:
            os.remove(RUTA_REGISTRO)
        except OSError:
            pass
        self.repintar()

    # ------------------------------------------------------------------
    def trabajo(self):
        import time
        fu = self.fuente()
        d = self.app.datos
        self.fase(0.0, 0.12)          # leer los XML (o sacarlos del cache)
        r = d.espectro()
        can, neto = r["canales"], r["neto"]
        manual = ((float(self.lo.get()), float(self.hi.get()))
                  if self.ventana.get() == "manual" else None)
        iv = d.ventana_ajuste(manual=manual)
        ncc = d.borde_estimado()
        B = float(ncc["edge_ch"])

        filas = []                       # (etiqueta, D, err, color)
        avisos = list(fu.avisos())
        res = {"fuente": fu, "iv": iv, "B": B, "rho": float(ncc["rho"]),
               "r": r, "filas": filas, "avisos": avisos, "presupuesto": [],
               "ppal": None, "fits": {}, "registro": {}}

        # Cuantos ajustes se van a correr. Se cuenta ANTES de empezar para que
        # la barra pueda decir "3 de 10" y estimar cuanto falta: es la unica
        # forma honesta de dar un ETA, porque cada ajuste cuesta parecido.
        ok_sim, msg_sim = d.compat_simulacion()
        usar_sim = d.hay_simulacion() and (ok_sim or d.forzar_sim)
        limites = [0.20 * B, 0.35 * B, 0.50 * B, 0.65 * B]
        n_total = 1 + (2 if usar_sim else 0)
        if self.sist.get():
            n_total += len(limites) + (3 if usar_sim else 0)
        hechos = [0]

        def paso(texto):
            self.progreso("[%d/%d] %s" % (hechos[0] + 1, n_total, texto),
                          hechos[0] / n_total)
            hechos[0] += 1

        self.fase(0.12, 1.0)

        # ---------------- 1) ruta analitica: siempre corre ----------------
        paso("ajuste analítico (Klein-Nishina)")
        if fu.doble:
            dob = ajustar_dos_bordes(can, neto, iv, E_gamma1=fu.lineas[0],
                                     E_gamma2=fu.lineas[1], ligar_origen=True,
                                     ligar_sigmas=False,
                                     B1_rango=(0.60 * B, 1.25 * B))
            D_kn, sig_kn, B_kn = dob["D"], dob["sigma1"], dob["B1"]
            etiq_kn = "Klein-Nishina, 2 bordes (método del paper)"
        else:
            uno = ajustar_borde_unico(can, neto, fu.lineas[0], iv, p0=ncc["p0"])
            D_kn, sig_kn, B_kn = (uno["edge_ch"] / fu.E_ref, uno["sigma"],
                                  uno["edge_ch"])
            etiq_kn = "Klein-Nishina, 1 borde (método del paper)"
        filas.append((etiq_kn, D_kn, None, "0.55"))
        res.update(D_kn=D_kn, sigma_kn=sig_kn, B_kn=B_kn)

        # los tres clasicos, que es contra lo que compara el paper (Fig. 4).
        # Van SIEMPRE, con o sin simulacion: son gratis y ubican al resto.
        clas = comparar_clasicos(can, neto, fu.E_ref, B_ref=B,
                                 errores=r["error"])
        res["clasicos"] = clas
        for nombre, color in (("50% de altura", "0.7"),
                              ("diferencia", "tab:olive"),
                              ("FWHM gaussiano", "tab:brown")):
            v = clas.get(nombre, {})
            if np.isfinite(v.get("D", float("nan"))):
                filas.append(("clásico · " + nombre, v["D"], None, color))

        # ---------------- 2) ruta con simulacion, si corresponde ----------
        res["msg_sim"] = msg_sim
        if d.hay_simulacion() and not ok_sim:
            avisos.append("simulación incompatible → " + msg_sim)
        if not usar_sim:
            avisos.append("sin plantilla Geant4 utilizable: el número que sale "
                          "es el analítico, que en el Co-60 de este proyecto "
                          "quedó +4.3% alto por no tener scattering múltiple "
                          "(§7.18). Es un límite del método, no del ajuste.")

        fits = {}
        if usar_sim:
            for rama in ("Lvis_keVee", "Edep_keV"):
                paso("modelo extendido, plantilla %s" % rama)
                fits[rama] = ajuste_extendido(d, r, iv, rama=rama)
            res["fits"] = fits
            filas.append(("Geant4 Edep → keV (Etapa 9)", fits["Edep_keV"]["D"],
                          fits["Edep_keV"]["err"].get("D"), "tab:green"))
            filas.append(("Geant4 Lvis → keVee (Etapa 9)",
                          fits["Lvis_keVee"]["D"],
                          fits["Lvis_keVee"]["err"].get("D"), "tab:blue"))

        ppal = fits.get(self.rama.get())
        res["ppal"] = ppal
        escala = ("keVee" if (ppal and self.rama.get() == "Lvis_keVee")
                  else "keV")
        D = float(ppal["D"]) if ppal else float(D_kn)
        err_stat = (float(ppal["err"].get("D", float("nan"))) if ppal
                    else float("nan"))

        # ---------------- 3) el presupuesto de error ----------------------
        presupuesto = []
        if np.isfinite(err_stat):
            presupuesto.append(("estadística (ajuste)", err_stat, False, True))
        if self.sist.get():
            lo_hi = iv[1]
            Ds = []
            for k, lo in enumerate(limites):
                paso("sistemático: ventana desde %.0f mV·ns" % lo)
                ivk = (float(lo), float(lo_hi))
                if ppal:
                    g = ajuste_extendido(d, r, ivk, rama=self.rama.get())
                    Ds.append(g["D"])
                elif fu.doble:
                    g = ajustar_dos_bordes(can, neto, ivk,
                                           E_gamma1=fu.lineas[0],
                                           E_gamma2=fu.lineas[1],
                                           ligar_origen=True,
                                           B1_rango=(0.60 * B, 1.25 * B))
                    Ds.append(g["D"])
                else:
                    g = ajustar_borde_unico(can, neto, fu.lineas[0], ivk,
                                            p0=ncc["p0"])
                    Ds.append(g["edge_ch"] / fu.E_ref)
            res["barrido_ventana"] = list(zip(limites, Ds))
            presupuesto.append(("ventana (límite inferior)",
                                0.5 * (max(Ds) - min(Ds)), False, False))

            if ppal:
                paso("sistemático: exceso como ley de potencia")
                g = ajuste_extendido(d, r, iv, rama=self.rama.get(),
                                     forma="pot")
                filas.append(("exceso como ley de potencia", g["D"],
                              g["err"].get("D"), "tab:orange"))
                presupuesto.append(("forma del exceso", abs(g["D"] - D),
                                    False, False))
                res["D_pot"] = g["D"]

                for etiqueta, w_fijo, color in (
                        ("w fijo = 1 (solo 1 scattering)", 1.0, "crimson"),
                        ("w fijo = %.3f (simulación)" % ppal["w_nat"],
                         ppal["w_nat"], "0.35")):
                    paso("sistemático: " + etiqueta)
                    g = ajuste_extendido(d, r, iv, rama=self.rama.get(),
                                         w_libre=False, w_fijo=w_fijo)
                    filas.append((etiqueta, g["D"], g["err"].get("D"), color))

            presupuesto.append(("definición del área (§7.35)",
                                self.ERR_AREA_REL * D, True, False))
        res["presupuesto"] = presupuesto

        sis = [v for _, v, _, est in presupuesto if not est]
        err_sist = float(np.sqrt(np.sum(np.square(sis)))) if sis else float("nan")
        err_tot = float(np.sqrt(np.nansum(np.square(
            [err_stat if np.isfinite(err_stat) else 0.0,
             err_sist if np.isfinite(err_sist) else 0.0]))))
        res.update(D=D, escala=escala, err_stat=err_stat, err_sist=err_sist,
                   err_total=err_tot,
                   ruta=("modelo extendido (Geant4 %s)" % self.rama.get()
                         if ppal else
                         ("Klein-Nishina (la simulación es de otra fuente)"
                          if d.hay_simulacion()
                          else "Klein-Nishina (sin simulación)")))

        if fu.nombre == "Co-60":
            filas.append(("referencia: Etapa 9 guardada (este SiPM)",
                          cal.D_LUZ, cal.D_LUZ_ERR_STAT, "0.75"))

        # ---------------- 4) al registro ----------------------------------
        entrada = {
            "fuente": fu.nombre, "lineas": list(fu.lineas),
            "E_ref": float(fu.E_ref), "bordes": [float(E) for E in fu.bordes],
            "D": D, "err_stat": err_stat, "err_sist": err_sist,
            "err_total": err_tot, "escala": escala, "ruta": res["ruta"],
            "D_kn": float(D_kn), "ventana": [float(iv[0]), float(iv[1])],
            "clasicos": {k: float(v.get("D", float("nan")))
                         for k, v in res.get("clasicos", {}).items()},
            "alpha": (float(ppal["alpha"]) if ppal else None),
            "fwhm_rel_ref": (float(ppal["fwhm_rel_ref"]) if ppal else
                             float(2.355 * sig_kn / max(B_kn, 1e-9))),
            "w": (float(ppal["w"]) if ppal else None),
            "chi2_red": (float(ppal["chi2_red"]) if ppal else None),
            "muestra": os.path.basename(d.archivo_muestra),
            "fondo": os.path.basename(d.archivo_fondo),
            "sim": (os.path.basename(d.archivo_sim) if usar_sim else ""),
            "huella": d.huella(),
            "fecha": time.strftime("%Y-%m-%d %H:%M"),
        }
        res["entrada"] = entrada
        res["registro"] = guardar_en_registro(
            "%s|%s" % (fu.nombre, d.huella()), entrada)
        return res

    # ------------------------------------------------------------------
    def pintar(self, res):
        fu = res["fuente"]
        ax1 = self.fig.add_subplot(121)
        ax2 = self.fig.add_subplot(122)

        # ---- izquierda: todas las estimaciones sobre los mismos datos ----
        filas = res["filas"]
        vals = np.array([f[1] for f in filas], float)
        errs = np.array([f[2] if f[2] else 0.0 for f in filas], float)
        ypos = np.arange(len(filas))
        ax1.errorbar(vals, ypos, xerr=errs, fmt="o", ecolor="0.4", capsize=3,
                     linestyle="none", markersize=0)
        for yy, v, c in zip(ypos, vals, [f[3] for f in filas]):
            ax1.plot(v, yy, "o", color=c, markersize=8)
        if np.isfinite(res["err_sist"]):
            ax1.axvspan(res["D"] - res["err_sist"], res["D"] + res["err_sist"],
                        color="tab:blue", alpha=0.10,
                        label="D ± sistemático")
            ax1.legend(fontsize=8)
        ax1.set_yticks(ypos)
        ax1.set_yticklabels([f[0] for f in filas], fontsize=8)
        ax1.invert_yaxis()
        ax1.set_xlabel("D [mV·ns por keV depositado, o por keVee de luz]")
        ax1.set_title("%s — todas las estimaciones, mismos datos" % fu.nombre,
                      fontsize=10)
        ax1.grid(alpha=0.3, axis="x")

        # ---- derecha: la recta, con las otras fuentes medidas ------------
        reg = res.get("registro", {})
        mejores = {}
        for e in reg.values():
            E = float(e.get("E_ref", 0.0))
            if E <= 0 or float(e.get("D", 0.0)) <= 0:
                continue
            clave = (e.get("fuente"), round(E, 1))
            previa = mejores.get(clave)
            if (previa is None
                    or (e.get("err_total") or 9e9) < (previa.get("err_total") or 9e9)):
                mejores[clave] = e
        puntos = []
        for e in sorted(mejores.values(), key=lambda e: e.get("E_ref", 0)):
            E = float(e["E_ref"])
            Dv = float(e["D"])
            sg = float(e.get("err_total") or 0.0)
            puntos.append((E, Dv * E, max(sg * E, 1e-9), e))
        res["puntos"] = puntos
        E_top = max([1.25 * max(fu.bordes)] + [1.15 * p[0] for p in puntos])
        E = np.linspace(0, E_top, 200)
        ax2.plot(E, res["D"] * E, "tab:blue", lw=2,
                 label="esta corrida: %.4f mV·ns/%s" % (res["D"], res["escala"]))
        if np.isfinite(res["err_total"]):
            ax2.fill_between(E, (res["D"] - res["err_total"]) * E,
                             (res["D"] + res["err_total"]) * E,
                             color="tab:blue", alpha=0.12)
        for E_b in fu.bordes:
            ax2.plot(E_b, res["D"] * E_b, "k*", markersize=13)
            ax2.annotate("borde %.0f" % E_b, (E_b, res["D"] * E_b),
                         textcoords="offset points", xytext=(-58, 6), fontsize=8)
        for Ep, Ap, sp, e in puntos:
            propio = e.get("huella") == res["entrada"].get("huella")
            ax2.errorbar(Ep, Ap, yerr=sp, fmt="s" if not propio else "o",
                         color="tab:red" if not propio else "tab:blue",
                         markersize=6, capsize=3,
                         label=None if propio else "%s (%s)" % (e.get("fuente"),
                                                                e.get("escala")))
        if len(puntos) >= 2:
            Dg, sg = recta_por_origen([(p[0], p[1], p[2]) for p in puntos])
            ax2.plot(E, Dg * E, "k--", lw=1.2,
                     label="todas las fuentes: %.4f ± %.4f" % (Dg, sg))
            res["D_global"] = (Dg, sg)
        ax2.set_xlabel("Energía del borde [%s]" % res["escala"])
        ax2.set_ylabel("Área del pulso [mV·ns]")
        ax2.set_title("Recta de calibración", fontsize=10)
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)

        self.escribir(self._informe(res) + "\n\n" + self.explicacion.strip())

    # ------------------------------------------------------------------
    def _informe(self, res):
        fu = res["fuente"]
        L = ["FUENTE: " + fu.descripcion()]
        if fu.nota():
            L.append("   " + fu.nota())
        L.append("ruta usada: " + res["ruta"])
        L.append("ventana [%.0f, %.0f] mV·ns   (borde NCC B = %.0f, ρ = %.4f)"
                 % (res["iv"][0], res["iv"][1], res["B"], res["rho"]))
        L.append("")
        L.append("D = %.4f mV·ns/%s   →  %.1f %s por 1000 mV·ns"
                 % (res["D"], res["escala"], 1000 / res["D"], res["escala"]))
        est = ("± %.4f (est)" % res["err_stat"] if np.isfinite(res["err_stat"])
               else "sin error estadístico (ajuste analítico)")
        sist = ("± %.4f (sist)" % res["err_sist"] if np.isfinite(res["err_sist"])
                else "sistemáticos SIN barrer (destildado)")
        L.append("    %s  %s   →  total ± %.4f (%.1f%%)"
                 % (est, sist, res["err_total"],
                    100 * res["err_total"] / res["D"]))
        ppal = res["ppal"]
        if ppal:
            L.append("α = %.3f √keV  →  σ(%.0f) = %.0f keV, FWHM = %.1f%%"
                     % (ppal["alpha"], ppal["E_ref"], ppal["sigma_E_ref"],
                        100 * ppal["fwhm_rel_ref"]))
            L.append("w = %.3f ± %.3f  (la simulación trae %.3f)"
                     % (ppal["w"], ppal["err"].get("w", float("nan")),
                        ppal["w_nat"]))
            L.append("χ²/ndf = %.3f  (%d gdl)"
                     % (ppal["chi2_red"], ppal["ndf"]))
            otra = ("Edep_keV" if self.rama.get() == "Lvis_keVee"
                    else "Lvis_keVee")
            if otra in res["fits"]:
                L.append("la otra escala (%s): D = %.4f   cociente %.3f"
                         % (otra, res["fits"][otra]["D"],
                            res["D"] / res["fits"][otra]["D"]))
        else:
            L.append("σ/B = %.3f  →  FWHM = %.1f%% en el borde"
                     % (res["sigma_kn"] / max(res["B_kn"], 1e-9),
                        235.5 * res["sigma_kn"] / max(res["B_kn"], 1e-9)))
        if res["ppal"]:
            L.append("Klein-Nishina daría D = %.4f  (%+.1f%% respecto del "
                     "elegido)" % (res["D_kn"],
                                   100 * (res["D_kn"] / res["D"] - 1)))

        if res.get("clasicos"):
            L.append("")
            L.append("MÉTODOS CLÁSICOS (la comparación de la Fig. 4 del paper)")
            for nombre, v in res["clasicos"].items():
                if not np.isfinite(v.get("D", float("nan"))):
                    L.append(f"   {nombre:<15} no ajustó")
                    continue
                L.append(f"   {nombre:<15} D = {v['D']:.4f}  "
                         f"({100*(v['D']/res['D'] - 1):+.1f}% respecto del "
                         f"elegido)")
            L.append("   Los tres sobrestiman por construcción: ninguno "
                     "contiene la acumulación de Klein-Nishina antes del "
                     "borde. Sobre sintético con verdad conocida el sesgo va "
                     "de +0.4% a +9% según la resolución.")

        if res["presupuesto"]:
            L.append("")
            L.append("PRESUPUESTO DE ERROR")
            for nombre, val, arrastrado, _est in res["presupuesto"]:
                L.append("   %-28s ± %.4f  (%.2f%%)%s"
                         % (nombre, val, 100 * val / res["D"],
                            "  ← arrastrado" if arrastrado else ""))
        if res.get("barrido_ventana"):
            L.append("")
            L.append("BARRIDO DE LA VENTANA (límite inferior)")
            for lo, Dv in res["barrido_ventana"]:
                L.append("   lo = %6.0f  →  D = %.4f" % (lo, Dv))

        reg = res.get("registro", {})
        if reg:
            L.append("")
            L.append("REGISTRO DE FUENTES MEDIDAS (.cache/calibraciones.json)")
            for clave, e in sorted(reg.items(),
                                   key=lambda kv: kv[1].get("E_ref", 0)):
                L.append("   %-8s borde %7.1f keV   D = %.4f ± %.4f  [%s]"
                         % (e.get("fuente", "?"), e.get("E_ref", 0),
                            e.get("D", 0), e.get("err_total") or 0,
                            e.get("escala", "?")))
                L.append("            %s · %s" % (e.get("fecha", ""),
                                                  e.get("ruta", "")))
            if "D_global" in res:
                Dg, sg = res["D_global"]
                L.append("   → todas juntas por el origen: D = %.4f ± %.4f"
                         % (Dg, sg))
                disp = max(abs(p[3]["D"] / Dg - 1)
                           for p in res.get("puntos", []))
                L.append("     dispersión máxima entre fuentes: %.1f%%  ← el "
                         "test de linealidad" % (100 * disp))
            else:
                L.append("   (una sola fuente: la linealidad se está "
                         "ASUMIENDO, no verificada)")

        if res["avisos"]:
            L.append("")
            L.append("AVISOS")
            for a in res["avisos"]:
                L.append("   ⚠ " + a)
        return "\n".join(L)


# ======================================================================
#  Glosario de parametros
# ======================================================================
# (simbolo, nombre, unidades, valor en este proyecto, donde vive, explicacion)
GLOSARIO = [
    ("canal / área", "Área del pulso integrada", "mV·ns", "0 – 1900",
     "spectrum._area_pulso",
     "El observable medido. NO es un canal de ADC: es la integral en el tiempo del "
     "pulso del SiPM, con la línea base restada. Es proporcional a la LUZ recogida, "
     "y por eso es proporcional a la energía en el régimen lineal. Todo el eje "
     "horizontal de este proyecto está en mV·ns."),

    ("E_γ", "Energía del fotón incidente", "keV", "1173.2 y 1332.5",
     "physics.compton_edge_energy",
     "Las dos líneas gamma del Co-60. Son datos conocidos de la fuente, NO se "
     "ajustan. Cada una produce su propio borde Compton."),

    ("E_max", "Energía del borde Compton", "keV", "963.4 y 1118.1",
     "physics.compton_edge_energy",
     "Máxima energía que el fotón puede transferir a un electrón, en "
     "retrodispersión (θ=180°):  E_max = 2α·E_γ/(1+2α).  Es el ANCLA de toda la "
     "calibración: una energía conocida sin ajustar ningún parámetro.\n\n"
     "⚠️ ERRATUM del paper: la Ec. 1 impresa da E_γ/(1+2α), que es la energía del "
     "fotón dispersado, no la del electrón."),

    ("α  (Klein-Nishina)", "Energía del fotón en unidades de m_e c²", "adimensional",
     "2.30 y 2.61", "physics.dsigma_dE",
     "α = E_γ / m_e c², con m_e c² = 510.999 keV. Aparece en la sección eficaz de "
     "Klein-Nishina.\n\n⚠️ OJO: en este proyecto la letra α se usa para TRES cosas "
     "distintas — esta, el coeficiente de enfriamiento del SA, y el parámetro de "
     "resolución α·√E. No son lo mismo."),

    ("A, A1, A2", "Amplitud de cada componente", "cuentas (escala)", "A2/A1 ≈ 0.75",
     "model.forward_model, fit_double",
     "Escalan la altura de la lineshape. Absorben el prefactor π r0²/(m_ec²) de la "
     "Ec. 10: para comparar formas solo importa la escala relativa.\n\n"
     "Son LINEALES en el modelo, así que se resuelven exacto por mínimos cuadrados "
     "(nnls) en vez de buscarlas — eso es lo que hace `_amplitudes_optimas`.\n\n"
     "⚠️ El coeficiente de Pearson (NCC) es invariante a la amplitud: con ese "
     "objetivo, A es INDETERMINABLE. Verificado: multiplicar A por 7.3 cambia el "
     "costo NCC en exactamente 0."),

    ("B  (edge_ch)", "Canal donde cae el borde Compton", "mV·ns", "1215",
     "model.forward_model, fit_single",
     "EL parámetro que se busca. Junto con E_max da el par (canal, energía) que "
     "calibra el detector. En el paper es 'B'; en el código, `edge_ch`.\n\n"
     "El mapa canal↔energía pasa por el ORIGEN: E_e = canal · E_max/B. Área cero "
     "corresponde a energía cero, sin offset."),

    ("B1, B2  (CE_ch,2)", "Canal de cada uno de los dos bordes", "mV·ns",
     "1215 y 1410", "fit_double.componentes",
     "Para el Co-60 hay dos bordes solapados. El segundo se amarra al primero con "
     "la Ec. 12:  CE_ch,2 = B1 + D·ΔE.\n\n"
     "Su separación (195 mV·ns) es ahora 1.8σ (σ1 = 109 mV·ns), contra 0.93σ con "
     "el área vieja. Mejoró, pero para tener un segundo punto REAL de la recta "
     "hace falta otra fuente."),

    ("σ, σ1, σ2", "Resolución del detector, en canal", "mV·ns", "109 y 118",
     "model._kernel_gaussiano",
     "Ancho de la gaussiana con que se convoluciona la lineshape teórica (Ec. 9). "
     "Es lo que convierte el escalón vertical de la teoría en la rampa que medís.\n\n"
     "El kernel está NORMALIZADO (suma 1), así que conserva el área: subir σ "
     "ensancha sin cambiar las cuentas totales, y A y σ quedan desacoplados.\n\n"
     "Físicamente σ2 > σ1 porque la resolución escala como √E; si un ajuste te da "
     "σ1 > σ2, está usando libertad que la física no le da."),

    ("σ/B  (σ_rel)", "Resolución relativa", "adimensional", "0.090  (21.1% FWHM)",
     "ncc.plantilla_referencia",
     "Resolución en unidades del borde. Se usa relativa porque así la plantilla del "
     "NCC no depende de la ganancia del detector.\n\n"
     "FWHM = 2.355·σ/B. El 42% medido es alto incluso para plástico+SiPM, y NO es "
     "artefacto del modelo: las tres plantillas (Klein-Nishina, Geant4 completo, "
     "Geant4 de scattering simple) piden 38-44%."),

    ("D", "Coeficiente de conversión energía→canal", "mV·ns/keV", "1.261 ± 0.036",
     "fit_double, sim_template, calibracion",
     "LA CALIBRACIÓN. Es la pendiente de la recta que convierte área de pulso en "
     "energía. D = 1/pendiente inversa: 4.23 mV·ns/keV = 236.5 keV por 1000 mV·ns.\n\n"
     "Con B1 y D ambos libres el mapa deja de pasar por el origen, lo que "
     "contradice la Ec. 8. Comparar D contra B1/E_max,1 es el TEST de esa "
     "hipótesis: salió compatible (offset -2.3%)."),

    ("D_luz vs D_edep", "Las dos escalas de calibración", "mV·ns/keVee | mV·ns/keV",
     "4.229 y 4.045", "calibracion.area_a_energia",
     "El SiPM ve LUZ, no energía. Según tu simulación, Lvis/Edep = 0.954 arriba de "
     "700 keV y 0.863 a 50-100 keV: eso es el quenching de Birks.\n\n"
     "• keVee (luz): para la rama de NEUTRONES, donde los protones de retroceso se "
     "comparan en keV electrón-equivalente.\n"
     "• keV (energía depositada): para gammas.\n"
     "Difieren 4.6%. No mezclarlas."),

    ("ΔE", "Separación entre los dos bordes", "keV", "154.7",
     "fit_double.delta_E",
     "Entra en la Ec. 12 para amarrar el segundo borde al primero.\n\n"
     "⚠️ ERRATUM del paper: dice 'energy difference between gamma lines' (159.3 "
     "keV). Está mal: D convierte energía del ELECTRÓN a canal, y lo que se separa "
     "en canal son los BORDES. Con el ΔE equivocado, D sale ~3% corrido."),

    ("Θ", "Vector de parámetros del modelo de dos bordes", "—",
     "{A1,A2,B1,σ1,σ2,D}", "fit_double.ajustar_dos_bordes",
     "Los 6 parámetros de la Ec. 11.\n\n"
     "⚠️ Está SUB-DETERMINADO cuando los bordes se solapan: el optimizador global "
     "encuentra soluciones con mejor χ² pero sin sentido físico (offset de -5928 "
     "mV·ns, o directamente A1=0). Hay que restringirlo con física."),

    ("ρ  (Pearson)", "Coeficiente de correlación normalizada", "−1 a 1", "0.9918",
     "ncc.pearson  (Ec. 15)",
     "Mide cuánto se parecen en FORMA los datos y la plantilla teórica. Resta la "
     "media y divide por la desviación de cada señal, así que ignora offset y "
     "amplitud; y como ambas se reinterpolan a 100 puntos, también ignora la escala "
     "horizontal.\n\n"
     "Esa triple invariancia es lo que le permite encontrar el borde sin conocer la "
     "ganancia — que es justo lo que queremos medir."),

    ("S1, S2", "Las dos señales que compara el NCC", "u.a. (100 puntos)", "—",
     "ncc.plantilla_referencia",
     "S2 = plantilla teórica (Klein-Nishina convolucionada), 100 puntos del plateau "
     "hasta pasado el borde. S1 = una ventana candidata del espectro medido, "
     "reinterpolada a los mismos 100 puntos.\n\n"
     "Ganá el intervalo con mayor ρ. Eso reemplaza al 'elegir el intervalo a ojo'."),

    ("inicio_rel, fin_rel", "Extensión de la plantilla, en unidades de E_max",
     "adimensional", "0.45 y ~1.0+4.5σ_rel", "ncc.plantilla_referencia",
     "Dónde empieza y termina S2. 0.45·E_max cae en el plateau; el final se pone "
     "4.5σ por encima del borde, donde la cola ya cayó a cero.\n\n"
     "Como el mapa pasa por el origen, fijado B la ventana queda DETERMINADA: "
     "[0.45·B, fin_rel·B]. Por eso el barrido correcto es 1D sobre B y no 2D sobre "
     "(inicio, fin) como propone el paper."),

    ("α_res, β, a_rel", "Parámetros del modelo de resolución", "√keV, keV, —",
     "2.79, 0 y 0",
     "sim_template.sigma_energia",
     "σ_E(E) = √(α²·E + β²).  El término α²E es estadística de fotoelectrones "
     "(σ ~ √N_pe, Tabla I del paper); β es el piso constante de ruido electrónico.\n\n"
     "Medido: α = 5.73 √keV, o sea σ(963 keV) = 178 keV. β salió compatible con 0."),

    ("χ²/ndf", "Bondad del ajuste por grado de libertad", "adimensional", "1.008",
     "sim_template.ajustar_forward_folding",
     "Suma de residuos al cuadrado PONDERADOS por la barra de error de cada bin, "
     "dividida por los grados de libertad. Cerca de 1 = el modelo describe los "
     "datos dentro del ruido estadístico.\n\n"
     "Sin pesos (mínimos cuadrados planos) los bins del plateau dominan y NO se "
     "puede poner barra de error a nada. El error del neto es Poisson propagado por "
     "la resta de fondo: √(N_m·k_m² + N_f·k_f²)."),

    ("T, T_init", "Temperatura del recocido simulado", "unidades del COSTO",
     "auto (2.1e9 ó 3.19)", "sa.recocido_simulado",
     "Controla cuánto tolera el algoritmo empeorar: acepta una subida ΔE con "
     "probabilidad exp(−ΔE/T).\n\n"
     "⚠️ El 'T_init = 2000' del paper no es transferible: T se compara contra ΔE, "
     "así que depende de las UNIDADES del costo. Medido acá: hace falta T ≈ 2.1e9 "
     "con el costo de mínimos cuadrados y ≈ 3.19 con el costo NCC. Nueve órdenes de "
     "magnitud. Por eso el código lo calibra solo."),

    ("α_enfriamiento", "Coeficiente de enfriamiento geométrico", "adimensional",
     "0.95 – 0.97", "sa.recocido_simulado",
     "T_{k+1} = α·T_k. El paper usa 0.95 y admite [0.85, 0.99]. Más cerca de 1 = "
     "enfría más lento = busca más (y tarda más).\n\n"
     "No confundir con el α de Klein-Nishina ni con el α de la resolución."),

    ("NCompton", "Número de dispersiones Compton por evento", "entero",
     "25% tiene >1", "simulación Geant4",
     "Rama de tu .root. 640082 eventos con 1 sola dispersión, 145497 con 2, 35072 "
     "con 3...\n\n"
     "El scattering MÚLTIPLE deposita más energía que el borde (hasta 1332 keV) y "
     "rellena la zona por encima. Klein-Nishina no lo contiene, así que un modelo "
     "analítico tiene que correr el borde hacia arriba para explicar esas cuentas: "
     "por eso sobreestima D un ~6%."),

    ("QF / Lvis / Edep", "Quenching y luz visible (simulación)", "— / keVee / keV",
     "QF medio 0.926", "simulación Geant4",
     "Edep = energía depositada. Lvis = luz producida, ya con el quenching de Birks "
     "aplicado. QF = Lvis/Edep evento a evento.\n\n"
     "Satura en 0.954 sobre 700 keV y baja a 0.863 en 50-100 keV. Es la razón de "
     "que existan dos escalas de calibración."),

    # ---------------- Etapa 8: MCMC ----------------
    ("posterior", "P(parámetros | datos)", "densidad de probabilidad",
     "3 parámetros: A, D, α", "mcmc.log_posterior_sim",
     "Toda la información que los datos dan sobre los parámetros, como una nube en "
     "el espacio de parámetros. Con priors planos es proporcional a exp(−χ²/2): la "
     "misma χ² del ajuste, leída como probabilidad en vez de como algo a minimizar."
     "\n\nEl ajuste da UN punto y una elipse aproximada; el posterior da la nube "
     "entera, incluso si está torcida o tiene dos máximos. La mediana y los "
     "percentiles 16/84 son el resultado a reportar."),

    ("R̂", "Estadístico de Gelman-Rubin", "adimensional", "1.000 – 1.005",
     "mcmc.gelman_rubin",
     "Compara la dispersión ENTRE cadenas con la dispersión DENTRO de cada cadena. "
     "Si varias cadenas que arrancaron en puntos distintos terminan describiendo la "
     "misma nube, las dos varianzas coinciden y R̂ → 1. Criterio usual: R̂ < 1.01."
     "\n\n⚠️ Sólo sirve si los arranques fueron DISPERSOS. Cuatro cadenas desde el "
     "mismo punto dan R̂ ≈ 1 aunque ninguna haya explorado nada."),

    ("τ_int / n_eff", "Autocorrelación y muestras independientes", "pasos / muestras",
     "τ ≈ 20, n_eff ≈ 2000", "mcmc.autocorr_integrada",
     "Los pasos de una cadena NO son independientes: cada uno arranca donde quedó el "
     "anterior. τ_int mide cuántos pasos hay que dar para olvidar de dónde se venía, "
     "y n_eff = n/τ es cuántas muestras independientes hay de verdad.\n\n"
     "Es lo que hace que 48000 pasos valgan como ~2000 muestras. Con n_eff de unos "
     "pocos cientos los percentiles ya son estables."),

    ("aceptación (MH)", "Fracción de propuestas aceptadas", "adimensional",
     "0.31", "mcmc.metropolis",
     "Debería caer en 0.2–0.5. Muy alta = los pasos son tan chicos que la cadena "
     "avanza a paso de hormiga; muy baja = son tan grandes que casi siempre caen en "
     "zonas de probabilidad despreciable y la cadena se queda quieta.\n\n"
     "Acá se logra sola porque la propuesta usa la covarianza del ajuste escalada "
     "por 2.4²/d, que es el óptimo teórico para un objetivo gaussiano."),

    ("w  (mezcla)", "Fracción de eventos de scattering SIMPLE", "adimensional",
     "0.796 ± 0.070", "modelo_ext.modelo_ext",
     "La plantilla es w·h_1scattering + (1−w)·h_múltiple. Con w = 1 es Klein-Nishina "
     "(un solo scattering, la plantilla del paper); con w = 0.751 es la simulación "
     "Geant4 tal cual sale.\n\n"
     "Antes esto era el 'sistemático de plantilla' (±2.7%), estimado comparando tres "
     "ajustes distintos. Ahora es un parámetro y los datos lo miden: piden 0.796 ± "
     "0.070, o sea COINCIDEN con la simulación (0.6σ).\n\n"
     "Que w esté realmente medido y no sea decorativo se ve fijándolo: en 0.751 el "
     "χ² no cambia (Δχ² = 0), en 1 (la plantilla analítica del paper) empeora "
     "Δχ² = +13 y D sube 3.4%.\n\n"
     "corr(D, w) = +0.96: por eso la ignorancia sobre la plantilla ya está dentro "
     "del error de D. Ese es todo el punto."),

    ("C, λ  (exceso)", "Componente empírica a bajo canal", "cuentas, mV·ns",
     "125 y 365", "modelo_ext.modelo_ext",
     "C·exp(−ch/λ): cuentas que los datos tienen y la simulación no, abajo de "
     "~800 mV·ns. Lo más probable es dispersión en el material de alrededor y el "
     "portafuente, que la simulación no modela.\n\n"
     "Modelarlo en vez de cortar la ventana permite ajustar desde 300 mV·ns y "
     "mata la deriva de D con el límite inferior (§7.32 tenía 3.5%; ahora 1.5% "
     "sin tendencia).\n\n"
     "⚠️ Su FORMA es el sistemático dominante que queda: pasar de exponencial a "
     "ley de potencia mueve D +2.35%, y las dos ajustan igual de bien (χ²/ndf "
     "1.008 contra 1.015). Eso no se cierra ajustando, hay que medirlo."),

    ("t_vivo", "Tiempo VIVO de la adquisición", "s", "283 s de 699 (Co-60)",
     "spectrum.tiempo_vivo",
     "El `<Scaler1>` de cada evento es la tasa de trigger REAL de la placa. En el "
     "Co-60 marca 170 Hz mientras solo se guardaron 71.5 ev/s: el DRSOSC vuelca "
     "~35 kB de XML por evento y no da abasto.\n\n"
     "Importa porque la resta de fondo se normaliza por tiempo: con tiempo de "
     "PARED el fondo entraba con k = 0.709 cuando corresponde 0.287, o sea 2.5x "
     "de más.\n\n"
     "Las marcas de <Time> NO sirven para esto: el DRSOSC las escribe en ráfagas "
     "(4124 eventos del Co-60 comparten timestamp al ms)."),

    ("compuerta / onset", "Cómo se integra el área del pulso", "ns",
     "[onset−6, onset+55]", "spectrum._area_pulso_gate",
     "Onset = cruce del 10% de la amplitud retrocediendo desde el pico. La línea "
     "base es la mediana hasta 4 ns ANTES del onset, y la integral va en una "
     "compuerta alrededor del pulso.\n\n"
     "El método viejo promediaba la línea base sobre los primeros 40 ns e "
     "integraba los 200 ns completos. Como el pulso arranca a los 25 ns, se "
     "restaba parte del propio pulso multiplicada por T/T_base ≈ 5: el 26% de los "
     "eventos salía con el área de signo invertido y el abs() los reciclaba.\n\n"
     "Verificado: el pulso medio devuelve el 100% de su área en 80 ns y "
     "área/amplitud es constante al 1% entre 20 y 220 mV. Ver la pestaña 9."),

    ("sesgo de Neyman", "χ² pesado con σ del dato vs del modelo", "% en D",
     "−0.5% a 25 cuentas/bin", "spectrum.varianza_predicha",
     "Si el peso de cada bin se calcula con σ = √N_observado, un bin que fluctuó "
     "BAJO recibe una barra más chica, o sea MÁS peso. El ajuste se deja arrastrar "
     "hacia abajo. El sesgo va como 1/N por bin.\n\n"
     "Medido con verdad conocida: −1.9% con 2.6 cuentas/bin, −0.7% con 26, −0.05% "
     "con 258. Este experimento tiene ~25 cuentas/bin, o sea el sesgo era comparable "
     "al error estadístico que reportábamos.\n\n"
     "La cura es el χ² de Pearson: calcular σ con el MODELO e iterar. Sobre 20 "
     "realizaciones la cobertura del 68% pasó de 8/20 a 15/20."),
]


# ojo: _normalizar pasa a minuscula ANTES de sustituir, asi que las mayusculas
# griegas hay que listarlas en su forma minuscula (Θ -> θ, Δ -> δ)
_ALIAS = {"σ": "sigma", "ρ": "rho", "α": "alpha", "β": "beta", "χ": "chi",
          "δ": "delta", "θ": "theta", "γ": "gamma", "·": " ", "²": "2",
          "á": "a", "é": "e", "í": "i", "ó": "o", "ú": "u", "ñ": "n",
          "τ": "tau", "̂": " hat"}       # "R̂" -> "r hat"; ojo: el acento va suelto


def _normalizar(texto):
    """Para que buscar 'sigma' encuentre las entradas escritas con σ, etc."""
    t = str(texto).lower()
    for k, v in _ALIAS.items():
        t = t.replace(k, v)
    return t


class TabParametros(ttk.Frame):
    titulo = "? Parámetros"
    auto = False

    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.calculada = True
        self._corriendo = False

        barra = ttk.Frame(self)
        barra.pack(side="top", fill="x", padx=6, pady=6)
        ttk.Label(barra, text="Filtrar:").pack(side="left")
        self.filtro = tk.StringVar()
        e = ttk.Entry(barra, textvariable=self.filtro, width=30)
        e.pack(side="left", padx=6)
        self.filtro.trace_add("write", lambda *_: self.llenar())
        ttk.Label(barra, text="(buscá 'sigma', 'D', 'paper', 'erratum'...)",
                  foreground="#666").pack(side="left")

        pw = ttk.PanedWindow(self, orient="vertical")
        pw.pack(fill="both", expand=True, padx=6, pady=(0, 6))

        marco = ttk.Frame(pw)
        cols = ("simbolo", "nombre", "unidades", "valor", "donde")
        self.tree = ttk.Treeview(marco, columns=cols, show="headings", height=14)
        for c, txt, anchoc in [("simbolo", "Símbolo", 150),
                               ("nombre", "Qué es", 300),
                               ("unidades", "Unidades", 130),
                               ("valor", "Valor acá", 150),
                               ("donde", "Dónde vive", 240)]:
            self.tree.heading(c, text=txt)
            self.tree.column(c, width=anchoc, anchor="w")
        sb = ttk.Scrollbar(marco, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        pw.add(marco, weight=3)

        det = ttk.LabelFrame(pw, text="Detalle")
        self.txt = tk.Text(det, wrap="word", relief="flat", height=10,
                           font=("Segoe UI", 10), background="#f7f7f7",
                           padx=8, pady=8)
        self.txt.pack(fill="both", expand=True)
        pw.add(det, weight=2)

        self.tree.bind("<<TreeviewSelect>>", self.al_seleccionar)
        self.llenar()

    def llenar(self):
        f = _normalizar(self.filtro.get().strip())
        self.tree.delete(*self.tree.get_children())
        self.filas = [g for g in GLOSARIO
                      if not f or any(f in _normalizar(campo) for campo in g)]
        for i, g in enumerate(self.filas):
            self.tree.insert("", "end", iid=str(i), values=g[:5])
        if self.filas:
            self.tree.selection_set("0")

    def al_seleccionar(self, _evt=None):
        sel = self.tree.selection()
        if not sel:
            return
        g = self.filas[int(sel[0])]
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        self.txt.insert("1.0", f"{g[0]}  —  {g[1]}\n", "titulo")
        self.txt.insert("end", f"unidades: {g[2]}     valor en este proyecto: {g[3]}"
                               f"     código: {g[4]}\n\n", "meta")
        self.txt.insert("end", g[5])
        self.txt.tag_configure("titulo", font=("Segoe UI", 12, "bold"))
        self.txt.tag_configure("meta", foreground="#555")
        self.txt.configure(state="disabled")

    def repintar(self):
        pass


# ======================================================================
#  aplicacion
# ======================================================================
class App(tk.Tk):
    def __init__(self):
        activar_dpi_awareness()          # antes de crear la ventana
        super().__init__()
        self.title("Calibración por Compton edge — visor por etapas")
        pref = cargar_preferencias()
        try:
            fuente = Fuente(pref["fuente"], pref["lineas"])
        except Exception:
            fuente = Fuente(FUENTE_DEF)
        self.datos = Datos(pref["muestra"], pref["fondo"], pref["sim"],
                           pref["arbol"], fuente)
        self.cola = queue.Queue()
        self.estado = tk.StringVar(value="listo")
        self.pct = tk.StringVar(value="")
        self._trabajos = set()
        self._t0 = 0.0
        self._modo = "indeterminado"
        self.configurar_nitidez()

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True)
        self.pestanas = []
        for clase in (TabArchivos, TabEspectro, TabFisica, TabForward,
                      TabUnBorde, TabNCC, TabAnalitico, TabSA, TabSimulacion,
                      TabAjusteGlobal, TabArea, TabModeloExt, TabCalibracion,
                      TabParametros):
            p = clase(self.nb, self)
            self.nb.add(p, text=clase.titulo)
            self.pestanas.append(p)
        self.actualizar_titulo()
        self.propagar_fuente()

        barra = ttk.Frame(self)
        barra.pack(side="bottom", fill="x")
        ttk.Label(barra, textvariable=self.estado, anchor="w",
                  relief="sunken").pack(side="left", fill="x", expand=True)
        self.barra_prog = ttk.Progressbar(barra, mode="determinate",
                                          maximum=100.0,
                                          length=int(180 * self.escala))
        self.barra_prog.pack(side="left", padx=(8, 4))
        ttk.Label(barra, textvariable=self.pct, width=20,
                  anchor="w").pack(side="left")
        ttk.Label(barra, text="nitidez").pack(side="left", padx=(10, 3))
        self.var_nitidez = tk.StringVar(value=f"{self.escala:.2f}")
        cb = ttk.Combobox(barra, textvariable=self.var_nitidez, width=6,
                          state="readonly",
                          values=["1.00", "1.25", "1.50", "1.75", "2.00", "2.50"])
        cb.pack(side="left", padx=(0, 8))
        cb.bind("<<ComboboxSelected>>", self.cambiar_nitidez)

        self.nb.bind("<<NotebookTabChanged>>", self.al_cambiar)
        self.after(100, self.revisar_cola)
        self.after(200, lambda: self.al_cambiar(None))

    # ------------------------------------------------------------------
    def progreso_iniciar(self, pest):
        """
        Arranca la barra de abajo cuando una pestana lanza su calculo.

        Empieza INDETERMINADA (la barrita que va y viene) porque la mayoria de
        los pasos no sabe de antemano cuanto falta, y para lo que se le pide
        —saber que el programa no se colgo— eso ya alcanza. En cuanto llega el
        primer `progreso(..., frac)` se pasa sola a determinada, con porcentaje
        y tiempo restante estimado.
        """
        import time
        if getattr(self, "barra_prog", None) is None:
            return
        self._trabajos.add(pest)
        self._t0 = time.time()
        self._modo = "indeterminado"
        self.barra_prog.configure(mode="indeterminate")
        self.barra_prog.start(50)
        self.pct.set("trabajando…")

    def progreso_actualizar(self, frac):
        """Pasa la barra a determinada y estima cuanto falta."""
        import time
        if getattr(self, "barra_prog", None) is None or not self._trabajos:
            return
        if frac is None:
            return
        frac = max(0.0, min(float(frac), 1.0))
        if self._modo != "determinado":
            self.barra_prog.stop()
            self.barra_prog.configure(mode="determinate", maximum=100.0)
            self._modo = "determinado"
        self.barra_prog["value"] = 100.0 * frac
        transcurrido = time.time() - self._t0
        # el ETA sale de una regla de tres, asi que antes de tener algo de
        # historia (3% y 2 s) es puro ruido y conviene no mostrarlo.
        if frac > 0.03 and transcurrido > 2.0:
            self.pct.set("%.0f%% · faltan ~%s"
                         % (100 * frac,
                            formato_tiempo(transcurrido * (1 - frac) / frac)))
        else:
            self.pct.set("%.0f%%" % (100 * frac))

    def progreso_terminar(self, pest):
        """Apaga la barra, pero solo cuando NO queda ninguna tarea corriendo."""
        if getattr(self, "barra_prog", None) is None:
            return
        self._trabajos.discard(pest)
        if self._trabajos:
            return
        self.barra_prog.stop()
        self.barra_prog.configure(mode="determinate", maximum=100.0)
        self.barra_prog["value"] = 0.0
        self._modo = "indeterminado"
        self.pct.set("")

    def actualizar_titulo(self):
        sim = (os.path.basename(self.datos.archivo_sim)
               if self.datos.hay_simulacion() else "sin simulación")
        self.title("Calibración por Compton edge — " + self.datos.fuente.nombre
                   + "   [" + os.path.basename(self.datos.archivo_muestra)
                   + "  −  " + os.path.basename(self.datos.archivo_fondo)
                   + "  ·  " + sim + "]")

    def propagar_fuente(self):
        """
        Avisa a las pestanas que cambio la fuente, en el hilo principal.

        Ademas refresca el texto de las solapas: la 6 dice "Dos bordes" o
        "Un borde" segun cuantas lineas tenga, porque es literalmente otro
        modelo el que se ajusta.
        """
        for p in self.pestanas:
            if not isinstance(p, Pestana):
                continue
            try:
                p.al_cambiar_fuente(self.datos.fuente)
                self.nb.tab(p, text=p.titulo_actual())
            except Exception:
                pass

    def reiniciar_pestanas(self, excepto=None):
        """
        Marca todo como pendiente despues de cambiar los archivos.

        No se recalcula solo a proposito: varias pestanas tardan minutos, y
        recalcular las 13 de una seria peor que dejar que el usuario abra las
        que le importan. Lo que si hay que hacer es BORRAR el grafico viejo: si
        no, queda en pantalla un resultado de otros datos.
        """
        for p in self.pestanas:
            if p is excepto or not isinstance(p, Pestana) or p._corriendo:
                continue
            p.calculada = False
            p._ultimo = None
            p.fig.clear()
            ax = p.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Cambiaron los archivos.\nApretá «Recalcular».",
                    ha="center", va="center", color="0.35")
            ax.set_axis_off()
            p.canvas.draw_idle()
        self.estado.set("archivos cambiados: cada pestaña se recalcula al abrirla")

    def configurar_nitidez(self):
        """
        Ajusta Tk y matplotlib al DPI real de la pantalla.

        Con DPI awareness activo, Windows ya no estira la ventana, pero entonces
        todo queda CHICO si no escalamos nosotros. Dos perillas:
          - `tk scaling`: Tk dibuja las fuentes (que estan en puntos) al DPI real.
          - figure.dpi de matplotlib: mas pixeles por pulgada = texto y lineas
            del grafico nitidos. Como el canvas se estira al widget, subir el dpi
            equivale a agrandar el contenido del grafico.
        """
        try:
            dpi = float(self.winfo_fpixels("1i"))       # px por pulgada reales
        except Exception:
            dpi = 96.0
        self.escala = max(1.0, min(dpi / 96.0, 3.0))
        try:
            self.tk.call("tk", "scaling", dpi / 72.0)
        except Exception:
            pass
        self.dpi_fig = 100.0 * self.escala

        matplotlib.rcParams.update({
            "figure.dpi": self.dpi_fig,
            "savefig.dpi": 150,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "lines.antialiased": True,
            "path.simplify": True,          # menos vertices en curvas densas
            "path.simplify_threshold": 0.8,
            "agg.path.chunksize": 10000,
        })
        # la ventana, proporcional al escalado (sin pasarse de la pantalla)
        ancho = min(int(1450 * self.escala), self.winfo_screenwidth() - 40)
        alto = min(int(880 * self.escala), self.winfo_screenheight() - 80)
        self.geometry(f"{ancho}x{alto}")
        self.minsize(1100, 680)

    def cambiar_nitidez(self, _evt=None):
        """Cambia el DPI de las figuras en caliente y repinta lo ya calculado."""
        try:
            f = float(self.var_nitidez.get())
        except ValueError:
            return
        self.escala = f
        self.dpi_fig = 100.0 * f
        matplotlib.rcParams["figure.dpi"] = self.dpi_fig
        for p in self.pestanas:
            fig = getattr(p, "fig", None)
            if fig is not None:
                fig.set_dpi(self.dpi_fig)
            p.repintar()
        self.estado.set(f"nitidez x{f:.2f}  (dpi de figura = {self.dpi_fig:.0f})")

    def al_cambiar(self, _evt):
        p = self.pestanas[self.nb.index(self.nb.select())]
        if not getattr(p, "calculada", True):
            if p.auto:
                p.lanzar()
            else:
                self.estado.set(f"{p.titulo}: apretá 'Recalcular' "
                                f"(este cálculo tarda)")

    def revisar_cola(self):
        try:
            while True:
                pest, res = self.cola.get_nowait()
                if res[0] == "progreso":            # aviso, no resultado final
                    self.estado.set(f"{pest.titulo}: {res[1]}")
                    self.progreso_actualizar(res[2] if len(res) > 2 else None)
                    continue
                pest.recibir(res)
        except queue.Empty:
            pass
        self.after(100, self.revisar_cola)


if __name__ == "__main__":
    activar_dpi_awareness()
    App().mainloop()
