"""
pipeline.py — la receta de calibracion, en UN solo lugar, con progreso.

Por que existe
--------------
La cadena completa (XML -> areas -> espectro -> plantilla -> NCC -> SA ->
calibracion) la necesitan la interfaz grafica y la linea de comandos. Si cada
una tuviera su copia, las dos divergirian: es exactamente lo que paso en el
proyecto anterior, donde la receta del ajuste vivia duplicada en una pestana de
la GUI y hubo que extraerla a una funcion cuando los numeros dejaron de
coincidir. Aca la receta es `calibrar_completo` y las dos interfaces la llaman.

Como se reparte la barra de progreso
------------------------------------
Medido en este equipo: parsear el XML de la muestra (1.76 GB) son ~40 s, el del
fondo (0.35 GB) ~8 s, y TODO el resto de la cadena junto son ~3 s. O sea que el
parseo es el 95% del trabajo. Dos consecuencias:

  * los dos XML se reparten la barra **por bytes**, no mitad y mitad. Con el
    reparto ingenuo la barra se clava en la primera mitad y el tiempo restante
    estimado sale el doble de lo que es.
  * si un archivo ya tiene las areas cacheadas, su etapa cuesta milisegundos y
    recibe peso ~0. Se consulta ANTES de empezar con `area.hay_cache`, que es lo
    unico que permite un ETA honesto desde el primer segundo.

El ETA no aparece hasta tener algo de historia (3% de avance y 2 s), porque
antes es una regla de tres sobre ruido.
"""

import os
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import (rutas, fisica, drs, area, espectro, simulacion,
                      calibracion)


class Cancelado(Exception):
    """El usuario pidio parar. La levantan los callbacks, no el codigo de calculo."""


class SimulacionIncompatible(Exception):
    """El .root no es de la fuente que se midio."""


# Orden y nombre de las etapas. La interfaz dibuja una fila por cada una.
ETAPAS = [
    ("muestra",   "Leer el XML de la muestra"),
    ("fondo",     "Leer el XML del fondo"),
    ("espectro",  "Espectro neto"),
    ("plantilla", "Plantilla Geant4"),
    ("ncc",       "NCC: intervalo y semilla"),
    ("sa",        "Recocido simulado"),
    ("calib",     "Calibracion"),
]
NOMBRES = dict(ETAPAS)

# Lo que cuestan las etapas que NO son parseo, en segundos medidos. Se usan como
# pesos relativos entre ellas; su peso total frente al parseo lo decide
# `repartir`.
COSTO_FIJO = {"espectro": 0.3, "plantilla": 0.6, "ncc": 0.9, "sa": 1.4,
              "calib": 0.2}


# --------------------------------------------------------------------------
# 1. Progreso
# --------------------------------------------------------------------------
class Progreso:
    """
    Lleva la cuenta del avance y se lo cuenta a quien escuche.

    `sink(mensaje)` recibe dicts. La GUI le pasa el `put` de una cola y los
    dibuja en el hilo de Tk; la linea de comandos le pasa una funcion que
    imprime. Nadie mas sabe de interfaces.

    `cancelar` es un objeto con `.is_set()` (un threading.Event). Se consulta en
    cada aviso, que es lo que permite cortar en medio de un parseo de 40 s: el
    bucle de lectura avisa cada 400 eventos, asi que la respuesta es inmediata
    sin tener que meter chequeos dentro del codigo de calculo.
    """

    def __init__(self, sink=None, cancelar=None, pesos=None):
        self.sink = sink or (lambda m: None)
        self.cancelar = cancelar
        self.pesos = dict(pesos or {})
        self.actual = None
        self.hechas = []
        self.t0 = time.time()
        self.t_etapa = None
        self._a, self._b = 0.0, 1.0     # sub-tramo activo dentro de la etapa

    # --- lo que llama el pipeline -----------------------------------------
    def etapa(self, clave, texto=""):
        """Marca el comienzo de una etapa. Idempotente."""
        if self.actual == clave:
            return
        self.terminar_etapa()
        self.actual = clave
        self.t_etapa = time.time()
        self._a, self._b = 0.0, 1.0
        self._avisar("etapa", clave=clave, texto=texto or NOMBRES.get(clave, clave))
        self.revisar()

    def terminar_etapa(self, nota=""):
        if self.actual is None:
            return
        dt = time.time() - self.t_etapa
        self.hechas.append(self.actual)
        self._avisar("fin_etapa", clave=self.actual, segundos=dt, nota=nota)
        self.actual = None

    def fase(self, a, b):
        """
        Mapea el 0-100% del sub-paso que viene al tramo [a, b] de la etapa.

        Sin esto, un sub-paso que reporta su propio 0-100% hace saltar la barra
        de la etapa entera cada vez que empieza uno nuevo.
        """
        self._a, self._b = float(a), float(b)

    def fraccion(self, f, texto=""):
        """Avance dentro de la etapa actual, 0..1. `None` = indeterminado."""
        self.revisar()
        if f is None:
            self._avisar("avance", frac=None, global_=None, texto=texto, eta=None)
            return
        f = self._a + float(np.clip(f, 0.0, 1.0)) * (self._b - self._a)
        g = self.global_(f)
        self._avisar("avance", frac=f, global_=g, texto=texto, eta=self.eta(g))

    def revisar(self):
        if self.cancelar is not None and self.cancelar.is_set():
            raise Cancelado()

    def log(self, texto):
        self._avisar("log", texto=texto)

    # --- cuentas -----------------------------------------------------------
    def global_(self, f_etapa):
        hecho = sum(self.pesos.get(k, 0.0) for k in self.hechas)
        return min(hecho + self.pesos.get(self.actual, 0.0) * f_etapa, 1.0)

    def eta(self, g):
        """Segundos restantes, o None si todavia no hay con que estimarlos."""
        t = time.time() - self.t0
        if g is None or g < 0.03 or t < 2.0:
            return None
        return max(t * (1.0 - g) / g, 0.0)

    def _avisar(self, tipo, **kw):
        kw["tipo"] = tipo
        kw["t"] = time.time() - self.t0
        self.sink(kw)


def repartir(corrida, opciones):
    """
    Pesos de cada etapa, medidos antes de empezar.

    El parseo domina, asi que su peso sale de los BYTES de cada XML; un archivo
    ya cacheado pesa casi nada. Lo que queda se reparte entre las otras cinco
    etapas segun `COSTO_FIJO`.
    """
    c = rutas.resolver(corrida)
    kw = dict(metodo=opciones.get("metodo", "global"),
              limite=opciones.get("limite"),
              plantilla=opciones.get("plantilla"))
    coste = {}
    for clave, ruta in (("muestra", c.muestra), ("fondo", c.fondo)):
        if not ruta or not os.path.exists(ruta):
            coste[clave] = 0.0
            continue
        tibio = area.hay_cache(ruta, **kw)
        # ~44 MB/s medidos; con cache es leer un npz de unos MB
        coste[clave] = 0.2 if tibio else os.path.getsize(ruta) / 44e6
    coste.update(COSTO_FIJO)
    total = sum(coste.values()) or 1.0
    return {k: v / total for k, v in coste.items()}


# --------------------------------------------------------------------------
# 2. La receta
# --------------------------------------------------------------------------
DEFAULTS = dict(metodo="global", rama="Lvis_keVee", bins=500, rango_max=None,
                inicio_rel=0.45, n_arranques=8, limite=None,
                forzar_simulacion=False, guardar=True,
                # el modelo: exceso de bajo canal adentro y ventana de ajuste
                # extendida. Los valores viven en calibracion.EXCESO /
                # INICIO_AJUSTE_REL, aca solo se los deja pasar (None = el
                # default del modulo, para no tener dos verdades).
                exceso=None, tipo_exceso="exp", inicio_ajuste_rel="auto",
                # `post_fijo` [ns] fuerza la MISMA compuerta en varias corridas.
                # Hace falta para compararlas: el `post` que mide cada archivo
                # depende de su propio pulso, y dos corridas con `post` distinto
                # recogen fracciones distintas de la carga, lo que entra derecho
                # en D. Medido en el Cs-137 del detector A, cuyas tres corridas
                # midieron 55, 98 y 124 ns: la ley de ganancia daba b = 10.4 con
                # residuos de +-4.9%; corrigiendo por la fraccion de carga baja
                # a b = 10.0 y +-2.3%. El Co-60, cuyas tres corridas quedaron
                # todas recortadas al mismo limite del registro (75 ns), da
                # b = 9.5 con residuos de +-0.5% sin corregir nada.
                post_fijo=None,
                # Solo con metodo="plantilla": que forma del pulso ajustar
                # ("junio", "ago-sep" o el dict ya cargado). None = la que le
                # toca a la MUESTRA por fecha (area.plantilla_para), y el fondo
                # usa esa misma, igual que usa la compuerta de la muestra.
                plantilla=None)


def calibrar_completo(corrida, opciones=None, progreso=None):
    """
    XML -> ... -> calibracion, avisando el avance en cada etapa.

    Devuelve un dict con el resultado de `calibracion.calibrar` y ademas todo lo
    que los diagnosticos necesitan para dibujarse sin recalcular: la ventana de
    integracion medida, el espectro, el barrido del NCC, las corridas del
    recocido con su semilla y su camino (`salida["sa"]`) y el `Problema` (que es
    lo que permite mapear el paisaje del chi2).
    """
    o = dict(DEFAULTS)
    o.update(opciones or {})
    c = rutas.resolver(corrida)
    pr = progreso or Progreso()
    if not pr.pesos:
        pr.pesos = repartir(c, o)

    if not c.fuente or c.fuente not in fisica.FUENTES:
        raise ValueError("fuente desconocida: %r. Conocidas: %s"
                         % (c.fuente, ", ".join(sorted(fisica.FUENTES))))

    salida = {"corrida": c, "opciones": o}

    # ---- 1. muestra: medir la ventana y leer las areas --------------------
    pr.etapa("muestra")
    pr.fase(0.0, 0.06)
    pr.fraccion(0.0, "midiendo la ventana de integracion")
    canal_m, _ = drs.detectar_canal(c.muestra)
    v = area.calibrar_ventana(c.muestra, canal=canal_m, verbose=False)
    if o["post_fijo"] is not None:
        v = dict(v)
        v["post_medido"] = v["post"]
        v["post"] = float(o["post_fijo"])
        v["fraccion"] = area.fraccion_de_carga(v, v["post"])
        v["frac_eventos_enteros"] = float((v["espacio"] >= v["post"]).mean())
    salida["ventana"] = v
    pr.log("ventana medida: pre = %.1f ns, post = %.1f ns (recoge %.2f%% de la "
           "carga); canal %s" % (v["pre"], v["post"], 100 * v["fraccion"], canal_m))

    pl = None
    if o["metodo"] == "plantilla":
        pl = o["plantilla"]
        if not isinstance(pl, dict):
            pl = area.cargar_plantilla(pl or area.plantilla_para(c.muestra))
        salida["plantilla_pulso"] = pl
        pr.log("area por AJUSTE de la plantilla '%s' (carga hasta %.0f ns del "
               "onset; hash %s)" % (pl.get("nombre", "?"), pl["t_carga"],
                                    pl["hash"]))

    pr.fase(0.06, 1.0)
    # se pregunta ANTES de leer: despues no hay manera de distinguir un parseo
    # rapido de una lectura de cache, y adivinarlo por el tiempo es una regla
    # que falla justo con los archivos chicos.
    tibio_m = area.hay_cache(c.muestra, metodo=o["metodo"], pre=v["pre"],
                             post=v["post"], canal=canal_m, limite=o["limite"],
                             plantilla=pl)
    d_m = area.areas_de_archivo(
        c.muestra, metodo=o["metodo"], pre=v["pre"], post=v["post"],
        canal=canal_m, limite=o["limite"], verbose=False, plantilla=pl,
        callback=lambda f, n: pr.fraccion(f, "%s eventos" % format(n, ",d")))
    salida["areas_muestra"] = d_m
    pr.terminar_etapa("%s pulsos%s" % (format(d_m["area"].size, ",d"),
                                       "  ·cache" if tibio_m else ""))

    # ---- 2. fondo ---------------------------------------------------------
    # La ventana es la MISMA que la de la muestra, a proposito: muestra y fondo
    # pueden venir de sesiones distintas, y si cada uno eligiera la suya
    # estarian midiendo observables diferentes y la resta no significaria nada.
    d_f = None
    pr.etapa("fondo")
    if c.fondo and os.path.exists(c.fondo):
        tibio_f = area.hay_cache(c.fondo, metodo=o["metodo"], pre=v["pre"],
                                 post=v["post"], limite=o["limite"],
                                 plantilla=pl)
        d_f = area.areas_de_archivo(
            c.fondo, metodo=o["metodo"], pre=v["pre"], post=v["post"],
            limite=o["limite"], verbose=False, plantilla=pl,
            callback=lambda f, n: pr.fraccion(f, "%s eventos" % format(n, ",d")))
        pr.terminar_etapa("%s pulsos%s" % (format(d_f["area"].size, ",d"),
                                           "  ·cache" if tibio_f else ""))
    else:
        pr.log("esta corrida no tiene fondo: no se resta nada")
        pr.terminar_etapa("sin fondo")
    salida["areas_fondo"] = d_f

    # ---- 3. espectro neto -------------------------------------------------
    pr.etapa("espectro")
    pr.fraccion(None, "histograma con binning comun y resta de fondo")
    e = espectro.espectro_neto(d_m, d_f, bins=o["bins"],
                               rango_max=o["rango_max"], verbose=False)
    e["corrida"], e["fuente"] = c.clave, c.fuente
    salida["espectro"] = e
    pr.log("espectro: %d bines de %.2f mV*ns; muestra %.0f s vivos de %.0f de "
           "pared (%.1f%%), k_fondo = %.4f"
           % (e["bins"], e["rango_max"] / e["bins"], e["t_vivo_muestra"],
              e["t_pared_muestra"], 100 * e["frac_viva_muestra"], e["k_f"]))
    pr.terminar_etapa()

    # ---- 4. plantilla Geant4 ---------------------------------------------
    pr.etapa("plantilla")
    if not c.sim or not os.path.exists(c.sim):
        raise SimulacionIncompatible(
            "no hay archivo de simulacion. La calibracion por forward folding "
            "necesita el .root de Geant4 de esta fuente.")
    pr.fase(0.0, 0.5)
    pr.fraccion(0.2, "leyendo los gammas primarios del .root")
    lineas = simulacion.lineas_primarias(c.sim)
    salida["lineas_sim"] = lineas
    coincide = simulacion.coincide_fuente(c.sim, c.fuente)
    pr.log("el .root simula %s keV; la fuente elegida es %s (%s)"
           % (", ".join("%.0f" % e0 for e0, _ in lineas), c.fuente,
              "coincide" if coincide else "NO COINCIDE"))
    if not coincide and not o["forzar_simulacion"]:
        raise SimulacionIncompatible(
            "la simulacion no es de %s: sus gammas primarios son %s keV.\n"
            "Calibrar contra una plantilla de otra fuente no da error en "
            "ningun lado, da un numero equivocado bien formateado."
            % (c.fuente, ", ".join("%.0f" % e0 for e0, _ in lineas)))
    pr.fase(0.5, 1.0)
    pr.fraccion(0.5, "cargando la plantilla %s" % o["rama"])
    plantilla = simulacion.cargar(c.sim, o["rama"])
    salida["plantilla"] = plantilla
    pr.terminar_etapa(o["rama"])

    # ---- 5, 6 y 7: NCC -> recocido -> calibracion -------------------------
    # Las tres viven adentro de `calibracion.calibrar`; las etapas se marcan
    # desde sus callbacks, que es lo que evita tener aca una segunda copia de
    # la receta.
    def cb_ncc(f, txt=""):
        pr.etapa("ncc")
        pr.fraccion(f, txt)

    def cb_sa(f, txt=""):
        pr.etapa("sa")
        pr.fraccion(f, txt)

    detalles = {}
    r, p, _, ncc_r = calibracion.calibrar(
        c, rama=o["rama"], metodo=o["metodo"], n_arranques=o["n_arranques"],
        verbose=False, esp=e, plantilla=plantilla,
        inicio_rel=o["inicio_rel"], cb_ncc=cb_ncc, cb_sa=cb_sa,
        exceso=o["exceso"], tipo_exceso=o["tipo_exceso"],
        inicio_ajuste_rel=o["inicio_ajuste_rel"], detalles=detalles)

    pr.etapa("calib")
    pr.fraccion(0.5, "errores por curvatura y bordes")
    salida["resultado"] = r
    salida["problema"] = p
    salida["ncc"] = ncc_r
    # el recocido, con su semilla y su camino: es lo que dibujan las pestanas
    # del SA, y sin esto habria que correrlo otra vez para verlo.
    salida["sa"] = detalles
    if o["guardar"]:
        salida["json"] = calibracion.guardar(r)
        espectro.guardar(e, "espectro_%s_%s" % (c.clave, o["metodo"]))
    pr.log("D = %.4f +- %.4f mV*ns/%s ; bordes en %s mV*ns ; chi2/ndf = %.3f"
           % (r["D"], r["err_D"],
              "keVee" if o["rama"].startswith("Lvis") else "keV",
              ", ".join("%.0f" % b for b in r["bordes_ch"]), r["chi2_ndf"]))
    pr.terminar_etapa()
    pr.fraccion(1.0, "listo")
    return salida


# --------------------------------------------------------------------------
# demo: la misma receta, por linea de comandos
# --------------------------------------------------------------------------
def _barra_consola():
    estado = {"etapa": ""}

    ancho = 96

    def limpiar():
        sys.stdout.write("\r" + " " * ancho + "\r")

    def sink(m):
        # Todo lo transitorio se escribe con \r sobre la misma linea y se limpia
        # antes de imprimir algo permanente. Si no, los restos de la barra
        # quedan pegados al final de la linea siguiente.
        if m["tipo"] == "etapa":
            estado["etapa"] = m["texto"]
        elif m["tipo"] == "fin_etapa":
            limpiar()
            sys.stdout.write("  %-30s %6.1f s   %s\n"
                             % (NOMBRES.get(m["clave"], ""), m["segundos"],
                                m.get("nota", "")))
        elif m["tipo"] == "avance" and m["frac"] is not None:
            n = int(24 * m["frac"])
            eta = ("  ~%3.0f s" % m["eta"]) if m["eta"] else "        "
            sys.stdout.write("\r  %-30s [%s%s] %3.0f%%%s %-22s"
                             % (estado["etapa"], "#" * n, "." * (24 - n),
                                100 * m["frac"], eta, m.get("texto", "")[:22]))
        elif m["tipo"] == "log":
            limpiar()
            sys.stdout.write("      %s\n" % m["texto"])
        sys.stdout.flush()
    return sink


def _demo():
    rutas.salida_utf8()
    # python analisis/pipeline.py [clave] [--area global|convergencia|plantilla]
    args = [a for a in sys.argv[1:]]
    opciones = {}
    if "--area" in args:
        i = args.index("--area")
        opciones["metodo"] = args[i + 1]
        del args[i:i + 2]
    clave = args[0] if args else "co60"
    c = rutas.corrida(clave)
    print("=" * 78)
    print("%s  -  %s%s" % (clave, c.descripcion,
                           ("   [area: %s]" % opciones["metodo"])
                           if opciones else ""))
    print("=" * 78)
    pr = Progreso(sink=_barra_consola())
    pr.pesos = repartir(c, dict(DEFAULTS, **opciones))
    print("  pesos de la barra: %s"
          % ", ".join("%s %.0f%%" % (k, 100 * v)
                      for k, v in sorted(pr.pesos.items(), key=lambda x: -x[1])))
    t0 = time.time()
    s = calibrar_completo(c, opciones, pr)
    print("\n")
    calibracion.informe(s["resultado"])
    print("  (%.1f s en total)" % (time.time() - t0))


if __name__ == "__main__":
    _demo()
