"""
calibrar_todo.py — calibra TODAS las corridas del catalogo y deja las figuras
ordenadas por detector, voltaje y muestra, mas los dos resumenes.

Que hace
--------
Corre `pipeline.calibrar_completo` sobre cada corrida disponible y guarda, para
cada una, las mismas paginas que dibuja la interfaz, en
`figuras/calibraciones/Detector X/<HV>/<fuente>/`. Al final arma los dos
resumenes que responden las dos preguntas distintas que se pueden hacer con
estos datos (CLAUDE.md 7.46):

  * **TEST B, la ganancia**: D contra alto voltaje, una curva por detector y
    fuente. Lo que se espera de un PMT no es una recta sino `D ~ V^b` con b
    entre 6 y 10, asi que se ajustan las dos y se comparan.
  * **TEST A, la linealidad del eje**: canal contra energia con DOS fuentes
    medidas en el MISMO punto de operacion. Es el unico test que valida la
    calibracion; con una sola fuente el mapa por el origen tiene un solo
    parametro y cualquier punto cae exacto sobre su propia recta.

Cuatro decisiones que valen para TODAS las corridas, y por que
-------------------------------------------------------------
  * **Escala `Edep_keV`** (pedido explicito): la plantilla es la energia
    DEPOSITADA, no la luz. La D que sale convierte area a keV de energia
    depositada, con el quenching de Birks adentro de D. No es la misma cantidad
    que la D en keVee de `calibracion.py`, y las dos NO se mezclan: difieren un
    ~4.5% en el Co-60 y ~7% en el Cs-137, que es justamente el quenching.
  * **Compuerta comun corta** (`POST_FIJO`, en ns despues del onset): el `post`
    que mide cada archivo depende de su propio pulso, y dos corridas con
    compuertas distintas recogen fracciones distintas de la carga, lo que entra
    derecho en D (CLAUDE.md 7.53). Para COMPARAR corridas -que es todo lo que
    hace este script- la compuerta tiene que ser la misma, y es la misma que
    usa la interfaz (`gui.COMPUERTA_DEFECTO`).
  * **Modelo sin el exceso de bajo canal** (`EXCESO`): la plantilla
    convolucionada y nada mas. Es la misma decision que toma la interfaz
    (`gui.EXCESO`), y tiene su contrapartida: sin esa componente la ventana de
    ajuste no puede entrar en la zona sucia, asi que `inicio_rel` pasa a ser un
    sistematico que se elige a mano.
  * **Sin fondo donde no lo hay**: solo dos corridas tienen fondo propio (el
    Co-60 de junio y el Cs-137 del detector A a 1600 V). Prestarle a las demas
    el fondo de otra sesion seria peor que no restar nada: el fondo depende del
    punto de operacion (a 1700 V el mismo cosmico deja el triple de area), asi
    que meteria una forma ajena en el espectro. Las figuras y la tabla dicen
    cuales llevan fondo restado y cuales no.

    python analisis/calibrar_todo.py                # todo lo que falte
    python analisis/calibrar_todo.py --rehacer      # ignora lo ya calculado
    python analisis/calibrar_todo.py --solo a_co60_1600 cs137_1600
    python analisis/calibrar_todo.py --resumen      # solo rehace los resumenes
    python analisis/calibrar_todo.py --area plantilla   # otro metodo de area,
                                     # en figuras/calibraciones_plantilla/
"""

import json
import os
import sys
import time
import traceback

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.ticker import ScalarFormatter

from analisis import (rutas, pipeline, registro, voltaje, calibracion, fisica,
                     area, drs)
from analisis import gui as g

# --------------------------------------------------------------------------
# Configuracion: lo que es igual para todas las corridas
# --------------------------------------------------------------------------
RAMA = "Edep_keV"            # energia depositada, no luz

# Compuerta y modelo: los mismos que la interfaz, importados de ahi.
#
# No se copian los valores a proposito. Son la definicion del observable, y con
# dos copias basta una tarde para que la interfaz y el lote midan cosas
# distintas sin que nada falle: los numeros seguirian saliendo, solo que no
# serian comparables. Si hay que cambiarlos, se cambian en `gui.py`, que es
# donde estan documentados.
#
# Lo que hay que recordar de la compuerta corta: recoge una fraccion de la
# carga del pulso (la columna «carga» del informe dice cuanta), asi que estas D
# NO son la ganancia de pulso entero. Son un observable comun, que es lo que
# hace falta para COMPARAR corridas (CLAUDE.md 7.52 y 7.53).
POST_FIJO = g.COMPUERTA_DEFECTO
EXCESO = g.EXCESO
METODO = "global"
OPCIONES = dict(rama=RAMA, post_fijo=POST_FIJO, metodo=METODO, exceso=EXCESO,
                guardar=True)

SALIDA = os.path.join(rutas.FIGURAS, "calibraciones")
RESULTADOS = os.path.join(SALIDA, "resultados.json")


def configurar_area(metodo):
    """
    `--area plantilla` (o `convergencia`): otro metodo de area para el lote
    ENTERO, sin tocar el default.

    Todo lo del metodo no global va a OTRA carpeta
    (`figuras/calibraciones_<metodo>/`) con su propio `resultados.json`. Si
    compartieran archivo, un lote con plantilla pisaria la tabla que leen
    `informe.py` y el visor, que es la del metodo vigente, y nada avisaria. La
    compuerta comun se sigue pasando: con "plantilla" define la BANDERA de cada
    evento (los mismos descartes que el global) y deja calculada al lado el
    area por compuerta, para compararlas evento a evento.
    """
    global METODO, OPCIONES, SALIDA, RESULTADOS
    METODO = metodo
    OPCIONES = dict(OPCIONES, metodo=metodo)
    if metodo != "global":
        SALIDA = os.path.join(rutas.FIGURAS, "calibraciones_%s" % metodo)
        RESULTADOS = os.path.join(SALIDA, "resultados.json")

# Las paginas de la interfaz que se guardan por corrida. Son las mismas
# funciones que dibuja `gui.py`, llamadas sin Tk (ver `_Stub`), asi que lo que
# queda en el PNG es exactamente lo que muestra la interfaz.
PAGINAS = [("1_calibracion", g.App._fig_resultado, ()),
           ("2_parametros", g.App._fig_parametros, ()),
           ("3_espectro_y_simulacion", g.App._fig_espectro, ()),
           ("4_ncc_e_intervalo", g.App._fig_ncc, ()),
           ("5_recocido", g.App._fig_recocido, ("sa",)),
           ("6_ventana_de_integracion", g.App._fig_ventana, ("ventana",))]


class _Stub:
    """Las figuras solo usan atributos de CLASE, asi que no hace falta una App
    de verdad (que pediria un root de Tk). Mismo truco que
    `verificaciones/_verif_gui_pestanas.py`."""
    _NOMBRE_LIBRE = g.App._NOMBRE_LIBRE
    _CORTO = g.App._CORTO


# --------------------------------------------------------------------------
# 1. Que corridas, y donde va cada una
# --------------------------------------------------------------------------
def descartada(c):
    """
    True si la corrida esta DESCARTADA como medicion valida.

    El unico caso hoy: las tres corridas de Cs-137 del detector D del 1-sep. En
    esa sesion el detector tenia la tierra mal conectada; las buenas son las
    del 8-sep. Se calibran igual y se guardan, porque miden el tamano de esa
    falla, pero no entran en ningun ajuste: ni en la ley de ganancia ni en el
    test de linealidad.

    Ojo con la tentacion de leer su comparacion con las buenas como una medida
    de reproducibilidad del montaje. No lo es: una de las dos mediciones esta
    mal por una causa identificada, asi que lo que se compara es "con falla" y
    "sin falla", no dos intentos equivalentes.
    """
    return rutas.descartada(c.clave)


def error_usable(f):
    """
    La barra de error de D que se puede usar para pesar un ajuste.

    El detector D da a veces una matriz de covarianza no definida positiva y
    `err_D` sale nan (ya estaba visto: CLAUDE.md 7.57). Un nan en los pesos
    contamina el ajuste entero y lo deja en nan sin decir por que, asi que se
    reemplaza por un piso del 1% de D -del orden de lo que dan las corridas
    sanas- y la fila queda marcada. Es un parche para poder seguir, no una
    medicion: lo que hay que mirar en esas corridas es el chi2, no el error.
    """
    e = f.get("err_D")
    if e is not None and np.isfinite(e) and e > 0:
        return float(e), False
    return 0.01 * float(f["D"]), True


def corridas():
    """Las corridas disponibles, ordenadas por detector, voltaje y fuente."""
    cs = [c for c in rutas.CORRIDAS.values() if c.disponible()]
    return sorted(cs, key=lambda c: (c.detector or "", c.hv or 0, c.fuente,
                                     descartada(c)))


def etiqueta(c):
    """Nombre corto y unico: 'A_1600V_Co-60'. Va en cada archivo, para que un
    PNG sacado de su carpeta siga diciendo de que es."""
    return "%s_%gV_%s%s" % (c.detector, c.hv, c.fuente,
                            "_descartada" if descartada(c) else "")


def carpeta(c):
    """`figuras/calibraciones/Detector A/1600 V/Co-60/`, creada si falta."""
    sub = c.fuente + (" (descartada)" if descartada(c) else "")
    d = os.path.join(SALIDA, "Detector %s" % c.detector, "%g V" % c.hv, sub)
    os.makedirs(d, exist_ok=True)
    return d


# --------------------------------------------------------------------------
# 2. Una corrida
# --------------------------------------------------------------------------
def _sink(prefijo="      "):
    """Imprime solo lo permanente. La barra de la consola de `pipeline` escribe
    con \\r, que en un archivo de log deja una sola linea ilegible de kilometros."""
    def sink(m):
        if m["tipo"] == "log":
            print("%s%s" % (prefijo, m["texto"]))
        elif m["tipo"] == "fin_etapa":
            print("%s%-26s %6.1f s   %s"
                  % (prefijo, pipeline.NOMBRES.get(m["clave"], m["clave"]),
                     m["segundos"], m.get("nota", "")))
        sys.stdout.flush()
    return sink


def figuras_de(salida, c):
    """Guarda las paginas de la interfaz de esta corrida. Devuelve los caminos."""
    stub = _Stub()
    destino, hechas = carpeta(c), []
    for nombre, metodo, necesita in PAGINAS:
        if any(k not in salida for k in necesita):
            print("      %-26s salteada (falta %s)"
                  % (nombre, ", ".join(k for k in necesita if k not in salida)))
            continue
        fig = Figure(figsize=(15, 8.6), dpi=100)
        FigureCanvasAgg(fig)
        try:
            metodo(stub, fig, salida)
        except Exception as e:                                  # noqa: BLE001
            print("      %-26s FALLA: %s: %s" % (nombre, type(e).__name__, e))
            continue
        ruta = os.path.join(destino, "%s_%s.png" % (etiqueta(c), nombre))
        fig.savefig(ruta)
        hechas.append(rutas.relativa(ruta))
    return hechas


def calibrar_una(c):
    """La cadena completa + las figuras. Devuelve la fila de resultados."""
    pr = pipeline.Progreso(sink=_sink())
    pr.pesos = pipeline.repartir(c, OPCIONES)
    t0 = time.time()
    salida = pipeline.calibrar_completo(c, OPCIONES, pr)
    r, v, e = salida["resultado"], salida["ventana"], salida["espectro"]
    print("      D = %.4f +- %.4f mV*ns/keV   chi2/ndf = %.3f   FWHM = %.1f%%"
          % (r["D"], r["err_D"], r["chi2_ndf"], 100 * r["fwhm_rel"]))

    fila = dict(
        clave=c.clave, detector=c.detector, hv=c.hv, fuente=c.fuente,
        descartada=descartada(c), descripcion=c.descripcion,
        archivo=rutas.relativa(c.muestra), fecha=rutas.fecha(c.muestra),
        fondo=rutas.relativa(c.fondo) if c.fondo else "",
        hay_fondo=bool(salida.get("areas_fondo") is not None),
        sim=rutas.relativa(c.sim), rama=RAMA, metodo=METODO,
        # el resultado
        D=r["D"], err_D=r["err_D"], E_ref=r["E_ref"], punto=r["punto"],
        inv_D=r["inv_D"], err_inv_D=r["err_inv_D"],
        sigma_rel=r["sigma_rel"], fwhm_rel=r["fwhm_rel"], sigma_E=r["sigma_E"],
        chi2_ndf=r["chi2_ndf"], ndf=r["ndf"],
        corr_D_sigma=r["corr_D_sigma"],
        dispersion_multiarranque=r["dispersion_multiarranque"],
        bordes_E=list(r["bordes_E"]), bordes_ch=list(r["bordes_ch"]),
        ventana=list(r["ventana"]), rho=r["rho"],
        lam_exceso=r.get("lam_exceso"),
        # la cadena de medicion, sin la cual D no significa nada.
        # Dos partes: la de SENAL (como llego el pulso al digitalizador) y la
        # de INTEGRACION (que pedazo del pulso se midio). `atenuador` importa
        # aparte de las demas: 20 dB dividen D por 10, asi que sin ese dato una
        # D no se puede comparar con la de otra sesion (CLAUDE.md 7.70).
        atenuador=r.get("atenuador"), cadena=r.get("cadena"),
        umbral=c.umbral, muestreo=c.muestreo,
        post=v["post"], post_medido=v.get("post_medido"),
        fraccion_carga=v["fraccion"],
        frac_eventos_enteros=v.get("frac_eventos_enteros"),
        pulsos=int(e["n_muestra"]) if "n_muestra" in e else None,
        pulsos_leidos=int(salida["areas_muestra"]["area"].size),
        t_vivo=e.get("t_vivo_muestra"), t_pared=e.get("t_pared_muestra"),
        frac_viva=e.get("frac_viva_muestra"), k_f=e.get("k_f"),
        segundos=time.time() - t0,
        figuras=figuras_de(salida, c),
    )
    d_m = salida["areas_muestra"]
    if "area_compuerta" in d_m:
        # metodo="plantilla": que forma se ajusto y como salio el ajuste
        ok = area.buenos(d_m)
        fila.update(
            plantilla=str(d_m.get("plantilla", "")),
            plantilla_hash=str(d_m.get("plantilla_hash", "")),
            frac_extrapolada=float(np.mean(d_m["extrapolada"][ok])),
            chi2_red_mediana=float(np.median(d_m["chi2_red"][ok])),
            cociente_plantilla_compuerta=float(np.median(
                d_m["area"][ok & (d_m["amplitud"] > 20)]
                / d_m["area_compuerta"][ok & (d_m["amplitud"] > 20)])))
    return fila


# --------------------------------------------------------------------------
# 3. Guardar y recuperar lo ya calculado
# --------------------------------------------------------------------------
def cargar():
    try:
        with open(RESULTADOS, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def escribir(res):
    os.makedirs(SALIDA, exist_ok=True)
    with open(RESULTADOS, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False, default=float)


def csv(res):
    """Tabla plana, para abrirla en cualquier lado."""
    cols = ["detector", "hv", "fuente", "descartada", "D", "err_D", "inv_D",
            "E_ref", "chi2_ndf", "ndf", "fwhm_rel", "sigma_E", "rho",
            "hay_fondo", "post", "fraccion_carga", "pulsos", "t_vivo",
            "clave", "archivo"]
    ruta = os.path.join(SALIDA, "resultados.csv")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for fila in ordenadas(res):
            f.write(",".join(_celda(fila.get(k)) for k in cols) + "\n")
    return ruta


def _celda(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "si" if v else "no"
    if isinstance(v, float):
        return "%.6g" % v
    s = str(v)
    return '"%s"' % s if "," in s else s


def ordenadas(res):
    """Las filas de calibracion, en orden. Ignora las entradas que no son una
    corrida (la clave empieza con '_'), como el control de compuerta."""
    filas = [f for k, f in res.items() if not k.startswith("_")]
    return sorted(filas, key=lambda f: (f["detector"], f["hv"], f["fuente"],
                                        f["descartada"]))


# --------------------------------------------------------------------------
# 4. Los dos tests
# --------------------------------------------------------------------------
COLOR = {"A": "#1f4e79", "B": "#c4630a", "C": "#1f7a3d", "D": "#8a2be2"}
MARCA = {"Co-60": "o", "Cs-137": "s"}


def test_b(res):
    """
    D contra alto voltaje, una entrada por (detector, fuente).

    Se excluyen las corridas descartadas (`descartada`): entrarian como un
    segundo valor en el mismo voltaje y el ajuste las promediaria con las
    buenas, que es exactamente lo que no hay que hacer con una medicion que se
    sabe mala. Su comparacion con las buenas va aparte (`pares_descartadas`).
    """
    grupos = {}
    for f in ordenadas(res):
        if f["descartada"]:
            continue
        grupos.setdefault((f["detector"], f["fuente"]), []).append(f)
    salida = {}
    for k, filas in sorted(grupos.items()):
        if len(filas) < 2:
            salida[k] = dict(filas=filas, ajuste=None, avisos=[],
                             fechas=sorted({f["fecha"] for f in filas}),
                             motivo="un solo voltaje")
            continue
        V = np.array([f["hv"] for f in filas], float)
        D = np.array([f["D"] for f in filas], float)
        e = np.array([error_usable(f)[0] for f in filas], float)
        salida[k] = dict(filas=filas, ajuste=voltaje.ajustar(V, D, e),
                         motivo="",
                         fechas=sorted({f["fecha"] for f in filas}),
                         avisos=rutas.misma_cadena(
                             [rutas.corrida(f["clave"]) for f in filas]))
    return salida


_QUENCH = {}


def quenching(fuente):
    """
    <Lvis/Edep> de la simulacion de esa fuente, en la zona del borde.

    Hace falta para leer el test A en escala Edep, y no es un detalle: el area
    medida es proporcional a la LUZ, y la luz por keV depositado no es la misma
    a 477 que a 963 keV porque el quenching de Birks todavia no saturo. O sea
    que en escala Edep **dos fuentes tienen que dar D distintas aunque la
    electronica sea perfectamente lineal** (CLAUDE.md 7.46). El cociente que
    corresponde esperar es q(Cs)/q(Co) = 0.980, no 1.
    """
    if fuente not in _QUENCH:
        sim = (rutas.SIM_CO60 if fuente == "Co-60" else rutas.SIM_CS137)
        _QUENCH[fuente] = calibracion.quenching_simulacion(
            sim, fisica.e_ref(fuente))
    return _QUENCH[fuente]


def cociente_esperado():
    """Lo que da D(Cs-137)/D(Co-60) con una cadena perfectamente lineal."""
    return quenching("Cs-137") / quenching("Co-60")


def test_a(res):
    """
    Canal contra energia en cada punto de operacion, con las fuentes que haya.

    Devuelve tambien los puntos que NO se pueden testear y por que: un detector
    con una sola fuente no falla el test, no lo puede correr, y eso hay que
    decirlo con todas las letras en vez de dejar el hueco.
    """
    grupos = {}
    for f in ordenadas(res):
        if f["descartada"]:
            continue
        grupos.setdefault((f["detector"], f["hv"]), []).append(f)
    salida = {}
    for k, filas in sorted(grupos.items()):
        # `linealidad` pesa con err_D; si viene nan hay que darle el piso
        entradas = [(f["clave"], dict(f, err_D=error_usable(f)[0]))
                    for f in filas]
        lin = registro.linealidad(entradas)
        razon = None
        if lin is None:
            razon = ("una sola fuente (%s): con el mapa por el origen el punto "
                     "cae exacto sobre su propia recta"
                     % ", ".join(sorted({f["fuente"] for f in filas})))
        # el cociente entre fuentes es el numero corto que resume el test
        porf = {f["fuente"]: f for f in filas}
        cociente = None
        if "Cs-137" in porf and "Co-60" in porf:
            a, b = porf["Cs-137"], porf["Co-60"]
            ea, eb = error_usable(a)[0], error_usable(b)[0]
            q = cociente_esperado()
            cociente = dict(
                valor=a["D"] / b["D"],
                err=(a["D"] / b["D"]) * np.hypot(ea / a["D"], eb / b["D"]),
                esperado=q, corregido=(a["D"] / b["D"]) / q)
        fechas = sorted({f["fecha"] for f in filas})
        salida[k] = dict(filas=filas, lin=lin, razon=razon, cociente=cociente,
                         fechas=fechas, misma_sesion=len(fechas) == 1,
                         avisos=rutas.misma_cadena(
                             [rutas.corrida(f["clave"]) for f in filas]))
    return salida


def pares_descartadas(res):
    """Cada corrida descartada junto a la buena del mismo punto y fuente."""
    pares = []
    for f in ordenadas(res):
        if not f["descartada"]:
            continue
        orig = [o for o in ordenadas(res)
                if not o["descartada"] and o["detector"] == f["detector"]
                and o["hv"] == f["hv"] and o["fuente"] == f["fuente"]]
        if orig:
            pares.append((orig[0], f))
    return pares


def control_compuerta(res, post=75.0):
    """
    El mismo test A con OTRA compuerta, donde el registro lo permite.

    Es el control de la unica decision de analisis que este script impone a
    todas las corridas por igual. Si el desacuerdo entre fuentes fuera un
    artefacto de integrar solo 36 ns -por ejemplo porque el pulso de una fuente
    es mas largo que el de la otra-, con 75 ns tendria que moverse. Solo se
    puede correr donde las DOS corridas del punto tienen registro de sobra: la
    sesion del 08-sep (detector D) no llega, y eso tambien es un resultado.
    """
    salida = {}
    for (det, hv), d in sorted(test_a(res).items()):
        if d["lin"] is None or not d["cociente"]:
            continue
        # el registro tiene que dar: `post_medido` es lo que la corrida podia
        # pagar antes de que se le impusiera la compuerta comun. Sin este
        # filtro, una corrida corta devuelve un numero -no falla- calculado
        # sobre el puñado de eventos que no salieron truncados, que son los de
        # onset temprano, o sea los de mayor amplitud: sesgo puro.
        corto = []
        for f in d["filas"]:
            c = rutas.corrida(f["clave"])
            canal, _ = drs.detectar_canal(c.muestra)
            v = area.calibrar_ventana(c.muestra, canal=canal, verbose=False)
            entero = float((np.asarray(v["espacio"], float) >= post).mean())
            if entero < 0.90:
                corto.append((f["fuente"], entero))
        if corto:
            print("  %s @ %g V: el registro no da para %g ns (%s)"
                  % (det, hv, post,
                     ", ".join("%s: solo el %.0f%% de los eventos entra"
                               % (n, 100 * e) for n, e in corto)))
            continue
        try:
            filas = {}
            for f in d["filas"]:
                c = rutas.corrida(f["clave"])
                o = dict(OPCIONES, post_fijo=post, guardar=False)
                r = pipeline.calibrar_completo(c, o, pipeline.Progreso())["resultado"]
                filas[f["fuente"]] = r["D"]
            q = filas["Cs-137"] / filas["Co-60"]
            salida["%s @ %g V" % (det, hv)] = dict(
                D_Cs=filas["Cs-137"], D_Co=filas["Co-60"], cociente=q,
                cociente_36=d["cociente"]["valor"],
                cambio=100 * (q / d["cociente"]["valor"] - 1))
            print("  %s @ %g V con compuerta de %g ns: Cs/Co = %.4f "
                  "(con %g ns daba %.4f, %+.2f%%)"
                  % (det, hv, post, q, POST_FIJO, d["cociente"]["valor"],
                     100 * (q / d["cociente"]["valor"] - 1)))
        except Exception as e:                                  # noqa: BLE001
            print("  %s @ %g V: no se pudo con %g ns (%s: %s)"
                  % (det, hv, post, type(e).__name__, e))
    return salida


# --------------------------------------------------------------------------
# 5. Las figuras de resumen
# --------------------------------------------------------------------------
def _leyenda_unica(ax, **kw):
    """Leyenda sin etiquetas repetidas (las descartadas comparten etiqueta)."""
    man, etq = ax.get_legend_handles_labels()
    vistos, m, e = set(), [], []
    for h, t in zip(man, etq):
        if t in vistos:
            continue
        vistos.add(t)
        m.append(h)
        e.append(t)
    ax.legend(m, e, **kw)


def _tabla(ax, filas, columnas, anchos, titulo=None, nota=None):
    ax.set_axis_off()
    if titulo:
        ax.set_title(titulo, fontsize=10, loc="left")
    alto = 0.90 if nota else 1.0
    t = ax.table(cellText=filas, colLabels=columnas, colWidths=list(anchos),
                 cellLoc="left", loc="upper left",
                 bbox=[0, 1.0 - alto, 1, alto])
    t.auto_set_font_size(False)
    t.set_fontsize(7.2)
    for (i, _), celda in t.get_celld().items():
        celda.set_linewidth(0.3)
        celda.set_edgecolor("#d5dce3")
        if i == 0:
            celda.set_facecolor("#1f4e79")
            celda.set_text_props(color="white", fontweight="bold")
        elif i % 2 == 0:
            celda.set_facecolor("#f3f6f9")
    if nota:
        ax.text(0, -0.01, nota, fontsize=7, color="#4a5560",
                va="top", transform=ax.transAxes, wrap=True)


def fig_ganancia(res, tb):
    fig = Figure(figsize=(15.5, 9.5), dpi=110)
    FigureCanvasAgg(fig)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.95))

    # ---- (a) D(V) en log-log, con la ley de potencia ajustada ------------
    ax = fig.add_subplot(gs[0, 0])
    for (det, fuente), d in sorted(tb.items()):
        filas, aj = d["filas"], d["ajuste"]
        V = np.array([f["hv"] for f in filas], float)
        D = np.array([f["D"] for f in filas], float)
        e = np.array([f["err_D"] for f in filas], float)
        col, m = COLOR.get(det, "0.4"), MARCA.get(fuente, "^")
        ax.errorbar(V, D, yerr=e, fmt=m, color=col, ms=6, capsize=3,
                    ls="none",
                    label="%s · %s%s" % (det, fuente,
                                         "" if aj is None
                                         else "  b = %.2f ± %.2f" % (aj["b"],
                                                                    aj["err_b"])))
        if aj is not None:
            vv = np.linspace(V.min() * 0.99, V.max() * 1.01, 60)
            ax.plot(vv, voltaje.evaluar(aj, vv), lw=1.1, color=col, alpha=.75)
    # las descartadas no entran en ningun ajuste, pero verlas es la mitad del
    # diagnostico: muestran cuanto se corrio la ganancia por la falla de tierra
    # del 1-sep, que es la escala de lo que un problema de montaje puede mover
    for f in ordenadas(res):
        if not f["descartada"]:
            continue
        ax.plot([f["hv"]], [f["D"]], "x", ms=9, mew=2,
                color=COLOR.get(f["detector"], "0.4"),
                label="%s · %s (%s, descartada: tierra)"
                      % (f["detector"], f["fuente"], f["fecha"]))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([1550, 1600, 1650, 1700])
    ax.set_xticks([], minor=True)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())
    ax.set_yticks([0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
    ax.set_title("ganancia contra alto voltaje (las lineas son D ∝ V^b)",
                 fontsize=10)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("D [mV·ns / keV depositado]")
    _leyenda_unica(ax, fontsize=7.5, ncol=2)
    ax.grid(alpha=.3, which="both")

    # ---- (b) residuos de la ley de potencia ------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="k", lw=.8)
    for (det, fuente), d in sorted(tb.items()):
        aj = d["ajuste"]
        if aj is None:
            continue
        r = 100 * aj["residuos_potencia"] / aj["D"]
        ax.plot(aj["V"], r, MARCA.get(fuente, "^") + "-", color=COLOR.get(det),
                ms=6, lw=1, label="%s · %s" % (det, fuente))
    ax.set_title("residuos de la ley de potencia", fontsize=10)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("(D − ajuste) / D  [%]")
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(alpha=.3)

    # ---- (c) el exponente b, con la banda tipica de un PMT ---------------
    ax = fig.add_subplot(gs[1, 0])
    ax.axhspan(voltaje.B_PMT[0], voltaje.B_PMT[1], color="#1f7a3d", alpha=.10,
               label="rango tipico de un PMT (%g a %g)" % voltaje.B_PMT)
    k = 0
    etiq = []
    for (det, fuente), d in sorted(tb.items()):
        aj = d["ajuste"]
        if aj is None:
            continue
        ax.errorbar([k], [aj["b"]], yerr=[aj["err_b"]], fmt=MARCA.get(fuente, "^"),
                    color=COLOR.get(det), ms=8, capsize=4)
        etiq.append("%s\n%s" % (det, fuente))
        k += 1
    ax.set_xticks(range(len(etiq)))
    ax.set_xticklabels(etiq, fontsize=8)
    ax.set_xlim(-0.6, len(etiq) - 0.4)
    ax.set_title("el exponente de la ley de ganancia", fontsize=10)
    ax.set_ylabel("b")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3, axis="y")

    # ---- (d) la tabla ----------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    filas = []
    for (det, fuente), d in sorted(tb.items()):
        aj = d["ajuste"]
        if aj is None:
            filas.append([det, fuente, "%d" % len(d["filas"]), "-", "-", "-",
                          d["motivo"]])
            continue
        veredicto = ("potencia" if aj["chi2_potencia"] < aj["chi2_recta"]
                     else "recta")
        if not aj["decide"]:
            veredicto += " (no decide)"
        filas.append([det, fuente, "%d" % aj["n"],
                      "%.2f ± %.2f" % (aj["b"], aj["err_b"]),
                      "%.0f" % aj["chi2_potencia"], "%.0f" % aj["chi2_recta"],
                      veredicto])
    _tabla(ax, filas,
           ["det", "fuente", "n V", "b", "χ² pot.", "χ² recta", "mejor"],
           (0.07, 0.13, 0.07, 0.20, 0.13, 0.13, 0.27),
           titulo="TEST B · la ganancia contra el voltaje",
           nota="'no decide' = el rango de voltajes medido no separa las dos "
                "hipotesis mas de lo que valen las barras de error.")
    fig.suptitle("Resumen · ganancia de los cuatro detectores  (escala %s, "
                 "compuerta comun de %g ns)" % (RAMA, POST_FIJO),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    ruta = os.path.join(SALIDA, "RESUMEN_ganancia.png")
    fig.savefig(ruta)
    return ruta


def fig_linealidad(res, ta):
    fig = Figure(figsize=(15.5, 9.5), dpi=110)
    FigureCanvasAgg(fig)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.95))

    # ---- (a) canal contra energia, un ajuste por punto de operacion ------
    ax = fig.add_subplot(gs[0, 0])
    for (det, hv), d in sorted(ta.items()):
        lin = d["lin"]
        col = COLOR.get(det, "0.4")
        E = np.array([f["E_ref"] for f in d["filas"]], float)
        ch = np.array([f["D"] * f["E_ref"] for f in d["filas"]], float)
        ec = np.array([f["err_D"] * f["E_ref"] for f in d["filas"]], float)
        o = np.argsort(E)
        ax.errorbar(E[o], ch[o], yerr=ec[o], fmt="o", ms=5, color=col,
                    capsize=3, ls="none", alpha=.9)
        if lin is not None:
            xx = np.array([0, 1.15 * E.max()])
            ax.plot(xx, lin["D_origen"] * xx, lw=1.0, color=col, alpha=.7,
                    label="%s @ %g V   D = %.4f" % (det, hv, lin["D_origen"]))
    ax.set_title("TEST A · canal contra energia, con las fuentes de cada punto",
                 fontsize=10)
    ax.set_xlabel("energia del borde Compton [keV depositados]")
    ax.set_ylabel("area del borde [mV·ns]")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=.3)

    # ---- (b) residuos de la recta por el origen --------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="k", lw=.8)
    etiquetas = {}
    for (det, hv), d in sorted(ta.items()):
        lin = d["lin"]
        if lin is None:
            continue
        r = 100 * lin["residuos_origen"] / lin["canal"]
        ax.plot(lin["E"], r, "o-", ms=6, lw=1, color=COLOR.get(det),
                label="%s @ %g V" % (det, hv))
        etiquetas = {float(E): f for E, f in zip(lin["E"], lin["fuentes"])}
    ax.set_title("cuanto se aparta cada fuente de la recta por el origen",
                 fontsize=10)
    ax.set_xlabel("borde Compton de cada fuente")
    ax.set_ylabel("(canal − D·E) / canal  [%]")
    if etiquetas:
        ax.set_xticks(sorted(etiquetas))
        ax.set_xticklabels(["%s\n%.0f keV" % (etiquetas[E], E)
                            for E in sorted(etiquetas)], fontsize=8)
    ax.legend(fontsize=7.5, ncol=3)
    ax.grid(alpha=.3)

    # ---- (c) el cociente entre fuentes, que es el test en un numero ------
    ax = fig.add_subplot(gs[1, 0])
    q = cociente_esperado()
    ax.axhline(q, color="k", ls="--", lw=1.2,
               label="cadena lineal: %.3f  (el quenching de Birks, no 1)" % q)
    ax.axhline(1.0, color="0.7", ls=":", lw=1, label="cociente 1")
    for (det, hv), d in sorted(ta.items()):
        q = d["cociente"]
        if q is None:
            continue
        ax.errorbar([hv], [q["valor"]], yerr=[q["err"]], fmt="o", ms=7,
                    color=COLOR.get(det), capsize=4)
    for det in sorted({k[0] for k in ta}):
        pares = [(hv, d["cociente"]["valor"]) for (dd, hv), d in sorted(ta.items())
                 if dd == det and d["cociente"]]
        if len(pares) > 1:
            ax.plot([p[0] for p in pares], [p[1] for p in pares], "-",
                    lw=1, color=COLOR.get(det), label="detector %s" % det)
    ax.set_title("D(Cs-137) / D(Co-60) en el mismo punto de operacion",
                 fontsize=10)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("cociente de las dos calibraciones")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3)

    # ---- (d) la tabla ----------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    filas = []
    for (det, hv), d in sorted(ta.items()):
        lin, q = d["lin"], d["cociente"]
        if lin is None:
            filas.append(["%s @ %g V" % (det, hv), "-", "-", "-", "-",
                          "1 sola fuente"])
            continue
        peor = 100 * np.max(np.abs(lin["residuos_origen"] / lin["canal"]))
        c0 = "%+.0f ± %.0f" % (lin["ordenada"], lin["err_ordenada"])
        veredicto = ("compatible con el origen"
                     if lin["compatible_con_origen"] else "NO pasa por el origen")
        filas.append(["%s @ %g V" % (det, hv),
                      "%.4f" % lin["D_origen"], "±%.1f %%" % peor, c0,
                      ("misma" if d["misma_sesion"] else "distintas"),
                      "-" if q is None else "%.3f (%+.0f %%)"
                      % (q["valor"], 100 * (q["corregido"] - 1))])
    _tabla(ax, filas,
           ["punto", "D por el origen", "residuo", "ordenada [mV·ns]",
            "sesiones", "Cs/Co (desacuerdo)"],
           (0.17, 0.19, 0.13, 0.20, 0.16, 0.15),
           titulo="TEST A · la linealidad del eje, detector por detector",
           nota="La ordenada sale de un ajuste de 2 puntos y 2 parametros, o "
                "sea que pasa exacto: lo unico que dice algo es su error.\n"
                "Un pedestal del digitalizador seria carga fija por evento, "
                "asi que tendria que dar lo MISMO a los tres voltajes.")
    fig.suptitle("Resumen · linealidad de la calibracion  (escala %s, "
                 "compuerta comun de %g ns)" % (RAMA, POST_FIJO),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    ruta = os.path.join(SALIDA, "RESUMEN_linealidad.png")
    fig.savefig(ruta)
    return ruta


# --------------------------------------------------------------------------
# 6. El informe de texto
# --------------------------------------------------------------------------
def informe(res, tb, ta, destino=None):
    L = []
    def P(s=""):
        L.append(s)

    P("=" * 78)
    P("CALIBRACION DE LAS %d CORRIDAS  ·  escala %s  ·  compuerta comun %g ns"
      % (len(ordenadas(res)), RAMA, POST_FIJO))
    P("=" * 78)
    P()
    P("  La escala es ENERGIA DEPOSITADA (plantilla Edep_keV de Geant4), no luz.")
    P("  D convierte area en keV depositados y se lleva el quenching de Birks")
    P("  adentro; no es comparable con las D en keVee del resto del proyecto.")
    P()
    P("-" * 78)
    P("1. TODAS LAS CALIBRACIONES")
    P("-" * 78)
    P("  %-4s %-6s %-16s %-18s %8s %7s %8s %6s %s"
      % ("det", "HV", "fuente", "D [mV*ns/keV]", "chi2/ndf", "FWHM", "pulsos",
         "fondo", "medida el"))
    for f in ordenadas(res):
        P("  %-4s %-6.0f %-16s %8.4f ± %-7.4f %8.3f %6.1f%% %8s %6s %s"
          % (f["detector"], f["hv"],
             f["fuente"] + (" (desc)" if f["descartada"] else ""),
             f["D"], f["err_D"], f["chi2_ndf"], 100 * f["fwhm_rel"],
             f["pulsos"] if f["pulsos"] else "?",
             "si" if f["hay_fondo"] else "no", f.get("fecha", "")))
    P()
    P("-" * 78)
    P("2. TEST B  ·  ganancia contra alto voltaje")
    P("-" * 78)
    for (det, fuente), d in sorted(tb.items()):
        aj = d["ajuste"]
        if aj is None:
            P("  %s · %s: %s" % (det, fuente, d["motivo"]))
            continue
        P("  %s · %s   (%d voltajes)" % (det, fuente, aj["n"]))
        P("     ley de potencia D ∝ V^b : b = %.2f ± %.2f   chi2 = %.1f (%d gdl)"
          % (aj["b"], aj["err_b"], aj["chi2_potencia"], aj["gdl"]))
        P("     recta                   : chi2 = %.1f" % aj["chi2_recta"])
        P("     residuos de la potencia : %s"
          % ", ".join("%+.2f%%" % v for v in
                      100 * aj["residuos_potencia"] / aj["D"]))
        P("     b dentro del rango de un PMT (%g a %g): %s"
          % (voltaje.B_PMT[0], voltaje.B_PMT[1],
             "SI" if aj["b_plausible"] else "NO"))
        if not aj["decide"]:
            P("     el rango de voltajes NO separa recta de potencia "
              "(separacion %.2f%% contra error %.2f%%)"
              % (100 * aj["separacion"], 100 * aj["err_rel"]))
        P("     medidas el: %s" % ", ".join(d["fechas"]))
        for av in d.get("avisos", []):
            P("     AVISO: %s" % av)
    P()
    P("-" * 78)
    P("3. TEST A  ·  linealidad del eje canal-energia")
    P("-" * 78)
    for (det, hv), d in sorted(ta.items()):
        lin = d["lin"]
        P("  %s @ %g V" % (det, hv))
        if lin is None:
            P("     no se puede correr: %s" % d["razon"])
            continue
        P("     fuentes: %s" % ", ".join("%s (%.1f keV, D = %.4f)"
                                         % (f["fuente"], f["E_ref"], f["D"])
                                         for f in d["filas"]))
        P("     recta por el origen : D = %.4f" % lin["D_origen"])
        P("     residuos            : %s"
          % ", ".join("%s %+.2f%%" % (n, 100 * r / c)
                      for n, r, c in zip(lin["fuentes"],
                                         lin["residuos_origen"],
                                         lin["canal"])))
        P("     ordenada libre      : c0 = %+.0f ± %.0f mV*ns  ->  %s"
          % (lin["ordenada"], lin["err_ordenada"],
             "compatible con cero" if lin["compatible_con_origen"]
             else "NO compatible con cero"))
        if d["cociente"]:
            q = d["cociente"]
            P("     D(Cs-137)/D(Co-60) = %.4f ± %.4f   (una cadena lineal "
              "daria %.4f en escala Edep)"
              % (q["valor"], q["err"], q["esperado"]))
            P("     o sea %+.1f%% de desacuerdo entre las dos fuentes"
              % (100 * (q["corregido"] - 1)))
        P("     medidas el %s  ->  %s"
          % (", ".join(d["fechas"]),
             "MISMA sesion" if d["misma_sesion"] else
             "SESIONES DISTINTAS: el punto de operacion real puede no ser el "
             "mismo aunque el voltaje nominal lo sea"))
        for av in d.get("avisos", []):
            P("     AVISO: %s" % av)
    P()
    rep = pares_descartadas(res)
    if rep:
        P("-" * 78)
        P("4. CORRIDAS DESCARTADAS  ·  el detector D con la tierra mal")
        P("-" * 78)
        P("  En la sesion del 1-sep el detector D tenia la tierra mal")
        P("  conectada. Las buenas son las del 8-sep; estas se calibran para")
        P("  medir el tamano de la falla y no entran en ningun ajuste.")
        for a, b in rep:
            P("  %s @ %g V · %s:  %.4f ± %.4f  contra  %.4f ± %.4f   "
              "(%+.2f%%)"
              % (a["detector"], a["hv"], a["fuente"], a["D"], a["err_D"],
                 b["D"], b["err_D"], 100 * (b["D"] - a["D"]) / a["D"]))
        P()
    P("-" * 78)
    P("5. LA COMPUERTA DE INTEGRACION")
    P("-" * 78)
    P("  Todas las corridas se integran en la MISMA ventana, %g ns despues del"
      % POST_FIJO)
    P("  onset de cada pulso, que es la que usa la interfaz. Es lo que las hace")
    P("  comparables entre si.")
    P("  La contra: esa ventana NO recoge el pulso entero. La columna 'carga'")
    P("  dice que fraccion del pulso registrado entra, o sea cuanto mas alta")
    P("  seria D con la ventana larga de cada corrida. Mientras la ventana sea")
    P("  la misma para todas, eso no afecta ninguna comparacion; si afecta")
    P("  citar una D como ganancia absoluta.")
    P()
    P("  %-4s %-6s %-16s %8s %10s %12s %s"
      % ("det", "HV", "fuente", "carga", "eventos", "usados", "medida el"))
    for f in ordenadas(res):
        leidos = f.get("pulsos_leidos") or 0
        P("  %-4s %-6.0f %-16s %7.1f%% %10s %11.1f%% %s"
          % (f["detector"], f["hv"],
             f["fuente"] + (" (desc)" if f["descartada"] else ""),
             100 * f["fraccion_carga"], format(leidos, ",d"),
             100 * (f["pulsos"] / leidos) if leidos else 0.0,
             f.get("fecha", "")))
    P()
    ctrl = res.get("_control_compuerta")
    if ctrl:
        P("-" * 78)
        P("6. CONTROL  ·  el mismo test A con otra compuerta")
        P("-" * 78)
        P("  Si el desacuerdo entre fuentes fuera un artefacto de integrar solo")
        P("  %g ns, cambiar la compuerta lo moveria." % POST_FIJO)
        for punto, d in sorted(ctrl.items()):
            P("  %-14s Cs/Co = %.4f con 75 ns  contra  %.4f con %g ns  "
              "(%+.2f%%)"
              % (punto, d["cociente"], d["cociente_36"], POST_FIJO,
                 d["cambio"]))
        P()
    P("-" * 78)
    P("7. LO QUE NO ENTRA, Y POR QUE")
    P("-" * 78)
    P("  · Na-22 (detector C a 1600 V): no hay simulacion de Geant4 de esa")
    P("    fuente con esta geometria, y `fisica.FUENTES` no la tiene. La unica")
    P("    ruta seria la plantilla analitica, que sesga D un +4% respecto de")
    P("    las demas: un punto asi en la recta mediria la plantilla, no el eje.")
    P("  · Las corridas sin fondo propio se calibran con el espectro CRUDO.")
    P("    Prestarles el fondo de otra sesion seria peor: depende del punto de")
    P("    operacion. Se indica corrida por corrida en la tabla 1.")
    P("  · Las corridas del 1-sep del detector D no entran en los tests A ni B:")
    P("    se tomaron con la tierra mal conectada. Las validas son las del")
    P("    8-sep.")
    P()
    texto = "\n".join(L)
    if destino:
        with open(destino, "w", encoding="utf-8") as f:
            f.write(texto + "\n")
    return texto


# --------------------------------------------------------------------------
# 7. main
# --------------------------------------------------------------------------
def main():
    rutas.salida_utf8()
    args = sys.argv[1:]
    rehacer = "--rehacer" in args
    solo_resumen = "--resumen" in args
    con_control = "--control" in args
    if "--area" in args:
        configurar_area(args[args.index("--area") + 1])
    solo = []
    if "--solo" in args:
        solo = []
        for a in args[args.index("--solo") + 1:]:
            if a.startswith("--"):
                break
            solo.append(a)

    os.makedirs(SALIDA, exist_ok=True)
    res = cargar()
    if rehacer:
        if solo:
            # Con `--solo` lo que se pide rehacer son ESAS corridas. Vaciar la
            # tabla entera borraria las otras 22 del archivo de resultados, que
            # es lo que hacia antes y cuesta un cuarto de hora recuperar.
            for k in solo:
                res.pop(k, None)
        else:
            res = {}

    if not solo_resumen:
        cs = [c for c in corridas() if not solo or c.clave in solo]
        print("=" * 78)
        print("%d corridas  ·  escala %s  ·  compuerta comun %g ns  ·  %s  ·  "
              "area: %s  ->  %s"
              % (len(cs), RAMA, POST_FIJO,
                 "con exceso" if EXCESO else "sin exceso de bajo canal",
                 METODO, rutas.relativa(SALIDA)))
        print("=" * 78)
        # Dos corridas del catalogo que apunten al mismo archivo son una
        # contradiccion: una medicion es de un detector. Ver `rutas.duplicados`.
        for grupo in rutas.duplicados(cs):
            print("  AVISO: estas corridas son EL MISMO archivo, byte a byte:")
            for c in grupo:
                print("     %-18s %s" % (c.clave, rutas.relativa(c.muestra)))
        t0 = time.time()
        for i, c in enumerate(cs, 1):
            if c.clave in res and not rehacer:
                print("[%2d/%d] %-18s ya calculada (D = %.4f)"
                      % (i, len(cs), c.clave, res[c.clave]["D"]))
                continue
            print()
            print("[%2d/%d] %-18s %s" % (i, len(cs), c.clave, c.descripcion))
            try:
                res[c.clave] = calibrar_una(c)
            except Exception as e:                              # noqa: BLE001
                print("      FALLA: %s: %s" % (type(e).__name__, e))
                traceback.print_exc()
                continue
            escribir(res)      # incremental: una falla no pierde lo anterior
        print()
        print("  %d de %d corridas calibradas en %.1f min"
              % (len(res), len(cs), (time.time() - t0) / 60))

    if not res:
        print("  no hay resultados que resumir")
        return 1

    if con_control:
        print()
        print("  control: el mismo test A con compuerta de 75 ns")
        res["_control_compuerta"] = control_compuerta(res)
    escribir(res)
    print()
    print("  tabla  -> %s" % rutas.relativa(csv(res)))
    tb, ta = test_b(res), test_a(res)
    print("  resumen-> %s" % rutas.relativa(fig_ganancia(res, tb)))
    print("  resumen-> %s" % rutas.relativa(fig_linealidad(res, ta)))
    texto = informe(res, tb, ta, os.path.join(SALIDA, "RESUMEN.txt"))
    print("  informe-> %s"
          % rutas.relativa(os.path.join(SALIDA, "RESUMEN.txt")))
    print()
    print(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
