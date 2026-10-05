"""
informe.py — arma el informe LaTeX de la campana de calibracion.

Que hace
--------
Lee `figuras/calibraciones/resultados.json` (lo que dejo `calibrar_todo.py`),
dibuja TODAS las figuras que el informe incluye y escribe las tablas y los
numeros en `informe/generado/`. Despues compila `informe.tex` con pdflatex,
si esta.

Las figuras se hacen DEL ANCHO QUE VAN A TENER EN LA PAGINA, con la letra del
documento y en PDF (ver la seccion 3). Es lo que hace que no se vean chicas:
una figura de 15 pulgadas reducida a `\\linewidth` (6.4) pierde el 57% de su
letra, y ninguna fuente grande lo compensa.

El reparto de responsabilidades es a proposito
----------------------------------------------
  * `informe/informe.tex` y `informe/preambulo.tex` son TUYOS: el texto, el
    orden de las secciones, el estilo. Este script no los toca nunca.
  * `informe/generado/` lo reescribe el script en cada corrida. Ahi van las
    tablas y un archivo de macros con cada numero que el texto cita, asi que
    en el .tex se escribe `\\Nc corridas` y no un 25 que manana es otro. Si
    editas algo ahi, se pierde.

O sea que se puede editar el informe a mano y volver a correr esto sin perder
los cambios, que es justamente lo que hace falta para que el documento crezca
con el proyecto.

    python analisis/informe.py             # figuras + tablas + PDF
    python analisis/informe.py --sin-pdf   # solo figuras y tablas
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from analisis import (rutas, espectro, fisica, voltaje, simulacion, modelo,
                      calibracion, area, drs, pulso)
from analisis import sa as sa_mod
from analisis import calibrar_todo as ct

DESTINO = os.path.join(rutas.RAIZ, "informe")
GENERADO = os.path.join(DESTINO, "generado")
FIGURAS = os.path.join(DESTINO, "figuras")

COLOR = ct.COLOR
COLOR_FUENTE = {"Co-60": "#1f4e79", "Cs-137": "#c4630a"}


# --------------------------------------------------------------------------
# 1. Utilidades de LaTeX
# --------------------------------------------------------------------------
def tex(s):
    """Escapa lo que LaTeX se toma a mal. Se usa en TODO lo que viene de datos."""
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def escribir(nombre, texto):
    os.makedirs(GENERADO, exist_ok=True)
    ruta = os.path.join(GENERADO, nombre)
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(texto if texto.endswith("\n") else texto + "\n")
    return ruta


def tabla(columnas, filas, alineacion, nota=None, largo=False):
    """
    Un `tabular` (o `longtable`) listo para `\\input`.

    Sin `\\begin{table}` ni `\\caption` a proposito: el epigrafe y la posicion
    de la tabla son decisiones del documento, no de los datos, asi que viven en
    `informe.tex` donde se los puede editar.
    """
    L = []
    entorno = "longtable" if largo else "tabular"
    L.append(r"\begin{%s}{%s}" % (entorno, alineacion))
    L.append(r"\toprule")
    L.append(" & ".join(r"\textbf{%s}" % c for c in columnas) + r" \\")
    L.append(r"\midrule")
    if largo:
        L.append(r"\endfirsthead")
        L.append(r"\toprule")
        L.append(" & ".join(r"\textbf{%s}" % c for c in columnas) + r" \\")
        L.append(r"\midrule")
        L.append(r"\endhead")
    for fila in filas:
        if fila is None:                       # separador
            L.append(r"\midrule")
            continue
        L.append(" & ".join(str(x) for x in fila) + r" \\")
    L.append(r"\bottomrule")
    L.append(r"\end{%s}" % entorno)
    if nota:
        # La nota va AFUERA de la tabla, no como una fila con `\multicolumn`.
        # Adentro habria que darle un ancho fijo (`p{0.95\linewidth}`) y ese
        # ancho no tiene como saber cuanto mide la tabla: en las anchas se
        # desborda y LaTeX solo avisa con un "Overfull \hbox" que es facil no
        # leer. Con un minipage detras, el ancho lo fija el documento.
        L.append(r"\par\smallskip")
        L.append(r"\begin{minipage}{0.95\linewidth}\footnotesize %s\end{minipage}"
                 % nota)
    return "\n".join(L)


def macro(nombre, valor):
    return r"\newcommand{\%s}{%s}" % (nombre, valor)


# --------------------------------------------------------------------------
# 2. Los datos
# --------------------------------------------------------------------------
def cargar():
    res = ct.cargar()
    if not res:
        raise SystemExit(
            "no hay resultados. Corre primero:  python analisis/calibrar_todo.py")
    return res


def por_clave(res, det, hv, fuente, descartada=False):
    for f in ct.ordenadas(res):
        if (f["detector"] == det and f["hv"] == hv and f["fuente"] == fuente
                and bool(f["descartada"]) == descartada):
            return f
    return None


def detectores(res):
    return sorted({f["detector"] for f in ct.ordenadas(res)})


def voltajes(res):
    return sorted({f["hv"] for f in ct.ordenadas(res)})


def corrida_ejemplo(res):
    """La corrida que ilustra el NCC y el ajuste: detector A a 1600 V con
    Cs-137, que es el punto que tienen los cuatro detectores."""
    return por_clave(res, "A", 1600, "Cs-137") or ct.ordenadas(res)[0]


def resultado_guardado(f):
    """
    El resultado COMPLETO de una fila, tal como lo dejo `calibracion.guardar`
    en `.cache/`, o None. Ahi esta lo que `calibrar_todo` no copia a su tabla:
    la semilla del NCC, su ventana y su rho. Se verifica que la D coincida con
    la de la tabla: en `.cache/` conviven resultados de corridas viejas y de
    la interfaz, y un numero de otra epoca no se notaria en ningun lado.
    """
    ruta = rutas.cache("calibracion_%s_%s.json"
                       % (f["clave"], "lvis" if f["rama"].startswith("Lvis")
                          else "edep"))
    try:
        with open(ruta, encoding="utf-8") as fh:
            r = json.load(fh)
    except (OSError, ValueError):
        return None
    return r if abs(r.get("D", 0) - f["D"]) < 1e-9 else None


def cuentas_por_bin(f):
    """
    Cuentas medianas por bin DENTRO de la ventana de ajuste de esa corrida.

    Es el numero que dice si el sesgo de pesar con las barras del dato (el
    chi2 de Neyman) importa: va como 1/N por bin, asi que con unos pocos
    eventos por bin es del orden del error estadistico y con cientos es
    despreciable.
    """
    e = espectro.cargar("espectro_%s_global" % f["clave"])
    lo, hi = f["ventana"]
    m = (e["canales"] >= lo) & (e["canales"] <= hi)
    y = np.asarray(e["neto"], float)[m]
    return float(np.median(y)) if y.size else float("nan")


def semillas(res):
    """
    (semilla del NCC, D final) de cada corrida.

    `calibrar_todo` no guarda la semilla en su tabla, pero `calibracion.guardar`
    escribe el resultado completo de cada corrida en `.cache/` y ahi si esta.
    Se verifica que el D del JSON coincida con el de la tabla antes de usarlo:
    en `.cache/` conviven resultados de corridas viejas y de la interfaz, y un
    numero de otra epoca metido en la estadistica no se notaria en ningun lado.
    """
    pares = []
    for f in ct.ordenadas(res):
        ruta = rutas.cache("calibracion_%s_%s.json"
                           % (f["clave"],
                              "lvis" if f["rama"].startswith("Lvis") else "edep"))
        if not os.path.exists(ruta):
            continue
        try:
            with open(ruta, encoding="utf-8") as fh:
                r = json.load(fh)
        except (OSError, ValueError):
            continue
        if r.get("semilla_D") and abs(r.get("D", 0) - f["D"]) < 1e-9:
            pares.append((float(r["semilla_D"]), float(f["D"])))
    return pares


# --------------------------------------------------------------------------
# 2-bis. La comparacion con el paper
# --------------------------------------------------------------------------
def paper():
    """
    Los numeros PUBLICADOS del paper, leidos de `replica_paper/auditoria.py`.

    No se tipean aca por lo mismo que no se tipea ningun otro numero del
    informe: ya estan transcritos, con la referencia a la figura de la que
    salen, en un solo lugar. Si algun dia se corrige una lectura, se corrige
    ahi y este informe la toma.

    Devuelve None si `replica_paper/` no esta: el informe se compila igual y
    la seccion de comparacion lo dice, en vez de romperse.
    """
    try:
        from replica_paper import auditoria as au
    except Exception:                                       # noqa: BLE001
        return None
    sds, r2s, cosenos = [], [], []
    for d in au.FIG3.values():
        imp = d["impreso"]
        sds += [imp["sd_azul"], imp["sd_rojo"]]
        r2s += [imp["r2_azul"], imp["r2_rojo"]]
        cosenos.append(d["coseno"])
    v = np.array([b for _, b in au.FIG5], float)
    return dict(
        sd_min=min(sds), sd_max=max(sds), r2_min=min(r2s), r2_max=max(r2s),
        coseno_min=min(cosenos), coseno_max=max(cosenos),
        # dispersion de la Fig. 5, definida como el rango sobre la media, que
        # es como la mide `auditoria.fig5`
        fig5_todas=float((v.max() - v.min()) / v.mean()),
        fig5_sin_angosta=float((v[1:].max() - v[1:].min()) / v[1:].mean()),
        reducciones=au.REDUCCIONES,
        intervalos=_intervalos_fig5())


def _intervalos_fig5():
    """Los seis intervalos de la Fig. 5, en unidades del canal del borde."""
    try:
        from replica_paper.correr import INTERVALOS_FIG5
        return list(INTERVALOS_FIG5)
    except Exception:                                       # noqa: BLE001
        return []


def robustez_intervalo(res, claves, min_bins=30):
    """
    El test de la Fig. 5 del paper, con NUESTRO modelo y NUESTROS datos.

    Ajusta la misma corrida en los seis intervalos anidados que usa el paper
    -importados de `replica_paper`, para que haya una sola definicion- y
    devuelve la dispersion de D. Es la comparacion mas directa que se puede
    hacer: el mismo test, sobre otro detector y con otro modelo.

    No pasa por `calibracion.calibrar` porque ahi la ventana la elige el NCC, y
    justamente lo que se quiere es imponerla. Se arma el `Problema` a mano con
    el espectro cacheado y la plantilla, que es lo mismo que hace la cadena
    despues del NCC.
    """
    from analisis import modelo, sa, simulacion
    intervalos = _intervalos_fig5()
    salida = {}
    for clave in claves:
        f = [x for x in ct.ordenadas(res) if x["clave"] == clave]
        if not f or not intervalos:
            continue
        f = f[0]
        c = rutas.corrida(clave)
        e = espectro.cargar("espectro_%s_global" % clave)
        cen, pesos = simulacion.cargar(c.sim, f["rama"])
        B = f["D"] * f["E_ref"]
        paso = float(np.median(np.diff(e["canales"])))
        filas = []
        for a_, b_ in intervalos:
            lo, hi = a_ * B, b_ * B
            if int(round((hi - lo) / paso)) < min_bins:
                filas.append((lo, hi, float("nan")))
                continue
            p = modelo.Problema(e["canales"], e["neto"], e["error"], cen, pesos,
                                f["E_ref"], lo, hi, exceso=ct.EXCESO)
            cotas = [(0.6 * f["D"], 1.6 * f["D"]), (0.02, 0.30)]
            corr = sa.multiarranque(p.costo, cotas, n=4, semilla=7)
            x = min(corr, key=lambda r: r["costo_pulido"])["x_pulido"]
            filas.append((lo, hi, float(x[0])))
        buenos = [d for _, _, d in filas if np.isfinite(d)]
        if len(buenos) < 2:
            continue
        b = np.array(buenos, float)
        # el mismo estimador que usa la auditoria sobre la Fig. 5: rango/media
        salida[clave] = dict(
            filas=filas, D=f["D"], detector=f["detector"], hv=f["hv"],
            fuente=f["fuente"],
            dispersion=float((b.max() - b.min()) / b.mean()),
            dispersion_sin_angosta=float((b[1:].max() - b[1:].min())
                                         / b[1:].mean()) if len(b) > 2 else
            float("nan"))
    return salida


def metricas_paper(res, ta):
    """
    Nuestros puntos, medidos con las metricas con las que el paper valida.

    La columna que importa es el coseno: entre dos rectas POR EL ORIGEN vale 1
    exactamente, sean cuales sean las pendientes, asi que esa metrica no puede
    ver un desacuerdo de pendiente. Es el motivo por el que este informe
    reporta residuos en keV y no similitudes.
    """
    E = np.linspace(50.0, 1400.0, 400)
    filas = []
    for (det, hv), d in sorted(ta.items()):
        lin, q = d["lin"], d["cociente"]
        if lin is None or not q:
            continue
        cs = por_clave(res, det, hv, "Cs-137")
        co = por_clave(res, det, hv, "Co-60")
        v1, v2 = cs["D"] * E, co["D"] * E
        coseno = float(v1 @ v2 / np.sqrt((v1 @ v1) * (v2 @ v2)))
        Ei = np.asarray(lin["E"], float)
        ch = np.asarray(lin["canal"], float)
        D0 = lin["D_origen"]
        r2 = 1.0 - (float(np.sum((ch - D0 * Ei) ** 2))
                    / float(np.sum((ch - ch.mean()) ** 2)))
        resid_keV = np.asarray(lin["residuos_origen"], float) / D0
        filas.append(dict(punto="%s @ %g V" % (det, hv),
                          desacuerdo=100 * (q["corregido"] - 1),
                          coseno=coseno, r2=r2,
                          sd_keV=float(np.sqrt(np.mean(resid_keV ** 2)))))
    return filas


def igualar_voltaje(f_ref, f, b):
    """
    A que voltaje habria que poner `f` para que tenga la ganancia de `f_ref`.

    Con `D ~ V^b`, `V_pedido = V * (D_ref/D)^(1/b)`. Es la traduccion util del
    cociente de ganancias: decir "el detector D da 3 veces menos que el A" no
    dice que hacer, y decir "hay que subirlo 180 V" si.
    """
    if not b or not np.isfinite(b) or b <= 0:
        return float("nan")
    return f["hv"] * (f_ref["D"] / f["D"]) ** (1.0 / b)


# --------------------------------------------------------------------------
# 3. El estilo de las figuras: al tamano del documento, con su letra
# --------------------------------------------------------------------------
# Lo que hacia que las figuras se vieran chicas en el PDF no era la letra que
# tenian sino el factor con que LaTeX las reducia. Se generaban de 15 in de
# ancho con letra de 8-10 pt y entraban a `\linewidth`, que en este documento
# mide 6.38 in: escala 0.43, o sea que 9 pt quedaban en 4. Subir la fuente no
# lo arregla, lo corre. La regla es la inversa de lo intuitivo: la figura se
# hace DEL ANCHO QUE VA A TENER EN LA PAGINA, con letra de tamano normal, y
# se incluye SIN escalar. Asi 1 pt en la figura es 1 pt en el documento.
#
# Los tres numeros salen del propio documento (informe.tex y preambulo.tex),
# y si alguno cambia alla hay que cambiarlo aca:
CUERPO_PT = 11              # \documentclass[11pt]
ANCHO_PAGINA_CM = 21.0      # a4paper
MARGEN_CM = 2.4             # \usepackage[a4paper,margin=2.4cm]{geometry}
ANCHO_TEXTO = (ANCHO_PAGINA_CM - 2 * MARGEN_CM) / 2.54   # \textwidth, en in

# Con cuerpo de 11 pt: \small = 10 (lo que usa \caption), \footnotesize = 9,
# \scriptsize = 8. Rotulos de los ejes al tamano del epigrafe; numeros de los
# ticks y leyendas un escalon menos; nada por debajo de 8 pt.
ESTILO = {
    "font.size": CUERPO_PT - 1,
    "axes.labelsize": CUERPO_PT - 1,
    "axes.titlesize": CUERPO_PT - 1,
    "xtick.labelsize": CUERPO_PT - 2,
    "ytick.labelsize": CUERPO_PT - 2,
    "legend.fontsize": CUERPO_PT - 2,
    "legend.title_fontsize": CUERPO_PT - 2,
    "figure.titlesize": CUERPO_PT,
    "lines.linewidth": 1.2,
    "lines.markersize": 5,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "grid.linewidth": 0.5,
    "legend.framealpha": 0.9, "legend.edgecolor": "0.85",
    "legend.handlelength": 1.6,
    "axes.titlepad": 4.0, "axes.labelpad": 2.5,
    # La tipografia del documento: `lmodern` en el .tex, Latin Modern aca
    # (ver `fuentes_del_documento`), y la matematica de Computer Modern, que
    # es la misma familia.
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "savefig.format": "pdf",
}
_ESTILO = {"puesto": False, "familia": None}


def fuentes_del_documento():
    """
    Registra en matplotlib la Latin Modern que usa el .tex (`lmodern`).

    MiKTeX la trae en su arbol (`fonts/opentype/public/lm/`) pero no la instala
    como fuente del sistema, asi que matplotlib no la ve: hay que pedirle la
    ruta a `kpsewhich` y agregarla a mano. Con eso la letra de las figuras es
    LA MISMA del texto, no una parecida. Sin MiKTeX en el PATH se queda con la
    serif que haya (DejaVu Serif), que es lo mas cercano. Devuelve el nombre
    de la familia a usar.
    """
    from matplotlib import font_manager
    kpse = shutil.which("kpsewhich")
    if not kpse:
        return "serif"
    n = 0
    for archivo in ("lmroman10-regular.otf", "lmroman10-bold.otf",
                    "lmroman10-italic.otf", "lmroman10-bolditalic.otf"):
        try:
            # kpsewhich escribe UTF-8 y no la pagina de codigos de la consola:
            # con `text=True` una ruta con acentos sale rota y "no existe".
            salida = subprocess.run([kpse, archivo], capture_output=True,
                                    timeout=30).stdout
        except (OSError, subprocess.SubprocessError):
            return "serif"
        ruta = salida.decode("utf-8", "replace").strip()
        if ruta and os.path.exists(ruta):
            font_manager.fontManager.addfont(ruta)
            n += 1
    return "Latin Modern Roman" if n == 4 else "serif"


def aplicar_estilo():
    """Pone el estilo del informe en los rcParams (una vez). Devuelve la
    familia tipografica que quedo."""
    if _ESTILO["puesto"]:
        return _ESTILO["familia"]
    familia = fuentes_del_documento()
    matplotlib.rcParams.update(ESTILO)
    if familia != "serif":
        matplotlib.rcParams["font.serif"] = (
            [familia] + list(matplotlib.rcParams["font.serif"]))
    _ESTILO.update(puesto=True, familia=familia)
    return familia


def figura(ancho=1.0, alto=0.62):
    """
    Una figura del ancho que va a tener en la pagina.

    ancho : fraccion de `\\textwidth` con que se la incluye. Tiene que ser LA
        MISMA que el `width=` del `\\includegraphics`: eso es lo que cierra el
        circulo y deja la letra a su tamano.
    alto : alto como fraccion del ancho (0.62 anda cerca de la razon aurea).
    """
    aplicar_estilo()
    w = ancho * ANCHO_TEXTO
    fig = Figure(figsize=(w, alto * w), layout="constrained")
    FigureCanvasAgg(fig)
    return fig


def guardar(fig, nombre):
    """
    PDF, y sin `bbox_inches="tight"`.

    Vector: la letra es texto y no depende de ningun dpi. Y sin recortar el
    margen a proposito: `tight` deja la figura un poco mas angosta que
    `\\textwidth`, LaTeX la vuelve a estirar, y la letra ya no mide lo que
    dice. Con `constrained` no hace falta.
    """
    os.makedirs(FIGURAS, exist_ok=True)
    ruta = os.path.join(FIGURAS, nombre + ".pdf")
    fig.savefig(ruta)
    return ruta


def etiqueta_panel(ax, letra):
    """
    El "(a)" de un panel, adentro y arriba a la izquierda.

    Reemplaza al titulo: en un informe el titulo es el epigrafe, y un titulo
    de 10 pt no entra en un panel de 2 pulgadas sin pisarse con el vecino.
    """
    ax.text(0.03, 0.96, letra, transform=ax.transAxes, ha="left", va="top",
            fontsize=ESTILO["font.size"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.85))


def pocos_ticks(ax, nx=5, ny=5):
    """En un panel angosto los ticks de fabrica se pisan."""
    from matplotlib.ticker import MaxNLocator
    if ax.get_xscale() == "linear":
        ax.xaxis.set_major_locator(MaxNLocator(nx))
    if ax.get_yscale() == "linear":
        ax.yaxis.set_major_locator(MaxNLocator(ny))


def leyenda_abajo(fig, manijas, etiquetas, ncol):
    """Una leyenda para toda la figura, debajo de los paneles: es lo que cabe
    cuando los paneles son chicos y las series muchas."""
    fig.legend(manijas, etiquetas, loc="outside lower center", ncol=ncol,
               frameon=False, columnspacing=1.2, handletextpad=0.5)


def _sin_notacion(ax):
    """Ticks legibles en un eje logaritmico corto: 0.5, 0.7, 1, ... y sin las
    etiquetas de los ticks menores, que en un decada y media son ruido."""
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    for eje in (ax.xaxis, ax.yaxis):
        eje.set_major_formatter(ScalarFormatter())
        eje.set_minor_formatter(NullFormatter())


# --------------------------------------------------------------------------
# 4. Las figuras
# --------------------------------------------------------------------------
def figura_por_punto(res, ta):
    """
    Los dos espectros de cada punto de operacion, en el MISMO eje de energia.

    La clave del dibujo es que las dos curvas se calibran con la MISMA D -la
    del ajuste conjunto por el origen-, no cada una con la suya. Si cada una
    usara la propia, su borde caeria exacto sobre la marca teorica por
    construccion y la figura no diria nada. Con la D comun, lo que se ve es
    cuanto se corre cada borde respecto de donde deberia estar, que es
    exactamente el desacuerdo que mide el test A.

    Grilla de detectores x voltajes a ancho de pagina, con ejes compartidos:
    los rotulos van una vez por fila y por columna, y lo que cada panel dice
    de si mismo (el cociente, o que tiene una sola fuente) va adentro.
    """
    # Los voltajes son los que tienen al menos un punto testeable; los
    # detectores, TODOS los que midieron en esos voltajes. Asi entra tambien el
    # detector B, que tiene una sola fuente: su panel no es un test -no hay con
    # que testearlo- pero deja ver su espectro y su borde al lado de los demas.
    hvs = sorted({k[1] for k, d in ta.items() if d["lin"] is not None})
    puntos = {k: d for k, d in ta.items() if k[1] in hvs}
    dets = sorted({k[0] for k in puntos})
    n_f, n_c = len(dets), len(hvs)
    fig = figura(1.0, 0.74 * n_f / n_c)
    ejes = fig.subplots(n_f, n_c, squeeze=False, sharex=True, sharey=True)
    chico = ESTILO["legend.fontsize"] - 1

    for i, det in enumerate(dets):
        for j, hv in enumerate(hvs):
            ax = ejes[i][j]
            d = puntos.get((det, hv))
            if d is None:
                ax.set_axis_off()
                continue
            # Con dos fuentes, el eje comun es la D del ajuste conjunto. Con
            # una sola no hay ajuste conjunto: se usa su propia D, y entonces
            # su borde cae sobre la marca por construccion. Eso NO es un test
            # que salio bien, y por eso el panel lo dice.
            solo_una = d["lin"] is None
            D_conj = (d["filas"][0]["D"] if solo_una
                      else d["lin"]["D_origen"])
            for f in d["filas"]:
                e = espectro.cargar("espectro_%s_%s"
                                    % (f["clave"], f.get("metodo", "global")))
                x = e["canales"] / D_conj
                y = np.asarray(e["neto"], float)
                # se normaliza por el plateau de CADA espectro (la mediana en
                # la mitad de abajo de su ventana de ajuste): sin eso el Co-60,
                # que reparte las mismas cuentas en el doble de energia, queda
                # aplastado contra el eje y no se puede comparar la posicion
                # de los bordes, que es lo unico que importa aca.
                lo, hi = f["ventana"]
                m = (e["canales"] > lo) & (e["canales"] < 0.5 * (lo + hi))
                escala = float(np.median(y[m])) if m.any() else float(np.max(y))
                col = COLOR_FUENTE.get(f["fuente"], "0.4")
                ax.step(x, y / max(escala, 1e-9), where="mid", lw=0.8,
                        color=col, alpha=.85, label=f["fuente"])
                # el borde MEDIDO por esa fuente, leido en el eje comun
                for E in fisica.bordes(f["fuente"]):
                    ax.axvline(E * f["D"] / D_conj, color=col, lw=1.2)
            for E in sorted({b for f in d["filas"]
                             for b in fisica.bordes(f["fuente"])}):
                ax.axvline(E, color="0.35", ls=":", lw=1.0)
                ax.annotate("%.0f" % E, xy=(E, 0.03), fontsize=chico,
                            color="0.35", ha="center", rotation=90,
                            xycoords=("data", "axes fraction"))
            # lo que el panel dice de si mismo va arriba y al CENTRO (500 a
            # 850 keV, encima del plateau): a la izquierda el Co-60 arranca
            # por encima del techo del eje, y a la derecha cruzan las lineas
            # de sus dos bordes
            q = d["cociente"]
            if q is None:
                ax.set_title("%s @ %g V" % (det, hv), color="#7a3d05")
                ax.text(0.55, 0.95, "una sola fuente:\nsin test",
                        transform=ax.transAxes, ha="center", va="top",
                        fontsize=chico, color="#7a3d05")
                ax.set_facecolor("#fdf6ef")
            else:
                ax.set_title("%s @ %g V" % (det, hv))
                # el segundo numero es el de la columna "Desac." de la tabla
                # del test A: el cociente con el quenching de Birks descontado
                ax.text(0.53, 0.95, "Cs/Co = %.3f\ndesac. %+.1f %%"
                        % (q["valor"], 100 * (q["corregido"] - 1)),
                        transform=ax.transAxes, ha="center", va="top",
                        fontsize=chico, zorder=5,
                        bbox=dict(boxstyle="round,pad=0.15", fc="w",
                                  ec="none", alpha=.85))
            ax.grid(alpha=.3)
            if i == 0 and j == n_c - 1:
                # a media altura y a la derecha: arriba esta el texto del
                # panel, y pasados 1150 keV solo queda la cola del Co-60
                ax.legend(loc="center right")

    ejes[0][0].set_xlim(0, 1.35 * max(fisica.bordes("Co-60")))
    ejes[0][0].set_ylim(0, 1.9)
    pocos_ticks(ejes[0][0], 5, 4)
    for j in range(n_c):
        ejes[-1][j].set_xlabel("$E$ con la $D$ conjunta [keV]")
    for i in range(n_f):
        ejes[i][0].set_ylabel("cuentas / plateau")
    return guardar(fig, "comparacion_por_punto")


def figura_lineshape(res):
    """
    De donde sale el borde y que ve el detector.

    Tres paneles: la lineshape de Klein-Nishina tal como la escribe el paper
    (Ec. 10), la misma despues de la resolucion real de una de nuestras
    corridas, y el barrido del NCC sobre esa corrida. Los tres con las mismas
    energias, para que se vea que el maximo del barrido cae donde esta el borde.
    Lo que cada panel es lo dice el epigrafe; aca solo llevan su letra.
    """
    from replica_paper import ecuaciones as ec
    from analisis import ncc as mod_ncc

    fig = figura(1.0, 0.36)
    ax = fig.subplots(1, 3)
    E = np.linspace(1.0, 1400.0, 1400)

    def lineshape(fuente):
        y = np.zeros_like(E)
        for eg in fisica.FUENTES[fuente]:
            y = y + ec.dsigma_de(E, eg)
        return y

    # ---- (a) la lineshape teorica ---------------------------------------
    for fuente in ("Cs-137", "Co-60"):
        col = COLOR_FUENTE[fuente]
        y = lineshape(fuente)
        ax[0].plot(E, y / y.max(), color=col, label=fuente)
        for b in fisica.bordes(fuente):
            ax[0].axvline(b, color=col, ls=":", lw=0.9)
    ax[0].set_xlabel("energía del electrón [keV]")
    ax[0].set_ylabel("$d\\sigma/dE$ (normalizada)")
    ax[0].legend(loc="lower left")

    # ---- (b) la misma, con la resolucion medida --------------------------
    f = corrida_ejemplo(res)
    sigma = f["sigma_rel"] * f["E_ref"]
    paso = float(E[1] - E[0])
    k = ec.kernel_gaussiano(sigma, paso, largo_senal=E.size)
    for fuente in ("Cs-137", "Co-60"):
        col = COLOR_FUENTE[fuente]
        z = np.convolve(lineshape(fuente), k, mode="same")
        ax[1].plot(E, z / z.max(), color=col)
        for b in fisica.bordes(fuente):
            ax[1].axvline(b, color=col, ls=":", lw=0.9)
    # abajo a la izquierda: entre 60 y 500 keV las dos curvas van por encima
    # de 0.3, y a la derecha cualquier banda horizontal cruza la caida del
    # Co-60, que baja de 1 a 0 entre 900 y 1250 keV.
    ax[1].text(0.04, 0.05, "$\\sigma$ = %.0f keV\n(%s @ %g V)"
               % (sigma, f["detector"], f["hv"]), transform=ax[1].transAxes,
               ha="left", va="bottom", fontsize=ESTILO["legend.fontsize"],
               color="0.3")
    ax[1].set_xlabel("energía del electrón [keV]")
    ax[1].set_ylabel("cuentas (normalizadas)")

    # ---- (c) el barrido del NCC sobre una corrida real -------------------
    c = rutas.corrida(f["clave"])
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    cen, pesos = simulacion.cargar(c.sim, f["rama"])
    r = mod_ncc.buscar_intervalo(e["canales"], e["neto"], cen, pesos,
                                 fisica.bordes(c.fuente))
    ax[2].plot(r["Bs"], r["rhos"], color="#1f7a3d")
    ax[2].plot([r["B"]], [r["rho"]], "o", ms=6, color="#e02020")
    ax[2].annotate("$\\rho$ = %.3f" % r["rho"], xy=(r["B"], r["rho"]),
                   xytext=(-6, -4), textcoords="offset points", ha="right",
                   va="top", fontsize=ESTILO["legend.fontsize"], color="#e02020")
    ax[2].axvline(f["D"] * f["E_ref"], color="0.4", ls="--", lw=1.0,
                  label="borde del ajuste final")
    ax[2].set_xlabel("$B$ de prueba [mV·ns]")
    ax[2].set_ylabel("$\\rho$ (Pearson)")
    ax[2].legend(loc="lower right")

    for a, letra in zip(ax, "abc"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3)
        pocos_ticks(a, 4, 5)
    return guardar(fig, "teoria_ncc")


def figura_recocido():
    """
    Los tres ingredientes del recocido simulado, sin datos de por medio.

    (a) por que hace falta: una funcion de juguete con varios valles, un
    optimizador local que se queda en el primero y el recocido que sale;
    (b) la regla de Metropolis para tres temperaturas;
    (c) los dos enfriamientos, el que se usa y el que tiene teorema.
    """
    from scipy.optimize import minimize

    fig = figura(1.0, 0.36)
    ax = fig.subplots(1, 3)

    f = lambda x: (x - 1.6) ** 2 / 6.0 - np.cos(3.4 * x) + 1.0
    xs = np.linspace(-4, 6, 900)
    ax[0].plot(xs, f(xs), lw=1.2, color="0.75", zorder=1)
    r = minimize(lambda v: f(v[0]), [-3.2], method="Nelder-Mead")
    ax[0].plot([-3.2], [f(-3.2)], "o", ms=6, color="#c4630a", zorder=4)
    ax[0].annotate("", xy=(r.x[0], f(r.x[0])), xytext=(-3.2, f(-3.2)),
                   arrowprops=dict(arrowstyle="->", color="#c4630a", lw=1.3))
    ax[0].plot([r.x[0]], [f(r.x[0])], "*", ms=12, color="#c4630a", zorder=4)
    rng = np.random.default_rng(4)
    cam, x, T = [], -3.2, 4.0
    c = f(x)
    for _ in range(240):
        p = x + rng.uniform(-1, 1) * 1.4
        cp = f(p)
        if cp <= c or rng.random() < np.exp(-(cp - c) / T):
            x, c = p, cp
        cam.append((x, c))
        T *= 0.975
    cam = np.array(cam)
    sc = ax[0].scatter(cam[:, 0], cam[:, 1], c=np.arange(len(cam)), s=7,
                       cmap="plasma", zorder=3)
    ax[0].plot([1.6], [f(1.6)], "v", ms=8, color="#1f7a3d", zorder=5)
    # a la izquierda del marcador y pegado al piso: arriba estan los puntos
    # del recocido, que convergen justo ahi
    ax[0].annotate("mínimo global", xy=(1.6, f(1.6)), xytext=(-10, 5),
                   textcoords="offset points", ha="right", va="bottom",
                   fontsize=ESTILO["legend.fontsize"], color="#1f7a3d")
    cb = fig.colorbar(sc, ax=ax[0], pad=0.02)
    cb.set_label("iteración")
    ax[0].set_xlabel("parámetro")
    ax[0].set_ylabel("costo")
    # sin leyenda: la estrella y el recorrido los explica el epigrafe, y en
    # este panel (angosto por la barra de color) cualquier leyenda tapa algo

    dE = np.linspace(0, 10, 400)
    for T_, col in ((4.0, "#1f4e79"), (1.0, "#c4630a"), (0.25, "#1f7a3d")):
        ax[1].plot(dE, np.exp(-dE / T_), color=col, label="$T$ = %.2f" % T_)
    ax[1].set_xlabel("empeoramiento del costo, $\\Delta E$")
    ax[1].set_ylabel("$P$ de aceptar")
    ax[1].legend(loc="upper right")

    k = np.arange(1, 400)
    ax[2].loglog(k, 4.0 * 0.95 ** k, color="#c4630a",
                 label="geométrico")
    ax[2].loglog(k, 4.0 / np.log(1.0 + k), color="#8a2be2",
                 label="logarítmico")
    ax[2].set_ylim(1e-8, 1e2)
    ax[2].set_xlabel("ciclo de temperatura")
    ax[2].set_ylabel("$T$")
    ax[2].legend(loc="lower left")

    for a, letra in zip(ax, "abc"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3, which="both" if a is ax[2] else "major")
        pocos_ticks(a, 5, 5)
    return guardar(fig, "teoria_sa")


def _ajuste_de(f):
    """
    `(espectro, problema, resultado)` de una fila: el modelo tal como quedo
    ajustado, reconstruido de lo que se guardo, SIN volver a calibrar.

    El pipeline deja el resultado completo en `.cache/calibracion_<clave>_
    <escala>.json` (ventana, D, sigma, exceso, amplitudes); con eso y el
    espectro cacheado se arma el mismo `Problema` y `prediccion` devuelve la
    curva del ajuste. Recalibrar daria un numero parecido pero no el mismo
    (el recocido es estocastico), y la figura tiene que mostrar EL ajuste que
    esta en las tablas, no otro.
    """
    ruta = rutas.cache("calibracion_%s_%s.json"
                       % (f["clave"], "lvis" if f["rama"].startswith("Lvis")
                          else "edep"))
    with open(ruta, encoding="utf-8") as h:
        r = json.load(h)
    c = rutas.corrida(f["clave"])
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    cen, pesos = simulacion.cargar(c.sim, r["rama"])
    lo, hi = r["ventana"]
    p = modelo.Problema(e["canales"], e["neto"], e["error"], cen, pesos,
                        r["E_ref"], lo, hi, exceso=bool(r.get("exceso")),
                        tipo_exceso=r.get("tipo_exceso") or "exp",
                        ordenada=bool(r.get("libre_ordenada")))
    return e, p, r


def figura_ajuste(f, nombre):
    """
    Un ajuste, en dos paneles: donde caen los bordes en canal, y el espectro
    con el eje ya calibrado. Es la pagina "Resultado" de la interfaz, redibujada
    a tamano de informe.

    Dos decisiones de dibujo que no estan en la interfaz:
    - El eje vertical se corta en 1.5 veces el plateau. El pico de bajo canal
      (el exceso, fuera de la ventana de ajuste) llega a 5 o 6 veces eso en el
      Co-60, y con la escala completa el ajuste -que es lo que se mira- queda
      aplastado en el quinto de abajo del panel. La banda de arriba que sobra
      es donde van las anotaciones de los bordes, que asi no pisan el dato.
    - Una sola leyenda debajo de los dos paneles. Adentro de un panel de 3
      pulgadas una leyenda de tres lineas tapa o el pico o el flanco.
    """
    e, p, r = _ajuste_de(f)
    u = "keVee" if r["rama"].startswith("Lvis") else "keV"
    y_mod = p.prediccion(calibracion.parametros(r))
    # la escala vertical la fija lo que hay DENTRO de la ventana de ajuste
    # (plateau y borde). Medirla desde un cuarto del borde, como hace la
    # interfaz, deja entrar la cola del exceso, que en el Co-60 a 1600 V vale
    # el doble del plateau y aplasta el flanco.
    tope = float(np.max(p.dato)) if p.dato.size else float(np.max(e["neto"]))
    x_max = 1.6 * max(r["bordes_ch"])
    chico = ESTILO["legend.fontsize"] - 1

    fig = figura(1.0, 0.46)
    a1, a2 = fig.subplots(1, 2)
    h_dato = a1.step(e["canales"], e["neto"], where="mid", lw=0.8,
                     color="0.55")[0]
    h_mod = a1.plot(p.canales, y_mod, lw=1.4, color="C1")[0]
    for k, (E, ch, ech) in enumerate(zip(r["bordes_E"], r["bordes_ch"],
                                         r["err_bordes_ch"])):
        h_borde = a1.axvline(ch, color="C3", ls="--", lw=1.0)
        a1.axvspan(ch - ech, ch + ech, color="C3", alpha=.25)
        # en la banda libre de arriba, y escalonadas: con dos bordes a ~190
        # mV*ns una al lado de la otra se pisan
        a1.annotate("%.0f mV·ns (%.0f %s)" % (ch, E, u), xy=(ch, 0),
                    xytext=(ch + 0.015 * x_max, (1.42 - 0.14 * k) * tope),
                    fontsize=chico, color="C3", va="center")
    a1.set_xlim(0, x_max)
    a1.set_ylim(0, 1.5 * tope)
    a1.set_xlabel("área del pulso [mV·ns]")
    a1.set_ylabel("cuentas")

    D = r["D"]
    a2.step(e["canales"] / D, e["neto"] * D, where="mid", lw=0.8, color="0.55")
    a2.plot(p.canales / D, y_mod * D, lw=1.4, color="C1")
    for E in r["bordes_E"]:
        a2.axvline(E, color="C3", ls="--", lw=1.0)
    h_sig = a2.errorbar([r["E_ref"]], [0.35 * tope * D],
                        xerr=[[r["sigma_E"]], [r["sigma_E"]]], fmt="o",
                        color="C0", capsize=3)
    a2.set_xlim(0, 1.6 * max(r["bordes_E"]))
    a2.set_ylim(0, 1.5 * tope * D)
    a2.set_xlabel("energía [%s]" % u)
    a2.set_ylabel("cuentas por %s" % u)
    # el numero del ajuste, arriba a la derecha: pasado el borde el espectro
    # ya cayo a cero y ese rincon esta libre en las dos fuentes
    a2.text(0.97, 0.95, "$D$ = %.4f ± %.4f mV·ns/%s\n$\\chi^2$/ndf = %.2f"
            % (r["D"], r["err_D"], u, r["chi2_ndf"]), transform=a2.transAxes,
            ha="right", va="top", fontsize=ESTILO["legend.fontsize"])
    for a in (a1, a2):
        a.grid(alpha=.3)
        pocos_ticks(a, 6, 5)
    leyenda_abajo(fig, [h_dato, h_mod, h_borde, h_sig],
                  ["espectro neto", "modelo ajustado",
                   "borde Compton ajustado (± error)",
                   "resolución, $\\pm1\\sigma$ = %.0f %s" % (r["sigma_E"], u)],
                  ncol=2)
    return guardar(fig, "ajuste_%s" % nombre)


def figura_cociente(ta):
    """El test A en un numero: D(Cs-137)/D(Co-60) por punto de operacion.
    Va a 0.72 del ancho del texto, que es como la incluye el .tex. La leyenda
    debajo: adentro, con tres detectores y dos lineas de referencia, tapa
    algun punto en cualquier rincon."""
    fig = figura(0.72, 0.70)
    ax = fig.subplots()
    q0 = ct.cociente_esperado()
    manijas = [ax.axhline(q0, color="k", ls="--", lw=1.0,
                          label="cadena lineal: %.3f (Birks)" % q0),
               ax.axhline(1.0, color="0.6", ls=":", lw=0.9, label="cociente 1")]
    con_cociente = set()
    for det in sorted({k[0] for k in ta}):
        pares = [(hv, d["cociente"]) for (dd, hv), d in sorted(ta.items())
                 if dd == det and d["cociente"]]
        if not pares:
            continue
        con_cociente.update(hv for hv, _ in pares)
        manijas.append(ax.errorbar(
            [p[0] for p in pares], [p[1]["valor"] for p in pares],
            yerr=[p[1]["err"] for p in pares], fmt="o-", ms=5, lw=1.0,
            color=COLOR.get(det, "0.4"), capsize=3, label="detector %s" % det))
    # solo los voltajes que tienen cociente: el 1550 V del Co-60 de junio no
    # tiene pareja y dejaba un tercio del eje vacio
    hvs = sorted(con_cociente)
    if hvs:
        ax.set_xticks(hvs)
        ax.set_xlim(hvs[0] - 15, hvs[-1] + 15)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("$D$(Cs-137) / $D$(Co-60)")
    ax.grid(alpha=.3)
    leyenda_abajo(fig, manijas, [m.get_label() for m in manijas], ncol=3)
    return guardar(fig, "cociente_fuentes")


def figura_ganancia(res, tb):
    """
    D(V) en log-log con la ley de potencia, y sus residuos. Antes era un
    recorte de la figura-resumen de `calibrar_todo`; ahora se dibuja al ancho
    del texto, con una sola leyenda debajo de los dos paneles.
    """
    fig = figura(1.0, 0.52)
    a1, a2 = fig.subplots(1, 2)
    manijas, etiquetas = [], []
    for (det, fuente), d in sorted(tb.items()):
        filas, aj = d["filas"], d["ajuste"]
        V = np.array([f["hv"] for f in filas], float)
        D = np.array([f["D"] for f in filas], float)
        err = np.array([f["err_D"] for f in filas], float)
        col, m = COLOR.get(det, "0.4"), ct.MARCA.get(fuente, "^")
        h = a1.errorbar(V, D, yerr=err, fmt=m, color=col, ms=5, capsize=2.5,
                        ls="none")
        manijas.append(h)
        etiquetas.append("%s · %s" % (det, fuente)
                         + ("" if aj is None
                            else "  ($b$ = %.2f ± %.2f)" % (aj["b"], aj["err_b"])))
        if aj is not None:
            vv = np.linspace(V.min() * 0.99, V.max() * 1.01, 60)
            a1.plot(vv, voltaje.evaluar(aj, vv), lw=1.0, color=col, alpha=.75)
            a2.plot(aj["V"], 100 * aj["residuos_potencia"] / aj["D"], m + "-",
                    color=col, ms=5, lw=1.0)
    # las descartadas no entran en ningun ajuste, pero verlas es la mitad del
    # diagnostico: muestran cuanto se corrio la ganancia por la falla de tierra
    vistas = set()
    for f in ct.ordenadas(res):
        if not f["descartada"]:
            continue
        h, = a1.plot([f["hv"]], [f["D"]], "x", ms=8, mew=1.8,
                     color=COLOR.get(f["detector"], "0.4"))
        etq = "%s · %s (descartada: tierra)" % (f["detector"], f["fuente"])
        if etq not in vistas:
            vistas.add(etq)
            manijas.append(h)
            etiquetas.append(etq)
    hvs = sorted({f["hv"] for f in ct.ordenadas(res)})
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xticks(hvs)
    a1.set_xticks([], minor=True)
    a1.set_yticks([0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
    _sin_notacion(a1)
    a1.set_xlabel("alto voltaje [V]")
    a1.set_ylabel("$D$ [mV·ns / keV]")
    a1.grid(alpha=.3, which="both")
    a2.axhline(0, color="k", lw=.7)
    a2.set_xticks(hvs)
    a2.set_xlabel("alto voltaje [V]")
    a2.set_ylabel("($D$ − ajuste) / $D$ [%]")
    a2.grid(alpha=.3)
    etiqueta_panel(a1, "(a)")
    etiqueta_panel(a2, "(b)")
    leyenda_abajo(fig, manijas, etiquetas, ncol=3)
    return guardar(fig, "ganancia_dv")


def figura_detectores(res, tb):
    """
    Comparacion entre detectores, 2x2: ganancia, ganancia relativa al A,
    resolucion y calidad de los ajustes. El color es el detector y el marcador
    la fuente en los cuatro paneles, asi que alcanza con UNA leyenda debajo.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig = figura(1.0, 0.80)
    ejes = fig.subplots(2, 2)
    dets = detectores(res)
    hvs = voltajes(res)
    marca = ct.MARCA

    # ---- (a) la ganancia de cada detector, fuente por fuente -------------
    ax = ejes[0][0]
    ancho = 0.8 / max(len(dets), 1)
    for k, det in enumerate(dets):
        for fuente, hatch in (("Cs-137", None), ("Co-60", "///")):
            xs, ys = [], []
            for i, hv in enumerate(hvs):
                f = por_clave(res, det, hv, fuente)
                if f is None:
                    continue
                xs.append(i + (k - (len(dets) - 1) / 2) * ancho
                          + (0.18 * ancho if fuente == "Co-60" else -0.18 * ancho))
                ys.append(f["D"])
            if xs:
                ax.bar(xs, ys, width=0.36 * ancho, color=COLOR.get(det),
                       hatch=hatch, edgecolor="w", linewidth=.4)
    ax.set_xticks(range(len(hvs)))
    ax.set_xticklabels(["%g V" % v for v in hvs])
    ax.set_ylabel("$D$ [mV·ns / keV]")
    ax.set_title("ganancia, voltaje por voltaje")
    ax.grid(alpha=.3, axis="y")

    # ---- (b) la misma ganancia, relativa al detector A -------------------
    ax = ejes[0][1]
    ax.axhline(1.0, color="k", ls="--", lw=1)      # el A contra si mismo
    for det in dets:
        for fuente in ("Cs-137", "Co-60"):
            xs, ys = [], []
            for hv in hvs:
                f, ref = por_clave(res, det, hv, fuente), \
                         por_clave(res, "A", hv, fuente)
                if f is None or ref is None:
                    continue
                xs.append(hv)
                ys.append(f["D"] / ref["D"])
            if xs:
                ax.plot(xs, ys, marca.get(fuente, "^") + "-",
                        color=COLOR.get(det), ms=5, lw=1)
    ax.set_yscale("log")
    ax.set_yticks([0.3, 0.5, 0.7, 1.0])
    _sin_notacion(ax)
    ax.set_xticks(hvs)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("$D$ / $D$(detector A)")
    ax.set_title("ganancia relativa al detector A")
    ax.grid(alpha=.3, which="both")

    # ---- (c) la resolucion ------------------------------------------------
    ax = ejes[1][0]
    for det in dets:
        for fuente in ("Cs-137", "Co-60"):
            xs, ys = [], []
            for hv in hvs:
                f = por_clave(res, det, hv, fuente)
                if f is None:
                    continue
                xs.append(hv)
                ys.append(100 * f["fwhm_rel"])
            if xs:
                ax.plot(xs, ys, marca.get(fuente, "^") + "-",
                        color=COLOR.get(det), ms=5, lw=1)
    ax.set_xticks(hvs)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("FWHM en $E_\\mathrm{ref}$ [%]")
    ax.set_title("resolución")
    ax.grid(alpha=.3)

    # ---- (d) la calidad del ajuste ---------------------------------------
    ax = ejes[1][1]
    ax.axhline(1.0, color="k", ls="--", lw=1)      # chi2/ndf = 1
    n = 0
    for f in ct.ordenadas(res):
        ax.plot([f["hv"] + (6 if f["fuente"] == "Co-60" else -6)],
                [f["chi2_ndf"]], marca.get(f["fuente"], "^"),
                color=COLOR.get(f["detector"]), ms=5)
        n += 1
    ax.set_xticks(hvs)
    ax.set_xlabel("alto voltaje [V]")
    ax.set_ylabel("$\\chi^2$ / ndf")
    ax.set_title("calidad de los %d ajustes" % n)
    ax.grid(alpha=.3)

    manijas = [Patch(facecolor=COLOR.get(d), label="detector %s" % d)
               for d in dets]
    manijas += [Line2D([], [], marker=marca["Co-60"], color="0.3", ls="none",
                       label="Co-60 (barras rayadas)"),
                Line2D([], [], marker=marca["Cs-137"], color="0.3", ls="none",
                       label="Cs-137")]
    leyenda_abajo(fig, manijas, [m.get_label() for m in manijas], ncol=3)
    return guardar(fig, "comparacion_detectores")


# --------------------------------------------------------------------------
# 4-bis. Los pulsos y su integracion
# --------------------------------------------------------------------------
# La seccion "Del pulso al espectro" del informe. Todo lo que muestra sale de
# los mismos archivos y las mismas funciones que uso la calibracion:
#   - `area.calibrar_ventana` mide la ventana sobre los primeros ~600 pulsos
#     (es determinista: lee siempre los mismos eventos) y devuelve tambien el
#     pulso promedio alineado, que es lo que se dibuja;
#   - `area.areas_de_archivo` con la MISMA compuerta que uso el lote sale del
#     cache `.npz` que la calibracion ya dejo, asi que no se re-parsea nada;
#   - los pulsos crudos se leen con `drs.leer_pulsos`, unos cientos.
# Nada se recalcula de otra manera: si un numero de aca no coincidiera con la
# tabla de calibraciones seria un bug, no una version alternativa.
# --------------------------------------------------------------------------
# 4-bis. Marco teorico, medido: la estadistica del ajuste, el NCC por dentro
#        y los metodos clasicos. Todo lo que estas funciones dibujan o tabulan
#        sale de correr la cadena, no de un numero recordado.
# --------------------------------------------------------------------------
SIGMAS_CLASICOS = (0.03, 0.06, 0.09, 0.12, 0.15, 0.18)
N_CIERRE = 30           # realizaciones Poisson del cierre con verdad conocida
SEMILLA_CIERRE = 11


def _clasicos():
    """
    Los tres metodos clasicos viven en `codigo/etapas/clasicos.py`, junto con
    el modelo analitico (Klein-Nishina convolucionada) con el que se les mide
    el sesgo. Se importan de ahi en vez de copiarlos: una copia se
    desincroniza sola.
    """
    raiz_codigo = os.path.join(rutas.RAIZ, "codigo")
    if raiz_codigo not in sys.path:
        sys.path.insert(0, raiz_codigo)
    from etapas import clasicos, model, physics
    from etapas.fit_single import _espectro_sintetico
    return clasicos, model, physics, _espectro_sintetico


def cierre_conocido(f, n=N_CIERRE, semilla=SEMILLA_CIERRE):
    """
    El cierre con verdad conocida de la cadena completa, sobre una corrida.

    Se toma el ajuste real de la corrida como VERDAD (D, sigma, A), se fabrica
    el espectro que ese modelo predice y se lo muestrea `n` veces con ruido de
    Poisson a la estadistica real. Cada realizacion pasa por la cadena entera
    -NCC, recocido multiarranque, pulido, curvatura- dos veces: pesando con
    las barras del dato (Neyman, que es como se calibra la campana) y con las
    del modelo (Pearson). Lo que sale es el sesgo de la cadena, si el error
    reportado es el error real (cobertura) y cuanto de eso es Neyman.

    Fuera de la ventana de ajuste el modelo no describe el espectro (abajo, el
    exceso de bajo canal y la subida del umbral; arriba, apilamiento y
    cosmicos), asi que ahi la esperanza de Poisson se toma del propio dato
    suavizado: la cadena tiene que elegir la ventana frente al mismo espectro
    que ve en la corrida real, o el cierre probaria el ajuste y no la cadena.

    Cuesta ~2.5 s por ajuste, o sea unos minutos por corrida, y se cachea en
    `.cache/informe_cierre_<clave>.json` con la verdad como clave: si cambia
    el ajuste de la corrida, se rehace solo.
    """
    r = resultado_guardado(f)
    if r is None or r.get("exceso"):
        return None
    D0, s0, A0 = float(r["D"]), float(r["sigma_rel"]), float(r["A"])
    firma = "%s|%d|%d|%.6f|%.6f|%g" % (f["clave"], n, semilla, D0, s0, ct.POST_FIJO)
    ruta = rutas.cache("informe_cierre_%s.json" % f["clave"])
    try:
        with open(ruta, encoding="utf-8") as h:
            viejo = json.load(h)
        if viejo.get("firma") == firma:
            return viejo
    except (OSError, ValueError):
        pass

    from analisis import ncc as mod_ncc
    c = rutas.corrida(f["clave"])
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    cen, pesos = simulacion.cargar(c.sim, r["rama"])
    canales = np.asarray(e["canales"], float)
    lo, hi = r["ventana"]
    mu = A0 * mod_ncc.plantilla_en_canal(canales, cen, pesos, D0, s0, r["E_ref"])
    fuera = (canales < lo) | (canales > hi)
    mu[fuera] = np.maximum(mod_ncc.suavizar(np.asarray(e["neto"], float), 21),
                           0.0)[fuera]
    mu = np.maximum(mu, 0.0)

    rng = np.random.default_rng(semilla)
    filas = []
    t0 = time.time()
    for i in range(n):
        y = rng.poisson(mu).astype(float)
        esp = dict(e)
        esp["neto"], esp["error"] = y, np.sqrt(np.maximum(y, 1.0))
        fila = {}
        for var in ("dato", "modelo"):
            rr, _, _, _ = calibracion.calibrar(
                c, rama=r["rama"], esp=esp, plantilla=(cen, pesos),
                exceso=False, verbose=False, varianza=var)
            fila[var] = dict(D=rr["D"], err_D=rr["err_D"],
                             sigma_rel=rr["sigma_rel"],
                             err_sigma_rel=rr["err_sigma_rel"],
                             chi2_ndf=rr["chi2_ndf"], semilla_D=rr["semilla_D"],
                             rho=rr["rho"], ventana=list(rr["ventana"]))
        filas.append(fila)
        print("    cierre %2d/%d  Neyman %.4f  Pearson %.4f  (%.0f s)"
              % (i + 1, n, fila["dato"]["D"], fila["modelo"]["D"],
                 time.time() - t0))

    m = (canales >= lo) & (canales <= hi)
    salida = dict(firma=firma, clave=f["clave"], n=n, D=D0, sigma_rel=s0,
                  cuentas_bin=float(np.median(mu[m])), filas=filas)
    for var in ("dato", "modelo"):
        Ds = np.array([x[var]["D"] for x in filas])
        errs = np.array([x[var]["err_D"] for x in filas])
        ss = np.array([x[var]["sigma_rel"] for x in filas])
        sem = np.array([x[var]["semilla_D"] for x in filas])
        dv = Ds / D0 - 1.0
        salida[var] = dict(
            sesgo=float(100 * dv.mean()),
            err_sesgo=float(100 * dv.std(ddof=1) / np.sqrt(n)),
            rms=float(100 * np.sqrt(np.mean(dv ** 2))),
            dispersion=float(100 * dv.std(ddof=1)),
            err_medio=float(100 * np.mean(errs) / D0),
            cobertura=int(np.sum(np.abs(Ds - D0) <= errs)),
            sesgo_sigma=float(100 * (ss / s0 - 1.0).mean()),
            semilla_desvio=float(100 * np.mean(np.abs(sem / D0 - 1.0))),
            chi2_ndf=float(np.mean([x[var]["chi2_ndf"] for x in filas])))
    with open(ruta, "w", encoding="utf-8") as h:
        json.dump(salida, h, indent=1)
    return salida


def perfil_chi2(f, n=41, ancho=3.0):
    """
    El perfil de Delta chi2 en D para una corrida, contra la parabola de la
    curvatura.

    Para cada D de una grilla alrededor del optimo se reoptimiza sigma (por
    eso es un PERFIL y no un corte) y se anota cuanto sube el chi2. La
    parabola es lo que la matriz de covarianza supone: que el minimo es
    cuadratico. Donde los dos coinciden, la barra de la curvatura es la barra
    de verdad; donde se separan, gana el perfil.
    """
    from scipy.optimize import minimize
    e, p, r = _ajuste_de(f)
    x = np.asarray(r["x_libres"], float)
    escala = float(r["escala_error"])
    err_curv = float(r["err_libres"][0]) / escala        # sin escalar
    chi2_min = p.costo(x)
    Ds = x[0] + np.linspace(-ancho, ancho, n) * err_curv
    prof = np.empty(n)
    for i, D in enumerate(Ds):
        res = minimize(lambda v: p.costo((D, v[0])), [x[1]],
                       method="Nelder-Mead",
                       options=dict(xatol=1e-7, fatol=1e-5, maxiter=800))
        prof[i] = res.fun - chi2_min
    # el error del perfil: donde Delta chi2 cruza 1, promedio de los dos lados
    izq = Ds[Ds < x[0]]; der = Ds[Ds > x[0]]
    pi = prof[Ds < x[0]]; pd_ = prof[Ds > x[0]]
    e_izq = x[0] - np.interp(1.0, pi[::-1], izq[::-1])
    e_der = np.interp(1.0, pd_, der) - x[0]
    return dict(Ds=Ds, prof=prof, D=float(x[0]), err_curv=err_curv,
                err_perfil=float(0.5 * (e_izq + e_der)), escala=escala,
                chi2_ndf=float(r["chi2_ndf"]), cov=np.asarray(r["cov_libres"]),
                sigma_rel=float(x[1]))


def _elipse(ax, centro, cov, n_sigma=1.0, **kw):
    """El contorno de `n_sigma` de una covarianza 2x2, como polilinea."""
    val, vec = np.linalg.eigh(cov)
    t = np.linspace(0, 2 * np.pi, 200)
    circ = np.stack([np.cos(t), np.sin(t)])
    pts = (vec @ (np.sqrt(np.maximum(val, 0))[:, None] * circ)) * n_sigma
    ax.plot(centro[0] + pts[0], centro[1] + pts[1], **kw)


def figura_estadistica(res, cie, perf):
    """
    De donde sale la barra de error, y si es la barra de verdad.

    (a) el perfil de Delta chi2 en D contra la parabola de la curvatura, con
        las dos lineas de corte: Delta chi2 = 1 (la barra estadistica) y
        Delta chi2 = chi2/ndf (la barra escalada, que es la que se reporta);
    (b) el plano (D, sigma): la elipse de Delta chi2 = 1 del ajuste real y las
        D, sigma recuperadas en las realizaciones del cierre;
    (c) la D recuperada en cada realizacion, con su barra, pesando con el dato
        (Neyman) y con el modelo (Pearson), contra la verdad.
    """
    fig = figura(1.0, 0.36)
    ax = fig.subplots(1, 3)
    azul, naranja, verde = "#1f4e79", "#c4630a", "#1f7a3d"

    # ---- (a) perfil contra parabola ---------------------------------------
    D, ec = perf["D"], perf["err_curv"]
    x = (perf["Ds"] - D) / ec
    # sin leyenda: el perfil ocupa el panel entero (sube hacia los dos lados)
    # y cualquier leyenda tapa justo el cruce con Delta chi2 = 1, que es lo
    # que hay que ver. Que es cada cosa lo dice el epigrafe; las dos lineas
    # horizontales llevan su rotulo pegado.
    ax[0].plot(x, ((perf["Ds"] - D) / ec) ** 2, color="0.6", lw=1.2)
    ax[0].plot(x, perf["prof"], "o", ms=3.2, color=azul)
    ax[0].axhline(1.0, color=naranja, ls="--", lw=1.0)
    ax[0].axhline(perf["chi2_ndf"], color=naranja, ls=":", lw=1.0)
    tope = min(max(perf["prof"].max(), 4.0), 6.0)
    # los rotulos van en el centro, sobre el fondo del valle: a los costados
    # las dos lineas cruzan el perfil y el texto se pisa con los puntos
    ax[0].text(0.0, 1.0, "$\\Delta\\chi^2 = 1$", ha="center", va="bottom",
               fontsize=ESTILO["legend.fontsize"] - 1, color=naranja)
    ax[0].text(0.0, perf["chi2_ndf"], "$\\chi^2/\\mathrm{ndf}$", ha="center",
               va="bottom", fontsize=ESTILO["legend.fontsize"] - 1, color=naranja)
    ax[0].set_xlim(-2.5, 2.5)
    ax[0].set_ylim(0, tope)
    ax[0].set_xlabel("$(D - \\hat D)\\,/\\,\\sigma_{\\mathrm{curv}}$")
    ax[0].set_ylabel("$\\Delta\\chi^2$")

    # ---- (b) el plano D-sigma ---------------------------------------------
    cov = perf["cov"] / perf["escala"] ** 2          # Delta chi2 = 1
    if cie:
        for var, col, mk in (("dato", azul, "o"), ("modelo", naranja, "s")):
            Ds = [f[var]["D"] for f in cie["filas"]]
            ss = [f[var]["sigma_rel"] for f in cie["filas"]]
            ax[1].plot(Ds, ss, mk, ms=3.4, color=col, alpha=0.75,
                       label="Neyman" if var == "dato" else "Pearson")
        ax[1].plot([cie["D"]], [cie["sigma_rel"]], "*", ms=11, color=verde,
                   zorder=5, label="verdad")
        # la dispersion REAL de la cadena: la covarianza empirica de lo que
        # recupero (Neyman), centrada en su media. Contra la elipse de la
        # curvatura mide lo que la barra reportada no contiene (por ejemplo,
        # que la ventana se elige de nuevo en cada realizacion).
        Ds = np.array([f["dato"]["D"] for f in cie["filas"]])
        ss = np.array([f["dato"]["sigma_rel"] for f in cie["filas"]])
        _elipse(ax[1], (Ds.mean(), ss.mean()), np.cov(Ds, ss), 1.0,
                color=azul, ls="--", lw=1.0, label="dispersión del cierre")
    _elipse(ax[1], (perf["D"], perf["sigma_rel"]), cov, 1.0, color="0.3",
            lw=1.1, label="$\\Delta\\chi^2 = 1$ (ajuste real)")
    # la region conjunta del 68 %: con dos parametros el contorno no es
    # Delta chi2 = 1 sino el cuantil 0.683 de una chi2 con 2 gdl (2.30)
    _elipse(ax[1], (perf["D"], perf["sigma_rel"]), cov,
            np.sqrt(DELTA_CHI2_68_2D), color="0.3", lw=1.0, ls=":")
    # los limites los fija la elipse (2.8 sigma a cada lado), no los puntos:
    # asi la elipse entra entera y la leyenda tiene una esquina libre
    sd, ss_ = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])
    ax[1].set_xlim(perf["D"] - 2.8 * sd, perf["D"] + 2.8 * sd)
    ax[1].set_ylim(perf["sigma_rel"] - 2.8 * ss_, perf["sigma_rel"] + 3.4 * ss_)
    ax[1].set_xlabel("$D$ [mV·ns/keV]")
    ax[1].set_ylabel("$\\sigma_{\\mathrm{rel}}$")
    # sin leyenda: cinco entradas no caben en un panel de 2 in sin tapar la
    # nube de puntos; el epigrafe dice que es cada marca y cada elipse

    # ---- (c) realizacion por realizacion ----------------------------------
    if cie:
        k = np.arange(1, cie["n"] + 1)
        for var, col, dx in (("dato", azul, -0.18), ("modelo", naranja, 0.18)):
            Ds = np.array([f[var]["D"] for f in cie["filas"]])
            er = np.array([f[var]["err_D"] for f in cie["filas"]])
            ax[2].errorbar(k + dx, 100 * (Ds / cie["D"] - 1), yerr=100 * er / cie["D"],
                           fmt="o" if var == "dato" else "s", ms=2.8, lw=0.8,
                           color=col, label="Neyman" if var == "dato" else "Pearson")
        ax[2].axhline(0, color=verde, lw=1.2)
        for var, col in (("dato", azul), ("modelo", naranja)):
            ax[2].axhline(cie[var]["sesgo"], color=col, ls="--", lw=0.9)
        ax[2].set_xlabel("realización")
        ax[2].set_ylabel("$D_{\\mathrm{rec}}/D_{\\mathrm{real}} - 1$ [%]")
        ax[2].legend(loc="upper right", fontsize=ESTILO["legend.fontsize"] - 1)

    for a, letra in zip(ax, "abc"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3)
        pocos_ticks(a, 4, 5)
    from matplotlib.ticker import MaxNLocator
    ax[2].xaxis.set_major_locator(MaxNLocator(5, integer=True))
    return guardar(fig, "teoria_estadistica")


def tabla_cierre(cie):
    if not cie:
        return "% sin cierre"
    filas = []
    for var, nombre in (("dato", "dato (Neyman)"), ("modelo", "modelo (Pearson)")):
        d = cie[var]
        filas.append([nombre,
                      "$%+.2f \\pm %.2f$" % (d["sesgo"], d["err_sesgo"]),
                      "%.2f" % d["dispersion"], "%.2f" % d["err_medio"],
                      "%d / %d" % (d["cobertura"], cie["n"]),
                      "$%+.1f$" % d["sesgo_sigma"], "%.2f" % d["chi2_ndf"]])
    # encabezados cortos: con siete columnas y los nombres largos la tabla se
    # salia del ancho del texto (y LaTeX solo avisa con un Overfull)
    return tabla(["barras", "sesgo $D$ [\\%]", "disp. [\\%]",
                  "error rep. [\\%]", "cobertura",
                  "sesgo $\\sigma$ [\\%]", "$\\chi^2$/ndf"],
                 filas, "lcccccc",
                 nota="Sesgo: media de $D_{\\text{rec}}/D_{\\text{real}} - 1$ "
                      "sobre las realizaciones, con el error de esa media. "
                      "Disp.: desviación estándar de las $D$ recuperadas. "
                      "Error rep.: la barra que la cadena reportó, en promedio. "
                      "Si la barra es correcta, coincide con la dispersión y la "
                      "cobertura queda cerca de 68\\,\\%.")


# --------------------------------------------------------------------------
# 4-bis-2. El chi2 y lo que dice (informe/PLAN_estadistica.md, bloque 1)
# --------------------------------------------------------------------------
# Cuanto se aleja un chi2/ndf de 1 no se juzga a ojo: si el modelo es correcto
# chi2_min sigue una distribucion chi2 con ndf grados de libertad, de media ndf
# y desviacion sqrt(2 ndf). Lo que sigue lo mide sobre la campana entera (con
# lo que ya esta en resultados.json) y, sobre la corrida de ejemplo, mira los
# residuos normalizados bin por bin y evalua en el MISMO minimo las tres formas
# de escribir el chi2 de un histograma (Neyman, Pearson y el cociente de
# verosimilitud de Poisson de Baker y Cousins).
PULL_LIM = 6.0          # los pulls de la figura se recortan a +-PULL_LIM


def _chi2_sf(chi2, ndf):
    """P(chi2 >= valor) para una chi2 con ndf grados de libertad."""
    from scipy.stats import chi2 as dist
    return float(dist.sf(chi2, ndf))


def estadistica_chi2(res, f):
    """
    El chi2 de la campana y los residuos de la corrida `f`.

    Los residuos normalizados ("pulls") son (dato - modelo)/sigma bin a bin. Si
    el modelo es correcto y las barras son las de Poisson, son variables
    N(0, 1) independientes: media 0, desviacion 1, un 0.27 % fuera de +-3, y
    sin correlacion entre bines vecinos (las cuentas de dos bines son
    independientes; la resolucion ensancha el MODELO, no correlaciona el
    DATO). Una correlacion entre vecinos es la firma de estructura que el
    modelo no tiene.
    """
    from analisis import ncc
    filas = [x for x in ct.ordenadas(res)
             if x.get("ndf") and np.isfinite(x.get("chi2_ndf") or np.nan)]
    nd = np.array([x["ndf"] for x in filas], float)
    ch = np.array([x["chi2_ndf"] * x["ndf"] for x in filas], float)
    sig = (ch - nd) / np.sqrt(2.0 * nd)

    e, p, r = _ajuste_de(f)
    x = calibracion.parametros(r)
    y = p.prediccion(x)
    pull = (p.dato - y) / p.error
    # las tres formas del chi2, en el mismo punto (el minimo de Neyman)
    chi2_n = float(np.sum(pull ** 2))
    mu = np.maximum(y, 1e-9)
    chi2_p = float(np.sum((p.dato - y) ** 2 / np.maximum(y, 1.0)))
    n = np.clip(p.dato, 0.0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(n > 0, n * np.log(n / mu), 0.0)
    chi2_bc = float(2.0 * np.sum(mu - n + t))
    if abs(chi2_n - float(r["chi2"])) > 1e-3 * float(r["chi2"]):
        # el modelo reconstruido no es el que se ajusto: los pulls no valen
        print("    AVISO: chi2 reconstruido %.2f contra %.2f guardado en %s"
              % (chi2_n, float(r["chi2"]), f["clave"]))
    # el modelo sobre toda la grilla, para ver los residuos fuera de la ventana
    _, A = p.chi2(r["D"], r["sigma_rel"])
    y_full = A * ncc.plantilla_en_canal(e["canales"], p.cen_E, p.pesos,
                                        r["D"], r["sigma_rel"], p.E_ref)
    err_full = np.maximum(np.asarray(e["error"], float), 1.0)
    pull_full = (np.asarray(e["neto"], float) - y_full) / err_full
    ndf = int(r["ndf"])
    return dict(
        ndf=ndf, chi2=float(r["chi2"]), chi2_ndf=float(r["chi2_ndf"]),
        sigma_ndf=float(np.sqrt(2.0 / ndf)),
        n_sigma=float((r["chi2"] - ndf) / np.sqrt(2.0 * ndf)),
        p=_chi2_sf(float(r["chi2"]), ndf),
        chi2_n=chi2_n, chi2_p=chi2_p, chi2_bc=chi2_bc,
        pull=pull, canales=p.canales, lo=p.lo, hi=p.hi,
        canales_full=np.asarray(e["canales"], float), pull_full=pull_full,
        bordes=r["bordes_ch"], umbral=r.get("umbral"),
        media=float(np.mean(pull)), std=float(np.std(pull, ddof=1)),
        fuera=float(100.0 * np.mean(np.abs(pull) > 3.0)),
        autocorr=float(np.corrcoef(pull[:-1], pull[1:])[0, 1]),
        ndf_min=int(nd.min()), ndf_max=int(nd.max()),
        sig_min=float(sig.min()), sig_max=float(sig.max()),
        n_corridas=len(filas))


def figura_pulls(st):
    """
    (a) los residuos normalizados contra el canal, dentro y fuera de la
        ventana; (b) su histograma dentro de la ventana contra N(0,1) y contra
        la gaussiana de ancho sqrt(chi2/ndf), que es lo que el escalado de la
        barra supone.
    """
    fig = figura(1.0, 0.36)
    a1, a2 = fig.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.6, 1.0]))
    azul, gris = "#1f4e79", "0.6"
    xs, ps = st["canales_full"], st["pull_full"]
    x_max = 1.4 * max(st["bordes"])
    m = xs <= x_max
    dentro = (xs >= st["lo"]) & (xs <= st["hi"])
    for sel, col in ((m & ~dentro, gris), (m & dentro, azul)):
        v = ps[sel]
        rec = np.clip(v, -PULL_LIM, PULL_LIM)
        norm = np.abs(v) <= PULL_LIM
        a1.plot(xs[sel][norm], rec[norm], ".", ms=2.2, color=col)
        arriba = ~norm & (v > 0)
        abajo = ~norm & (v < 0)
        a1.plot(xs[sel][arriba], rec[arriba], "^", ms=3, color=col, ls="none")
        a1.plot(xs[sel][abajo], rec[abajo], "v", ms=3, color=col, ls="none")
    a1.axvspan(st["lo"], st["hi"], color=azul, alpha=0.07, lw=0)
    for yv, ls in ((0, "-"), (1, ":"), (-1, ":"), (3, "--"), (-3, "--")):
        a1.axhline(yv, color="k", lw=0.6, ls=ls)
    a1.set_xlim(0, x_max)
    a1.set_ylim(-PULL_LIM - 0.5, PULL_LIM + 0.5)
    a1.set_xlabel("área del pulso [mV·ns]")
    a1.set_ylabel("(dato $-$ modelo) / $\\sigma$")
    etiqueta_panel(a1, "a")

    v = st["pull"]
    bins = np.linspace(-PULL_LIM, PULL_LIM, 37)
    a2.hist(np.clip(v, -PULL_LIM, PULL_LIM), bins=bins, density=True,
            color=azul, alpha=0.55)
    z = np.linspace(-PULL_LIM, PULL_LIM, 300)
    s = np.sqrt(st["chi2_ndf"])
    a2.plot(z, np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi), color="k", lw=1.1)
    a2.plot(z, np.exp(-z ** 2 / (2 * s * s)) / (s * np.sqrt(2 * np.pi)),
            color="k", lw=1.0, ls="--")
    a2.set_xlabel("(dato $-$ modelo) / $\\sigma$, dentro de la ventana")
    a2.set_ylabel("densidad")
    etiqueta_panel(a2, "b")
    for a in (a1, a2):
        a.grid(alpha=.3)
        pocos_ticks(a, 6, 5)
    return guardar(fig, "teoria_pulls")


# --------------------------------------------------------------------------
# 4-bis-3. La incertidumbre de los parametros (PLAN_estadistica, bloque 2)
# --------------------------------------------------------------------------
PROB_1SIGMA = 0.6827                     # la probabilidad de +-1 sigma gaussiano
DELTA_CHI2_68_2D = -2.0 * np.log(1.0 - PROB_1SIGMA)   # cuantil de chi2(2 gdl)


def macros_incertidumbre(res, cie=None, perf=None):
    """
    Los numeros del bloque 2: el paso del hessiano, los cuantiles de Delta chi2
    para una y dos dimensiones, cuantas corridas necesitaron el parche de
    `calibrar_todo.error_usable`, y cuantas realizaciones del cierre caen
    dentro de la region conjunta del 68 % alrededor de la verdad.
    """
    from scipy.stats import chi2 as dist
    filas = ct.ordenadas(res)
    parche = sum(1 for f in filas if ct.error_usable(f)[1])
    L = [macro("PasoCurvRel", _cientifica(modelo.PASO_CURV_REL, 0)),
         macro("DeltaChiUnoD", "%.2f" % dist.ppf(PROB_1SIGMA, 1)),
         macro("DeltaChiDosD", "%.2f" % dist.ppf(PROB_1SIGMA, 2)),
         macro("DeltaChiTresD", "%.2f" % dist.ppf(PROB_1SIGMA, 3)),
         macro("ProbElipseUno", "%.0f" % (100 * dist.cdf(1.0, 2))),
         macro("NErrorParche", "%d" % parche)]
    if cie and perf:
        C = np.asarray(perf["cov"], float) / perf["escala"] ** 2
        Ci = np.linalg.inv(C)
        d = np.array([[x["dato"]["D"] - cie["D"],
                       x["dato"]["sigma_rel"] - cie["sigma_rel"]]
                      for x in cie["filas"]])
        q = np.einsum("ij,jk,ik->i", d, Ci, d)
        L.append(macro("CierreDentroElipse",
                       "%d" % int(np.sum(q <= dist.ppf(PROB_1SIGMA, 2)))))
    return L


def macros_chi2(st):
    if not st:
        return []
    return [macro("EjNdf", "%d" % st["ndf"]),
            macro("EjChiDos", "%.0f" % st["chi2"]),
            macro("EjSigmaNdf", "%.2f" % st["sigma_ndf"]),
            macro("EjDosSigmaNdf", "%.2f" % (2 * st["sigma_ndf"])),
            macro("EjChiSigmas", "%.1f" % st["n_sigma"]),
            macro("EjChiPvalor", _cientifica(st["p"], 0) if st["p"] > 0
                  else "10^{-300}"),
            macro("EjChiNeyman", "%.1f" % st["chi2_n"]),
            macro("EjChiPearson", "%.1f" % st["chi2_p"]),
            macro("EjChiBC", "%.1f" % st["chi2_bc"]),
            macro("EjPullMedia", "%+.2f" % st["media"]),
            macro("EjPullSigma", "%.2f" % st["std"]),
            macro("EjPullFuera", "%.1f" % st["fuera"]),
            macro("EjPullAutocorr", "%+.2f" % st["autocorr"]),
            macro("NdfMin", "%d" % st["ndf_min"]),
            macro("NdfMax", "%d" % st["ndf_max"]),
            macro("ChiSigmasMin", "%.0f" % st["sig_min"]),
            macro("ChiSigmasMax", "%.0f" % st["sig_max"])]


def ncc_detalle(f, n_lista=(12, 25, 50, 100, 200, 400)):
    """
    Lo que el NCC compara en la ventana ganadora, y el efecto del numero de
    puntos de remuestreo.
    """
    from analisis import ncc as mod_ncc
    c = rutas.corrida(f["clave"])
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    cen, pesos = simulacion.cargar(c.sim, f["rama"])
    bordes = fisica.bordes(c.fuente)
    r = mod_ncc.buscar_intervalo(e["canales"], e["neto"], cen, pesos, bordes)
    canales = np.asarray(e["canales"], float)
    y = mod_ncc.suavizar(e["neto"], 21)
    m = (canales >= r["lo"]) & (canales <= r["hi"])
    p = mod_ncc.plantilla_en_canal(canales, cen, pesos, r["D"], r["sigma_rel"],
                                   r["E_ref"])
    u = np.linspace(r["lo"], r["hi"], r["n_pts"])
    s1 = np.interp(u, canales[m], y[m])
    s2 = np.interp(u, canales[m], p[m])
    z = lambda s: (s - s.mean()) / s.std()
    barrido = []
    for n in list(n_lista) + [r["n_pts"]]:
        rn = mod_ncc.buscar_intervalo(e["canales"], e["neto"], cen, pesos,
                                      bordes, n_pts=n)
        barrido.append((n, rn["rho"], rn["B"]))
    barrido.sort()
    return dict(r=r, u=u, z1=z(s1), z2=z(s2), rho=mod_ncc.pearson(s1, s2),
                barrido=barrido, n_auto=r["n_pts"], canales=canales, y=y,
                p=p, m=m)


def figura_ncc_detalle(d):
    """
    (a) las dos curvas estandarizadas en la ventana ganadora: eso es lo unico
        que Pearson ve; (b) una contra otra, con rho; (c) rho maximo y B
        contra el numero de puntos de remuestreo.
    """
    fig = figura(1.0, 0.36)
    ax = fig.subplots(1, 3)
    azul, naranja, verde = "#1f4e79", "#c4630a", "#1f7a3d"
    r = d["r"]

    ax[0].plot(d["u"], d["z1"], color="0.35", lw=1.2, label="espectro suavizado")
    ax[0].plot(d["u"], d["z2"], color=naranja, lw=1.4, label="plantilla en $B$")
    ax[0].axvline(r["B"], color=verde, ls="--", lw=1.0, label="$B$ ganador")
    ax[0].set_xlabel("área [mV·ns]")
    ax[0].set_ylabel("curva estandarizada")
    ax[0].legend(loc="lower left", fontsize=ESTILO["legend.fontsize"] - 1)

    ax[1].plot([-3, 3], [-3, 3], color="0.6", lw=0.9)
    sc = ax[1].scatter(d["z2"], d["z1"], c=d["u"], s=9, cmap="viridis")
    ax[1].text(0.95, 0.07, "$\\rho$ = %.4f" % d["rho"], transform=ax[1].transAxes,
               ha="right", va="bottom", fontsize=ESTILO["legend.fontsize"])
    ax[1].set_xlabel("plantilla (estandarizada)")
    ax[1].set_ylabel("espectro (estandarizado)")
    lim = 1.1 * max(np.abs(d["z1"]).max(), np.abs(d["z2"]).max())
    ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
    cb = fig.colorbar(sc, ax=ax[1], pad=0.02)
    cb.set_label("área [mV·ns]")

    ns = np.array([b[0] for b in d["barrido"]], float)
    rhos = np.array([b[1] for b in d["barrido"]])
    Bs = np.array([b[2] for b in d["barrido"]])
    ax[2].semilogx(ns, rhos, "o-", ms=3.5, color=verde)
    ax[2].axvline(d["n_auto"], color="0.5", ls=":", lw=1.0)
    ax[2].set_xlabel("puntos de remuestreo")
    ax[2].set_ylabel("$\\rho$ máximo", color=verde)
    ax[2].tick_params(axis="y", colors=verde)
    ax2 = ax[2].twinx()
    ax2.semilogx(ns, 100 * (Bs / Bs[-1] - 1), "s--", ms=3.5, color=azul)
    ax2.set_ylabel("$B$ respecto del auto [%]", color=azul)
    ax2.tick_params(axis="y", colors=azul)
    ax2.set_ylim(-1, 1)
    # el "auto" va como linea punteada y no como rotulo: cae pegado al 200
    marcas = [v for v in ns if int(v) != d["n_auto"]]
    ax[2].set_xticks(marcas)
    ax[2].set_xticklabels(["%d" % v for v in marcas])
    ax[2].minorticks_off()

    for a, letra in zip(ax, "abc"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3)
    pocos_ticks(ax[0], 4, 5); pocos_ticks(ax[1], 4, 4)
    return guardar(fig, "teoria_ncc_detalle")


def sesgo_clasicos(sigmas=SIGMAS_CLASICOS, B=600.0, semilla=3):
    """
    Los tres metodos clasicos contra un borde de posicion CONOCIDA.

    Un borde de Klein-Nishina (la linea del Cs-137) convolucionado con una
    gaussiana de ancho sigma = s*B, con y sin ruido de Poisson. Devuelve, por
    sigma, el sesgo de cada metodo en las dos versiones: el "sin ruido" es el
    sesgo del METODO solo.
    """
    clasicos, model, physics, sintetico = _clasicos()
    Eg = float(fisica.FUENTES["Cs-137"][0])
    E_max = float(physics.compton_edge_energy(Eg))
    canales = np.linspace(0, 1000, 1024)
    salida = []
    for s in sigmas:
        sig = s * B
        y = sintetico(canales, [(Eg, 5000.0, B, sig)], semilla=semilla).astype(float)
        ideal = model.forward_model(canales, Eg, 5000.0, B, sig)
        con = clasicos.comparar(canales, y, E_max, B_ref=B)
        sin = clasicos.comparar(canales, ideal, E_max, B_ref=B)
        fila = dict(sigma=s)
        for nombre, k in (("50% de altura", "mitad"), ("diferencia", "dif"),
                          ("FWHM gaussiano", "gauss")):
            fila[k] = 100 * (con[nombre]["borde"] / B - 1)
            fila[k + "_sin"] = 100 * (sin[nombre]["borde"] / B - 1)
        salida.append(fila)
    return dict(filas=salida, B=B, canales=canales, Eg=Eg, E_max=E_max)


def clasicos_reales(res):
    """Los tres clasicos sobre cada corrida, contra la D del ajuste."""
    clasicos, _, _, _ = _clasicos()
    salida = {}
    for f in ct.ordenadas(res):
        e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
        try:
            r = clasicos.comparar(np.asarray(e["canales"], float),
                                  np.asarray(e["neto"], float), f["E_ref"],
                                  errores=np.asarray(e["error"], float),
                                  B_ref=f["D"] * f["E_ref"])
        except Exception:                                  # noqa: BLE001
            continue
        salida[f["clave"]] = {k: 100 * (r[n]["D"] / f["D"] - 1)
                              for n, k in (("50% de altura", "mitad"),
                                           ("diferencia", "dif"),
                                           ("FWHM gaussiano", "gauss"))}
    return salida


def figura_clasicos(res, sc, f_ej):
    """
    (a) un borde sintetico con la resolucion de la corrida de ejemplo, el
        borde verdadero y donde lo pone cada metodo; (b) el sesgo de cada
        metodo contra sigma/B, con ruido (lleno) y sin ruido (trazos), con la
        banda de resoluciones de la campana sombreada.
    """
    clasicos, model, physics, sintetico = _clasicos()
    fig = figura(1.0, 0.40)
    ax = fig.subplots(1, 2)
    B, canales, Eg, E_max = sc["B"], sc["canales"], sc["Eg"], sc["E_max"]
    col = {"mitad": "#1f7a3d", "dif": "#c4630a", "gauss": "#8a2be2"}
    nombre = {"mitad": "50 % de altura", "dif": "mínimo de la derivada",
              "gauss": "escalón gaussiano"}
    clave = {"mitad": "50% de altura", "dif": "diferencia", "gauss": "FWHM gaussiano"}

    s = float(f_ej["sigma_rel"])
    y = sintetico(canales, [(Eg, 5000.0, B, s * B)], semilla=3).astype(float)
    r = clasicos.comparar(canales, y, E_max, B_ref=B)
    # solo el entorno del borde, que es donde los tres metodos miran; con el
    # espectro entero el plateau se lleva el panel y la leyenda lo tapa
    z = (canales > 0.4 * B) & (canales < 1.55 * B)
    ax[0].step(canales[z], y[z], where="mid", color="0.75", lw=0.8)
    v = clasicos._ventana_auto(canales, y, (0.55 * B, 1.45 * B), 3)
    ax[0].plot(canales[z], clasicos._suavizar(y, v, 3)[z], color="0.3", lw=1.1)
    ax[0].axvline(B, color="k", lw=1.5, label="borde verdadero")
    for k in ("mitad", "dif", "gauss"):
        b = r[clave[k]]["borde"]
        if np.isfinite(b):
            ax[0].axvline(b, color=col[k], ls="--", lw=1.1,
                          label="%s (%+.1f %%)" % (nombre[k], 100 * (b / B - 1)))
    ax[0].text(0.04, 0.84, "$\\sigma/B$ = %.3f" % s, transform=ax[0].transAxes,
               va="top", fontsize=ESTILO["legend.fontsize"], color="0.3")
    ax[0].set_xlabel("canal [u.a.]")
    ax[0].set_ylabel("cuentas")
    # abajo a la izquierda: el plateau esta arriba y la caida a la derecha
    ax[0].legend(loc="lower left", fontsize=ESTILO["legend.fontsize"] - 1)

    sig = [x["sigma"] for x in sc["filas"]]
    for k in ("mitad", "dif", "gauss"):
        ax[1].plot(sig, [x[k] for x in sc["filas"]], "o-", ms=3.5, color=col[k],
                   label=nombre[k])
        ax[1].plot(sig, [x[k + "_sin"] for x in sc["filas"]], "--", lw=0.9,
                   color=col[k])
    filas = ct.ordenadas(res)
    smin, smax = min(x["sigma_rel"] for x in filas), max(x["sigma_rel"] for x in filas)
    ax[1].axvspan(smin, smax, color="0.85", alpha=0.6, label="esta campaña")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xlabel("$\\sigma/B$")
    ax[1].set_ylabel("sesgo del borde [%]")
    ax[1].legend(loc="lower right", fontsize=ESTILO["legend.fontsize"] - 1)

    for a, letra in zip(ax, "ab"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3)
        pocos_ticks(a, 5, 5)
    return guardar(fig, "teoria_clasicos")


def tabla_clasicos(sc):
    filas = []
    for x in sc["filas"]:
        filas.append(["%.2f" % x["sigma"]] +
                     ["$%+.1f$" % x[k] for k in ("mitad", "dif", "gauss")] +
                     ["$%+.1f$" % x[k + "_sin"] for k in ("mitad", "dif", "gauss")])
    # "s/r" = sin ruido; el epigrafe del .tex lo dice. Con los nombres largos
    # la tabla se salia 2 cm del ancho del texto.
    return tabla(["$\\sigma/B$", "50\\,\\%", "derivada", "escalón",
                  "50\\,\\% s/r", "derivada s/r", "escalón s/r"],
                 filas, "lcccccc",
                 # sin nota: lo que diria ya lo dice el epigrafe en el .tex
                 nota=None)


def macros_teoria(res, cie, perf, nd, sc, cr):
    L = []
    if perf:
        L += [macro("EjErrCurv", "%.4f" % perf["err_curv"]),
              macro("EjErrPerfil", "%.4f" % perf["err_perfil"]),
              macro("EjErrEscalado", "%.4f" % (perf["err_curv"] * perf["escala"])),
              macro("EjEscala", "%.2f" % perf["escala"]),
              macro("EjChiNdf", "%.2f" % perf["chi2_ndf"])]
    if cie:
        L += [macro("CierreN", "%d" % cie["n"]),
              macro("CierreCuentasBin", "%.0f" % cie["cuentas_bin"]),
              macro("CierreD", "%.4f" % cie["D"])]
        for var, suf in (("dato", "Neyman"), ("modelo", "Pearson")):
            d = cie[var]
            L += [macro("CierreSesgo%s" % suf, "%+.2f" % d["sesgo"]),
                  macro("CierreErrSesgo%s" % suf, "%.2f" % d["err_sesgo"]),
                  macro("CierreDisp%s" % suf, "%.2f" % d["dispersion"]),
                  macro("CierreErrRep%s" % suf, "%.2f" % d["err_medio"]),
                  macro("CierreCob%s" % suf, "%d" % d["cobertura"]),
                  macro("CierreSesgoSigma%s" % suf, "%+.1f" % d["sesgo_sigma"]),
                  macro("CierreSemilla%s" % suf, "%.2f" % d["semilla_desvio"])]
    if nd:
        b = nd["barrido"]
        rho_auto = [x[1] for x in b if x[0] == nd["n_auto"]][0]
        L += [macro("EjNptsAuto", "%d" % nd["n_auto"]),
              macro("EjRhoRecalc", "%.4f" % nd["rho"]),
              macro("NptsMin", "%d" % b[0][0]),
              macro("NptsMax", "%d" % b[-1][0]),
              macro("NptsRhoPocos", "%.4f" % b[0][1]),
              macro("NptsRhoAuto", "%.4f" % rho_auto),
              macro("NptsBDesvMax", "%.1f" % max(
                  100 * abs(x[2] / nd["r"]["B"] - 1) for x in b))]
    if sc:
        f = {x["sigma"]: x for x in sc["filas"]}
        lo, hi = min(f), max(f)
        L += [macro("ClasSigmaMin", "%.2f" % lo), macro("ClasSigmaMax", "%.2f" % hi),
              macro("ClasEgamma", "%.1f" % sc["Eg"]),
              macro("ClasEmax", "%.1f" % sc["E_max"])]
        for k, suf in (("mitad", "Mitad"), ("dif", "Dif"), ("gauss", "Gauss")):
            L += [macro("Clas%sLo" % suf, "%+.1f" % f[lo][k]),
                  macro("Clas%sHi" % suf, "%+.1f" % f[hi][k]),
                  macro("Clas%sSinLo" % suf, "%+.1f" % f[lo][k + "_sin"]),
                  macro("Clas%sSinHi" % suf, "%+.1f" % f[hi][k + "_sin"])]
    if cr:
        L.append(macro("ClasNReales", "%d" % len(cr)))
        for k, suf in (("mitad", "Mitad"), ("dif", "Dif"), ("gauss", "Gauss")):
            v = np.array([x[k] for x in cr.values()])
            v = v[np.isfinite(v)]
            L += [macro("ClasReal%sMed" % suf, "%+.1f" % np.median(v)),
                  macro("ClasReal%sMin" % suf, "%+.1f" % v.min()),
                  macro("ClasReal%sMax" % suf, "%+.1f" % v.max())]
    return L


# --------------------------------------------------------------------------
# 4-ter. La fisica de base: interaccion de fotones, cinematica Compton y
#        Klein-Nishina. Los datos de secciones eficaces son de NIST XCOM
#        (datos/xcom/, con su procedencia en la cabecera de cada archivo); lo
#        demas es analitico y se verifica contra esos datos.
# --------------------------------------------------------------------------
N_AVOGADRO = 6.02214076e23          # 1/mol (valor exacto del SI desde 2019)
R0_CM = 2.818e-13                   # radio clasico del electron [cm], el de ecuaciones.py
LINEAS_CAMPANA = (("Cs-137", 661.657), ("Co-60", 1173.228), ("Co-60", 1332.492))
LADO_DETECTOR_CM = 2 * 2.54         # el largo del cilindro simulado (2 in), ver geant4_geometria


def xcom(nombre):
    """
    Una tabla de NIST XCOM guardada en `datos/xcom/xcom_<nombre>.txt`.

    Devuelve un dict de arrays con la energia en keV y cada proceso en cm2/g,
    mas los constituyentes (Z, fraccion en peso) que XCOM dedujo de la
    formula. La cabecera del archivo dice de donde salio y cuando.
    """
    ruta = os.path.join(rutas.DATOS, "xcom", "xcom_%s.txt" % nombre)
    consti = []
    filas = []
    with open(ruta, encoding="utf-8") as h:
        for linea in h:
            if linea.startswith("#"):
                m = re.search(r"Constituyentes.*?: (.*)$", linea)
                if m:
                    consti = [(int(z), float(w)) for z, w in
                              re.findall(r"(\d+) : ([0-9.]+)", m.group(1))]
                continue
            v = [float(x) for x in linea.split()]
            if len(v) == 8:
                filas.append(v)
    a = np.array(filas)
    return dict(E=a[:, 0] * 1000.0, coherente=a[:, 1], compton=a[:, 2],
                foto=a[:, 3], pares=a[:, 4] + a[:, 5], total=a[:, 6],
                constituyentes=consti, nombre=nombre)


def _interp_log(tabla, clave, E):
    """Interpola una columna de XCOM en log-log, que es como esta tabulada."""
    x, y = np.log(tabla["E"]), np.log(np.maximum(tabla[clave], 1e-300))
    return float(np.exp(np.interp(np.log(E), x, y)))


def electrones_por_gramo(constituyentes):
    """N_e por gramo de un compuesto, de las fracciones en peso de XCOM."""
    # masas atomicas de H y C (g/mol); las demas no hacen falta aca
    masa = {1: 1.00794, 6: 12.0107, 11: 22.98977, 53: 126.90447}
    return N_AVOGADRO * sum(w * z / masa[z] for z, w in constituyentes)


def kn_total(E, mec2=None):
    """
    Seccion eficaz TOTAL de Klein-Nishina por electron [cm2], la integral
    analitica de la Ec. 2 sobre angulo solido. Se verifica contra la
    integral numerica en `figura_cinematica` (coinciden a 1e-8).
    """
    from replica_paper import ecuaciones as ec
    a = np.asarray(E, float) / (ec.M_E if mec2 is None else mec2)
    return 2 * np.pi * R0_CM ** 2 * (
        (1 + a) / a ** 2 * (2 * (1 + a) / (1 + 2 * a) - np.log(1 + 2 * a) / a)
        + np.log(1 + 2 * a) / (2 * a) - (1 + 3 * a) / (1 + 2 * a) ** 2)


def sigma_thomson():
    """El limite de baja energia de Klein-Nishina, 8 pi r0^2 / 3 [cm2]."""
    return 8 * np.pi * R0_CM ** 2 / 3


def energia_dispersada(E, theta):
    """El foton despues de dispersar a `theta` (rad): E / (1 + a (1 - cos))."""
    from replica_paper import ecuaciones as ec
    a = E / ec.M_E
    return E / (1 + a * (1 - np.cos(theta)))


SIM_DIR = os.path.join(rutas.RAIZ, "simulacion")   # la fuente que produjo los .root


def _fuente_sim(nombre):
    """El texto de un archivo de `simulacion/`, o "" si no esta."""
    try:
        with open(os.path.join(SIM_DIR, nombre), encoding="utf-8",
                  errors="replace") as h:
            return h.read()
    except OSError:
        return ""


def geant4_geometria():
    """
    Lo que dice la fuente de la simulacion (`simulacion/`, la copia de
    `ComptonSpectrum` de WSL que produjo los .root de `datos/`): radio y
    semialtura del cilindro, grosores de la carcasa de Al, largo del stub de
    PMT, constante de Birks, distancia de la fuente, tipo de haz, lista de
    fisica, corte de produccion y version de Geant4. Se lee del codigo para
    no tipearlo; lo que no se encuentra queda en None.

    Ojo: en `JcV/ComptonSimulation/` (Windows) hay OTRA simulacion, mas vieja
    y sin carcasa, que NO es la de los .root. Hasta el 21-sep-2026 se leia esa;
    los cuatro numeros que compartian coinciden, la descripcion no.
    """
    det = _fuente_sim(os.path.join("src", "DetectorConstruction.cc"))
    if not det:
        return None
    phys = _fuente_sim(os.path.join("src", "PhysicsList.cc"))
    readme = _fuente_sim("README.md")

    def num(patron, texto=det):
        m = re.search(patron, texto)
        return float(m.group(1)) if m else None

    def cad(patron, texto=det):
        m = re.search(patron, texto)
        return m.group(1) if m else None

    return dict(
        radio=num(r'fScintRadius\(\s*([0-9.]+)\s*\*\s*cm'),
        semialtura=num(r'fScintHalfZ\(\s*([0-9.]+)\s*\*\s*cm'),
        pared=num(r'fHousingWall\(\s*([0-9.]+)\s*\*\s*mm'),
        tapa=num(r'fHousingFront\(\s*([0-9.]+)\s*\*\s*mm'),
        birks=num(r'fBirks\(\s*([0-9.]+)\s*\*\s*mm\s*/\s*MeV'),
        distancia=num(r'fSourceDistance\(\s*([0-9.]+)\s*\*\s*cm'),
        pmt_semilargo=num(r'pmtHalfZ\s*=\s*([0-9.]+)\s*\*\s*cm'),
        haz=cad(r'fBeamType\("(\w+)"\)'),
        material=cad(r'FindOrBuildMaterial\("(G4_PLASTIC_SC_\w+)"\)'),
        # la densidad del material NIST esta en el comentario que lo presenta
        densidad=num(r'([0-9.]+)\s*g/cm3'),
        ventana_pmt=cad(r'FindOrBuildMaterial\("(G4_Pyrex_Glass)"\)'),
        fisica=cad(r'RegisterPhysics\(new (G4EmStandardPhysics_option\d)', phys),
        corte=num(r'SetDefaultCutValue\(\s*([0-9.]+)\s*\*\s*mm', phys),
        version=cad(r'Geant4\s+(\d+\.\d+)', readme))


def figura_interaccion(plastico, nai):
    """
    Los procesos por los que un foton interactua, en el plastico y en un
    inorganico de alto Z, de la misma base de datos: (a) BC-404 (C9H10),
    (b) NaI. Sobre el plastico va ademas Klein-Nishina por electron libre,
    multiplicada por los electrones por gramo, para ver que el "incoherente"
    de XCOM es eso mismo por encima de ~100 keV.
    """
    fig = figura(1.0, 0.50)          # 0.08 mas de alto para la leyenda de abajo (dos filas)
    ax = fig.subplots(1, 2, sharey=True)
    col = {"foto": "#c4630a", "compton": "#1f4e79", "pares": "#8a2be2",
           "coherente": "0.55", "total": "k"}
    nombre = {"foto": "fotoeléctrico", "compton": "Compton", "pares": "pares",
              "coherente": "coherente (Rayleigh)", "total": "total"}
    for a, t in zip(ax, (plastico, nai)):
        m = (t["E"] >= 10) & (t["E"] <= 10000)
        for k in ("total", "compton", "foto", "pares", "coherente"):
            y = t[k][m]
            ok = y > 0
            a.loglog(t["E"][m][ok], y[ok], color=col[k], lw=1.6 if k == "total" else 1.2,
                     ls="-" if k != "total" else "-", label=nombre[k])
        for fuente, E in LINEAS_CAMPANA:
            a.axvline(E, color=COLOR_FUENTE[fuente], ls=":", lw=0.9)
        a.set_xlabel("energía del fotón [keV]")
        a.set_ylim(1e-5, 1e2)
        a.grid(alpha=.3, which="both")
    ne = electrones_por_gramo(plastico["constituyentes"])
    E = np.logspace(1, 4, 200)
    ax[0].loglog(E, kn_total(E) * ne, color=col["compton"], ls="--", lw=1.0,
                 label="Klein-Nishina, e$^-$ libre")
    ax[0].set_ylabel("coeficiente másico [cm$^2$/g]")
    # una sola leyenda debajo de los dos paneles: adentro de cualquiera tapa
    # o la caida del fotoelectrico o su cruce con Compton, que es lo que hay
    # que ver
    manijas, rotulos = ax[0].get_legend_handles_labels()
    leyenda_abajo(fig, manijas, rotulos, 3)
    ax[0].text(0.97, 0.95, "BC-404 (C$_9$H$_{10}$)", transform=ax[0].transAxes,
               ha="right", va="top", fontsize=ESTILO["legend.fontsize"])
    ax[1].text(0.97, 0.95, "NaI", transform=ax[1].transAxes, ha="right", va="top",
               fontsize=ESTILO["legend.fontsize"])
    for a, letra in zip(ax, "ab"):
        etiqueta_panel(a, "(%s)" % letra)
    return guardar(fig, "teoria_interaccion")


def figura_cinematica(plastico):
    """
    (a) la energia del electron de retroceso contra el angulo del foton,
        para las tres lineas: el maximo esta en 180 grados y es el borde;
    (b) Klein-Nishina en angulo, contra el limite de Thomson;
    (c) la seccion eficaz total por electron contra energia: la formula
        analitica, la integral numerica de la Ec. 2 y el "incoherente" de
        XCOM para el plastico llevado a por-electron.
    """
    from replica_paper import ecuaciones as ec
    from scipy.integrate import quad
    fig = figura(1.0, 0.36)
    ax = fig.subplots(1, 3)
    th = np.linspace(0, np.pi, 400)
    grados = np.degrees(th)

    for fuente, E in LINEAS_CAMPANA:
        T = E - energia_dispersada(E, th)
        ax[0].plot(grados, T, color=COLOR_FUENTE[fuente], lw=1.3,
                   ls="-" if E < 1200 else "--")
        ax[0].plot([180], [T[-1]], "o", ms=4, color=COLOR_FUENTE[fuente])
        ax[0].annotate("%.0f" % T[-1], xy=(180, T[-1]), xytext=(-4, 3),
                       textcoords="offset points", ha="right", va="bottom",
                       fontsize=ESTILO["legend.fontsize"] - 1,
                       color=COLOR_FUENTE[fuente])
    ax[0].set_xlabel("ángulo del fotón dispersado [°]")
    ax[0].set_ylabel("energía del electrón [keV]")
    ax[0].set_xlim(0, 185)
    ax[0].set_ylim(0, 1.12 * LINEAS_CAMPANA[2][1])
    # sin rotulo de fuentes: los colores son los de todo el informe y la
    # leyenda del panel (b) los nombra

    for fuente, E in LINEAS_CAMPANA:
        y = ec.dsigma_domega(th, E) / R0_CM ** 2
        ax[1].plot(grados, y, color=COLOR_FUENTE[fuente], lw=1.3,
                   ls="-" if E < 1200 else "--",
                   label="%s, %.0f keV" % (fuente, E))
    ax[1].plot(grados, (1 + np.cos(th) ** 2) / 2, color="0.4", lw=1.0, ls=":",
               label="Thomson ($E \\to 0$)")
    ax[1].set_xlabel("ángulo del fotón dispersado [°]")
    ax[1].set_ylabel("$d\\sigma/d\\Omega$ [$r_0^2$]")
    ax[1].set_xlim(0, 180)
    ax[1].legend(loc="upper right", fontsize=ESTILO["legend.fontsize"] - 1)

    E = np.logspace(1, 4, 60)
    barn = 1e24
    ax[2].loglog(E, kn_total(E) * barn, color="#1f4e79", lw=1.3,
                 label="Klein-Nishina (analítica)")
    Ei = np.logspace(1, 4, 12)
    num = np.array([quad(lambda t: ec.dsigma_domega(t, e) * 2 * np.pi * np.sin(t),
                         0, np.pi)[0] for e in Ei])
    ax[2].loglog(Ei, num * barn, "o", ms=3.5, color="#1f4e79", mfc="none",
                 label="integral numérica de la Ec. 2")
    ne = electrones_por_gramo(plastico["constituyentes"])
    m = (plastico["E"] >= 10) & (plastico["E"] <= 10000)
    ax[2].loglog(plastico["E"][m], plastico["compton"][m] / ne * barn, "s", ms=3,
                 color="#c4630a", mfc="none", label="XCOM incoherente, C$_9$H$_{10}$")
    ax[2].axhline(sigma_thomson() * barn, color="0.4", ls=":", lw=1.0,
                  label="Thomson")
    ax[2].set_xlabel("energía del fotón [keV]")
    ax[2].set_ylabel("$\\sigma$ por electrón [barn]")
    ax[2].legend(loc="lower left", fontsize=ESTILO["legend.fontsize"] - 2)

    for a, letra in zip(ax, "abc"):
        etiqueta_panel(a, "(%s)" % letra)
        a.grid(alpha=.3, which="both" if a is ax[2] else "major")
    pocos_ticks(ax[0], 5, 5); pocos_ticks(ax[1], 5, 5)
    return guardar(fig, "teoria_cinematica")


def tabla_cinematica():
    from replica_paper import ecuaciones as ec
    filas = []
    for fuente, E in LINEAS_CAMPANA:
        a = E / fisica.M_E
        emax = fisica.borde_compton(E)
        filas.append([fuente, "%.1f" % E, "%.3f" % a, "%.1f" % emax,
                      "%.1f" % (100 * emax / E), "%.1f" % (E - emax),
                      "%.1f" % energia_dispersada(E, np.pi / 2)])
    return tabla(["fuente", "$E_\\gamma$ [keV]", "$\\alpha$", "$E_{\\max}$ [keV]",
                  "$E_{\\max}/E_\\gamma$ [\\%]", "retrodispersión [keV]",
                  "$E'$ a 90° [keV]"],
                 filas, "lcccccc", nota=None)


def _cientifica(x, cifras=1):
    """`4\\times10^{-5}`, para usar DENTRO de modo matematico en el .tex."""
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10 ** exp
    return "%.*f\\times10^{%d}" % (cifras, mant, exp)


def macros_fisica(plastico, nai):
    L = []
    g = geant4_geometria() or {}
    rho = g.get("densidad") or 1.032          # G4_PLASTIC_SC_VINYLTOLUENE
    E_cs = LINEAS_CAMPANA[0][1]
    E_co = LINEAS_CAMPANA[2][1]
    for nombre, E in (("Cs", E_cs), ("CoB", E_co)):
        f, c = _interp_log(plastico, "foto", E), _interp_log(plastico, "compton", E)
        L.append(macro("FotoSobreCompton%s" % nombre, _cientifica(f / c, 0)))
        L.append(macro("FotoSobreComptonNaI%s" % nombre,
                       "%.2f" % (_interp_log(nai, "foto", E) / _interp_log(nai, "compton", E))))
        mu = _interp_log(plastico, "total", E) * rho
        L.append(macro("LibreCamino%s" % nombre, "%.1f" % (1 / mu)))
        L.append(macro("ProbInteraccion%s" % nombre,
                       "%.0f" % (100 * (1 - np.exp(-mu * LADO_DETECTOR_CM)))))
    # a que energia el fotoelectrico deja de dominar en NaI y en el plastico
    for t, suf in ((nai, "NaI"), (plastico, "Plastico")):
        E = np.logspace(1, 4, 2000)
        f = np.array([_interp_log(t, "foto", e) for e in E])
        c = np.array([_interp_log(t, "compton", e) for e in E])
        cruce = E[np.argmax(f < c)] if np.any(f < c) else float("nan")
        L.append(macro("CruceFotoCompton%s" % suf, "%.0f" % cruce))
    L += [macro("DensidadBC", "%.3f" % rho),
          macro("LadoDetector", "%.2f" % LADO_DETECTOR_CM),
          macro("SigmaThomson", "%.3f" % (sigma_thomson() * 1e24)),
          macro("SigmaKNCs", "%.3f" % (kn_total(E_cs) * 1e24)),
          macro("ElectronesPorGramo",
                _cientifica(electrones_por_gramo(plastico["constituyentes"]), 2)),
          macro("XcomFraccionH", "%.4f" % dict(plastico["constituyentes"])[1]),
          macro("XcomFraccionC", "%.4f" % dict(plastico["constituyentes"])[6])]
    if g.get("birks"):
        L.append(macro("BirksSim", "%.3f" % g["birks"]))
    # Klein-Nishina contra XCOM: desvio maximo por encima de 100 keV
    ne = electrones_por_gramo(plastico["constituyentes"])
    m = (plastico["E"] >= 100) & (plastico["E"] <= 3000)
    desv = np.abs(kn_total(plastico["E"][m]) * ne / plastico["compton"][m] - 1)
    L.append(macro("KNvsXcom", "%.1f" % (100 * desv.max())))
    m2 = (plastico["E"] >= 10) & (plastico["E"] < 100)
    desv2 = np.abs(kn_total(plastico["E"][m2]) * ne / plastico["compton"][m2] - 1)
    L.append(macro("KNvsXcomBajo", "%.0f" % (100 * desv2.max())))
    # el paper imprime E_gamma/(1+2a) como Ec. 1: cuanto se equivoca
    for nombre, E in (("Cs", E_cs), ("CoA", LINEAS_CAMPANA[1][1])):
        L.append(macro("EcUnoImpresa%s" % nombre, "%.1f" % energia_dispersada(E, np.pi)))
    return L


# --------------------------------------------------------------------------
# 4-quater. La simulacion de Geant4: que hay en los .root, de donde sale la
#           plantilla y que le falta. Los numeros se leen de los archivos de
#           `datos/` y de la fuente en `simulacion/`; el binning se MIDE
#           reajustando la corrida de ejemplo.
# --------------------------------------------------------------------------
BINS_PLANTILLA = (100, 175, 350, 700)     # 350 es el default de simulacion.cargar
QF_BAJO = (50.0, 100.0)                   # [keV] la zona donde Birks muerde mas
QF_ALTO = 700.0                           # [keV] desde donde el quenching satura
_SIM_EXTRA = {}                           # lo que mide figura_simulacion al dibujar


def _primarios_del_nombre(ruta):
    """El `beamOn` que el programa escribio en el nombre: `_3M`, `_10M`, `_20k`."""
    m = re.search(r"_(\d+(?:\.\d+)?)([kM])\.root$", os.path.basename(ruta))
    if not m:
        return None
    return int(float(m.group(1)) * {"k": 1e3, "M": 1e6}[m.group(2)])


def _perfil_qf(edep, qf, n=40):
    """Mediana y percentiles 16/84 de Lvis/Edep en bines logaritmicos de Edep."""
    bordes = np.logspace(np.log10(10.0), np.log10(edep.max()), n + 1)
    filas = []
    for a, b in zip(bordes[:-1], bordes[1:]):
        m = (edep >= a) & (edep < b)
        if m.sum() < 50:
            continue
        p16, p50, p84 = np.percentile(qf[m], (16, 50, 84))
        filas.append((float(np.sqrt(a * b)), float(p50), float(p16), float(p84)))
    return filas


def simulacion_resumen():
    """
    Que tiene cada .root de la campana, leido del archivo: primarios (del
    nombre), filas con deposito, la reparticion por numero de scatterings, las
    absorciones fotoelectricas, el quenching en la zona baja y en la alta y
    el perfil Lvis/Edep contra Edep. Se cachea por tamano + fecha del .root
    porque el de Cs-137 tiene 3.5 millones de filas.
    """
    import uproot
    salida = {}
    for fuente, ruta in (("Co-60", rutas.SIM_CO60), ("Cs-137", rutas.SIM_CS137)):
        if not os.path.exists(ruta):
            continue
        st = os.stat(ruta)
        firma = "%s|%d|%d|v2" % (os.path.basename(ruta), st.st_size, int(st.st_mtime))
        cache = rutas.cache("informe_sim_%s.json" % fuente.replace("-", "").lower())
        try:
            with open(cache, encoding="utf-8") as h:
                viejo = json.load(h)
            if viejo.get("firma") == firma:
                salida[fuente] = viejo
                continue
        except (OSError, ValueError):
            pass
        d = uproot.open(ruta)[simulacion.ARBOL].arrays(
            ["Edep_keV", "Lvis_keVee", "E0_keV", "NCompton", "NPhoto", "FullAbs"],
            library="np")
        edep, lvis, nc = d["Edep_keV"], d["Lvis_keVee"], d["NCompton"]
        qf = lvis / edep
        n = len(edep)
        e_ref = fisica.e_ref(fuente)
        lineas = [(float(v), int(c)) for v, c in
                  zip(*np.unique(np.round(d["E0_keV"], 1), return_counts=True))]
        r = dict(firma=firma, archivo=os.path.basename(ruta),
                 primarios=_primarios_del_nombre(ruta), filas=n, lineas=lineas,
                 nc0=int((nc == 0).sum()), nc1=int((nc == 1).sum()),
                 # sin Compton en el centellador y sin fotoelectrico: el
                 # deposito lo hizo un electron nacido afuera (en la tapa)
                 nc0_electron=int(((nc == 0) & (d["NPhoto"] == 0)).sum()),
                 nc0_edep_mediana=float(np.median(edep[(nc == 0) & (d["NPhoto"] == 0)])),
                 nc2=int((nc == 2).sum()), nc3=int((nc >= 3).sum()),
                 nc_max=int(nc.max()),
                 nphoto=int((d["NPhoto"] > 0).sum()),
                 fullabs=int((d["FullAbs"] == 1).sum()),
                 edep_max=float(edep.max()), lvis_max=float(lvis.max()),
                 qf_bajo=float(np.mean(qf[(edep >= QF_BAJO[0]) & (edep < QF_BAJO[1])])),
                 qf_alto=(float(np.mean(qf[edep >= QF_ALTO]))
                          if (edep >= QF_ALTO).sum() > 100 else None),
                 qf_ajuste=float(np.mean(qf[edep > 0.7 * e_ref])),
                 perfil=_perfil_qf(edep, qf))
        # el borde en LUZ: E_max por el quenching que tiene un electron de esa
        # energia (leido del perfil, no de una formula)
        ex = np.array([p[0] for p in r["perfil"]])
        qm = np.array([p[1] for p in r["perfil"]])
        r["bordes_lvis"] = [float(b * np.interp(b, ex, qm)) for b in fisica.bordes(fuente)]
        with open(cache, "w", encoding="utf-8") as h:
            json.dump(r, h, indent=1)
        salida[fuente] = r
    return salida


def binning_plantilla(f, bins_lista=BINS_PLANTILLA):
    """
    Cuanto cambia D con el binning de la plantilla, MEDIDO reajustando la
    corrida de ejemplo con la misma cadena (NCC, recocido, pulido) y distintos
    numeros de bines. El kernel de la resolucion es decenas de veces mas ancho
    que un bin, asi que no deberia cambiar nada; el numero lo dice.
    Se cachea con la D guardada como firma.
    """
    r = resultado_guardado(f)
    if r is None:
        return None
    firma = "%s|%.6f|%s|%g" % (f["clave"], r["D"], ",".join(map(str, bins_lista)),
                               ct.POST_FIJO)
    ruta = rutas.cache("informe_binning_%s.json" % f["clave"])
    try:
        with open(ruta, encoding="utf-8") as h:
            viejo = json.load(h)
        if viejo.get("firma") == firma:
            return viejo
    except (OSError, ValueError):
        pass
    c = rutas.corrida(f["clave"])
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    filas = []
    for nb in bins_lista:
        cen, pesos = simulacion.cargar(c.sim, r["rama"], bins=nb)
        rr, _, _, _ = calibracion.calibrar(c, rama=r["rama"], esp=e,
                                           plantilla=(cen, pesos), exceso=False,
                                           verbose=False)
        filas.append(dict(bins=nb, ancho=float(cen[1] - cen[0]), D=rr["D"],
                          err_D=rr["err_D"], sigma_E=rr["sigma_E"],
                          chi2_ndf=rr["chi2_ndf"]))
        print("    binning %4d bines (%.1f keV): D = %.4f ± %.4f, chi2/ndf %.3f"
              % (nb, filas[-1]["ancho"], rr["D"], rr["err_D"], rr["chi2_ndf"]))
    salida = dict(firma=firma, clave=f["clave"], D_ref=r["D"], filas=filas)
    with open(ruta, "w", encoding="utf-8") as h:
        json.dump(salida, h, indent=1)
    return salida


def figura_simulacion(sim, f_ej):
    """
    Que hay en la simulacion y que le falta, en cuatro paneles:
    (a) Co-60: energia depositada y luz, con los bordes en las dos escalas;
    (b) Cs-137: la luz separada por numero de scatterings Compton, que es lo
        que rellena la zona por encima del borde;
    (c) el quenching de Birks, Lvis/Edep contra Edep, leido evento a evento;
    (d) la razon dato/modelo del ajuste de la corrida de ejemplo en todo el
        rango: dentro de la ventana es 1, por debajo crece (el exceso que la
        simulacion no tiene) y en el umbral cae a cero.
    """
    import uproot
    fig = figura(1.0, 0.72)
    ax = fig.subplots(2, 2)
    chico = ESTILO["legend.fontsize"] - 1

    # (a) Co-60: Edep contra Lvis
    a = ax[0][0]
    if "Co-60" in sim:
        d = uproot.open(rutas.SIM_CO60)[simulacion.ARBOL].arrays(
            ["Edep_keV", "Lvis_keVee"], library="np")
        bordes = np.linspace(0, 1400, 281)
        cen = 0.5 * (bordes[1:] + bordes[:-1])
        he, _ = np.histogram(d["Edep_keV"], bordes)
        hl, _ = np.histogram(d["Lvis_keVee"], bordes)
        a.step(cen, he / 1e3, where="mid", color="0.45", lw=1.0, label="$E_{dep}$ [keV]")
        a.step(cen, hl / 1e3, where="mid", color=COLOR_FUENTE["Co-60"], lw=1.2,
               label="$L$ [keVee]")
        for b in fisica.bordes("Co-60"):
            a.axvline(b, color="0.45", ls="--", lw=0.8)
        for b in sim["Co-60"]["bordes_lvis"]:
            a.axvline(b, color=COLOR_FUENTE["Co-60"], ls=":", lw=1.0)
        a.set_xlim(0, 1400)
        # el primer bin (depositos de casi cero) es 3 veces el plateau y
        # aplasta el borde, que es lo que hay que ver
        a.set_ylim(0, 1.3 * he[cen > 100].max() / 1e3)
        a.set_ylabel("eventos por 5 keV [$\\times 10^3$]")
        a.legend(loc="upper left", fontsize=chico, frameon=False,
                 bbox_to_anchor=(0.08, 1.0))
    a.set_xlabel("energía [keV o keVee]")
    etiqueta_panel(a, "(a)")

    # (b) Cs-137 por NCompton
    a = ax[0][1]
    if "Cs-137" in sim:
        d = uproot.open(rutas.SIM_CS137)[simulacion.ARBOL].arrays(
            ["Lvis_keVee", "NCompton"], library="np")
        bordes = np.linspace(0, 700, 141)
        cen = 0.5 * (bordes[1:] + bordes[:-1])
        L, nc = d["Lvis_keVee"], d["NCompton"]
        ht, _ = np.histogram(L, bordes)
        h1, _ = np.histogram(L[nc == 1], bordes)
        h2, _ = np.histogram(L[nc == 2], bordes)
        h3, _ = np.histogram(L[nc >= 3], bordes)
        a.step(cen, ht / 1e3, where="mid", color="k", lw=1.2, label="todos")
        a.step(cen, h1 / 1e3, where="mid", color=COLOR_FUENTE["Cs-137"], lw=1.0,
               label="1 Compton")
        a.step(cen, h2 / 1e3, where="mid", color="C0", lw=1.0, label="2 Compton")
        a.step(cen, h3 / 1e3, where="mid", color="C2", lw=1.0, label="$\\geq 3$")
        a.axvline(fisica.bordes("Cs-137")[0], color="0.45", ls="--", lw=0.8)
        a.axvline(sim["Cs-137"]["bordes_lvis"][0], color=COLOR_FUENTE["Cs-137"],
                  ls=":", lw=1.0)
        a.set_xlim(0, 700)
        a.set_ylim(0, 1.3 * ht[cen > 50].max() / 1e3)
        a.set_ylabel("eventos por 5 keVee [$\\times 10^3$]")
        a.legend(loc="upper left", fontsize=chico, frameon=False,
                 bbox_to_anchor=(0.08, 1.0))
    a.set_xlabel("luz [keVee]")
    etiqueta_panel(a, "(b)")

    # (c) el quenching
    a = ax[1][0]
    for fuente, color in (("Co-60", COLOR_FUENTE["Co-60"]),
                          ("Cs-137", COLOR_FUENTE["Cs-137"])):
        if fuente not in sim:
            continue
        p = np.array(sim[fuente]["perfil"])
        a.fill_between(p[:, 0], p[:, 2], p[:, 3], color=color, alpha=0.15, lw=0)
        a.plot(p[:, 0], p[:, 1], color=color, lw=1.2, label=fuente)
    for fuente in ("Cs-137", "Co-60"):
        for b in fisica.bordes(fuente):
            a.axvline(b, color="0.45", ls="--", lw=0.8)
    a.set_xscale("log")
    a.set_xlim(10, 1500)
    a.set_ylim(0.5, 1.0)
    a.set_xlabel("energía depositada $E_{dep}$ [keV]")
    a.set_ylabel("$L / E_{dep}$")
    a.legend(loc="lower right", fontsize=chico, frameon=False)
    _sin_notacion(a)
    etiqueta_panel(a, "(c)")

    # (d) dato/modelo de la corrida de ejemplo
    a = ax[1][1]
    # el modelo en TODO el rango, no solo en la ventana: fuera de ella es
    # donde se ve lo que la simulacion no tiene (la campana no usa el exceso,
    # asi que el modelo es A por la plantilla convolucionada)
    from analisis import ncc as mod_ncc
    e, p, r = _ajuste_de(f_ej)
    cen, pesos = simulacion.cargar(rutas.corrida(f_ej["clave"]).sim, r["rama"])
    can = np.asarray(e["canales"], float)
    y_mod = r["A"] * mod_ncc.plantilla_en_canal(can, cen, pesos, r["D"],
                                                r["sigma_rel"], r["E_ref"])
    dato = np.asarray(e["neto"], float)
    err = np.asarray(e["error"], float)
    lo, hi = r["ventana"]
    m = y_mod > 0.05 * y_mod.max()
    # cuanto supera el dato al modelo por debajo de la ventana (y por encima
    # del umbral): el numero que la seccion cita como "el exceso"
    piso = espectro.fin_del_umbral(e["canales"], e["neto"])
    mm = m & (can >= piso) & (can < lo)
    _SIM_EXTRA["exceso_max"] = float((dato[mm] / y_mod[mm]).max()) if mm.any() else None
    a.axvspan(lo, hi, color="0.9", lw=0)
    a.errorbar(can[m], dato[m] / y_mod[m], yerr=err[m] / y_mod[m], fmt=".",
               ms=2.5, color="0.35", ecolor="0.7", elinewidth=0.5)
    a.axhline(1.0, color="C1", lw=1.2)
    a.axvline(piso, color="C3", ls=":", lw=1.0)
    a.set_xlim(0, hi * 1.05)
    a.set_ylim(0, 3.0)
    a.set_xlabel("área del pulso [mV·ns]")
    a.set_ylabel("dato / modelo")
    a.text(0.5 * (lo + hi), 2.7, "ventana de ajuste", ha="center", fontsize=chico)
    a.text(piso, 2.7, " umbral", ha="left", fontsize=chico, color="C3")
    etiqueta_panel(a, "(d)")

    for fila in ax:
        for a in fila:
            a.grid(alpha=.3)
    return guardar(fig, "simulacion_geant4")


def tabla_simulacion(sim):
    """Una columna por fuente: que tiene cada .root."""
    fuentes = [f for f in ("Co-60", "Cs-137") if f in sim]
    def pct(k, s):
        return "%.1f" % (100.0 * s[k] / s["filas"])
    filas = [
        # el nombre corto: el completo (ComptonSpectrum_<fuente>_<N>.root) no
        # entra en la columna y el texto ya dice como se arma
        ["archivo"] + [r"\texttt{%s}" % tex(sim[f]["archivo"]
                                              .replace("ComptonSpectrum_", "")
                                              .replace(".root", ""))
                       for f in fuentes],
        ["gammas primarios"] + ["%.0f M" % (sim[f]["primarios"] / 1e6) for f in fuentes],
        [r"eventos con $E_{dep} > 0$"] + ["%d" % sim[f]["filas"] for f in fuentes],
        [r"fracción que deposita [\%]"]
        + ["%.1f" % (100.0 * sim[f]["filas"] / sim[f]["primarios"]) for f in fuentes],
        None,
        [r"1 scattering Compton [\%]"] + [pct("nc1", sim[f]) for f in fuentes],
        [r"2 scatterings [\%]"] + [pct("nc2", sim[f]) for f in fuentes],
        [r"3 o más [\%]"] + [pct("nc3", sim[f]) for f in fuentes],
        [r"ningún Compton en el centellador [\%]"] + [pct("nc0", sim[f]) for f in fuentes],
        [r"con absorción fotoeléctrica [\%]"] + [pct("nphoto", sim[f]) for f in fuentes],
        [r"absorción total ($E_{dep} \approx E_0$)"] + ["%d" % sim[f]["fullabs"] for f in fuentes],
        None,
        [r"$\langle L/E_{dep}\rangle$, %.0f--%.0f keV" % QF_BAJO]
        + ["%.3f" % sim[f]["qf_bajo"] for f in fuentes],
        [r"$\langle L/E_{dep}\rangle$, $E_{dep} > %.0f$ keV" % QF_ALTO]
        + [("%.3f" % sim[f]["qf_alto"]) if sim[f]["qf_alto"] else "--" for f in fuentes],
        [r"$\langle L/E_{dep}\rangle$ en la zona del ajuste"]
        + ["%.3f" % sim[f]["qf_ajuste"] for f in fuentes],
        ["bordes en luz [keVee]"]
        + [", ".join("%.0f" % b for b in sim[f]["bordes_lvis"]) for f in fuentes],
    ]
    return tabla(["", *fuentes], filas, "lrr",
                 nota="Las fracciones por número de scatterings son sobre los "
                      "eventos con depósito. «Ningún Compton en el centellador» "
                      "son casi todos eventos cuyo depósito lo hace un electrón "
                      "Compton nacido en la tapa de aluminio que entra al "
                      "plástico (mediana de %.0f keV en el Cs-137); el resto "
                      "son fotoeléctricos directos. La "
                      % sim[fuentes[-1]]["nc0_edep_mediana"] +
                      "zona del ajuste es $E_{dep} > 0{,}7\\,E_{\\text{ref}}$, "
                      "que es donde se promedia el quenching que convierte "
                      "$D_{\\text{luz}}$ en $D_{\\text{dep}}$.")


def macros_simulacion(g, sim, binning):
    L = []
    g = g or {}
    for k, nombre, fmt in (("radio", "SimRadio", "%.2f"),
                           ("semialtura", "SimSemialtura", "%.2f"),
                           ("pared", "SimPared", "%.0f"), ("tapa", "SimTapa", "%.0f"),
                           ("birks", "SimBirks", "%.3f"),
                           ("distancia", "SimDistancia", "%.0f"),
                           ("pmt_semilargo", "SimPmtSemilargo", "%.1f"),
                           ("corte", "SimCorte", "%.1f"),
                           ("densidad", "SimDensidad", "%.3f")):
        if g.get(k) is not None:
            L.append(macro(nombre, fmt % g[k]))
    for k, nombre in (("haz", "SimHaz"), ("material", "SimMaterial"),
                      ("fisica", "SimFisica"), ("version", "SimVersion"),
                      ("ventana_pmt", "SimVentanaPmt")):
        if g.get(k):
            L.append(macro(nombre, r"\texttt{%s}" % tex(g[k])))
    if g.get("semialtura"):
        L.append(macro("SimLargo", "%.2f" % (2 * g["semialtura"])))
    # con cuantos bines se histograma la plantilla en la campana: el default
    # de simulacion.cargar, que es lo que modelo.armar usa
    L.append(macro("SimBinsPlantilla", "%d" % _defaults(simulacion.cargar)["bins"]))
    for fuente, suf in (("Co-60", "Co"), ("Cs-137", "Cs")):
        if fuente not in sim:
            continue
        s = sim[fuente]
        L += [macro("SimPrimarios%s" % suf, "%.0f" % (s["primarios"] / 1e6)),
              macro("SimFilas%s" % suf, "%d" % s["filas"]),
              macro("SimFracDeposita%s" % suf, "%.1f" % (100.0 * s["filas"] / s["primarios"])),
              macro("SimFracSimple%s" % suf, "%.1f" % (100.0 * s["nc1"] / s["filas"])),
              macro("SimFracMultiple%s" % suf,
                    "%.1f" % (100.0 * (s["nc2"] + s["nc3"]) / s["filas"])),
              macro("SimFracFoto%s" % suf, "%.1f" % (100.0 * s["nphoto"] / s["filas"])),
              macro("SimNcMax%s" % suf, "%d" % s["nc_max"]),
              macro("SimQfBajo%s" % suf, "%.3f" % s["qf_bajo"]),
              macro("SimQfAjuste%s" % suf, "%.3f" % s["qf_ajuste"]),
              macro("SimEdepMax%s" % suf, "%.0f" % s["edep_max"]),
              macro("SimBordeLuz%s" % suf, "%.0f" % s["bordes_lvis"][0])]
        if s["qf_alto"]:
            L.append(macro("SimQfAlto%s" % suf, "%.3f" % s["qf_alto"]))
        if len(s["lineas"]) == 2:
            L.append(macro("SimRazonLineas%s" % suf,
                           "%.3f" % (s["lineas"][0][1] / s["lineas"][1][1])))
    if "Co-60" in sim:
        L.append(macro("SimBordeLuzCoB", "%.0f" % sim["Co-60"]["bordes_lvis"][1]))
    if _SIM_EXTRA.get("exceso_max"):
        L.append(macro("SimExcesoMax", "%.1f" % _SIM_EXTRA["exceso_max"]))
    if binning:
        fl = binning["filas"]
        Ds = np.array([x["D"] for x in fl])
        L += [macro("BinsMin", "%d" % fl[0]["bins"]), macro("BinsMax", "%d" % fl[-1]["bins"]),
              macro("BinsAnchoMin", "%.1f" % fl[-1]["ancho"]),
              macro("BinsAnchoMax", "%.1f" % fl[0]["ancho"]),
              macro("BinsDesvMax", "%.2f" % (100 * np.max(np.abs(Ds / binning["D_ref"] - 1)))),
              macro("BinsSigmaE", "%.0f" % np.mean([x["sigma_E"] for x in fl])),
              macro("BinsSigmaSobreBin", "%.0f" % (np.mean([x["sigma_E"] for x in fl])
                                                   / fl[-1]["ancho"])),
              macro("BinsErrD", "%.4f" % np.mean([x["err_D"] for x in fl]))]
        # los nombres de macro no admiten digitos: una letra por binning
        for letra, x in zip("abcdefgh", fl):
            L += [macro("BinsN%s" % letra, "%d" % x["bins"]),
                  macro("BinsD%s" % letra, "%.4f" % x["D"])]
    return L


_PULSOS = {}
# Hasta donde se compara la carga acumulada entre corridas [ns tras el onset].
# La fraccion "de la carga del pulso registrado" (columna Carga de la tabla de
# calibraciones) depende de cuanto registro hay despues del onset, que va de
# ~47 ns (detector D, 8-sep) a ~160 (Co-60 de junio): el "100%" de cada
# corrida es donde se le acaba el registro. Para ver si la FORMA es la misma
# hay que normalizar a un tiempo que todos los registros contengan; 40 ns
# queda dentro del mas corto. Medido: Q(12)/Q(40) vale 0.82-0.83 en las 25
# corridas, mientras que Q(12)/Q(registro) va de 0.72 a 0.82.
T_COMUN = 40.0
_BANDERAS = (("saturados", pulso.MAL_SATURA), ("truncados", pulso.MAL_TRUNCADO),
             ("sin_pre", pulso.MAL_SIN_PRE), ("sin_pico", pulso.MAL_SIN_PICO))


def _defaults(fn):
    """Los valores por defecto de la firma de `fn`, para citarlos como macro
    en vez de copiarlos."""
    import inspect
    return {k: v.default for k, v in inspect.signature(fn).parameters.items()}


def pulsos_de(f):
    """
    Lo que hace falta saber de los pulsos de una corrida, medido sobre el
    propio archivo. Se calcula una vez por corrida (dict `_PULSOS`).

    Devuelve: la corrida, el canal, la ventana medida (`area.calibrar_ventana`),
    las areas de todos los eventos con la compuerta del lote, la mascara de
    eventos buenos, la fraccion marcada por cada bandera, la carga que recoge
    la compuerta, los percentiles del onset y de la amplitud, y el perfil de
    area/amplitud (el test de que la forma del pulso no depende del tamano).
    """
    if f["clave"] in _PULSOS:
        return _PULSOS[f["clave"]]
    c = rutas.corrida(f["clave"])
    canal, _ = drs.detectar_canal(c.muestra)
    v = area.calibrar_ventana(c.muestra, canal=canal, verbose=False)
    d = area.areas_de_archivo(c.muestra, metodo=f["metodo"], pre=v["pre"],
                              post=f["post"], canal=canal, verbose=False)
    ok = area.buenos(d)
    b = d["bandera"]
    a, A = d["area"][ok], d["amplitud"][ok]
    m = A > 10.0
    cx, cy, ce = area.perfil(A[m], a[m] / A[m], bins=12)
    Q = pulso.integral_acumulada(v["tau_rel"], v["forma"])
    espacio = float(np.median(v["espacio"]))      # registro tras el onset
    p = dict(
        corrida=c, canal=canal, ventana=v, areas=d, buenos=ok,
        espacio=espacio, Q=Q,
        # la carga en la compuerta relativa a la acumulada hasta T_COMUN: es
        # la que se puede comparar entre corridas con registros distintos
        carga_comun=float(_Q_en(v["tau_rel"], Q, f["post"])
                          / _Q_en(v["tau_rel"], Q, T_COMUN)),
        formato=drs.formato(c.muestra),
        n_eventos=int(d["n_eventos"]), n_pulsos=int(d["area"].size),
        n_buenos=int(ok.sum()),
        marcados={n: float((b & bit).astype(bool).mean()) for n, bit in _BANDERAS},
        marcados_total=float((b != 0).mean()),
        descartados=float(1.0 - ok.mean()),
        carga=area.fraccion_de_carga(v, f["post"]),
        onset=np.percentile(d["t_onset"][ok], [5, 50, 95]),
        amplitud=np.percentile(A, [50, 99]),
        area=np.percentile(a, [50, 99]),
        perfil=(cx, cy, ce),
        # cuanto se mueve area/amplitud de punta a punta del rango de
        # amplitudes, relativo a su media: es el tamano de lo que la forma del
        # pulso cambia con el tamano del pulso, visto por esta compuerta
        variacion_forma=float((cy.max() - cy.min()) / cy.mean()) if cy.size
        else float("nan"),
        # lo mismo, solo sobre los tramos por encima de la amplitud mediana:
        # abajo la amplitud (una sola muestra, con ~1 mV de ruido) y el onset
        # al 10% se miden peor, y el cociente cae por eso y no por la forma
        variacion_alta=_variacion(cx, cy, float(np.median(A))))
    _PULSOS[f["clave"]] = p
    return p


def _variacion(cx, cy, desde):
    """Rango sobre media de `cy` en los tramos con `cx >= desde`."""
    m = np.asarray(cx) >= desde
    if m.sum() < 2:
        return float("nan")
    y = np.asarray(cy)[m]
    return float((y.max() - y.min()) / y.mean())


def _Q_en(tau, Q, t):
    """La integral acumulada `Q` (sobre la grilla `tau`) evaluada en `t`."""
    return float(Q[min(int(np.searchsorted(tau, float(t))), Q.size - 1)])


def pulso_tipico(p, n=300):
    """
    Un pulso REAL de la corrida para dibujar la integracion paso a paso: sin
    banderas, con la amplitud mas cercana a la mediana de la corrida y el
    onset dentro del cuarto central. Devuelve `(tau, v, r)` con `r` lo que
    `pulso.area_evento` midio sobre el, con la compuerta del lote.
    """
    T, V, _, _ = drs.leer_pulsos(p["corrida"].muestra, n=n, canal=p["canal"])
    d = p["areas"]
    amp_med = float(np.median(d["amplitud"][p["buenos"]]))
    on_lo, on_hi = np.percentile(d["t_onset"][p["buenos"]], [25, 75])
    mejor = None
    for t, v in zip(T, V):
        r = pulso.area_evento(t, v, pre=p["ventana"]["pre"], post=d["post"])
        if r["bandera"] != 0 or not (on_lo <= r["t_onset"] <= on_hi):
            continue
        costo = abs(r["amplitud"] - amp_med)
        if mejor is None or costo < mejor[0]:
            mejor = (costo, t, v, r)
    if mejor is None:
        t, v = T[0], V[0]
        mejor = (0.0, t, v, pulso.area_evento(t, v, pre=p["ventana"]["pre"],
                                              post=d["post"]))
    _, t, v, r = mejor
    tau, vv = pulso.preparar(t, v)
    return tau, vv, r


def figura_integracion(res):
    """
    Como un registro del DRS4 se convierte en un area, sobre la corrida de
    ejemplo. Cuatro paneles: (a) pulsos crudos tal como vienen, (b) UN pulso
    con cada paso de la integracion marcado, (c) el pulso promedio alineado
    por onset con la compuerta, (d) la carga acumulada y que fraccion entra.
    """
    f = corrida_ejemplo(res)
    p = pulsos_de(f)
    v, d = p["ventana"], p["areas"]
    pre, post = v["pre"], float(d["post"])
    guarda = _defaults(pulso.area_evento)["guarda_base"]
    chico = ESTILO["legend.fontsize"] - 1
    verde, naranja, rojo = "#2e8b57", "#e08a1e", "#c0392b"

    fig = figura(1.0, 0.78)
    ejes = fig.subplots(2, 2)
    a, b, c_, d_ = ejes[0][0], ejes[0][1], ejes[1][0], ejes[1][1]

    # ---- (a) pulsos crudos -----------------------------------------------
    T, V, _, _ = drs.leer_pulsos(p["corrida"].muestra, n=40, canal=p["canal"])
    tau = None
    for t, vv in zip(T[:30], V[:30]):
        tau, y = pulso.preparar(t, vv)
        a.plot(tau, y, lw=0.5, alpha=0.55, color="#1f4e79")
    on5, on50, on95 = p["onset"]
    a.axvspan(on5, on95, color=naranja, alpha=0.18, lw=0)
    a.axvline(on50, color=naranja, lw=1.0, ls="--")
    a.text(on95 + 3, 0.05, "onset: mediana %.0f ns,\np5–p95 de %.0f a %.0f ns"
           % (on50, on5, on95), transform=a.get_xaxis_transform(),
           fontsize=chico, color="#7a3d05", ha="left", va="bottom")
    a.set_xlim(0, float(tau[-1]))
    a.set_xlabel("tiempo desde el inicio del registro [ns]")
    a.set_ylabel("señal [mV]")

    # ---- (b) un pulso, paso a paso ----------------------------------------
    tau, y, r = pulso_tipico(p)
    t_on, base, amp = r["t_onset"], r["base"], r["amplitud"]
    i_ext = int(np.argmax(np.abs(y - base)))
    signo = -1.0 if (y[i_ext] - base) < 0 else 1.0
    lo, hi = r["t_ini"], r["t_fin"]
    x0, x1 = t_on - 28.0, t_on + 40.0
    zona = (tau >= x0) & (tau <= x1)
    b.plot(tau[zona], y[zona], ".-", ms=2.2, lw=0.7, color="#1f4e79",
           zorder=3)
    b.axvspan(x0, t_on - guarda, color=verde, alpha=0.15, lw=0)
    b.axhline(base, color=verde, lw=1.0)
    b.axhline(base + signo * pulso.FRAC_ONSET * amp, color="0.4", lw=0.8,
              ls=":")
    b.axvline(t_on, color=naranja, lw=1.0, ls="--")
    b.axvspan(lo, hi, color=naranja, alpha=0.18, lw=0)
    g = (tau >= lo) & (tau <= hi)
    b.fill_between(tau[g], y[g], base, color=rojo, alpha=0.35, lw=0)
    # las anotaciones van en la mitad del panel que el pulso deja libre: la
    # linea base a la izquierda, sobre la zona verde, y lo demas a la derecha
    # del pulso, que a 16 ns del onset ya volvio a la linea base
    b.text(x0 + 1.5, base + signo * 0.22 * amp,
           "línea base: mediana\nhasta %.0f ns antes\ndel onset" % guarda,
           fontsize=chico, color=verde, ha="left",
           va="top" if signo < 0 else "bottom")
    b.annotate("onset: cruce del\n%.0f %% de la amplitud"
               % (100 * pulso.FRAC_ONSET), xy=(t_on, base + signo * 0.1 * amp),
               xytext=(t_on + 16, base + signo * 0.32 * amp), fontsize=chico,
               color="#7a3d05", ha="left", va="center",
               arrowprops=dict(arrowstyle="->", color="#7a3d05", lw=0.8))
    b.text(t_on + 16, base + signo * 0.68 * amp,
           "compuerta [−%.0f, +%.0f] ns\nárea = %.0f mV·ns\namplitud = %.0f mV"
           % (pre, post, r["area"], amp), fontsize=chico, color=rojo,
           ha="left", va="center")
    b.set_xlim(x0, x1)
    b.set_xlabel("tiempo desde el inicio del registro [ns]")
    b.set_ylabel("señal [mV]")

    # ---- (c) el pulso promedio alineado ----------------------------------
    tau_r, forma = v["tau_rel"], v["forma"]
    c_.plot(tau_r, forma, lw=1.2, color="#1f4e79")
    c_.axhline(0, color="k", lw=0.5)
    c_.axvspan(-pre, post, color=naranja, alpha=0.18, lw=0)
    c_.axvline(-pre, color=verde, lw=1.0, ls="--")
    c_.axvline(post, color=rojo, lw=1.0, ls="--")
    c_.axvline(v["post"], color="0.4", lw=0.9, ls=":")
    c_.text(v["post"] + 2, 0.55, "el pulso se hunde\nen el ruido: %.0f ns" % v["post"],
            fontsize=chico, color="0.3", ha="left", va="center")
    c_.text(post + 2, 0.92, "compuerta:\n%.0f ns" % post, fontsize=chico,
            color=rojo, ha="left", va="top")
    c_.set_xlim(-25, float(tau_r[-1]))
    c_.set_ylim(-0.06, 1.05)
    c_.set_xlabel("tiempo desde el onset [ns]")
    c_.set_ylabel("pulso promedio (%d pulsos)" % v["n_usados"])

    # ---- (d) la carga acumulada -------------------------------------------
    _, frac, _ = pulso.curva_convergencia(tau_r, forma)
    d_.plot(tau_r, 100 * frac, lw=1.4, color="#1f4e79")
    d_.axhline(100, color="k", lw=0.6, ls=":")
    d_.axvline(-pre, color=verde, lw=1.0, ls="--")
    d_.axvline(post, color=rojo, lw=1.0, ls="--")
    d_.axvline(v["post"], color="0.4", lw=0.9, ls=":")
    d_.plot([post], [100 * p["carga"]], "o", ms=5, color=rojo)
    d_.annotate("%.1f %% de la carga\nen %.0f ns" % (100 * p["carga"], post),
                xy=(post, 100 * p["carga"]), xytext=(post + 14, 45),
                fontsize=chico, color=rojo, ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color=rojo, lw=0.8))
    d_.text(v["post"] + 2, 20, "%.1f %% en %.0f ns"
            % (100 * v["fraccion"], v["post"]),
            fontsize=chico, color="0.3", ha="left", va="center")
    d_.set_xlim(-25, float(tau_r[-1]))
    d_.set_ylim(-5, 112)
    d_.set_xlabel("tiempo desde el onset [ns]")
    d_.set_ylabel("carga acumulada [%]")

    for ax, letra in zip((a, b, c_, d_), "abcd"):
        etiqueta_panel(ax, "(%s)" % letra)
        ax.grid(alpha=.3)
        pocos_ticks(ax, 6, 5)
    return guardar(fig, "integracion_pulso")


def figura_pulsos_detectores(res):
    """
    El pulso promedio de cada detector a 1600 V con Cs-137, y su carga
    acumulada, con la compuerta comun marcada.

    La carga se normaliza a la acumulada hasta `T_COMUN` y no al "total",
    porque el total de cada corrida es lo que su REGISTRO contiene: en las
    sesiones de septiembre el trigger quedo tarde y el registro termina a
    47-80 ns del onset, asi que normalizar ahi hace que el D "converja" a los
    45 ns y el A a los 130 aunque la forma sea la misma. Cada curva se dibuja
    hasta donde su registro llega, y ahi se marca.
    """
    filas = [f for f in (por_clave(res, det, 1600, "Cs-137")
                         for det in detectores(res)) if f]
    if not filas:
        return None
    post = float(filas[0]["post"])
    pre = pulsos_de(filas[0])["ventana"]["pre"]
    chico = ESTILO["legend.fontsize"] - 1
    fig = figura(1.0, 0.44)
    a, b = fig.subplots(1, 2)
    manijas, etiquetas = [], []
    for f in filas:
        p = pulsos_de(f)
        v, Q = p["ventana"], p["Q"]
        tau = v["tau_rel"]
        col = COLOR.get(f["detector"], "0.4")
        h, = a.plot(tau, v["forma"], lw=1.1, color=col)
        q = 100 * Q / _Q_en(tau, Q, T_COMUN)
        m = tau <= p["espacio"]
        b.plot(tau[m], q[m], lw=1.1, color=col)
        b.plot([tau[m][-1]], [q[m][-1]], "|", ms=9, mew=1.4, color=col)
        b.plot([post], [100 * p["carga_comun"]], "o", ms=4, color=col)
        manijas.append(h)
        etiquetas.append("detector %s: %.0f %% en %.0f ns; registro hasta %.0f ns"
                         % (f["detector"], 100 * p["carga_comun"], post,
                            p["espacio"]))
    for ax in (a, b):
        ax.axvspan(-pre, post, color="#e08a1e", alpha=0.15, lw=0)
        ax.axvline(post, color="#c0392b", lw=1.0, ls="--")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlim(-6, 60)
    a.set_ylim(-0.06, 1.05)
    a.set_xlabel("tiempo desde el onset [ns]")
    a.set_ylabel("pulso promedio (amplitud 1)")
    b.axhline(100, color="k", lw=0.6, ls=":")
    b.axvline(T_COMUN, color="0.4", lw=0.8, ls=":")
    b.text(T_COMUN + 2, 30, "referencia:\ncarga a %.0f ns" % T_COMUN,
           fontsize=chico, color="0.3", ha="left", va="center")
    b.set_xlim(-6, 140)
    b.set_ylim(-5, 122)
    b.set_xlabel("tiempo desde el onset [ns]")
    b.set_ylabel("$Q(t)\\,/\\,Q(%.0f\\,\\mathrm{ns})$ [%%]" % T_COMUN)
    b.text(post + 2, 6, "compuerta\ncomún", fontsize=chico, color="#c0392b",
           ha="left", va="bottom")
    for ax, letra in zip((a, b), "ab"):
        etiqueta_panel(ax, "(%s)" % letra)
        ax.grid(alpha=.3)
        pocos_ticks(ax, 6, 5)
    leyenda_abajo(fig, manijas, etiquetas, ncol=2)
    return guardar(fig, "pulsos_detectores")


def figura_area_a_espectro(res):
    """
    De las areas al espectro, sobre la corrida de ejemplo. (a) area contra
    amplitud, evento por evento, con la mediana de area/amplitud por tramo:
    si la forma del pulso no dependiera del tamano, el cociente seria una
    constante. (b) los tres histogramas -muestra, fondo normalizado y neto-
    con el binning comun, el fin de la subida del trigger y la ventana de
    ajuste.
    """
    f = corrida_ejemplo(res)
    p = pulsos_de(f)
    d, ok = p["areas"], p["buenos"]
    chico = ESTILO["legend.fontsize"] - 1
    fig = figura(1.0, 0.44)
    a, b = fig.subplots(1, 2)

    # ---- (a) area / amplitud ---------------------------------------------
    A, ar = d["amplitud"][ok], d["area"][ok]
    m = (A > 10.0) & (A < np.percentile(A, 99.7))
    q = ar[m] / A[m]
    y_lo, y_hi = np.percentile(q, [0.5, 99.5])
    a.hist2d(A[m], q, bins=[60, 60],
             range=[[float(A[m].min()), float(A[m].max())], [y_lo, y_hi]],
             cmap="Blues", cmin=1, rasterized=True)
    cx, cy, ce = p["perfil"]
    a.errorbar(cx, cy, yerr=ce, fmt="o-", ms=3.5, lw=1.0, color="#c0392b",
               capsize=2)
    a.text(0.97, 0.05, "mediana por tramo de amplitud:\nvaría %.1f %% de "
           "punta a punta" % (100 * p["variacion_forma"]),
           transform=a.transAxes, fontsize=chico, color="#c0392b", ha="right",
           va="bottom",
           bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=.85))
    a.set_xlabel("amplitud del pulso [mV]")
    a.set_ylabel("área / amplitud [ns]")

    # ---- (b) los histogramas ---------------------------------------------
    e = espectro.cargar("espectro_%s_%s" % (f["clave"], f.get("metodo", "global")))
    umbral = espectro.fin_del_umbral(e["canales"], e["neto"])
    x_max = float(e["bordes"][-1])
    # la muestra va ancha y clara DEBAJO del neto: donde el fondo pesa poco
    # las dos curvas coinciden y una tapa a la otra
    manijas = [b.step(e["canales"], e["muestra"], where="mid", lw=2.2,
                      color="#9dbbd8",
                      label="muestra (%s pulsos)"
                      % format(e["n_muestra"], ",d"))[0]]
    if e["hay_fondo"]:
        manijas.append(b.step(e["canales"], e["fondo"], where="mid", lw=0.9,
                              color="0.45",
                              label="fondo × %.2g (%s pulsos)"
                              % (e["k_f"], format(e["n_fondo"], ",d")))[0])
    manijas.append(b.step(e["canales"], np.maximum(e["neto"], 0.0),
                          where="mid", lw=0.9, color="#2e8b57",
                          label="neto")[0])
    lo, hi = f["ventana"]
    manijas.append(b.axvspan(lo, hi, color="#e08a1e", alpha=0.15, lw=0,
                             label="ventana de ajuste"))
    b.axvline(umbral, color="0.3", lw=0.9, ls=":")
    b.text(umbral + 0.01 * x_max, 0.60, "fin de la subida\ndel trigger",
           fontsize=chico, color="0.3", ha="left", va="top",
           transform=b.get_xaxis_transform())
    b.set_yscale("log")
    b.set_xlim(0, x_max)
    y_max = float(np.max(e["muestra"][e["canales"] > umbral]))
    # el piso deja ver el fondo normalizado aunque k_f sea de 1e-3: que sea
    # despreciable es parte de lo que el panel muestra
    piso = 0.5
    if e["hay_fondo"] and np.any(e["fondo"] > 0):
        piso = min(piso, 0.5 * float(np.percentile(e["fondo"][e["fondo"] > 0],
                                                   20)))
    b.set_ylim(max(piso, 1e-3), 4.0 * y_max)
    b.set_xlabel("área del pulso [mV·ns]")
    b.set_ylabel("cuentas por bin de %.1f mV·ns" % (e["rango_max"] / e["bins"]))

    for ax, letra in zip((a, b), "ab"):
        etiqueta_panel(ax, "(%s)" % letra)
        ax.grid(alpha=.3)
    pocos_ticks(a, 6, 5)
    # la leyenda del panel (b) va debajo: adentro, en cualquier rincon, tapa
    # o el plateau o el flanco o la cola
    leyenda_abajo(fig, manijas, [m.get_label() for m in manijas], ncol=2)
    return guardar(fig, "area_a_espectro")


def resumen_pulsos(res):
    """
    Lo que la seccion de pulsos cita sobre la campana entera: los pulsos de
    cada corrida (via `pulsos_de`), los metadatos del registro de la corrida
    de ejemplo y el paso temporal real de un archivo binario.
    """
    filas = ct.ordenadas(res)
    por = {f["clave"]: pulsos_de(f) for f in filas}
    f_ej = corrida_ejemplo(res)
    meta = drs.metadatos(rutas.corrida(f_ej["clave"]).muestra)
    # el paso entre celdas, en un archivo que traiga la calibracion por celda
    # (los binarios la traen; el XML del 26-ago viene con eje uniforme)
    clave_dt = next((f["clave"] for f in filas
                     if por[f["clave"]]["formato"] == "binario"), f_ej["clave"])
    T, _, _, _ = drs.leer_pulsos(rutas.corrida(clave_dt).muestra, n=1)
    return dict(por=por, meta=meta, ejemplo=por[f_ej["clave"]],
                dt=np.diff(np.asarray(T[0], float)), clave_dt=clave_dt)


# --------------------------------------------------------------------------
# 5. Las tablas
# --------------------------------------------------------------------------
def tabla_calibraciones(res):
    filas, det_previo = [], None
    for f in ct.ordenadas(res):
        if det_previo is not None and f["detector"] != det_previo:
            filas.append(None)
        det_previo = f["detector"]
        err = ("%.4f" % f["err_D"]) if np.isfinite(f["err_D"] or np.nan) \
            else r"\textit{nd}"
        filas.append([
            f["detector"], "%g" % f["hv"],
            tex(f["fuente"] + (" (descartada)" if f["descartada"] else "")),
            tex(f["fecha"]),
            "$%.4f \\pm %s$" % (f["D"], err),
            "%.1f" % f["inv_D"],
            "%.2f" % f["chi2_ndf"],
            "%.1f" % (100 * f["fwhm_rel"]),
            "%.0f" % (100 * f["fraccion_carga"]),
            "si" if f["hay_fondo"] else "no"])
    # Los encabezados van sin unidades a proposito: con diez columnas, el ancho
    # de la tabla lo fijan los rotulos y no los numeros, y "$1/D$ [keV/mV·ns]"
    # es tres veces mas ancho que cualquier dato de esa columna. Las unidades
    # van en la nota.
    return tabla(
        ["Det.", "HV", "Fuente", "Medida", r"$D$", r"$1/D$",
         r"$\chi^2$/ndf", "FWHM", "Carga", "Fondo"],
        filas, "llllrrrrrc", largo=True,
        nota="Unidades: HV en V, $D$ en mV$\\cdot$ns/keV, $1/D$ en "
             "keV/mV$\\cdot$ns, FWHM y carga en \\%. $1/D$ es lo que se usa "
             "para convertir: energia depositada por unidad de area. FWHM "
             "esta evaluada en el borde de referencia de cada fuente "
             "(477.3 keV el Cs-137, 963.4 keV el Co-60). «Carga» es la "
             "fraccion del pulso registrado que entra en la compuerta comun. "
             "«Fondo» dice si la corrida tiene fondo propio restado.")


def tabla_pulsos(res, pul):
    """
    Una fila por corrida con lo que paso ANTES del ajuste: cuantos eventos
    trae el archivo, cuantos sobreviven a las banderas, cuanto tiempo estuvo
    viva la adquisicion, hasta donde llega el registro despues del onset y
    que fraccion de la carga entra en la compuerta, en las dos normalizaciones
    (la del pulso registrado, que es la columna Carga de la tabla de
    calibraciones, y la de la carga a T_COMUN, que es comparable entre
    corridas).
    """
    filas, det_previo = [], None
    for f in ct.ordenadas(res):
        p = pul["por"][f["clave"]]
        if det_previo is not None and f["detector"] != det_previo:
            filas.append(None)
        det_previo = f["detector"]
        filas.append([
            f["detector"], "%g" % f["hv"],
            tex(f["fuente"] + (" (desc.)" if f["descartada"] else "")),
            "XML" if p["formato"] == "xml" else "bin.",
            format(p["n_eventos"], ",d").replace(",", "\\,"),
            "%.1f" % (100 * p["descartados"]),
            "%.0f" % f["t_pared"],
            "%.2f" % (100 * f["frac_viva"]),
            "%.0f" % p["espacio"],
            "%.0f" % (100 * p["carga"]),
            "%.0f" % (100 * p["carga_comun"])])
    return tabla(
        ["Det.", "HV", "Fuente", "Form.", "Eventos", "Desc.", "Pared",
         "Vivo", "Registro", "Carga",
         r"$Q_{%.0f}/Q_{%.0f}$" % (ct.POST_FIJO, T_COMUN)],
        filas, "lllcrrrrrrr", largo=True,
        nota="HV en V. «Form.»: formato del archivo (XML del DRSOSC o binario "
             "del DRS4). «Eventos»: registros en el archivo. «Desc.»: "
             "porcentaje descartado por las banderas de calidad (saturacion, "
             "compuerta fuera del registro, sin pulso). «Pared»: duracion de "
             "la corrida en s; «Vivo»: fraccion del tiempo de pared en que la "
             "adquisicion estuvo viva, en \\%%, medida con el scaler de la "
             "placa. «Registro»: cuanto registro queda despues del onset, en "
             "ns (mediana). «Carga»: fraccion de la carga del pulso "
             "\\emph{registrado} que entra en la compuerta de %g ns, en "
             "\\%% (la misma columna de la tabla~\\ref{tab:calib}); "
             "$Q_{%.0f}/Q_{%.0f}$: la misma compuerta, relativa a la carga "
             "acumulada hasta %.0f ns del onset, que todos los registros "
             "contienen." % (ct.POST_FIJO, ct.POST_FIJO, T_COMUN, T_COMUN))


def tabla_test_a(res, ta):
    filas = []
    for (det, hv), d in sorted(ta.items()):
        if d["lin"] is None:
            # 1 columna + un multicolumn de 7 = las 8 columnas de la tabla.
            # Ojo con rellenar con celdas vacias despues de un multicolumn:
            # es la forma tipica de romper un tabular sin que el error de
            # LaTeX diga donde.
            filas.append([r"%s @ %g V" % (det, hv),
                          r"\multicolumn{7}{l}{\textit{%s}}"
                          % tex("una sola fuente: el test no se puede correr")])
            continue
        cs = por_clave(res, det, hv, "Cs-137")
        co = por_clave(res, det, hv, "Co-60")
        q = d["cociente"]
        filas.append([
            "%s @ %g V" % (det, hv),
            "%.4f" % cs["D"], "%.4f" % co["D"],
            "%.4f" % d["lin"]["D_origen"],
            "%.3f" % q["valor"],
            "$%+.1f$" % (100 * (q["corregido"] - 1)),
            "$%+.0f \\pm %.0f$" % (d["lin"]["ordenada"],
                                   d["lin"]["err_ordenada"]),
            "misma" if d["misma_sesion"] else "distintas"])
    return tabla(
        ["Punto", r"$D$(Cs-137)", r"$D$(Co-60)", r"$D$ conjunta",
         "Cs/Co", "Desac. [\\%]", r"$c_0$ [mV$\cdot$ns]", "Sesiones"],
        filas, "lrrrrrrc",
        nota="«Desac.» es el desacuerdo entre las dos fuentes una vez "
             "descontado el quenching de Birks: en escala de energia "
             "depositada una cadena perfectamente lineal da "
             "Cs/Co $= 0.980$, no 1. $c_0$ es la ordenada de un ajuste de dos "
             "puntos y dos parametros, asi que pasa exacto: lo unico que dice "
             "algo es su error.")


def tabla_test_b(res, tb):
    filas = []
    for (det, fuente), d in sorted(tb.items()):
        aj = d["ajuste"]
        if aj is None:
            filas.append([det, tex(fuente), "%d" % len(d["filas"]),
                          r"\textit{%s}" % tex(d["motivo"]), "", "", ""])
            continue
        peor = 100 * np.max(np.abs(aj["residuos_potencia"] / aj["D"]))
        filas.append([det, tex(fuente), "%d" % aj["n"],
                      "$%.2f \\pm %.2f$" % (aj["b"], aj["err_b"]),
                      "%.0f" % aj["chi2_potencia"],
                      "%.0f" % aj["chi2_recta"],
                      "%.1f" % peor])
    return tabla(
        ["Det.", "Fuente", r"$n_V$", r"$b$", r"$\chi^2$ potencia",
         r"$\chi^2$ recta", "Peor residuo [\\%]"],
        filas, "llrrrrr",
        nota="Ajuste de $D(V)$ con las dos hipotesis. El $\\chi^2$ de la recta "
             "es una a dos ordenes de magnitud mayor en todos los casos, o sea "
             "que la ley de potencia gana sin ambiguedad. Un PMT tipico tiene "
             "$b$ entre 6 y 10.")


def tabla_detectores(res, tb):
    """La comparacion entre detectores, con la traduccion a voltaje."""
    hvs = voltajes(res)
    filas = []
    for fuente in ("Cs-137", "Co-60"):
        hay = [d for d in detectores(res)
               if any(por_clave(res, d, v, fuente) for v in hvs)]
        if not hay:
            continue
        filas.append(None)
        filas.append([r"\multicolumn{%d}{l}{\textbf{%s}}"
                      % (len(hvs) + 3, tex(fuente))])
        for det in hay:
            fila = [det]
            for hv in hvs:
                f = por_clave(res, det, hv, fuente)
                fila.append("%.4f" % f["D"] if f else "--")
            ref = por_clave(res, "A", 1600, fuente)
            f16 = por_clave(res, det, 1600, fuente)
            aj = (tb.get((det, fuente)) or {}).get("ajuste")
            if f16 and ref:
                fila.append("%.3f" % (f16["D"] / ref["D"]))
                v = igualar_voltaje(ref, f16, aj["b"] if aj else None)
                fila.append("--" if not np.isfinite(v) else "%.0f" % v)
            else:
                fila += ["--", "--"]
            filas.append(fila)
    return tabla(
        ["Det."] + ["%g V" % v for v in hvs] + [r"$D/D_A$", r"$V_{eq}$"],
        filas, "l" + "r" * (len(hvs) + 2),
        nota=r"Las columnas de voltaje son $D$ en mV$\cdot$ns/keV a ese alto "
             "voltaje; $D/D_A$ compara con el detector A a 1600 V. "
             "$V_{eq}$ responde la pregunta practica: a que voltaje "
             "hay que poner cada detector para que tenga la ganancia que el "
             "detector A tiene a 1600 V, usando su propia ley "
             "$D \\propto V^{b}$. Extrapolar fuera de 1600--1700 V es "
             "justamente eso, extrapolar.")


def tabla_descartadas(res):
    filas = []
    for a, b in ct.pares_descartadas(res):
        eb = ("%.4f" % b["err_D"]) if np.isfinite(b["err_D"] or np.nan) \
            else r"\textit{nd}"
        filas.append(["%s @ %g V" % (a["detector"], a["hv"]), tex(a["fuente"]),
                      "$%.4f \\pm %.4f$" % (a["D"], a["err_D"]),
                      tex(a["fecha"]),
                      "$%.4f \\pm %s$" % (b["D"], eb), tex(b["fecha"]),
                      "$%+.1f$" % (100 * (b["D"] - a["D"]) / a["D"])])
    return tabla(
        ["Punto", "Fuente", r"$D$ (sesion valida)", "Medida",
         r"$D$ (tierra mal)", "Medida", "Cambio [\\%]"],
        filas, "llrlrlr",
        nota="Mismo detector, mismo voltaje nominal y misma fuente. La sesion "
             "del 1-sep se tomo con la tierra del detector mal conectada y "
             "esta descartada; la del 8-sep es la valida. La columna del "
             "cambio mide, entonces, cuanto movio la ganancia esa falla, no la "
             "reproducibilidad del montaje.")


def tabla_metricas_paper(res, ta):
    """Nuestros puntos con las metricas de validacion del paper."""
    filas = [[tex(f["punto"]), "$%+.1f$" % f["desacuerdo"], "%.1f" % f["sd_keV"],
              "%.8f" % f["coseno"], "%.4f" % f["r2"]]
             for f in metricas_paper(res, ta)]
    return tabla(
        ["Punto", "Desacuerdo [\\%]", "SD [keV]", "Similitud coseno",
         r"$R^2$"],
        filas, "lrrrr",
        nota="Las dos ultimas columnas son las metricas con las que el paper "
             "valida su calibracion (su Fig. 3). La similitud coseno da 1 "
             "EXACTAMENTE en todos los puntos, incluido el que discrepa mas: "
             "entre dos rectas por el origen es 1 sean cuales sean las "
             "pendientes, asi que no puede ver este modo de falla. El $R^2$ "
             "se calcula sobre dos puntos y con la recta forzada por el "
             "origen; con ordenada libre y dos puntos daria 1 exacto.")


def tabla_robustez(rob, p):
    """Nuestra dispersion frente al intervalo, contra la del paper."""
    filas = []
    for clave, d in sorted(rob.items(), key=lambda kv: kv[1]["detector"]):
        filas.append(["%s @ %g V" % (d["detector"], d["hv"]), tex(d["fuente"]),
                      "%d de %d" % (sum(1 for _, _, x in d["filas"]
                                        if np.isfinite(x)), len(d["filas"])),
                      "%.1f" % (100 * d["dispersion"]),
                      "--" if not np.isfinite(d["dispersion_sin_angosta"])
                      else "%.1f" % (100 * d["dispersion_sin_angosta"])])
    if p:
        filas.append(None)
        filas.append([r"\textit{paper, Fig. 5}", r"\textit{Cs-137}",
                      "%d de %d" % (6, 6),
                      "%.1f" % (100 * p["fig5_todas"]),
                      "%.1f" % (100 * p["fig5_sin_angosta"])])
    return tabla(
        ["Punto", "Fuente", "Intervalos", "Dispersion [\\%]",
         "Sin el mas angosto [\\%]"],
        filas, "llrrr",
        nota="El mismo test de la Fig. 5 del paper: la misma corrida ajustada "
             "en seis intervalos anidados, definidos en unidades del canal del "
             "borde y tomados de sus propios paneles. Un intervalo con menos "
             "de 30 bines se saltea, igual que en la replica. La fila del "
             "paper sale de los seis valores que rotula su figura.")


def tabla_control(res):
    ctrl = res.get("_control_compuerta") or {}
    if not ctrl:
        return ("\\textit{Este control no esta calculado. Se genera con "
                "\\texttt{python analisis/calibrar\\_todo.py -{}-resumen "
                "-{}-control}.}")
    filas = [[tex(p), "%.4f" % d["cociente_36"], "%.4f" % d["cociente"],
              "$%+.2f$" % d["cambio"]] for p, d in sorted(ctrl.items())]
    return tabla(
        ["Punto", "Cs/Co con %g ns" % ct.POST_FIJO, "Cs/Co con 75 ns",
         "Cambio [\\%]"],
        filas, "lrrr",
        nota="Donde el registro alcanza para las dos compuertas, el cociente "
             "entre fuentes se mueve menos del 1\\%: la eleccion de la "
             "compuerta no explica el desacuerdo del test A. En el detector D "
             "el control no se puede correr, porque con 75 ns entra el 2\\% de "
             "los eventos.")


# --------------------------------------------------------------------------
# 6. Los numeros que cita el texto
# --------------------------------------------------------------------------
def numeros(res, ta, tb, rob=None, pap=None, pul=None):
    filas = ct.ordenadas(res)
    bs = [d["ajuste"]["b"] for d in tb.values() if d["ajuste"]]
    desac = [100 * abs(d["cociente"]["corregido"] - 1)
             for d in ta.values() if d["cociente"]]
    err_rel = [100 * f["err_D"] / f["D"] for f in filas
               if np.isfinite(f["err_D"] or np.nan)]
    rep = [100 * (b["D"] - a["D"]) / a["D"] for a, b in ct.pares_descartadas(res)]
    ctrl = res.get("_control_compuerta") or {}
    puntos_ta = [d for d in ta.values() if d["lin"] is not None]

    L = [r"% Numeros del informe. Los escribe analisis/informe.py: no editar.",
         r"% Cada uno se usa en el texto como \NombreDeLaMacro.",
         macro("NCorridas", "%d" % len(filas)),
         macro("NDetectores", "%d" % len(detectores(res))),
         macro("NPuntosTestA", "%d" % len(puntos_ta)),
         macro("NGruposTestB", "%d" % len(bs)),
         macro("Rama", tex(ct.RAMA)),
         macro("PostFijo", "%g" % ct.POST_FIJO),
         macro("CocEsperado", "%.3f" % ct.cociente_esperado()),
         macro("QuenchCo", "%.4f" % ct.quenching("Co-60")),
         macro("QuenchCs", "%.4f" % ct.quenching("Cs-137")),
         macro("BMin", "%.2f" % min(bs)),
         macro("BMax", "%.2f" % max(bs)),
         macro("DesacMin", "%.1f" % min(desac)),
         macro("DesacMax", "%.1f" % max(desac)),
         macro("ErrRelMediano", "%.2f" % float(np.median(err_rel))),
         macro("RepMin", "%.0f" % min(abs(r) for r in rep)),
         macro("RepMax", "%.0f" % max(abs(r) for r in rep)),
         macro("ControlMax", "%.2f"
               % (max(abs(d["cambio"]) for d in ctrl.values()) if ctrl else 0)),
         macro("ChiMin", "%.2f" % min(f["chi2_ndf"] for f in filas)),
         macro("ChiMax", "%.2f" % max(f["chi2_ndf"] for f in filas)),
         macro("FwhmMin", "%.0f" % (100 * min(f["fwhm_rel"] for f in filas))),
         macro("FwhmMax", "%.0f" % (100 * max(f["fwhm_rel"] for f in filas))),
         # cuanta carga del pulso entra en la compuerta: es el precio de que
         # sea corta, y el texto lo cita
         macro("CargaMin", "%.0f"
               % (100 * min(f["fraccion_carga"] for f in filas))),
         macro("CargaMax", "%.0f"
               % (100 * max(f["fraccion_carga"] for f in filas))),
         macro("NSinFondo", "%d" % sum(1 for f in filas if not f["hay_fondo"])),
         # el NCC: cuanto correlaciona la ventana ganadora y cuanto se corre el
         # ajuste respecto de la semilla que esa ventana entrega
         macro("RhoMin", "%.3f" % min(f["rho"] for f in filas)),
         macro("RhoMax", "%.3f" % max(f["rho"] for f in filas)),
         # --- constantes fisicas y de configuracion, para que el texto no las
         # tipee: si manana cambia una linea gamma del catalogo, cambia sola
         macro("MasaElectron", "%.3f" % fisica.M_E),
         macro("EgammaCs", "%.3f" % fisica.FUENTES["Cs-137"][0]),
         macro("EgammaCoA", "%.3f" % fisica.FUENTES["Co-60"][0]),
         macro("EgammaCoB", "%.3f" % fisica.FUENTES["Co-60"][1]),
         macro("BordeCs", "%.1f" % fisica.bordes("Cs-137")[0]),
         macro("BordeCoA", "%.1f" % fisica.bordes("Co-60")[0]),
         macro("BordeCoB", "%.1f" % fisica.bordes("Co-60")[1]),
         macro("NArranques", "%d" % ct.OPCIONES.get(
             "n_arranques", ct.pipeline.DEFAULTS["n_arranques"])),
         macro("NBins", "%d" % ct.OPCIONES.get(
             "bins", ct.pipeline.DEFAULTS["bins"])),
         # cuanto aporta el optimizador a la incerteza: nada, y conviene poder
         # decirlo con el numero al lado
         macro("DispersionMax", "%.2f" % max(
             100 * f["dispersion_multiarranque"] / f["D"] for f in filas)),
         # con el modelo de dos parametros la mediana da cero a cuatro
         # decimales; escribir "0.0000 %" en el texto se lee como un redondeo
         # sospechoso, asi que por debajo del umbral se dice que esta debajo
         macro("DispersionMediana", (lambda v: "menos de 0.001"
                                     if v < 0.001 else "%.3f" % v)(
             float(np.median([100 * f["dispersion_multiarranque"] / f["D"]
                              for f in filas])))),
         macro("CorrMin", "%+.2f" % min(f["corr_D_sigma"] for f in filas)),
         macro("CorrMax", "%+.2f" % max(f["corr_D_sigma"] for f in filas)),
         macro("NFigurasDiag", "%d" % sum(len(f.get("figuras") or ())
                                          for f in filas)),
         # cuentas por bin adentro de la ventana: es lo que decide si el sesgo
         # de usar las barras del dato (chi2 de Neyman) importa o no
         macro("CuentasBin", "%.0f" % float(np.median(
             [cuentas_por_bin(f) for f in filas]))),
         # el arranque de la ventana de ajuste sale de la configuracion del
         # lote, no de las filas: es el mismo para las 25 corridas
         macro("InicioRel", "%.2f" % ct.OPCIONES.get(
             "inicio_rel", ct.pipeline.DEFAULTS["inicio_rel"])),
         macro("FechaInforme", time.strftime("%d-%m-%Y")),
         ]
    # ---- lo que el NCC entrega y como lo usa el recocido -----------------
    # Los defaults se leen de las firmas y de las constantes del modulo, no se
    # copian: si alguien cambia el barrido o las cotas, el texto lo sigue.
    import inspect
    from analisis import ncc as mod_ncc
    d_ncc = {k: v.default for k, v in
             inspect.signature(mod_ncc.buscar_intervalo).parameters.items()}
    d_sa = {k: v.default for k, v in
            inspect.signature(sa_mod.recocido).parameters.items()}
    L += [macro("SigmaRelNCC", "%.2f" % d_ncc["sigma_rel"]),
          macro("NBarrido", "%d" % d_ncc["n_B"]),
          macro("CotaDMin", "%.1f" % calibracion.COTA_D_REL[0]),
          macro("CotaDMax", "%.1f" % calibracion.COTA_D_REL[1]),
          macro("CotaSigmaMin", "%.2f" % calibracion.COTA_SIGMA_REL[0]),
          macro("CotaSigmaMax", "%.2f" % calibracion.COTA_SIGMA_REL[1]),
          macro("CotaLambdaMin", "%.0f" % calibracion.COTA_LAMBDA_MIN),
          macro("LambdaSemillaRel", "%.0f" % (100 * calibracion.LAMBDA_SEMILLA_REL)),
          macro("InicioAjusteRel", "%.2f" % calibracion.INICIO_AJUSTE_REL),
          macro("NMarkov", "%d" % d_sa["n_markov"]),
          macro("Alfa", "%.2f" % d_sa["alfa"])]
    # la corrida de ejemplo, con lo que el barrido le entrego al recocido
    f_ej = corrida_ejemplo(res)
    r_ej = resultado_guardado(f_ej)
    if r_ej is not None and r_ej.get("semilla_D"):
        D0 = float(r_ej["semilla_D"])
        vn = r_ej.get("ventana_ncc", r_ej["ventana"])
        L += [macro("EjCorrida", tex("detector %s a %g V con %s"
                                     % (f_ej["detector"], f_ej["hv"],
                                        f_ej["fuente"]))),
              macro("EjB", "%.0f" % (D0 * r_ej["E_ref"])),
              macro("EjRho", "%.4f" % r_ej["rho"]),
              macro("EjVentanaLo", "%.0f" % vn[0]),
              macro("EjVentanaHi", "%.0f" % vn[1]),
              macro("EjSemillaD", "%.4f" % D0),
              macro("EjCotaDMin", "%.3f" % (calibracion.COTA_D_REL[0] * D0)),
              macro("EjCotaDMax", "%.3f" % (calibracion.COTA_D_REL[1] * D0)),
              macro("EjD", "%.4f" % r_ej["D"]),
              macro("EjDesvSemilla", "%+.1f" % (100 * (r_ej["D"] / D0 - 1)))]
    sem = semillas(res)
    if sem:
        desvio = [abs(D / s - 1) for s, D in sem]
        L += [macro("NSemillas", "%d" % len(sem)),
              macro("SemillaMediana", "%.2f" % (100 * float(np.median(desvio)))),
              macro("SemillaMax", "%.1f" % (100 * max(desvio)))]

    if pul:
        L += macros_pulsos(res, pul)

    # la ganancia de cada detector a 1600 V con Cs-137, que es el punto que
    # todos los detectores tienen medido
    for det in detectores(res):
        f = por_clave(res, det, 1600, "Cs-137")
        if f:
            L.append(macro("Dmilseis%s" % det, "%.4f" % f["D"]))
    razon = [por_clave(res, d, 1600, "Cs-137") for d in detectores(res)]
    razon = [f["D"] for f in razon if f]
    if razon:
        L.append(macro("RazonGanancias", "%.1f" % (max(razon) / min(razon))))

    # El ejemplo de la seccion "como usar estos numeros": a cuanta energia
    # corresponde un pulso de 1000 mV*ns en un caso concreto. Va como macro y
    # no escrito en el parrafo por lo mismo que todo el resto.
    ref = por_clave(res, "A", 1600, "Cs-137")
    if ref:
        L.append(macro("EjemploKeV", "%.0f" % (1000.0 / ref["D"])))

    # --- la comparacion con el paper --------------------------------------
    met = metricas_paper(res, ta)
    if met:
        L += [macro("SDNuestroMin", "%.0f" % min(f["sd_keV"] for f in met)),
              macro("SDNuestroMax", "%.0f" % max(f["sd_keV"] for f in met)),
              macro("CosenoNuestro",
                    "%.6f" % max(f["coseno"] for f in met)),
              macro("RdosNuestroMin", "%.3f" % min(f["r2"] for f in met))]
    if pap:
        L += [macro("SDPaperMin", "%.1f" % pap["sd_min"]),
              macro("SDPaperMax", "%.1f" % pap["sd_max"]),
              macro("RdosPaperMin", "%.6f" % pap["r2_min"]),
              macro("RdosPaperMax", "%.6f" % pap["r2_max"]),
              macro("CosenoPaperMin", "%.6f" % pap["coseno_min"]),
              macro("CosenoPaperMax", "%.6f" % pap["coseno_max"]),
              macro("RobustezPaper", "%.1f" % (100 * pap["fig5_todas"])),
              macro("RobustezPaperSinAngosta",
                    "%.1f" % (100 * pap["fig5_sin_angosta"])),
              macro("ReduccionDif", "%d" % pap["reducciones"]["BC408"]["diferencia"]),
              macro("ReduccionGauss", "%d" % pap["reducciones"]["BC408"]["gauss"])]
    if rob:
        ds = [100 * d["dispersion_sin_angosta"] for d in rob.values()
              if np.isfinite(d["dispersion_sin_angosta"])]
        if ds:
            L += [macro("RobustezNuestraMin", "%.1f" % min(ds)),
                  macro("RobustezNuestraMax", "%.1f" % max(ds))]
    # cuanto sube la ganancia al subir 50 V, con el b mediano de la campana
    bs_ = [d["ajuste"]["b"] for d in tb.values() if d["ajuste"]]
    if bs_:
        b_med = float(np.median(bs_))
        L.append(macro("SubidaCincuenta",
                       "%.0f" % (100 * ((1650.0 / 1600.0) ** b_med - 1))))
    return "\n".join(L)


def _kf_pared(e):
    """El factor del fondo si se normalizara por tiempo de pared: la misma
    regla de `espectro.espectro_neto` (al menor de los dos tiempos) con los
    tiempos de reloj en lugar de los vivos."""
    if not e["hay_fondo"] or e["t_pared_fondo"] <= 0:
        return 0.0
    t_norm = min(e["t_pared_muestra"], e["t_pared_fondo"])
    return float(t_norm / e["t_pared_fondo"])


def macros_pulsos(res, pul):
    """
    Los numeros de la seccion "Del pulso al espectro". Las constantes del
    analisis se leen de los modulos y de las firmas (`_defaults`), no se
    copian; lo medido sale de `resumen_pulsos`.
    """
    filas = ct.ordenadas(res)
    por, ej, meta = pul["por"], pul["ejemplo"], pul["meta"]
    d_ev = _defaults(pulso.area_evento)
    d_lb = _defaults(pulso.linea_base)
    d_cv = _defaults(area.calibrar_ventana)
    d_pp = _defaults(pulso.pulso_promedio)
    dt = pul["dt"]
    dt = dt[dt > 0]
    f_ej = corrida_ejemplo(res)
    e_ej = espectro.cargar("espectro_%s_%s"
                           % (f_ej["clave"], f_ej.get("metodo", "global")))
    umbral = espectro.fin_del_umbral(e_ej["canales"], e_ej["neto"])

    def pct(x):
        return "%.1f" % (100 * x)

    def rango(vals, fmt):
        return fmt % min(vals), fmt % max(vals)

    viva_cs = [f["frac_viva"] for f in filas if f["fuente"] == "Cs-137"]
    viva_co = [f["frac_viva"] for f in filas if f["fuente"] == "Co-60"]
    pres = [por[f["clave"]]["ventana"]["pre"] for f in filas]
    espacios = [por[f["clave"]]["espacio"] for f in filas]
    cargas_comun = [por[f["clave"]]["carga_comun"] for f in filas]
    variaciones = [por[f["clave"]]["variacion_forma"] for f in filas]
    var_altas = [por[f["clave"]]["variacion_alta"] for f in filas]
    marcados = [por[f["clave"]]["descartados"] for f in filas]
    sat = [por[f["clave"]]["marcados"]["saturados"] for f in filas]
    trunc = [por[f["clave"]]["marcados"]["truncados"] for f in filas]
    sinpre = [por[f["clave"]]["marcados"]["sin_pre"] for f in filas]
    tau_r, forma = ej["ventana"]["tau_rel"], ej["ventana"]["forma"]
    i_pico = int(np.argmax(forma))
    arriba = np.flatnonzero(forma > 0.5)
    _, fr, _ = pulso.curva_convergencia(tau_r, forma)

    L = [r"% --- la seccion 'Del pulso al espectro' ---",
         # el registro
         macro("Muestras", "%d" % meta["muestras"]),
         macro("Registro", "%.0f" % meta["duracion_registro"]),
         macro("Muestreo", "%.2f" % meta["muestreo_GSs"]),
         macro("DtMedio", "%.3f" % meta["dt_medio"]),
         macro("DtMin", "%.2f" % np.percentile(dt, 1)),
         macro("DtMax", "%.2f" % np.percentile(dt, 99)),
         macro("EventosTipico", format(int(np.median(
             [por[f["clave"]]["n_eventos"] for f in filas])), ",d")
             .replace(",", r"\,")),
         # las constantes de la integracion
         macro("RecorteBordes", "%d" % pulso.RECORTE_BORDES),
         macro("UmbralRuido", "%g" % area.UMBRAL_RUIDO),
         macro("FracOnset", "%.0f" % (100 * pulso.FRAC_ONSET)),
         macro("GuardaBase", "%g" % d_ev["guarda_base"]),
         macro("MinBase", "%d" % d_lb["minimo"]),
         macro("VSatura", "%.0f" % pulso.V_SATURA),
         macro("NVentana", "%d" % d_cv["n"]),
         macro("AmpMinimaForma", "%.0f" % d_pp["amp_minima"]),
         macro("KConv", "%g" % d_cv["k"]),
         macro("TramoConv", "%g" % d_cv["tramo"]),
         macro("MargenRegistro", "%g" % d_cv["margen_registro"]),
         macro("PreMin", "%.1f" % min(pres)),
         macro("PreMax", "%.1f" % max(pres)),
         macro("TComun", "%.0f" % T_COMUN),
         # la corrida de ejemplo
         macro("EjEventos", format(ej["n_eventos"], ",d").replace(",", r"\,")),
         macro("EjBuenos", format(ej["n_buenos"], ",d").replace(",", r"\,")),
         macro("EjDescartados", "%.2f" % (100 * ej["descartados"])),
         macro("EjOnsetMediano", "%.0f" % ej["onset"][1]),
         macro("EjOnsetCinco", "%.0f" % ej["onset"][0]),
         macro("EjOnsetNoventaycinco", "%.0f" % ej["onset"][2]),
         macro("EjAmpMediana", "%.0f" % ej["amplitud"][0]),
         macro("EjAmpNoventaynueve", "%.0f" % ej["amplitud"][1]),
         macro("EjAreaMediana", "%.0f" % ej["area"][0]),
         macro("EjNForma", "%d" % ej["ventana"]["n_usados"]),
         macro("EjPico", "%.1f" % tau_r[i_pico]),
         macro("EjFwhmPulso", "%.1f" % (tau_r[arriba[-1]] - tau_r[arriba[0]])),
         macro("EjCarga", pct(ej["carga"])),
         macro("EjCargaComun", pct(ej["carga_comun"])),
         macro("EjPostMedido", "%.0f" % ej["ventana"]["post"]),
         macro("EjCargaPostMedido", pct(ej["ventana"]["fraccion"])),
         macro("EjCargaPico", pct(fr[i_pico])),
         macro("EjRegistro", "%.0f" % ej["espacio"]),
         macro("EjVariacionForma", pct(ej["variacion_forma"])),
         macro("EjAnchoBin", "%.1f" % (e_ej["rango_max"] / e_ej["bins"])),
         macro("EjRangoMax", "%.0f" % e_ej["rango_max"]),
         macro("EjUmbral", "%.0f" % umbral),
         macro("EjTPared", "%.0f" % e_ej["t_pared_muestra"]),
         macro("EjTVivo", "%.1f" % e_ej["t_vivo_muestra"]),
         macro("EjViva", "%.2f" % (100 * e_ej["frac_viva_muestra"])),
         macro("EjScaler", format(int(round(e_ej["scaler_muestra"] / 1000.0)),
                                  "d")),
         macro("EjTasaGuardada", "%.0f" % (e_ej["n_muestra"]
                                            / e_ej["t_pared_muestra"])),
         macro("EjKf", "%.2g" % e_ej["k_f"]),
         macro("EjPre", "%.1f" % ej["ventana"]["pre"]),
         macro("EjNFondo", format(int(e_ej["n_fondo"]), ",d")
               .replace(",", r"\,")),
         macro("EjTVivoFondo", "%.0f" % e_ej["t_vivo_fondo"]),
         # lo que valdria k_f normalizando por tiempo de PARED, y cuantas
         # veces mas grande es que el k_f con tiempo vivo
         macro("EjKfPared", "%.2f" % _kf_pared(e_ej)),
         macro("EjFactorVivo", "%.0f" % (_kf_pared(e_ej) / e_ej["k_f"]
                                          if e_ej["k_f"] > 0 else 0)),
         # la campana
         macro("RangoFactor", "%g" % espectro.RANGO_FACTOR),
         macro("RangoPercentil", "%g" % espectro.RANGO_PERCENTIL),
         macro("VivaCsMin", "%.2f" % (100 * min(viva_cs))),
         macro("VivaCsMax", "%.2f" % (100 * max(viva_cs))),
         macro("VivaCoMin", "%.0f" % (100 * min(viva_co))),
         macro("VivaCoMax", "%.0f" % (100 * max(viva_co))),
         macro("RegistroMin", "%.0f" % min(espacios)),
         macro("RegistroMax", "%.0f" % max(espacios)),
         macro("CargaComunMin", "%.0f" % (100 * min(cargas_comun))),
         macro("CargaComunMax", "%.0f" % (100 * max(cargas_comun))),
         macro("VariacionFormaMax", pct(max(variaciones))),
         macro("VariacionFormaMediana", pct(float(np.median(variaciones)))),
         macro("VariacionAltaMax", pct(max(var_altas))),
         macro("VariacionAltaMediana", pct(float(np.median(var_altas)))),
         macro("EjVariacionAlta", pct(ej["variacion_alta"])),
         macro("DescartadosMax", "%.1f" % (100 * max(marcados))),
         macro("DescartadosMediana", "%.1f" % (100 * float(np.median(marcados)))),
         macro("SaturadosMax", "%.2f" % (100 * max(sat))),
         macro("TruncadosMax", "%.1f" % (100 * max(trunc))),
         macro("SinPreMax", "%.1f" % (100 * max(sinpre))),
         ]
    return L


# --------------------------------------------------------------------------
# 7. Figuras que se copian tal cual
# --------------------------------------------------------------------------
def copiar_figuras(res):
    """
    Trae al informe las dos figuras-resumen de `calibrar_todo`, tal cual.

    No entran en el PDF: sus paneles ya estan como tablas y como figuras a
    tamano de informe (`ganancia_dv`, `cociente_fuentes`). Se copian porque el
    texto las cita por nombre como "la figura completa", y para que la
    carpeta `informe/` se pueda mandar entera sin que le falte nada.
    """
    os.makedirs(FIGURAS, exist_ok=True)
    copiadas = []
    for origen, nombre in (
            (os.path.join(ct.SALIDA, "RESUMEN_ganancia.png"),
             "resumen_ganancia.png"),
            (os.path.join(ct.SALIDA, "RESUMEN_linealidad.png"),
             "resumen_linealidad.png")):
        if os.path.exists(origen):
            shutil.copy2(origen, os.path.join(FIGURAS, nombre))
            copiadas.append(nombre)
    return copiadas


# --------------------------------------------------------------------------
# 8. Compilar
# --------------------------------------------------------------------------
def compilar():
    """
    pdflatex dos veces (la segunda resuelve el indice y las referencias).

    Se corre con `cwd` en la carpeta del informe y el nombre RELATIVO: la ruta
    absoluta del proyecto tiene espacios y acentos, y pasarsela a pdflatex es
    pedir problemas.
    """
    exe = shutil.which("pdflatex")
    if not exe:
        print("  no hay pdflatex en el PATH: queda el .tex sin compilar")
        return None
    bib = shutil.which("bibtex")
    pasadas = ("latex", "bibtex", "latex", "latex") if bib else ("latex", "latex")
    for i, paso in enumerate(pasadas):
        if paso == "bibtex":
            # la bibliografia: pdflatex deja las citas en el .aux, bibtex arma
            # el .bbl desde referencias.bib, y hacen falta dos pasadas mas
            # para que las citas y el indice queden resueltos
            subprocess.run([bib, "informe"], cwd=DESTINO, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
            continue
        r = subprocess.run([exe, "-interaction=nonstopmode", "informe.tex"],
                           cwd=DESTINO, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        if r.returncode != 0 and i == len(pasadas) - 1:
            print("  pdflatex termino con error; el log esta en "
                  "informe/informe.log")
            for linea in (r.stdout or "").splitlines():
                if linea.startswith("!"):
                    print("     %s" % linea)
            return None
    pdf = os.path.join(DESTINO, "informe.pdf")
    return pdf if os.path.exists(pdf) else None


def main():
    rutas.salida_utf8()
    os.makedirs(GENERADO, exist_ok=True)
    os.makedirs(FIGURAS, exist_ok=True)
    res = cargar()
    ta, tb = ct.test_a(res), ct.test_b(res)

    print("=" * 74)
    print("informe de calibracion  ·  %d corridas" % len(ct.ordenadas(res)))
    print("=" * 74)

    familia = aplicar_estilo()
    print("  figuras: %.2f in de ancho (= \\textwidth), letra de %d/%d pt, %s"
          % (ANCHO_TEXTO, ESTILO["axes.labelsize"], ESTILO["xtick.labelsize"],
             familia))
    print("  los pulsos de cada corrida (ventana medida + areas del cache)...")
    pul = resumen_pulsos(res)
    hechas = [figura_integracion(res), figura_pulsos_detectores(res),
              figura_area_a_espectro(res),
              figura_lineshape(res), figura_recocido(),
              figura_por_punto(res, ta), figura_detectores(res, tb),
              figura_ganancia(res, tb), figura_cociente(ta)]
    # el marco teorico medido: cierre con verdad conocida (cacheado), perfil
    # del chi2, el NCC por dentro y los metodos clasicos
    f_ej = corrida_ejemplo(res)
    print("  el cierre con verdad conocida sobre %s (%d realizaciones, "
          "Neyman y Pearson; se cachea)..." % (f_ej["clave"], N_CIERRE))
    cie = cierre_conocido(f_ej)
    perf = perfil_chi2(f_ej)
    nd = ncc_detalle(f_ej)
    sc = sesgo_clasicos()
    cr = clasicos_reales(res)
    hechas += [figura_estadistica(res, cie, perf), figura_ncc_detalle(nd),
               figura_clasicos(res, sc, f_ej)]
    # el chi2 y lo que dice: su distribucion en la campana, y los residuos
    # normalizados de la corrida de ejemplo (PLAN_estadistica, bloque 1)
    st = estadistica_chi2(res, f_ej)
    hechas.append(figura_pulls(st))
    # la fisica de base, con las secciones eficaces de NIST XCOM
    plastico, nai = xcom("C9H10"), xcom("NaI")
    hechas += [figura_interaccion(plastico, nai), figura_cinematica(plastico)]
    # la simulacion: que hay en los .root y cuanto importa el binning
    print("  la simulacion de Geant4 (los .root, y el binning reajustando %s)..."
          % f_ej["clave"])
    sim = simulacion_resumen()
    binning = binning_plantilla(f_ej)
    hechas.append(figura_simulacion(sim, f_ej))
    # un ajuste de ejemplo por detector, a 1600 V (el punto que tienen todos)
    for det in detectores(res):
        f = (por_clave(res, det, 1600, "Cs-137")
             or por_clave(res, det, 1600, "Co-60"))
        if f is not None:
            hechas.append(figura_ajuste(f, det))
    for ruta in hechas:
        print("    %s" % rutas.relativa(ruta))
        viejo = ruta[:-4] + ".png"        # la version anterior, PNG de 15 in
        if os.path.exists(viejo):
            os.remove(viejo)
    for nombre in copiar_figuras(res):
        print("    %s (copiada)" % nombre)

    # El test de la Fig. 5 del paper sobre nuestros datos: un detector por
    # corrida de Cs-137 a 1600 V, que es el punto que los cuatro tienen.
    p = paper()
    print("  el test del intervalo (Fig. 5 del paper) sobre nuestras corridas...")
    claves = [f["clave"] for f in
              (por_clave(res, d, 1600, "Cs-137") for d in detectores(res)) if f]
    rob = robustez_intervalo(res, claves) if p else {}
    for clave, d in sorted(rob.items()):
        print("    %-16s dispersion %.1f %%  (sin el mas angosto, %.1f %%)"
              % (clave, 100 * d["dispersion"],
                 100 * d["dispersion_sin_angosta"]))

    print("  tablas y numeros...")
    for nombre, texto in (
            ("numeros.tex", numeros(res, ta, tb, rob, p, pul)
             + "\n" + "\n".join(macros_teoria(res, cie, perf, nd, sc, cr))
             + "\n" + "\n".join(macros_chi2(st))
             + "\n" + "\n".join(macros_incertidumbre(res, cie, perf))
             + "\n" + "\n".join(macros_fisica(plastico, nai))
             + "\n" + "\n".join(macros_simulacion(geant4_geometria(), sim, binning))),
            ("tabla_simulacion.tex", tabla_simulacion(sim)),
            ("tabla_cinematica.tex", tabla_cinematica()),
            ("tabla_cierre.tex", tabla_cierre(cie)),
            ("tabla_clasicos.tex", tabla_clasicos(sc)),
            ("tabla_pulsos.tex", tabla_pulsos(res, pul)),
            ("tabla_calibraciones.tex", tabla_calibraciones(res)),
            ("tabla_testA.tex", tabla_test_a(res, ta)),
            ("tabla_testB.tex", tabla_test_b(res, tb)),
            ("tabla_detectores.tex", tabla_detectores(res, tb)),
            ("tabla_descartadas.tex", tabla_descartadas(res)),
            ("tabla_control.tex", tabla_control(res)),
            ("tabla_metricas_paper.tex", tabla_metricas_paper(res, ta)),
            ("tabla_robustez.tex", tabla_robustez(rob, p))):
        escribir(nombre, texto)
        print("    generado/%s" % nombre)

    if "--sin-pdf" not in sys.argv:
        print("  compilando...")
        pdf = compilar()
        if pdf:
            print()
            print("PDF -> %s" % rutas.relativa(pdf))
    print()
    print("El texto del informe se edita en informe/informe.tex; "
          "lo de generado/ se reescribe solo.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
