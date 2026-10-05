"""
doc_sa_codigo.py — genera `docs/recocido_simulado_codigo.pdf`.

Que es
------
El recorrido del CODIGO del recocido simulado, linea por linea, escrito para
alguien que nunca vio el algoritmo. El documento hermano
(`docs/recocido_simulado.pdf`, lo genera `doc_sa.py`) discute el metodo y la
evidencia: si el optimo es global, que constantes del paper no son
transferibles, cuanto vale el resultado. Este explica QUE HACE CADA LINEA.

Dos reglas de construccion, que son las que hacen que el PDF no pueda quedar
desincronizado del programa:

  * **el codigo no se tipea**: cada recuadro se extrae del fuente con
    `inspect.getsourcelines` al generar el PDF, con los numeros de linea reales
    del archivo. Si alguien edita `sa.py`, el proximo PDF muestra lo editado.
  * **los numeros tampoco**: todos salen de correr el ajuste sobre el Co-60 al
    momento de generarlo, incluida la traza de las primeras propuestas, que se
    reproduce con la misma semilla del generador aleatorio y se VERIFICA contra
    lo que devolvio `sa.recocido`.

    python analisis/doc_sa_codigo.py            # completo (~1 min)
    python analisis/doc_sa_codigo.py --rapido   # menos variantes, para iterar
"""

import ast
import inspect
import os
import sys
import tempfile
import textwrap
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analisis import rutas, modelo, sa, calibracion, gui, pipeline

TMP = tempfile.mkdtemp(prefix="doc_sa_cod_")
SALIDA = os.path.join(rutas.RAIZ, "docs", "recocido_simulado_codigo.pdf")

AZUL = "#1f4e79"
NARANJA = "#c4630a"
VERDE = "#1f7a3d"
GRIS = "#4a5560"
ROJO = "#e02020"

CLAVE = "co60"


# ==========================================================================
# 1. EXTRACCION DEL CODIGO
# ==========================================================================
# La idea de todo este bloque: el PDF no lleva una COPIA del codigo sino una
# VISTA del archivo. `inspect.getsourcelines` devuelve las lineas y el numero
# de la primera, asi que cada recuadro puede decir "sa.py, lineas 151-168" y
# esos numeros son ciertos aunque el archivo cambie manana.
def lineas_de(obj, sin_doc=True):
    """
    (numero_de_linea, texto) por cada linea del fuente de `obj`, sin sangria.

    `sin_doc` reemplaza el docstring por un marcador de una linea. No es por
    ahorrar espacio: los docstrings de este proyecto son el argumento de por
    que el codigo es asi, y ese argumento se cuenta en el texto del documento.
    Repetirlo adentro del recuadro haria leer dos veces lo mismo.
    """
    crudas, n0 = inspect.getsourcelines(obj)
    texto = textwrap.dedent("".join(crudas))
    filas = [(n0 + i, s) for i, s in enumerate(texto.splitlines())]
    if sin_doc:
        try:
            fn = ast.parse(texto).body[0]
            if ast.get_docstring(fn) is not None:
                d = fn.body[0]
                linea = filas[d.lineno - 1][1]
                sangria = " " * (len(linea) - len(linea.lstrip()))
                filas = (filas[:d.lineno - 1]
                         + [(n0 + d.lineno - 1, sangria + '"""..."""')]
                         + filas[d.end_lineno:])
        except (SyntaxError, IndexError):
            pass
    return filas


def recorte(filas, desde=None, hasta=None, quitar=()):
    """
    Sub-bloque de un fuente, ubicado por SUBCADENAS y no por numeros.

    Buscar por texto es lo que hace que el documento sobreviva a una edicion de
    `sa.py`: si alguien agrega tres lineas mas arriba, un recorte por numero
    mostraria otra cosa sin avisar; uno por texto sigue apuntando al mismo
    lugar, o falla ruidosamente al generar el PDF.
    """
    i0 = 0
    if desde is not None:
        i0 = next(i for i, (_, s) in enumerate(filas) if desde in s)
    i1 = len(filas)
    if hasta is not None:
        i1 = next(i for i, (_, s) in enumerate(filas) if hasta in s and i > i0)
    sub = filas[i0:i1]
    if quitar:
        sub = [f for f in sub if not any(q in f[1] for q in quitar)]
    # sin lineas en blanco al principio ni al final: el recuadro tendria un
    # renglon vacio y el pie diria un rango de lineas mas largo que lo que se ve
    while sub and not sub[0][1].strip():
        sub = sub[1:]
    while sub and _sobra_al_final(sub[-1][1]):
        sub = sub[:-1]
    return sub


def _sobra_al_final(linea):
    """Blanco, o una linea de guiones: no es codigo, es separacion."""
    t = linea.strip()
    return (not t) or set(t) <= set("#- ")


def lineas_modulo(mod):
    """Igual que `lineas_de`, pero para un archivo entero.

    No se puede usar `inspect.getsourcelines` con un modulo y numerar desde
    ahi: para un modulo devuelve 0 como primera linea, asi que todo el listado
    quedaria corrido en uno respecto de lo que muestra un editor.
    """
    with open(inspect.getsourcefile(mod), encoding="utf-8") as f:
        return [(i + 1, s) for i, s in enumerate(f.read().splitlines())]


def rango(filas):
    return filas[0][0], filas[-1][0]


def archivo_de(obj):
    return os.path.basename(inspect.getsourcefile(obj))


# ==========================================================================
# 2. CALCULOS  (todo numero que despues aparece en el texto)
# ==========================================================================
def calcular(rapido=False):
    R = {}
    print("[1/7] armando el problema (espectro + plantilla + NCC)...")
    t0 = time.time()
    p, semilla, ncc_r, esp, corr = modelo.armar(CLAVE, verbose=False)
    R.update(p=p, semilla=np.asarray(semilla, float), ncc=ncc_r, esp=esp,
             corrida=corr, t_armar=time.time() - t0)
    R["nbins"] = int(p.m.sum())
    R["ndf"] = p.ndf
    cotas = [(0.6 * semilla[0], 1.6 * semilla[0]), (0.02, 0.30)]
    R["cotas"] = np.asarray(cotas, float)
    R["delta"] = 0.1 * (R["cotas"][:, 1] - R["cotas"][:, 0])
    R["chi2_semilla"] = p.costo(semilla)

    # cuanto cuesta UNA evaluacion del costo: es el numero que decide si un
    # optimizador que pide 20 000 evaluaciones es viable o no
    t0 = time.time()
    for _ in range(200):
        p.costo(semilla)
    R["t_costo"] = (time.time() - t0) / 200

    print("[2/7] una corrida del recocido, con semilla fija...")
    t0 = time.time()
    res = sa.recocido(p.costo, semilla, cotas, semilla=1)
    R["t_sa"] = time.time() - t0
    xp, cp = sa.pulir(p.costo, res["x"])
    R.update(res=res, x=xp, chi2=cp, chi2_ndf=cp / p.ndf)
    R["A"] = p.chi2(*xp)[1]
    R["chi2_sin_pulir"] = res["costo"]

    print("[3/7] la traza de las primeras propuestas...")
    R["traza"], R["traza_ok"] = _traza(p, R["semilla"], R["cotas"], res)

    print("[4/7] errores por curvatura...")
    err, C = sa.errores_por_curvatura(p.costo, xp)
    R.update(err=err, C=C, corr_Dsigma=C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]))

    print("[5/7] multiarranque...")
    n = 4 if rapido else 12
    t0 = time.time()
    ms = sa.multiarranque(p.costo, cotas, n=n, semilla=7)
    R["t_ms"] = time.time() - t0
    R["ms"] = ms
    Ds = np.array([m["x_pulido"][0] for m in ms])
    R["ms_D"] = Ds
    R["ms_chi2"] = np.array([m["costo_pulido"] for m in ms])
    R["ms_dispersion"] = float(Ds.std())
    R["ms_iguales"] = int(np.sum(R["ms_chi2"] < R["ms_chi2"].min() + 0.1))

    print("[6/7] que pasa si se toca cada constante...")
    R["variantes"] = _variantes(p, R["semilla"], R["cotas"], rapido)
    R["curvas_T"] = _curvas_T(p, R["semilla"], R["cotas"], res["T_init"])
    R["malo"] = _arranque_malo(p, R["cotas"], xp)
    R["encerrado"] = _cotas_malas(p)
    R["local"] = _local_vs_sa(p, R["cotas"], n=(12 if rapido else 40))

    print("[7/7] la calibracion completa, como la corre la interfaz...")
    det = {}
    t0 = time.time()
    r_cal, p_cal, _, _ = calibracion.calibrar(
        corr, verbose=False, esp=esp, plantilla=(p.cen_E, p.pesos),
        detalles=det)
    R["t_cal"] = time.time() - t0
    R.update(r_cal=r_cal, det_cal=det, p_cal=p_cal)
    return R


def _traza(p, semilla, cotas, res):
    """
    Reproduce a mano las primeras propuestas del recocido de `res`.

    Es la MISMA cuenta que hace el bucle interno de `sa.recocido`, escrita
    aparte para poder mirar los numeros de a uno. Que sea la misma se verifica:
    al terminar las 100 propuestas del primer ciclo, el costo, la temperatura y
    la tasa de aceptacion tienen que coincidir con la primera fila de
    `res["historia"]`.

    Ojo con un detalle que parece cosmetico y no lo es: `dE <= 0 or rng.random()
    < ...` NO sortea cuando la propuesta mejora, porque Python corta la
    evaluacion del `or`. Replicar eso mal desincroniza el generador y la traza
    deja de ser la del recocido real.
    """
    lo, hi = cotas[:, 0], cotas[:, 1]
    delta = 0.1 * (hi - lo)
    rng = np.random.default_rng(1)
    x = np.asarray(semilla, float).copy()
    T = sa.temperatura_inicial(p.costo, x, delta, rng=rng)
    c = p.costo(x)
    filas, aceptados = [], 0
    for k in range(100):
        prop = np.clip(x + rng.uniform(-1, 1, size=x.size) * delta, lo, hi)
        cp = p.costo(prop)
        dE = cp - c
        if dE <= 0:
            prob, sorteo, acepta = 1.0, None, True
        else:
            prob = float(np.exp(-dE / max(T, 1e-300)))
            sorteo = float(rng.random())
            acepta = sorteo < prob
        if k < 12:
            filas.append(dict(k=k + 1, x=x.copy(), prop=prop.copy(), c=c,
                              cp=cp, dE=dE, T=T, prob=prob, sorteo=sorteo,
                              acepta=acepta))
        if acepta:
            x, c = prop, cp
            aceptados += 1
    h0 = res["historia"][0]
    ok = (abs(c - h0[1]) < 1e-9 and abs(aceptados / 100.0 - h0[3]) < 1e-12
          and abs(T - h0[0]) < 1e-9)
    return filas, bool(ok)


def _variantes(p, semilla, cotas, rapido):
    """Cada constante del algoritmo, movida a proposito, con n semillas."""
    n = 3 if rapido else 6
    casos = [
        ("T_init calibrada (el default)", {}),
        ("T_init = 2000, el numero del paper", dict(T_init=2000.0)),
        ("T_init 100 veces mas chica", dict(T_init="x0.01")),
        ("T_init 100 veces mas grande", dict(T_init="x100")),
        ("alfa = 0.85, enfria rapido", dict(alfa=0.85)),
        ("alfa = 0.99, enfria lento", dict(alfa=0.99)),
        ("paso 10 veces mas chico", dict(delta="x0.1")),
        ("paso = el rango entero de las cotas", dict(delta="x10")),
        ("10 propuestas por temperatura, no 100", dict(n_markov=10)),
    ]
    T_auto = sa.temperatura_inicial(p.costo, semilla,
                                    0.1 * (cotas[:, 1] - cotas[:, 0]),
                                    rng=np.random.default_rng(0))
    base_delta = 0.1 * (cotas[:, 1] - cotas[:, 0])
    salida = []
    for nombre, kw in casos:
        kw = dict(kw)
        if isinstance(kw.get("T_init"), str):
            kw["T_init"] = T_auto * float(kw["T_init"][1:])
        if isinstance(kw.get("delta"), str):
            kw["delta"] = base_delta * float(kw["delta"][1:])
        Ds, evs, cic, chis = [], [], [], []
        for s in range(1, n + 1):
            r = sa.recocido(p.costo, semilla, cotas, semilla=s, **kw)
            xp, cp = sa.pulir(p.costo, r["x"])
            Ds.append(xp[0]); evs.append(r["n_eval"]); cic.append(r["ciclos"])
            chis.append(cp)
        salida.append(dict(nombre=nombre, D=float(np.median(Ds)),
                           sd=float(np.std(Ds)), n_eval=float(np.mean(evs)),
                           ciclos=float(np.mean(cic)),
                           chi2=float(np.median(chis)), n=n))
    return salida


def _curvas_T(p, semilla, cotas, T_auto):
    """Aceptacion y costo ciclo a ciclo para tres temperaturas iniciales."""
    salida = []
    for etiqueta, T in (("T_init / 100", T_auto / 100),
                        ("T_init calibrada", T_auto),
                        ("T_init x 100", T_auto * 100)):
        r = sa.recocido(p.costo, semilla, cotas, T_init=T, semilla=11)
        salida.append((etiqueta, T, r))
    return salida


def _arranque_malo(p, cotas, x_bueno):
    """El experimento que justifica el algoritmo: local contra recocido."""
    from scipy.optimize import minimize
    x0 = (0.70 * x_bueno[0], 0.28)
    loc = minimize(p.costo, x0, method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-4, maxiter=4000))
    r = sa.recocido(p.costo, x0, cotas, semilla=3)
    xs, cs = sa.pulir(p.costo, r["x"])
    return dict(x0=np.asarray(x0, float), loc_x=loc.x, loc_chi2=float(loc.fun),
                sa_x=xs, sa_chi2=float(cs), n_eval_loc=int(loc.nfev),
                n_eval_sa=int(r["n_eval"]))


def _cotas_malas(p):
    """Cotas que dejan el optimo afuera: el recocido devuelve la pared."""
    cotas = [(1.35, 1.75), (0.02, 0.30)]
    r = sa.recocido(p.costo, (1.4, 0.1), cotas, semilla=5)
    xp, cp = sa.pulir(p.costo, r["x"])
    return dict(cotas=np.asarray(cotas, float), x=xp, chi2=float(cp),
                chi2_ndf=float(cp) / p.ndf)


def _local_vs_sa(p, cotas, n=40, semilla=3):
    """
    N arranques al azar, resueltos por el optimizador LOCAL y por el recocido.

    Es el experimento que dice si en ESTE paisaje hace falta el recocido, en
    vez de suponerlo. Ojo con la asimetria de la comparacion, que es real y
    conviene decirla: `minimize` de scipy corre SIN cotas, mientras que el
    recocido recorta cada propuesta con `np.clip`. O sea que no se comparan
    solo dos estrategias de busqueda: tambien se compara buscar adentro de lo
    fisicamente posible contra buscar en todo el plano.
    """
    from scipy.optimize import minimize
    rng = np.random.default_rng(semilla)
    filas = []
    for _ in range(n):
        x0 = rng.uniform(cotas[:, 0], cotas[:, 1])
        m = minimize(p.costo, x0, method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-4, maxiter=4000))
        r = sa.recocido(p.costo, x0, cotas, semilla=int(rng.integers(1 << 30)))
        xs, cs = sa.pulir(p.costo, r["x"])
        filas.append(dict(x0=x0, loc_x=m.x, loc_chi2=float(m.fun),
                          sa_x=xs, sa_chi2=float(cs)))
    mejor = min(min(f["loc_chi2"] for f in filas),
                min(f["sa_chi2"] for f in filas))
    return dict(filas=filas, n=n, mejor=mejor,
                loc_ok=sum(f["loc_chi2"] < mejor + 0.1 for f in filas),
                sa_ok=sum(f["sa_chi2"] < mejor + 0.1 for f in filas),
                loc_peor=max(filas, key=lambda f: f["loc_chi2"]))


# ==========================================================================
# 3. FIGURAS
# ==========================================================================
def _guardar(fig, nombre):
    ruta = os.path.join(TMP, nombre)
    fig.tight_layout()
    fig.savefig(ruta, dpi=170)
    plt.close(fig)
    return ruta


def ec(texto, tam=15):
    """Renderiza una formula con mathtext y la devuelve como PNG."""
    fig = plt.figure(figsize=(0.1, 0.1))
    t = fig.text(0, 0, texto, fontsize=tam, color="#101820")
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    fig.set_size_inches(bb.width / 170 + 0.12, bb.height / 170 + 0.10)
    ruta = os.path.join(TMP, "ec_%d.png" % (abs(hash(texto)) % 10 ** 9))
    fig.savefig(ruta, dpi=170, bbox_inches="tight", pad_inches=0.03,
                transparent=True)
    plt.close(fig)
    return ruta


def fig_cadena():
    """Donde vive el recocido dentro del analisis."""
    fig, ax = plt.subplots(figsize=(12.6, 2.4))
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 2)
    pasos = [("archivo\ndel DRS4", "drs.py"),
             ("area de\ncada pulso", "area.py"),
             ("espectro\nneto", "espectro.py"),
             ("ventana y\nsemilla", "ncc.py"),
             ("RECOCIDO\nD y sigma", "sa.recocido"),
             ("pulido\nlocal", "sa.pulir"),
             ("errores por\ncurvatura", "sa.py"),
             ("D y su error", "calibracion.py")]
    an, sep = 1.06, 0.19
    for i, (texto, mod) in enumerate(pasos):
        x = 0.06 + i * (an + sep)
        clave = texto.startswith("RECOCIDO")
        ax.add_patch(plt.Rectangle(
            (x, 0.72), an, 0.72,
            facecolor=("#fde7cf" if clave else "#eef2f6"),
            edgecolor=(NARANJA if clave else "#c3ccd5"),
            lw=(1.8 if clave else 1.0), zorder=2))
        ax.text(x + an / 2, 1.08, texto, ha="center", va="center",
                fontsize=8.2, color=("#7a3d05" if clave else "#20303c"),
                fontweight=("bold" if clave else "normal"), zorder=3)
        ax.text(x + an / 2, 0.55, mod, ha="center", va="center", fontsize=7,
                color=GRIS, style="italic", zorder=3)
        if i:
            ax.annotate("", xy=(x - 0.02, 1.08),
                        xytext=(x - sep + 0.02, 1.08),
                        arrowprops=dict(arrowstyle="->", color="#8894a0",
                                        lw=1.2))
    ax.text(5.0, 0.12, "los tres pasos de la derecha viven en el mismo modulo "
            "sa.py; la interfaz los muestra juntos en la pestana 3",
            ha="center", fontsize=7.4, color=GRIS)
    return _guardar(fig, "cadena.png")


def fig_juguete():
    """Por que hace falta: una funcion de juguete con varios valles."""
    f = lambda x: (x - 1.6) ** 2 / 6.0 - np.cos(3.4 * x) + 1.0
    xs = np.linspace(-4, 6, 900)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.5))

    from scipy.optimize import minimize
    a = axes[0]
    a.plot(xs, f(xs), lw=1.6, color=AZUL)
    for x0, col in ((-3.2, NARANJA), (4.6, "#8a2be2")):
        r = minimize(lambda v: f(v[0]), [x0], method="Nelder-Mead")
        a.plot([x0], [f(x0)], "o", color=col, ms=7)
        a.annotate("", xy=(r.x[0], f(r.x[0])), xytext=(x0, f(x0)),
                   arrowprops=dict(arrowstyle="->", color=col, lw=1.4))
        a.plot([r.x[0]], [f(r.x[0])], "*", color=col, ms=14)
    a.plot([1.6], [f(1.6)], "v", color=VERDE, ms=9)
    a.text(1.6, f(1.6) - 1.1, "minimo global", ha="center", fontsize=8,
           color=VERDE)
    a.set_title("bajar siempre termina en el valle de al lado", fontsize=10)
    a.set_xlabel("parametro")
    a.set_ylabel("costo")
    a.grid(alpha=.3)

    # la MISMA mecanica del recocido, en una dimension y a la vista
    rng = np.random.default_rng(4)
    cam, subidas = [], []
    x, T = -3.2, 4.0
    c = f(x)
    for k in range(240):
        prop = x + rng.uniform(-1, 1) * 1.4
        cp = f(prop)
        sube = cp > c
        if not sube or rng.random() < np.exp(-(cp - c) / T):
            if sube:
                subidas.append((prop, cp))
            x, c = prop, cp
        cam.append((k, x, c, T))
        T *= 0.975
    cam = np.array(cam)

    a = axes[1]
    a.plot(xs, f(xs), lw=1.2, color="0.78", zorder=1)
    sc = a.scatter(cam[:, 1], cam[:, 2], c=cam[:, 0], s=10, cmap="plasma",
                   zorder=3)
    if subidas:
        s = np.array(subidas)
        a.plot(s[:, 0], s[:, 1], "x", color=ROJO, ms=5, mew=1.2, zorder=4,
               label="%d saltos cuesta ARRIBA" % len(subidas))
        a.legend(fontsize=7.5, loc="upper center")
    a.plot([1.6], [f(1.6)], "v", color=VERDE, ms=9, zorder=5)
    cb = fig.colorbar(sc, ax=a)
    cb.set_label("iteracion", fontsize=8)
    a.set_title("el recocido acepta subir, y por eso escapa", fontsize=10)
    a.set_xlabel("parametro")
    a.set_ylabel("costo")
    a.grid(alpha=.3)

    a = axes[2]
    a.plot(cam[:, 0], cam[:, 1], lw=1.1, color=AZUL, label="donde esta")
    a.axhline(1.6, color=VERDE, ls="--", lw=1.2, label="minimo global")
    a2 = a.twinx()
    a2.semilogy(cam[:, 0], cam[:, 3], lw=1.2, color=NARANJA)
    a2.set_ylabel("temperatura", color=NARANJA, fontsize=9)
    a.set_title("al enfriarse deja de pasear y se queda", fontsize=10)
    a.set_xlabel("iteracion")
    a.set_ylabel("parametro")
    a.legend(fontsize=8, loc="lower right")
    a.grid(alpha=.3)
    return _guardar(fig, "juguete.png")


def fig_costo(R):
    """Que es exactamente el numero que el algoritmo minimiza."""
    p, x, A = R["p"], R["x"], R["A"]
    e = R["esp"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))

    a = axes[0]
    a.step(e["canales"], e["neto"], where="mid", lw=.9, color="0.55",
           label="espectro neto medido")
    a.axvspan(p.lo, p.hi, color=NARANJA, alpha=.12,
              label="ventana de ajuste\n[%.0f, %.0f]" % (p.lo, p.hi))
    A0 = p.chi2(*R["semilla"])[1]
    a.plot(p.canales, A0 * p.forma(*R["semilla"]), lw=1.1, ls="--",
           color=AZUL, label="modelo con la semilla del NCC")
    a.plot(p.canales, A * p.forma(*x), lw=1.7, color=ROJO,
           label="modelo con el optimo")
    a.set_xlim(0, 1.15 * p.hi)
    a.set_title("el dato, la ventana y el modelo", fontsize=10)
    a.set_xlabel("area [mV*ns]")
    a.set_ylabel("cuentas")
    a.legend(fontsize=7)
    a.grid(alpha=.3)

    a = axes[1]
    r0 = (p.dato - A0 * p.forma(*R["semilla"])) / p.error
    r1 = (p.dato - A * p.forma(*x)) / p.error
    a.axhline(0, color="k", lw=.8)
    a.step(p.canales, r0, where="mid", lw=.8, color=AZUL, alpha=.7,
           label="semilla: $\\chi^2$/ndf = %.2f" % (R["chi2_semilla"] / p.ndf))
    a.step(p.canales, r1, where="mid", lw=1.1, color=ROJO,
           label="optimo: $\\chi^2$/ndf = %.2f" % R["chi2_ndf"])
    a.set_title("los residuos, que es lo que el $\\chi^2$ suma", fontsize=10)
    a.set_xlabel("area [mV*ns]")
    a.set_ylabel("(dato - modelo) / error")
    a.legend(fontsize=8)
    a.grid(alpha=.3)

    a = axes[2]
    Ds = np.linspace(R["cotas"][0, 0], R["cotas"][0, 1], 160)
    a.plot(Ds, [p.costo((d, x[1])) / p.ndf for d in Ds], lw=1.5, color=AZUL)
    a.axvline(x[0], color=ROJO, ls="--", lw=1.2,
              label="optimo D = %.4f" % x[0])
    a.axvline(R["semilla"][0], color=NARANJA, ls=":", lw=1.4,
              label="semilla D = %.4f" % R["semilla"][0])
    a.set_yscale("log")
    a.set_title("el costo contra D, con $\\sigma$ fijo en el optimo",
                fontsize=10)
    a.set_xlabel("D [mV*ns/keVee]")
    a.set_ylabel("$\\chi^2$/ndf")
    a.legend(fontsize=8)
    a.grid(alpha=.3, which="both")
    return _guardar(fig, "costo.png")


def fig_metropolis(R):
    """Las dos formulas del algoritmo, con la T real de esta corrida."""
    T0 = R["res"]["T_init"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.4))

    a = axes[0]
    dE = np.linspace(0, 4 * T0, 400)
    for frac, col in ((1.0, AZUL), (0.1, NARANJA), (0.01, VERDE)):
        a.plot(dE, np.exp(-dE / (T0 * frac)), lw=1.6, color=col,
               label="T = %.3g" % (T0 * frac))
    a.set_title("P(aceptar algo peor) = exp($-\\Delta E$/T)", fontsize=10)
    a.set_xlabel("cuanto EMPEORA el costo, $\\Delta E$ [unidades de $\\chi^2$]")
    a.set_ylabel("probabilidad")
    a.legend(fontsize=8)
    a.grid(alpha=.3)

    a = axes[1]
    k = np.arange(200)
    for al, col in ((0.85, AZUL), (0.95, NARANJA), (0.99, VERDE)):
        a.semilogy(k, T0 * al ** k, lw=1.5, color=col,
                   label="$\\alpha$ = %.2f" % al)
    a.axhline(1e-6 * T0, color="k", ls=":", lw=1,
              label="corte en $10^{-6}\\,T_0$")
    a.set_title("enfriamiento geometrico $T_k = \\alpha^k T_0$", fontsize=10)
    a.set_xlabel("ciclo de temperatura")
    a.set_ylabel("T")
    a.legend(fontsize=8)
    a.grid(alpha=.3)

    a = axes[2]
    h = R["res"]["historia"]
    a.plot(h[:, 3], "o-", lw=1.4, ms=4, color=AZUL, label="fraccion aceptada")
    a.set_ylim(-0.03, 1.03)
    a.set_xlabel("ciclo de temperatura")
    a.set_ylabel("fraccion aceptada", color=AZUL)
    a2 = a.twinx()
    a2.semilogy(h[:, 0], lw=1.3, color=NARANJA)
    a2.set_ylabel("T", color=NARANJA)
    a.set_title("y lo que paso de verdad en esta corrida", fontsize=10)
    a.legend(fontsize=8, loc="center right")
    a.grid(alpha=.3)
    return _guardar(fig, "metropolis.png")


def fig_corrida(R):
    """La corrida completa: como baja el costo y por donde paso."""
    p, h = R["p"], R["res"]["historia"]
    hx = R["res"]["historia_x"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))

    a = axes[0]
    a.semilogy(h[:, 1] / p.ndf, lw=1.0, color=AZUL, alpha=.8,
               label="costo donde esta la cadena")
    a.semilogy(h[:, 2] / p.ndf, lw=1.9, color=ROJO,
               label="mejor costo visto hasta aqui")
    a.axhline(R["chi2_semilla"] / p.ndf, color=NARANJA, ls=":", lw=1.3,
              label="semilla del NCC (%.2f)" % (R["chi2_semilla"] / p.ndf))
    a.axhline(R["chi2_ndf"], color="k", ls="--", lw=1,
              label="final, ya pulido (%.3f)" % R["chi2_ndf"])
    a.set_title("el costo, ciclo a ciclo", fontsize=10)
    a.set_xlabel("ciclo de temperatura")
    a.set_ylabel("$\\chi^2$/ndf")
    a.legend(fontsize=7.5)
    a.grid(alpha=.3, which="both")

    a = axes[1]
    Ds = np.linspace(R["cotas"][0, 0], R["cotas"][0, 1], 130)
    Ss = np.linspace(R["cotas"][1, 0], R["cotas"][1, 1], 110)
    Z = p.mapa(Ds, Ss) / p.ndf
    im = a.pcolormesh(Ds, Ss, np.log10(Z), shading="auto", cmap="viridis")
    a.contour(Ds, Ss, Z * p.ndf,
              levels=[R["chi2"] + 2.30, R["chi2"] + 6.17, R["chi2"] + 11.8],
              colors="w", linewidths=.9)
    a.plot(hx[:, 0], hx[:, 1], "-", lw=.9, color="0.85", zorder=3)
    a.scatter(hx[:, 0], hx[:, 1], c=np.arange(len(hx)), s=22, cmap="autumn",
              zorder=4, edgecolor="k", linewidth=.3)
    a.plot([R["semilla"][0]], [R["semilla"][1]], "*", ms=13, color="w",
           mec="k", zorder=5, label="semilla del NCC")
    a.plot([R["x"][0]], [R["x"][1]], "o", ms=8, color=ROJO, mec="k",
           zorder=6, label="optimo")
    cb = fig.colorbar(im, ax=a)
    cb.set_label("$\\log_{10}(\\chi^2$/ndf)", fontsize=8)
    a.set_title("el paisaje y el camino de la cadena", fontsize=10)
    a.set_xlabel("D [mV*ns/keVee]")
    a.set_ylabel("$\\sigma_{rel}$")
    a.legend(fontsize=7.5, loc="upper left")

    a = axes[2]
    for etiqueta, T, r in R["curvas_T"]:
        a.plot(r["historia"][:, 3], lw=1.4,
               label="%s (T = %.3g)" % (etiqueta, T))
    a.set_ylim(-0.03, 1.03)
    a.set_title("si T_init esta mal, se nota en la aceptacion", fontsize=10)
    a.set_xlabel("ciclo de temperatura")
    a.set_ylabel("fraccion aceptada")
    a.legend(fontsize=7.5)
    a.grid(alpha=.3)
    return _guardar(fig, "corrida.png")


def fig_multiarranque(R):
    """La evidencia de que el optimo es el mismo desde cualquier lado."""
    p, ms = R["p"], R["ms"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.9))

    a = axes[0]
    Ds = np.linspace(R["cotas"][0, 0], R["cotas"][0, 1], 130)
    Ss = np.linspace(R["cotas"][1, 0], R["cotas"][1, 1], 110)
    Z = p.mapa(Ds, Ss) / p.ndf
    a.pcolormesh(Ds, Ss, np.log10(Z), shading="auto", cmap="viridis",
                 alpha=.85)
    for m in ms:
        hx = m["historia_x"]
        if hx.size:
            a.plot(hx[:, 0], hx[:, 1], lw=.7, color="0.85", alpha=.8, zorder=2)
        a.plot([m["x0"][0]], [m["x0"][1]], "s", ms=5, color="w", mec="k",
               mew=.5, zorder=3)
        a.plot([m["x_pulido"][0]], [m["x_pulido"][1]], "o", ms=5, color=ROJO,
               zorder=4)
    a.set_xlim(R["cotas"][0, 0], R["cotas"][0, 1])
    a.set_ylim(R["cotas"][1, 0], R["cotas"][1, 1])
    a.set_title("%d arranques al azar (cuadrados) y sus finales (circulos)"
                % len(ms), fontsize=10)
    a.set_xlabel("D [mV*ns/keVee]")
    a.set_ylabel("$\\sigma_{rel}$")

    a = axes[1]
    a.scatter([m["x0"][0] for m in ms], R["ms_D"], s=48, color=AZUL, zorder=3,
              label="cada arranque")
    a.axhline(R["x"][0], color=ROJO, ls="--", lw=1.2,
              label="D = %.4f" % R["x"][0])
    rango_y = max(3e-4, 6 * R["ms_dispersion"])
    a.set_ylim(R["x"][0] - rango_y, R["x"][0] + rango_y)
    a.set_title("D final contra D inicial\n(dispersion %.1e)"
                % R["ms_dispersion"], fontsize=10)
    a.set_xlabel("D del punto INICIAL [mV*ns/keVee]")
    a.set_ylabel("D FINAL")
    a.legend(fontsize=8)
    a.grid(alpha=.3)
    return _guardar(fig, "multiarranque.png")


# ==========================================================================
# 4. EL PDF
# ==========================================================================
# Que hace cada argumento de `sa.recocido`. La tabla del documento se arma
# cruzando esto con `inspect.signature`, asi que los VALORES POR DEFECTO nunca
# se tipean: si alguien cambia uno en sa.py, el PDF lo muestra cambiado.
ARGUMENTOS = {
    "costo": "la funcion a minimizar. Recibe el vector de parametros y "
             "devuelve un numero. Al algoritmo no le importa que hay adentro.",
    "x0": "punto de partida. Aca, la semilla que entrego el NCC.",
    "cotas": "(minimo, maximo) de cada parametro. Acotan la busqueda a lo "
             "fisicamente posible y, de paso, fijan la escala del paso.",
    "delta": "cuanto se mueve como maximo cada parametro en una propuesta. "
             "Por defecto, el 10% del ancho de sus cotas: es el 0.1*d_i "
             "del paper.",
    "T_init": "temperatura inicial. 'auto' la calibra midiendo el costo, que "
              "es lo unico que la hace independiente de las unidades.",
    "alfa": "cuanto se enfria por ciclo: T_(k+1) = alfa * T_k. El paper usa "
            "0.95 y recomienda el rango 0.85 a 0.99.",
    "n_markov": "propuestas por temperatura. Las 100 del paper.",
    "T_min_rel": "corte por temperatura, como fraccion de la inicial.",
    "n_estable": "cuantos ciclos seguidos sin mejorar hacen falta para cortar.",
    "tol_rel": "que tan poco es 'no mejorar', en relativo. Relativo y no "
               "absoluto porque un corte absoluto depende de las unidades del "
               "costo, igual que T.",
    "paso_adaptativo": "agranda o achica el paso para mantener la aceptacion "
                       "cerca de un objetivo. Apagado: sobre este problema no "
                       "cambia nada y el default del paper es paso fijo.",
    "objetivo_aceptacion": "la aceptacion que persigue el paso adaptativo.",
    "semilla": "semilla del generador aleatorio. Fijarla hace que la corrida "
               "sea reproducible bit a bit.",
    "max_ciclos": "tope duro de ciclos de temperatura, por si ningun corte "
                  "se cumple.",
}

PANELES_GUI = [
    ("1 - la tabla",
     "de donde arranca cada corrida y donde termina",
     "muestra la semilla del NCC, las cotas y los N arranques con su D final. "
     "Los arranques NO salen de la semilla: salen de puntos al azar adentro "
     "de las cotas."),
    ("2 - el costo por ciclo",
     "como baja, en las N corridas a la vez",
     "la linea gruesa es la corrida ganadora y las finas las otras. Que "
     "terminen apiladas es la evidencia practica de que el minimo es global."),
    ("3 - T y aceptacion",
     "la mecanica del enfriamiento",
     "al principio se acepta casi todo (la cadena pasea), al final casi nada "
     "(se comporta como un optimizador local). Si arranca en ~1 o en ~0, "
     "T_init esta mal calibrada."),
    ("4 - el camino",
     "por donde paseo la busqueda en el plano (D, sigma)",
     "es lo unico que distingue 'bajo en linea recta' de 'cruzo de un valle a "
     "otro', que es la unica razon por la que se usa recocido."),
]


def nl(filas, sub):
    """Numero de linea real de la primera fila que contiene `sub`."""
    return next(n for n, s in filas if sub in s)


def construir(R):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                    Image, Table, TableStyle, PageBreak,
                                    KeepTogether, XPreformatted)
    from analisis.doc_sa import _registrar_fuentes

    F, FB, FI, FM = _registrar_fuentes()
    ANCHO = A4[0] - 4 * cm
    az = colors.HexColor(AZUL)
    gr = colors.HexColor("#555f6a")

    est = dict(
        titulo=ParagraphStyle("t", fontName=FB, fontSize=21, leading=25,
                              textColor=az, spaceAfter=6),
        sub=ParagraphStyle("s", fontName=F, fontSize=11.5, leading=15,
                           textColor=gr, spaceAfter=16),
        h1=ParagraphStyle("h1", fontName=FB, fontSize=14, leading=17,
                          textColor=az, spaceBefore=16, spaceAfter=7),
        h2=ParagraphStyle("h2", fontName=FB, fontSize=11.5, leading=14,
                          textColor=colors.HexColor("#101820"),
                          spaceBefore=11, spaceAfter=4),
        p=ParagraphStyle("p", fontName=F, fontSize=9.7, leading=13.6,
                         spaceAfter=7, alignment=4),
        li=ParagraphStyle("li", fontName=F, fontSize=9.7, leading=13.4,
                          leftIndent=14, bulletIndent=4, spaceAfter=4,
                          alignment=4),
        cap=ParagraphStyle("cap", fontName=FI, fontSize=8.3, leading=11,
                           textColor=gr, spaceBefore=3, spaceAfter=12,
                           alignment=4),
        nota=ParagraphStyle("n", fontName=F, fontSize=9.3, leading=12.6,
                            leftIndent=10, textColor=gr, spaceAfter=8,
                            alignment=4),
        ref=ParagraphStyle("r", fontName=FI, fontSize=7.8, leading=10,
                           textColor=gr, spaceBefore=1, spaceAfter=10),
        cod=ParagraphStyle("c", fontName=FM, fontSize=7.5, leading=9.6,
                           textColor=colors.HexColor("#16232e"),
                           backColor=colors.HexColor("#f3f6f9"),
                           borderPadding=(5, 6, 5, 6), spaceBefore=4,
                           spaceAfter=2),
    )

    def P(t, e="p"):
        return Paragraph(t, est[e])

    def B(t):
        return Paragraph(t, est["li"], bulletText="•")

    def N(t):
        return Paragraph(t, est["li"], bulletText="–")

    def IMG(ruta, ancho=ANCHO):
        from PIL import Image as PILImage
        w, h = PILImage.open(ruta).size
        return Image(ruta, width=ancho, height=ancho * h / w)

    def FIG(ruta, epigrafe):
        return KeepTogether([IMG(ruta), Paragraph(epigrafe, est["cap"])])

    def ECU(texto, escala=1.0):
        r = ec(texto)
        from PIL import Image as PILImage
        w, h = PILImage.open(r).size
        an = min(w / 170 * 72 * escala, ANCHO)
        return Image(r, width=an, height=an * h / w, hAlign="CENTER")

    def CODIGO(filas, obj, titulo=None):
        """Un recuadro de codigo REAL, con los numeros de linea del archivo."""
        def esc(s):
            return (s.replace("&", "&amp;").replace("<", "&lt;")
                     .replace(">", "&gt;"))
        cuerpo = "\n".join(
            '<font color="#9fb0bf">%4d</font>  %s' % (n, esc(s))
            for n, s in filas)
        a, b = rango(filas)
        pie_txt = "%s &middot; lineas %d a %d%s" % (
            archivo_de(obj), a, b,
            ("" if titulo is None else " &middot; " + titulo))
        return [XPreformatted(cuerpo, est["cod"]), Paragraph(pie_txt,
                                                             est["ref"])]

    def TABLA(datos, anchos=None, alinear_num=True, cuerpo=8.6,
              izquierda=()):
        t = Table(datos, colWidths=anchos, hAlign="LEFT")
        estilo = [
            ("FONT", (0, 0), (-1, 0), FB, cuerpo),
            ("FONT", (0, 1), (-1, -1), F, cuerpo),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("BACKGROUND", (0, 0), (-1, 0), az),
            ("LINEBELOW", (0, 0), (-1, -1), 0.4, colors.HexColor("#d5dce3")),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 3.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ]
        if alinear_num:
            estilo.append(("ALIGN", (1, 1), (-1, -1), "RIGHT"))
        for k in izquierda:
            estilo.append(("ALIGN", (k, 1), (k, -1), "LEFT"))
        t.setStyle(TableStyle(estilo))
        return t

    def CELDA(t, negrita=False, tam=8.4):
        return Paragraph(t, ParagraphStyle(
            "celda", fontName=(FB if negrita else F), fontSize=tam,
            leading=tam * 1.25))

    def pie(canv, doc):
        canv.saveState()
        canv.setFont(F, 7.5)
        canv.setFillColor(gr)
        canv.drawString(2 * cm, 1.3 * cm,
                        "Calibracion por borde Compton  -  el recocido "
                        "simulado, linea por linea")
        canv.drawRightString(A4[0] - 2 * cm, 1.3 * cm, "%d" % doc.page)
        canv.setStrokeColor(colors.HexColor("#dde3ea"))
        canv.line(2 * cm, 1.6 * cm, A4[0] - 2 * cm, 1.6 * cm)
        canv.restoreState()

    # ---- los fuentes que se van a mostrar --------------------------------
    L_costo = lineas_de(modelo.Problema.costo)
    L_chi2 = lineas_de(modelo.Problema.chi2)
    L_amp = lineas_de(modelo.amplitud_optima)
    L_temp = lineas_de(sa.temperatura_inicial)
    L_rec = lineas_de(sa.recocido)
    L_pul = lineas_de(sa.pulir)
    L_multi = lineas_de(sa.multiarranque)
    L_err = lineas_de(sa.errores_por_curvatura)
    L_cal = lineas_de(calibracion.calibrar)
    L_pipe = lineas_de(pipeline.calibrar_completo)

    rec_prep = recorte(L_rec, desde="def recocido", hasta="for ciclo in")
    rec_bucle = recorte(L_rec, desde="for ciclo in", hasta="# paso adaptativo")
    rec_cortes = recorte(L_rec, desde="# corte por estabilizacion",
                         hasta="return dict")
    rec_ret = recorte(L_rec, desde="return dict")
    cal_frag = recorte(L_cal, desde="cotas = [", hasta="if detalles is not")
    pipe_frag = recorte(L_pipe, desde="detalles = {}",
                        hasta='pr.etapa("calib")')

    x, err, ms = R["x"], R["err"], R["ms"]
    p, res, rc = R["p"], R["res"], R["r_cal"]
    E = []

    # ------------------------------------------------------------ portada -
    E += [
        P("El recocido simulado, linea por linea", "titulo"),
        P("Que hace cada instruccion del optimizador que la interfaz corre en "
          "la etapa &laquo;3 &middot; Semilla y recocido&raquo;, explicado "
          "desde cero: primero para que sirve, despues el algoritmo con una "
          "funcion de juguete, y despues el codigo, bloque por bloque, con "
          "una corrida real al lado.", "sub"),
    ]
    E += [P("Como esta hecho este documento", "h2")]
    E += [B("<b>El codigo no esta copiado</b>: cada recuadro se extrae del "
            "archivo fuente al generar el PDF, con los numeros de linea "
            "reales. Si manana alguien edita <font face='%s'>sa.py</font>, el "
            "proximo PDF muestra lo editado, no esto." % FM)]
    E += [B("<b>Los numeros tampoco</b>: salen de correr el ajuste sobre la "
            "corrida <b>%s</b> (%s, %s) mientras se arma el documento. La "
            "traza de las primeras propuestas se reproduce con la misma "
            "semilla del generador aleatorio y se verifica contra lo que "
            "devolvio el recocido de verdad: <b>%s</b>."
            % (R["corrida"].clave, R["corrida"].fuente,
               R["corrida"].punto(),
               "coincide" if R["traza_ok"] else "NO COINCIDE"))]
    E += [B("El documento hermano, <font face='%s'>docs/recocido_simulado.pdf"
            "</font>, discute otra cosa: el metodo, que dice el paper y que "
            "evidencia hay de que el resultado es el optimo global. Este "
            "explica el programa." % FM)]

    E += [P("Un aviso sobre el modelo que se usa para explicar", "h2")]
    E += [P(
        "Todo el documento trabaja con el modelo <b>desnudo de dos "
        "parametros</b> (D y &sigma;<sub>rel</sub>), porque un paisaje de dos "
        "dimensiones se puede dibujar y mirar. La calibracion de verdad agrega "
        "una tercera componente -el exceso de cuentas a bajo canal- y da otro "
        "numero: <b>D = %.4f</b> en vez de <b>%.4f</b>. El algoritmo es el "
        "mismo, con un parametro mas; la seccion 8 lo muestra tal como lo "
        "corre la interfaz." % (rc["D"], x[0]))]

    # ------------------------------------------------------------- 1 -----
    E += [P("1. Donde vive el recocido", "h1")]
    E += [P(
        "El recocido no lee archivos ni sabe de fisica. Es el anteultimo "
        "eslabon de una cadena que ya hizo todo el trabajo sucio: leer el "
        "archivo del digitalizador, integrar el area de cada pulso, armar el "
        "histograma con resta de fondo, y elegir automaticamente en que tramo "
        "del espectro conviene ajustar. Lo que recibe es un problema "
        "matematico limpio: <b>una funcion que devuelve un numero, y hay que "
        "encontrar donde ese numero es minimo</b>.")]
    E += [FIG(R["fig_cadena"],
              "La cadena completa. El recuadro naranja es lo que explica este "
              "documento; los dos que le siguen viven en el mismo archivo y "
              "la interfaz los muestra en la misma pestana.")]
    E += [P(
        "En la interfaz (<font face='%s'>analisis/gui.py</font>) esto es la "
        "pestana <b>3 &middot; Semilla y recocido</b>. La etapa anterior, "
        "<b>2 &middot; NCC e intervalo</b>, es la que entrega las dos cosas "
        "que el recocido necesita para arrancar: la <b>ventana</b> de ajuste "
        "y la <b>semilla</b>. Para esta corrida: ventana <b>[%.0f, %.0f] "
        "mV&middot;ns</b> con &rho; = %.5f, y semilla <b>D = %.4f</b>, "
        "&sigma;<sub>rel</sub> = %.3f."
        % (FM, R["ncc"]["lo"], R["ncc"]["hi"], R["ncc"]["rho"],
           R["semilla"][0], R["semilla"][1]))]

    E += [P("1.1 Los archivos que intervienen", "h2")]
    E += [TABLA(
        [["archivo", "que aporta"],
         [CELDA("<font face='%s'>analisis/modelo.py</font>" % FM),
          CELDA("arma el <b>problema</b>: junta datos, plantilla y ventana, y "
                "expone la funcion de costo")],
         [CELDA("<font face='%s'>analisis/sa.py</font>" % FM),
          CELDA("el <b>algoritmo</b>: temperatura inicial, recocido, pulido, "
                "multiarranque y barras de error")],
         [CELDA("<font face='%s'>analisis/calibracion.py</font>" % FM),
          CELDA("la <b>receta</b>: define las cotas, corre N recocidos, se "
                "queda con el mejor y traduce a calibracion")],
         [CELDA("<font face='%s'>analisis/pipeline.py</font>" % FM),
          CELDA("el <b>orquestador</b> que llaman las dos interfaces, para "
                "que no haya dos copias de la receta")],
         [CELDA("<font face='%s'>analisis/gui.py</font>" % FM),
          CELDA("los <b>cuatro paneles</b> de la pestana 3, que dibujan lo "
                "que el recocido dejo grabado")]],
        anchos=(5.2 * cm, ANCHO - 5.2 * cm), alinear_num=False)]

    # ------------------------------------------------------------- 2 -----
    E += [P("2. El vocabulario minimo", "h1")]
    E += [P(
        "Seis palabras alcanzan para leer todo lo que sigue. Vale la pena "
        "fijarlas antes, porque el algoritmo es corto y lo dificil no es el "
        "codigo sino el vocabulario.")]
    E += [TABLA(
        [["palabra", "que significa aca"],
         [CELDA("<b>parametro</b>"),
          CELDA("un numero que no conocemos y queremos determinar. Aca son "
                "dos: <b>D</b>, la ganancia en mV&middot;ns por keV, que es "
                "el resultado buscado; y <b>&sigma;<sub>rel</sub></b>, cuanto "
                "borronea el detector, como fraccion del canal del borde.")],
         [CELDA("<b>vector x</b>"),
          CELDA("los parametros juntos, en orden fijo. Para el algoritmo el "
                "punto x = (%.4f, %.4f) es una posicion, nada mas."
                % (x[0], x[1]))],
         [CELDA("<b>costo</b>"),
          CELDA("un numero que mide cuan mal le va a ese x contra los datos. "
                "Aca es el &chi;<sup>2</sup>. Cuanto mas chico, mejor "
                "describe el modelo lo medido.")],
         [CELDA("<b>paisaje</b>"),
          CELDA("el valor del costo en cada punto del plano (D, "
                "&sigma;<sub>rel</sub>). Es una superficie con valles y lomas: "
                "buscar el minimo es buscar el fondo del valle mas hondo.")],
         [CELDA("<b>minimo local</b>"),
          CELDA("el fondo de UN valle. Un optimizador que solo baja termina "
                "ahi y no tiene forma de saber, desde adentro, si hay un "
                "valle mas hondo del otro lado de la loma.")],
         [CELDA("<b>minimo global</b>"),
          CELDA("el fondo del valle mas hondo de todos. Es el que se busca.")]],
        anchos=(3.4 * cm, ANCHO - 3.4 * cm), alinear_num=False)]

    # ------------------------------------------------------------- 3 -----
    E += [P("3. Que numero se esta minimizando", "h1")]
    E += [P(
        "Antes del algoritmo hay que ver que es la funcion de costo, porque "
        "es lo unico que el recocido conoce del problema. El modelo dice que "
        "el espectro medido tiene que verse asi:")]
    E += [ECU(r"$S(ch)\;=\;A\;\cdot\;\left[\,P\!\left(\frac{ch}{D}\right)"
              r"\;\otimes\;G(\sigma)\,\right]$")]
    E += [P(
        "En palabras: se toma la <b>plantilla</b> P, que es el espectro que "
        "predice la simulacion de Geant4 en energia; se la lleva al eje "
        "medido dividiendo por <b>D</b> (que es justamente la ganancia que "
        "buscamos, en mV&middot;ns por keV); se la borronea con una gaussiana "
        "de ancho <b>&sigma;</b>, porque el detector no mide con precision "
        "infinita; y se la escala con una amplitud <b>A</b> para que tenga "
        "tantas cuentas como el histograma medido.")]
    E += [P(
        "Comparar ese modelo con el dato bin por bin, pesando cada bin por su "
        "barra de error, da el &chi;<sup>2</sup>:")]
    E += [ECU(r"$\chi^2(\theta)\;=\;\sum_i\left(\frac{dato_i - "
              r"modelo_i(\theta)}{error_i}\right)^{\!2}$")]
    E += [P(
        "Sobre esta corrida la ventana tiene <b>%d bins</b> y el modelo "
        "%d parametros, o sea <b>%d grados de libertad</b>. Un ajuste que "
        "describe bien los datos da &chi;<sup>2</sup> parecido al numero de "
        "grados de libertad; por eso se reporta el cociente "
        "&chi;<sup>2</sup>/ndf, que deberia dar alrededor de 1."
        % (R["nbins"], R["nbins"] - R["ndf"], R["ndf"]))]
    E += [FIG(R["fig_costo"],
              "Izquierda: el espectro medido, la ventana que eligio el NCC y "
              "el modelo antes y despues de ajustar. Centro: los residuos, "
              "que es lo que el &chi;<sup>2</sup> eleva al cuadrado y suma. "
              "Derecha: el costo como funcion de D, con &sigma; congelado. "
              "Ese corte es una rebanada del paisaje que el recocido recorre.")]

    E += [P("3.1 El costo, en el codigo", "h2")]
    E += [P(
        "El optimizador llama siempre a la misma funcion, y es de tres "
        "lineas utiles. Lo que hace primero es rechazar lo imposible: un "
        "ancho negativo no significa nada, asi que en vez de calcular devuelve "
        "un numero enorme. Ese 1e30 es un <b>muro</b>: el algoritmo lo ve como "
        "un costo altisimo y no vuelve.")]
    E += CODIGO(L_costo, modelo.Problema)
    E += [P(
        "Y el &chi;<sup>2</sup> propiamente dicho. La linea "
        "<font face='%s'>r = (self.dato - A * f) / self.error</font> es la "
        "formula de arriba escrita en numpy: resta vectorial, division "
        "vectorial, y despues la suma de los cuadrados." % FM)]
    E += CODIGO(recorte(L_chi2, desde="self.n_eval += 1",
                        hasta="g = self.exceso_col"), modelo.Problema)

    E += [P("3.2 Un parametro que no hace falta buscar", "h2")]
    E += [P(
        "El modelo tiene tres parametros -A, D y &sigma;- pero el recocido "
        "busca en dos. La razon es que <b>A entra linealmente</b>: fijados D y "
        "&sigma;, la amplitud que minimiza el &chi;<sup>2</sup> se despeja "
        "exacto, sin buscarla. Se llama proyeccion, o VARPRO, y en la practica "
        "compra dos cosas:")]
    E += [B("una dimension menos que explorar, que es exponencialmente mas "
            "barato;")]
    E += [B("no hay que elegir un paso comun para parametros de escalas "
            "distintas. Medido en esta corrida: A vale %.4g y D vale %.4f, o "
            "sea que un mismo paso no le sirve a ninguno de los dos."
            % (R["A"], x[0]))]
    E += CODIGO(L_amp, modelo)
    E += [P(
        "Con eso, el problema que ve el algoritmo es el mas simple posible: "
        "<b>una funcion de dos numeros que devuelve un numero</b>. Cuesta "
        "%.0f microsegundos evaluarla, que es el dato que decide si un metodo "
        "que pide decenas de miles de evaluaciones es viable."
        % (R["t_costo"] * 1e6))]

    # ------------------------------------------------------------- 4 -----
    E += [P("4. La idea, antes del codigo", "h1")]
    E += [P(
        "El metodo obvio para minimizar es bajar: mirar alrededor, ir hacia "
        "donde el costo disminuye, repetir. Eso es un <b>optimizador local</b> "
        "y funciona perfecto mientras haya un solo valle. Cuando hay varios, "
        "termina en el que le toco por donde arranco, y no hay manera de "
        "darse cuenta desde adentro.")]
    E += [P(
        "El recocido simulado cambia una sola regla: <b>a veces acepta "
        "moverse a un punto peor</b>. Esa es la unica forma de subir una loma "
        "y ver que hay del otro lado. Cuanto acepta lo controla un numero "
        "llamado <b>temperatura</b>: alta al principio, cuando conviene "
        "pasear; baja despues, cuando conviene afinar. El nombre viene de la "
        "metalurgia, donde enfriar despacio deja que los atomos encuentren la "
        "red de minima energia y enfriar de golpe los congela desordenados.")]
    E += [FIG(R["fig_juguete"],
              "El algoritmo entero sobre una funcion de juguete de una "
              "variable. Izquierda: dos arranques, dos valles distintos, y "
              "ninguno es el bueno. Centro: el recorrido del recocido desde "
              "el mismo arranque malo; las cruces rojas son los movimientos "
              "que EMPEORARON el costo y aun asi fueron aceptados. Derecha: "
              "al bajar la temperatura la cadena deja de saltar y se queda.")]

    E += [P("4.1 La regla de aceptacion", "h2")]
    E += [P("Cada paso propone un punto nuevo y decide si se muda o se queda:")]
    E += [ECU(r"$P(aceptar) = 1$  si el costo baja;      "
              r"$P(aceptar) = e^{-\Delta E/T}$  si empeora en $\Delta E$",
              escala=0.68)]
    E += [P(
        "Vale la pena leer la exponencial despacio, porque toda la mecanica "
        "del algoritmo esta ahi:")]
    E += [B("si el empeoramiento &Delta;E es chico comparado con T, el "
            "exponente es casi cero y la probabilidad casi 1: se acepta casi "
            "siempre;")]
    E += [B("si &Delta;E es grande comparado con T, la probabilidad se "
            "desploma: un salto muy malo casi nunca pasa;")]
    E += [B("cuando T tiende a cero, ninguna subida se acepta y el metodo se "
            "vuelve un optimizador local corriente.")]
    E += [P(
        "Lo importante de la comparacion es que <b>&Delta;E y T se miden en "
        "las mismas unidades</b>, las del costo. Un mismo valor de T puede ser "
        "enorme o ridiculo segun en que este medido el costo, y esa es la "
        "razon por la que el codigo calibra T en vez de copiarle el numero al "
        "paper (seccion 5.1).")]
    E += [P(
        "El enfriamiento es geometrico: <b>T se multiplica por &alpha; "
        "= %.2f</b> despues de cada tanda de propuestas. Con eso T cae "
        "exponencialmente y el algoritmo pasa de explorar a afinar sin que "
        "haya que decirle cuando."
        % inspect.signature(sa.recocido).parameters["alfa"].default)]
    E += [FIG(R["fig_metropolis"],
              "Izquierda: la regla de aceptacion para tres temperaturas, con "
              "la T que esta corrida calibro de verdad. Centro: el "
              "enfriamiento para los tres valores de &alpha; que menciona el "
              "paper. Derecha: la aceptacion medida ciclo a ciclo en la "
              "corrida real, que arranca en %.0f%% y termina en %.0f%%."
              % (100 * res["historia"][0, 3], 100 * res["historia"][-1, 3]))]

    # ------------------------------------------------------------- 5 -----
    E += [P("5. El codigo, bloque por bloque", "h1")]
    E += [P(
        "Ahora si, <font face='%s'>analisis/sa.py</font> de arriba abajo. Son "
        "cuatro funciones: calibrar la temperatura, recocer, pulir y repetir "
        "desde varios lados. Ninguna pasa de cuarenta lineas." % FM)]

    E += [P("5.1 Elegir la temperatura inicial", "h2")]
    E += [P(
        "El paper dice T_init = 2000 y no explica de donde sale. El problema "
        "es que T se compara contra &Delta;E, o sea que <b>2000 no significa "
        "nada hasta saber en que unidades esta el costo</b>. Con un "
        "&chi;<sup>2</sup> que vale ~%.0f, T = 2000 anda bien; con un costo "
        "normalizado a 1, el mismo 2000 acepta absolutamente todo y la cadena "
        "se vuelve una caminata al azar que nunca converge."
        % R["chi2_semilla"])]
    E += [P(
        "La solucion del codigo es no elegir T sino <b>medirla</b>: se tiran "
        "unas cuantas propuestas al azar alrededor del punto de partida, se "
        "mira cuanto empeoran en promedio, y se despeja la T que hace que un "
        "empeoramiento tipico se acepte con la probabilidad que uno quiera.")]
    E += CODIGO(L_temp, sa)
    E += [N("<b>linea %d</b>: se evalua el costo en el punto de partida, para "
            "tener contra que comparar." % nl(L_temp, "c0 = costo(x0)"))]
    E += [N("<b>lineas %d a %d</b>: %s propuestas al azar; de cada una se "
            "guarda cuanto EMPEORO, si es que empeoro. Las que mejoran no "
            "sirven para esto: se aceptan siempre, sin importar T."
            % (nl(L_temp, "for _ in range(n)"), nl(L_temp, "subidas.append"),
               inspect.signature(sa.temperatura_inicial)
               .parameters["n"].default))]
    E += [N("<b>linea %d</b>: se despeja T de exp(-&lt;&Delta;E&gt;/T) = "
            "aceptacion. Con la aceptacion objetivo del %.0f%% que trae por "
            "defecto, esta corrida dio <b>T_init = %.4g</b>, en unidades de "
            "&chi;<sup>2</sup>."
            % (nl(L_temp, "return float(-np.mean"),
               100 * inspect.signature(sa.temperatura_inicial)
               .parameters["aceptacion"].default, res["T_init"]))]

    E += [P("5.2 Los argumentos del recocido", "h2")]
    E += [P(
        "La funcion principal tiene muchos argumentos, pero casi todos son "
        "constantes del algoritmo con un valor por defecto sensato. Los tres "
        "primeros son el problema; el resto es como buscarlo.")]
    filas = [["argumento", "por defecto", "que es"]]
    for nombre, par in inspect.signature(sa.recocido).parameters.items():
        pordef = ("requerido" if par.default is inspect.Parameter.empty
                  else repr(par.default))
        filas.append([CELDA("<font face='%s'>%s</font>" % (FM, nombre), tam=7.8),
                      CELDA("<font face='%s'>%s</font>" % (FM, pordef), tam=7.6),
                      CELDA(ARGUMENTOS.get(nombre, ""), tam=8.2)])
    E += [TABLA(filas, anchos=(3.5 * cm, 1.9 * cm, ANCHO - 5.4 * cm),
                alinear_num=False)]

    E += [P("5.3 La preparacion", "h2")]
    E += CODIGO(rec_prep, sa)
    E += [N("<b>linea %d</b>: el generador aleatorio. Al fijarle una semilla, "
            "la corrida entera es <b>reproducible</b>: los mismos numeros al "
            "azar, en el mismo orden, siempre. Es lo que permite volver a mirar "
            "una corrida que dio raro." % nl(rec_prep, "rng = np.random"))]
    E += [N("<b>lineas %d y %d</b>: el paso. Si no se lo dan, cada parametro "
            "se mueve como maximo un 10%% del ancho de sus cotas. Para esta "
            "corrida eso es %.4f en D y %.4f en &sigma;<sub>rel</sub>."
            % (nl(rec_prep, "if delta is None"),
               nl(rec_prep, "delta = 0.1 *"), R["delta"][0], R["delta"][1]))]
    E += [N("<b>linea %d</b>: la temperatura inicial, calibrada como se vio "
            "recien salvo que se la pasen a mano."
            % nl(rec_prep, "T = (temperatura_inicial"))]
    E += [N("<b>lineas %d a %d</b>: la contabilidad. <font face='%s'>x</font> "
            "y <font face='%s'>c</font> son donde esta la cadena y cuanto vale "
            "ahi; <font face='%s'>mejor_x</font> y "
            "<font face='%s'>mejor_c</font> son el mejor punto visto hasta el "
            "momento. Esa separacion es clave y se explica sola en cuanto uno "
            "ve el bucle: la cadena puede irse a un lugar peor, pero el mejor "
            "punto no se pierde nunca."
            % (nl(rec_prep, "c = costo(x)"),
               nl(rec_prep, "mejor_x, mejor_c"), FM, FM, FM, FM))]

    E += [P("5.4 El corazon: los dos bucles", "h2")]
    E += [P(
        "Estas dieciocho lineas son el algoritmo. El bucle de afuera cuenta "
        "<b>ciclos de temperatura</b>; el de adentro hace las %d propuestas "
        "que corresponden a esa temperatura."
        % inspect.signature(sa.recocido).parameters["n_markov"].default)]
    E += CODIGO(rec_bucle, sa)
    E += [N("<b>linea %d</b> - <b>la propuesta</b>. Se toma el punto actual y "
            "se le suma un desplazamiento al azar, uniforme entre "
            "-&delta; y +&delta; en cada parametro. Es un salto a ciegas: "
            "nadie mira la pendiente, no hay derivadas en ningun lado. Por eso "
            "el metodo funciona con costos que no son derivables."
            % nl(rec_bucle, "prop = x + rng.uniform"))]
    E += [N("<b>linea %d</b> - <b>las cotas</b>. "
            "<font face='%s'>np.clip</font> recorta la propuesta al rectangulo "
            "permitido. Es la unica barrera fisica del algoritmo: sin esto la "
            "busqueda se puede ir a un &sigma; negativo o a una D absurda."
            % (nl(rec_bucle, "prop = np.clip"), FM))]
    E += [N("<b>lineas %d y %d</b> - <b>cuanto cambia</b>. Se evalua el costo "
            "en la propuesta y se calcula &Delta;E. Negativo quiere decir que "
            "la propuesta es MEJOR."
            % (nl(rec_bucle, "cp = costo(prop)"), nl(rec_bucle, "dE = cp - c")))]
    E += [N("<b>linea %d</b> - <b>la regla de Metropolis</b>, la unica linea "
            "que distingue este metodo de una busqueda al azar. "
            "<font face='%s'>dE &lt;= 0</font> acepta toda mejora; si no, se "
            "sortea un numero entre 0 y 1 y se lo compara con exp(-&Delta;E/T). "
            "Python evalua el <font face='%s'>or</font> de izquierda a derecha "
            "y corta apenas puede, asi que <b>cuando la propuesta mejora ni "
            "siquiera se sortea</b>."
            % (nl(rec_bucle, "if dE <= 0 or"), FM, FM))]
    E += [N("<b>linea %d</b> - el <font face='%s'>max(T, 1e-300)</font> evita "
            "una division por cero cuando la temperatura ya se enfrio del "
            "todo. Detalle de implementacion, no de metodo."
            % (nl(rec_bucle, "if dE <= 0 or"), FM))]
    E += [N("<b>lineas %d y %d</b> - <b>mudarse</b>. Si se acepta, la cadena "
            "pasa a vivir en la propuesta. Notar que esto pasa tambien cuando "
            "el costo empeoro: es exactamente lo que se queria."
            % (nl(rec_bucle, "x, c = prop, cp"),
               nl(rec_bucle, "aceptados += 1")))]
    E += [N("<b>lineas %d y %d</b> - <b>el mejor no se pierde</b>. Antes de "
            "seguir se compara con el record y, si es mejor, se guarda una "
            "COPIA. La copia importa: sin "
            "<font face='%s'>.copy()</font> se guardaria una referencia al "
            "mismo array que el bucle sigue modificando."
            % (nl(rec_bucle, "if c < mejor_c"),
               nl(rec_bucle, "mejor_x, mejor_c = x.copy()"), FM))]
    E += [N("<b>lineas %d a %d</b> - <b>la libreta</b>. Al cerrar cada ciclo "
            "se anota la temperatura, el costo actual, el mejor costo, la "
            "fraccion aceptada, y ademas DONDE estaba la cadena. Esa ultima "
            "columna es la que despues dibuja el camino en la pestana 3: sin "
            "ella se sabria que el costo bajo, pero no si para eso hubo que "
            "cruzar de un valle a otro."
            % (nl(rec_bucle, "tasa = aceptados"),
               nl(rec_bucle, "historia_mejor_x.append")))]

    E += [P("5.5 Cuando parar", "h2")]
    E += [P(
        "Hay dos cortes y hacen falta los dos. El primero mira si el resultado "
        "dejo de mejorar; el segundo, si la temperatura ya es tan baja que no "
        "puede pasar nada nuevo.")]
    E += CODIGO(rec_cortes, sa)
    E += [N("<b>lineas %d a %d</b> - <b>corte por estabilizacion</b>. Si el "
            "mejor costo mejoro menos que <font face='%s'>tol_rel</font> "
            "respecto del ciclo anterior, se cuenta un ciclo estable; si "
            "mejoro, el contador vuelve a cero. Con %d ciclos estables "
            "seguidos, se corta. Es RELATIVO y no absoluto por lo mismo que la "
            "temperatura: un umbral absoluto depende de las unidades del "
            "costo, y el 0.001 del paper no sobrevive a cambiar de costo."
            % (nl(rec_cortes, "if len(historia) > 1"),
               nl(rec_cortes, "break"), FM,
               inspect.signature(sa.recocido).parameters["n_estable"].default))]
    E += [N("<b>lineas %d y %d</b> - <b>enfriar y corte por temperatura</b>. "
            "T se multiplica por &alpha; y, si cayo por debajo de una "
            "millonesima de la inicial, se termina. En esta corrida gano el "
            "corte por <b>%s</b>, a los <b>%d ciclos</b> y <b>%d "
            "evaluaciones</b> del costo, en <b>%.2f segundos</b>."
            % (nl(rec_cortes, "T *= alfa"), nl(rec_cortes, "if T < T_min_rel"),
               res["motivo"], res["ciclos"], res["n_eval"], R["t_sa"]))]

    E += [P("5.6 Lo que devuelve", "h2")]
    E += [P(
        "Devuelve <b>mejor_x</b>, no donde quedo la cadena. Es una distincion "
        "que se pasa por alto y cambia el resultado: la ultima propuesta "
        "aceptada puede ser peor que el record, porque el metodo acepta "
        "empeorar.")]
    E += CODIGO(rec_ret, sa)

    E += [P("5.7 Pulir", "h2")]
    E += [P(
        "El recocido encuentra el valle correcto pero se detiene con la "
        "precision de su ultimo paso, que no es cero. El pulido es un "
        "optimizador local corriente que baja al fondo exacto desde ahi. Los "
        "dos hacen falta: el recocido solo da los primeros decimales, el "
        "pulido solo da el valle equivocado. En esta corrida el pulido movio "
        "el &chi;<sup>2</sup> de %.4f a %.4f y D de %.5f a %.5f."
        % (R["chi2_sin_pulir"], R["chi2"], res["x"][0], x[0]))]
    E += CODIGO(L_pul, sa)
    E += [P(
        "Un detalle que conviene tener presente: <font face='%s'>minimize</font> "
        "corre <b>sin cotas</b>. El argumento <font face='%s'>cotas</font> de "
        "esta funcion se acepta y no se usa, asi que el pulido puede terminar "
        "afuera del rectangulo en el que el recocido tenia prohibido entrar. "
        "En la practica no molesta -el pulido arranca ya adentro del valle "
        "bueno y se mueve poco- pero es una asimetria real entre las dos "
        "etapas." % (FM, FM))]

    E += [P("5.8 Correr muchas veces: el multiarranque", "h2")]
    E += [P(
        "El recocido usa numeros al azar, asi que <b>una corrida no es un "
        "resultado, es una muestra</b>. Lo que se hace es correrlo N veces "
        "desde puntos iniciales <b>dispersos</b> y comparar. Si N busquedas "
        "independientes, arrancando de lugares distintos, terminan todas en el "
        "mismo punto, eso es evidencia practica de que ese punto es el minimo "
        "global. No es una demostracion, pero es lo mejor que se puede tener "
        "sin explorar todo el espacio.")]
    E += CODIGO(L_multi, sa)
    E += [N("<b>linea %d</b>: sin <font face='%s'>x0</font>, cada arranque "
            "sale de un punto <b>uniforme dentro de las cotas</b>. Es "
            "deliberado: si todos arrancaran de la semilla, que todos terminen "
            "igual no probaria nada."
            % (nl(L_multi, "xi = rng.uniform"), FM))]
    E += [N("<b>linea %d</b>: cada corrida recibe su propia semilla, sorteada "
            "por el generador de arriba. Asi las N son independientes entre si "
            "y el conjunto sigue siendo reproducible."
            % nl(L_multi, "r = recocido(costo, xi"))]
    E += [N("<b>linea %d</b>: cada corrida se pule antes de compararla, si no "
            "se estarian comparando resultados a medio terminar."
            % nl(L_multi, 'r["x_pulido"]'))]

    E += [P("5.9 De donde salen las barras de error", "h2")]
    E += [P(
        "El optimo por si solo no dice nada sobre su incerteza. Esa sale de la "
        "<b>curvatura</b> del costo en el minimo: si el valle es angosto, "
        "moverse poco empeora mucho el ajuste y el parametro esta bien "
        "determinado; si es ancho y plano, no. La funcion calcula esa "
        "curvatura numericamente -evaluando el costo alrededor del optimo- y "
        "la invierte para obtener la matriz de covarianza.")]
    E += CODIGO(recorte(L_err, desde="x = np.asarray(x, float)"), sa)
    E += [P(
        "Sobre esta corrida da <b>D = %.4f &plusmn; %.4f</b> (%.1f%%) y "
        "<b>&sigma;<sub>rel</sub> = %.4f &plusmn; %.4f</b>, con una "
        "correlacion entre los dos de <b>%+.2f</b>. Esa correlacion no es un "
        "detalle: dice que D y la resolucion se compensan entre si, asi que "
        "hay que reportarlas juntas."
        % (x[0], err[0], 100 * err[0] / x[0], x[1], err[1], R["corr_Dsigma"]))]

    # ------------------------------------------------------------- 6 -----
    E += [P("6. Una corrida de verdad, numero por numero", "h1")]
    E += [P(
        "Hasta aca el algoritmo en abstracto. Esta seccion es la misma corrida "
        "mirada por dentro: las primeras doce propuestas con todas sus cuentas "
        "a la vista. Los numeros no estan tipeados ni son un ejemplo inventado: "
        "salen de repetir las mismas lineas del bucle con la misma semilla del "
        "generador, y se verifican contra lo que devolvio el recocido "
        "(coincide el costo, la temperatura y la tasa de aceptacion al cerrar "
        "el primer ciclo).")]
    filas = [["", "propuesta (D, sigma)", "chi2", "dE", "P = exp(-dE/T)",
              "sorteo", "que hace"]]
    for f in R["traza"]:
        filas.append([
            "%d" % f["k"],
            "(%.4f, %.4f)" % (f["prop"][0], f["prop"][1]),
            "%.1f" % f["cp"],
            "%+.1f" % f["dE"],
            ("mejora" if f["dE"] <= 0 else "%.3f" % f["prob"]),
            ("no sortea" if f["sorteo"] is None else "%.3f" % f["sorteo"]),
            ("acepta" if f["acepta"] else "rechaza")])
    E += [TABLA(filas, anchos=(0.8 * cm, 3.7 * cm, 2.0 * cm, 2.0 * cm,
                               2.6 * cm, 2.1 * cm, ANCHO - 13.2 * cm),
                cuerpo=8.0, izquierda=(4, 5, 6))]
    E += [P(
        "Vale la pena leer tres filas de esa tabla:", "p")]
    E += [B("una propuesta que <b>mejora</b> se acepta sin sortear nada, y en "
            "la columna del sorteo dice justamente eso;")]
    E += [B("una propuesta que <b>empeora</b> tiene una probabilidad calculada "
            "y se la compara con un numero al azar. A esta temperatura "
            "(T = %.4g, del orden del propio &chi;<sup>2</sup>) esa "
            "probabilidad es alta y varias subidas pasan;"
            % res["T_init"])]
    E += [B("las que se rechazan no mueven la cadena: el punto actual sigue "
            "siendo el de antes, y en la fila siguiente la propuesta se "
            "construye otra vez desde ahi.")]
    E += [P(
        "Sumadas las %d propuestas del primer ciclo, se aceptaron el "
        "<b>%.0f%%</b>. Esa es la fraccion de aceptacion que la interfaz "
        "grafica, y es el termometro del algoritmo: alta al principio, casi "
        "cero al final."
        % (inspect.signature(sa.recocido).parameters["n_markov"].default,
           100 * res["historia"][0, 3]))]

    E += [P("6.1 Los ciclos, de principio a fin", "h2")]
    filas = [["ciclo", "T", "costo actual", "mejor costo", "aceptacion"]]
    for i, fila in enumerate(res["historia"]):
        filas.append(["%d" % (i + 1), "%.4g" % fila[0], "%.2f" % fila[1],
                      "%.2f" % fila[2], "%.0f %%" % (100 * fila[3])])
    E += [TABLA(filas, anchos=(2.0 * cm, 3.0 * cm, 3.4 * cm, 3.4 * cm,
                               ANCHO - 11.8 * cm), cuerpo=8.2)]
    E += [P(
        "Se ve el algoritmo entero en una tabla: la temperatura baja un %.0f%% "
        "por fila, la aceptacion cae con ella, el costo actual deja de saltar "
        "y el mejor costo se aplana. Cuando se aplana durante %d ciclos "
        "seguidos, corta."
        % (100 * (1 - inspect.signature(sa.recocido).parameters["alfa"].default),
           inspect.signature(sa.recocido).parameters["n_estable"].default))]
    E += [FIG(R["fig_corrida"],
              "Izquierda: el costo ciclo a ciclo. La linea fina es donde esta "
              "la cadena -sube y baja, porque acepta empeorar- y la gruesa es "
              "el record, que solo puede bajar. Centro: el paisaje real del "
              "&chi;<sup>2</sup> y el camino que recorrio la cadena; las "
              "lineas blancas son los contornos de 1, 2 y 3 sigma. Derecha: "
              "que pasa con la aceptacion si la temperatura inicial esta cien "
              "veces mal en cada direccion.")]

    # ------------------------------------------------------------- 7 -----
    E += [P("7. Como se sabe que ese es el minimo", "h1")]
    E += [P(
        "Un optimizador siempre devuelve algo. La pregunta es si ese algo es "
        "el minimo global o el fondo del primer pozo que encontro, y hay tres "
        "controles independientes.")]

    E += [P("7.1 Multiarranque", "h2")]
    E += [P(
        "%d recocidos desde puntos al azar dentro de las cotas. Los %d "
        "terminan en el mismo &chi;<sup>2</sup> (dentro de 0.1) y la "
        "dispersion de las D finales es <b>%.1e</b> mV&middot;ns/keVee. Para "
        "comparar, la barra de error estadistica de D es %.4f, o sea "
        "<b>%.0e veces mas grande</b>: la incerteza del resultado no viene "
        "del optimizador."
        % (len(ms), R["ms_iguales"], R["ms_dispersion"], err[0],
           err[0] / max(R["ms_dispersion"], 1e-12)))]
    E += [FIG(R["fig_multiarranque"],
              "Izquierda: los %d arranques (cuadrados blancos), sus caminos y "
              "sus finales (circulos rojos), sobre el paisaje. Derecha: la D "
              "final contra la D del punto de partida. Que la nube sea plana "
              "es el resultado: de donde se arranque no cambia donde se "
              "termina." % len(ms))]

    E += [P("7.2 Contra un optimizador local", "h2")]
    E += [P(
        "El otro control es comparar con el metodo simple. Se resolvio el "
        "mismo problema desde <b>%d puntos de partida al azar</b>, una vez con "
        "un optimizador local corriente (Nelder-Mead, el mismo que usa el "
        "pulido) y otra con el recocido:" % R["local"]["n"])]
    E += [TABLA(
        [["metodo", "llegan al minimo", "peor resultado"],
         [CELDA("optimizador local, sin cotas"),
          CELDA("%d de %d" % (R["local"]["loc_ok"], R["local"]["n"])),
          CELDA("D = %.3f con &chi;<sup>2</sup> = %.0f"
                % (R["local"]["loc_peor"]["loc_x"][0],
                   R["local"]["loc_peor"]["loc_chi2"]))],
         [CELDA("recocido + pulido"),
          CELDA("%d de %d" % (R["local"]["sa_ok"], R["local"]["n"])),
          CELDA("-")]],
        anchos=(6.4 * cm, 4.2 * cm, ANCHO - 10.6 * cm), alinear_num=False)]
    E += [P(
        "El resultado hay que leerlo con honestidad: en <b>este</b> paisaje el "
        "optimizador local casi siempre acierta, porque la region de interes "
        "tiene un solo valle. O sea que el recocido <b>no cambia el numero</b>, "
        "y eso es justamente lo que se queria saber. Lo que aporta es no "
        "depender de la suerte: el caso que fallo se escapo a D = %.2f con un "
        "&chi;<sup>2</sup> %.0f veces peor, y nada en el codigo lo hubiera "
        "marcado como sospechoso si se corriera una sola vez."
        % (R["local"]["loc_peor"]["loc_x"][0],
           R["local"]["loc_peor"]["loc_chi2"] / R["chi2"]))]
    E += [P(
        "Notar la asimetria: el optimizador local corre sin cotas y el "
        "recocido recorta cada propuesta. Parte de la ventaja es del "
        "algoritmo y parte de haberle dicho donde puede buscar.", "nota")]

    E += [P("7.3 El paisaje", "h2")]
    E += [P(
        "El tercer control es mirar. El panel del centro de la figura anterior "
        "es el &chi;<sup>2</sup> en todo el rectangulo de busqueda: si "
        "hubiera un segundo valle comparable, se veria. Y el minimo tiene que "
        "verse aproximadamente <b>cuadratico</b>, porque de esa suposicion "
        "salen las barras de error de la seccion 5.9. Si el valle fuera "
        "torcido o tuviera una meseta, la curvatura local mentiria y habria "
        "que reportar el error de otra manera.")]

    # ------------------------------------------------------------- 8 -----
    E += [P("8. Como lo llama la interfaz", "h1")]
    E += [P(
        "El modulo <font face='%s'>sa.py</font> no sabe nada de bordes "
        "Compton. Quien lo conecta con la fisica es "
        "<font face='%s'>calibracion.calibrar</font>, que decide las cotas, "
        "corre el multiarranque y se queda con la mejor corrida." % (FM, FM))]
    E += CODIGO(cal_frag, calibracion)
    E += [N("<b>linea %d</b> - <b>las cotas</b>, y de donde salen. D se "
            "busca entre el 60%% y el 160%% de la semilla del NCC, que es una "
            "ventana generosa alrededor de una estimacion que ya es buena. "
            "&sigma;<sub>rel</sub> entre 0.02 y 0.30: menos del 2%% seria un "
            "detector que no borronea nada y mas del 30%% no seria un borde."
            % nl(cal_frag, "cotas = ["))]
    E += [N("<b>lineas %d a %d</b> - si el modelo incluye el exceso de bajo "
            "canal, aparece un tercer parametro (&lambda;) con su propia cota "
            "y su propia semilla. El algoritmo no cambia: busca en tres "
            "dimensiones en vez de dos."
            % (nl(cal_frag, "if exceso:"), nl(cal_frag, "semilla = tuple(semilla) + (0.15")))]
    E += [N("<b>linea %d</b> - el multiarranque, con <b>%d</b> corridas por "
            "defecto."
            % (nl(cal_frag, "corridas = sa.multiarranque"),
               inspect.signature(calibracion.calibrar)
               .parameters["n_arranques"].default))]
    E += [N("<b>lineas %d a %d</b> - se elige la corrida de menor "
            "&chi;<sup>2</sup>. Las otras no se tiran: su dispersion se "
            "reporta como diagnostico."
            % (nl(cal_frag, "j = int(np.argmin"), nl(cal_frag, "Ds_ms = ")))]

    E += [P("8.1 El detalle que la interfaz necesita para dibujar", "h2")]
    E += [P(
        "Un ajuste devuelve numeros; una pestana que explica el ajuste "
        "necesita ademas el <b>camino</b>. Por eso <font face='%s'>calibrar"
        "</font> acepta un diccionario opcional y lo llena con la semilla, "
        "las cotas y las N corridas completas. No entra en el resultado "
        "normal a proposito: el resultado se guarda como JSON y esto son "
        "arrays de numpy." % FM)]
    E += CODIGO(pipe_frag, pipeline)
    E += [P(
        "Con eso, la pestana <b>3 &middot; Semilla y recocido</b> puede "
        "mostrar todo sin volver a correr nada:")]
    E += [TABLA(
        [["panel", "que muestra", "que hay que mirar"]] +
        [[CELDA("<b>%s</b>" % a), CELDA(b), CELDA(c)]
         for a, b, c in PANELES_GUI],
        anchos=(3.2 * cm, 4.6 * cm, ANCHO - 7.8 * cm), alinear_num=False)]

    E += [P("8.2 El resultado, con el modelo completo", "h2")]
    E += [P(
        "Corrido tal como lo corre la interfaz -o sea con el exceso de bajo "
        "canal adentro, que agrega el tercer parametro- el ajuste de esta "
        "corrida da:")]
    E += [TABLA(
        [["", "valor", "que es"],
         [CELDA("<b>D</b>"), CELDA("%.4f &plusmn; %.4f" % (rc["D"], rc["err_D"])),
          CELDA("mV&middot;ns por keVee: la calibracion buscada")],
         [CELDA("<b>&sigma;<sub>rel</sub></b>"),
          CELDA("%.4f &plusmn; %.4f" % (rc["sigma_rel"], rc["err_sigma_rel"])),
          CELDA("resolucion, como fraccion del canal del borde "
                "(%.1f %% FWHM)" % (100 * rc["fwhm_rel"]))],
         [CELDA("<b>&chi;<sup>2</sup>/ndf</b>"),
          CELDA("%.3f" % rc["chi2_ndf"]),
          CELDA("%d grados de libertad" % rc["ndf"])],
         [CELDA("<b>bordes</b>"),
          CELDA(", ".join("%.0f" % b for b in rc["bordes_ch"])),
          CELDA("donde caen, en mV&middot;ns, los bordes Compton de "
                "%s" % rc["fuente"])],
         [CELDA("<b>dispersion</b>"),
          CELDA("%.1e" % rc["dispersion_multiarranque"]),
          CELDA("entre los %d arranques: el optimizador no aporta incerteza"
                % rc["n_arranques"])]],
        anchos=(2.6 * cm, 4.0 * cm, ANCHO - 6.6 * cm), alinear_num=False)]
    E += [P(
        "Tardo <b>%.1f segundos</b> en total, de los cuales el recocido son "
        "unas decimas: lo caro de la cadena esta antes, en leer el archivo y "
        "armar el espectro." % R["t_cal"])]

    # ------------------------------------------------------------- 9 -----
    E += [P("9. Que pasa si se mueve cada constante", "h1")]
    E += [P(
        "La forma honesta de justificar una constante es moverla y ver que "
        "cambia. Cada fila es el mismo problema resuelto %d veces con semillas "
        "distintas, cambiando una sola cosa respecto del default."
        % R["variantes"][0]["n"])]
    filas = [["variante", "D (mediana)", "dispersion", "ciclos",
              "evaluaciones"]]
    for v in R["variantes"]:
        filas.append([CELDA(v["nombre"], tam=8.2), "%.4f" % v["D"],
                      "%.1e" % v["sd"], "%.0f" % v["ciclos"],
                      "%.0f" % v["n_eval"]])
    E += [TABLA(filas, anchos=(7.4 * cm, 2.4 * cm, 2.2 * cm, 1.6 * cm,
                               ANCHO - 13.6 * cm), cuerpo=8.2)]
    E += [P(
        "La lectura: sobre este problema <b>ninguna de las constantes cambia "
        "el resultado</b>, solo cuanto tarda. Eso es lo que uno quiere de un "
        "optimizador -que el numero sea del problema y no del metodo- y "
        "tambien es la razon por la que no hay que gastar tiempo afinandolas. "
        "Lo que si cambia es el precio: enfriar mas lento o proponer poco "
        "cuesta mas evaluaciones para llegar al mismo lugar.")]
    E += [P(
        "El caso interesante es el de las <b>10 propuestas por temperatura</b> "
        "en vez de 100: llega al mismo sitio con una decima parte del trabajo. "
        "En un paisaje con varios valles no seria asi -la cadena necesita "
        "propuestas suficientes por temperatura para explorar antes de que se "
        "enfrie- y por eso el default sigue siendo el del paper.")]

    E += [P("9.1 Lo unico que si rompe: cotas que dejan el optimo afuera",
            "h2")]
    E += [P(
        "Se corrio a proposito con las cotas de D puestas en [%.2f, %.2f], que "
        "no contienen el optimo (%.4f). El recocido, que recorta cada "
        "propuesta, no puede hacer otra cosa que apoyarse contra la pared. Lo "
        "llamativo es lo que pasa despues: el <b>pulido corre sin cotas</b>, "
        "se escapa de la caja y termina en D = %.4f, el valor correcto. "
        "Salio bien por casualidad, y muestra el punto: el algoritmo no tiene "
        "forma de avisar que las cotas estaban mal, porque un numero devuelto "
        "contra la pared se ve igual que cualquier otro."
        % (R["encerrado"]["cotas"][0, 0], R["encerrado"]["cotas"][0, 1],
           x[0], R["encerrado"]["x"][0]))]
    E += [P(
        "La defensa contra eso no esta en el optimizador sino aguas arriba: "
        "las cotas se derivan de la semilla del NCC -del 60 al 160 por ciento "
        "de ella- y no se eligen a mano, y el multiarranque cubre todo el "
        "rectangulo. Si el "
        "optimo estuviera contra un borde, los N arranques terminarian todos "
        "ahi y se veria en el panel del camino.")]

    # ------------------------------------------------------------ 10 -----
    E += [P("10. Preguntas que aparecen la primera vez", "h1")]
    for preg, resp in [
        ("Si acepta cosas peores, como puede terminar en un lugar bueno?",
         "Porque el mejor punto se guarda aparte y no se pierde nunca "
         "(<font face='%s'>mejor_x</font>). La cadena pasea; el resultado es "
         "el record. Y ademas la probabilidad de aceptar algo peor se "
         "desploma a medida que baja la temperatura, asi que el paseo se "
         "apaga solo." % FM),
        ("El resultado es aleatorio?",
         "La busqueda si; el resultado, en la practica, no. Con la semilla "
         "fija la corrida es reproducible bit a bit, y el multiarranque mide "
         "cuanto cambia el numero entre corridas independientes: %.1e, o sea "
         "%.0f ordenes de magnitud menos que la barra de error."
         % (R["ms_dispersion"],
            np.log10(err[0] / max(R["ms_dispersion"], 1e-12)))),
        ("Garantiza encontrar el minimo global?",
         "No. El teorema que garantiza convergencia pide un enfriamiento "
         "logaritmico (T ~ c/log k), que es inaplicablemente lento. Con "
         "enfriamiento geometrico, que es lo que usa todo el mundo y lo que "
         "dice el paper, no hay garantia: hay evidencia, y es la de la "
         "seccion 7."),
        ("Por que no usar una funcion de scipy?",
         "Porque el paper especifica este algoritmo con estas constantes, y el "
         "objetivo era replicarlo y poder auditarlo. "
         "<font face='%s'>scipy.optimize.dual_annealing</font> existe y "
         "resuelve un problema parecido, pero es otro algoritmo (recocido "
         "generalizado, con saltos de Cauchy y pulido local automatico): "
         "usarlo seria comparar contra otra cosa." % FM),
        ("Cuanto tarda?",
         "Una corrida, %.2f segundos (%d evaluaciones del costo a %.0f "
         "microsegundos cada una). El multiarranque de %d, %.1f segundos. Es "
         "despreciable al lado de leer el archivo del digitalizador."
         % (R["t_sa"], res["n_eval"], R["t_costo"] * 1e6, len(ms), R["t_ms"])),
        ("Y si el ajuste da mal?",
         "Antes de tocar el optimizador conviene mirar las cuatro cosas de la "
         "pestana 3, en este orden: la aceptacion (si arranca en 100% o en 0%, "
         "la temperatura esta mal), la dispersion entre arranques (si es "
         "grande, hay mas de un valle), el camino (si termina contra una "
         "pared, las cotas cortan), y el paisaje (si el valle no es un pozo "
         "cerrado, las barras de error mienten)."),
    ]:
        E += [P("<b>%s</b>" % preg, "h2"), P(resp)]

    # ------------------------------------------------------------ 11 -----
    E += [P("11. Donde esta cada cosa", "h1")]
    filas = [["funcion", "archivo y lineas", "que hace"]]
    for obj, que in [
        (modelo.armar, "arma el problema: espectro + plantilla + ventana"),
        (modelo.Problema.costo, "el numero que el optimizador minimiza"),
        (modelo.amplitud_optima, "despeja la amplitud sin buscarla (VARPRO)"),
        (sa.temperatura_inicial, "calibra T_init midiendo el costo"),
        (sa.recocido, "el algoritmo: Metropolis + enfriamiento"),
        (sa.pulir, "optimizador local que baja al fondo exacto"),
        (sa.multiarranque, "N recocidos desde puntos dispersos"),
        (sa.errores_por_curvatura, "barras de error por curvatura del chi2"),
        (calibracion.calibrar, "la receta: cotas, multiarranque, resultado"),
        (pipeline.calibrar_completo, "la cadena entera, que llaman las dos "
                                     "interfaces"),
        (gui.App._fig_recocido, "los cuatro paneles de la pestana 3"),
    ]:
        a, b = rango(lineas_de(obj, sin_doc=False))
        filas.append([
            CELDA("<font face='%s'>%s</font>" % (FM, obj.__qualname__), tam=8.0),
            CELDA("%s: %d-%d" % (archivo_de(obj), a, b), tam=8.0),
            CELDA(que, tam=8.2)])
    E += [TABLA(filas, anchos=(5.6 * cm, 4.0 * cm, ANCHO - 9.6 * cm),
                alinear_num=False)]
    E += [P(
        "Para verlo corriendo: <font face='%s'>python analisis/sa.py</font> "
        "deja una figura con el algoritmo sobre esta misma corrida, "
        "<font face='%s'>python analisis/pipeline.py co60</font> hace la "
        "cadena completa por consola, y "
        "<font face='%s'>python analisis/gui.py</font> abre la interfaz."
        % (FM, FM, FM))]

    # --------------------------------------------------------- apendice --
    E += [P("Apendice. El modulo completo", "h1")]
    E += [P(
        "Lo que sigue es <font face='%s'>analisis/sa.py</font> entero -sin la "
        "demo del final- con sus comentarios y docstrings, que es donde el "
        "codigo explica sus propias decisiones. Los recuadros de las secciones "
        "anteriores son fragmentos de esto." % FM)]
    E += CODIGO(recorte(lineas_modulo(sa), desde="import os",
                        hasta="# demo"), sa, "el modulo, sin la demo")

    doc = SimpleDocTemplate(SALIDA, pagesize=A4,
                            leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=2 * cm, bottomMargin=2 * cm,
                            title="El recocido simulado, linea por linea",
                            author="Calibracion por borde Compton")
    doc.build(E, onFirstPage=pie, onLaterPages=pie)
    return SALIDA


# ==========================================================================
# 5. MAIN
# ==========================================================================
def main():
    rutas.salida_utf8()
    rapido = "--rapido" in sys.argv
    t0 = time.time()
    R = calcular(rapido=rapido)

    print("dibujando...")
    R["fig_cadena"] = fig_cadena()
    R["fig_juguete"] = fig_juguete()
    R["fig_costo"] = fig_costo(R)
    R["fig_metropolis"] = fig_metropolis(R)
    R["fig_corrida"] = fig_corrida(R)
    R["fig_multiarranque"] = fig_multiarranque(R)

    print("armando el PDF...")
    salida = construir(R)
    print()
    print("  corrida        : %s (%s, %s)"
          % (R["corrida"].clave, R["corrida"].fuente, R["corrida"].punto()))
    print("  modelo desnudo : D = %.4f +- %.4f, chi2/ndf = %.3f"
          % (R["x"][0], R["err"][0], R["chi2_ndf"]))
    print("  modelo completo: D = %.4f +- %.4f, chi2/ndf = %.3f"
          % (R["r_cal"]["D"], R["r_cal"]["err_D"], R["r_cal"]["chi2_ndf"]))
    print("  traza verificada contra el recocido: %s"
          % ("SI" if R["traza_ok"] else "NO"))
    print("  %.1f s en total" % (time.time() - t0))
    print()
    print("PDF -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    main()
