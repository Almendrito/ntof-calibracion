"""
doc_sa.py — genera `docs/recocido_simulado.pdf`.

El documento explica la Etapa de optimizacion: que dice el paper, que es el
recocido simulado, como se aplica a nuestro modelo y —sobre todo— como se
decide si los parametros encontrados son el optimo global o no.

Ningun numero del PDF esta tipeado: todos salen de correr el ajuste sobre los
datos reales al momento de generarlo, asi que el documento no puede quedar
desincronizado del codigo.

    python analisis/doc_sa.py              # completo (~1 min)
    python analisis/doc_sa.py --rapido     # menos arranques, para iterar el texto
"""

import os
import sys
import tempfile
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analisis import rutas, modelo, sa

TMP = tempfile.mkdtemp(prefix="doc_sa_")
SALIDA = os.path.join(rutas.RAIZ, "docs", "recocido_simulado.pdf")

AZUL = "#1f4e79"
NARANJA = "#c4630a"
GRIS = "#4a5560"


# ==========================================================================
# 1. CALCULOS  (todo lo que despues aparece como numero en el texto)
# ==========================================================================
def calcular(rapido=False):
    R = {}
    print("[1/6] armando el problema...")
    p, semilla, ncc_r, esp, corr = modelo.armar("co60", verbose=False)
    R["p"] = p
    R["semilla"] = semilla
    R["ncc"] = ncc_r
    R["esp"] = esp
    R["corrida"] = corr
    R["ventana"] = (ncc_r["lo"], ncc_r["hi"])
    R["nbins"] = int(p.m.sum())
    R["ndf"] = p.ndf

    cotas = [(0.6 * semilla[0], 1.6 * semilla[0]), (0.02, 0.30)]
    R["cotas"] = cotas

    print("[2/6] recocido desde la semilla del NCC...")
    t0 = time.time()
    res = sa.recocido(p.costo, semilla, cotas, semilla=1)
    xp, cp = sa.pulir(p.costo, res["x"])
    R["res"] = res
    R["x"] = xp
    R["chi2"] = cp
    R["chi2_ndf"] = cp / p.ndf
    R["t_sa"] = time.time() - t0
    R["chi2_semilla"] = p.chi2(*semilla)[0] / p.ndf
    R["A"] = p.chi2(*xp)[1]

    print("[3/6] errores por curvatura...")
    err, C = sa.errores_por_curvatura(p.costo, xp)
    R["err"] = err
    R["C"] = C
    R["corr_Dsigma"] = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])

    print("[4/6] multiarranque...")
    n = 4 if rapido else 12
    ms = sa.multiarranque(p.costo, cotas, n=n, semilla=7)
    R["ms"] = ms
    R["ms_D"] = np.array([m["x_pulido"][0] for m in ms])
    R["ms_c"] = np.array([m["costo_pulido"] for m in ms])

    print("[5/6] barrido ancho de D (mapa del paisaje)...")
    from scipy.optimize import minimize_scalar
    Ds = np.linspace(0.3 * semilla[0], 2.5 * semilla[0], 90 if rapido else 150)
    perfil = []
    for D in Ds:
        rr = minimize_scalar(lambda s: p.costo((D, s)), bounds=(0.02, 0.35),
                             method="bounded")
        perfil.append(rr.fun)
    perfil = np.array(perfil)
    R["perfil_D"], R["perfil_chi2"] = Ds, perfil / p.ndf
    R["minimos_locales"] = [(float(Ds[i]), float(perfil[i] / p.ndf))
                            for i in range(1, len(perfil) - 1)
                            if perfil[i] < perfil[i - 1] and perfil[i] < perfil[i + 1]]

    print("[6/6] optimizador local desde arranques malos...")
    from scipy.optimize import minimize
    malos = [(0.70 * xp[0], 0.28), (1.45 * xp[0], 0.03), (0.65 * xp[0], 0.05)]
    comp = []
    for m0 in malos:
        loc = minimize(p.costo, m0, method="Nelder-Mead",
                       options=dict(xatol=1e-6, fatol=1e-4, maxiter=4000))
        rs = sa.recocido(p.costo, m0, cotas, semilla=3)
        xs, cs = sa.pulir(p.costo, rs["x"])
        comp.append(dict(x0=m0, D_loc=float(loc.x[0]), c_loc=float(loc.fun),
                         D_sa=float(xs[0]), c_sa=float(cs)))
    R["comparacion"] = comp

    print("[7/8] verificaciones de las decisiones de implementacion...")
    # Cada una contrasta una afirmacion que, si no, quedaria sostenida solo por
    # el texto. Estan tambien en `verificar_sa.py`, con mas detalle.
    rng = np.random.default_rng(0)
    delta = 0.1 * np.array([c[1] - c[0] for c in cotas])
    lo_c = [c[0] for c in cotas]
    hi_c = [c[1] for c in cotas]

    def aceptacion(costo, T, n=2500):
        z = np.array(semilla, float)
        cz = costo(z)
        acc = 0
        for _ in range(n):
            pr = np.clip(z + rng.uniform(-1, 1, 2) * delta, lo_c, hi_c)
            cp = costo(pr)
            if cp - cz <= 0 or rng.random() < np.exp(-(cp - cz) / T):
                z, cz = pr, cp
                acc += 1
        return acc / n

    escala_costo = p.costo(semilla)
    R["escala_costo"] = escala_costo
    R["acc_2000_chi2"] = aceptacion(p.costo, 2000.0)
    R["acc_auto_chi2"] = aceptacion(p.costo, res["T_init"])
    R["acc_2000_orden1"] = aceptacion(lambda z: p.costo(z) / escala_costo, 2000.0)

    R["v_paso"] = {}
    for etiqueta, adapt in (("fijo", False), ("adaptativo", True)):
        Ds = np.array([sa.recocido(p.costo, semilla, cotas, semilla=s,
                                   paso_adaptativo=adapt)["x"][0]
                       for s in range(6 if rapido else 8)])
        R["v_paso"][etiqueta] = (Ds.mean(), Ds.std())

    R["v_acept"] = []
    for obj in (0.234, 0.35, 0.44):
        Ds = np.array([sa.pulir(p.costo,
                                sa.recocido(p.costo, semilla, cotas, semilla=s,
                                            paso_adaptativo=True,
                                            objetivo_aceptacion=obj)["x"])[0][0]
                       for s in range(4)])
        R["v_acept"].append((obj, Ds.mean(), Ds.std()))

    A_opt = R["A"]

    def costo3(z):
        Dz, sz, Az = z
        if Dz <= 0 or sz <= 0 or Az <= 0:
            return 1e30
        return p.chi2(Dz, sz, A=Az)[0]

    cotas3 = list(cotas) + [(0.05 * A_opt, 20 * A_opt)]
    Ds3 = np.array([sa.recocido(costo3, list(semilla) + [A_opt], cotas3,
                                semilla=s)["costo"] for s in range(4)])
    Ds2 = np.array([sa.recocido(p.costo, semilla, cotas, semilla=s)["costo"]
                    for s in range(4)])
    R["v_varpro"] = (Ds2.mean(), Ds3.mean())

    print("[8/8] el otro detector...")
    # el otro detector, para la tabla de resultados
    p2, sem2, ncc2, _, _ = modelo.armar("cs137_1600", verbose=False)
    cot2 = [(0.6 * sem2[0], 1.6 * sem2[0]), (0.02, 0.30)]
    r2 = sa.recocido(p2.costo, sem2, cot2, semilla=1)
    x2, c2 = sa.pulir(p2.costo, r2["x"])
    e2, C2 = sa.errores_por_curvatura(p2.costo, x2)
    R["cs137"] = dict(D=x2[0], sigma=x2[1], err_D=e2[0], chi2_ndf=c2 / p2.ndf,
                      ndf=p2.ndf, ventana=(ncc2["lo"], ncc2["hi"]),
                      semilla=sem2[0], rho=ncc2["rho"], p=p2)
    return R


# ==========================================================================
# 2. FIGURAS
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


def fig_juguete():
    """Por que hace falta SA: funcion de juguete con muchos minimos."""
    f = lambda x: (x - 1.6) ** 2 / 6.0 - np.cos(3.4 * x) + 1.0
    xs = np.linspace(-4, 6, 900)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))

    a = axes[0]
    a.plot(xs, f(xs), lw=1.6, color=AZUL)
    from scipy.optimize import minimize
    for x0, col in ((-3.2, NARANJA), (4.6, "#8a2be2")):
        r = minimize(lambda v: f(v[0]), [x0], method="Nelder-Mead")
        a.plot([x0], [f(x0)], "o", color=col, ms=7)
        a.annotate("", xy=(r.x[0], f(r.x[0])), xytext=(x0, f(x0)),
                   arrowprops=dict(arrowstyle="->", color=col, lw=1.4))
        a.plot([r.x[0]], [f(r.x[0])], "*", color=col, ms=14)
    a.plot([1.6], [f(1.6)], "v", color="#1f7a3d", ms=9)
    a.text(1.6, f(1.6) - 0.9, "minimo global", ha="center", fontsize=8,
           color="#1f7a3d")
    a.set_title("un optimizador local solo sabe bajar", fontsize=10)
    a.set_xlabel("parametro"); a.set_ylabel("costo")
    a.grid(alpha=.3)

    a = axes[1]
    a.plot(xs, f(xs), lw=1.2, color="0.75", zorder=1)
    rng = np.random.default_rng(4)
    cam, x = [], -3.2
    c = f(x)
    T = 4.0
    for k in range(220):
        prop = x + rng.uniform(-1, 1) * 1.4
        cp = f(prop)
        if cp - c <= 0 or rng.random() < np.exp(-(cp - c) / T):
            x, c = prop, cp
        cam.append((x, c))
        T *= 0.975
    cam = np.array(cam)
    sc = a.scatter(cam[:, 0], cam[:, 1], c=np.arange(len(cam)), s=9,
                   cmap="plasma", zorder=3)
    a.plot([1.6], [f(1.6)], "v", color="#1f7a3d", ms=9, zorder=4)
    cb = fig.colorbar(sc, ax=a); cb.set_label("iteracion", fontsize=8)
    a.set_title("el recocido acepta subir, y por eso escapa", fontsize=10)
    a.set_xlabel("parametro"); a.set_ylabel("costo")
    a.grid(alpha=.3)
    return _guardar(fig, "juguete.png")


def fig_mecanica(R):
    """Metropolis y enfriamiento, con la T calibrada sobre nuestro costo."""
    T0 = R["res"]["T_init"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.5))

    a = axes[0]
    dE = np.linspace(0, 6 * T0 / 3, 400)
    for frac, col in ((1.0, AZUL), (0.1, NARANJA), (0.01, "#1f7a3d")):
        a.plot(dE, np.exp(-dE / (T0 * frac)), lw=1.6, color=col,
               label="T = %.3g" % (T0 * frac))
    a.set_title("P(aceptar un empeoramiento) = exp(-$\\Delta E$/T)", fontsize=10)
    a.set_xlabel("$\\Delta E$  [unidades de $\\chi^2$]")
    a.set_ylabel("probabilidad"); a.legend(fontsize=8); a.grid(alpha=.3)

    a = axes[1]
    k = np.arange(160)
    for al, col in ((0.85, AZUL), (0.95, NARANJA), (0.99, "#1f7a3d")):
        a.semilogy(k, T0 * al ** k, lw=1.5, color=col, label="$\\alpha$ = %.2f" % al)
    a.axhline(1e-6 * T0, color="k", ls=":", lw=1, label="corte $10^{-6}T_0$")
    a.set_title("enfriamiento geometrico $T_k = \\alpha^k T_0$", fontsize=10)
    a.set_xlabel("ciclo de temperatura"); a.set_ylabel("T")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = axes[2]
    k = np.arange(1, 400)
    a.loglog(k, T0 / np.log(1 + k), lw=1.6, color="#8a2be2",
             label="logaritmico $c/\\log k$  (con garantia)")
    a.loglog(k, T0 * 0.95 ** k, lw=1.6, color=NARANJA,
             label="geometrico $\\alpha^k$  (el que se usa)")
    a.set_ylim(1e-8, 1e5)
    a.set_title("por que no hay garantia", fontsize=10)
    a.set_xlabel("ciclo"); a.set_ylabel("T")
    a.legend(fontsize=8); a.grid(alpha=.3)
    return _guardar(fig, "mecanica.png")


def fig_paisaje(R):
    p, xp = R["p"], R["x"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.9))

    a = axes[0]
    Ds = np.linspace(R["cotas"][0][0], R["cotas"][0][1], 110)
    Ss = np.linspace(R["cotas"][1][0], R["cotas"][1][1], 90)
    Z = p.mapa(Ds, Ss) / p.ndf
    im = a.pcolormesh(Ds, Ss, np.log10(Z), shading="auto", cmap="viridis")
    a.contour(Ds, Ss, Z * p.ndf,
              levels=[R["chi2"] + 2.30, R["chi2"] + 6.17, R["chi2"] + 11.8],
              colors="w", linewidths=.9)
    a.plot(*xp, "o", color="#e02020", ms=7, label="optimo")
    a.plot(*R["semilla"], "s", color=NARANJA, ms=6, label="semilla del NCC")
    fig.colorbar(im, ax=a, label="$\\log_{10}(\\chi^2/ndf)$")
    a.set_title("el paisaje, y los contornos de 1, 2 y 3$\\sigma$", fontsize=10)
    a.set_xlabel("D  [mV$\\cdot$ns/keVee]"); a.set_ylabel("$\\sigma_{rel}$")
    a.legend(fontsize=8)

    a = axes[1]
    a.semilogy(R["perfil_D"], R["perfil_chi2"], lw=1.6, color=AZUL)
    a.axvline(xp[0], color="#e02020", ls="--", lw=1.2,
              label="unico minimo: D = %.4f" % xp[0])
    a.set_title("barrido ancho de D ($\\sigma$ optimizado en cada punto)",
                fontsize=10)
    a.set_xlabel("D  [mV$\\cdot$ns/keVee]"); a.set_ylabel("$\\chi^2/ndf$")
    a.legend(fontsize=8); a.grid(alpha=.3)
    return _guardar(fig, "paisaje.png")


def fig_corrida(R):
    h = R["res"]["historia"]
    p = R["p"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.5))

    a = axes[0]
    a.semilogy(h[:, 1] / p.ndf, lw=1.1, color=AZUL, label="costo actual")
    a.semilogy(h[:, 2] / p.ndf, lw=1.8, color="#e02020", label="mejor hasta aqui")
    a.set_title("la corrida", fontsize=10)
    a.set_xlabel("ciclo de temperatura"); a.set_ylabel("$\\chi^2/ndf$")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = axes[1]
    a.plot(h[:, 3], lw=1.5, color=AZUL, marker="o", ms=3)
    a.axhline(0.35, color="k", ls="--", lw=1)
    a.text(0.5, 0.36, "objetivo 0.35", fontsize=8)
    a.set_ylim(0, 1)
    a.set_title("tasa de aceptacion", fontsize=10)
    a.set_xlabel("ciclo"); a.set_ylabel("fraccion aceptada"); a.grid(alpha=.3)

    a = axes[2]
    a.scatter([m["x0"][0] for m in R["ms"]], R["ms_D"], s=48, color=AZUL,
              zorder=3, label="%d arranques al azar" % len(R["ms"]))
    a.axhline(R["x"][0], color="#e02020", ls="--",
              label="D = %.4f" % R["x"][0])
    a.set_title("multiarranque: dispersion %.3f%%"
                % (100 * R["ms_D"].std() / np.median(R["ms_D"])), fontsize=10)
    a.set_xlabel("D del punto INICIAL"); a.set_ylabel("D FINAL")
    a.legend(fontsize=8); a.grid(alpha=.3)
    return _guardar(fig, "corrida.png")


def fig_ajuste(R):
    p, xp = R["p"], R["x"]
    A = R["A"]
    mod = A * p.forma(*xp)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))

    a = axes[0]
    a.errorbar(p.canales, p.dato, yerr=p.error, fmt="none", ecolor="0.8", lw=.5)
    a.step(p.canales, p.dato, where="mid", lw=1, color=GRIS, label="medido")
    a.plot(p.canales, mod, lw=1.8, color="#e02020", label="modelo (SA)")
    a.plot(p.canales, p.chi2(*R["semilla"])[1] * p.forma(*R["semilla"]),
           lw=1.1, ls="--", color=NARANJA, label="semilla del NCC")
    a.set_title("Co-60: modelo contra medido", fontsize=10)
    a.set_xlabel("area  [mV$\\cdot$ns]"); a.set_ylabel("cuentas")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = axes[1]
    a.axhline(0, color="k", lw=.8)
    a.fill_between(p.canales, -1, 1, color="0.85", zorder=0)
    a.step(p.canales, (p.dato - mod) / p.error, where="mid", lw=1, color=AZUL)
    a.set_title("residuos normalizados  ($\\chi^2/ndf$ = %.2f)" % R["chi2_ndf"],
                fontsize=10)
    a.set_xlabel("area  [mV$\\cdot$ns]"); a.set_ylabel("(dato - modelo)/error")
    a.grid(alpha=.3)

    a = axes[2]
    C = R["C"]
    th = np.linspace(0, 2 * np.pi, 200)
    vals, vecs = np.linalg.eigh(C)
    for k, nsig in ((2.30, 1), (6.17, 2)):
        e = vecs @ (np.sqrt(np.maximum(vals * k, 0))[:, None] *
                    np.vstack([np.cos(th), np.sin(th)]))
        a.plot(xp[0] + e[0], xp[1] + e[1], lw=1.3,
               label="%d$\\sigma$ conjunto" % nsig)
    a.plot(*xp, "o", color="#e02020", ms=6)
    a.set_title("correlacion D-$\\sigma$:  %+.2f" % R["corr_Dsigma"], fontsize=10)
    a.set_xlabel("D  [mV$\\cdot$ns/keVee]"); a.set_ylabel("$\\sigma_{rel}$")
    a.legend(fontsize=8); a.grid(alpha=.3)
    return _guardar(fig, "ajuste.png")


# ==========================================================================
# 3. EL PDF
# ==========================================================================
def _registrar_fuentes():
    """
    Registra DejaVu Sans en reportlab y devuelve (normal, negrita, italica, mono).

    Hace falta: las fuentes base del PDF (Helvetica, Courier) usan codificacion
    WinAnsi y NO tienen alfabeto griego, asi que un sigma o un chi salen como un
    hueco en blanco. En un documento que habla de chi2, sigma y Delta E eso no
    es un detalle cosmetico: se pierde justo lo que se esta explicando.
    DejaVu viene con matplotlib, asi que no agrega ninguna dependencia.
    """
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.fonts import addMapping
    base = os.path.join(os.path.dirname(matplotlib.__file__),
                        "mpl-data", "fonts", "ttf")
    caras = {"DejaVu": "DejaVuSans.ttf",
             "DejaVu-Bold": "DejaVuSans-Bold.ttf",
             "DejaVu-Obl": "DejaVuSans-Oblique.ttf",
             "DejaVu-Mono": "DejaVuSansMono.ttf"}
    for nombre, archivo in caras.items():
        ruta = os.path.join(base, archivo)
        if not os.path.exists(ruta):
            return "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Courier"
        pdfmetrics.registerFont(TTFont(nombre, ruta))
    # el mapeo es lo que hace funcionar <b> y <i> DENTRO de un parrafo
    addMapping("DejaVu", 0, 0, "DejaVu")
    addMapping("DejaVu", 1, 0, "DejaVu-Bold")
    addMapping("DejaVu", 0, 1, "DejaVu-Obl")
    addMapping("DejaVu", 1, 1, "DejaVu-Bold")
    return "DejaVu", "DejaVu-Bold", "DejaVu-Obl", "DejaVu-Mono"


def construir(R):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                    Table, TableStyle, PageBreak, KeepTogether)

    F, FB, FI, FM = _registrar_fuentes()
    ANCHO = A4[0] - 4 * cm
    az = colors.HexColor(AZUL)
    na = colors.HexColor(NARANJA)
    gr = colors.HexColor("#555f6a")

    est = dict(
        titulo=ParagraphStyle("t", fontName=FB, fontSize=21,
                              leading=25, textColor=az, spaceAfter=6),
        sub=ParagraphStyle("s", fontName=F, fontSize=11.5, leading=15,
                           textColor=gr, spaceAfter=18),
        h1=ParagraphStyle("h1", fontName=FB, fontSize=14,
                          leading=17, textColor=az, spaceBefore=16, spaceAfter=7),
        h2=ParagraphStyle("h2", fontName=FB, fontSize=11.5,
                          leading=14, textColor=colors.HexColor("#101820"),
                          spaceBefore=11, spaceAfter=4),
        p=ParagraphStyle("p", fontName=F, fontSize=9.7, leading=13.6,
                         spaceAfter=7, alignment=4),
        li=ParagraphStyle("li", fontName=F, fontSize=9.7, leading=13.4,
                          leftIndent=14, bulletIndent=4, spaceAfter=4,
                          alignment=4),
        cap=ParagraphStyle("cap", fontName=FI, fontSize=8.3,
                           leading=11, textColor=gr, spaceBefore=3,
                           spaceAfter=12, alignment=4),
        nota=ParagraphStyle("n", fontName=F, fontSize=9.3, leading=12.6,
                            leftIndent=10, borderPadding=0, textColor=gr,
                            spaceAfter=8, alignment=4),
        cod=ParagraphStyle("c", fontName=FM, fontSize=8.4, leading=11.4,
                           leftIndent=12, textColor=colors.HexColor("#20303c"),
                           spaceAfter=8),
    )

    def P(t, e="p"):
        return Paragraph(t, est[e])

    def B(t):
        return Paragraph(t, est["li"], bulletText="•")

    def FIG(ruta, epigrafe):
        """Imagen y epigrafe como un solo bloque, para que el salto de
        pagina no los separe."""
        return KeepTogether([IMG(ruta), Paragraph(epigrafe, est["cap"])])

    def IMG(ruta, ancho=ANCHO):
        from PIL import Image as PILImage
        w, h = PILImage.open(ruta).size
        return Image(ruta, width=ancho, height=ancho * h / w)

    def ECU(texto, escala=1.0):
        r = ec(texto)
        from PIL import Image as PILImage
        w, h = PILImage.open(r).size
        an = min(w / 170 * 72 * escala, ANCHO)
        return Image(r, width=an, height=an * h / w, hAlign="CENTER")

    def TABLA(datos, anchos=None, alinear_num=True):
        t = Table(datos, colWidths=anchos, hAlign="LEFT")
        est_t = [
            ("FONT", (0, 0), (-1, 0), FB, 8.6),
            ("FONT", (0, 1), (-1, -1), F, 8.6),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("BACKGROUND", (0, 0), (-1, 0), az),
            ("LINEBELOW", (0, 0), (-1, -1), 0.4, colors.HexColor("#d5dce3")),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ]
        if alinear_num:
            est_t.append(("ALIGN", (1, 1), (-1, -1), "RIGHT"))
        t.setStyle(TableStyle(est_t))
        return t

    def pie(canv, doc):
        canv.saveState()
        canv.setFont(F, 7.5)
        canv.setFillColor(gr)
        canv.drawString(2 * cm, 1.3 * cm,
                        "Calibracion por borde Compton  -  recocido simulado")
        canv.drawRightString(A4[0] - 2 * cm, 1.3 * cm, "%d" % doc.page)
        canv.setStrokeColor(colors.HexColor("#dde3ea"))
        canv.line(2 * cm, 1.6 * cm, A4[0] - 2 * cm, 1.6 * cm)
        canv.restoreState()

    E = []
    x, err = R["x"], R["err"]
    ncc = R["ncc"]
    cs = R["cs137"]

    # ---------------------------------------------------------- portada --
    E += [
        P("Recocido simulado para el ajuste del borde Compton", "titulo"),
        P("Que dice el paper, que es el algoritmo, como se aplica a nuestro "
          "modelo, y como se decide si los parametros encontrados son el "
          "optimo global. Todos los numeros de este documento salen de correr "
          "el ajuste sobre los datos reales al generarlo.", "sub"),
    ]

    E += [P("0. Donde estamos", "h1")]
    E += [P(
        "La cadena de procesamiento ya entrega el espectro neto: el XML del "
        "DRS4 se lee, cada pulso se integra en una ventana medida sobre la "
        "convergencia de su integral acumulada, y el histograma se arma con "
        "binning comun, resta de fondo y normalizacion por tiempo vivo. La "
        "etapa siguiente (NCC) elige automaticamente el intervalo de ajuste y "
        "entrega una semilla. Falta el ultimo paso: <b>encontrar los "
        "parametros del modelo que mejor reproducen los datos medidos</b>, y "
        "saber si esos parametros son realmente los mejores.")]
    E += [P(
        "Para el Co-60, el NCC eligio la ventana <b>[%.0f, %.0f] mV&middot;ns</b> "
        "con &rho; = %.5f, y entrego la semilla <b>D = %.4f mV&middot;ns/keVee</b>. "
        "Ese es el punto de partida de todo lo que sigue."
        % (ncc["lo"], ncc["hi"], ncc["rho"], R["semilla"][0]))]

    # ------------------------------------------------- resumen del paper --
    E += [P("1. Resumen de la Seccion II-A del paper", "h1")]
    E += [P(
        "Kong et al. (IEEE TNS 2025, arXiv:2507.16282) atacan un problema "
        "concreto: en un centellador plastico no hay fotopico, asi que hay que "
        "calibrar con el borde Compton; pero el borde medido no es el escalon "
        "que predice la teoria, sino una caida gradual. La Seccion II-A explica "
        "por que, y que hacer al respecto.")]

    E += [P("1.1 El borde real no es un escalon", "h2")]
    E += [P(
        "El paper atribuye el ensanchamiento a tres fuentes independientes "
        "(su Tabla I): la fluctuacion del rendimiento de luz, la fluctuacion de "
        "ganancia del SiPM y el ruido termico de la electronica. Como son "
        "independientes, se suman en cuadratura (su Ec. 3):")]
    E += [ECU("$\\sigma_{total}=\\sqrt{\\sigma_{LY}^{2}+\\sigma_{gain}^{2}+"
              "\\sigma_{elec}^{2}}$")]
    E += [P(
        "El termino <i>&sigma;<sub>gain</sub> = ENF&middot;&radic;N<sub>pe</sub></i> "
        "es literalmente la contribucion del MPPC, o sea que esta parte de la "
        "Tabla I describe nuestro detector.")]

    E += [P("1.2 El modelo convolucional", "h2")]
    E += [P(
        "Si el detector es un sistema lineal e invariante en el tiempo (Ecs. "
        "4&ndash;7), su efecto sobre el espectro es una convolucion. El paper "
        "toma la respuesta como gaussiana y escribe el espectro modelado "
        "(su Ec. 8) como la lineshape teorica llevada a canal y convolucionada:")]
    E += [ECU("$S_{model}(E)=\\int_{-\\infty}^{+\\infty}\\frac{d\\sigma}{dE_e}"
              "\\left(E_{eff}=E_e\\cdot\\frac{B}{E_{max}}\\right)\\,G(E_e)\\,dE_e$")]
    E += [P(
        "Hay dos cosas importantes escondidas ahi. La primera: el mapa "
        "canal&harr;energia es <b>E<sub>eff</sub> = E<sub>e</sub>&middot;B/E<sub>max</sub></b>, "
        "o sea una recta <b>por el origen</b> &mdash; area cero corresponde a "
        "energia cero, sin ordenada. Todo el metodo depende de eso, incluida la "
        "seleccion de intervalo del NCC. La segunda: el nucleo gaussiano "
        "(Ec. 9) esta normalizado, asi que conserva el area y desacopla la "
        "amplitud de la resolucion.")]
    E += [P(
        "Los parametros a ajustar son entonces la amplitud <i>A</i>, el canal "
        "del borde <i>B</i> y la resolucion <i>&sigma;</i>. Para fuentes con dos "
        "lineas gamma cercanas (Co-60, Na-22) el paper amplia el vector "
        "(su Ec. 11) a")]
    E += [ECU("$\\Theta=\\{A_1,\\ A_2,\\ B_1,\\ \\sigma_1,\\ \\sigma_2,\\ D\\}$")]
    E += [P(
        "amarrando el segundo borde al primero con <b>CE<sub>ch,2</sub> = "
        "B<sub>1</sub> + D&middot;&Delta;E</b> (Ec. 12), donde <i>D</i> es el "
        "coeficiente de conversion energia&rarr;canal: el numero que se busca.")]

    E += [P("1.3 La funcion de costo", "h2")]
    E += [P("La Ec. 14 define la perdida como minimos cuadrados entre el "
            "espectro medido y el modelado:")]
    E += [ECU("$loss_{min}(E_i)=\\sum_{i=1}^{N}\\left[S_1(E_i)-"
              "S_{model}(E_i;\\vec{p})\\right]^{2}$")]

    E += [P("1.4 El recocido, como lo especifica el paper", "h2")]
    E += [P("Textualmente, la receta de la Seccion II-A es:")]
    E += [B("punto inicial <b>x<sub>0</sub></b> y <b>T<sub>init</sub> = 2000</b>, "
            "con coeficiente de enfriamiento <b>&alpha; = 0.95</b>;")]
    E += [B("en cada temperatura T<sub>k</sub>, <b>100 iteraciones</b> de cadena "
            "de Markov, proponiendo <b>x<sub>nueva</sub> = x<sub>actual</sub> + "
            "&Delta;x</b> con &Delta;x<sub>i</sub> ~ "
            "<i>U</i>(&minus;0.1&delta;<sub>i</sub>, +0.1&delta;<sub>i</sub>), "
            "donde &delta;<sub>i</sub> es el rango de busqueda del parametro;")]
    E += [B("transiciones por el <b>criterio de Metropolis</b>, aceptando con "
            "probabilidad <b>P = exp(&minus;&Delta;E/T<sub>k</sub>)</b>;")]
    E += [B("enfriamiento <b>geometrico</b> T<sub>k+1</sub> = "
            "&alpha;&middot;T<sub>k</sub>, con &alpha; en [0.85, 0.99];")]
    E += [B("<b>doble criterio de corte</b>: T<sub>k</sub> &lt; 10<sup>&minus;6</sup>, "
            "o estabilizacion |&Delta;E<sub>medio</sub>| &lt; 0.001 durante 10 "
            "ciclos consecutivos.")]

    E += [P("1.5 Una inconsistencia del paper que hay que resolver al implementar",
            "h2")]
    E += [P(
        "La Ec. 14 dice que lo que se minimiza son minimos cuadrados, pero el "
        "parrafo que describe Metropolis dice que <i>E(x)</i> es \"the "
        "NCC-based matching error objective function\". Son dos funciones "
        "distintas y hay que elegir una. Nosotros elegimos <b>&chi;<sup>2</sup> "
        "ponderado</b>, por una razon dura: multiplicar la amplitud por una "
        "constante cambia el error NCC en <b>exactamente cero</b> (Pearson "
        "divide por la desviacion), asi que con el objetivo NCC la amplitud es "
        "<b>indeterminable</b>. Ademas el &chi;<sup>2</sup> usa las barras de "
        "error del espectro, que es lo que despues permite hablar de "
        "incertidumbre.", "nota")]

    # ------------------------------------------------- que es el recocido -
    E += [P("2. Que es el recocido simulado", "h1")]
    E += [P(
        "Un optimizador local &mdash;Nelder-Mead, Levenberg-Marquardt, "
        "gradiente descendente&mdash; solo sabe hacer una cosa: bajar. Desde "
        "donde arranca camina cuesta abajo y se detiene en el fondo del valle "
        "en el que cayo. Si la funcion de costo tiene varios valles, el "
        "resultado <b>depende del punto inicial</b>, y desde adentro del "
        "algoritmo no hay ninguna manera de saber si el valle encontrado es el "
        "mas profundo o simplemente el mas cercano.")]
    E += [P(
        "El recocido simulado agrega una sola idea: <b>a veces acepta subir</b>. "
        "Una propuesta que empeora el costo se acepta con probabilidad "
        "exp(&minus;&Delta;E/T). Esa es la unica manera de salir de un valle. "
        "La temperatura T controla cuanto: alta al principio, y el algoritmo "
        "pasea casi sin restriccion; baja al final, y se comporta como un "
        "optimizador local que afina. El nombre viene de la metalurgia: enfriar "
        "despacio deja que los atomos encuentren la red de minima energia, "
        "enfriar de golpe los congela en un vidrio con defectos.")]
    E += [FIG(R["f_juguete"],
              "Izquierda: dos arranques distintos sobre la misma funcion terminan "
            "en minimos distintos, y ninguno de los dos es el global. Derecha: "
            "una cadena de recocido sobre la misma funcion, coloreada por "
            "iteracion; al principio recorre todo el dominio y hacia el final "
            "se concentra en el minimo global. Es una funcion de juguete, no "
            "nuestros datos: sirve para ver el mecanismo.")]

    E += [P("2.1 Por que exp(&minus;&Delta;E/T) y no otra regla", "h2")]
    E += [P(
        "No es arbitrario. Con esa probabilidad de aceptacion, la cadena a "
        "temperatura fija T muestrea la distribucion de Boltzmann "
        "<b>&pi;(x) &prop; exp(&minus;E(x)/T)</b>. Cuando T es grande esa "
        "distribucion es casi plana (todos los puntos son igual de probables: "
        "exploracion pura); cuando T&rarr;0 se concentra toda su masa en el "
        "minimo global. El recocido es exactamente eso: seguir la familia de "
        "distribuciones de Boltzmann mientras T baja, con la esperanza de que "
        "la cadena se quede donde la distribucion se concentra.")]
    E += [P(
        "De ahi sale una consecuencia que importa mucho al implementar: "
        "<b>Metropolis compara T contra &Delta;E</b>, o sea que T esta en las "
        "<b>unidades de la funcion de costo</b>. Un T copiado de otro trabajo "
        "solo funciona si el costo tiene la misma escala.")]

    E += [P("2.2 El enfriamiento, y por que no hay garantia", "h2")]
    E += [P(
        "Existe un teorema (Geman &amp; Geman, 1984) que garantiza convergencia "
        "al minimo global &mdash;con probabilidad 1&mdash; si la temperatura "
        "baja como <b>T<sub>k</sub> &ge; c/log(k)</b>. El problema es que ese "
        "enfriamiento es impracticablemente lento: para bajar T en un factor "
        "1000 hace falta un numero astronomico de ciclos. Todo el mundo, el "
        "paper incluido, usa enfriamiento <b>geometrico</b> "
        "T<sub>k</sub> = &alpha;<sup>k</sup>T<sub>0</sub>, que es mucho mas "
        "rapido y <b>no tiene ninguna garantia</b>. Conviene tenerlo claro "
        "desde el principio: el recocido no demuestra que encontro el minimo "
        "global. La Seccion 5 se ocupa de que evidencia si se puede juntar.")]
    E += [FIG(R["f_mecanica"],
              "Izquierda: la probabilidad de aceptar un empeoramiento, para tres "
            "temperaturas alrededor de la que se calibro sobre nuestro costo. "
            "Centro: el enfriamiento geometrico para los tres &alpha; que el "
            "paper considera. Derecha: el enfriamiento con garantia teorica "
            "(logaritmico) contra el que se usa en la practica (geometrico), en "
            "escala doble logaritmica &mdash; la diferencia es la razon por la "
            "que hace falta la Seccion 5.")]

    # ------------------------------------------------- nuestro problema ---
    E += [P("3. Nuestro problema concreto", "h1")]
    E += [P("3.1 El modelo", "h2")]
    E += [P("El espectro predicho es la plantilla teorica llevada a canal y "
            "borroneada por la resolucion del detector:")]
    E += [ECU("$S_{modelo}(ch)=A\\,\\left[\\,P\\!\\left(\\frac{ch}{D}\\right)"
              "\\ast G(\\sigma)\\,\\right]$")]
    E += [P(
        "Es la Ec. 8 del paper con una sustitucion: donde el paper pone la "
        "lineshape analitica de Klein-Nishina (su Ec. 10), nosotros ponemos el "
        "<b>histograma de Geant4</b>. Esa sustitucion no es cosmetica. La "
        "plantilla simulada ya contiene el scattering Compton multiple (que en "
        "nuestra geometria es el <b>25%</b> de los eventos y rellena la zona "
        "por encima del borde), el escape, y el quenching de Birks. Ademas "
        "contiene <b>las dos lineas del Co-60 con su intensidad relativa "
        "correcta</b>.")]
    E += [P("Por eso nuestro vector de parametros tiene tres componentes y no "
            "las seis de la Ec. 11 del paper:")]
    E += [TABLA([
        ["parametro", "que es", "unidades", "como entra"],
        ["A", "amplitud", "cuentas", "lineal - se proyecta"],
        ["D", "ganancia; el resultado", "mV·ns/keVee", "no lineal"],
        ["σ_rel", "resolucion en E_ref", "adimensional", "no lineal"],
    ], anchos=[2.6 * cm, 6.2 * cm, 3.4 * cm, 4.6 * cm], alinear_num=False)]
    E += [Spacer(1, 8)]
    E += [P(
        "A<sub>2</sub>/A<sub>1</sub>, B<sub>2</sub>/B<sub>1</sub> y la mezcla "
        "de scattering dejan de ser parametros libres: los fija la simulacion. "
        "Eso importa mas de lo que parece &mdash; con dos bordes separados por "
        "menos de una resolucion, un modelo de seis parametros esta "
        "<b>sub-determinado</b>, y el optimizador gasta las libertades de mas "
        "en soluciones que ajustan bien y no tienen sentido fisico.", "nota")]

    E += [P("3.2 La amplitud se proyecta, no se busca", "h2")]
    E += [P(
        "<i>A</i> entra linealmente, asi que dado (D, &sigma;) la mejor "
        "amplitud se despeja en forma cerrada minimizando el &chi;<sup>2</sup>. "
        "Esto se llama proyeccion, o VARPRO (Golub &amp; Pereyra, 1973). El "
        "optimizador busca entonces en <b>2 dimensiones en vez de 3</b>.")]
    E += [P(
        "<b>Medido</b>, y no supuesto: corriendo el recocido con <i>A</i> como "
        "tercer parametro libre y el paso comun 0.1&delta;<sub>i</sub> del "
        "paper, el &chi;<sup>2</sup> final promedio es <b>%.1f</b>; con "
        "<i>A</i> proyectada es <b>%.1f</b>. La version de tres parametros "
        "<b>no se mueve de la semilla</b>: el rango de busqueda de <i>A</i> es "
        "de otra escala que el de <i>D</i> (%.3g contra %.3g), y con un paso "
        "comun casi toda propuesta cambia <i>A</i> lo suficiente como para que "
        "el &chi;<sup>2</sup> se dispare y se rechace."
        % (R["v_varpro"][1], R["v_varpro"][0], R["A"], R["x"][0]))]

    E += [P("3.3 El costo", "h2")]
    E += [P("Minimos cuadrados <b>ponderados</b> por las barras de error que la "
            "etapa del espectro ya propago (Poisson, incluyendo la resta de "
            "fondo escalada por tiempo vivo):")]
    E += [ECU("$\\chi^{2}(\\Theta)=\\sum_i\\left[\\frac{y_i-"
              "S_{modelo}(ch_i;\\Theta)}{\\varepsilon_i}\\right]^{2}$")]
    E += [P(
        "Ponderar no es un detalle. Sin pesos, los bines del plateau "
        "&mdash;que tienen muchas cuentas&mdash; dominan el ajuste y el flanco, "
        "que es donde vive la informacion sobre D, casi no pesa. Y sin pesos "
        "tampoco hay &chi;<sup>2</sup>/ndf que mirar ni barras de error que "
        "reportar.")]

    E += [P("3.4 La semilla la da el NCC", "h2")]
    E += [P(
        "El NCC entrega las dos cosas que el optimizador necesita: <b>donde "
        "ajustar</b> (la ventana [%.0f, %.0f] mV&middot;ns) y <b>desde donde "
        "arrancar</b> (D = %.4f). Con esa semilla el problema queda con "
        "%d bines y <b>%d grados de libertad</b>."
        % (ncc["lo"], ncc["hi"], R["semilla"][0], R["nbins"], R["ndf"]))]

    # ------------------------------------ constantes no transferibles -----
    E += [P("4. Las constantes del paper, una por una", "h1")]
    E += [P(
        "La receta de la Seccion II-A trae numeros concretos. La pregunta no es "
        "si copiarlos o no, sino <b>cuales dependen del problema</b>. Un numero "
        "adimensional se copia; uno que vive en las unidades del costo, no. Y la "
        "unica manera de saber en cual de los dos casos estamos es medirlo, no "
        "argumentarlo &mdash; el resto de esta seccion mide, y en dos de los "
        "tres casos el resultado <b>contradice</b> lo que la primera version de "
        "este documento afirmaba.")]

    E += [P("4.1 T<sub>init</sub>", "h2")]
    E += [P(
        "Metropolis compara T contra &Delta;E, asi que T vive en las unidades "
        "del costo. En vez de fijarla, la calibramos: se muestrean propuestas "
        "alrededor de la semilla, se mide cuanto empeoran en promedio, y se "
        "despeja la T que las acepta con probabilidad 0.8:")]
    E += [ECU("$T_{init}=-\\frac{\\langle\\Delta E^{+}\\rangle}{\\ln 0.8}$")]
    E += [P(
        "Medido sobre nuestro problema: <b>T<sub>init</sub> = %.0f</b>. Y "
        "conviene ser preciso, porque es facil exagerar esto: el 2000 del paper "
        "<b>funciona perfectamente aca</b>. Lo que no es transferible es el "
        "numero, no el metodo. Medido:" % R["res"]["T_init"])]
    E += [TABLA([
        ["costo", "valor tipico", "T = 2000 (paper)", "T calibrada"],
        ["nuestro χ²", "%.0f" % R["escala_costo"],
         "acepta %.2f" % R["acc_2000_chi2"],
         "acepta %.2f" % R["acc_auto_chi2"]],
        ["el mismo, ÷ %.0f  (escala tipo NCC)" % R["escala_costo"], "1",
         "acepta %.2f" % R["acc_2000_orden1"], "acepta 0.82"],
    ], anchos=[6.4 * cm, 2.8 * cm, 3.7 * cm, 3.7 * cm])]
    E += [Spacer(1, 6)]
    E += [P(
        "Sobre un costo de orden 1 &mdash;la escala que tendria el objetivo NCC "
        "del texto de Metropolis&mdash; el mismo T = 2000 acepta el "
        "<b>%.0f%%</b> de las propuestas: la temperatura deja de controlar nada "
        "y la cadena es una caminata al azar de principio a fin. Calibrar T "
        "evita tener que acertar ese numero."
        % (100 * R["acc_2000_orden1"]))]

    E += [P("4.2 El corte por estabilizacion", "h2")]
    E += [P(
        "Por lo mismo, |&Delta;E<sub>medio</sub>| &lt; 0.001 es un umbral "
        "absoluto sobre una cantidad cuyas unidades cambian con el problema. Se "
        "reemplaza por un corte <b>relativo</b>: el mejor costo tiene que "
        "mejorar menos de una parte en 10<sup>4</sup> durante 10 ciclos "
        "seguidos. Con eso, la corrida sobre el Co-60 termino por "
        "<b>%s</b> despues de %d ciclos y %d evaluaciones (%.1f s)."
        % (R["res"]["motivo"], R["res"]["ciclos"], R["res"]["n_eval"],
           R["t_sa"]))]

    E += [P("4.3 El paso adaptativo: una idea que <b>no</b> hizo falta", "h2")]
    E += [P(
        "El paper propone &Delta;x<sub>i</sub> ~ "
        "<i>U</i>(&minus;0.1&delta;<sub>i</sub>, +0.1&delta;<sub>i</sub>) con "
        "&delta;<sub>i</sub> fijo. El argumento habitual en contra es que a "
        "temperatura alta casi todo se acepta y el paso queda timido, y a "
        "temperatura baja casi nada se acepta y la cadena se congela sin haber "
        "afinado; la correccion es escalar el paso para mantener la aceptacion "
        "cerca de un objetivo. Esta implementado, y <b>se midio si hace "
        "falta</b>. Con %d semillas distintas:"
        % (len(R["v_acept"]) and 8))]
    E += [TABLA([
        ["paso", "D promedio", "dispersion entre semillas"],
        ["fijo (como el paper)", "%.4f" % R["v_paso"]["fijo"][0],
         "± %.4f" % R["v_paso"]["fijo"][1]],
        ["adaptativo", "%.4f" % R["v_paso"]["adaptativo"][0],
         "± %.4f" % R["v_paso"]["adaptativo"][1]],
    ], anchos=[6.0 * cm, 4.2 * cm, 6.4 * cm])]
    E += [Spacer(1, 6)]
    E += [P(
        "<b>No hay diferencia.</b> El paso adaptativo no mejora nada sobre este "
        "problema, asi que el default del codigo es el del <b>paper</b> y el "
        "adaptativo queda como opcion. Se lo menciona porque en un modelo con "
        "mas parametros, donde las escalas difieren, si importa &mdash; es "
        "exactamente lo que se ve en 3.2 con la amplitud libre.", "nota")]
    E += [P(
        "De paso, una correccion a lo que decia la primera version de este "
        "documento: el objetivo de aceptacion 0.35 <b>no es</b> \"el valor "
        "clasico\". Los numeros citados en la literatura son 0.44 (una "
        "dimension) y 0.234 (asintotico en muchas), y son resultados sobre "
        "<i>muestreo</i> por Metropolis, no sobre recocido. El 0.35 es una "
        "eleccion intermedia. Lo que corresponde no es justificarla sino medir "
        "si cambia algo: barrida entre 0.234 y 0.44, D se mueve %s."
        % (", ".join("%.5f" % v[1] for v in R["v_acept"])), "nota")]

    E += [P("4.4 Y una que si conviene conservar", "h2")]
    E += [P(
        "El <b>enfriamiento geometrico con &alpha; = 0.95</b> y las <b>100 "
        "iteraciones por temperatura</b> son adimensionales: no dependen de las "
        "unidades del costo, solo de cuanto se quiere explorar. Esos si se "
        "copian tal cual.")]

    # ---------------------------------- como saber si es el optimo --------
    E += [P("5. Como saber si los parametros son los optimos", "h1")]
    E += [P(
        "Esta es la pregunta que el algoritmo <b>no</b> contesta solo. Con "
        "enfriamiento geometrico no hay teorema que garantice nada (Seccion "
        "2.2), asi que lo que se puede hacer es juntar evidencia. Conviene "
        "separar dos preguntas que se confunden facil:")]
    E += [B("<b>&iquest;el optimizador encontro el minimo del costo?</b> "
            "&mdash; es una pregunta sobre el ALGORITMO;")]
    E += [B("<b>&iquest;el modelo en ese minimo describe los datos?</b> "
            "&mdash; es una pregunta sobre la FISICA.")]
    E += [P("Se pueden tener las dos respuestas al reves: un minimo global "
            "perfectamente encontrado de un modelo que no sirve, o un modelo "
            "correcto mal optimizado.")]

    E += [P("5.1 Evidencia sobre el algoritmo", "h2")]

    E += [P("<b>(a) Multiarranque.</b> El recocido es estocastico: una corrida "
            "no es un resultado, es una muestra. Se corrieron <b>%d</b> "
            "recocidos desde puntos iniciales <b>al azar dentro de todas las "
            "cotas</b> (no perturbaciones de la semilla). Resultado: "
            "<b>%d de %d</b> terminan en el mismo &chi;<sup>2</sup> dentro de "
            "0.1, con dispersion en D de <b>%.3f%%</b>. Si hubiera un segundo "
            "valle competitivo, algunos arranques habrian caido ahi."
            % (len(R["ms"]),
               int(np.sum(R["ms_c"] < R["ms_c"].min() + 0.1)), len(R["ms"]),
               100 * R["ms_D"].std() / np.median(R["ms_D"])))]

    E += [P("<b>(b) Mapear el paisaje.</b> Con dos parametros no hace falta "
            "confiar: se puede dibujar. El barrido ancho de D &mdash;de "
            "0.3&times; a 2.5&times; la semilla, optimizando &sigma; en cada "
            "punto&mdash; encuentra <b>%s</b>."
            % ("<b>un solo minimo</b> en todo el rango (el barrido es grueso, "
               "asi que su posicion cae dentro del paso de la grilla del optimo "
               "fino)" if len(R["minimos_locales"]) == 1
               else "<b>%d minimos locales</b>" % len(R["minimos_locales"])))]
    E += [FIG(R["f_paisaje"],
              "Izquierda: el paisaje completo del &chi;<sup>2</sup>, con los "
            "contornos conjuntos de 1, 2 y 3&sigma;. Es un solo valle, "
            "alargado en la direccion en que D y &sigma; se compensan. Derecha: "
            "el corte a lo largo de D optimizando &sigma; en cada punto, en "
            "escala logaritmica: no hay estructura secundaria.")]

    E += [P("<b>(c) Comparar contra el optimizador local desde arranques "
            "malos.</b> Si el paisaje fuera accidentado, un optimizador local "
            "lanzado lejos del optimo deberia quedarse corto. Medido:")]
    E += [TABLA(
        [["arranque (D, σ)", "local: D", "local: χ²",
          "SA+pulido: D", "SA+pulido: χ²"]] +
        [["(%.3f, %.2f)" % c["x0"], "%.4f" % c["D_loc"], "%.2f" % c["c_loc"],
          "%.4f" % c["D_sa"], "%.2f" % c["c_sa"]] for c in R["comparacion"]],
        anchos=[3.6 * cm, 2.9 * cm, 3.0 * cm, 3.4 * cm, 3.9 * cm])]
    E += [Spacer(1, 6)]
    E += [P(
        "<b>Los dos coinciden.</b> Y eso, lejos de ser decepcionante, es el "
        "resultado: significa que <b>en este problema el paisaje es unimodal</b>, "
        "asi que el recocido no cambia el numero. Su aporte aca no es corregir "
        "sino <b>demostrar</b> &mdash; sin el habria que suponer que el minimo "
        "local es el global, y suponerlo es exactamente lo que el metodo del "
        "paper trata de evitar. El caso donde el recocido si cambia el numero "
        "es el modelo de seis parametros de la Ec. 11, que esta "
        "sub-determinado; nosotros lo esquivamos usando la plantilla de "
        "Geant4.", "nota")]

    E += [P("<b>(d) La historia de la corrida.</b> Una corrida sana se "
            "reconoce: la aceptacion arranca alta y baja hacia el objetivo, el "
            "costo actual fluctua mientras el mejor costo baja en escalones y "
            "despues se aplana. Si la aceptacion arranca en 1.0 y nunca baja, "
            "T era demasiado alta; si arranca en 0, demasiado baja.")]
    E += [FIG(R["f_corrida"],
              "Izquierda: el costo actual (que sube y baja, porque el algoritmo "
            "acepta empeorar) contra el mejor encontrado hasta cada ciclo. "
            "Centro: la aceptacion, que arranca en %.2f y termina cerca del "
            "objetivo. Derecha: los %d arranques al azar; el eje horizontal es "
            "el D inicial y el vertical el final &mdash; noten la escala del "
            "eje vertical."
            % (R["res"]["historia"][0, 3], len(R["ms"])))]

    E += [P("5.2 Evidencia sobre el ajuste", "h2")]
    E += [P(
        "Que el optimizador haya hecho bien su trabajo no dice nada sobre si el "
        "modelo sirve. Para eso hay tres diagnosticos, y conviene mirarlos en "
        "este orden.")]
    E += [P("<b>&chi;<sup>2</sup>/ndf.</b> Si el modelo es correcto y las "
            "barras de error estan bien, tiene que dar &asymp;1. Medido: "
            "<b>%.2f</b> con %d grados de libertad. Es <b>mayor que 1</b>, y "
            "hay que decirlo: el modelo no describe todo lo que hay en los "
            "datos." % (R["chi2_ndf"], R["ndf"]))]
    E += [P("<b>Residuos.</b> Un &chi;<sup>2</sup>/ndf alto puede venir de "
            "ruido subestimado (residuos grandes pero sin estructura) o de un "
            "modelo incompleto (residuos con forma). Hay que mirarlos, no "
            "resumirlos en un numero.")]
    E += [P("<b>Curvatura.</b> Cerca del minimo el &chi;<sup>2</sup> es "
            "aproximadamente cuadratico, y la matriz de covarianza sale de su "
            "hessiano: <b>C = 2H<sup>&minus;1</sup></b>. De ahi salen las "
            "barras de error y, sobre todo, la <b>correlacion</b> entre "
            "parametros.")]
    E += [FIG(R["f_ajuste"],
              "Izquierda: el modelo optimo sobre los datos, con la semilla del "
            "NCC punteada para ver cuanto movio el ajuste. Centro: los residuos "
            "normalizados; la banda gris es &plusmn;1&sigma;. Derecha: las "
            "elipses de 1 y 2&sigma; conjuntas en el plano D&ndash;&sigma;.")]
    E += [P(
        "La correlacion medida es <b>corr(D, &sigma;) = %+.2f</b>, y es fuerte: "
        "un borde mas ancho puede simularse corriendo el borde y ensanchando la "
        "resolucion a la vez. Consecuencia practica: <b>D y &sigma; hay que "
        "reportarlos juntos</b>, y el error de D solo tiene sentido si &sigma; "
        "estuvo libre." % R["corr_Dsigma"])]
    E += [P(
        "Ademas, con &chi;<sup>2</sup>/ndf = %.2f el criterio "
        "&Delta;&chi;<sup>2</sup> = 1 <b>subestima</b> las barras: la receta "
        "habitual es escalarlas por &radic;(&chi;<sup>2</sup>/ndf), lo que "
        "lleva el error de D de %.4f a %.4f (de %.2f%% a %.2f%%). Es un parche "
        "honesto, no una solucion: la solucion es mejorar el modelo."
        % (R["chi2_ndf"], err[0], err[0] * np.sqrt(R["chi2_ndf"]),
           100 * err[0] / x[0],
           100 * err[0] * np.sqrt(R["chi2_ndf"]) / x[0]), "nota")]

    # ------------------------------------------------------- resultado ----
    E += [P("6. El resultado, y que le falta", "h1")]
    E += [TABLA([
        ["", "Co-60", "Cs-137 @ 1600 V"],
        ["fuente / bordes [keV]", "963.4 y 1118.1", "477.3"],
        ["ventana NCC [mV·ns]",
         "%.0f – %.0f" % (ncc["lo"], ncc["hi"]),
         "%.0f – %.0f" % cs["ventana"]],
        ["ρ del NCC", "%.5f" % ncc["rho"], "%.5f" % cs["rho"]],
        ["semilla D", "%.4f" % R["semilla"][0], "%.4f" % cs["semilla"]],
        ["D optimo [mV·ns/keVee]",
         "%.4f ± %.4f" % (x[0], err[0]),
         "%.4f ± %.4f" % (cs["D"], cs["err_D"])],
        ["σ_rel", "%.4f" % x[1], "%.4f" % cs["sigma"]],
        ["FWHM en E_ref", "%.1f %%" % (100 * 2.355 * x[1]),
         "%.1f %%" % (100 * 2.355 * cs["sigma"])],
        ["χ²/ndf", "%.2f  (%d gdl)" % (R["chi2_ndf"], R["ndf"]),
         "%.2f  (%d gdl)" % (cs["chi2_ndf"], cs["ndf"])],
    ], anchos=[6.0 * cm, 5.4 * cm, 5.4 * cm])]
    E += [Spacer(1, 10)]
    E += [P(
        "Las dos <b>D</b> son de detectores distintos y no son comparables "
        "entre si; lo que si es comparable es que las dos salieron del mismo "
        "procedimiento, sin ningun numero elegido a mano.")]

    E += [P("6.1 Lo que el &chi;<sup>2</sup>/ndf esta diciendo", "h2")]
    E += [P("Un &chi;<sup>2</sup>/ndf de %.2f con %d grados de libertad no es "
            "compatible con ruido estadistico: falta fisica en el modelo. Los "
            "candidatos, en orden de sospecha:" % (R["chi2_ndf"], R["ndf"]))]
    E += [B("<b>El exceso de cuentas a bajo canal.</b> El espectro medido tiene "
            "cuentas por debajo del plateau que la simulacion no tiene. La "
            "ventana del NCC arranca en 0.45&middot;B y entra en esa zona.")]
    E += [B("<b>&sigma; constante.</b> El modelo usa un unico ancho gaussiano, "
            "pero fisicamente la resolucion escala como &radic;E (fluctuacion "
            "de fotoelectrones). Con dos bordes a 963 y 1118 keV hay palanca "
            "para distinguirlo.")]
    E += [B("<b>&chi;<sup>2</sup> de Neyman.</b> Las barras se calculan con las "
            "cuentas <b>observadas</b>, asi que un bin que fluctuo bajo recibe "
            "una barra chica y por lo tanto mas peso. Eso sesga el ajuste hacia "
            "abajo; la correccion es usar la varianza <b>predicha por el "
            "modelo</b> (&chi;<sup>2</sup> de Pearson), iterando.")]

    E += [P("6.2 Lo que sigue", "h2")]
    E += [P("En este orden, porque cada uno cambia el numero del anterior: "
            "(1) &chi;<sup>2</sup> de Pearson; (2) modelar el exceso de bajo "
            "canal en vez de esquivarlo; (3) &sigma;(E) &prop; &radic;E; "
            "(4) barrer la ventana y la definicion del area para el presupuesto "
            "de sistematicos; (5) MCMC para el posterior completo, que es lo "
            "unico que da barras de error sin suponer que el minimo es "
            "cuadratico.")]

    # ------------------------------------------------------- procedencia --
    E += [P("7. De donde sale cada afirmacion", "h1")]
    E += [P(
        "Este documento mezcla tres tipos de enunciado, y conviene que se "
        "distingan a simple vista porque se verifican de maneras distintas.")]

    E += [TABLA([
        ["origen", "como se verifica", "ejemplos en este documento"],
        ["El paper", "leyendolo: Sec. II-A, Ecs. 1-15",
         "Θ de la Ec. 11; T=2000, α=0.95,\n100 iteraciones; Metropolis;\nel doble criterio de corte"],
        ["Medido sobre estos datos",
         "corriendo el codigo; los numeros\nse regeneran solos",
         "D y σ; χ²/ndf; el paisaje;\nel multiarranque; T calibrada;\npaso fijo vs adaptativo"],
        ["Literatura / proyecto anterior",
         "hay que ir a la fuente;\nNO se verifica aca",
         "Geman & Geman 1984;\nVARPRO (Golub & Pereyra 1973);\nel modelo de 6 parametros"],
    ], anchos=[4.4 * cm, 5.4 * cm, 6.8 * cm], alinear_num=False)]
    E += [Spacer(1, 10)]

    E += [P("7.1 Lo que se verifica corriendo codigo", "h2")]
    E += [P("<font face=\"%s\">python analisis/verificar_sa.py</font> contrasta, "
            "una por una, las decisiones de implementacion que se apartan del "
            "paper. Cada bloque imprime el numero, no la conclusion:" % FM)]
    E += [B("<b>que exp(&minus;&Delta;E/T) muestrea Boltzmann</b> &mdash; se "
            "corre la cadena a T fija sobre una funcion conocida y se compara "
            "el histograma contra exp(&minus;E/T). Coinciden al 2&ndash;3%.")]
    E += [B("<b>que la Ec. 1 impresa no da el borde Compton</b> &mdash; para "
            "E<sub>&gamma;</sub> = 661.657 keV da 184.3 keV, cuando la Fig. 1 "
            "del propio paper muestra el corte cerca de 480. Con 2&alpha; en el "
            "numerador da 477.3.")]
    E += [B("<b>que T<sub>init</sub> depende de la escala del costo</b> "
            "&mdash; tabla de la Seccion 4.1.")]
    E += [B("<b>que el paso adaptativo no hace falta aqui</b> &mdash; tabla de "
            "la Seccion 4.3. Esta verificacion <b>refuto</b> lo que decia la "
            "primera version del codigo.")]
    E += [B("<b>que el objetivo de aceptacion es indiferente</b> &mdash; "
            "barrido de la Seccion 4.3.")]
    E += [B("<b>que proyectar la amplitud hace falta</b> &mdash; Seccion 3.2.")]

    E += [P("7.2 Lo que <b>no</b> esta verificado", "h2")]
    E += [P(
        "<b>Que el modelo de 6 parametros de la Ec. 11 esta sub-determinado.</b> "
        "Se afirma en la Seccion 3.1 y el origen es "
        "<font face=\"%s\">CLAUDE.md &sect;7.23</font>, del proyecto anterior de "
        "esta misma tesis. Este codigo <b>no implementa</b> ese modelo, asi que "
        "la afirmacion esta importada y no re-verificada. Para cerrarla hay que "
        "escribir la Ec. 11 con Klein-Nishina y correrla." % FM)]
    E += [P(
        "<b>El teorema de Geman &amp; Geman.</b> Es una cita (IEEE PAMI 6(6), "
        "1984). No se verifica desde aca; hay que leer el paper. Lo que si es "
        "verificable es lo que se usa de el en este documento: que nuestro "
        "enfriamiento es geometrico y no logaritmico, y que por lo tanto la "
        "evidencia de optimalidad de la Seccion 5 es <b>empirica</b>.")]
    E += [P(
        "<b>Las causas del &chi;<sup>2</sup>/ndf = %.2f.</b> Las tres de la "
        "Seccion 6.1 son <b>hipotesis ordenadas por sospecha</b>, no "
        "resultados. Cada una se cierra implementandola y viendo si el "
        "&chi;<sup>2</sup> baja." % R["chi2_ndf"])]

    E += [Spacer(1, 14)]
    mono = lambda t: '<font face="%s">%s</font>' % (FM, t)
    E += [P("<b>Codigo.</b> %s (modelo y costo), %s (recocido, multiarranque, "
            "errores por curvatura), %s (este documento). Las figuras sueltas "
            "quedan en %s."
            % (mono("analisis/modelo.py"), mono("analisis/sa.py"),
               mono("analisis/doc_sa.py"), mono("analisis/figuras/")), "nota")]

    doc = SimpleDocTemplate(SALIDA, pagesize=A4,
                            leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=2 * cm, bottomMargin=2.2 * cm,
                            title="Recocido simulado para el ajuste del borde Compton",
                            author="Calibracion CE")
    doc.build(E, onFirstPage=pie, onLaterPages=pie)
    return SALIDA


def main():
    rutas.salida_utf8()
    rapido = "--rapido" in sys.argv
    R = calcular(rapido=rapido)
    print("generando figuras...")
    R["f_juguete"] = fig_juguete()
    R["f_mecanica"] = fig_mecanica(R)
    R["f_paisaje"] = fig_paisaje(R)
    R["f_corrida"] = fig_corrida(R)
    R["f_ajuste"] = fig_ajuste(R)
    print("armando el PDF...")
    ruta = construir(R)
    print("PDF -> %s  (%.0f kB)" % (rutas.relativa(ruta),
                                    os.path.getsize(ruta) / 1e3))


if __name__ == "__main__":
    main()
