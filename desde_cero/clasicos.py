"""
Etapa 10 (opcional) — los metodos convencionales, y la robustez frente al
intervalo (paper, Seccion III, Fig. 4 y 5).

    python desde_cero/clasicos.py

Que resuelve esta etapa
-----------------------
El paper no solo propone NCC-SA: mide cuanto mejor es que los metodos de
siempre. Aca estan esos metodos, para medirles el sesgo con verdad conocida
y sobre nuestros datos, contra la D del ajuste con plantilla (Etapa 8).

Los cuatro metodos
------------------
1. **50 % de altura**: el canal donde el espectro suavizado cae a la mitad de
   su maximo, pasado el maximo. El mas ingenuo.
2. **Diferencia**: el canal donde la derivada del espectro suavizado es mas
   negativa (la caida mas empinada), refinado sub-bin con una parabola.
3. **Escalon gaussiano**: se ajusta C/2 erfc((ch - B)/(sqrt2 s)) + a + b ch.
   El erfc es EXACTAMENTE un escalon convolucionado con una gaussiana, asi
   que B es su punto de media altura y s la resolucion; el termino lineal
   absorbe que el plateau no es plano. Es la lectura de `codigo/etapas`.
4. **Gaussiana del paper**: lo que el paper escribe literalmente ("Gaussian
   FWHM method"): una gaussiana ajustada al PICO de la acumulacion antes del
   borde, y el borde en mu + s sqrt(2 ln 2), el lado descendente a media
   altura. Es la lectura de `replica_paper`.

Si el borde fuera un escalon perfecto borroneado, 1, 2 y 3 darian lo mismo.
No lo es: Klein-Nishina tiene una acumulacion creciente ANTES del borde, y
eso corre la media altura y el minimo de la derivada hacia ARRIBA. Todos
sobrestiman el borde, y mas cuanto peor la resolucion (§7.44). Esa es la
justificacion cuantitativa de ajustar la forma entera.

Dos decisiones que importan (§7.44)
-----------------------------------
* **La ventana de suavizado sale del ancho del flanco**, no es una constante:
  la derivada del borde es un pozo de ancho ~ sigma, y con 21 bines fijos, a
  resolucion mala aparecen varios minimos casi iguales y el ruido elige. El
  ancho se mide sobre el ESPECTRO (caida del 80 % al 20 %, ~1.68 sigma), no
  sobre la derivada, y la ventana es 0.75 de eso.
* **La region se ancla al borde del NCC** (B +- 45 %), no a la ventana de
  ajuste: debajo del borde esta el exceso de bajo canal, cuya caida es mucho
  mas empinada que el borde, y el minimo de la derivada se va ahi (§7.44:
  D = 0.36 en vez de 1.27).

Robustez (Fig. 5)
-----------------
Seis intervalos anidados alrededor del borde, en unidades de B, los mismos
que `replica_paper` saco de los seis paneles del paper: del mas angosto
[0.977, 1.042] B al mas ancho [0.521, 1.433] B. Cada metodo se corre SOLO
dentro de cada intervalo; un metodo robusto da el mismo borde en todos. Con
menos de 30 bines un intervalo se saltea para todos igual.
"""

import os
import sys

import numpy as np
from scipy.optimize import curve_fit, least_squares
from scipy.signal import savgol_filter
from scipy.special import erfc

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ORDEN = 3
FACTOR_VENTANA = 0.75
ANCHO_REGION = 0.45
INTERVALOS_FIG5 = [(0.977, 1.042), (0.912, 1.042), (0.782, 1.042),
                   (0.782, 1.173), (0.651, 1.303), (0.521, 1.433)]
MIN_BINES = 30
METODOS = ("50 % de altura", "diferencia", "escalon gaussiano", "gaussiana del paper")


def suavizar(y, ventana, orden=ORDEN):
    v = int(ventana) | 1
    v = max(v, orden + 2 + (orden % 2))
    v = min(v, (len(y) - 1) | 1)
    return savgol_filter(np.asarray(y, float), v, orden)


def _mascara(x, region):
    return (x >= region[0]) & (x <= region[1])


def ventana_auto(x, y, region):
    """
    Bines del Savitzky-Golay: 0.75 x la caida 80 % -> 20 % del flanco,
    medida sobre el espectro apenas suavizado (5 % de la region).
    """
    idx = np.flatnonzero(_mascara(x, region))
    n = max(idx.size, 8)
    v_max = max(int(0.5 * n) | 1, ORDEN + 3)
    v0 = max(int(0.05 * n) | 1, ORDEN + 3)
    s = suavizar(y, v0)
    tope = int(idx[np.argmax(s[idx])])
    cola = idx[idx > tope]
    if cola.size < 5:
        return v0
    piso, alto = float(np.percentile(s[cola], 10)), float(s[tope])
    if alto <= piso:
        return v0
    b80 = cola[s[cola] <= piso + 0.8 * (alto - piso)]
    b20 = cola[s[cola] <= piso + 0.2 * (alto - piso)]
    if b80.size == 0 or b20.size == 0:
        return v0
    ancho = max(int(b20[0]) - int(b80[0]), 3)          # ~1.68 sigma en bines
    return int(np.clip(int(FACTOR_VENTANA * ancho) | 1, ORDEN + 3, v_max))


# --------------------------------------------------------------------------
def mitad_altura(x, y, region, ventana="auto"):
    v = ventana_auto(x, y, region) if ventana == "auto" else int(ventana)
    s = suavizar(y, v)
    m = _mascara(x, region)
    xs, ys = x[m], s[m]
    i = int(np.argmax(ys))
    j = np.flatnonzero(ys[i:] <= 0.5 * ys[i])
    return float(xs[i + j[0]]) if j.size else float("nan")


def diferencia(x, y, region, ventana="auto"):
    v = ventana_auto(x, y, region) if ventana == "auto" else int(ventana)
    s = suavizar(y, v)
    paso = float(x[1] - x[0])
    d = np.gradient(s, paso)
    idx = np.flatnonzero(_mascara(x, region))
    i = int(idx[np.argmin(d[idx])])
    borde = x[i]
    if 0 < i < x.size - 1:
        # vertice de la parabola por los tres puntos: saca al metodo del piso
        # de un bin, que si no seria un sesgo artificial en la comparacion
        den = d[i - 1] - 2 * d[i] + d[i + 1]
        if abs(den) > 1e-12:
            borde += paso * float(np.clip(0.5 * (d[i - 1] - d[i + 1]) / den, -1, 1))
    return float(borde)


def escalon(ch, B, s, C, a, b):
    ch = np.asarray(ch, float)
    return 0.5 * C * erfc((ch - B) / (np.sqrt(2) * max(s, 1e-9))) + a + b * (ch - B)


def escalon_gaussiano(x, y, region, errores=None):
    m = _mascara(x, region)
    xx, yy = x[m], y[m]
    B0 = 0.5 * (region[0] + region[1])
    alto = float(np.percentile(yy, 90) - np.percentile(yy, 10))
    p0 = [B0, max(0.08 * B0, 1e-3), max(alto, 1.0), float(np.percentile(yy, 10)), 0.0]
    cotas = ([xx[0], 1e-6, 0, -np.inf, -np.inf], [xx[-1], xx[-1] - xx[0], np.inf,
                                                  np.inf, np.inf])
    sy = None if errores is None else errores[m]
    p, _ = curve_fit(escalon, xx, yy, p0=p0, sigma=sy, absolute_sigma=sy is not None,
                     bounds=cotas, maxfev=20000)
    return float(p[0])


def gaussiana_del_paper(x, y, region):
    m = _mascara(x, region)
    xx, yy = x[m], y[m]
    i = int(np.argmax(yy))
    p0 = [max(yy[i], 1.0), xx[i], max((xx[-1] - xx[0]) / 4.0, 1.0)]
    r = least_squares(lambda p: p[0] * np.exp(-(xx - p[1]) ** 2 / (2 * p[2] ** 2)) - yy,
                      p0, method="lm", max_nfev=20000)
    A, mu, s = r.x[0], r.x[1], abs(r.x[2])
    borde = mu + s * np.sqrt(2 * np.log(2))
    # sin maximo adentro de la region (intervalo solo sobre el flanco) la
    # gaussiana se va al infinito: eso no es un borde
    return float(borde) if A > 0 and xx[0] <= borde <= xx[-1] else float("nan")


def los_cuatro(x, y, region, errores=None):
    """{metodo: canal del borde} dentro de `region`. nan si el metodo no pudo."""
    out = {}
    for nombre, f in zip(METODOS, (mitad_altura, diferencia, None, gaussiana_del_paper)):
        try:
            out[nombre] = (escalon_gaussiano(x, y, region, errores) if f is None
                           else f(x, y, region))
        except Exception:                                   # noqa: BLE001
            out[nombre] = float("nan")
    return out


def region_anclada(B, x):
    return (max((1 - ANCHO_REGION) * B, x[0]), min((1 + ANCHO_REGION) * B, x[-1]))


# --------------------------------------------------------------------------
def robustez(x, y, errores, B, fuente, plantilla=None):
    """
    Fig. 5: los cuatro metodos y nuestro ajuste, cada uno solo dentro de cada
    intervalo. Devuelve [(intervalo, bines, {metodo: borde})], y B el canal
    del borde en el que se anclan los intervalos.
    """
    from desde_cero import ajuste, ncc
    E_ref = float(ncc.bordes_de(fuente)[0])
    paso = float(x[1] - x[0])
    filas = []
    for a, b in INTERVALOS_FIG5:
        iv = (a * B, b * B)
        nb = int(round((iv[1] - iv[0]) / paso))
        if nb < MIN_BINES:
            filas.append((iv, nb, None))
            continue
        r = los_cuatro(x, y, iv, errores)
        w = dict(lo=iv[0], hi=iv[1], D0=B / E_ref, sigma_rel=0.1,
                 cota_D=(0.6 * B / E_ref, 1.6 * B / E_ref), cota_sigma=(0.02, 0.30))
        q = ajuste.calibrar(x, y, errores, fuente, ventana=w, plantilla=plantilla,
                            n_arranques=4)
        r["ajuste Geant4"] = q["D"] * E_ref
        filas.append((iv, nb, r))
    return filas


# --------------------------------------------------------------------------
def sesgo_sintetico(sigmas=(0.03, 0.06, 0.09, 0.12, 0.18), ruido=True, n=10,
                    semilla=3, A=2e5, D=1.0, paso=2.5):
    """
    Borde de Klein-Nishina (Cs-137) de posicion CONOCIDA, a varias
    resoluciones, con y sin ruido Poisson. Devuelve {metodo: [sesgo relativo
    medio por sigma]}. La region se ancla en el B verdadero.
    """
    from desde_cero import fisica, modelo
    rng = np.random.default_rng(semilla)
    B = D * float(fisica.e_max(661.657))
    x = np.arange(0.5 * paso, 2.2 * B, paso)
    out = {m: [] for m in METODOS}
    for s in sigmas:
        mu = modelo.modelo(x, A, D, s, "Cs-137")
        tiradas = [rng.poisson(mu).astype(float) for _ in range(n)] if ruido else [mu]
        acum = {m: [] for m in METODOS}
        for y in tiradas:
            e = np.sqrt(np.maximum(y, 1.0))
            for m, v in los_cuatro(x, y, region_anclada(B, x), e).items():
                acum[m].append(v / B - 1)
        for m in METODOS:
            out[m].append(float(np.nanmean(acum[m])))
    return out


# --------------------------------------------------------------------------
def _demo():
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero import espectro, lote, ncc, plantilla
    from desde_cero.lector import CORRIDAS_TRABAJO

    sig = (0.03, 0.06, 0.09, 0.12, 0.15, 0.18)
    con, sin = sesgo_sintetico(sig), sesgo_sintetico(sig, ruido=False)
    print("  sesgo del borde, Klein-Nishina de posicion conocida (con ruido / sin ruido):")
    print("  %-20s" % "sigma/B" + "".join("%13.2f" % s for s in sig))
    for m in METODOS:
        print("  %-20s" % m + "".join("%+6.1f/%+5.1f%%" % (100 * a, 100 * b)
                                      for a, b in zip(con[m], sin[m])))

    fuentes = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
    lot = json.load(open(lote.SALIDA, encoding="utf-8")) if os.path.exists(lote.SALIDA) else {}
    fig, ejes = plt.subplots(1, 3, figsize=(18, 5))
    ax = ejes[0]
    for m in METODOS:
        l, = ax.plot(sig, 100 * np.array(con[m]), "o-", ms=4, label=m)
        ax.plot(sig, 100 * np.array(sin[m]), ":", color=l.get_color())
    ax.axhline(0, color="k", lw=0.8)
    ax.set(xlabel="sigma / B", ylabel="borde medido / verdadero - 1 [%]",
           title="sesgo con verdad conocida (lleno: con ruido; punteado: sin)")
    ax.legend(fontsize=8)

    print()
    print("  sobre los datos, D de cada metodo contra la del ajuste con Geant4 (Etapa 8):")
    rob_datos = {}
    for clave, f in fuentes.items():
        if not os.path.exists(CORRIDAS_TRABAJO[clave][0]):
            continue
        e, _, _ = espectro.de_corrida(clave)
        x, y, err = e["canales"], e["neto"], e["error"]
        E_ref = float(ncc.bordes_de(f)[0])
        pl = plantilla.Plantilla(f)
        D_aj = lot[clave]["D"] if clave in lot else None
        if D_aj is None:
            from desde_cero import ajuste
            D_aj = ajuste.calibrar(x, y, err, f, plantilla=pl)["D"]
        B_ncc = ncc.barrido_B(x, y, f, plantilla=pl)["B"]
        r = los_cuatro(x, y, region_anclada(B_ncc, x), err)
        print("  %-12s ajuste %.4f | " % (clave, D_aj)
              + " | ".join("%s %+.1f %%" % (m, 100 * (v / E_ref / D_aj - 1))
                           for m, v in r.items()))
        rob_datos[clave] = robustez(x, y, err, D_aj * E_ref, f, pl)

    for ax, (clave, filas) in zip(ejes[1:], list(rob_datos.items())[:2]):
        ok = [(i, fila) for i, fila in enumerate(filas) if fila[2] is not None]
        B = lot.get(clave, {}).get("D", np.nan) * float(ncc.bordes_de(fuentes[clave])[0])
        for m in METODOS + ("ajuste Geant4",):
            ax.plot([i + 1 for i, _ in ok], [100 * (f[2][m] / B - 1) for _, f in ok],
                    "o-", ms=4, label=m)
        ax.set(xlabel="intervalo de la Fig. 5 (1 = el mas angosto)",
               ylabel="borde / B del ajuste completo - 1 [%]",
               title="%s: robustez frente al intervalo" % clave, xticks=range(1, 7))
        ax.legend(fontsize=7)
        print("  %-12s robustez (dispersion max-min entre intervalos):" % clave)
        for m in METODOS + ("ajuste Geant4",):
            v = np.array([f[2][m] for _, f in ok])
            v = v[np.isfinite(v)]
            print("      %-20s %5.1f %% (%d intervalos)"
                  % (m, 100 * np.ptp(v) / np.mean(v) if v.size else np.nan, v.size))
    for ax in ejes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa10_clasicos.png")
    fig.savefig(salida, dpi=110)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
