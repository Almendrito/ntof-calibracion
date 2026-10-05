"""
Etapa 6 — la ventana automatica: correlacion normalizada (paper, Seccion II-B,
Ec. 15).

    python desde_cero/ncc.py

Que resuelve esta etapa
-----------------------
Para ajustar el borde (Etapa 7) hay que decir EN QUE TRAMO del espectro.
Elegirlo a mano es subjetivo, y el resultado depende de la eleccion (§7.32:
3.5 % de D moviendo el limite inferior). El paper lo decide solo: compara la
FORMA del dato contra la plantilla teorica en muchos tramos candidatos y se
queda con el que mejor correlaciona. De paso entrega una primera D, la
SEMILLA del ajuste, y de ahi las cotas donde el ajuste va a buscar.

Pearson (Ec. 15)
----------------
    rho = sum (a - <a>)(b - <b>) / sqrt( sum (a - <a>)^2 * sum (b - <b>)^2 )

Restar la media borra cualquier pedestal; dividir por la norma borra
cualquier escala. Queda la FORMA y nada mas: rho = 1 si una curva es la otra
por una constante mas otra constante. Por eso sirve aca, donde la plantilla
no sabe cuantas cuentas hay (A) y el dato puede tener un fondo plano. Y por
eso mismo es peligroso: a un tramo de ruido puro tambien se lo puede llevar a
parecerse a cualquier cosa (§7.11).

Los dos barridos
----------------
* `barrido_libre` — el del PAPER, literal: la plantilla fija de 100 puntos
  entre 0.45 y ~1.3 veces el borde; cada bin no vacio es un inicio posible y
  cada bin posterior un fin posible; el tramo de dato se reinterpola a 100
  puntos. Dos parametros (inicio, fin).
* `barrido_B` — el NUESTRO: un solo parametro, el canal del borde B. Por que
  alcanza con uno: la Ec. 8 dice que el mapa pasa por el ORIGEN, ch = D*E.
  Fijado B (y con el, D = B/E_ref), la ventana ya esta determinada:
  [0.45 B, f B]. Barrer inicio y fin por separado le da a la busqueda un
  grado de libertad que el modelo no tiene, y en datos reales se lo gasto en
  ventanas absurdas (§7.10: una de 32 bines en la cola vacia con rho = 0.984).
  `f` deja entrar el ULTIMO borde mas 3 sigma: en el Co-60 una ventana que
  corte en el primer borde perderia justo lo que distingue a los dos.

Un punto por bin, no 100
------------------------
El paper reinterpola a 100 puntos. Con pocos puntos la interpolacion promedia
bines vecinos y eso PROMEDIA RUIDO: rho sube sin que la forma se parezca mas
(§7.27: rho 0.994 con 50 puntos contra 0.992 con uno por bin). Y rho es
justamente el criterio con el que se comparan ventanas. Con MAS puntos que
bines no se agrega informacion, solo puntos correlacionados. Un punto por bin
es "no inventar ni tirar nada". El numero se fija UNA vez para todo el
barrido: si cambiara con B, los rho de ventanas distintas no se podrian
comparar.

Suavizado Savitzky-Golay
------------------------
El barrido nuestro correlaciona contra el dato suavizado (ventana 21 bines,
polinomio de grado 3). Un promedio movil de 21 bines achataria el flanco y
correria el borde; Savitzky-Golay ajusta un polinomio local y conserva la
posicion y la pendiente de un flanco suave. Solo se usa para ELEGIR la
ventana: el ajuste de la Etapa 7 va sobre el dato crudo.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import fisica, modelo                       # noqa: E402

INICIO_REL = 0.45        # inicio de la ventana en unidades de B (el del paper)
SIGMA_REL = 0.09         # ancho de la plantilla (el paper no lo dice; ver abajo)
N_B = 240                # cuantos B se prueban
SUAVE = 21               # ventana del Savitzky-Golay [bines]
COTA_D = (0.6, 1.6)      # la busqueda de la Etapa 7, en veces D0
COTA_SIGMA = (0.02, 0.30)
N_PAPER = 100            # puntos del paper
MIN_BINES = 25           # largo minimo de un tramo en el barrido libre


def pearson(a, b):
    """Ec. 15."""
    a = np.asarray(a, float) - np.mean(a)
    b = np.asarray(b, float) - np.mean(b)
    d = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / d) if d > 0 else 0.0


def suavizar(y, ventana=SUAVE, orden=3):
    from scipy.signal import savgol_filter
    return savgol_filter(np.asarray(y, float), ventana, orden)


def bordes_de(fuente):
    E, _ = modelo.lineas_de(fuente)
    return np.sort(fisica.e_max(E))


def borde_aproximado(canales, y):
    """
    Una primera idea de donde esta el borde, SOLO para acotar el barrido:
    el primer canal, pasada la mitad de las cuentas, donde el espectro
    suavizado cae a la mitad de su nivel tipico en esa zona.
    """
    pos = np.maximum(y, 0.0)
    cum = np.cumsum(pos) / max(pos.sum(), 1e-9)
    zona = (cum > 0.5) & (cum < 0.9)
    if zona.sum() < 5:
        return float(canales[len(canales) // 2])
    alto = float(np.median(y[zona]))
    i0 = int(np.flatnonzero(zona)[-1])
    j = np.flatnonzero((y < 0.5 * alto) & (np.arange(y.size) > i0))
    return float(canales[j[0]]) if j.size else float(canales[i0])


# --------------------------------------------------------------------------
# El nuestro: un parametro
# --------------------------------------------------------------------------
def barrido_B(canales, cuentas, fuente, sigma_rel=SIGMA_REL, inicio_rel=INICIO_REL,
              n_B=N_B, suave=SUAVE, B_rango=None, n_pts=None, por_bin=True,
              plantilla=None):
    """
    Barre el canal del borde B. Para cada B: D = B/E_ref, la plantilla es
    `modelo.modelo(canales, 1, D, sigma_rel)` y la ventana [inicio_rel B,
    fin_rel B]; rho entre dato suavizado y plantilla en la ventana.

    Devuelve un dict: B, D0 (= semilla), rho, lo, hi, cotas de D y de
    sigma_rel para la Etapa 7, y las curvas del barrido (Bs, rhos).

    `por_bin`: plantilla integrada por bin (default, la misma del ajuste) o
    evaluada en el centro del bin (la de analisis/, para compararse con el).
    Mueve B ~0.3-0.6 % en los datos reales: es semilla, no calibracion.
    `plantilla`: None = Klein-Nishina; o la de Geant4 (`plantilla.py`).
    """
    ch = np.asarray(canales, float)
    bordes = bordes_de(fuente)
    E_ref = float(bordes[0])
    fin_rel = float(bordes[-1] / E_ref + 3.0 * sigma_rel)
    y = suavizar(cuentas, suave)
    if B_rango is None:
        B0 = borde_aproximado(ch, y)
        B_rango = (0.40 * B0, 1.10 * B0)
    Bs = np.linspace(B_rango[0], B_rango[1], n_B)
    paso = float(ch[1] - ch[0])
    if n_pts is None:
        # un punto por bin de la ventana MEDIA del barrido, fijo para todos
        n_pts = max(int((fin_rel - inicio_rel) * Bs.mean() / paso), 24)

    rhos = np.full(n_B, -1.0)
    for i, B in enumerate(Bs):
        lo, hi = inicio_rel * B, fin_rel * B
        m = (ch >= lo) & (ch <= hi)
        if m.sum() < 12:
            continue
        p = modelo.modelo(ch, 1.0, B / E_ref, sigma_rel, fuente, por_bin=por_bin,
                          plantilla=plantilla)
        u = np.linspace(lo, hi, n_pts)
        rhos[i] = pearson(np.interp(u, ch[m], y[m]), np.interp(u, ch[m], p[m]))

    j = int(np.argmax(rhos))
    B = float(Bs[j])
    D0 = B / E_ref
    return dict(B=B, D0=D0, rho=float(rhos[j]), lo=inicio_rel * B,
                hi=fin_rel * B, E_ref=E_ref, sigma_rel=sigma_rel,
                inicio_rel=inicio_rel, fin_rel=fin_rel, n_pts=int(n_pts),
                cota_D=(COTA_D[0] * D0, COTA_D[1] * D0), cota_sigma=COTA_SIGMA,
                Bs=Bs, rhos=rhos, suavizado=y)


# --------------------------------------------------------------------------
# El del paper: dos parametros
# --------------------------------------------------------------------------
def plantilla_paper(fuente, sigma_rel=0.06, n=N_PAPER, ini_rel=INICIO_REL,
                    fin_rel=None):
    """
    S2 del paper: `n` puntos uniformes de la forma entre ini_rel y fin_rel
    veces el borde, con fin_rel = 1 + 4.5 sigma_rel por defecto (el paper
    dice "hasta el minimo tras el borde" y no da el numero; replica_paper
    usa este). Se construye con D = 1 (canal = keV): Pearson no ve escalas.
    """
    E_ref = float(bordes_de(fuente)[0])
    if fin_rel is None:
        fin_rel = 1.0 + 4.5 * sigma_rel
    x = np.linspace(ini_rel, fin_rel, int(n))
    # por_bin=False: el paper evalua la forma en cada punto, no la integra
    return x, modelo.modelo(x * E_ref, 1.0, 1.0, sigma_rel, fuente,
                            por_bin=False), fin_rel


def barrido_libre(canales, cuentas, fuente, sigma_rel=0.06, n=N_PAPER,
                  ini_rel=INICIO_REL, fin_rel=None, min_bines=MIN_BINES, n_max=220):
    """
    El barrido de la Seccion II-B, literal: cada bin NO NULO es un inicio
    posible y cada bin a >= `min_bines` un fin posible; el tramo se
    reinterpola a `n` puntos y se correlaciona con S2. Sin suavizar (el paper
    no suaviza).

    Si el espectro tiene mas de `n_max` bines, inicios y fines se prueban
    cada ceil(N/n_max) bines (O(N^2) si no). El borde se deduce del span de
    la plantilla: si S2 va de ini a fin veces el borde y el tramo ganador es
    [lo, hi], el borde esta en lo + (1 - ini)/(fin - ini) * (hi - lo).
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    _, s2, fin_rel = plantilla_paper(fuente, sigma_rel, n, ini_rel, fin_rel)
    s2c = s2 - s2.mean()
    n2 = np.sqrt(np.sum(s2c * s2c))
    u = np.linspace(0.0, 1.0, int(n))
    salto = max(1, int(np.ceil(ch.size / float(n_max))))
    mejor, mapa = (-2.0, -1, -1), []
    for i in np.flatnonzero(y != 0.0)[::salto]:
        for j in range(i + min_bines, y.size, salto):
            t = y[i:j + 1]
            s1 = np.interp(u, np.linspace(0.0, 1.0, t.size), t)
            s1c = s1 - s1.mean()
            n1 = np.sqrt(np.sum(s1c * s1c))
            if n1 <= 0:
                continue
            rho = float(s1c @ s2c / (n1 * n2))
            mapa.append((ch[i], ch[j], rho))
            if rho > mejor[0]:
                mejor = (rho, int(i), int(j))
    rho, i, j = mejor
    lo, hi = float(ch[i]), float(ch[j])
    B = lo + (1.0 - ini_rel) / (fin_rel - ini_rel) * (hi - lo)
    return dict(lo=lo, hi=hi, rho=rho, B=float(B),
                D0=float(B) / float(bordes_de(fuente)[0]),
                mapa=np.array(mapa), salto=salto, fin_rel=fin_rel)


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero import espectro
    from desde_cero.lector import CORRIDAS_TRABAJO

    fuentes = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
    claves = [k for k in fuentes if os.path.exists(CORRIDAS_TRABAJO[k][0])]
    fig, ejes = plt.subplots(3, len(claves), figsize=(6 * len(claves), 11),
                             squeeze=False)
    for col, clave in enumerate(claves):
        e, _, _ = espectro.de_corrida(clave)
        x, y = e["canales"], e["neto"]
        r = barrido_B(x, y, fuentes[clave])
        lib = barrido_libre(x, y, fuentes[clave])
        print("  %-12s nuestro: B = %.1f, D0 = %.4f mV*ns/keV, rho = %.5f, "
              "ventana [%.0f, %.0f]" % (clave, r["B"], r["D0"], r["rho"],
                                        r["lo"], r["hi"]))
        print("  %-12s paper:   B = %.1f, D0 = %.4f,          rho = %.5f, "
              "ventana [%.0f, %.0f] (%d bines)"
              % ("", lib["B"], lib["D0"], lib["rho"], lib["lo"], lib["hi"],
                 int(round((lib["hi"] - lib["lo"]) / (x[1] - x[0])))))

        ax = ejes[0, col]
        ax.plot(r["Bs"], r["rhos"], lw=1)
        ax.axvline(r["B"], color="tab:red", lw=0.8)
        ax.set(xlabel="B [mV*ns]", ylabel="rho", title="%s: rho(B)" % clave,
               ylim=(max(r["rhos"].min(), 0.5), 1.0))

        ax = ejes[1, col]
        ax.step(x, y, where="mid", color="0.6", lw=0.7, label="neto")
        ax.plot(x, r["suavizado"], color="k", lw=1, label="suavizado (SG 21)")
        ax.axvspan(r["lo"], r["hi"], color="tab:blue", alpha=0.12,
                   label="ventana nuestra")
        ax.axvspan(lib["lo"], lib["hi"], color="tab:orange", alpha=0.20,
                   label="ventana del paper")
        ax.axvline(r["B"], color="tab:blue", ls="--", lw=1)
        ax.axvline(lib["B"], color="tab:orange", ls="--", lw=1)
        ax.set(xlabel="area [mV*ns]", ylabel="cuentas", title="ventanas elegidas",
               xlim=(0, 1.6 * r["hi"]))
        ax.legend(fontsize=7)

        ax = ejes[2, col]
        m = (x >= r["lo"]) & (x <= r["hi"])
        p = modelo.modelo(x, 1.0, r["D0"], r["sigma_rel"], fuentes[clave])
        z = lambda v: (v - v.mean()) / v.std()              # noqa: E731
        ax.plot(x[m], z(r["suavizado"][m]), lw=1, label="dato suavizado")
        ax.plot(x[m], z(p[m]), lw=1, label="plantilla (Klein-Nishina)")
        ax.set(xlabel="area [mV*ns]", ylabel="estandarizado",
               title="lo que compara Pearson (rho = %.4f)" % r["rho"])
        ax.legend(fontsize=7)
    for ax in ejes.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa6_ncc.png")
    fig.savefig(salida, dpi=110)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
