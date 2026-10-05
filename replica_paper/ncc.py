"""
ncc.py — Seccion II-B del paper: seleccion automatica del intervalo del borde.

El paper lo dice asi, y esta implementado literal:

    "we first generate 100 reference points (S2) by uniformly sampling the
     predicted lineshape from the Compton plateau to the edge minimum. For
     experimental data, each non-zero point serves as a potential edge start
     position, with subsequent points forming test intervals that are
     interpolated to 100 points (S1). The interval showing maximum Pearson
     correlation (rho) between S1 and S2 is identified as the optimal Compton
     edge region."

O sea: barrido de DOS parametros (inicio y fin del intervalo), plantilla fija de
100 puntos, Pearson (Ec. 15) como criterio.

Lo que el paper NO dice, y aca queda como parametro explicito
-------------------------------------------------------------
1. Con que sigma se construye S2. El paper habla del "predicted lineshape ...
   to the edge minimum": una Klein-Nishina SIN ensanchar no tiene descenso, se
   corta de golpe, asi que S2 tiene que estar convolucionada con ALGO. Pero
   sigma es justo lo que el SA ajusta despues. -> `sigma_rel`, y
   `barrer_sigma_rel` mide cuanto depende el resultado.
2. Donde empieza "the Compton plateau" y donde termina "the edge minimum". Esto
   NO es cosmetico: fija el mapa. Si S2 abarca [ini_rel*E_max, fin_rel*E_max] y
   el intervalo ganador es [lo, hi], entonces el borde cae en
       B = lo + (1 - ini_rel)/(fin_rel - ini_rel) * (hi - lo),
   asi que mover el span de la plantilla mueve B en proporcion.
3. El largo minimo de un intervalo candidato. Sin un minimo, Pearson (que ignora
   escala y offset) puede encontrar correlacion alta en un puñado de bins de
   cola vacia, donde solo hay ruido.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import ecuaciones as ec

N_PUNTOS = 100          # "interpolated to 100 points": el numero es del paper
INI_REL = 0.45          # inicio de S2, en unidades de E_max  [NO esta en el paper]
FIN_REL = None          # fin de S2; None = borde + 4.5*sigma  [NO esta en el paper]
MIN_BINS = 25           # largo minimo de un intervalo         [NO esta en el paper]


def plantilla(e_gamma, sigma_rel=0.06, n=N_PUNTOS, ini_rel=INI_REL,
              fin_rel=FIN_REL, mec2=ec.M_E):
    """
    S2: `n` puntos uniformes de la lineshape predicha (Ec. 8-10).

    Se construye en el propio eje de energia (B = E_max, o sea canal == keV),
    porque Pearson es invariante a escala y offset: lo unico que importa es la
    FORMA y el span relativo al borde.

    Devuelve (x_rel, s2) con x_rel en unidades de E_max, para que quien la use
    sepa donde cae el borde dentro de la plantilla.
    """
    emax = float(ec.e_max(e_gamma, mec2))
    if fin_rel is None:
        fin_rel = 1.0 + 4.5 * sigma_rel
    if not (0.0 < ini_rel < 1.0 < fin_rel):
        raise ValueError("se espera 0 < ini_rel < 1 < fin_rel")
    x_rel = np.linspace(ini_rel, fin_rel, int(n))
    y = ec.modelo(x_rel * emax, e_gamma, A=1.0, B=emax,
                  sigma=sigma_rel * emax, mec2=mec2)
    return x_rel, y


def buscar_intervalo(canales, cuentas, e_gamma, sigma_rel=0.06, n=N_PUNTOS,
                     ini_rel=INI_REL, fin_rel=FIN_REL, min_bins=MIN_BINS,
                     paso_inicio="auto", paso_fin="auto", n_max=220,
                     mec2=ec.M_E):
    """
    El barrido de la Seccion II-B, literal.

    Recorre todos los pares (inicio, fin) con inicio en un bin de cuentas NO
    NULAS (el unico filtro que pone el paper) y fin - inicio >= `min_bins`,
    reinterpola el tramo de datos a `n` puntos y se queda con el de mayor rho.

    Devuelve un dict con lo, hi, rho, B (canal del borde deducido del span de la
    plantilla), sigma (semilla, = sigma_rel*B) y el mapa rho(lo, hi) para poder
    dibujarlo.
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    if ch.size != y.size:
        raise ValueError("canales y cuentas tienen que tener el mismo largo")

    x_rel, s2 = plantilla(e_gamma, sigma_rel, n, ini_rel, fin_rel, mec2)
    s2c = s2 - s2.mean()
    n2 = np.sqrt(np.sum(s2c * s2c))
    u = np.linspace(0.0, 1.0, int(n))

    # El barrido es O(N^2) en bins. Con binning fino eso se va de las manos, y
    # muestrear el barrido no cambia el resultado porque rho varia suave con los
    # extremos. `n_max` fija cuantos inicios y cuantos finales se prueban.
    if paso_inicio == "auto":
        paso_inicio = max(1, int(np.ceil(ch.size / float(n_max))))
    if paso_fin == "auto":
        paso_fin = max(1, int(np.ceil(ch.size / float(n_max))))
    inicios = np.nonzero(y != 0.0)[0][::int(paso_inicio)]
    mejor = dict(rho=-2.0, i=-1, j=-1)
    mapa = []
    for i in inicios:
        jmax = y.size - 1
        for j in range(i + int(min_bins), jmax + 1, int(paso_fin)):
            tramo = y[i:j + 1]
            s1 = np.interp(u, np.linspace(0.0, 1.0, tramo.size), tramo)
            s1c = s1 - s1.mean()
            n1 = np.sqrt(np.sum(s1c * s1c))
            if n1 <= 0.0:
                continue
            rho = float(np.dot(s1c, s2c) / (n1 * n2))
            mapa.append((ch[i], ch[j], rho))
            if rho > mejor["rho"]:
                mejor = dict(rho=rho, i=int(i), j=int(j))

    if mejor["i"] < 0:
        raise ValueError("ningun intervalo candidato: revisar min_bins y los datos")

    lo, hi = float(ch[mejor["i"]]), float(ch[mejor["j"]])
    if fin_rel is None:
        fin_rel = 1.0 + 4.5 * sigma_rel
    frac_borde = (1.0 - ini_rel) / (fin_rel - ini_rel)
    B = lo + frac_borde * (hi - lo)
    return dict(lo=lo, hi=hi, i=mejor["i"], j=mejor["j"], rho=mejor["rho"],
                paso_inicio=int(paso_inicio), paso_fin=int(paso_fin),
                B=float(B), sigma=float(sigma_rel * B),
                sigma_rel=float(sigma_rel), ini_rel=float(ini_rel),
                fin_rel=float(fin_rel), n_puntos=int(n),
                x_rel=x_rel, s2=s2, mapa=np.array(mapa, float))


def barrer_sigma_rel(canales, cuentas, e_gamma, sigmas, **kw):
    """
    Cuanto se mueve el resultado con el sigma de la plantilla, que el paper no
    especifica. Devuelve una lista de dicts de `buscar_intervalo`.

    OJO CON COMPARAR LOS rho ENTRE SI: `fin_rel` depende de `sigma_rel`, asi que
    cambiar sigma_rel cambia el span de S2 y por lo tanto el conjunto de
    ventanas candidatas. Los rho de dos sigma_rel distintos NO se midieron sobre
    el mismo problema, y rho tiende a crecer con sigma_rel simplemente porque la
    plantilla se vuelve mas suave. Quedarse con el de mayor rho es el unico
    criterio que el paper ofrece, pero es un criterio confundido: por eso el
    informe imprime el barrido entero y no solo el ganador.
    """
    return [buscar_intervalo(canales, cuentas, e_gamma, sigma_rel=float(s), **kw)
            for s in sigmas]


def _demo():
    """Cierre con VERDAD CONOCIDA: espectro sintetico con B y sigma dados."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from replica_paper import rutas
    rutas.salida_utf8()
    rng = np.random.default_rng(7)

    eg = 661.657
    B_real, sigma_real, A_real = 600.0, 40.0, 1.0
    ch = np.linspace(2.5, 1200.0, 240)          # bins de 5
    limpio = ec.modelo(ch, eg, A_real, B_real, sigma_real)
    limpio = limpio / limpio.max() * 800.0
    datos = rng.poisson(np.maximum(limpio, 0.0)).astype(float)

    print("=" * 74)
    print("ncc.py - Seccion II-B sobre un espectro sintetico")
    print("=" * 74)
    print("  verdad: B = %.1f, sigma = %.1f (sigma/B = %.3f)"
          % (B_real, sigma_real, sigma_real / B_real))
    print()
    print("  sigma_rel de S2   rho        intervalo          B        error")
    res = None
    for s in (0.03, 0.05, 0.0667, 0.09, 0.12):
        r = buscar_intervalo(ch, datos, eg, sigma_rel=s)
        print("       %.4f       %.5f   [%6.0f, %6.0f]   %7.1f   %+6.1f"
              % (s, r["rho"], r["lo"], r["hi"], r["B"], r["B"] - B_real))
        if res is None or r["rho"] > res["rho"]:
            res = r
    print()
    print("  gana sigma_rel = %.4f con rho = %.5f -> B = %.1f (%+.2f%%)"
          % (res["sigma_rel"], res["rho"], res["B"],
             100 * (res["B"] / B_real - 1)))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].step(ch, datos, where="mid", color="0.5", lw=0.9, label="datos")
    ax[0].axvspan(res["lo"], res["hi"], color="tab:orange", alpha=0.18,
                  label="intervalo NCC")
    ax[0].axvline(B_real, color="k", ls="--", lw=1, label="B verdadero")
    ax[0].axvline(res["B"], color="tab:red", lw=1.2, label="B del NCC")
    ax[0].set_xlabel("canal")
    ax[0].set_ylabel("cuentas")
    ax[0].legend(fontsize=8)
    ax[0].set_title("intervalo elegido (rho = %.5f)" % res["rho"])

    tramo = datos[res["i"]:res["j"] + 1]
    u = np.linspace(0, 1, res["n_puntos"])
    s1 = np.interp(u, np.linspace(0, 1, tramo.size), tramo)
    z = lambda v: (v - v.mean()) / v.std()
    ax[1].plot(u, z(s1), lw=1.1, label="S1 (datos)")
    ax[1].plot(u, z(res["s2"]), lw=1.1, label="S2 (Ec. 8-10)")
    ax[1].set_xlabel("punto reinterpolado")
    ax[1].set_ylabel("estandarizado")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Ec. 15: S1 contra S2")
    fig.tight_layout()
    salida = rutas.figura("ncc_sintetico.png")
    fig.savefig(salida, dpi=130)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
