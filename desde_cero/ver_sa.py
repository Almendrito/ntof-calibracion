"""
Ver el recocido simulado trabajando sobre los datos reales.

    python desde_cero/ver_sa.py

Esto NO es una etapa de la construccion desde cero (el SA vendria en la etapa 5,
despues del area y del espectro). Es una MIRADA al SA que ya usa `analisis/`,
para entender que hace antes de reimplementarlo.

Muestra cuatro cosas, en orden:

  1. QUE se optimiza: el vector Theta y la funcion de costo.
  2. COMO se ve el paisaje de esa funcion.
  3. POR QUE hace falta el recocido: un optimizador local desde un mal arranque
     contra el recocido desde EL MISMO arranque.
  4. COMO se sabe que el resultado es bueno: multiarranque.

El documento largo, con las 7 verificaciones, esta en `docs/recocido_simulado.pdf`
(lo genera `analisis/doc_sa.py`, y las pruebas estan en `analisis/verificar_sa.py`).
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scipy.optimize import minimize                          # noqa: E402
from analisis import modelo, rutas, sa                       # noqa: E402


def _titulo(n, t):
    print()
    print("=" * 74)
    print("%d. %s" % (n, t))
    print("=" * 74)


def main():
    # El problema: Co-60, sin el termino de exceso, para que Theta tenga solo
    # DOS parametros y el paisaje se pueda dibujar en un plano. En produccion
    # son tres (se agrega `lam`, la caida del exceso de bajo canal).
    p, semilla, _, _, _ = modelo.armar(rutas.corrida("co60"), verbose=False,
                                       exceso=False)
    cotas = [(0.6 * semilla[0], 1.6 * semilla[0]), (0.02, 0.30)]

    # ---------------------------------------------------------------------
    _titulo(1, "QUE se optimiza")
    print("  Theta      = %s" % p.libres)
    print("  semilla    = D = %.4f mV*ns/keVee, sigma_rel = %.4f   (la da el NCC)"
          % tuple(semilla))
    print("  cotas      = %s" % [tuple(round(v, 3) for v in c) for c in cotas])
    print("  ventana    = [%.0f, %.0f] mV*ns, %d bins, %d gdl"
          % (p.lo, p.hi, p.dato.size, p.ndf))
    print()
    print("  costo(x) = chi2 = suma[ (dato - A*modelo(x))^2 / error^2 ]")
    print("  La AMPLITUD A no esta en Theta: para D y sigma fijos, el chi2 es")
    print("  cuadratico en A y su optimo se despeja de una (minimos cuadrados")
    print("  lineales). Eso se llama VARPRO y saca una dimension del problema")
    print("  gratis: el recocido busca en 2, no en 3.")
    print("  costo(semilla) = %.2f" % p.costo(semilla))

    # ---------------------------------------------------------------------
    _titulo(2, "COMO se ve el paisaje")
    Ds = np.linspace(cotas[0][0], cotas[0][1], 90)
    Ss = np.linspace(cotas[1][0], cotas[1][1], 70)
    t0 = time.time()
    Z = p.mapa(Ds, Ss)
    print("  %d x %d evaluaciones del modelo en %.1f s" % (len(Ds), len(Ss),
                                                           time.time() - t0))
    i, j = np.unravel_index(np.argmin(Z), Z.shape)
    print("  minimo de la grilla: D = %.4f, sigma_rel = %.4f, chi2 = %.2f"
          % (Ds[j], Ss[i], Z[i, j]))
    print("  chi2 va de %.0f a %.0f: el valle es angosto y el resto es pared."
          % (Z.min(), Z.max()))

    # ---------------------------------------------------------------------
    _titulo(3, "HACE FALTA el recocido en ESTE problema?")
    # La pregunta no se contesta con una anecdota: hay que barrer arranques.
    # Se prueban 25 puntos de partida repartidos por toda la caja de cotas y se
    # cuenta desde cuantos llega cada metodo al minimo global.
    ref = float(min(Z.min(), sa.pulir(p.costo, [Ds[j], Ss[i]])[1]))
    partidas = [(d, sg) for d in np.linspace(cotas[0][0], cotas[0][1], 5)
                for sg in np.linspace(cotas[1][0], cotas[1][1], 5)]
    fallos_loc, fallos_sa, t_loc, t_sa = [], [], 0.0, 0.0
    for x0 in partidas:
        t = time.time()
        l = minimize(p.costo, np.array(x0), method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-4, maxiter=4000))
        t_loc += time.time() - t
        t = time.time()
        rr = sa.recocido(p.costo, x0, cotas, semilla=1)
        _, cc = sa.pulir(p.costo, rr["x"])
        t_sa += time.time() - t
        # "llego" = quedo a menos de 1 unidad de chi2 del mejor conocido
        if l.fun > ref + 1.0:
            fallos_loc.append((x0, float(l.fun), float(l.x[0])))
        if cc > ref + 1.0:
            fallos_sa.append((x0, cc))

    n = len(partidas)
    print("  mejor chi2 conocido: %.2f" % ref)
    print("  %d arranques repartidos por toda la caja de cotas:" % n)
    print("    Nelder-Mead solo : llega desde %d/%d   (%.2f s en total)"
          % (n - len(fallos_loc), n, t_loc))
    print("    recocido+pulido  : llega desde %d/%d   (%.2f s en total, %.0fx mas)"
          % (n - len(fallos_sa), n, t_sa, t_sa / max(t_loc, 1e-9)))
    if fallos_loc:
        print()
        print("  Donde falla el local (arranque -> a donde se fue):")
        for x0, c_, d_ in fallos_loc[:6]:
            fuera = "  <- FUERA de las cotas" if not (cotas[0][0] <= d_ <= cotas[0][1]) else ""
            print("    D0 = %.3f, sigma0 = %.3f  ->  chi2 = %.0f, D = %.4f%s"
                  % (x0[0], x0[1], c_, d_, fuera))
        print()
        print("  Fijate a donde se va: D fuera de la caja. Nelder-Mead NO conoce")
        print("  las cotas — se le da un punto inicial y camina hasta donde sea.")
        print("  El recocido si: cada propuesta pasa por np.clip(prop, lo, hi),")
        print("  asi que no puede salirse de lo fisicamente posible. Esa es una")
        print("  ventaja concreta suya en este problema, aparte de los valles.")
    print()
    print("  Y hay que decir lo que la medicion NO muestra: con Theta de 2")
    print("  parametros y VARPRO sacando la amplitud, el paisaje tiene un solo")
    print("  valle, asi que el recocido casi no esta corrigiendo minimos locales")
    print("  — es un seguro. Se gano su lugar en el modelo de SEIS parametros de")
    print("  bordes solapados (CLAUDE.md 7.21-7.23): ahi no mejoro el numero,")
    print("  DEMOSTRO que el modelo estaba sub-determinado, porque el optimo")
    print("  global que encontro no tenia sentido fisico.")
    print()

    # una corrida sola, para mirarle el enfriamiento por dentro
    malo = np.array([0.55 * semilla[0], 3.0 * semilla[1]])
    r = sa.recocido(p.costo, malo, cotas, semilla=1)
    x_sa, c_sa = sa.pulir(p.costo, r["x"])
    loc = minimize(p.costo, malo, method="Nelder-Mead",
                   options=dict(xatol=1e-7, fatol=1e-4, maxiter=4000))
    print("  Una corrida por dentro, desde D0 = %.3f (%.0f%% del real), sigma0 = %.3f:"
          % (malo[0], 100 * malo[0] / semilla[0], malo[1]))
    print("  El recocido hizo %d evaluaciones en %d ciclos de temperatura,"
          % (r["n_eval"], r["ciclos"]))
    print("  arrancando en T = %.1f y cortando por %s." % (r["T_init"], r["motivo"]))
    print()
    print("  Como cae la temperatura y con ella la tasa de aceptacion:")
    h = r["historia"]
    print("    %-8s %-12s %-12s %s" % ("ciclo", "T", "mejor chi2", "aceptacion"))
    for k in np.linspace(0, len(h) - 1, min(8, len(h))).astype(int):
        print("    %-8d %-12.4g %-12.2f %.2f" % (k, h[k][0], h[k][2], h[k][3]))
    print("  ^ arranca aceptando ~%.0f%% (explora) y termina en ~%.0f%% (afina)."
          % (100 * h[0][3], 100 * h[-1][3]))

    # ---------------------------------------------------------------------
    _titulo(4, "COMO se sabe que el resultado es bueno")
    t0 = time.time()
    corridas = sa.multiarranque(p.costo, cotas, n=8, semilla=7)
    Ds_ms = np.array([c["x_pulido"][0] for c in corridas])
    cs_ms = np.array([c["costo_pulido"] for c in corridas])
    print("  8 arranques INDEPENDIENTES, desde puntos al azar dentro de las cotas")
    print("  (%.0f s):" % (time.time() - t0))
    for k in range(len(corridas)):
        print("    arranque %d: partio de D = %.3f  ->  D = %.4f   chi2 = %.2f"
              % (k + 1, corridas[k]["x0"][0], Ds_ms[k], cs_ms[k]))
    print()
    print("  D = %.4f +- %.4f   (dispersion relativa %.3f %%)"
          % (Ds_ms.mean(), Ds_ms.std(), 100 * Ds_ms.std() / Ds_ms.mean()))
    print("  Que 8 busquedas independientes caigan en el mismo punto es la")
    print("  evidencia practica de que ese punto es el minimo GLOBAL. No es una")
    print("  demostracion, pero es lo unico que se puede afirmar desde adentro.")

    # ---------------------------------------------------------------------
    _figura(p, Ds, Ss, Z, r, malo, loc.x, x_sa, corridas)


def _figura(p, Ds, Ss, Z, r, malo, x_loc, x_sa, corridas):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15.5, 4.3))

    # --- el paisaje, en escala log para que se vea el valle ---
    Zr = np.log10(Z - Z.min() + 1.0)
    im = a1.contourf(Ds, Ss, Zr, levels=30, cmap="viridis")
    fig.colorbar(im, ax=a1, label="log10(chi2 - chi2_min + 1)")
    a1.plot(*malo, "*", color="w", mec="k", mew=1.2, ms=17,
            label="arranque malo")
    a1.plot(*x_loc, "X", color="#d62728", mec="w", mew=1.0, ms=12,
            label="Nelder-Mead")
    a1.plot(*x_sa, "o", color="#ff7f0e", mec="w", mew=1.0, ms=8,
            label="recocido + pulido")
    a1.set_xlabel("D [mV·ns/keVee]")
    a1.set_ylabel("sigma_rel")
    a1.set_title("el paisaje del chi2: UN valle", fontsize=10)
    a1.legend(fontsize=8, loc="lower right", framealpha=.9)

    # --- temperatura y aceptacion ---
    h = r["historia"]
    ciclos = np.arange(len(h))
    a2.semilogy(ciclos, h[:, 0], color="C3", label="temperatura")
    a2.semilogy(ciclos, h[:, 2] - h[:, 2].min() + 1e-3, color="C0",
                label="mejor chi2 - min")
    a2.set_xlabel("ciclo de temperatura")
    a2.set_title("enfriamiento geometrico  T <- 0.95 T", fontsize=10)
    a2.legend(fontsize=8)
    a2.grid(alpha=.3)

    a3b = a2.twinx()
    a3b.plot(ciclos, h[:, 3], color="0.5", lw=1, ls="--")
    a3b.set_ylabel("tasa de aceptacion", color="0.4")
    a3b.set_ylim(0, 1)

    # --- multiarranque: DE DONDE partio cada uno y A DONDE llego ---
    # Un histograma no sirve aca: los 8 dan el mismo D y el eje queda en
    # notacion cientifica sobre un ancho de 1e-10. Lo que hay que ver es que el
    # punto de llegada NO depende del de partida.
    x0 = np.array([c["x0"][0] for c in corridas])
    xf = np.array([c["x_pulido"][0] for c in corridas])
    a3.scatter(x0, xf, s=55, color="C0", zorder=3, label="cada arranque")
    a3.axhline(xf.mean(), color="C3", ls="--",
               label="llegada comun = %.4f" % xf.mean())
    lim = (min(x0.min(), xf.min()) * 0.95, max(x0.max(), xf.max()) * 1.05)
    a3.plot(lim, lim, color="0.8", lw=1, zorder=1, label="y = x (no se movio)")
    a3.set_xlim(*lim)
    a3.set_xlabel("D del punto de PARTIDA")
    a3.set_ylabel("D de LLEGADA")
    a3.set_title("8 arranques: la llegada no depende de la partida",
                 fontsize=10)
    a3.legend(fontsize=8)
    a3.grid(alpha=.3)

    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras")
    os.makedirs(salida, exist_ok=True)
    salida = os.path.join(salida, "sa_en_accion.png")
    fig.savefig(salida, dpi=130)
    print()
    print("  figura: %s" % salida)


if __name__ == "__main__":
    main()
