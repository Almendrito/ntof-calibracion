"""
ver_espectros.py — dibuja los espectros ya calculados, sin recalcular nada.

Lee los `.npz` que dejo `hacer_espectro.py` en `analisis/.cache/` y los grafica.
Es barato (milisegundos) porque no vuelve a tocar los XML: sirve para mirar el
resultado, cambiar escalas o comparar corridas sin pagar el parseo de 4 GB.

    python analisis/ver_espectros.py                  # todos los que haya
    python analisis/ver_espectros.py --metodo global  # solo un metodo
    python analisis/ver_espectros.py --log            # eje y logaritmico
"""

import argparse
import glob
import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, espectro


def disponibles(metodo=None):
    """Nombres de los espectros guardados, sin la extension."""
    patron = os.path.join(rutas.CACHE, "espectro_*.npz")
    nombres = [os.path.splitext(os.path.basename(f))[0]
               for f in sorted(glob.glob(patron))]
    if metodo:
        nombres = [n for n in nombres if n.endswith("_" + metodo)]
    return nombres


def _borde_50(x, y, suave=9, lo=0.35, hi=0.92):
    """
    Posicion del flanco al 50% de la altura del plateau.

    NO es la calibracion: es el metodo clasico, que **sobrestima** el borde
    Compton porque la lineshape de Klein-Nishina no es un escalon (tiene una
    acumulacion creciente antes del corte). Se dibuja solo como referencia
    visual de donde cae el flanco.

    La region del plateau se define por la ACUMULADA de cuentas, entre el `lo` y
    el `hi` del total, y no por posicion en el eje. Definirla por posicion no
    funciona: el Co-60 tiene un exceso enorme a bajo canal que corre el maximo,
    y el Cs-137 arranca en el umbral de trigger, asi que cualquier fraccion fija
    del eje cae en un lugar distinto en cada corrida (la primera version de esto
    daba un "plateau" de 6 cuentas en el Cs-137, que era la cola vacia).
    """
    k = np.ones(suave) / suave
    ys = np.convolve(y, k, mode="same")
    pos = np.maximum(ys, 0.0)
    total = pos.sum()
    if total <= 0:
        return np.nan, np.nan
    cum = np.cumsum(pos) / total
    region = (cum > lo) & (cum < hi)
    if region.sum() < 5:
        region = np.ones_like(cum, bool)
    # MEDIANA y no maximo: el maximo se lo lleva el pico mas alto que caiga
    # dentro de la region, que en el Co-60 es la cola del exceso de bajo canal
    # y no el plateau Compton (daba 223 cuentas donde el plateau son 120).
    alto = float(np.median(ys[region]))
    i_alto = int(np.flatnonzero(region)[-1])
    j = np.flatnonzero((ys < 0.5 * alto) & (np.arange(ys.size) > i_alto))
    return (float(x[j[0]]), alto) if j.size else (np.nan, alto)


def dibujar(nombres, log=False, salida="espectros_medidos.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    espectros = [(n, espectro.cargar(n)) for n in nombres]
    # se agrupa por corrida: cada corrida es un panel, los metodos se superponen
    corridas = []
    for n, e in espectros:
        c = str(e.get("corrida", n))
        if c not in corridas:
            corridas.append(c)

    fig, axes = plt.subplots(1, len(corridas), figsize=(7.5 * len(corridas), 5.2),
                             squeeze=False)
    for j, clave in enumerate(corridas):
        ax = axes[0][j]
        for n, e in espectros:
            if str(e.get("corrida", n)) != clave:
                continue
            met = str(e["metodo"])
            ax.errorbar(e["canales"], e["neto"], yerr=e["error"], fmt="none",
                        ecolor="0.85", lw=0.5, zorder=1)
            ax.step(e["canales"], e["neto"], where="mid", lw=1.3, zorder=2,
                    label="%s  (%d pulsos, %.0f s vivos)"
                          % (met, e["n_muestra"], e["t_vivo_muestra"]))
            if met == "global":
                b50, alto = _borde_50(e["canales"], e["neto"])
                if np.isfinite(b50):
                    ax.axvline(b50, color="C3", ls=":", lw=1.2, zorder=3,
                               label="50%% de altura: %.0f mV*ns" % b50)
                    print("  %-14s %-13s plateau %.0f cuentas, flanco al 50%% "
                          "en %.0f mV*ns" % (clave, met, alto, b50))
        e0 = dict(espectros[[str(x[1].get("corrida", x[0])) for x in espectros]
                            .index(clave)][1])
        ax.axhline(0, color="k", lw=0.8)
        if log:
            ax.set_yscale("log")
        ax.set_title("%s  —  %s\nventana %.1f + %.1f ns, %d bins de %.2f mV*ns"
                     % (clave, e0.get("fuente", ""), e0["pre"], e0["post"],
                        e0["bins"], e0["rango_max"] / e0["bins"]))
        ax.set_xlabel("area del pulso [mV*ns]")
        ax.set_ylabel("cuentas (fondo restado, normalizado por tiempo vivo)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.tight_layout()
    ruta = rutas.figura(salida)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metodo", default=None,
                   choices=["global", "convergencia"],
                   help="dibuja solo ese metodo (por defecto, todos)")
    p.add_argument("--log", action="store_true", help="eje y logaritmico")
    p.add_argument("--salida", default="espectros_medidos.png")
    a = p.parse_args(argv)

    rutas.salida_utf8()
    nombres = [n for n in disponibles(a.metodo) if "_lim" not in n]
    if not nombres:
        print("No hay espectros guardados. Corre primero:")
        print("   python analisis/hacer_espectro.py --corrida co60")
        return
    print("Espectros guardados que voy a dibujar:")
    for n in nombres:
        print("   %s" % n)
    print()
    ruta = dibujar(nombres, log=a.log, salida=a.salida)
    print()
    print("figura -> %s" % rutas.relativa(ruta))


if __name__ == "__main__":
    main()
