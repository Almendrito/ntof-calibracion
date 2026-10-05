"""
simulacion.py — leer los .root de Geant4 y armar la plantilla de energia.

La plantilla es un histograma de energia depositada (o de luz visible) que ya
trae metida la fisica que la formula analitica de Klein-Nishina no tiene:
scattering multiple, escape, y el quenching de Birks. Se usa como S2 en el NCC.

El TTree `ComptonSpectrum` de los dos archivos tiene las mismas ramas:

    Edep_keV     energia depositada en el centellador
    Lvis_keVee   luz visible, con Birks ya aplicado  [keVee]
    E0_keV       energia del gamma primario -> dice de QUE fuente es el archivo
    NCompton     cuantos scatterings Compton hubo (==1 es scattering simple)
    NPhoto, FullAbs, QF

Correr `python analisis/simulacion.py` compara las dos plantillas con los bordes
teoricos y deja `figuras/04_simulacion.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, fisica

ARBOL = "ComptonSpectrum"


def lineas_primarias(ruta, arbol=ARBOL, tol=1.0):
    """
    Que gammas simula el archivo, leyendo `E0_keV`. Devuelve [(energia, cuentas)].

    Sirve para verificar que la plantilla sea de la MISMA fuente que los datos:
    usar el .root de Co-60 contra un espectro de Cs-137 no da error en ningun
    lado, da un numero equivocado bien formateado.
    """
    import uproot
    e0 = uproot.open(ruta)[arbol]["E0_keV"].array(library="np")
    vals, cuentas = np.unique(np.round(e0 / tol) * tol, return_counts=True)
    return [(float(v), int(c)) for v, c in zip(vals, cuentas)]


def coincide_fuente(ruta, fuente, tol=2.0):
    """True si las lineas del .root son las de `fuente` (fisica.FUENTES)."""
    sim = [e for e, _ in lineas_primarias(ruta)]
    esperadas = fisica.FUENTES[fuente]
    if len(sim) != len(esperadas):
        return False
    return all(abs(a - b) < tol for a, b in zip(sorted(sim), sorted(esperadas)))


def cargar(ruta, rama="Lvis_keVee", bins=350, rango=None, solo_simple=False,
           arbol=ARBOL):
    """
    Devuelve `(centros, pesos)`: la plantilla en ENERGIA, ya histogrameada.

    rama : "Lvis_keVee" (luz, escala keVee) o "Edep_keV" (energia depositada).
        Las dos difieren por el quenching de Birks; hay que decir cual se usa
        porque la D que sale despues esta en las unidades de la que se eligio.
    solo_simple : quedarse con `NCompton == 1`. Eso es exactamente lo que
        describe Klein-Nishina, asi que sirve para medir cuanto aporta el
        scattering multiple comparando las dos plantillas.
    rango : (min, max) en keV. Por defecto de 0 al maximo de la rama.
    """
    import uproot
    t = uproot.open(ruta)[arbol]
    ramas = [rama] + (["NCompton"] if solo_simple else [])
    d = t.arrays(ramas, library="np")
    e = d[rama]
    if solo_simple:
        e = e[d["NCompton"] == 1]
    e = e[e > 0]
    if rango is None:
        rango = (0.0, float(e.max()))
    pesos, bordes = np.histogram(e, bins=bins, range=rango)
    return 0.5 * (bordes[:-1] + bordes[1:]), pesos.astype(float)


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES
              if rutas.corrida(k).sim and os.path.exists(rutas.corrida(k).sim)]
    fig, axes = plt.subplots(1, len(claves), figsize=(7.5 * len(claves), 5),
                             squeeze=False)

    for j, clave in enumerate(claves):
        c = rutas.corrida(clave)
        print("=" * 72)
        print("%s  -  %s" % (clave, os.path.basename(c.sim)))
        print("=" * 72)

        print("  gammas primarios del .root:")
        for e0, n in lineas_primarias(c.sim):
            print("     %9.1f keV   %8d eventos   -> borde Compton %7.1f keV"
                  % (e0, n, fisica.borde_compton(e0)))
        ok = coincide_fuente(c.sim, c.fuente)
        print("  fuente medida: %s  ->  %s"
              % (c.fuente, "coincide" if ok else "NO COINCIDE, no usar"))

        cen, tot = cargar(c.sim, "Lvis_keVee")
        _, simple = cargar(c.sim, "Lvis_keVee", solo_simple=True,
                           rango=(cen[0] - 0.5 * (cen[1] - cen[0]),
                                  cen[-1] + 0.5 * (cen[1] - cen[0])))
        _, edep = cargar(c.sim, "Edep_keV",
                         rango=(cen[0] - 0.5 * (cen[1] - cen[0]),
                                cen[-1] + 0.5 * (cen[1] - cen[0])))
        f_simple = simple.sum() / tot.sum()
        print("  %d eventos con energia > 0; %.1f%% son scattering SIMPLE"
              % (tot.sum(), 100 * f_simple))
        print("  bordes teoricos: %s keV"
              % ", ".join("%.1f" % b for b in fisica.bordes(c.fuente)))

        ax = axes[0][j]
        ax.step(cen, tot / tot.max(), where="mid", lw=1.2,
                label="Lvis (luz, keVee)")
        ax.step(cen, edep / edep.max(), where="mid", lw=1.0, alpha=0.8,
                label="Edep (energia depositada, keV)")
        ax.step(cen, simple / tot.max(), where="mid", lw=1.0,
                label="Lvis, solo 1 scattering (%.0f%%)" % (100 * f_simple))
        for b in fisica.bordes(c.fuente):
            ax.axvline(b, color="C3", ls="--", lw=1)
            ax.annotate("%.0f keV" % b, xy=(b, 1.0), xytext=(b - 8, 1.03),
                        rotation=90, fontsize=8, color="C3", ha="right")
        ax.set_title("%s  -  %s" % (clave, c.fuente))
        ax.set_xlabel("energia [keV o keVee]")
        ax.set_ylabel("cuentas (normalizado)")
        ax.set_ylim(0, 1.18)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        print()

    fig.tight_layout()
    salida = rutas.figura("04_simulacion.png")
    fig.savefig(salida, dpi=130)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
