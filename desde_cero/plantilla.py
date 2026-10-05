"""
Etapa 8 — de Klein-Nishina a Geant4: la plantilla con la fisica del detector.

    python desde_cero/plantilla.py

Que resuelve esta etapa
-----------------------
La plantilla del paper es Klein-Nishina: UN scattering Compton sobre un
electron libre, y el gamma se va. En un centellador de 5 cm no es asi. El
gamma dispersado puede volver a interactuar adentro (scattering MULTIPLE) y
depositar mas energia; un cuarto de los eventos del Co-60 simulado son asi
(§7.17). Esos eventos caen por ENCIMA del borde de un solo scattering y
rellenan el hueco que Klein-Nishina deja entre el borde y cero. Un modelo
que no los tiene explica esas cuentas de mas corriendo el borde hacia arriba
en canal: D sale grande (+4 % en §7.18). Las Etapas 6 y 7 lo mostraron en
los residuos: la Klein-Nishina ajusta con chi2/ndf 1.4-2.1 y residuos con
estructura.

La salida es no escribir la forma a mano sino SIMULARLA: Geant4 tira gammas
de la fuente contra un cilindro de plastico con su carcasa y cuenta cuanta
energia queda adentro evento por evento. Eso trae el scattering multiple,
los electrones que se escapan, la carcasa de aluminio, todo lo que la
geometria de la simulacion tenga (ver `simulacion/`, CLAUDE.md §7.69).

Las decisiones, y por que
-------------------------
1. **Rama `Edep_keV`, energia DEPOSITADA** (la del informe). La otra,
   `Lvis_keVee`, es luz con el quenching de Birks aplicado: da una D en
   keVee, 4.5 % mas grande (§7.17). Las dos son validas; hay que decir cual.

2. **Se verifica que la simulacion sea de la fuente.** Un .root de Co-60
   contra un espectro de Cs-137 no da error en ningun lado: da una D
   perfectamente formateada y equivocada (§7.41). Se leen los gammas
   primarios (`E0_keV`) y si no son las lineas de la fuente, se levanta un
   error en vez de seguir.

3. **Histograma de 350 bines, INTEGRADO por bin del espectro.** La
   simulacion llega como eventos; se histogramean en energia y la fraccion de
   cuentas en un intervalo [E1, E2] se saca de la ACUMULADA (dentro de cada
   bin de 350 se supone uniforme). Es el mismo remedio de la Etapa 7:
   evaluar la plantilla en el centro de cada bin del espectro haria el
   modelo discontinuo en D (y con un histograma encima, escalonado: el
   "peine" de §7.25). Con la acumulada el modelo es continuo. 350 bines de
   ~3-4 keV son 30-40 veces mas finos que la resolucion (sigma ~ 60-120
   keV), asi que el numero de bines no mueve D (§7.26, §7.69).

4. **Nada mas cambia.** La plantilla entra por `modelo.modelo(...,
   plantilla=)` y el mapa, la resolucion, sigma_rel (referido al borde
   Compton teorico mas bajo), el NCC y el ajuste son los de las Etapas 5-7.
   Asi la diferencia de D entre las dos plantillas es SOLO la plantilla.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import fisica                               # noqa: E402
from desde_cero.lector import RAIZ                          # noqa: E402

SIMULACIONES = {
    "Cs-137": os.path.join(RAIZ, "datos", "ComptonSpectrum_Cs137_10M.root"),
    "Co-60": os.path.join(RAIZ, "datos", "ComptonSpectrum_Co60_3M.root"),
}
ARBOL = "ComptonSpectrum"
RAMA = "Edep_keV"
BINES = 350
TOL_LINEA = 2.0           # [keV] para reconocer una linea primaria


def lineas_primarias(ruta):
    """[(energia, eventos)] de los gammas que simula el archivo (rama E0_keV)."""
    import uproot
    e0 = uproot.open(ruta)[ARBOL]["E0_keV"].array(library="np")
    v, c = np.unique(np.round(e0, 1), return_counts=True)
    return [(float(a), int(b)) for a, b in zip(v, c)]


class Plantilla:
    """
    La simulacion como plantilla de `modelo.modelo`.

    `fraccion(E1, E2)`: fraccion de los eventos (con deposito > 0) entre E1 y
    E2 keV. `densidad(E)`: la misma en 1/keV, en el centro de los bines.
    `e_tope`: la energia mas alta del histograma.
    """

    def __init__(self, fuente, rama=RAMA, bines=BINES, ruta=None,
                 n_compton=None):
        self.fuente, self.rama = fuente, rama
        self.ruta = ruta or SIMULACIONES[fuente]
        esperadas = sorted(e for e, _ in fisica.FUENTES[fuente])
        sim = sorted(e for e, _ in lineas_primarias(self.ruta))
        if len(sim) != len(esperadas) or any(abs(a - b) > TOL_LINEA
                                             for a, b in zip(sim, esperadas)):
            raise ValueError(
                "%s simula gammas de %s keV y la fuente %s tiene %s keV: esta "
                "plantilla no es de esta fuente" % (os.path.basename(self.ruta),
                                                     sim, fuente, esperadas))
        self.lineas_simuladas = sim

        import uproot
        d = uproot.open(self.ruta)[ARBOL].arrays([rama, "NCompton"], library="np")
        e, nc = d[rama], d["NCompton"]
        ok = e > 0      # el .root ya guarda solo eventos con deposito; por las dudas
        e, nc = e[ok], nc[ok]
        if n_compton is not None:
            e = e[np.isin(nc, np.atleast_1d(n_compton))]
        self.energias, self.n_compton = e, nc
        self.bordes = np.linspace(0.0, float(d[rama][ok].max()), bines + 1)
        h, _ = np.histogram(e, bins=self.bordes)
        self.cuentas = h.astype(float)
        self.acum = np.concatenate(([0.0], np.cumsum(h))) / h.sum()
        self.e_tope = float(self.bordes[-1])
        ancho = np.diff(self.bordes)
        self.centros = 0.5 * (self.bordes[1:] + self.bordes[:-1])
        self._dens = self.cuentas / (h.sum() * ancho)

    def fraccion(self, E_inf, E_sup):
        # acumulada lineal dentro de cada bin = densidad uniforme en el bin
        return (np.interp(E_sup, self.bordes, self.acum, left=0.0, right=1.0)
                - np.interp(E_inf, self.bordes, self.acum, left=0.0, right=1.0))

    def densidad(self, E):
        return np.interp(E, self.centros, self._dens, left=0.0, right=0.0)


class Mezcla:
    """
    w * (eventos de UN scattering) + (1 - w) * (el resto), cada parte con su
    propia forma normalizada. Con w = fraccion simulada de NCompton == 1 es la
    plantilla completa; con w = 1, la de un solo scattering.

    Existe para MEDIR cuanto scattering multiple pide el dato. Resultado
    (27-sep-2026, compuerta de 12 ns, sin componente de exceso; ver PLAN.md):
    - con la ventana del NCC (desde 0.45 B) el dato "pide" w ~ 0.9-1.1 contra
      0.73-0.75 de la simulacion (Delta chi2 de 90 a 250). NO es el
      scattering: es el exceso de bajo canal, que el modelo no tiene;
    - arrancando la ventana en 0.80 B, lejos del exceso, el chi2 ya no
      distingue w entre 0.75 y 1 (diferencias <= 1 en 3 de 3 corridas) pero D
      se mueve 4-6 % en ese rango: el dato no fija w, lo fija la simulacion,
      y lo que la geometria simulada tenga mal entra derecho en D.
    """

    def __init__(self, fuente, w, rama=RAMA, bines=BINES):
        self.uno = Plantilla(fuente, rama, bines, n_compton=1)
        self.resto = Plantilla(fuente, rama, bines, n_compton=[0] + list(range(2, 100)))
        self.w = float(w)
        self.e_tope = max(self.uno.e_tope, self.resto.e_tope)

    def fraccion(self, E_inf, E_sup):
        return (self.w * self.uno.fraccion(E_inf, E_sup)
                + (1 - self.w) * self.resto.fraccion(E_inf, E_sup))

    def densidad(self, E):
        return self.w * self.uno.densidad(E) + (1 - self.w) * self.resto.densidad(E)


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def comparar(clave, fuente, e, plantillas):
    """Ajusta el mismo espectro con cada plantilla. Devuelve {nombre: res}."""
    from desde_cero import ajuste
    return {n: ajuste.calibrar(e["canales"], e["neto"], e["error"], fuente,
                               plantilla=p) for n, p in plantillas.items()}


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero import espectro, modelo, ncc
    from desde_cero.lector import CORRIDAS_TRABAJO

    fuentes = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
    geant = {f: Plantilla(f) for f in ("Cs-137", "Co-60")}

    fig, ejes = plt.subplots(3, 3, figsize=(18, 13))
    # fila 1: las plantillas en energia, y la separacion por dispersiones
    for col, f in enumerate(("Cs-137", "Co-60")):
        p = geant[f]
        ax = ejes[0, col]
        E = np.linspace(0, 1.15 * p.e_tope, 3000)
        paso = E[1] - E[0]
        kn = modelo.cuentas_por_bin(E - paso / 2, E + paso / 2,
                                    *modelo.lineas_de(f)) / paso
        g = p.fraccion(E - paso / 2, E + paso / 2) / paso
        ax.plot(E, kn, lw=1.2, label="Klein-Nishina (1 scattering, sin detector)")
        ax.plot(E, g, lw=1.2, label="Geant4, Edep")
        # por numero de dispersiones: apilado
        h_tot = np.zeros(len(p.centros))
        nc = p.n_compton
        for nombre, sel in (("NCompton = 1", nc == 1), ("NCompton = 2", nc == 2),
                            ("NCompton >= 3", nc >= 3), ("NCompton = 0", nc == 0)):
            h, _ = np.histogram(p.energias[sel], bins=p.bordes)
            d = h / (p.cuentas.sum() * np.diff(p.bordes))
            ax.fill_between(p.centros, h_tot, h_tot + d, step="mid", alpha=0.35,
                            label="%s (%.1f %%)" % (nombre, 100 * sel.mean()))
            h_tot = h_tot + d
        for b in fisica.e_max(np.array([x for x, _ in fisica.FUENTES[f]])):
            ax.axvline(b, color="k", ls=":", lw=0.8)
        ax.set(xlabel="energia depositada [keV]", ylabel="densidad [1/keV]",
               title="%s: %d eventos con deposito (%s)"
               % (f, p.energias.size, os.path.basename(p.ruta)),
               xlim=(0, 1.1 * p.e_tope))
        ax.legend(fontsize=7)
    ejes[0, 2].axis("off")

    print("  %-12s %-14s %9s %8s %9s %8s   %s"
          % ("corrida", "plantilla", "D", "+-", "sigma_rel", "chi2/ndf", "ventana"))
    for col, clave in enumerate(fuentes):
        if not os.path.exists(CORRIDAS_TRABAJO[clave][0]):
            continue
        f = fuentes[clave]
        e, _, _ = espectro.de_corrida(clave)
        res = comparar(clave, f, e, {"Klein-Nishina": None, "Geant4": geant[f]})
        for n, r in res.items():
            print("  %-12s %-14s %9.4f %8.4f %9.4f %8.2f   [%.0f, %.0f]"
                  % (clave, n, r["D"], r["err_D"], r["sigma_rel"],
                     r["chi2"] / r["ndf"], r["ventana"]["lo"], r["ventana"]["hi"]))
        kn, g4 = res["Klein-Nishina"], res["Geant4"]
        print("  %-12s Klein-Nishina / Geant4 - 1 = %+.2f %%"
              % ("", 100 * (kn["D"] / g4["D"] - 1)))

        ax = ejes[1, col]
        x = e["canales"]
        m = (x > 0.3 * g4["ventana"]["lo"]) & (x < 1.4 * g4["ventana"]["hi"])
        ax.errorbar(x[m], e["neto"][m], e["error"][m], fmt=".", ms=2, lw=0.4,
                    color="0.55", label="neto")
        for n, r, c in (("Klein-Nishina", kn, "tab:blue"), ("Geant4", g4, "tab:red")):
            ax.plot(r["problema"].x, r["problema"].prediccion(r["p"]), color=c,
                    lw=1.4, label="%s: D = %.4f, chi2/ndf = %.2f"
                    % (n, r["D"], r["chi2"] / r["ndf"]))
        ax.set(title=clave, xlabel="area [mV*ns]", ylabel="cuentas")
        ax.legend(fontsize=7)

        ax = ejes[2, col]
        for n, r, c in (("Klein-Nishina", kn, "tab:blue"), ("Geant4", g4, "tab:red")):
            P = r["problema"]
            ax.plot(P.x, (P.y - P.prediccion(r["p"])) / P.e, ".", ms=3, color=c,
                    alpha=0.7, label=n)
        ax.axhline(0, color="k", lw=0.8)
        ax.set(xlabel="area [mV*ns]", ylabel="(dato - modelo) / error",
               title="residuos")
        ax.legend(fontsize=7)
    for ax in ejes.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa8_plantilla.png")
    fig.savefig(salida, dpi=100)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
