"""
Comprobar la Etapa 8.

    python desde_cero/verificar_plantilla.py            (~20 min)

1. La plantilla se NIEGA a usar una simulacion de otra fuente (el .root de
   Co-60 declarado como Cs-137).

2. CONTRA EL ORACULO: con la plantilla de Geant4 (Edep, 350 bines), la D de
   `cs137_1600`, `co60` y `a_co60_1600` tiene que coincidir con la de
   `analisis/figuras/calibraciones/resultados.json` (el resultado del
   informe) dentro de la barra combinada. No es la misma cuenta al bit: el
   oraculo mide `pre` en cada archivo (2.0-2.3 ns contra 2 fijo), recorta 2
   muestras del registro en vez de 4, toma la base de los eventos SIN_PRE
   desde +16 ns en vez de +84, y evalua la plantilla por interpolacion en vez
   de integrarla por bin.

3. EL CAMBIO: D con Klein-Nishina contra D con Geant4, sobre las mismas
   corridas. §7.18 predice del orden de +4 %.

4. CIERRE con la plantilla de Geant4: 8 realizaciones Poisson del ajuste del
   Cs-137, por la cadena entera (NCC con la plantilla de Geant4 -> recocido ->
   barra). Sesgo compatible con cero.

5. INFORMATIVO, sin criterio: que fraccion de scattering simple w pide el
   dato, en la ventana de la plantilla completa. Es el hallazgo de esta
   etapa y no tiene "respuesta correcta" contra la cual comparar.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import ajuste, espectro, lector, modelo, ncc, plantilla  # noqa: E402

FUENTE = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
ORACULO = os.path.join(lector.RAIZ, "analisis", "figuras", "calibraciones",
                       "resultados.json")


def main():
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-38s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    print("=" * 78)
    print("1. Rechazo de una simulacion de otra fuente")
    print("=" * 78)
    try:
        plantilla.Plantilla("Cs-137", ruta=plantilla.SIMULACIONES["Co-60"])
        check("Co-60.root declarado como Cs-137", False, "NO se nego")
    except ValueError as e:
        check("Co-60.root declarado como Cs-137", True, "se niega: %s..." % str(e)[:60])

    geant = {f: plantilla.Plantilla(f) for f in ("Cs-137", "Co-60")}
    ora = json.load(open(ORACULO, encoding="utf-8"))
    print()
    print("=" * 78)
    print("2 y 3. Geant4 contra el oraculo, y Klein-Nishina contra Geant4")
    print("=" * 78)
    reales = {}
    for clave, f in FUENTE.items():
        if not os.path.exists(lector.CORRIDAS_TRABAJO[clave][0]):
            continue
        e, _, _ = espectro.de_corrida(clave)
        g4 = ajuste.calibrar(e["canales"], e["neto"], e["error"], f, plantilla=geant[f])
        kn = ajuste.calibrar(e["canales"], e["neto"], e["error"], f)
        reales[clave] = (e, g4)
        o = ora[clave]
        err = np.hypot(g4["err_D"], o["err_D"])
        check("%s: Geant4 contra el oraculo" % clave, abs(g4["D"] - o["D"]) < err,
              "%.4f +- %.4f / %.4f +- %.4f (%+.2f %%, %.1f barras) ventana "
              "[%.0f, %.0f] / [%.0f, %.0f]"
              % (g4["D"], g4["err_D"], o["D"], o["err_D"], 100 * (g4["D"] / o["D"] - 1),
                 abs(g4["D"] - o["D"]) / err, g4["ventana"]["lo"], g4["ventana"]["hi"],
                 *o["ventana"]))
        print("  %-38s Klein-Nishina %.4f / Geant4 %.4f - 1 = %+.2f %%   "
              "(chi2/ndf %.2f / %.2f, ventanas distintas)"
              % ("", kn["D"], g4["D"], 100 * (kn["D"] / g4["D"] - 1),
                 kn["chi2"] / kn["ndf"], g4["chi2"] / g4["ndf"]))

    print()
    print("=" * 78)
    print("4. Cierre con la plantilla de Geant4 (Cs-137, 8 realizaciones)")
    print("=" * 78)
    if "cs137_1600" in reales:
        e, r = reales["cs137_1600"]
        x = e["canales"]
        mu = modelo.modelo(x, r["A"], r["D"], r["sigma_rel"], "Cs-137",
                           plantilla=geant["Cs-137"])
        rng = np.random.default_rng(21)
        dif, errs = [], []
        for k in range(8):
            y = rng.poisson(mu).astype(float)
            q = ajuste.calibrar(x, y, np.sqrt(np.maximum(y, 1.0)), "Cs-137",
                                plantilla=geant["Cs-137"], semilla=200 + k)
            dif.append(q["D"] - r["D"])
            errs.append(q["err_D"])
        dif = np.array(dif)
        s, es = dif.mean(), dif.std(ddof=1) / np.sqrt(dif.size)
        check("sesgo de D", abs(s) < 3 * es + 1e-12,
              "%+.5f +- %.5f (%+.3f %%); dispersion %.5f contra barra media %.5f"
              % (s, es, 100 * s / r["D"], dif.std(ddof=1), np.mean(errs)))

    print()
    print("=" * 78)
    print("5. Que fraccion de scattering simple pide el dato (informativo)")
    print("=" * 78)
    # Con la ventana del NCC (desde 0.45 B) y otra vez desde 0.80 B, lejos del
    # exceso de bajo canal que el modelo no tiene. Si la preferencia por w
    # desaparece al alejarse del exceso, era el exceso y no el scattering.
    for clave, (e, g4) in reales.items():
        f = FUENTE[clave]
        w_sim = float(np.mean(geant[f].n_compton == 1))
        for ini in (None, 0.80):
            v = g4["ventana"] if ini is None else dict(g4["ventana"],
                                                       lo=ini * g4["ventana"]["B"])
            filas, ndf = [], 0
            for w in (w_sim, 0.85, 1.0):
                q = ajuste.calibrar(e["canales"], e["neto"], e["error"], f,
                                    ventana=v, plantilla=plantilla.Mezcla(f, w),
                                    n_arranques=3)
                ndf = q["ndf"]
                filas.append("w = %.3f: D %.4f, chi2 %.0f" % (w, q["D"], q["chi2"]))
            print("  %-12s desde %s (ndf %d) %s"
                  % (clave, "NCC " if ini is None else "%.2fB" % ini, ndf,
                     " | ".join(filas)))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
