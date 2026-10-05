"""
Comprobar la Etapa 7.

    python desde_cero/verificar_ajuste.py               (~15 min)
    python desde_cero/verificar_ajuste.py --rapido      (8 realizaciones en vez de 30)

1. DATOS REALES (cs137_1600 y co60, plantilla Klein-Nishina):
   a. los 8 arranques, pulidos, caen en el MISMO minimo (dispersion de D
      menor que un decimo de su barra);
   b. la barra por curvatura coincide con el chi2 PERFILADO: donde el chi2,
      re-optimizando sigma y A en cada D, sube 1 (§7.49: cuando no coinciden,
      gana el perfil y la diferencia es informacion).

2. CIERRE CON VERDAD CONOCIDA. Se toma el ajuste real del Cs-137 como verdad,
   se fabrica el espectro que predice (en todo el rango, no solo la ventana:
   el NCC tiene que encontrar la ventana solo) y se muestrea 30 veces con
   Poisson. Cada realizacion pasa por la cadena ENTERA: NCC -> 8 recocidos ->
   pulido -> barra. Se mide el sesgo de D (con su error), la dispersion de
   las D contra la barra que el ajuste dice, y la COBERTURA: en cuantas la
   verdad cae dentro de +-1 barra. Si la barra esta bien, ~68 %.

3. CONTRA replica_paper/metodo.py sobre NUESTROS espectros (misma compuerta de
   12 ns; los numeros de replica_paper/informe.txt son con la compuerta que
   mide cada archivo y no se pueden comparar). NO es la misma cuenta: el
   oraculo minimiza la Ec. 14 sin pesos, mete A en el recocido, elige el
   sigma de la plantilla del NCC por barrido, usa el barrido libre, y en el
   Co-60 tiene sigma1 y sigma2. Se reporta la diferencia en unidades de
   nuestra barra; no se le exige nada.
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from replica_paper import metodo as oraculo                 # noqa: E402
from desde_cero import ajuste, espectro, lector, modelo, ncc  # noqa: E402

FUENTE = {"cs137_1600": "Cs-137", "co60": "Co-60"}
LINEAS = {"Cs-137": (661.657,), "Co-60": (1173.228, 1332.492)}


def cruce_perfil(res):
    """Semiancho donde el chi2 perfilado sube 1, promediando los dos lados."""
    Ds, dchi = ajuste.perfil_D(res, n=41, ancho=2.0)
    D0 = res["D"]
    izq = Ds <= D0
    der = Ds >= D0
    a = np.interp(1.0, dchi[izq][::-1], D0 - Ds[izq][::-1])
    b = np.interp(1.0, dchi[der], Ds[der] - D0)
    return 0.5 * (a + b)


def main():
    rapido = "--rapido" in sys.argv
    n_real = 8 if rapido else 30
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-40s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    print("=" * 78)
    print("1. Datos reales, plantilla Klein-Nishina")
    print("=" * 78)
    reales = {}
    for clave, f in FUENTE.items():
        if not os.path.exists(lector.CORRIDAS_TRABAJO[clave][0]):
            continue
        e, _, _ = espectro.de_corrida(clave)
        r = ajuste.calibrar(e["canales"], e["neto"], e["error"], f)
        reales[clave] = (e, r)
        print("  %-12s D = %.4f +- %.4f (sin escalar %.4f)  sigma_rel = %.4f  "
              "chi2/ndf = %.2f  ventana [%.0f, %.0f]"
              % (clave, r["D"], r["err_D"], r["err_D_crudo"], r["sigma_rel"],
                 r["chi2"] / r["ndf"], r["ventana"]["lo"], r["ventana"]["hi"]))
        disp = np.ptp(r["D_arranques"])
        check("   a. 8 arranques pulidos", disp < 0.1 * r["err_D_crudo"],
              "rango de D %.1e (barra %.1e); sin pulir %.1e"
              % (disp, r["err_D_crudo"], np.ptp(r["D_recocido"])))
        perf = cruce_perfil(r)
        check("   b. perfil chi2 + 1 contra curvatura",
              abs(perf / r["err_D_crudo"] - 1) < 0.15,
              "%.5f contra %.5f" % (perf, r["err_D_crudo"]))

    print()
    print("=" * 78)
    print("2. Cierre con verdad conocida (Cs-137, %d realizaciones)" % n_real)
    print("=" * 78)
    if "cs137_1600" in reales:
        e, r = reales["cs137_1600"]
        x = e["canales"]
        verdad = dict(D=r["D"], s=r["sigma_rel"], A=r["A"])
        mu = modelo.modelo(x, r["A"], r["D"], r["sigma_rel"], "Cs-137")
        rng = np.random.default_rng(11)
        Ds, errs, t0 = [], [], time.time()
        for k in range(n_real):
            y = rng.poisson(mu).astype(float)
            q = ajuste.calibrar(x, y, np.sqrt(np.maximum(y, 1.0)), "Cs-137",
                                semilla=100 + k)
            Ds.append(q["D"])
            errs.append(q["err_D"])
        Ds, errs = np.array(Ds), np.array(errs)
        dif = Ds - verdad["D"]
        sesgo, err_sesgo = dif.mean(), dif.std(ddof=1) / np.sqrt(n_real)
        cob = np.mean(np.abs(dif) <= errs)
        print("  verdad D = %.4f, sigma_rel = %.4f, A = %.0f cuentas; ~%.0f cuentas "
              "por bin en la ventana; %.0f s"
              % (verdad["D"], verdad["s"], verdad["A"],
                 np.median(mu[(x >= r["ventana"]["lo"]) & (x <= r["ventana"]["hi"])]),
                 time.time() - t0))
        check("sesgo de D", abs(sesgo) < 3 * err_sesgo,
              "%+.5f +- %.5f (%+.3f %%)" % (sesgo, err_sesgo, 100 * sesgo / verdad["D"]))
        check("dispersion de D / barra media", 0.7 < dif.std(ddof=1) / errs.mean() < 1.4,
              "%.5f / %.5f = %.2f" % (dif.std(ddof=1), errs.mean(),
                                      dif.std(ddof=1) / errs.mean()))
        # con n realizaciones la cobertura binomial fluctua sqrt(.68*.32/n)
        tol = 2.5 * np.sqrt(0.68 * 0.32 / n_real)
        check("cobertura de +-1 barra", abs(cob - 0.68) < tol,
              "%d/%d = %.0f %% (esperado 68 +- %.0f %%)"
              % (int(round(cob * n_real)), n_real, 100 * cob, 100 * tol / 2.5))

    print()
    print("=" * 78)
    print("3. Contra replica_paper/metodo.py (sobre nuestros espectros)")
    print("=" * 78)
    for clave, (e, r) in reales.items():
        f = FUENTE[clave]
        kw = dict(ligar_origen=True) if f == "Co-60" else {}
        t0 = time.time()
        o = oraculo.ajustar(e["canales"], e["neto"], LINEAS[f], **kw)
        print("  %-12s nuestro D = %.4f +- %.4f   replica_paper D = %.4f "
              "(%+.2f %%, %+.1f barras; ventana [%.0f, %.0f], chi2/ndf segun "
              "ella %.2f; %.0f s)"
              % (clave, r["D"], r["err_D"], o["D"], 100 * (o["D"] / r["D"] - 1),
                 (o["D"] - r["D"]) / r["err_D"], o["intervalo"][0],
                 o["intervalo"][1], o["chi2_ndf"], time.time() - t0))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
