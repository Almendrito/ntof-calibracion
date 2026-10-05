"""
calibrar_v2.py - LA CADENA COMPLETA con el observable arreglado (Etapa 9).

    XML crudo -> area con compuerta -> normalizacion por tiempo vivo
              -> modelo extendido (mezcla de scattering + exceso)
              -> posterior por MCMC -> D y su error

Que lo separa de `calibrar_ncc.py` (Etapas 4-8)
-----------------------------------------------
Aquel script daba D = 4.32 mV*ns/keVee. Ese numero salia de un AREA mal
calculada: la linea base se promediaba sobre los primeros 40 ns del registro y
el pulso empieza a los ~25 ns, asi que la linea base se comia la subida del
propio pulso (CLAUDE.md 7.33). Ademas la resta de fondo usaba tiempo de PARED
cuando la corrida de Co-60 estuvo viva solo el 40% del tiempo (7.34).

Los dos arreglos viven en `spectrum.py` y cambian el resultado por un factor
~3.4, asi que TODO lo que dependa del numero viejo hay que rehacerlo.

Corre solo:  python calibrar_v2.py            (con MCMC, ~3 min)
             python calibrar_v2.py --rapido   (sin MCMC, ~10 s)
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/scripts/calibrar_v2.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas


import sys
import time
import numpy as np

from datos_io import spectrum
from ajuste import modelo_ext as me
from ajuste.mcmc import correr_cadenas, resumen, contornos_2d

MUESTRA = rutas.MUESTRA
FONDO = rutas.FONDO

N_CANALES, RANGO_MAX = 500, 2500.0      # 5 mV*ns por canal
VENTANA = (300.0, 1900.0)               # el borde vive en ~1200-1400
VENTANA_LIMPIA = (900.0, 1900.0)        # arriba del exceso, para sembrar
LIBRES = ["A", "D", "alpha", "w", "C", "lambda"]
E_REF = 963.4                           # borde Compton de la linea de 1173 keV

RAMAS = [
    ("Lvis (keVee)", "Lvis_keVee"),
    ("Edep (keV)", "Edep_keV"),
]


def preparar():
    """Espectro neto con el area de compuerta y el tiempo vivo."""
    r = spectrum.espectro_neto(MUESTRA, FONDO, num_canales=N_CANALES,
                               rango_max=RANGO_MAX, recortar_negativos=False)
    print(f"  tiempo vivo : muestra {r['t_vivo_muestra']:.0f} s de "
          f"{r['t_muestra']:.0f} s de pared  ({100*r['t_vivo_muestra']/r['t_muestra']:.0f}%)")
    print(f"                fondo   {r['t_vivo_fondo']:.0f} s de "
          f"{r['t_fondo']:.0f} s        -> k_fondo = {r['k_f']:.3f}")
    print(f"  descartados por calidad del pulso: {r['n_descartados'][0]} muestra, "
          f"{r['n_descartados'][1]} fondo")
    return r


def var_de(r, iv):
    m = (r["canales"] >= iv[0]) & (r["canales"] <= iv[1])
    return spectrum.varianza_predicha(r["fondo"][m], r["k_m"], r["k_f"])


def ajuste_completo(r, centros, h1, hm, w_nat):
    """Semilla de 3 parametros en la zona limpia -> modelo completo."""
    can, neto, err = r["canales"], r["neto"], r["error"]
    p0 = np.array([1.0, 1.25, 3.5, w_nat, 0.0, 250.0, 0.0, 0.0])
    z = (can > 1000) & (can < 1100)
    p0[me.I_A] = neto[z].mean() / me.modelo_ext(can, centros, h1, hm, p0)[z].mean()

    f0 = me.ajustar(can, neto, centros, h1, hm, VENTANA_LIMPIA, p0,
                    ["A", "D", "alpha"], var_fn=var_de(r, VENTANA_LIMPIA),
                    errores=err)
    p = f0["p"].copy()
    p[me.I_C], p[me.I_LAM] = 80.0, 350.0
    f = me.ajustar(can, neto, centros, h1, hm, VENTANA, p, LIBRES,
                   var_fn=var_de(r, VENTANA), errores=err)
    f["semilla"] = f0
    return f


def sistematicos(r, centros, h1, hm, f, w_nat):
    """
    Lo que el modelo viejo llamaba 'sistematico de plantilla' es ahora un
    parametro: se ve fijando w a mano y mirando cuanto se mueve D.
    """
    can, neto, err = r["canales"], r["neto"], r["error"]
    salida = {}
    for etiqueta, w_fijo in [("w = 1 (solo 1 scattering)", 1.0),
                             (f"w = {w_nat:.3f} (simulacion completa)", w_nat)]:
        p = f["p"].copy()
        p[me.I_W] = w_fijo
        g = me.ajustar(can, neto, centros, h1, hm, VENTANA, p,
                       ["A", "D", "alpha", "C", "lambda"],
                       var_fn=var_de(r, VENTANA), errores=err)
        salida[etiqueta] = g
    return salida


def barrido_ventana(r, centros, h1, hm, f, limites=(250, 300, 400, 500, 700)):
    can, neto, err = r["canales"], r["neto"], r["error"]
    out = []
    for lo in limites:
        iv = (float(lo), VENTANA[1])
        g = me.ajustar(can, neto, centros, h1, hm, iv, f["p"], LIBRES,
                       var_fn=var_de(r, iv), errores=err)
        out.append((lo, g["D"], g["chi2_red"]))
    return out


def posterior(r, centros, h1, hm, f, n_pasos=9000, n_cadenas=4):
    can = r["canales"]
    lp, expandir, n_bins = me.log_posterior(
        can, r["neto"], centros, h1, hm, VENTANA, f["p"], LIBRES,
        var_fn=var_de(r, VENTANA), errores=r["error"])
    cov = np.asarray(f["cov"], float)
    if cov.ndim != 2 or not np.all(np.isfinite(cov)):
        cov = np.diag((0.02 * np.abs(f["q"]) + 1e-9) ** 2)
    cad = correr_cadenas(lp, f["q"], cov, n_cadenas=n_cadenas, n_pasos=n_pasos,
                         dispersion=1.0, semilla=1)
    res = resumen(cad, nombres=LIBRES)
    res["n_bins"] = n_bins
    return res


def main(rapido=False):
    print("=" * 76)
    print("CALIBRACION v2 - area con compuerta + tiempo vivo + modelo extendido")
    print("=" * 76)
    r = preparar()
    can, neto, err = r["canales"], r["neto"], r["error"]

    resultados = {}
    for etiqueta, rama in RAMAS:
        print("\n" + "-" * 76)
        print(f"  PLANTILLA: {etiqueta}")
        print("-" * 76)
        centros, h1, hm, w_nat = me.cargar_componentes(rama=rama, bins=175)
        f = ajuste_completo(r, centros, h1, hm, w_nat)

        print(f"  semilla (3 par, ventana limpia): D = {f['semilla']['D']:.4f}"
              f"   chi2/ndf = {f['semilla']['chi2_red']:.3f}")
        print(f"  modelo completo, ventana {VENTANA[0]:.0f}-{VENTANA[1]:.0f}:")
        print(f"    D     = {f['D']:.4f} +- {f['err']['D']:.4f} mV*ns/keV"
              f"   ({1000/f['D']:.1f} keV por 1000 mV*ns)")
        print(f"    alpha = {f['alpha']:.3f} +- {f['err']['alpha']:.3f} sqrt(keV)"
              f"   -> FWHM({E_REF:.0f}) = {100*f['fwhm_rel_ref']:.1f}%")
        print(f"    w     = {f['w']:.3f} +- {f['err']['w']:.3f}"
              f"   (la simulacion trae {w_nat:.3f})")
        print(f"    exceso: C = {f['C']:.1f}  lambda = {f['lambda']:.0f} mV*ns"
              f"   ({100*f['C']*f['lambda']/max(np.trapezoid(neto[(can>=VENTANA[0])&(can<=VENTANA[1])], can[(can>=VENTANA[0])&(can<=VENTANA[1])]),1):.1f}% del area de la ventana)")
        print(f"    chi2/ndf = {f['chi2_red']:.3f}  ({f['ndf']} gdl)")

        print("\n  w FIJO (lo que hacia el analisis viejo con cada plantilla):")
        for k, g in sistematicos(r, centros, h1, hm, f, w_nat).items():
            print(f"    {k:34s} D = {g['D']:.4f}   chi2/ndf = {g['chi2_red']:.3f}"
                  f"   (dchi2 = {g['chi2']-f['chi2']:+.0f})")

        print("\n  estabilidad frente al limite inferior de la ventana:")
        for lo, D, c2 in barrido_ventana(r, centros, h1, hm, f):
            print(f"    lo = {lo:4d}   D = {D:.4f}   chi2/ndf = {c2:.3f}")

        if not rapido:
            t0 = time.time()
            res = posterior(r, centros, h1, hm, f)
            iD = LIBRES.index("D")
            print(f"\n  POSTERIOR ({len(res['aceptacion'])} cadenas, "
                  f"aceptacion {res['aceptacion'].mean():.2f}, {time.time()-t0:.0f} s)")
            for j, n in enumerate(LIBRES):
                print(f"    {n:7s} = {res['mediana'][j]:10.4f}"
                      f"  -{res['err_menos'][j]:.4f} +{res['err_mas'][j]:.4f}"
                      f"   R_hat = {res['R_hat'][j]:.4f}  n_eff = {res['n_eff'][j]:.0f}")
            print(f"    corr(D, w) = {res['corr'][iD, LIBRES.index('w')]:+.3f}"
                  f"   corr(D, alpha) = {res['corr'][iD, LIBRES.index('alpha')]:+.3f}")
            f["post"] = res
        f["centros"], f["h1"], f["hm"], f["w_nat"] = centros, h1, hm, w_nat
        resultados[etiqueta] = f

    # ------------------------------ resumen ------------------------------
    fl = resultados["Lvis (keVee)"]
    fe = resultados["Edep (keV)"]
    print("\n" + "=" * 76)
    print("  RESULTADO")
    print("=" * 76)
    for nom, f, u in [("D_LUZ ", fl, "mV*ns/keVee"), ("D_EDEP", fe, "mV*ns/keV")]:
        p = f.get("post")
        if p is not None:
            j = LIBRES.index("D")
            print(f"  {nom} = {p['mediana'][j]:.4f} -{p['err_menos'][j]:.4f} "
                  f"+{p['err_mas'][j]:.4f} {u}")
        else:
            print(f"  {nom} = {f['D']:.4f} +- {f['err']['D']:.4f} {u}")
    print(f"  cociente D_LUZ/D_EDEP = {fl['D']/fe['D']:.4f}"
          f"   (Birks de la simulacion: 1/0.954 = 1.048)")
    print(f"  resolucion: sigma_E = {fl['alpha']:.2f}*sqrt(E) keV"
          f"  -> {100*fl['fwhm_rel_ref']:.1f}% FWHM a {E_REF:.0f} keV")

    figura(r, resultados)
    return resultados


def figura(r, resultados):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    can, neto, err = r["canales"], r["neto"], r["error"]
    f = resultados["Lvis (keVee)"]
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))

    z = (can >= 100) & (can <= 2100)
    a = ax[0, 0]
    a.errorbar(can[z], neto[z], yerr=err[z], fmt=".", ms=2, lw=0.5, color="0.65",
               label="datos (neto)")
    for et, g in resultados.items():
        a.plot(can[z], me.modelo_ext(can[z], g["centros"], g["h1"], g["hm"], g["p"]),
               lw=1.6, label=f"{et}: D = {g['D']:.3f}")
    p = f["p"].copy()
    a.plot(can[z], p[me.I_C] * np.exp(-can[z] / p[me.I_LAM]), "k--", lw=1,
           label="exceso empirico")
    q = p.copy()
    q[me.I_C] = 0.0
    a.plot(can[z], me.modelo_ext(can[z], f["centros"], f["h1"], f["hm"], q),
           color="C3", ls=":", lw=1.2, label="solo Compton")
    a.axvspan(*VENTANA, color="green", alpha=0.06)
    a.set_xlabel("Area del pulso [mV*ns]")
    a.set_ylabel("Cuentas netas")
    a.set_title("Espectro con el area corregida")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[0, 1]
    mod = me.modelo_ext(can[z], f["centros"], f["h1"], f["hm"], f["p"])
    a.plot(can[z], (neto[z] - mod) / np.maximum(err[z], 1), ".", ms=3)
    a.axhline(0, color="k", lw=1)
    a.axhspan(-1, 1, color="green", alpha=0.12)
    a.axvspan(*VENTANA, color="green", alpha=0.06)
    a.set_ylim(-6, 6)
    a.set_xlabel("Area del pulso [mV*ns]")
    a.set_ylabel("(dato - modelo)/sigma")
    a.set_title(f"Residuos   chi2/ndf = {f['chi2_red']:.3f} ({f['ndf']} gdl)")
    a.grid(alpha=0.3)

    a = ax[1, 0]
    post = f.get("post")
    if post is not None:
        iD, iw = LIBRES.index("D"), LIBRES.index("w")
        x, y = post["muestras"][:, iD], post["muestras"][:, iw]
        H, cx, cy, niv = contornos_2d(x, y)
        a.plot(x[::20], y[::20], ".", ms=1, alpha=0.25, color="0.6")
        if niv:
            a.contour(cx, cy, H, levels=niv, colors="C0")
        a.axhline(f["w_nat"], color="C3", ls="--", lw=1,
                  label=f"w de la simulacion = {f['w_nat']:.3f}")
        a.set_xlabel("D [mV*ns/keVee]")
        a.set_ylabel("w (fraccion de scattering simple)")
        a.set_title("Posterior: D contra la mezcla de scattering")
        a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "sin MCMC (--rapido)", ha="center", va="center")
    a.grid(alpha=0.3)

    a = ax[1, 1]
    d = np.load(MUESTRA + ".areas_gate.npz")
    viejo = np.load(MUESTRA + ".areas.npz")["areas"]
    a.hist(viejo, bins=200, range=(0, 8000), histtype="step", color="C3",
           label="area vieja (linea base contaminada)")
    a.hist(d["areas"][d["bandera"] == 0], bins=200, range=(0, 8000),
           histtype="step", color="C0", label="area con compuerta")
    a.set_yscale("log")
    a.set_xlabel("Area del pulso [mV*ns]")
    a.set_ylabel("Cuentas (crudas, Co-60)")
    a.set_title("El observable, antes y despues")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("calibracion_v2.png"), dpi=130)
    print("\n  Figura guardada: figuras/calibracion_v2.png")


if __name__ == "__main__":
    main(rapido="--rapido" in sys.argv)
