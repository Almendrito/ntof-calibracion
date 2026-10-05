"""
Etapa 9 — todas las corridas, y las dos pruebas del informe.

    python desde_cero/lote.py              calibra lo que falte (~40 min la 1a vez)
    python desde_cero/lote.py --rehacer    recalibra todo
    python desde_cero/lote.py --tests      solo los tests, sobre lo ya calibrado

Que resuelve esta etapa
-----------------------
Pasa las 28 corridas del catalogo por la cadena de desde_cero (Etapas 1-8:
compuerta comun de 12 ns, espectro neto, NCC, recocido, plantilla de Geant4
en energia depositada) y con las D hace las dos pruebas:

TEST B — la ganancia contra el alto voltaje. Un PMT multiplica los
electrones en cada dinodo por k ~ V^a, y con n dinodos la ganancia total va
como V^(a n) = V^b. Se ajusta D(V) de cada detector y cada fuente con dos
hipotesis, una RECTA y una LEY DE POTENCIA, y se comparan los chi2. La ley
de potencia se ajusta en log-log, donde es una recta:

    ln D = ln D0 + b ln(V / V0)

con V0 la media GEOMETRICA de los voltajes. Por que ahi: el error de la
ordenada de una recta crece con la distancia del pivote al centro de los
datos, y la ordenada y la pendiente salen correlacionadas. Centrando el
pivote, ln D0 y b salen independientes y el error de b es solo el de b.
Los errores de D pasan a log como err/D (d ln D = dD / D).

TEST A — el eje: las dos fuentes al mismo punto de operacion (mismo detector
y mismo HV) tienen que dar la misma D si el area es proporcional a la
energia. Pero en escala de energia DEPOSITADA no tiene que dar 1: el area es
proporcional a la LUZ, y la luz por keV (quenching de Birks) no es igual a
477 keV que a 963 keV. Lo que se espera es

    D(Cs)/D(Co) = q(Cs)/q(Co) ~ 0.980,   q = <Lvis/Edep> de la simulacion

y el desacuerdo es (D(Cs)/D(Co)) / 0.980 - 1. Con dos fuentes hay ademas una
recta canal = D E + c0 que pasa exacta por los dos bordes; su ordenada c0 es
lo que diria un pedestal, y lo unico que dice algo es su ERROR (dos puntos,
dos parametros: siempre pasa). La D conjunta es la media de las dos pesada
por 1/sigma^2.

Propagacion de errores (informe, apendice ap:propagacion): para una funcion
f de variables independientes, sigma_f^2 = sum (df/dx_i)^2 sigma_i^2. Asi
    cociente r = D1/D2:      sigma_r / r = sqrt((s1/D1)^2 + (s2/D2)^2)
    media pesada:            sigma = 1 / sqrt(sum 1/s_i^2)
    c0 = E1 E2 (D1 - D2) / (E2 - E1):  sigma = E1 E2 / (E2 - E1) * sqrt(s1^2 + s2^2)

Del catalogo `analisis/rutas.py` se usan SOLO las rutas y los metadatos
(detector, HV, fuente, descartada): ningun calculo.
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import ajuste, area, espectro, fisica, lector, plantilla  # noqa: E402

SALIDA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "resultados_lote.json")
FRAC_QUENCH = 0.7        # q = <Lvis/Edep> para Edep > 0.7 * E_ref


def catalogo():
    """Las corridas: rutas y metadatos, sin calculo."""
    from analisis import rutas
    return [dict(clave=k, muestra=c.muestra, fondo=c.fondo, fuente=c.fuente,
                 detector=c.detector, hv=c.hv, descartada=rutas.descartada(k))
            for k, c in rutas.CORRIDAS.items()]


# --------------------------------------------------------------------------
# Calibrar
# --------------------------------------------------------------------------
def calibrar_corrida(c, plantillas):
    t0 = time.time()
    d_m = area.areas_de_archivo(c["muestra"])
    d_f = area.areas_de_archivo(c["fondo"]) if c["fondo"] else None
    e = espectro.espectro_neto(d_m, d_f)
    r = ajuste.calibrar(e["canales"], e["neto"], e["error"], c["fuente"],
                        plantilla=plantillas[c["fuente"]])
    w = r["ventana"]
    return dict(c, D=r["D"], err_D=r["err_D"], err_D_crudo=r["err_D_crudo"],
                sigma_rel=r["sigma_rel"], chi2=r["chi2"], ndf=r["ndf"],
                corr=r["corr"], ventana=[w["lo"], w["hi"]], D0=w["D0"],
                rho=w["rho"], hay_fondo=e["hay_fondo"], frac_viva=e["frac_viva_m"],
                pulsos=int(e["n_m"]), fecha=d_m["reloj"][0][:10],
                dispersion_arranques=float(np.ptp(r["D_arranques"])),
                segundos=time.time() - t0)


def correr(rehacer=False):
    res = {}
    if os.path.exists(SALIDA) and not rehacer:
        res = json.load(open(SALIDA, encoding="utf-8"))
    plantillas = {f: plantilla.Plantilla(f) for f in ("Cs-137", "Co-60")}
    cat = catalogo()
    for i, c in enumerate(cat):
        if c["clave"] in res:
            continue
        if not os.path.exists(c["muestra"]):
            print("  [%2d/%d] %-18s no esta en disco" % (i + 1, len(cat), c["clave"]))
            continue
        r = calibrar_corrida(c, plantillas)
        res[c["clave"]] = r
        print("  [%2d/%d] %-18s D = %.4f +- %.4f  chi2/ndf %.2f  (%.0f s)"
              % (i + 1, len(cat), c["clave"], r["D"], r["err_D"],
                 r["chi2"] / r["ndf"], r["segundos"]), flush=True)
        # se guarda despues de CADA corrida: si algo se corta a la corrida 20,
        # las 19 anteriores no se pierden
        json.dump(res, open(SALIDA, "w", encoding="utf-8"), indent=1)
    return res


# --------------------------------------------------------------------------
# Test B
# --------------------------------------------------------------------------
def recta_ponderada(x, y, err):
    """y = a + m x por minimos cuadrados con pesos 1/err^2 -> a, m, ea, em."""
    w = 1.0 / err ** 2
    S, Sx, Sxx = w.sum(), (w * x).sum(), (w * x * x).sum()
    Sy, Sxy = (w * y).sum(), (w * x * y).sum()
    det = S * Sxx - Sx * Sx
    a = (Sxx * Sy - Sx * Sxy) / det
    m = (S * Sxy - Sx * Sy) / det
    return a, m, np.sqrt(Sxx / det), np.sqrt(S / det)


def ley_de_ganancia(V, D, err):
    """Recta contra ley de potencia. Devuelve b, su error, y los dos chi2."""
    V, D, err = (np.asarray(v, float) for v in (V, D, err))
    a, m, _, _ = recta_ponderada(V, D, err)
    chi2_recta = float(np.sum(((D - a - m * V) / err) ** 2))
    V0 = float(np.exp(np.mean(np.log(V))))
    lD0, b, _, eb = recta_ponderada(np.log(V / V0), np.log(D), err / D)
    pred = np.exp(lD0) * (V / V0) ** b
    chi2_pot = float(np.sum(((D - pred) / err) ** 2))
    gdl = V.size - 2
    return dict(n=int(V.size), b=float(b), err_b=float(eb),
                # el error de b de arriba NO esta escalado por sqrt(chi2/gdl)
                # (igual que analisis/voltaje.py); cuando chi2 >> gdl es corto
                err_b_escalado=float(eb * np.sqrt(max(chi2_pot / gdl, 1.0)))
                if gdl > 0 else float("nan"),
                chi2_potencia=chi2_pot, chi2_recta=chi2_recta, gdl=int(gdl),
                # relativo a la D MEDIDA, como la tabla del informe
                peor_residuo=float(np.max(np.abs(D - pred) / D)),
                V0=V0, D0=float(np.exp(lD0)))


def test_B(res):
    grupos = {}
    for r in res.values():
        if r["descartada"]:
            continue
        grupos.setdefault((r["detector"], r["fuente"]), []).append(r)
    out = {}
    for (det, f), rs in sorted(grupos.items()):
        rs = sorted(rs, key=lambda r: r["hv"])
        if len(rs) < 2:
            continue
        out[(det, f)] = dict(ley_de_ganancia([r["hv"] for r in rs],
                                             [r["D"] for r in rs],
                                             [r["err_D"] for r in rs]),
                             claves=[r["clave"] for r in rs])
    return out


# --------------------------------------------------------------------------
# Test A
# --------------------------------------------------------------------------
_Q = {}


def quenching(fuente):
    """<Lvis/Edep> de la simulacion por encima de 0.7 E_ref."""
    if fuente not in _Q:
        import uproot
        E_ref = float(min(fisica.e_max(e) for e, _ in fisica.FUENTES[fuente]))
        d = uproot.open(plantilla.SIMULACIONES[fuente])[plantilla.ARBOL].arrays(
            ["Edep_keV", "Lvis_keVee"], library="np")
        m = d["Edep_keV"] > FRAC_QUENCH * E_ref
        _Q[fuente] = float(np.mean(d["Lvis_keVee"][m] / d["Edep_keV"][m]))
    return _Q[fuente]


def test_A(res):
    esperado = quenching("Cs-137") / quenching("Co-60")
    E_cs = float(fisica.e_max(661.657))
    E_co = float(fisica.e_max(1173.228))
    puntos = {}
    for r in res.values():
        if not r["descartada"]:
            puntos.setdefault((r["detector"], r["hv"]), {})[r["fuente"]] = r
    out = {}
    for (det, hv), fs in sorted(puntos.items()):
        if len(fs) < 2:
            out[(det, hv)] = None
            continue
        cs, co = fs["Cs-137"], fs["Co-60"]
        D1, s1, D2, s2 = cs["D"], cs["err_D"], co["D"], co["err_D"]
        w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
        razon = D1 / D2
        err_razon = razon * np.hypot(s1 / D1, s2 / D2)
        k = E_cs * E_co / (E_co - E_cs)
        out[(det, hv)] = dict(
            D_cs=D1, err_cs=s1, D_co=D2, err_co=s2,
            D_conj=(w1 * D1 + w2 * D2) / (w1 + w2), err_conj=1 / np.sqrt(w1 + w2),
            razon=razon, err_razon=err_razon,
            desacuerdo=razon / esperado - 1, err_desacuerdo=err_razon / esperado,
            c0=k * (D1 - D2), err_c0=k * np.hypot(s1, s2),
            sesiones="misma" if cs["fecha"] == co["fecha"] else "distintas",
            fechas=(cs["fecha"], co["fecha"]))
    return out, esperado


# --------------------------------------------------------------------------
# Salida
# --------------------------------------------------------------------------
def imprimir(res):
    b, (a, esperado) = test_B(res), test_A(res)
    print()
    print("  TEST B: ley de ganancia D = D0 (V/V0)^b")
    print("  %-4s %-7s %3s %16s %12s %12s %8s" % ("det", "fuente", "nV", "b",
                                                  "chi2 pot.", "chi2 recta",
                                                  "peor res."))
    for (det, f), g in b.items():
        print("  %-4s %-7s %3d %8.2f +- %.2f %12.1f %12.1f %7.1f %%"
              % (det, f, g["n"], g["b"], g["err_b"], g["chi2_potencia"],
                 g["chi2_recta"], 100 * g["peor_residuo"]))
    print()
    print("  TEST A: las dos fuentes al mismo punto (esperado Cs/Co = %.4f por "
          "quenching)" % esperado)
    print("  %-10s %8s %8s %8s %13s %15s %12s %s"
          % ("punto", "D(Cs)", "D(Co)", "D conj", "Cs/Co", "desacuerdo", "c0", "sesiones"))
    for (det, hv), t in a.items():
        if t is None:
            print("  %s @ %d V   una sola fuente" % (det, hv))
            continue
        print("  %s @ %d V %8.4f %8.4f %8.4f %6.3f+-%.3f %+6.1f +- %.1f %% %+5.0f +- %2.0f  %s"
              % (det, hv, t["D_cs"], t["D_co"], t["D_conj"], t["razon"],
                 t["err_razon"], 100 * t["desacuerdo"], 100 * t["err_desacuerdo"],
                 t["c0"], t["err_c0"], t["sesiones"]))
    return b, a, esperado


def figura(res, b, a, esperado):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ejes = plt.subplots(1, 2, figsize=(14, 5.2))
    ax = ejes[0]
    colores = dict(A="tab:blue", B="tab:orange", C="tab:green", D="tab:red")
    for (det, f), g in b.items():
        rs = [res[k] for k in g["claves"]]
        V = np.array([r["hv"] for r in rs])
        D = np.array([r["D"] for r in rs])
        E = np.array([r["err_D"] for r in rs])
        mk = "o" if f == "Co-60" else "s"
        ax.errorbar(V, D, E, fmt=mk, color=colores[det], ms=5, capsize=2,
                    mfc="white" if f == "Cs-137" else colores[det],
                    label="%s %s: b = %.2f" % (det, f, g["b"]))
        vv = np.linspace(V.min() - 20, V.max() + 20, 50)
        ax.plot(vv, g["D0"] * (vv / g["V0"]) ** g["b"], color=colores[det], lw=0.8,
                ls="-" if f == "Co-60" else "--")
    ax.set(xscale="log", yscale="log", xlabel="alto voltaje [V]",
           ylabel="D [mV*ns/keV]", title="test B: D(V) y ley de potencia")
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    ax.xaxis.set_major_locator(FixedLocator([1550, 1600, 1650, 1700]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    ax.legend(fontsize=7, ncol=2)
    ax = ejes[1]
    for (det, hv), t in a.items():
        if t is None:
            continue
        ax.errorbar(hv + {"A": -9, "B": -3, "C": 3, "D": 9}[det], t["razon"],
                    t["err_razon"], fmt="o", color=colores[det], capsize=2,
                    label=det if hv == 1600 else None)
    ax.axhline(esperado, color="k", ls=":", label="quenching: %.3f" % esperado)
    ax.set(xlabel="alto voltaje [V]", ylabel="D(Cs-137) / D(Co-60)",
           title="test A: cociente entre fuentes al mismo punto")
    ax.legend(fontsize=8)
    for ax in ejes:
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa9_lote.png")
    fig.savefig(salida, dpi=120)
    print("  figura: %s" % salida)


def main():
    if "--tests" in sys.argv:
        res = json.load(open(SALIDA, encoding="utf-8"))
    else:
        res = correr(rehacer="--rehacer" in sys.argv)
    b, a, esperado = imprimir(res)
    figura(res, b, a, esperado)


if __name__ == "__main__":
    main()
