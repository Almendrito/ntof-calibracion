"""
correr.py — el analisis de resultados del paper, aplicado a nuestras corridas.

    python replica_paper/correr.py --auditoria   solo los numeros del paper
    python replica_paper/correr.py --rapido      sin los 10 bloques (~1 min)
    python replica_paper/correr.py --cierre      cierre con verdad conocida
    python replica_paper/correr.py               todo (~6 min la primera vez)

Que hace, en el orden del paper:

  Seccion II    NCC -> SA -> borde, sobre cada corrida (`metodo.py`), y los dos
                metodos convencionales sobre la MISMA ventana.
  Seccion II-C  "Ten independent measurements": la corrida partida en 10 bloques
                (`datos.bloques`), de donde salen las barras de error.
  Seccion III   recta canal-energia por punto de operacion, con SD y R^2;
                similitud coseno; errores relativos por punto; reduccion del
                error medio; robustez frente al intervalo.

Lo que NO se puede replicar con lo medido esta listado al final del informe, con
el motivo. No se rellena con nada.
"""

import argparse
import os
import sys
import time

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import auditoria
from replica_paper import convencionales as conv
from replica_paper import datos
from replica_paper import ecuaciones as ec
from replica_paper import metodo
from replica_paper import resultados as R
from replica_paper import rutas

# Los intervalos de la Fig. 5, en unidades del canal del borde. Salen de los seis
# paneles del paper con su borde en 1535: [1500,1600] -> [0.977, 1.042], etc.
INTERVALOS_FIG5 = [(0.977, 1.042), (0.912, 1.042), (0.782, 1.042),
                   (0.782, 1.173), (0.651, 1.303), (0.521, 1.433)]

METODOS = ("NCC-SA", "Gauss FWHM", "Diferencia")

# Bins minimos para que un intervalo de la Fig. 5 sea evaluable. Por debajo de
# esto no hay con que ajustar 3 parametros ni con que derivar, y el numero que
# saliera seria ruido. Se aplica IGUAL a los tres metodos, para no favorecer a
# ninguno.
MIN_BINS_ROBUSTEZ = 30


# --------------------------------------------------------------------------
def _borde_de(res):
    """Canal del borde que reporta un ajuste NCC-SA."""
    return res["B"] if res["tipo"] == "un_borde" else res["B1"]


def tres_metodos(canales, cuentas, lineas, intervalo=None, sigma_rel=None,
                 ligar_origen=True, verbose=False):
    """
    Los tres metodos sobre el mismo espectro y la misma ventana.

    La ventana sale del NCC (Seccion II-B). Darsela tambien a los dos metodos
    convencionales es lo unico que hace la comparacion justa: el paper no dice
    que ventana usan, y si cada uno eligiera la suya se estaria comparando la
    eleccion de ventana y no el metodo.
    """
    r = metodo.ajustar(canales, cuentas, lineas, intervalo=intervalo,
                       sigma_rel=sigma_rel,
                       **({"ligar_origen": ligar_origen} if len(lineas) == 2 else {}))
    iv = r["intervalo"]
    try:
        g = conv.metodo_gauss(canales, cuentas, iv)["borde"]
    except Exception:                                    # noqa: BLE001
        g = float("nan")
    try:
        d = conv.metodo_diferencia(canales, cuentas, iv)["borde"]
    except Exception:                                    # noqa: BLE001
        d = float("nan")
    return {"NCC-SA": _borde_de(r), "Gauss FWHM": g, "Diferencia": d}, r


def diagnostico_dos_bordes(e, r, verbose=True):
    """
    Por que el Co-60 no entrega dos puntos de calibracion, con numeros.

    Dos medidas, las dos dentro del marco del paper:
      1. Que pasa con D libre (Ec. 11) contra D = B1/E_max,1 (que es lo que
         impone la Ec. 8 del mismo paper).
      2. A cuantos sigma quedan los dos bordes. Es la cantidad que decide si son
         separables, y explica por que al paper si le funciona: su BC408 tiene
         sigma/B ~ 0.06 y este detector ~0.13.
    Mas la Ec. 12 impresa contra la corregida, sobre los mismos datos.
    """
    eg1, eg2 = e["lineas"]
    if verbose:
        print("      diagnostico del solapamiento (Ec. 11-13):")
    salida = {}
    for lig in (False, True):
        f = metodo.ajustar_dos_bordes(e["canales"], e["cuentas"], eg1, eg2,
                                      intervalo=r["intervalo"],
                                      sigma_rel=r["sigma_rel_s2"],
                                      ligar_origen=lig)
        salida[lig] = f
        if verbose:
            sep = (f["B2"] - f["B1"]) / max(f["sigma1"], 1e-9)
            print("        D %-8s B1 = %7.1f  B2 = %7.1f  D = %.4f  "
                  "A2/A1 = %7.3f  sep = %.2f sigma  chi2/ndf = %.2f"
                  % ("libre" if not lig else "= B1/Emax", f["B1"], f["B2"],
                     f["D"], f["A2"] / max(f["A1"], 1e-30), sep, f["chi2_ndf"]))
    if verbose:
        print("        (con Ec. 12 impresa, dE de las lineas gamma:)")
    g = metodo.ajustar_dos_bordes(e["canales"], e["cuentas"], eg1, eg2,
                                  intervalo=r["intervalo"],
                                  sigma_rel=r["sigma_rel_s2"],
                                  ligar_origen=True, delta_e="gammas")
    salida["gammas"] = g
    if verbose:
        print("        D = B1/Emax  B1 = %7.1f  B2 = %7.1f  D = %.4f  "
              "chi2/ndf = %.2f" % (g["B1"], g["B2"], g["D"], g["chi2_ndf"]))
    return salida


def analizar(clave, bloques=10, verbose=True):
    """Todo lo que se le puede sacar a una corrida."""
    e = datos.espectro(clave)
    t0 = time.time()
    bordes, r = tres_metodos(e["canales"], e["cuentas"], e["lineas"])
    if verbose:
        print("  %-12s %s" % (clave, e["corrida"].descripcion))
        print("      fuente %s  ->  bordes teoricos %s keV"
              % (e["fuente"], ", ".join("%.1f" % ec.e_max(x) for x in e["lineas"])))
        print("      NCC: intervalo [%.0f, %.0f] mV*ns   rho = %.5f   "
              "(sigma_rel de S2 = %.2f)"
              % (r["intervalo"][0], r["intervalo"][1], r["rho"], r["sigma_rel_s2"]))
        if r["tipo"] == "un_borde":
            print("      SA : B = %.1f   sigma = %.1f (sigma/B = %.3f)   "
                  "chi2/ndf = %.2f"
                  % (r["B"], r["sigma"], r["sigma"] / r["B"], r["chi2_ndf"]))
        else:
            print("      SA : B1 = %.1f  B2 = %.1f  D = %.4f  sigma = (%.0f, %.0f)"
                  % (r["B1"], r["B2"], r["D"], r["sigma1"], r["sigma2"]))
            print("           A2/A1 = %.3f   chi2/ndf = %.2f   (Ec.12 con dE de %s%s)"
                  % (r["A2"] / max(r["A1"], 1e-30), r["chi2_ndf"], r["delta_e"],
                     ", D ligada a la Ec. 8" if r.get("ligar_origen") else ""))
        for m in METODOS:
            print("      %-11s borde = %8.1f mV*ns" % (m, bordes[m]))
        print("      (%.0f s)" % (time.time() - t0))

    if verbose:
        print("      barrido del sigma de S2 (el paper no lo especifica; gana el")
        print("      de mayor rho, pero ojo: fin_rel depende de sigma_rel, asi que")
        print("      los rho no se midieron sobre el mismo conjunto de ventanas):")
        print("        sigma_rel   rho        B (NCC)")
        for sr, rho, Bn in r["barrido_sigma_rel"]:
            print("          %.2f      %.5f   %8.1f%s"
                  % (sr, rho, Bn, "   <- elegido" if sr == r["sigma_rel_s2"] else ""))

    diag = None
    if len(e["lineas"]) == 2:
        diag = diagnostico_dos_bordes(e, r, verbose=verbose)

    # --- las "diez mediciones independientes" -----------------------------
    disp = None
    if bloques:
        if verbose:
            print("      %d bloques..." % bloques, end="", flush=True)
        vals = {m: [] for m in METODOS}
        for b in datos.bloques(clave, n=bloques):
            try:
                bb, _ = tres_metodos(b["canales"], b["cuentas"], b["lineas"],
                                     intervalo=r["intervalo"],
                                     sigma_rel=r["sigma_rel_s2"])
            except Exception:                            # noqa: BLE001
                bb = {m: float("nan") for m in METODOS}
            for m in METODOS:
                vals[m].append(bb[m])
        disp = {}
        for m in METODOS:
            v = np.array(vals[m], float)
            ok = np.isfinite(v)
            disp[m] = dict(valores=v, media=float(np.mean(v[ok])) if ok.any() else np.nan,
                           sd=float(np.std(v[ok], ddof=1)) if ok.sum() > 1 else np.nan,
                           n=int(ok.sum()))
        if verbose:
            print(" listo (%.0f s)" % (time.time() - t0))
            npb = e["n_pulsos"] / float(bloques)
            print("      cada bloque tiene ~%.0f pulsos (1/%d de la corrida). Las"
                  % (npb, bloques))
            print("      diez mediciones del paper tienen cada una la estadistica")
            print("      COMPLETA, asi que estos SD son mas grandes que los suyos por")
            print("      un factor ~sqrt(%d) = %.1f si la dispersion fuera solo de"
                  % (bloques, np.sqrt(bloques)))
            print("      conteo. La columna escalada se da con esa salvedad.")
            for m in METODOS:
                sd = disp[m]["sd"]
                b0 = abs(bordes[m]) if np.isfinite(bordes[m]) else float("nan")
                print("      %-11s media %8.1f   SD %7.2f (%5.2f%%)   "
                      "SD/sqrt(%d) = %6.2f (%5.2f%%)   (%d/%d bloques)"
                      % (m, disp[m]["media"], sd, 100 * sd / b0, bloques,
                         sd / np.sqrt(bloques), 100 * sd / np.sqrt(bloques) / b0,
                         disp[m]["n"], bloques))

    # --- robustez frente al intervalo (Fig. 5) ----------------------------
    B = bordes["NCC-SA"]
    paso = float(e["canales"][1] - e["canales"][0])
    ivs, por_intervalo, saltados = [], [], []
    for a_, b_ in INTERVALOS_FIG5:
        iv = (a_ * B, b_ * B)
        nbins = int(round((iv[1] - iv[0]) / paso))
        if nbins < MIN_BINS_ROBUSTEZ:
            saltados.append((iv, nbins))
            continue
        ivs.append(iv)
        try:
            bb, _ = tres_metodos(e["canales"], e["cuentas"], e["lineas"],
                                 intervalo=iv, sigma_rel=r["sigma_rel_s2"])
        except Exception:                                # noqa: BLE001
            bb = {m: float("nan") for m in METODOS}
        por_intervalo.append((iv, bb, nbins))
    rob = {}
    for m in METODOS:
        tabla = {tuple(iv): bb[m] for iv, bb, _ in por_intervalo}
        rob[m] = R.robustez_intervalo(lambda iv, t=tabla: t[tuple(iv)], ivs)
        rob[m]["nbins"] = [n for _, _, n in por_intervalo]
    if verbose:
        print("      robustez de intervalo (los intervalos relativos de la Fig. 5):")
        print("        intervalo [mV*ns]     bins    " +
              "".join("%-13s" % m for m in METODOS))
        for iv, bb, nb in por_intervalo:
            print("        [%7.0f, %7.0f]   %4d    " % (iv[0], iv[1], nb) +
                  "".join("%-13.1f" % bb[m] for m in METODOS))
        for iv, nb in saltados:
            print("        [%7.0f, %7.0f]   %4d    salteado: menos de %d bins"
                  % (iv[0], iv[1], nb, MIN_BINS_ROBUSTEZ))
        print("        %-24s" % "dispersion (todos)" +
              "".join("%-13s" % ("%.2f%%" % (100 * rob[m]["dispersion_rel"]))
                      for m in METODOS))
        # El paper tiene el mismo efecto: su intervalo mas angosto ([1500,1600])
        # es el unico que se aparta (1521 contra 1535-1540). Se reporta la
        # dispersion sin el mas angosto, igual que en `auditoria.fig5`.
        for m in METODOS:
            v = rob[m]["valores"][1:]
            ok = np.isfinite(v)
            rob[m]["dispersion_sin_el_mas_angosto"] = (
                float((v[ok].max() - v[ok].min()) / np.mean(v[ok]))
                if ok.sum() > 1 else float("nan"))
        print("        %-24s" % "  sin el mas angosto" +
              "".join("%-13s" % ("%.2f%%" % (100 * rob[m]["dispersion_sin_el_mas_angosto"]))
                      for m in METODOS))
        print()

    return dict(clave=clave, espectro=e, ajuste=r, bordes=bordes,
                dispersion=disp, robustez=rob, punto=e["punto"],
                lineas=e["lineas"], diagnostico=diag)


# --------------------------------------------------------------------------
def puntos_de_calibracion(analisis):
    """
    Los pares (canal, energia) que cada metodo aporta, por punto de operacion.

    Con una fuente de DOS lineas los tres metodos entregan un solo canal: el
    paper mismo dice que los convencionales "failed to resolve overlapping
    edges", y en nuestros datos al NCC-SA le pasa lo mismo (una de las dos
    amplitudes de la Ec. 11 se va a cero). Ese unico canal se asigna al borde de
    MENOR energia, que es la componente que el ajuste conserva, y el otro borde
    queda sin punto. Es la misma regla para los tres metodos.
    """
    puntos = []
    for a in analisis:
        E1 = float(ec.e_max(a["lineas"][0]))
        sd = a.get("dispersion")
        for m in METODOS:
            puntos.append(dict(clave=a["clave"], punto=a["punto"], metodo=m,
                               canal=float(a["bordes"][m]), energia=E1,
                               sd_canal=(float(sd[m]["sd"]) if sd else float("nan")),
                               fuente=a["espectro"]["fuente"]))
    return puntos


def seccion_iii(puntos, verbose=True):
    """Seccion III: recta, SD, R^2, similitud coseno, errores relativos."""
    grupos = R.rectas_por_punto(puntos)
    if verbose:
        print("=" * 78)
        print("Seccion III - recta canal-energia, por punto de operacion")
        print("=" * 78)
    rectas = {}
    for (punto, met), g in sorted(grupos.items()):
        if g["recta"] is None:
            if verbose:
                print("  %-14s  %-11s  %d punto(s): %s"
                      % (punto, met, len(g["puntos"]), g["motivo"]))
            continue
        r = g["recta"]
        rectas[(punto, met)] = r
        if verbose:
            print("  %-14s  %-11s  n = %d   y = %.4f x %+.2f   SD = %.4f   R2 = %.6f"
                  % (punto, met, r["n"], r["m"], r["b"], r["sd"], r["r2"]))
            if g["motivo"]:
                print("                 %s" % g["motivo"])
    if verbose and not rectas:
        print()
        print("  NINGUN punto de operacion tiene dos o mas bordes medidos, asi que")
        print("  no hay recta que ajustar. Es la Fig. 3 del paper, y no se puede")
        print("  hacer con lo medido: ver la lista del final.")
    return rectas


def fig4(analisis, verbose=True):
    """
    Fig. 4: error relativo por punto y por metodo, con barra = SD entre bloques.

    Sin recta de calibracion (hace falta mas de un borde por punto de operacion)
    el error relativo no esta definido: la energia reconstruida saldria de la
    misma medicion que se quiere evaluar. Se dice y no se inventa.
    """
    if verbose:
        print()
        print("=" * 78)
        print("Fig. 4 - error relativo por metodo")
        print("=" * 78)
    puntos = puntos_de_calibracion(analisis)
    grupos = R.rectas_por_punto(puntos)
    utiles = {k: g for k, g in grupos.items() if g["recta"] is not None
              and g["recta"]["gdl"] > 0}
    if not utiles:
        if verbose:
            print("  No calculable: hace falta una recta de calibracion con al menos")
            print("  3 bordes en un mismo punto de operacion (con 2 la recta pasa")
            print("  exacta y todos los errores relativos dan 0 por construccion).")
            print("  Lo que SI se puede reportar con lo medido es la DISPERSION entre")
            print("  bloques, que es la barra de error de la Fig. 4 sin el valor")
            print("  central:")
            print()
            print("  corrida       metodo        borde [mV*ns]   SD entre bloques")
            for a in analisis:
                if not a.get("dispersion"):
                    continue
                for m in METODOS:
                    d = a["dispersion"][m]
                    print("  %-12s  %-11s   %10.1f      %7.2f  (%.2f%%)"
                          % (a["clave"], m, a["bordes"][m], d["sd"],
                             100 * d["sd"] / max(abs(a["bordes"][m]), 1e-9)))
        return None
    salida = {}
    for (punto, met), g in sorted(utiles.items()):
        lista = g["puntos"]
        e = R.errores_de_metodo([p["canal"] for p in lista],
                                [p["sd_canal"] for p in lista],
                                [p["energia"] for p in lista], g["recta"])
        salida[(punto, met)] = e
        if verbose:
            print("  %-14s %-11s" % (punto, met))
            for c, et, er, b in zip(e["canal"], e["e_teorica"], e["error_rel"],
                                    e["barra_rel"]):
                print("      canal %8.1f   E_teo %7.1f keV   error %+6.2f%% "
                      "+- %.2f%%" % (c, et, er, b))
    # reducciones del error medio
    if verbose:
        print()
        print("  Reduccion del error medio de NCC-SA:")
    for punto in sorted({p for p, _ in salida}):
        base = salida.get((punto, "NCC-SA"))
        if base is None:
            continue
        for otro in ("Diferencia", "Gauss FWHM"):
            o = salida.get((punto, otro))
            if o is None:
                continue
            red = R.reduccion_media(o["error_rel"], base["error_rel"])
            if verbose:
                print("      %-14s contra %-11s : %+.1f%%" % (punto, otro, red))
    return salida


# --------------------------------------------------------------------------
def figuras(analisis, rectas):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(analisis)
    fig, ax = plt.subplots(2, n, figsize=(5.0 * n, 7.4), squeeze=False)
    for k, a in enumerate(analisis):
        e, r = a["espectro"], a["ajuste"]
        a0 = ax[0][k]
        a0.step(e["canales"], e["cuentas"], where="mid", color="0.55", lw=0.8,
                label="neto")
        a0.axvspan(r["intervalo"][0], r["intervalo"][1], color="tab:orange",
                   alpha=0.13, label="intervalo NCC")
        a0.plot(r["canales"], r["modelo"], color="tab:red", lw=1.4,
                label="Ec. 8-13 (SA)")
        for m, col, ls in zip(METODOS, ("tab:red", "tab:green", "tab:blue"),
                              ("-", "--", ":")):
            a0.axvline(a["bordes"][m], color=col, ls=ls, lw=1.2, label=m)
        a0.set_title("%s  (%s)" % (a["clave"], a["punto"]), fontsize=10)
        a0.set_xlabel("canal [mV*ns]")
        a0.set_ylabel("cuentas netas")
        a0.legend(fontsize=7)
        a0.grid(alpha=0.3)

        a1 = ax[1][k]
        for m, col in zip(METODOS, ("tab:red", "tab:green", "tab:blue")):
            v = np.array([b for _, b in a["robustez"][m]["filas"]], float)
            a1.plot(range(1, len(v) + 1), v, "o-", color=col, lw=1.1, ms=4,
                    label="%s (SD %.1f)" % (m, a["robustez"][m]["sd"]))
        a1.set_xlabel("intervalo de la Fig. 5 (1 = el mas angosto)")
        a1.set_ylabel("borde [mV*ns]")
        a1.set_title("robustez frente al intervalo", fontsize=10)
        a1.legend(fontsize=7)
        a1.grid(alpha=0.3)
    fig.tight_layout()
    f1 = rutas.figura("replica_ajustes.png")
    fig.savefig(f1, dpi=130)

    fig2, ax2 = plt.subplots(figsize=(6.4, 5.0))
    hay = False
    for (punto, met), r in sorted(rectas.items()):
        ax2.plot(r["canal"], r["energia"], "o", label="%s / %s" % (punto, met))
        x = np.linspace(0, max(r["canal"]) * 1.1, 50)
        ax2.plot(x, r["m"] * x + r["b"], "--", lw=1.0)
        hay = True
    if not hay:
        for a in analisis:
            E1 = float(ec.e_max(a["lineas"][0]))
            ax2.plot([a["bordes"]["NCC-SA"]], [E1], "o", label=a["punto"])
        ax2.text(0.03, 0.95, "un solo borde por punto de operacion:\nno hay recta"
                 " que ajustar (ver el informe)", transform=ax2.transAxes,
                 va="top", fontsize=9,
                 bbox=dict(fc="w", ec="0.7"))
    ax2.set_xlabel("canal [mV*ns]")
    ax2.set_ylabel("energia del borde [keV]")
    ax2.set_title("Fig. 3 del paper, con nuestras corridas")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    f2 = rutas.figura("replica_fig3.png")
    fig2.savefig(f2, dpi=130)
    return f1, f2


# --------------------------------------------------------------------------
def no_replicable(analisis):
    print()
    print("=" * 78)
    print("Lo que NO se puede replicar con lo medido, y por que")
    print("=" * 78)
    puntos = {}
    for a in analisis:
        puntos.setdefault(a["punto"], []).append(a["espectro"]["fuente"])
    print("  1. Fig. 3 (recta canal-energia con SD y R^2).")
    print("     El paper la arma con 4 fuentes medidas con el MISMO detector y la")
    print("     MISMA ganancia: 5 o 6 bordes por recta. Aca hay:")
    for p, fs in sorted(puntos.items()):
        print("       %-16s %s" % (p, ", ".join(sorted(set(fs)))))
    print("     o sea una fuente por punto de operacion. Con un solo borde no hay")
    print("     recta; con dos, la recta pasa exacta y SD = 0, R^2 = 1 no miden")
    print("     nada. Hace falta medir Cs-137 y Co-60 (y ojala Na-22 y Mn-54) al")
    print("     mismo alto voltaje.")
    print("     Ojo con un detalle estructural: la Ec. 8 mapea canal a energia POR")
    print("     EL ORIGEN, asi que todos los bordes que salen de UNA fuente caen")
    print("     exactos sobre una recta por el origen. La recta solo mide algo si")
    print("     los puntos vienen de fuentes DISTINTAS.")
    print()
    print("  2. Fig. 4 (error relativo por metodo).")
    print("     Depende de la recta del punto 1: la energia reconstruida sale de")
    print("     ella. Sin recta, el error relativo no esta definido.")
    print("     La BARRA de error si se puede: es la dispersion entre mediciones,")
    print("     y esta calculada arriba con los 10 bloques.")
    print()
    print("  3. La comparacion contra la calibracion por fotopico (y su similitud")
    print("     coseno). Un centellador plastico no tiene fotopico. El paper")
    print("     tampoco la hace para BC408 y EJ309, por lo mismo.")
    print()
    print("  4. Las 10 mediciones INDEPENDIENTES.")
    print("     Aca son 10 bloques de una misma corrida, y eso las cambia en dos")
    print("     direcciones opuestas:")
    print("       - hacia abajo: entre bloques solo varia la estadistica de")
    print("         conteo, mientras que entre las diez corridas del paper varian")
    print("         ademas el montaje, la posicion de la fuente y la ganancia;")
    print("       - hacia arriba: cada bloque tiene la DECIMA parte de las cuentas,")
    print("         asi que su dispersion estadistica es ~sqrt(10) veces la de una")
    print("         medicion completa.")
    print("     Por eso el SD se reporta crudo y escalado por sqrt(10), y ninguno")
    print("     de los dos es directamente la barra del paper.")
    print()
    print("  5. El segundo borde del Co-60 (1118.1 keV).")
    for a in analisis:
        d = a.get("diagnostico")
        if not d:
            continue
        libre, ligada = d[False], d[True]
        print("     Con Theta = {A1,A2,B1,s1,s2,D} (Ec. 11) D queda LIBRE ademas de")
        print("     B1, y el modelo se sobre-parametriza respecto de su propia")
        print("     Ec. 8, que ya fija el mapa por el origen. Medido en %s:"
              % a["clave"])
        print("       D libre         : D = %.4f  A2/A1 = %.2f  chi2/ndf = %.2f"
              % (libre["D"], libre["A2"] / max(libre["A1"], 1e-30),
                 libre["chi2_ndf"]))
        print("       D = B1/E_max,1  : D = %.4f  A2/A1 = %.2f  chi2/ndf = %.2f"
              % (ligada["D"], ligada["A2"] / max(ligada["A1"], 1e-30),
                 ligada["chi2_ndf"]))
        print("     o sea %.0f%% de diferencia en D con chi2 casi igual, y un"
              % abs(100 * (libre["D"] / ligada["D"] - 1)))
        print("     cociente de amplitudes sin sentido fisico (el Co-60 emite las")
        print("     dos lineas con intensidades practicamente iguales).")
        sep = (ligada["B2"] - ligada["B1"]) / max(ligada["sigma1"], 1e-9)
        print("     Los dos bordes quedan a %.1f sigma1. El cierre con verdad" % sep)
        print("     conocida de `metodo.py` dice que a 2 sigma D sale sesgada ~-14%")
        print("     con +-15% de dispersion, y que recien a 4 sigma es utilizable.")
        if a.get("dispersion"):
            print("     Consistente con eso, el SD de B1 entre los 10 bloques es del")
            print("     %.0f%%." % (100 * a["dispersion"]["NCC-SA"]["sd"]
                                    / max(abs(a["bordes"]["NCC-SA"]), 1e-9)))


# --------------------------------------------------------------------------
def cierre(D_real=1.30, sigma_rel=0.07, semilla=2026):
    """
    Cierre con VERDAD CONOCIDA de todo el analisis de resultados.

    Genera dos espectros sinteticos -Cs-137 y Co-60- con la MISMA D, que es la
    situacion que hoy no tenemos medida, y corre la Seccion III completa. Sirve
    para dos cosas: verificar que el camino de la recta, el SD, el R^2, la
    similitud coseno y los errores relativos funciona, y dejarlo listo para
    cuando existan las dos corridas al mismo voltaje.

    D_real    : mV*ns por keV, comun a las dos fuentes.
    sigma_rel : resolucion relativa en el borde.

    Devuelve (analisis, rectas) para que `gui.py` pueda dibujar lo mismo que
    imprime, sin volver a generar los espectros por su cuenta.
    """
    rutas.salida_utf8()
    rng = np.random.default_rng(semilla)
    print("=" * 78)
    print("Cierre con verdad conocida: dos fuentes al MISMO punto de operacion")
    print("=" * 78)
    print("  D verdadera = %.4f mV*ns/keV,  sigma/B = %.3f, recta por el origen"
          % (D_real, sigma_rel))
    print()

    fuentes = [("Cs-137", (661.657,)), ("Co-60", (1173.228, 1332.492))]
    analisis = []
    for nombre, lineas in fuentes:
        emax = [float(ec.e_max(x)) for x in lineas]
        ch = np.linspace(5.0, 1.9 * D_real * max(emax), 380)
        if len(lineas) == 1:
            limpio = ec.modelo(ch, lineas[0], 1.0, D_real * emax[0],
                               sigma_rel * D_real * emax[0])
        else:
            limpio = ec.modelo_dos_bordes(
                ch, lineas[0], lineas[1], 1.0, 0.95, D_real * emax[0],
                sigma_rel * D_real * emax[0], sigma_rel * D_real * emax[1], D_real)
        limpio = limpio / limpio.max() * 900.0
        cuentas = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
        bordes, r = tres_metodos(ch, cuentas, lineas)
        # "diez mediciones": diez realizaciones de Poisson del mismo espectro
        vals = {m: [] for m in METODOS}
        for _ in range(10):
            c = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
            bb, _ = tres_metodos(ch, c, lineas, intervalo=r["intervalo"],
                                 sigma_rel=r["sigma_rel_s2"])
            for m in METODOS:
                vals[m].append(bb[m])
        disp = {m: dict(valores=np.array(vals[m], float),
                        media=float(np.nanmean(vals[m])),
                        sd=float(np.nanstd(vals[m], ddof=1)),
                        n=int(np.isfinite(vals[m]).sum())) for m in METODOS}
        print("  %s  bordes teoricos %s keV"
              % (nombre, ", ".join("%.1f" % x for x in emax)))
        for m in METODOS:
            print("    %-11s borde = %8.1f   canal verdadero = %8.1f   "
                  "(%+.2f%%)   SD entre 10 = %.1f"
                  % (m, bordes[m], D_real * emax[0],
                     100 * (bordes[m] / (D_real * emax[0]) - 1), disp[m]["sd"]))
        analisis.append(dict(clave=nombre,
                             espectro=dict(fuente=nombre, canales=ch,
                                           cuentas=cuentas, limpio=limpio),
                             ajuste=r, bordes=bordes, dispersion=disp,
                             punto="sintetico @ D = %.2f" % D_real,
                             lineas=lineas, robustez={},
                             canal_verdadero=D_real * emax[0]))
    print()
    puntos = puntos_de_calibracion(analisis)
    rectas = seccion_iii(puntos)
    print()
    print("  Contra la verdad (m = 1/D = %.4f keV/canal, b = 0):" % (1.0 / D_real))
    for (p, met), r in sorted(rectas.items()):
        print("    %-11s m = %.4f (%+.2f%%)   b = %+.2f keV   SD = %.3f   R2 = %.6f"
              % (met, r["m"], 100 * (r["m"] * D_real - 1), r["b"], r["sd"], r["r2"]))
    print()
    print("  Similitud coseno entre las relaciones energia-canal de cada par:")
    claves = sorted(rectas)
    for i in range(len(claves)):
        for j in range(i + 1, len(claves)):
            s = R.similitudes(rectas[claves[i]], rectas[claves[j]])
            print("    %-11s vs %-11s  en los canales: %.8f   coeficientes: %.8f"
                  % (claves[i][1], claves[j][1], s["en_los_canales"],
                     s["coeficientes"]))
    print()
    fig4(analisis)
    print()
    print("  Con 2 puntos la recta pasa exacta. El cierre demuestra que el camino")
    print("  corre; para que SD, R^2 y los errores relativos midan algo hacen")
    print("  falta 3 o mas bordes en el mismo punto de operacion.")
    return analisis, rectas


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--auditoria", action="store_true",
                    help="solo la auditoria de los numeros publicados")
    ap.add_argument("--rapido", action="store_true",
                    help="sin los 10 bloques")
    ap.add_argument("--cierre", action="store_true",
                    help="cierre con verdad conocida del analisis de resultados")
    ap.add_argument("--corridas", nargs="*", default=None,
                    help="claves del catalogo (por defecto, las que tengan fondo)")
    ap.add_argument("--bloques", type=int, default=10)
    args = ap.parse_args()

    rutas.salida_utf8()
    if args.cierre:
        cierre()
        return
    auditoria.todo()
    if args.auditoria:
        return

    claves = args.corridas
    if claves is None:
        claves = datos.con_fondo()
    print()
    print("=" * 78)
    print("Seccion II - NCC -> SA sobre nuestras corridas")
    print("=" * 78)
    analisis = [analizar(k, bloques=0 if args.rapido else args.bloques)
                for k in claves]

    puntos = puntos_de_calibracion(analisis)
    print()
    rectas = seccion_iii(puntos)
    fig4(analisis)
    no_replicable(analisis)
    f1, f2 = figuras(analisis, rectas)
    print()
    print("  figuras: %s" % rutas.relativa_o(f1))
    print("           %s" % rutas.relativa_o(f2))


if __name__ == "__main__":
    main()
