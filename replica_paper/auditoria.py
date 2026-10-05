"""
auditoria.py — reproducir los resultados publicados del paper con sus propios datos.

Todo lo que esta aca sale de las figuras y del texto del paper, nada de nuestras
mediciones. Sirve para dos cosas: fijar las definiciones que el paper no escribe
(que es "SD", por ejemplo) y ver si sus numeros cierran entre si.

Los puntos de la Fig. 3 estan rotulados uno por uno en la propia figura, asi que
no hay que digitalizar pixeles: se leen. La prueba de que se leyeron bien es que
el ajuste devuelve las rectas impresas a 4-5 cifras.

`python replica_paper/auditoria.py`
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import ecuaciones as ec
from replica_paper import resultados as R
from replica_paper import rutas

# --------------------------------------------------------------------------
# Fig. 3: puntos rotulados en la propia figura, (canal, energia [keV]).
# El nombre de la fuente es el que IMPRIME el paper; `fuente_real` es la que
# corresponde a esa energia. En el panel de NaI no coinciden.
# --------------------------------------------------------------------------
FIG3 = {
    "LaBr3:Ce": dict(
        fotopico=[("Na22", 249, 511), ("Cs137", 322, 662), ("Mn54", 406, 834),
                  ("Co60", 570, 1173), ("Na22", 617, 1275), ("Co60", 647, 1333)],
        borde=[("Cs137", 239, 477), ("Mn54", 314, 639), ("Co60", 469, 963),
               ("Na22", 518, 1061), ("Co60", 543, 1117)],
        impreso=dict(fotopico=(2.07, -4.76), borde=(2.10, -21.73),
                     sd_azul=2.183036, r2_azul=0.999924,
                     sd_rojo=1.609198, r2_rojo=0.999974),
        coseno=0.999993),
    "NaI:Tl": dict(
        fotopico=[("Na22", 667, 511), ("Na22", 846, 662), ("Cs137", 1050, 834),
                  ("Mn54", 1442, 1173), ("Na22", 1568, 1275), ("Co60", 1631, 1333)],
        borde=[("Co60", 620, 477), ("Na22", 801, 639), ("Mn54", 1207, 963),
               ("Co60", 1303, 1061)],
        impreso=dict(fotopico=(0.85, -58.88), borde=(0.84, -40.18),
                     sd_azul=7.274524, r2_azul=0.999055,
                     sd_rojo=1.981123, r2_rojo=0.999960),
        coseno=0.999999),
}

# Que fuente emite cada energia, para detectar rotulos cruzados.
FUENTE_DE_ENERGIA = {
    511: "Na22", 662: "Cs137", 834: "Mn54", 1173: "Co60", 1275: "Na22",
    1333: "Co60",                                   # fotopicos
    341: "Na22", 477: "Cs137", 639: "Mn54", 963: "Co60", 1061: "Na22",
    1117: "Co60",                                   # bordes Compton
}

# Fig. 5: intervalo de ajuste -> canal del borde que reporta cada panel.
FIG5 = [((1500, 1600), 1521), ((1400, 1600), 1538), ((1200, 1600), 1540),
        ((1200, 1800), 1537), ((1000, 2000), 1540), ((800, 2200), 1535)]

# Seccion III / Conclusion: reducciones de error medio que afirma el texto.
REDUCCIONES = {"BC408": dict(diferencia=80, gauss=85),
               "EJ309": dict(diferencia=44, gauss=66)}


def _linea(nombre, puntos):
    ch = [c for _, c, _ in puntos]
    en = [e for _, _, e in puntos]
    return R.ajuste_lineal(ch, en)


def fig3():
    """Reproduce los dos ajustes de cada panel y los compara con lo impreso."""
    print("=" * 78)
    print("Fig. 3 - recta canal-energia (ajuste de los puntos que el paper rotula)")
    print("=" * 78)
    salida = {}
    for det, d in FIG3.items():
        rf = _linea("fotopico", d["fotopico"])
        rb = _linea("borde", d["borde"])
        imp = d["impreso"]
        print("  %s" % det)
        print("    fotopico (azul)  n=%d  y = %.4f x %+.3f   RMS = %.6f  R2 = %.6f"
              % (rf["n"], rf["m"], rf["b"], rf["sd"], rf["r2"]))
        print("                          impreso: y = %.2f x %+.2f" % imp["fotopico"])
        print("    borde   (rojo)   n=%d  y = %.4f x %+.3f   RMS = %.6f  R2 = %.6f"
              % (rb["n"], rb["m"], rb["b"], rb["sd"], rb["r2"]))
        print("                          impreso: y = %.2f x %+.2f" % imp["borde"])
        print("    recuadro del paper:  SD_azul = %.6f  R2_azul = %.6f"
              % (imp["sd_azul"], imp["r2_azul"]))
        print("                         SD_rojo = %.6f  R2_rojo = %.6f"
              % (imp["sd_rojo"], imp["r2_rojo"]))
        cruzado = (abs(rf["sd"] - imp["sd_rojo"]) < 1e-3
                   and abs(rb["sd"] - imp["sd_azul"]) < 1e-3)
        derecho = (abs(rf["sd"] - imp["sd_azul"]) < 1e-3
                   and abs(rb["sd"] - imp["sd_rojo"]) < 1e-3)
        if cruzado:
            print("    -> SD y R2 estan INTERCAMBIADOS entre azul y rojo en el")
            print("       recuadro: el SD del ajuste azul es %.6f (impreso en rojo)"
                  % rf["sd"])
            print("       y el del rojo es %.6f (impreso en azul)." % rb["sd"])
        elif derecho:
            print("    -> el recuadro coincide con los ajustes.")
        else:
            print("    -> no coincide de ninguna de las dos formas; revisar la lectura")
        salida[det] = dict(fotopico=rf, borde=rb)
    print()
    print("  Consecuencia sobre lo que AFIRMA la Seccion III:")
    print('    "the standard deviation associated with the NCC-SA calibration')
    print('     consistently remains below that of the full energy peak fitting')
    print('     method, indicating higher precision"')
    for det, r in salida.items():
        print("    %-10s  SD(NCC-SA) = %.3f   SD(fotopico) = %.3f   -> %s"
              % (det, r["borde"]["sd"], r["fotopico"]["sd"],
                 "se cumple" if r["borde"]["sd"] < r["fotopico"]["sd"]
                 else "NO se cumple"))
    return salida


def rotulos():
    """Chequea que el nombre de la fuente de cada punto corresponda a su energia."""
    print()
    print("=" * 78)
    print("Fig. 3 - coherencia de los rotulos (fuente contra energia)")
    print("=" * 78)
    malos = 0
    for det, d in FIG3.items():
        for serie in ("fotopico", "borde"):
            for nombre, ch, en in d[serie]:
                esperado = FUENTE_DE_ENERGIA.get(int(en))
                if esperado and esperado != nombre:
                    print("  %-9s %-8s (%4d, %4d) rotulado %-6s -> deberia ser %s"
                          % (det, serie, ch, en, nombre, esperado))
                    malos += 1
    print("  %d rotulos cruzados." % malos)
    if malos:
        print("  No cambian ningun ajuste (los pares canal-energia estan bien),")
        print("  pero si se leen los rotulos se le atribuyen lineas a la fuente")
        print("  equivocada.")


def ec1():
    """La Ec. 1 impresa contra las energias que la propia Fig. 3 rotula."""
    print()
    print("=" * 78)
    print("Ec. 1 - la formula impresa contra los bordes que rotula el paper")
    print("=" * 78)
    print("  gamma [keV]   rotulo   Ec.1 impresa   Ec.1 corregida")
    for eg, rot in ((661.657, 477), (834.848, 639), (1173.228, 963),
                    (1274.537, 1061), (1332.492, 1117)):
        print("  %10.3f   %6d   %12.1f   %13.1f"
              % (eg, rot, ec.e_max_impresa(eg), ec.e_max(eg)))
    print("  -> los rotulos siguen la CORREGIDA. La Ec. 1 impresa da la energia")
    print("     del foton retrodispersado y esta mal por un factor ~2.6 a 5.2.")


def ec12():
    """La Ec. 12 impresa contra la separacion de los dos bordes del Co-60."""
    print()
    print("=" * 78)
    print("Ec. 12 - dE = 'energy difference between gamma lines'")
    print("=" * 78)
    d = FIG3["LaBr3:Ce"]
    borde = {int(e): c for _, c, e in d["borde"]}
    if 963 in borde and 1117 in borde:
        dB = borde[1117] - borde[963]
        rb = _linea("borde", d["borde"])
        D = 1.0 / rb["m"]                      # canal por keV
        dg = ec.delta_e_impresa(1173.228, 1332.492)
        db = ec.delta_e_bordes(1173.228, 1332.492)
        print("  LaBr3: los dos bordes del Co-60 caen en los canales %d y %d"
              % (borde[963], borde[1117]))
        print("         separacion medida            : %5.1f canales" % dB)
        print("         D*dE con dE de GAMMAS  (%.1f keV): %5.1f canales" % (dg, D * dg))
        print("         D*dE con dE de BORDES  (%.1f keV): %5.1f canales" % (db, D * db))
        print("  -> la separacion medida coincide con la version de BORDES.")
        print("     Con la formula impresa el segundo borde queda %.1f canales "
              "mas arriba," % (D * dg - dB))
        print("     o sea %.1f keV, contra un SD del ajuste de %.1f keV."
              % ((D * dg - dB) * rb["m"], rb["sd"]))
    print("  (El panel de NaI no sirve para este chequeo: no tiene el borde de")
    print("   1117 keV.)")


def coseno():
    """Las lecturas posibles de la similitud coseno que reporta el paper."""
    print()
    print("=" * 78)
    print("Similitud coseno - el paper reporta 99.9993% (LaBr3) y 99.9999% (NaI)")
    print("=" * 78)
    for det, d in FIG3.items():
        rf = _linea("fotopico", d["fotopico"])
        rb = _linea("borde", d["borde"])
        ch = np.array(sorted([c for _, c, _ in d["fotopico"]]
                             + [c for _, c, _ in d["borde"]]), float)
        s = R.similitudes(rf, rb, ch)
        et_f = np.array([e for _, _, e in d["fotopico"]], float)
        et_b = np.array([e for _, _, e in d["borde"]], float)
        s["recon_vs_teorico_fotopico"] = R.similitud_coseno(
            R.energia_reconstruida(rf, [c for _, c, _ in d["fotopico"]]), et_f)
        s["recon_vs_teorico_borde"] = R.similitud_coseno(
            R.energia_reconstruida(rb, [c for _, c, _ in d["borde"]]), et_b)
        print("  %-10s (objetivo %.6f)" % (det, d["coseno"]))
        for k, v in s.items():
            print("      %-28s %.8f" % (k, v))
    print("  -> ninguna lectura reproduce los DOS valores impresos. El paper no")
    print("     define entre que vectores calcula la similitud, asi que este")
    print("     numero no es verificable.")
    print("  Nota aparte: entre dos rectas casi proporcionales esta metrica es")
    print("     insensible por construccion, asi que 99.999% no acota gran cosa.")


def fig5():
    """La robustez frente al intervalo, resumida con numeros."""
    print()
    print("=" * 78)
    print("Fig. 5 - robustez frente al intervalo de ajuste")
    print("=" * 78)
    v = np.array([b for _, b in FIG5], float)
    print("  intervalo          borde")
    for (lo, hi), b in FIG5:
        print("  [%4d, %4d]        %4d" % (lo, hi, b))
    print("  media %.1f   SD %.2f   rango %.0f canales   dispersion %.2f%%"
          % (v.mean(), v.std(ddof=0), v.max() - v.min(),
             100 * (v.max() - v.min()) / v.mean()))
    print("  -> el paper lo describe como 'highly consistent' y no da un numero;")
    print("     el que sale de sus propios paneles es 1.2% de dispersion total.")
    print("     El intervalo mas angosto ([1500,1600], que no alcanza a contener")
    print("     el flanco entero) es el unico que se aparta: sacandolo, la")
    print("     dispersion baja a %.2f%%." % (100 * (v[1:].max() - v[1:].min())
                                              / v[1:].mean()))


def fig4():
    """Que se puede y que no se puede auditar de la Fig. 4."""
    print()
    print("=" * 78)
    print("Fig. 4 - lo que el texto afirma, y por que no es auditable")
    print("=" * 78)
    for det, d in REDUCCIONES.items():
        print("  %-7s reduccion del error medio: %d%% contra el metodo de la"
              % (det, d["diferencia"]))
        print("          diferencia y %d%% contra el gaussiano." % d["gauss"])
    print("  La Fig. 4 es el unico soporte de esos numeros y no publica la tabla")
    print("  subyacente: los errores relativos hay que leerlos del grafico, con")
    print("  una precision (~0.1%) del mismo orden que los propios valores del")
    print("  NCC-SA (0.2-0.6%). Digitalizarlos daria una cifra inventada, asi")
    print("  que este bloque se deja SIN auditar a proposito.")
    print("  Lo que si se puede decir: la Fig. 4 tiene 6 puntos por metodo, o sea")
    print("  6 bordes Compton (341, 477, 639, 963, 1061 y 1117 keV), que son los")
    print("  de las cuatro fuentes. El panel de LaBr3 usa 5 y el de NaI 4, asi")
    print("  que el conjunto de bordes NO es el mismo en las dos figuras.")


def inconsistencias():
    """Discrepancias internas del texto, sin interpretacion."""
    print()
    print("=" * 78)
    print("Otras inconsistencias internas del paper")
    print("=" * 78)
    print("  1. Ec. 14 define la loss como minimos cuadrados; el parrafo de")
    print("     Metropolis dice que E(x) es 'the NCC-based matching error'. Son")
    print("     objetivos distintos y hay que elegir uno.")
    print("  2. Seccion III dice que la robustez se probo 'calibrating EJ309")
    print("     spectra (137Cs source)'; el epigrafe de la Fig. 5 dice BC408.")
    print("  3. El borde del Co-60 de 1332.5 keV se rotula 1117 keV en la Fig. 3;")
    print("     la Ec. 1 corregida da 1118.1 keV.")
    print("  4. T_init = 2000 y |dE_mean| < 0.001 son constantes con unidades del")
    print("     COSTO, y el costo no esta fijado (punto 1). No son transferibles.")


def todo():
    rutas.salida_utf8()
    fig3()
    rotulos()
    ec1()
    ec12()
    coseno()
    fig5()
    fig4()
    inconsistencias()


if __name__ == "__main__":
    todo()
