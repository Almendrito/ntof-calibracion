"""
Comprobar la Etapa 4: contra replica_paper y contra verdades conocidas.

    python desde_cero/verificar_fisica.py

1. ORACULO. `replica_paper/ecuaciones.py` implementa las mismas ecuaciones
   con las constantes del paper (511 keV, 2.818e-13 cm). Se le pasan las
   NUESTRAS (mec2=, r0=) para comparar la formula y no las constantes: tiene
   que dar lo mismo a redondeo de maquina.

2. VERDAD CONOCIDA, sin oraculo:
   a. los bordes del PLAN: 477.3, 963.4 y 1118.1 keV;
   b. la Ec. 10 escrita (`dsigma_de`) contra la Ec. 2 pasada a energia por el
      jacobiano (`dsigma_de_jacobiano`): el jacobiano del paper, verificado;
   c. la seccion eficaz total por tres caminos independientes (integral en
      angulo de la Ec. 2, integral en energia de la Ec. 10, formula cerrada)
      a 6 decimales para 500, 662, 1000 y 2000 keV;
   d. el limite de Thomson: sigma -> 8 pi r0^2/3 = 0.665 b a baja energia, con
      la correccion de primer orden conocida, sigma ~ sigma_T (1 - 2 alpha);
   e. los valores que se leen en la Fig. 1 del paper (con sus constantes).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from replica_paper import ecuaciones as oraculo             # noqa: E402
from desde_cero import fisica as f                          # noqa: E402


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.abs(b) > 0
    return float(np.max(np.abs(a[m] - b[m]) / np.abs(b[m])))


def main():
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-52s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    print("=" * 78)
    print("1. Contra replica_paper/ecuaciones.py (con nuestras constantes)")
    print("=" * 78)
    kw = dict(mec2=f.M_E, r0=f.R0)
    eg = np.array([200.0, 500.0, 661.657, 1173.228, 1332.492, 2000.0, 5000.0])
    check("Ec. 1 corregida", rel(f.e_max(eg), oraculo.e_max(eg, f.M_E)) < 1e-14,
          "max rel %.1e" % rel(f.e_max(eg), oraculo.e_max(eg, f.M_E)))
    check("Ec. 1 impresa",
          rel(f.e_max_impresa(eg), oraculo.e_max_impresa(eg, f.M_E)) < 1e-14,
          "max rel %.1e" % rel(f.e_max_impresa(eg), oraculo.e_max_impresa(eg, f.M_E)))
    th = np.linspace(0.0, np.pi, 2001)
    d2 = max(rel(f.dsigma_domega(th, e, **kw), oraculo.dsigma_domega(th, e, **kw))
             for e in eg)
    check("Ec. 2 (en angulo)", d2 < 1e-12, "max rel %.1e" % d2)
    d10 = max(rel(f.dsigma_de(np.linspace(0, 1.2 * f.e_max(e), 3001), e, **kw),
                  oraculo.dsigma_de(np.linspace(0, 1.2 * f.e_max(e), 3001), e,
                                    A=1.0, **kw)) for e in eg)
    check("Ec. 10 (en energia del electron)", d10 < 1e-12, "max rel %.1e" % d10)
    d12 = abs((f.e_max(1332.492) - f.e_max(1173.228))
              - oraculo.delta_e_bordes(1173.228, 1332.492, f.M_E))
    check("Ec. 12 corregida (entre bordes)", d12 < 1e-10, "|d| %.1e keV" % d12)

    print()
    print("=" * 78)
    print("2. Verdades conocidas")
    print("=" * 78)
    for fuente, lineas in f.FUENTES.items():
        for (e, _), esperado in zip(lineas, {"Cs-137": (477.3,),
                                             "Co-60": (963.4, 1118.1)}[fuente]):
            check("a. borde de %s %.3f keV" % (fuente, e),
                  round(float(f.e_max(e)), 1) == esperado,
                  "%.2f keV (PLAN: %.1f)" % (f.e_max(e), esperado))

    dj = max(rel(f.dsigma_de_jacobiano(np.linspace(1e-3, f.e_max(e), 5001), e),
                 f.dsigma_de(np.linspace(1e-3, f.e_max(e), 5001), e)) for e in eg)
    check("b. Ec. 10 == jacobiano(Ec. 2)", dj < 1e-9, "max rel %.1e" % dj)

    for e in (500.0, 662.0, 1000.0, 2000.0):
        sa = f.sigma_total_angulo(e) / f.BARN
        se = f.sigma_total_energia(e) / f.BARN
        sc = f.sigma_total_cerrada(e) / f.BARN
        igual = round(sa, 6) == round(se, 6) == round(sc, 6)
        check("c. sigma total a %4.0f keV [b]" % e, igual,
              "%.6f / %.6f / %.6f" % (sa, se, sc))

    st = f.sigma_thomson() / f.BARN
    for e in (0.01, 0.1, 1.0):
        a = e / f.M_E
        s = f.sigma_total_angulo(e) / f.BARN
        # a primer orden en alpha, sigma = sigma_T (1 - 2 alpha); lo que queda
        # tiene que ser de segundo orden (el siguiente termino es 26/5 alpha^2)
        resto = abs(s / st - (1 - 2 * a))
        check("d. Thomson, E = %5.2f keV (alpha = %.1e)" % (e, a),
              resto < 6 * a ** 2 + 1e-12,
              "sigma/sigma_T = %.8f, resto tras 1er orden %.1e" % (s / st, resto))
    # la formula cerrada a energia muy baja pierde precision por cancelacion:
    # se muestra, no se exige (para eso esta el camino por angulo)
    print("  %-52s cerrada / por angulo - 1: %.0e a 10 eV, %.0e a 1 eV, "
          "%.0e a 0.01 eV" % ("   (nota: cancelacion en la cerrada)",
                              *[f.sigma_total_cerrada(e) / f.sigma_total_angulo(e) - 1
                                for e in (1e-2, 1e-3, 1e-5)]))

    # e. la Fig. 1 del paper, leida del eje: valor en E -> 0 y en el borde
    fig1 = {500: (1.02, 1.68), 662: (0.58, 1.12), 1000: (0.255, 0.645),
            2000: (0.064, 0.283), 5000: (0.010, 0.11)}
    for e0, (y0, ypico) in fig1.items():
        E = np.linspace(1e-6, float(f.e_max(e0, 511.0)), 20001)
        y = f.dsigma_de(E, e0, mec2=511.0, r0=2.818e-13) / 1e-27
        # lo leido del grafico tiene 2-3 cifras: tolerancia del 3 % (5 % en
        # el pico de 5000, que se lee "0.11")
        tol = 0.05 if e0 == 5000 else 0.03
        bien = abs(y[0] / y0 - 1) < tol and abs(y.max() / ypico - 1) < tol
        check("e. Fig. 1, E0 = %4d keV (x1e-27 cm^2/keV)" % e0, bien,
              "E->0 %.3f (fig %.3f)  borde %.3f (fig %.3f)"
              % (y[0], y0, y.max(), ypico))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
