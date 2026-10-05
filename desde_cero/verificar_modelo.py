"""
Comprobar la Etapa 5: contra replica_paper y con verdades conocidas.

    python desde_cero/verificar_modelo.py

1. ORACULO. `replica_paper.ecuaciones.modelo` (un borde) y
   `modelo_dos_bordes` (Co-60, con la Ec. 12 corregida). Las convenciones
   son distintas y se convierten, no se esconden:
     - el paper: S = A_r * dsigma/dE(ch * E_max / B), A_r sin significado;
     - nosotros: S = A * forma(ch/D) * paso/D, A = cuentas, forma de integral 1.
   O sea A_r = A * w_i * paso / (D * sum_j w_j sigma_j), B = D * E_max y
   sigma = sigma_rel * D * E_ref. Con las mismas constantes (m_e c^2, r0)
   tienen que dar lo mismo a redondeo de maquina.

2. VERDAD CONOCIDA.
   a. Conservacion: la convolucion no crea ni destruye cuentas. Sumando en
      una grilla que incluya los canales negativos (la resolucion lleva
      cuentas ahi), el total es el mismo para cualquier sigma, e igual a A
      salvo el error de sumar la forma en los centros de bin.
   b. La grilla extendida importa: con la convolucion solo en la ventana el
      modelo sale bajo en el borde izquierdo de la ventana (§7.14).
   c. El 50 % de altura se corre por encima del borde verdadero, cada vez
      mas con sigma; se compara con la tabla de §7.44.
   d. Co-60: el segundo borde cae en D * E_max2, o sea a D * 154.7 keV del
      primero, sin haberlo escrito en ningun lado.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from replica_paper import ecuaciones as oraculo             # noqa: E402
from desde_cero import fisica, modelo                       # noqa: E402

K = dict(mec2=fisica.M_E, r0=fisica.R0)


def rel_max(a, b, piso=1e-12):
    m = np.abs(b) > piso * np.abs(b).max()
    return float(np.max(np.abs(a[m] - b[m]) / np.abs(b[m])))


def main():
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-50s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    A, paso = 1e5, 2.5
    print("=" * 78)
    print("1. Contra replica_paper/ecuaciones.py")
    print("=" * 78)
    for D, s in ((1.9, 0.03), (1.9, 0.10), (1.25, 0.08)):
        # un borde: Cs-137
        eg = 661.657
        B = D * float(fisica.e_max(eg))
        ch = np.arange(0.5 * paso, 1.6 * B, paso)
        # por_bin=False: el oraculo evalua la forma en el centro del bin
        nos = modelo.modelo(ch, A, D, s, "Cs-137", por_bin=False)
        A_r = A * paso / (D * fisica.sigma_total_cerrada(eg))
        ora = oraculo.modelo(ch, eg, A_r, B, s * B, **K)
        check("Cs-137, D = %.2f, sigma_rel = %.2f" % (D, s), rel_max(nos, ora) < 1e-10,
              "max rel %.1e" % rel_max(nos, ora))

        # dos bordes: Co-60, un solo sigma para los dos (ver modelo.py)
        (e1, e2), (w1, w2) = modelo.lineas_de("Co-60")
        tot = w1 * fisica.sigma_total_cerrada(e1) + w2 * fisica.sigma_total_cerrada(e2)
        B1 = D * float(fisica.e_max(e1))
        ch = np.arange(0.5 * paso, 1.2 * D * float(fisica.e_max(e2)) + 5 * s * B1, paso)
        nos = modelo.modelo(ch, A, D, s, "Co-60", por_bin=False)
        ora = oraculo.modelo_dos_bordes(ch, e1, e2, A * w1 * paso / (D * tot),
                                        A * w2 * paso / (D * tot), B1, s * B1,
                                        s * B1, D, delta_e="bordes", **K)
        check("Co-60,  D = %.2f, sigma_rel = %.2f" % (D, s), rel_max(nos, ora) < 1e-10,
              "max rel %.1e" % rel_max(nos, ora))

    print()
    print("=" * 78)
    print("2. Verdades conocidas")
    print("=" * 78)
    D = 1.9
    B = D * float(fisica.e_max(661.657))
    # a. conservacion. La grilla arranca bien por debajo de cero: ahi caen
    # las cuentas que la resolucion saca de la zona de energia ~0.
    ch = np.arange(-600.0 + 0.5 * paso, 1800.0, paso)
    E, w = modelo.lineas_de("Cs-137")
    crudo = A * modelo.forma_energia(ch / D, E, w) * paso / D   # sin resolucion
    # misma forma (en centros de bin) a los dos lados: lo que se prueba es la
    # convolucion, no la manera de llenar los bines
    sumas = [modelo.modelo(ch, A, D, s, "Cs-137", por_bin=False).sum()
             for s in (0.02, 0.08, 0.15)]
    d = max(abs(x / crudo.sum() - 1) for x in sumas)
    check("a. suma(modelo) = suma(sin resolucion), 3 sigmas", d < 1e-12,
          "max rel %.1e" % d)
    check("a. suma(en centros de bin) = A", abs(crudo.sum() / A - 1) < 5e-3,
          "%.1f contra %.0f (%+.2f %%: forma sumada en centros de bin)"
          % (crudo.sum(), A, 100 * (crudo.sum() / A - 1)))
    total = modelo.modelo(ch, A, D, 0.08, "Cs-137").sum()
    check("a. suma(integrado por bin) = A", abs(total / A - 1) < 1e-9,
          "%.4f contra %.0f" % (total, A))

    # a'. continuidad en D: el modelo integrado por bin no salta cuando el
    # borde cruza un centro de bin; el evaluado en el centro, si.
    chv = np.arange(0.5 * paso, 1800.0, paso)
    Ds = np.linspace(1.9, 1.9 + 2 * paso / 477.0, 400)     # el borde cruza 2 bines
    salto = lambda pb: max(np.abs(modelo.modelo(chv, A, b, 0.08, "Cs-137", por_bin=pb)
                                  - modelo.modelo(chv, A, a, 0.08, "Cs-137", por_bin=pb)).max()
                           for a, b in zip(Ds[:-1], Ds[1:]))   # noqa: E731
    s_bin, s_cen = salto(True), salto(False)
    check("a'. continuidad en D (mayor cambio entre D vecinos)", s_bin < 0.2 * s_cen,
          "integrado %.3f cuentas, en el centro %.3f" % (s_bin, s_cen))
    ch_pos = np.arange(0.5 * paso, 1800.0, paso)
    perdido = 1 - modelo.modelo(ch_pos, A, D, 0.15, "Cs-137").sum() / crudo.sum()
    print("  %-50s %.1f %% de las cuentas (sigma_rel 0.15) caen en canal < 0"
          % ("   (nota)", 100 * perdido))

    # b. la grilla extendida
    ch = np.arange(0.45 * B, 1.3 * B, paso)
    con = modelo.modelo(ch, A, D, 0.08, "Cs-137")
    sin = modelo.modelo(ch, A, D, 0.08, "Cs-137", extender=False)
    falta = 1 - sin[0] / con[0]
    check("b. sin grilla extendida, borde izquierdo de la ventana", falta > 0.3,
          "sale %.0f %% bajo (el ajuste lo pagaria achicando sigma)" % (100 * falta))
    medio = np.abs(sin / con - 1)[ch.size // 2]
    check("b.   ... y en el medio de la ventana, igual", medio < 1e-6,
          "diferencia %.1e" % medio)

    # c. el 50 % de altura (sin ruido) contra la columna "sin ruido" de §7.44
    # no se puede comparar directo (aquella es del metodo de la derivada); se
    # exige que sea positivo, crezca con sigma, y quede en el orden de §7.44
    # (+0.8 / +2.2 / +3.8 / +5.4 / +8.5 % a 0.03 / .06 / .09 / .12 / .18, con
    # ruido).
    ch = np.arange(0.5 * paso, 3000.0, paso)
    ss = (0.03, 0.06, 0.09, 0.12, 0.18)
    sesgo = [modelo.mitad_de_altura(ch, modelo.modelo(ch, A, D, s, "Cs-137"), 0.5 * B)
             / B - 1 for s in ss]
    crece = all(b > a for a, b in zip(sesgo, sesgo[1:])) and sesgo[0] > 0
    check("c. 50 % de altura: sesgo positivo y creciente", crece,
          " ".join("%.2f:%+.1f%%" % (s, 100 * b) for s, b in zip(ss, sesgo)))

    # d. el segundo borde del Co-60. Con sigma muy chico el modelo es el
    # escalon: el ultimo canal con cuentas de la linea 2 esta en D*E_max2.
    paso_f = 0.1
    ch = np.arange(0.5 * paso_f, 2500.0, paso_f)
    y = modelo.modelo(ch, A, D, 1e-4, "Co-60")
    (e1, e2), _ = modelo.lineas_de("Co-60")
    b2 = D * float(fisica.e_max(e2))
    ultimo = ch[np.flatnonzero(y > 1e-6 * y.max())[-1]]
    # el kernel, por chico que sea sigma, tiene radio 4 sigma: el ultimo canal
    # con cuentas puede estar hasta ahi (mas un paso) por encima del escalon
    tol = modelo.N_SIGMAS * 1e-4 * D * float(fisica.e_max(e1)) + paso_f
    check("d. Co-60: fin del espectro en D*E_max2", 0 <= ultimo - b2 < tol,
          "%.2f contra %.2f mV*ns (D*154.7 = %.1f entre bordes)"
          % (ultimo, b2, D * (fisica.e_max(e2) - fisica.e_max(e1))))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
