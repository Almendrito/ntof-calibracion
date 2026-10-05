"""
Comprobar la Etapa 6: con B conocido y contra los dos oraculos.

    python desde_cero/verificar_ncc.py          (~5 min: arma los espectros reales)

1. VERDAD CONOCIDA. Espectros fabricados con `modelo.modelo` (D y sigma
   conocidos, Poisson, la grilla y el numero de cuentas de la corrida real):
   el B elegido tiene que caer a menos de un bin del verdadero cuando la
   plantilla tiene el sigma correcto. Y se mide cuanto se corre con un sigma
   equivocado, en los dos barridos: el paper no dice con que sigma construir
   la plantilla, asi que eso es parte del metodo.

2. `replica_paper/ncc.py` (el barrido del paper) contra `barrido_libre`: la
   misma cuenta, tienen que dar la MISMA ventana y el mismo rho. Solo con
   Cs-137: replica_paper arma la plantilla con UNA linea, y para el Co-60 la
   nuestra usa las dos, asi que ahi no son la misma cuenta.

3. `analisis/ncc.py:buscar_intervalo` contra `barrido_B`, con la MISMA
   plantilla (Klein-Nishina de `modelo.forma_energia`, pasada como tabla, y
   `por_bin=False`: evaluada en el centro del bin, como el oraculo).
   No es bit a bit: el oraculo interpola la plantilla de una tabla y usa un
   kernel de 6 sigma contra 4 del nuestro. Se exige el mismo B (o el vecino
   de la grilla del barrido) y rho a 1e-4.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import ncc as oraculo                          # noqa: E402
from replica_paper import ncc as oraculo_paper               # noqa: E402
from desde_cero import espectro, fisica, lector, modelo, ncc  # noqa: E402

FUENTE = {"cs137_1600": "Cs-137", "co60": "Co-60"}
D_VERDAD = {"cs137_1600": 1.40, "co60": 1.00}   # cerca de lo que dan los datos


def main():
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-44s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    espectros = {}
    for clave in FUENTE:
        if os.path.exists(lector.CORRIDAS_TRABAJO[clave][0]):
            e, _, _ = espectro.de_corrida(clave)
            espectros[clave] = (e["canales"], e["neto"])

    print("=" * 78)
    print("1. Espectros fabricados con B conocido (sigma verdadero 0.09, 10 "
          "realizaciones)")
    print("=" * 78)
    rng = np.random.default_rng(5)
    for clave, (x, y_real) in espectros.items():
        f, D = FUENTE[clave], D_VERDAD[clave]
        B = D * float(ncc.bordes_de(f)[0])
        paso = x[1] - x[0]
        mu = modelo.modelo(x, y_real.sum(), D, 0.09, f)
        for s in (0.06, 0.09, 0.12):
            nos, pap = [], []
            for r in range(10):
                y = rng.poisson(mu).astype(float)
                nos.append((ncc.barrido_B(x, y, f, sigma_rel=s)["B"] - B) / paso)
                if r < 3:
                    pap.append((ncc.barrido_libre(x, y, f, sigma_rel=s)["B"] - B) / paso)
            detalle = ("nuestro %+.2f +- %.2f bines   paper %+.2f bines"
                       % (np.mean(nos), np.std(nos), np.mean(pap)))
            if s == 0.09:
                check("%s, plantilla sigma = %.2f (la verdadera)" % (f, s),
                      abs(np.mean(nos)) < 1.0, detalle)
            else:
                print("  %-44s %s" % ("%s, plantilla sigma = %.2f" % (f, s), detalle))

    print()
    print("=" * 78)
    print("2. barrido_libre contra replica_paper/ncc.py (Cs-137)")
    print("=" * 78)
    casos = []
    if "cs137_1600" in espectros:
        casos.append(("real", *espectros["cs137_1600"]))
        x = espectros["cs137_1600"][0]
        casos.append(("fabricado", x, rng.poisson(
            modelo.modelo(x, 1e5, 1.4, 0.09, "Cs-137")).astype(float)))
    for nombre, x, y in casos:
        a = ncc.barrido_libre(x, y, "Cs-137", sigma_rel=0.06)
        b = oraculo_paper.buscar_intervalo(x, y, 661.657, sigma_rel=0.06,
                                           mec2=fisica.M_E)
        igual = a["lo"] == b["lo"] and a["hi"] == b["hi"] and abs(a["rho"] - b["rho"]) < 1e-12
        check("Cs-137 %s" % nombre, igual,
              "[%.1f, %.1f] rho %.6f / [%.1f, %.1f] rho %.6f"
              % (a["lo"], a["hi"], a["rho"], b["lo"], b["hi"], b["rho"]))

    print()
    print("=" * 78)
    print("3. barrido_B contra analisis/ncc.py (misma plantilla Klein-Nishina)")
    print("=" * 78)
    for clave, (x, y) in espectros.items():
        f = FUENTE[clave]
        E, w = modelo.lineas_de(f)
        bordes = ncc.bordes_de(f)
        cen = np.arange(0.0, 1.3 * bordes[-1], 0.02)
        # por_bin=False: el oraculo evalua la plantilla en el centro del bin
        a = ncc.barrido_B(x, y, f, por_bin=False)
        b = oraculo.buscar_intervalo(x, y, cen, modelo.forma_energia(cen, E, w),
                                     bordes, sigma_rel=ncc.SIGMA_REL)
        dB = abs(a["B"] - b["B"])
        paso_B = a["Bs"][1] - a["Bs"][0]
        bien = dB <= paso_B * 1.0001 and abs(a["rho"] - b["rho"]) < 1e-4
        check(clave, bien, "B %.2f / %.2f (paso %.2f)  rho %.6f / %.6f"
              % (a["B"], b["B"], paso_B, a["rho"], b["rho"]))
        c = ncc.barrido_B(x, y, f)
        print("  %-44s con la plantilla integrada por bin: B %.2f (%+.2f %%)"
              % ("", c["B"], 100 * (c["B"] / a["B"] - 1)))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
