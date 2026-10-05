"""
Comprobar la Etapa 10.

    python desde_cero/verificar_clasicos.py        (~5 min: arma los espectros reales)

1. CONTRA `codigo/etapas/clasicos.py`: 50 % de altura, diferencia y escalon
   gaussiano, con la MISMA region, sobre espectros sinteticos y reales:
   iguales a redondeo de maquina (la ventana de suavizado automatica incluida).

2. CONTRA `replica_paper/convencionales.py`: la gaussiana del paper
   (`metodo_gauss`, fraccion 0.5), misma region: igual.

3. SESGO CON VERDAD CONOCIDA, contra la tabla de §7.44 (informativo: el
   espectro sintetico de §7.44 no es el mismo, B = 600 canales y otro modelo
   directo). Se exige solo lo que es propiedad del METODO: sesgo positivo y
   creciente con sigma para 50 % de altura y escalon.

4. DATOS REALES: D de cada metodo contra la del ajuste con Geant4 (Etapa 8).
   En el Co-60 se reporta contra los DOS bordes: los metodos clasicos
   detectan la caida del conjunto, que la domina el segundo borde (1118 keV),
   y dividir por el primero (963 keV, como hace el informe) mete un +16 %
   que es solo 1118/963.
"""

import json
import os
import sys

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)
sys.path.insert(0, os.path.join(RAIZ, "codigo"))       # codigo/ importa por paquete

from etapas import clasicos as oraculo                      # noqa: E402
from replica_paper import convencionales as oraculo_paper   # noqa: E402
from desde_cero import clasicos, espectro, fisica, lector, lote, modelo, ncc  # noqa: E402

FUENTE = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
TABLA_744 = {  # §7.44, con ruido, sigma/B = 0.03, 0.06, 0.09, 0.12, 0.18
    "50 % de altura": (0.8, 2.2, 3.8, 5.4, 8.5),
    "diferencia": (0.7, 1.5, 2.8, 4.4, 7.4),
    "escalon gaussiano": (2.6, 4.7, 6.3, 7.5, 9.2),
}


def del_oraculo(x, y, region, e):
    """Los tres de codigo/ con la misma region y ventana automatica."""
    v = oraculo._ventana_auto(x, y, region, 3)
    s = oraculo._suavizar(y, v, 3)
    m = oraculo._mascara(x, region)
    from etapas.model import punto_mitad_altura
    return {"50 % de altura": float(punto_mitad_altura(x[m], s[m])),
            "diferencia": oraculo.metodo_diferencia(x, y, region=region)["borde"],
            "escalon gaussiano": oraculo.metodo_fwhm_gauss(x, y, region=region,
                                                            errores=e)["borde"],
            "gaussiana del paper": oraculo_paper.metodo_gauss(x, y, region)["borde"]}


def main():
    ok = True

    def check(nombre, bien, detalle):
        nonlocal ok
        ok &= bool(bien)
        print("  %-34s %s   %s" % (nombre, detalle, "ok" if bien else "FALLA"))

    casos = []
    rng = np.random.default_rng(9)
    B = float(fisica.e_max(661.657))
    x = np.arange(1.25, 2.2 * B, 2.5)
    for s in (0.05, 0.12):
        y = rng.poisson(modelo.modelo(x, 2e5, 1.0, s, "Cs-137")).astype(float)
        casos.append(("sintetico sigma %.2f" % s, x, y, np.sqrt(np.maximum(y, 1)), B, None))
    reales = {}
    for clave, f in FUENTE.items():
        if os.path.exists(lector.CORRIDAS_TRABAJO[clave][0]):
            e, _, _ = espectro.de_corrida(clave)
            B_ncc = ncc.barrido_B(e["canales"], e["neto"], f)["B"]
            reales[clave] = (e, B_ncc)
            casos.append((clave, e["canales"], e["neto"], e["error"], B_ncc, f))

    print("=" * 78)
    print("1 y 2. Contra codigo/etapas/clasicos.py y replica_paper/convencionales.py")
    print("=" * 78)
    for nombre, x, y, err, Bc, _ in casos:
        region = clasicos.region_anclada(Bc, x)
        nos = clasicos.los_cuatro(x, y, region, err)
        ora = del_oraculo(x, y, region, err)
        for m in clasicos.METODOS:
            a, b = nos[m], ora[m]
            igual = (np.isnan(a) and np.isnan(b)) or abs(a - b) <= 1e-6 * abs(b)
            check("%s: %s" % (nombre, m), igual, "%.4f / %.4f" % (a, b))

    print()
    print("=" * 78)
    print("3. Sesgo con verdad conocida (Klein-Nishina, B conocido, 10 realizaciones)")
    print("=" * 78)
    sig = (0.03, 0.06, 0.09, 0.12, 0.18)
    con = clasicos.sesgo_sintetico(sig)
    for m in clasicos.METODOS:
        v = 100 * np.array(con[m])
        ref = TABLA_744.get(m)
        print("  %-20s %s%s" % (m, " ".join("%+5.1f" % a for a in v),
                                "   (§7.44: %s)" % " ".join("%+5.1f" % a for a in ref)
                                if ref else ""))
    for m in ("50 % de altura", "escalon gaussiano"):
        v = np.array(con[m])
        check("%s: positivo y creciente" % m, v[0] > 0 and np.all(np.diff(v) > 0),
              " ".join("%+.1f" % (100 * a) for a in v))

    print()
    print("=" * 78)
    print("4. Datos reales contra el ajuste con Geant4")
    print("=" * 78)
    lot = json.load(open(lote.SALIDA, encoding="utf-8"))
    for clave, (e, B_ncc) in reales.items():
        f = FUENTE[clave]
        bordes = ncc.bordes_de(f)
        D = lot[clave]["D"]
        r = clasicos.los_cuatro(e["canales"], e["neto"],
                                clasicos.region_anclada(B_ncc, e["canales"]), e["error"])
        for m, v in r.items():
            txt = "contra %.1f keV: %+.1f %%" % (bordes[0], 100 * (v / bordes[0] / D - 1))
            if len(bordes) > 1:
                txt += "   contra %.1f keV: %+.1f %%" % (bordes[-1],
                                                      100 * (v / bordes[-1] / D - 1))
            print("  %-12s %-20s %s" % (clave, m, txt))

    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
