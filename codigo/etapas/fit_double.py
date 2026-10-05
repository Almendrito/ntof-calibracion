"""
fit_double.py - Etapa 6: modelo de DOS bordes solapados (Ec. 11-13).

Por que hace falta
------------------
El Co-60 emite dos gammas (1173.2 y 1332.5 keV) -> DOS bordes Compton, en
963.4 y 1118.1 keV, separados apenas un 16%. Con el broadening del detector se
funden en una sola bajada. Ajustar UN borde a esa bajada da un B sesgado, que
cae entre ambos, y ademas infla sigma para poder cubrir los dos (CLAUDE.md 7.5
y 7.13: sobre los datos reales el ajuste de un borde pedia sigma/B = 0.185, o
sea 43% FWHM, absurdo para plastico+SiPM).

El modelo (Ec. 11-13)
---------------------
    Theta = {A1, A2, B1, sigma1, sigma2, D}                        (Ec. 11)
    CE_ch,2 = B1 + D * dE                                          (Ec. 12)
    S_comp(ch) = f1(ch; A1, B1, sigma1) + f2(ch; A2, CE_ch,2, sigma2)  (Ec. 13)

donde cada f_i es el forward model de la Etapa 2 (Klein-Nishina de ESA linea
gamma, convolucionada con SU gaussiana), y D es el coeficiente de conversion
energia->canal: la PENDIENTE DE CALIBRACION que buscamos.

Dos correcciones al paper
-------------------------
1. dE: el paper dice "energy difference between gamma lines" (1332.5 - 1173.2 =
   159.3 keV). Esta mal. D convierte energia del ELECTRON a canal, y lo que se
   separa en canal son los BORDES, asi que dE = E_max,2 - E_max,1 = 154.7 keV.
   Con el dE equivocado, D sale ~3% corrido.
2. Con B1 y D ambos libres el mapa canal<->energia deja de pasar por el origen
   (queda afin: ch = D*E + c, con c = B1 - D*E_max,1). Eso contradice la Ec. 8,
   que si pasa por el origen. NO lo amarramos: dejarlo libre convierte la
   comparacion D vs B1/E_max,1 en un TEST de la hipotesis del origen. Si se
   quiere amarrar, esta la opcion ligar_origen=True.

Truco numerico [fuera del paper]: VARPRO
----------------------------------------
S_comp es LINEAL en A1 y A2. Entonces, fijados {B1, sigma1, sigma2, D}, las
amplitudes optimas se resuelven EXACTO con minimos cuadrados no negativos
(nnls), sin iterar. Eso baja la busqueda no lineal de 6 a 4 parametros y hace
barato un barrido grueso de semillas, que es lo que le falta al ajuste local
(CLAUDE.md 7.4). El SA de la Etapa 5 despues ataca lo mismo por la via del paper.

Corre solo:  python fit_double.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/fit_double.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas


import numpy as np
from scipy.optimize import least_squares, nnls
from etapas.model import forward_model
from etapas.physics import compton_edge_energy

# lineas del Co-60 y sus bordes Compton
E_GAMMA_1, E_GAMMA_2 = 1173.2, 1332.5


def delta_E(E_gamma1=E_GAMMA_1, E_gamma2=E_GAMMA_2):
    """dE de la Ec. 12: diferencia de los BORDES, no de las lineas gamma."""
    return compton_edge_energy(E_gamma2) - compton_edge_energy(E_gamma1)


# ----------------------------------------------------------------------
# 1) el modelo compuesto (Ec. 13)
# ----------------------------------------------------------------------
def componentes(canales, B1, sigma1, sigma2, D,
                E_gamma1=E_GAMMA_1, E_gamma2=E_GAMMA_2):
    """Las dos formas con amplitud 1 (lo que multiplica a A1 y A2)."""
    ce2 = B1 + D * delta_E(E_gamma1, E_gamma2)          # Ec. 12
    f1 = forward_model(canales, E_gamma1, 1.0, max(B1, 1e-9), sigma1)
    f2 = forward_model(canales, E_gamma2, 1.0, max(ce2, 1e-9), sigma2)
    return f1, f2, ce2


def modelo_dos_bordes(canales, A1, A2, B1, sigma1, sigma2, D,
                      E_gamma1=E_GAMMA_1, E_gamma2=E_GAMMA_2):
    """S_comp(ch) de la Ec. 13."""
    f1, f2, _ = componentes(canales, B1, sigma1, sigma2, D, E_gamma1, E_gamma2)
    return A1 * f1 + A2 * f2


def _amplitudes_optimas(f1, f2, y):
    """A1, A2 exactas por minimos cuadrados NO NEGATIVOS (VARPRO)."""
    M = np.vstack([f1, f2]).T
    try:
        coef, _ = nnls(M, y)
    except Exception:
        coef = np.zeros(2)
    return float(coef[0]), float(coef[1])


# ----------------------------------------------------------------------
# 2) semilla por barrido grueso (global de pobre)
# ----------------------------------------------------------------------
def semilla_por_barrido(canales, cuentas, intervalo, B1_rango, sigma_rel_rango=(0.03, 0.15),
                        n_B=60, n_s=20, E_gamma1=E_GAMMA_1, E_gamma2=E_GAMMA_2,
                        ligar_sigmas=True):
    """
    Barre (B1, sigma_rel) con D amarrado al origen (D = B1/E_max,1) y resuelve
    A1,A2 exacto en cada punto. Devuelve la mejor semilla [A1,A2,B1,s1,s2,D].

    ligar_sigmas: sigma2 = sigma1 * sqrt(E_max2/E_max1). La resolucion de un
    centellador escala ~sqrt(E) (estadistica de fotoelectrones), asi que amarrar
    asi quita una degeneracion del barrido sin mentirle a la fisica. En el ajuste
    final sigma1 y sigma2 quedan libres (como en el paper).
    """
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]
    Emax1 = compton_edge_energy(E_gamma1)
    razon_s = np.sqrt(compton_edge_energy(E_gamma2) / Emax1)

    mejor = (np.inf, None)
    for B1 in np.linspace(B1_rango[0], B1_rango[1], n_B):
        D = B1 / Emax1                                   # arranque por el origen
        for s_rel in np.linspace(sigma_rel_rango[0], sigma_rel_rango[1], n_s):
            s1 = s_rel * B1
            s2 = s1 * razon_s if ligar_sigmas else s1
            f1, f2, _ = componentes(x, B1, s1, s2, D, E_gamma1, E_gamma2)
            A1, A2 = _amplitudes_optimas(f1, f2, y)
            res = A1 * f1 + A2 * f2 - y
            sse = float(res @ res)
            if sse < mejor[0]:
                mejor = (sse, [A1, A2, B1, s1, s2, D])
    return mejor[1]


# ----------------------------------------------------------------------
# 3) el ajuste de los 6 parametros
# ----------------------------------------------------------------------
def ajustar_dos_bordes(canales, cuentas, intervalo, p0=None,
                       E_gamma1=E_GAMMA_1, E_gamma2=E_GAMMA_2,
                       ligar_origen=False, ligar_sigmas=False, B1_rango=None):
    """
    Ajusta Theta = {A1, A2, B1, sigma1, sigma2, D} (Ec. 11-13) por minimos
    cuadrados (Ec. 14) sobre el intervalo dado.

    ligar_origen : si True, D = B1/E_max,1 (mapa por el origen, 5 parametros).
    ligar_sigmas : si True, sigma2 = sigma1*sqrt(E_max2/E_max1) (5 parametros).
    B1_rango     : (min,max) para el barrido de semilla; por defecto se deduce
                   del intervalo.

    Devuelve dict con los parametros, los dos puntos de calibracion, y D.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]

    Emax1 = compton_edge_energy(E_gamma1)
    Emax2 = compton_edge_energy(E_gamma2)
    dE = delta_E(E_gamma1, E_gamma2)
    razon_s = np.sqrt(Emax2 / Emax1)

    if p0 is None:
        if B1_rango is None:
            B1_rango = (lo + 0.15 * (hi - lo), hi)
        p0 = semilla_por_barrido(canales, cuentas, intervalo, B1_rango,
                                 E_gamma1=E_gamma1, E_gamma2=E_gamma2,
                                 ligar_sigmas=True)

    A1_0, A2_0, B1_0, s1_0, s2_0, D_0 = p0

    # --- empaquetado: que parametros son libres segun las ligaduras ---
    def desempaquetar(p):
        if ligar_origen and ligar_sigmas:
            A1, A2, B1, s1 = p
            return A1, A2, B1, s1, s1 * razon_s, B1 / Emax1
        if ligar_origen:
            A1, A2, B1, s1, s2 = p
            return A1, A2, B1, s1, s2, B1 / Emax1
        if ligar_sigmas:
            A1, A2, B1, s1, D = p
            return A1, A2, B1, s1, s1 * razon_s, D
        return tuple(p)

    if ligar_origen and ligar_sigmas:
        q0 = [A1_0, A2_0, B1_0, s1_0]
        cotas = ([0, 0, lo, 1e-3], [np.inf, np.inf, hi, hi - lo])
    elif ligar_origen:
        q0 = [A1_0, A2_0, B1_0, s1_0, s2_0]
        cotas = ([0, 0, lo, 1e-3, 1e-3], [np.inf, np.inf, hi, hi - lo, hi - lo])
    elif ligar_sigmas:
        q0 = [A1_0, A2_0, B1_0, s1_0, D_0]
        cotas = ([0, 0, lo, 1e-3, 0.2 * D_0],
                 [np.inf, np.inf, hi, hi - lo, 3.0 * D_0])
    else:
        q0 = [A1_0, A2_0, B1_0, s1_0, s2_0, D_0]
        cotas = ([0, 0, lo, 1e-3, 1e-3, 0.2 * D_0],
                 [np.inf, np.inf, hi, hi - lo, hi - lo, 3.0 * D_0])

    q0 = [min(max(v, a), b) for v, a, b in zip(q0, *cotas)]

    def residuo(p):
        A1, A2, B1, s1, s2, D = desempaquetar(p)
        return modelo_dos_bordes(x, A1, A2, B1, s1, s2, D, E_gamma1, E_gamma2) - y

    sol = least_squares(residuo, q0, bounds=cotas)
    A1, A2, B1, s1, s2, D = desempaquetar(sol.x)
    B2 = B1 + D * dE

    return {
        "A1": A1, "A2": A2, "B1": B1, "B2": B2,
        "sigma1": s1, "sigma2": s2, "D": D, "dE": dE,
        "E_max1": Emax1, "E_max2": Emax2,
        "puntos_calib": [(B1, Emax1), (B2, Emax2)],
        # test de la hipotesis del origen: D vs la pendiente que da B1 solo
        "D_origen": B1 / Emax1,
        "offset_ch": B1 - D * Emax1,       # canal en E=0; deberia ser ~0
        "coste": sol.cost, "exito": sol.success,
        "intervalo": (lo, hi), "p0": p0,
    }


# ======================================================================
#  DEMO: primero verdad conocida, despues se usa en calibrar_ncc.py
# ======================================================================
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from etapas.fit_single import _espectro_sintetico, ajustar_borde_unico

    Emax1 = compton_edge_energy(E_GAMMA_1)
    Emax2 = compton_edge_energy(E_GAMMA_2)

    # --- verdad conocida: dos bordes por el origen, sigma ~ sqrt(E) ---
    B1_true = 600.0
    D_true = B1_true / Emax1
    B2_true = B1_true + D_true * delta_E()
    s1_true = 45.0
    s2_true = s1_true * np.sqrt(Emax2 / Emax1)
    A1_true, A2_true = 5000.0, 4500.0

    canales = np.linspace(0, 1200, 1200)
    y = _espectro_sintetico(canales,
                            [(E_GAMMA_1, A1_true, B1_true, s1_true),
                             (E_GAMMA_2, A2_true, B2_true, s2_true)],
                            semilla=7).astype(float)

    print("=" * 70)
    print("ETAPA 6 - dos bordes solapados (Ec. 11-13), verdad conocida")
    print("=" * 70)
    print(f"  bordes teoricos : {Emax1:.1f} y {Emax2:.1f} keV  (dE = {delta_E():.1f} keV)")
    print(f"  VERDAD  B1={B1_true:.1f}  B2={B2_true:.1f}  s1={s1_true:.1f}  "
          f"s2={s2_true:.1f}  D={D_true:.4f} canal/keV")

    intervalo = (0.45 * B1_true, 1.55 * B1_true)
    r2 = ajustar_dos_bordes(canales, y, intervalo)
    r1 = ajustar_borde_unico(canales, y, E_GAMMA_1, intervalo)

    print(f"\n  --- ajuste de DOS bordes (6 parametros) ---")
    print(f"  B1 = {r2['B1']:7.1f}  ({r2['B1']-B1_true:+.1f})     "
          f"B2 = {r2['B2']:7.1f}  ({r2['B2']-B2_true:+.1f})")
    print(f"  s1 = {r2['sigma1']:7.1f}  ({r2['sigma1']-s1_true:+.1f})     "
          f"s2 = {r2['sigma2']:7.1f}  ({r2['sigma2']-s2_true:+.1f})")
    print(f"  A1 = {r2['A1']:7.0f}  (real {A1_true:.0f})     "
          f"A2 = {r2['A2']:7.0f}  (real {A2_true:.0f})")
    print(f"  D  = {r2['D']:.4f} canal/keV   (real {D_true:.4f}, "
          f"error {100*(r2['D']/D_true-1):+.2f}%)")
    print(f"  test del origen: D_origen = {r2['D_origen']:.4f}, "
          f"offset = {r2['offset_ch']:+.1f} canales (deberia ser ~0)")

    print(f"\n  --- el mismo espectro con UN solo borde (Etapa 3) ---")
    print(f"  B  = {r1['edge_ch']:7.1f}   -> sesgo de {r1['edge_ch']-B1_true:+.1f} canales "
          f"({100*(r1['edge_ch']/B1_true-1):+.1f}%) respecto de B1")
    print(f"  s  = {r1['sigma']:7.1f}   -> s/B = {r1['sigma']/r1['edge_ch']:.3f} "
          f"(inflado; el real es {s1_true/B1_true:.3f})")

    # ------------------------------ figura ------------------------------
    xf = np.linspace(intervalo[0], intervalo[1], 800)
    f1, f2, _ = componentes(xf, r2["B1"], r2["sigma1"], r2["sigma2"], r2["D"])

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))

    ax[0].step(canales, y, where="mid", color="0.75", lw=1, label="datos sinteticos")
    ax[0].plot(xf, r2["A1"] * f1, "g--", lw=1.6, label=f"f1 (963 keV), B1={r2['B1']:.0f}")
    ax[0].plot(xf, r2["A2"] * f2, "m--", lw=1.6, label=f"f2 (1118 keV), B2={r2['B2']:.0f}")
    ax[0].plot(xf, r2["A1"] * f1 + r2["A2"] * f2, "b", lw=2.2, label="S_comp (Ec. 13)")
    ax[0].axvline(B1_true, color="k", ls="--", alpha=0.7, label=f"B1 real = {B1_true:.0f}")
    ax[0].axvline(B2_true, color="k", ls=":", alpha=0.7, label=f"B2 real = {B2_true:.0f}")
    ax[0].set_xlim(150, 1000)
    ax[0].set_title("Dos bordes: el modelo los separa")
    ax[0].set_xlabel("Canal"); ax[0].set_ylabel("Cuentas")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    ax[1].step(canales, y, where="mid", color="0.75", lw=1, label="datos sinteticos")
    ax[1].plot(xf, r2["A1"] * f1 + r2["A2"] * f2, "b", lw=2.2, label="dos bordes")
    ax[1].plot(xf, forward_model(xf, E_GAMMA_1, r1["A"], r1["edge_ch"], r1["sigma"]),
               "r--", lw=2, label="un borde (sesgado)")
    ax[1].axvline(B1_true, color="k", ls="--", alpha=0.7, label=f"B1 real = {B1_true:.0f}")
    ax[1].axvline(r1["edge_ch"], color="r", ls=":", label=f"B un borde = {r1['edge_ch']:.0f}")
    ax[1].set_xlim(150, 1000)
    ax[1].set_title("Por que importa: el ajuste de un borde queda corrido")
    ax[1].set_xlabel("Canal"); ax[1].set_ylabel("Cuentas")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("fit_double_demo.png"), dpi=130)
    print("\nFigura guardada: figuras/fit_double_demo.png")
