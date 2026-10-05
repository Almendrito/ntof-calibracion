"""
fit_single.py — Etapa 3: primer ajuste sobre datos, UN solo borde.

Ajusta el forward model S_model(A, edge_ch, sigma) de la Etapa 2 a un espectro
(canales, cuentas) sobre un intervalo elegido a mano, con least_squares.
E_gamma es FIJO (fuente conocida). El parametro edge_ch (=B) es la posicion
del borde -> junto con E_max da tu primer par (canal, energia).

Este es el metodo "ingenuo": optimizador local. En la Etapa 5 lo comparamos
con simulated annealing para ver cuando se queda pegado en minimos locales.
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/fit_single.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
from scipy.optimize import least_squares
from etapas.model import forward_model, punto_mitad_altura
from etapas.physics import compton_edge_energy


def ajustar_borde_unico(canales, cuentas, E_gamma, intervalo, p0=None):
    """
    Devuelve dict con A, edge_ch, sigma, la energia del borde (E_max) y el par
    (canal, energia) para la calibracion.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]

    # --- estimaciones iniciales razonables ---
    if p0 is None:
        edge0 = 0.5 * (lo + hi)
        sigma0 = 0.12 * (hi - lo)
        # A escala linealmente: fijamos A para que el pico del modelo ~ pico de datos
        forma = forward_model(x, E_gamma, A=1.0, edge_ch=edge0, sigma=sigma0)
        A0 = (np.max(y) / np.max(forma)) if np.max(forma) > 0 else 1.0
        p0 = [A0, edge0, sigma0]

    # --- residuo a minimizar: modelo - datos (Ec. 14, minimos cuadrados) ---
    def residuo(p):
        A, edge_ch, sigma = p
        return forward_model(x, E_gamma, A, edge_ch, sigma) - y

    # cotas: A>0 ; el borde dentro del intervalo ; sigma acotado y > 0
    # el tope de sigma se abre si la semilla lo pide (una semilla del NCC puede
    # traer sigma > ancho del intervalo; least_squares aborta si p0 viola cotas)
    sig_max = max(hi - lo, 2.0 * p0[2])
    bounds = ([0.0, lo, 1e-3], [np.inf, hi, sig_max])
    p0 = [min(max(v, b0), b1) for v, b0, b1 in zip(p0, *bounds)]
    sol = least_squares(residuo, p0, bounds=bounds)

    A, edge_ch, sigma = sol.x
    E_max = compton_edge_energy(E_gamma)
    return {
        "A": A, "edge_ch": edge_ch, "sigma": sigma,
        "E_max": E_max, "punto_calib": (edge_ch, E_max),
        "coste": sol.cost, "exito": sol.success,
    }


# ---------- utilidades para generar datos SINTETICOS de prueba ----------
def _espectro_sintetico(canales, edges, semilla=0):
    """Suma de bordes (cada uno = forward_model) + ruido de Poisson. 'edges' = lista
    de (E_gamma, A, edge_ch, sigma). Devuelve cuentas ~ enteras con ruido."""
    rng = np.random.default_rng(semilla)
    limpio = np.zeros_like(canales, float)
    for Eg, A, ech, s in edges:
        limpio += forward_model(canales, Eg, A, ech, s)
    limpio = np.clip(limpio, 0, None)
    return rng.poisson(limpio)  # ruido de Poisson -> como cuentas reales


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    canales = np.linspace(0, 1000, 1024)
    Eg1, Eg2 = 1173.2, 1332.5           # lineas del Co-60

    # =========================================================
    # PRUEBA 1: UN borde sintetico con B=600, sigma=40 conocidos
    # =========================================================
    B_true, s_true = 600.0, 40.0
    esc = 5000.0                        # escala de cuentas
    y1 = _espectro_sintetico(canales, [(Eg1, esc, B_true, s_true)], semilla=1)

    r1 = ajustar_borde_unico(canales, y1, Eg1, intervalo=(400, 800))
    ch50 = punto_mitad_altura(canales, y1.astype(float))

    print("=== PRUEBA 1: un borde, verdad conocida ===")
    print(f"  B verdadero        = {B_true:.1f}")
    print(f"  B ajustado (fit)   = {r1['edge_ch']:.1f}   (error {r1['edge_ch']-B_true:+.1f})")
    print(f"  B metodo 50%       = {ch50:.1f}   (error {ch50-B_true:+.1f})")
    print(f"  sigma verdadero    = {s_true:.1f}  ->  ajustado {r1['sigma']:.1f}")

    # =========================================================
    # PRUEBA 2: DOS bordes solapados (Co-60), ajuste de UN borde
    # =========================================================
    B1 = 600.0                          # borde 963 keV
    B2 = B1 * (compton_edge_energy(Eg2) / compton_edge_energy(Eg1))  # 1118 keV, lineal
    y2 = _espectro_sintetico(canales,
                             [(Eg1, esc, B1, 40.0), (Eg2, 0.9*esc, B2, 40.0)],
                             semilla=2)
    r2 = ajustar_borde_unico(canales, y2, Eg1, intervalo=(400, 850))

    print("\n=== PRUEBA 2: dos bordes solapados, ajuste de UNO ===")
    print(f"  bordes verdaderos  = {B1:.1f} (963 keV) y {B2:.1f} (1118 keV)")
    print(f"  medio entre bordes = {0.5*(B1+B2):.1f}")
    print(f"  B ajustado (uno)   = {r2['edge_ch']:.1f}")
    print(f"  sigma ajustado     = {r2['sigma']:.1f}   (vs 40 real de cada borde)")

    # ---- figura ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    ax1.step(canales, y1, where="mid", color="0.6", label="datos sinteticos")
    xf = np.linspace(400, 800, 500)
    ax1.plot(xf, forward_model(xf, Eg1, r1["A"], r1["edge_ch"], r1["sigma"]),
             "b", lw=2, label="ajuste (modelo)")
    ax1.axvline(B_true, color="k", ls="--", label=f"B real = {B_true:.0f}")
    ax1.axvline(r1["edge_ch"], color="b", ls=":", label=f"B fit = {r1['edge_ch']:.0f}")
    ax1.axvline(ch50, color="r", ls=":", label=f"50% = {ch50:.0f}")
    ax1.set_title("Prueba 1: un borde (el fit recupera B sin sesgo)")
    ax1.set_xlabel("Canal"); ax1.set_ylabel("Cuentas"); ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    ax1.set_xlim(300, 900)

    ax2.step(canales, y2, where="mid", color="0.6", label="Co-60 sintetico (2 bordes)")
    xf2 = np.linspace(400, 850, 500)
    ax2.plot(xf2, forward_model(xf2, Eg1, r2["A"], r2["edge_ch"], r2["sigma"]),
             "b", lw=2, label="ajuste de UN borde")
    ax2.axvline(B1, color="k", ls="--", label=f"borde 963 = {B1:.0f}")
    ax2.axvline(B2, color="g", ls="--", label=f"borde 1118 = {B2:.0f}")
    ax2.axvline(r2["edge_ch"], color="b", ls=":", label=f"B fit = {r2['edge_ch']:.0f}")
    ax2.set_title("Prueba 2: dos bordes -> el fit de uno queda en el medio")
    ax2.set_xlabel("Canal"); ax2.set_ylabel("Cuentas"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
    ax2.set_xlim(300, 950)

    plt.tight_layout()
    plt.savefig(rutas.figura("fit_single_demo.png"), dpi=130)
    print("\nFigura guardada: figuras/fit_single_demo.png")
