"""
model.py — Etapa 2: el forward model (Seccion II-A, mitad "modelo").

Toma la lineshape teorica de physics.py (Etapa 1) y hace dos cosas (Ec. 8-9):
  1. la lleva de espacio de ENERGIA a espacio de CANAL  (el parametro B)
  2. la convoluciona con una gaussiana de ancho sigma   (el broadening)

Resultado: S_model(canal; A, edge_ch, sigma), la curva que en la Etapa 3
ajustaremos contra tu espectro medido. Un solo borde por ahora.
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/model.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
from etapas.physics import dsigma_dE, compton_edge_energy


def _kernel_gaussiano(sigma_muestras, truncate=4.0, max_radio=None):
    """Gaussiana discreta NORMALIZADA (suma 1 -> conserva area al convolucionar).
    max_radio: cota superior al radio, para que el kernel nunca sea mas largo
    que la senal (si no, np.convolve mode='same' devuelve el largo del kernel)."""
    if sigma_muestras <= 0:
        return np.array([1.0])
    radio = int(truncate * sigma_muestras + 0.5)
    if max_radio is not None:
        radio = min(radio, max_radio)
    radio = max(radio, 1)
    x = np.arange(-radio, radio + 1)
    k = np.exp(-0.5 * (x / sigma_muestras) ** 2)
    return k / k.sum()


def forward_model(canales, E_gamma, A, edge_ch, sigma, pad=True, max_pad=20000):
    """
    S_model(canal): lineshape ideal en canal, convolucionada con gaussiana.

    canales : array de canales, UNIFORMEMENTE espaciado (los bins de tu histograma)
    E_gamma : energia del gamma [keV] (fija la FORMA via E_max)
    A       : amplitud
    edge_ch : canal donde cae el Compton edge  (el 'B' del paper)
    sigma   : ancho gaussiano en unidades de CANAL (resolucion del detector)
    pad     : evalua la convolucion sobre una grilla EXTENDIDA (hasta el canal 0
              por la izquierda y 4 sigma por la derecha) y despues recorta.

    Mapa canal->energia: lineal por el origen, con el edge_ch <-> E_max.

    Por que el pad importa: si se convoluciona solo dentro de la ventana de
    ajuste, el modelo pierde los eventos de MAS ABAJO del limite izquierdo que
    el broadening deberia meter hacia adentro (np.convolve rellena con ceros).
    El modelo queda flojo cerca del borde izquierdo y el ajuste lo compensa
    ACHICANDO sigma: en la prueba sintetica de dos bordes eso daba sigma un 18%
    bajo. Con pad=True la ventana de ajuste ya no muerde la convolucion.
    """
    canales = np.asarray(canales, dtype=float)
    E_max = compton_edge_energy(E_gamma)
    n = len(canales)

    # 0) grilla extendida (misma resolucion) para que la convolucion no se trunque
    n_izq = n_der = 0
    grid = canales
    if pad and n > 1:
        paso_g = (canales[-1] - canales[0]) / (n - 1)
        if paso_g > 0:
            n_izq = int(np.ceil(max(canales[0], 0.0) / paso_g))
            n_der = int(np.ceil(4.0 * abs(sigma) / paso_g))
            n_izq = min(n_izq, max_pad)
            n_der = min(n_der, max_pad)
            if n_izq or n_der:
                izq = canales[0] - paso_g * np.arange(n_izq, 0, -1)
                der = canales[-1] + paso_g * np.arange(1, n_der + 1)
                grid = np.concatenate([izq, canales, der])

    # 1) canal -> energia:  el canal 'edge_ch' corresponde a la energia E_max
    E_e = grid * (E_max / edge_ch)
    ideal = dsigma_dE(E_e, E_gamma, A=A)     # lineshape ideal, ya en espacio de canal

    # 2) convolucion con la gaussiana (sigma en canales -> en muestras del grid)
    paso = np.mean(np.diff(grid))
    max_radio = (len(grid) - 1) // 2             # kernel nunca mas largo que la senal
    k = _kernel_gaussiano(sigma / paso, max_radio=max_radio)
    borroso = np.convolve(ideal, k, mode="same")  # ahora siempre len == len(grid)
    return borroso[n_izq:n_izq + n] if grid is not canales else borroso


def punto_mitad_altura(canales, y):
    """Canal del punto al 50% del maximo sobre el flanco DESCENDENTE (metodo clasico)."""
    i_max = int(np.argmax(y))
    mitad = 0.5 * y[i_max]
    post = y[i_max:]
    cruces = np.where(post <= mitad)[0]
    return canales[i_max + cruces[0]] if len(cruces) else np.nan


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Eg = 1173.2                     # una linea del Co-60 (borde teorico 963 keV)
    edge_ch = 600.0                 # <- canal ARBITRARIO del borde (unidades de tu eje)
    canales = np.linspace(0, 1000, 1024)

    plt.figure(figsize=(9, 5.5))
    for sigma in [0, 15, 40, 80]:
        y = forward_model(canales, Eg, A=1.0, edge_ch=edge_ch, sigma=sigma)
        etiqueta = "ideal (sigma=0)" if sigma == 0 else f"sigma = {sigma} canales"
        plt.plot(canales, y, lw=1.7, label=etiqueta)

    # el borde teorico verdadero (donde deberia estar) y el 50% del caso mas borroso
    y80 = forward_model(canales, Eg, A=1.0, edge_ch=edge_ch, sigma=80)
    ch50 = punto_mitad_altura(canales, y80)
    plt.axvline(edge_ch, color="k", ls="--", alpha=0.7, label=f"borde real B = {edge_ch:.0f}")
    plt.axvline(ch50, color="red", ls=":", alpha=0.8, label=f"50% de max (sigma=80) = {ch50:.0f}")

    plt.title("Etapa 2 - Forward model: la gaussiana derrite el escalon")
    plt.xlabel("Canal  (area del pulso, u.a.)")
    plt.ylabel("S_model  (u.a.)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(rutas.figura("forward_model_demo.png"), dpi=130)
    print(f"borde real B      = {edge_ch:.1f} canales")
    print(f"punto 50% (s=80)  = {ch50:.1f} canales  ->  sesgo = {ch50-edge_ch:+.1f} canales")
    print("Figura guardada: figuras/forward_model_demo.png")
