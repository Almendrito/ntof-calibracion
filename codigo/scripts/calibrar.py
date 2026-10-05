"""
calibrar.py — corre TODO de una sola vez:
  lee tus XML  ->  muestra el espectro  ->  te pide el intervalo del borde
  ->  ajusta  ->  te da el punto de calibracion.

Ejecuta:   python3 calibrar.py

Debe estar en la MISMA carpeta que physics.py, model.py, fit_single.py,
spectrum.py y tus .xml.
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/scripts/calibrar.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
import matplotlib.pyplot as plt
from datos_io import spectrum
from etapas.fit_single import ajustar_borde_unico
from etapas.model import forward_model

# ---------- tus archivos (edita solo si cambian de nombre/carpeta) ----------
ARCHIVO_MUESTRA = rutas.MUESTRA
ARCHIVO_FONDO   = rutas.FONDO
UMBRAL_RUIDO = 0.05      # umbral peak-to-peak [mV]
NUM_CANALES  = 500
E_GAMMA      = 1173.2    # linea del Co-60 a ajustar (borde = 963 keV)
RANGO_MAX    = 2500.0    # [mV*ns]; con el area corregida la fisica vive bajo 1900
# ---------------------------------------------------------------------------

# 1) espectro neto
r = spectrum.espectro_neto(ARCHIVO_MUESTRA, ARCHIVO_FONDO,
                           umbral_ruido=UMBRAL_RUIDO, num_canales=NUM_CANALES,
                           rango_max=RANGO_MAX, recortar_negativos=False)
can, neto = r["canales"], r["neto"]
print(f"rango_max = {r['rango_max']:.3f} mV*ns | t_muestra={r['t_muestra']:.0f}s "
      f"| t_fondo={r['t_fondo']:.0f}s | cuentas netas={neto.sum():.0f}")

# 2) mostrar el espectro para ubicar el borde
plt.figure(figsize=(10, 6))
plt.step(can, np.maximum(neto, 0), where="mid", color="red",  lw=1.6, label="neto")
plt.step(can, r["muestra"],        where="mid", color="blue", alpha=0.35, label="muestra")
plt.step(can, r["fondo"],          where="mid", color="gray", alpha=0.35, label="fondo")
plt.xlabel("Canal (area del pulso) [mV*ns]"); plt.ylabel("Cuentas normalizadas")
plt.title("Ubica el borde (la bajada). CIERRA esta ventana para continuar.")
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig(rutas.figura("espectro_neto.png"), dpi=130)
print("\nSe abrio el espectro (y se guardo en 'espectro_neto.png').")
print("Mira la bajada del espectro NETO (rojo) y CIERRA la ventana para seguir.")
plt.show()

# 3) pedir el intervalo por teclado
print("\nEscribe el intervalo del borde como   lo,hi   (en mV*ns).")
print("Regla: 'lo' un poco antes de la bajada, 'hi' un poco despues, centrado en el borde.")
lo, hi = [float(s) for s in input("Intervalo lo,hi = ").replace(" ", "").split(",")]

# 4) ajustar y reportar
fit = ajustar_borde_unico(can, neto, E_gamma=E_GAMMA, intervalo=(lo, hi))
B, sigma, Emax = fit["edge_ch"], fit["sigma"], fit["E_max"]
print("\n===== RESULTADO =====")
print(f"  B (canal del borde) = {B:.3f} mV*ns")
print(f"  sigma (resolucion)  = {sigma:.3f} mV*ns")
print(f"  E_max (energia)     = {Emax:.1f} keV")
print(f"  PUNTO DE CALIBRACION -> ({B:.3f} mV*ns , {Emax:.1f} keV)")
print("  (Co-60 tiene 2 bordes solapados -> este punto es APROXIMADO;")
print("   el punto sin sesgo llega con el modelo de dos bordes, Etapa 6.)")

# 5) mostrar el ajuste
plt.figure(figsize=(10, 6))
plt.step(can, np.maximum(neto, 0), where="mid", color="0.6", label="neto")
xf = np.linspace(lo, hi, 500)
plt.plot(xf, forward_model(xf, E_GAMMA, fit["A"], B, sigma), "b", lw=2, label="ajuste")
plt.axvline(B, color="b", ls=":", label=f"B = {B:.2f} mV*ns")
plt.xlim(0, hi * 1.4)
plt.xlabel("Canal (area del pulso) [mV*ns]"); plt.ylabel("Cuentas")
plt.title("Ajuste del borde (la curva azul debe pegarse a la bajada)")
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.savefig(rutas.figura("ajuste_borde.png"), dpi=130)
plt.show()
print("\nGuardado: 'ajuste_borde.png'. Si la curva azul NO se pega a la bajada,")
print("vuelve a correr y da un intervalo mas centrado en el borde.")
