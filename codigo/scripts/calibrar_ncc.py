"""
calibrar_ncc.py - Etapa 4 sobre DATOS REALES: la logica del paper, sin simulacion.

Cadena completa, igual que el paper (Fig. 2 del paper):
    XML -> espectro neto -> NCC contra la plantilla Klein-Nishina -> ajuste
    -> punto de calibracion (canal del borde, E_max)

La unica "teoria" que entra es la lineshape de Klein-Nishina (Ec. 10): no hay
simulacion Monte Carlo en ningun lado. El espectro medido se compara contra esa
prediccion analitica convolucionada con la gaussiana del detector (Ec. 8-9).

A diferencia de calibrar.py, aqui NADIE elige el intervalo a mano: lo elige el
NCC (Ec. 15).

Uso:   python calibrar_ncc.py            (abre las figuras)
       python calibrar_ncc.py --no-show  (solo las guarda)
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/scripts/calibrar_ncc.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas


import sys
import numpy as np

MOSTRAR = "--no-show" not in sys.argv
import matplotlib
if not MOSTRAR:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

from datos_io import spectrum
from etapas.ncc import buscar_intervalo, plantilla_referencia
from etapas.fit_single import ajustar_borde_unico
from etapas.fit_double import ajustar_dos_bordes, componentes
from etapas.model import forward_model, punto_mitad_altura
from etapas.physics import compton_edge_energy

# ----------------------------- configuracion -----------------------------
ARCHIVO_MUESTRA = rutas.MUESTRA
ARCHIVO_FONDO   = rutas.FONDO
UMBRAL_RUIDO = 0.05        # umbral peak-to-peak [mV]
NUM_CANALES  = 500
RANGO_MAX    = 2500.0      # [mV*ns] FIJO: con None, unos cosmicos lo estiran y el
                           # borde queda sub-muestreado (CLAUDE.md 7.7). El valor
                           # era 8000 con el area vieja; con el area corregida
                           # (CLAUDE.md 7.33) toda la fisica vive bajo 1900.
REGION_NCC   = (250.0, 2500.0)    # zona de busqueda: excluye el pico de bajo canal
                                  # (umbral de trigger), que no es Compton
E_GAMMA      = 1173.2      # linea del Co-60 que se ajusta (borde teorico 963.4 keV)
VENTANA_SUAV = 21          # bins del Savitzky-Golay para el NCC
# -------------------------------------------------------------------------

print("=" * 72)
print("CALIBRACION POR COMPTON EDGE - NCC sobre datos reales de Co-60")
print("=" * 72)
print("  AVISO: esta es la cadena de las Etapas 4-6 (plantilla analitica).")
print("  El resultado VIGENTE sale de `calibrar_v2.py` (Etapa 9), que usa la")
print("  simulacion con la mezcla de scattering ajustable y modela el exceso")
print("  a bajo canal. Aca la plantilla es Klein-Nishina = 1 solo scattering,")
print("  asi que D sale ~2% alto.")

# ---------------------------------------------------------------- 1) espectro
r = spectrum.espectro_neto(ARCHIVO_MUESTRA, ARCHIVO_FONDO,
                           umbral_ruido=UMBRAL_RUIDO, num_canales=NUM_CANALES,
                           rango_max=RANGO_MAX, recortar_negativos=False)
can, neto = r["canales"], r["neto"]
ancho_bin = can[1] - can[0]
print(f"\n  bins        : {NUM_CANALES} de {ancho_bin:.1f} mV*ns  (rango 0-{RANGO_MAX:.0f})")
print(f"  tiempos     : muestra {r['t_muestra']:.0f} s | fondo {r['t_fondo']:.0f} s")
print(f"  tasas       : muestra {50000/r['t_muestra']:.1f} cps | "
      f"fondo {10000/r['t_fondo']:.1f} cps")
print(f"  cuentas netas en el rango: {neto.sum():.0f}")

# suavizado SOLO para el NCC: el Pearson compara forma y el ruido de Poisson
# bin a bin le baja el rho sin aportar informacion (el ajuste usa el neto crudo)
neto_suave = spectrum.suavizar(neto, ventana=VENTANA_SUAV, orden=3)

# --------------------------------------------------- 2) NCC: intervalo automatico
print("\n  --- NCC (Ec. 15): buscando el intervalo del borde ---")
ncc = buscar_intervalo(can, neto_suave, E_GAMMA, region=REGION_NCC, paso=1)
lo, hi = ncc["intervalo"]
print(f"  intervalo elegido : [{lo:.0f}, {hi:.0f}] mV*ns   (rho = {ncc['rho']:.4f})")
print(f"  sigma_rel ganador : {ncc['sigma_rel']:.3f}")
print(f"  B estimado (NCC)  : {ncc['edge_ch']:.0f} mV*ns")
print(f"  sigma estimado    : {ncc['sigma']:.0f} mV*ns")

# el barrido por sigma_rel, para ver cuanto depende el resultado de esa eleccion
print("\n  --- sensibilidad a la plantilla (sigma_rel forzado) ---")
for s_rel in (0.03, 0.06, 0.09, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25):
    rr = buscar_intervalo(can, neto_suave, E_GAMMA, sigma_rel=s_rel,
                          region=REGION_NCC, paso=2)
    print(f"    sigma_rel={s_rel:5.3f} -> rho={rr['rho']:.4f}  "
          f"intervalo=[{rr['lo']:.0f},{rr['hi']:.0f}]  B={rr['edge_ch']:.0f} mV*ns")

# ------------------------------------------- 3) ajuste sembrado por el NCC
fit = ajustar_borde_unico(can, neto, E_GAMMA, (lo, hi), p0=ncc["p0"])
B, sigma = fit["edge_ch"], fit["sigma"]
E_max = compton_edge_energy(E_GAMMA)

# metodo clasico del 50%, para contraste (sesgado por el broadening)
m = (can >= lo) & (can <= hi)
ch50 = punto_mitad_altura(can[m], neto_suave[m])

print("\n" + "=" * 72)
print("  RESULTADO")
print("=" * 72)
print(f"  B (canal del borde)   = {B:.0f} mV*ns      [ajuste sembrado por NCC]")
print(f"  B (solo NCC)          = {ncc['edge_ch']:.0f} mV*ns")
print(f"  B (metodo 50%)        = {ch50:.0f} mV*ns      [clasico, sesgado]")
print(f"  sigma (resolucion)    = {sigma:.0f} mV*ns   ->  sigma/B = {sigma/B:.3f}")
print(f"  E_max (teorica)       = {E_max:.1f} keV   (de E_gamma = {E_GAMMA} keV)")
print(f"\n  PUNTO DE CALIBRACION  ->  ({B:.0f} mV*ns , {E_max:.1f} keV)")
print(f"  ganancia implicada    =  {B/E_max:.3f} mV*ns/keV  "
      f"( = {E_max/B*1000:.2f} keV por 1000 mV*ns )")
print("\n  (el Co-60 tiene DOS bordes solapados: este numero esta sesgado, CLAUDE.md 7.5)")

# --------------------------------------- 3b) Etapa 6: modelo de DOS bordes
print("\n" + "=" * 72)
print("  ETAPA 6 - modelo de dos bordes solapados (Ec. 11-13)")
print("=" * 72)
dob = ajustar_dos_bordes(can, neto, (lo, hi),
                         B1_rango=(0.60 * ncc["edge_ch"], 1.05 * ncc["edge_ch"]))
print(f"  B1 = {dob['B1']:6.0f} mV*ns  <->  {dob['E_max1']:.1f} keV   (borde de 1173 keV)")
print(f"  B2 = {dob['B2']:6.0f} mV*ns  <->  {dob['E_max2']:.1f} keV   (borde de 1333 keV)")
print(f"  s1 = {dob['sigma1']:6.0f} mV*ns   -> s1/B1 = {dob['sigma1']/dob['B1']:.3f}"
      f"  ({2.355*dob['sigma1']/dob['B1']*100:.0f}% FWHM)")
print(f"  s2 = {dob['sigma2']:6.0f} mV*ns   -> s2/B2 = {dob['sigma2']/dob['B2']:.3f}")
print(f"  A1 = {dob['A1']:6.0f}   A2 = {dob['A2']:6.0f}   (A2/A1 = {dob['A2']/dob['A1']:.2f})")
print(f"\n  D (pendiente de calibracion) = {dob['D']:.4f} mV*ns/keV")
print(f"  test del origen: D vs B1/E_max1 = {dob['D']:.4f} vs {dob['D_origen']:.4f}"
      f"  ({100*(dob['D']/dob['D_origen']-1):+.1f}%)")
print(f"  offset en E=0 : {dob['offset_ch']:+.0f} mV*ns  "
      f"({100*dob['offset_ch']/dob['B1']:+.1f}% de B1)")

# --- cuanto se pueden separar realmente los dos bordes ---
sep = dob["B2"] - dob["B1"]
print(f"\n  separacion B2-B1 = {sep:.0f} mV*ns   vs   sigma1 = {dob['sigma1']:.0f} mV*ns")
print(f"  -> los bordes estan a {sep/dob['sigma1']:.2f} sigma: por debajo de 1 sigma la")
print(f"     descomposicion interna (y por lo tanto D) queda mal condicionada, aunque")
print(f"     la POSICION del conjunto (B1) si este bien determinada.")

# ajuste con el mapa amarrado al origen: D deja de ser libre y se estabiliza
dob_o = ajustar_dos_bordes(can, neto, (lo, hi), ligar_origen=True,
                           B1_rango=(0.60 * ncc["edge_ch"], 1.05 * ncc["edge_ch"]))
print(f"\n  --- variante con el mapa amarrado al ORIGEN (5 parametros) ---")
print(f"  B1 = {dob_o['B1']:.0f} mV*ns   B2 = {dob_o['B2']:.0f} mV*ns   "
      f"D = {dob_o['D']:.4f} mV*ns/keV")

print("\n" + "=" * 72)
print("  CALIBRACION RECOMENDADA")
print("=" * 72)
print(f"  D = {dob_o['D']:.3f} mV*ns/keV     ( = {1000/dob_o['D']:.1f} keV por 1000 mV*ns )")
print(f"  puntos: ({dob_o['B1']:.0f} mV*ns, {dob_o['E_max1']:.1f} keV) y "
      f"({dob_o['B2']:.0f} mV*ns, {dob_o['E_max2']:.1f} keV)")
print(f"\n  Corrimiento respecto del ajuste de UN borde: "
      f"{100*(dob['B1']/B - 1):+.1f}% en el punto de 963 keV")
print(f"  sigma/B: {sigma/B:.3f} (un borde) -> {dob['sigma1']/dob['B1']:.3f} (dos bordes)")
print(f"  Ese {2.355*dob['sigma1']/dob['B1']*100:.0f}% de FWHM NO es solo resolucion: el modelo")
print(f"  gaussiano del paper absorbe ahi tambien el scattering Compton multiple, que en")
print(f"  un centellador organico rellena la zona sobre el borde. [fuera del paper]")

# --------------------------------------------------------------- 4) figuras
fig, ax = plt.subplots(2, 2, figsize=(14, 9))

# (a) espectro completo con muestra/fondo/neto
ax[0, 0].step(can, r["muestra"], where="mid", color="blue", alpha=0.35, label="muestra")
ax[0, 0].step(can, r["fondo"], where="mid", color="gray", alpha=0.5, label="fondo")
ax[0, 0].step(can, neto, where="mid", color="red", lw=1.2, label="neto")
ax[0, 0].set_yscale("log")
ax[0, 0].set_title("Espectro medido (log)")
ax[0, 0].set_xlabel("Area del pulso [mV*ns]"); ax[0, 0].set_ylabel("Cuentas norm.")
ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)

# (b) zoom a la zona Compton con el intervalo elegido por el NCC
z = (can >= 1000) & (can <= 7000)
ax[0, 1].step(can[z], neto[z], where="mid", color="0.7", lw=1, label="neto")
ax[0, 1].plot(can[z], neto_suave[z], color="crimson", lw=1.8, label="neto suavizado (NCC)")
ax[0, 1].axvspan(lo, hi, color="green", alpha=0.13,
                 label=f"intervalo NCC (rho={ncc['rho']:.3f})")
ax[0, 1].axvline(ncc["edge_ch"], color="green", ls=":", lw=2,
                 label=f"B por NCC = {ncc['edge_ch']:.0f}")
ax[0, 1].set_title("Intervalo elegido AUTOMATICAMENTE (Ec. 15)")
ax[0, 1].set_xlabel("Area del pulso [mV*ns]"); ax[0, 1].set_ylabel("Cuentas norm.")
ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)

# (c) S1 vs S2 estandarizados: lo que realmente compara Pearson
i, w = ncc["i_ini"], ncc["ancho_bins"]
xq = np.linspace(can[i], can[i + w], 100)
S1 = np.interp(xq, can, neto_suave)
_, S2 = plantilla_referencia(E_GAMMA, sigma_rel=ncc["sigma_rel"],
                             inicio_rel=ncc["inicio_rel"], fin_rel=ncc["fin_rel"])
zz = lambda v: (v - np.mean(v)) / np.std(v)
ax[1, 0].plot(np.arange(100), zz(S1), color="0.35", lw=1.4, label="S1 (datos reales)")
ax[1, 0].plot(np.arange(100), zz(S2), "g", lw=2, label="S2 (Klein-Nishina)")
ax[1, 0].set_title(f"Datos vs teoria, estandarizados: rho = {ncc['rho']:.4f}")
ax[1, 0].set_xlabel("punto reinterpolado (1..100)"); ax[1, 0].set_ylabel("z-score")
ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)

# (d) el ajuste de DOS bordes sobre los datos, con sus componentes
xf = np.linspace(lo, hi, 600)
f1, f2, _ = componentes(xf, dob["B1"], dob["sigma1"], dob["sigma2"], dob["D"])
ax[1, 1].step(can[z], neto[z], where="mid", color="0.7", lw=1, label="neto")
ax[1, 1].plot(xf, dob["A1"] * f1, "g--", lw=1.5,
              label=f"f1 (963 keV), B1={dob['B1']:.0f}")
ax[1, 1].plot(xf, dob["A2"] * f2, "m--", lw=1.5,
              label=f"f2 (1118 keV), B2={dob['B2']:.0f}")
ax[1, 1].plot(xf, dob["A1"] * f1 + dob["A2"] * f2, "b", lw=2.2, label="S_comp (Ec. 13)")
ax[1, 1].plot(xf, forward_model(xf, E_GAMMA, fit["A"], B, sigma), "r:", lw=1.8,
              label=f"un borde (B={B:.0f}, sesgado)")
ax[1, 1].axvline(dob["B1"], color="g", ls=":", lw=1.5)
ax[1, 1].axvline(dob["B2"], color="m", ls=":", lw=1.5)
ax[1, 1].set_title(f"Etapa 6: D = {dob['D']:.3f} mV*ns/keV")
ax[1, 1].set_xlabel("Area del pulso [mV*ns]"); ax[1, 1].set_ylabel("Cuentas norm.")
ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(rutas.figura("calibracion_ncc.png"), dpi=130)
print("\n  Figura guardada: figuras/calibracion_ncc.png")
if MOSTRAR:
    plt.show()
