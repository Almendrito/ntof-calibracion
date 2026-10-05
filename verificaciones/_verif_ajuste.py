"""Etapa 3/6: del intervalo del NCC al ajuste de la convolucion (B y sigma)."""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python verificaciones/_verif_ajuste.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.join(_os.path.dirname(_os.path.dirname(
    _os.path.abspath(__file__))), "codigo")
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas

import numpy as np
from datos_io import spectrum
from etapas.ncc import buscar_intervalo
from etapas.physics import compton_edge_energy
from etapas.model import forward_model, punto_mitad_altura
from etapas.fit_single import ajustar_borde_unico
from etapas.fit_double import ajustar_dos_bordes, modelo_dos_bordes, delta_E

M = rutas.MUESTRA
F = rutas.FONDO
EG1, EG2 = 1173.2, 1332.5
E1, E2 = compton_edge_energy(EG1), compton_edge_energy(EG2)

r = spectrum.espectro_neto(M, F, num_canales=500, rango_max=2500.0,
                           recortar_negativos=False)
can, neto, err = r["canales"], r["neto"], r["error"]
suave = spectrum.suavizar(neto, ventana=21, orden=3)
k = buscar_intervalo(can, suave, EG1, region=(0.08*can[-1], can[-1]), paso=1,
                     n_puntos="auto", cuentas_crudas=neto)
B_ncc = k["edge_ch"]
V_NCC = (k["lo"], k["hi"])
V_AJ = (0.25*B_ncc, min(1.55*B_ncc, can[-1]))
print(f"NCC: B = {B_ncc:.1f}   ventana NCC = [{V_NCC[0]:.0f},{V_NCC[1]:.0f}]"
      f"   ventana de ajuste = [{V_AJ[0]:.0f},{V_AJ[1]:.0f}]")
print(f"semilla p0 = A={k['p0'][0]:.1f}  B={k['p0'][1]:.1f}  sigma={k['p0'][2]:.1f}")

def chi2(mod, m):
    return float(np.sum(((neto[m]-mod)/err[m])**2))

# ---------------- 1) UN borde: Theta = {A, B, sigma} ----------------
print("\n[1] UN borde (Ec. 8-10), sembrado por NCC")
for et, V in [("ventana NCC   ", V_NCC), ("ventana ajuste", V_AJ)]:
    f = ajustar_borde_unico(can, neto, EG1, V, p0=k["p0"])
    m = (can >= V[0]) & (can <= V[1])
    mod = forward_model(can[m], EG1, f["A"], f["edge_ch"], f["sigma"])
    nd = m.sum() - 3
    print(f"  {et}: B = {f['edge_ch']:7.1f}  sigma = {f['sigma']:6.1f}"
          f"  sigma/B = {f['sigma']/f['edge_ch']:.3f}"
          f"  D = B/E_max = {f['edge_ch']/E1:.4f}"
          f"  chi2/ndf = {chi2(mod, m)/nd:.3f}")

# el metodo clasico del 50%, para tener la referencia del sesgo
xs = np.linspace(V_AJ[0], V_AJ[1], 2000)
f = ajustar_borde_unico(can, neto, EG1, V_AJ, p0=k["p0"])
ys = forward_model(xs, EG1, f["A"], f["edge_ch"], f["sigma"])
c50 = punto_mitad_altura(xs, ys)
print(f"  metodo clasico del 50%: canal = {c50:.1f} -> D = {c50/E1:.4f}"
      f"  ({100*(c50/f['edge_ch']-1):+.1f}% respecto de B del ajuste)")

# ---------------- 2) DOS bordes: Theta = {A1,A2,B1,s1,s2,D} ----------------
print("\n[2] DOS bordes solapados (Ec. 11-13), ventana de ajuste")
dEb = delta_E(EG1, EG2)
print(f"  Delta E de los BORDES = {dEb:.1f} keV (no 159.3 de las lineas gamma)")
B1r = (0.75*B_ncc, 1.15*B_ncc)
for et, lo_, ls_ in [("6p, D libre        ", False, False),
                     ("5p, origen ligado  ", True,  False),
                     ("5p, sigmas ligadas ", False, True),
                     ("4p, ambas ligadas  ", True,  True)]:
    d = ajustar_dos_bordes(can, neto, V_AJ, ligar_origen=lo_, ligar_sigmas=ls_,
                           B1_rango=B1r)
    m = (can >= V_AJ[0]) & (can <= V_AJ[1])
    mod = modelo_dos_bordes(can[m], d["A1"], d["A2"], d["B1"], d["sigma1"],
                            d["sigma2"], d["D"], EG1, EG2)
    npar = 6 - int(lo_) - int(ls_)
    print(f"  {et}: B1={d['B1']:7.1f}  B2={d['B2']:7.1f}  s1/B1={d['sigma1']/d['B1']:.3f}"
          f"  s2/B2={d['sigma2']/d['B2']:.3f}  A2/A1={d['A2']/max(d['A1'],1e-9):5.2f}"
          f"  D={d['D']:.4f}  offset={d['offset_ch']:+7.1f}"
          f"  chi2/ndf={chi2(mod,m)/(m.sum()-npar):.3f}")

# ---------------- 3) estabilidad de D con la ventana ----------------
print("\n[3] estabilidad con la ventana (origen ligado, 5p)")
for lo_ch in (250, 300, 400, 500, 600, 800):
    d = ajustar_dos_bordes(can, neto, (lo_ch, V_AJ[1]), ligar_origen=True,
                           B1_rango=B1r)
    print(f"  lo={lo_ch:4d} -> D={d['D']:.4f}  B1={d['B1']:7.1f}"
          f"  s1/B1={d['sigma1']/d['B1']:.3f}")

# ---------------- 4) que significa sigma como RESOLUCION ----------------
print("\n[4] sigma -> resolucion")
d = ajustar_dos_bordes(can, neto, V_AJ, ligar_origen=True, B1_rango=B1r)
for et, s, Bx, Ex in [("borde 963.4", d["sigma1"], d["B1"], E1),
                      ("borde 1118.1", d["sigma2"], d["B2"], E2)]:
    sE = s / d["D"]
    print(f"  {et} keV: sigma = {s:.1f} mV*ns = {sE:.1f} keV"
          f"   FWHM/E = {2.355*sE/Ex*100:.1f}%   alpha = sigma/sqrt(E) = {sE/np.sqrt(Ex):.2f}")
print(f"\n  referencia Etapa 9 (modelo_ext + MCMC): D_luz = 1.261 +- 0.014,"
      f" 2.79*sqrt(E) keV -> 21.1% FWHM a 963 keV")
