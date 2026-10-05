"""Figura de verificacion de la Etapa 4 (NCC) sobre los datos reales de Co-60."""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python verificaciones/_verif_ncc_fig.py` pone en sys.path la carpeta DEL ARCHIVO, no
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

import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datos_io import spectrum
from etapas.ncc import (buscar_intervalo, pearson, _z, plantilla_referencia, fin_rel_kn)
from etapas.physics import compton_edge_energy

M = rutas.MUESTRA
F = rutas.FONDO
EG1 = 1173.2
E1 = compton_edge_energy(EG1)

r = spectrum.espectro_neto(M, F, num_canales=500, rango_max=2500.0,
                           recortar_negativos=False)
can, neto = r["canales"], r["neto"]
suave = spectrum.suavizar(neto, ventana=21, orden=3)
region = (0.08 * can[-1], can[-1])
k = buscar_intervalo(can, suave, EG1, region=region, paso=1,
                     n_puntos="auto", cuentas_crudas=neto)

s_rel, ini = k["sigma_rel"], k["inicio_rel"]
fin, n = fin_rel_kn(s_rel), k["n_puntos"]
S2 = plantilla_referencia(EG1, sigma_rel=s_rel, n_puntos=n,
                          inicio_rel=ini, fin_rel=fin)[1]
z2 = _z(S2)
B = k["edge_ch"]

def ventana(f):
    lo, hi = ini * B * f, fin * B * f
    s1 = np.interp(np.linspace(lo, hi, n), can, suave)
    return lo, hi, s1, pearson(s1, z2)

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))

# (a) rho(B): el criterio de la Etapa 4
Bc, rc = k["B_cands"], k["rho_de_B"]
ok = rc > -1.5
ax[0, 0].plot(Bc[ok], rc[ok], "b", lw=1.6)
ax[0, 0].axvline(B, color="g", ls="--", label=f"B* = {B:.0f} mV·ns  (ρ = {k['rho']:.4f})")
ax[0, 0].set_xlabel("B candidato [mV·ns]"); ax[0, 0].set_ylabel("ρ de Pearson")
ax[0, 0].set_title("(a) El barrido: ρ vs B — máximo ÚNICO e interior")
ax[0, 0].legend(fontsize=9); ax[0, 0].grid(alpha=0.3)

# (b) espectro con la ventana ganadora y dos malas
ax[0, 1].step(can, neto, where="mid", color="0.75", lw=0.9, label="neto (crudo)")
ax[0, 1].plot(can, suave, "k", lw=1.3, label="suavizado (S-G 21)")
for f, col, et in [(1.00, "green", "ganadora"), (0.85, "orange", "B −15%"),
                   (0.70, "red", "B −30%")]:
    lo, hi, _, rho = ventana(f)
    ax[0, 1].axvspan(lo, hi, color=col, alpha=0.12)
    ax[0, 1].plot([], [], color=col, lw=6, alpha=0.35,
                  label=f"{et}: [{lo:.0f},{hi:.0f}]  ρ={rho:.3f}")
ax[0, 1].axvline(B, color="g", ls="--")
ax[0, 1].set_xlim(0, 2500); ax[0, 1].set_xlabel("Área del pulso [mV·ns]")
ax[0, 1].set_ylabel("Cuentas netas")
ax[0, 1].set_title("(b) Ventanas candidatas sobre el espectro real")
ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)

# (c) lo que realmente compara Pearson: S1 y S2 estandarizados
for f, col, et in [(1.00, "green", "ganadora"), (0.85, "orange", "B −15%"),
                   (0.70, "red", "B −30%")]:
    _, _, s1, rho = ventana(f)
    ax[1, 0].plot(np.linspace(0, 1, n), _z(s1), color=col, lw=1.2,
                  label=f"S1 {et}  (ρ={rho:.3f})")
ax[1, 0].plot(np.linspace(0, 1, n), z2, "k", lw=2.4, label="S2 (Klein-Nishina)")
ax[1, 0].set_xlabel(f"punto reinterpolado (0..1, n = {n})")
ax[1, 0].set_ylabel("z-score")
ax[1, 0].set_title("(c) Lo que compara Pearson: forma contra forma")
ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)

# (d) rho vs sigma_rel: la plantilla tambien estima la resolucion
ss = np.array([0.03, 0.06, 0.09, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.21])
rr = [buscar_intervalo(can, suave, EG1, sigma_rel=s, region=region, paso=1,
                       n_puntos="auto", cuentas_crudas=neto) for s in ss]
ax[1, 1].plot(ss, [x["rho"] for x in rr], "o-", color="purple")
ax[1, 1].axvline(s_rel, color="g", ls="--", label=f"σ_rel* = {s_rel:.3f}")
ax[1, 1].set_xlabel("σ_rel de la plantilla"); ax[1, 1].set_ylabel("ρ máximo")
ax[1, 1].set_title("(d) ρ vs σ_rel: el máximo es interior a la grilla")
ax[1, 1].legend(fontsize=9); ax[1, 1].grid(alpha=0.3)
ax2 = ax[1, 1].twinx()
ax2.plot(ss, [x["edge_ch"] for x in rr], "s--", color="0.5", ms=4)
ax2.set_ylabel("B [mV·ns] (gris)", color="0.5")

plt.tight_layout()
plt.savefig(rutas.figura("ncc_verificacion.png"), dpi=130)
print("figuras/ncc_verificacion.png")
print(f"B={B:.1f}  rho={k['rho']:.5f}  sigma_rel={s_rel}  n={n}  "
      f"ventana=[{k['lo']:.0f},{k['hi']:.0f}]  D_seed={k['D_seed']:.4f}")
print("B spread sobre sigma_rel:", f"{min(x['edge_ch'] for x in rr):.0f}-{max(x['edge_ch'] for x in rr):.0f}")
