"""Del intervalo del NCC al ajuste: que dan B y sigma segun la plantilla."""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python verificaciones/_verif_convolucion.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
from ajuste import modelo_ext as me
from etapas.ncc import buscar_intervalo
from etapas.physics import compton_edge_energy
from etapas.model import forward_model
from etapas.fit_single import ajustar_borde_unico
from etapas.fit_double import ajustar_dos_bordes, modelo_dos_bordes, componentes
from ajuste.sim_template import cargar_simulacion, ajustar_pearson, modelo_simulado

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
B0 = k["edge_ch"]
V = (0.25*B0, min(1.55*B0, can[-1]))
m = (can >= V[0]) & (can <= V[1])
x = can[m]
def c2(mod, npar):  return float(np.sum(((neto[m]-mod)/err[m])**2))/(m.sum()-npar)

filas = []

# M1: un borde analitico
f1 = ajustar_borde_unico(can, neto, EG1, V, p0=k["p0"])
y1 = forward_model(x, EG1, f1["A"], f1["edge_ch"], f1["sigma"])
sE = f1["sigma"]/(f1["edge_ch"]/E1)
filas.append(("M1  Klein-Nishina, 1 borde", f1["edge_ch"]/E1, 2.355*sE/E1*100,
              c2(y1, 3), y1, "tab:red"))

# M2: dos bordes analiticos, origen ligado
d = ajustar_dos_bordes(can, neto, V, ligar_origen=True,
                       B1_rango=(0.75*B0, 1.15*B0))
y2 = modelo_dos_bordes(x, d["A1"], d["A2"], d["B1"], d["sigma1"], d["sigma2"],
                       d["D"], EG1, EG2)
sE2 = d["sigma1"]/d["D"]
filas.append(("M2  Klein-Nishina, 2 bordes", d["D"], 2.355*sE2/E1*100,
              c2(y2, 5), y2, "tab:orange"))

# M3: plantilla Geant4, {A, D, alpha}
cen, pes = cargar_simulacion(rama="Lvis_keVee", bins=175)
var_fn = spectrum.varianza_predicha(r["fondo"][m], r["k_m"], r["k_f"])
g = ajustar_pearson(can, neto, cen, pes, V, var_fn, errores=err, E_ref=E1)
y3 = modelo_simulado(x, cen, pes, g["A"], g["D"], g["alpha"])
filas.append(("M3  Geant4, sin exceso", g["D"], g["fwhm_rel_ref"]*100,
              c2(y3, 3), y3, "tab:blue"))

# M4: modelo extendido (Etapa 9)
cen4, h1, hm, w_nat = me.cargar_componentes(rama="Lvis_keVee", bins=175)
p0 = np.array([1.0, 1.25, 3.5, w_nat, 0.0, 250.0, 0.0, 0.0])
z = (can > 1000) & (can < 1100)
p0[me.I_A] = neto[z].mean()/me.modelo_ext(can, cen4, h1, hm, p0)[z].mean()
f0 = me.ajustar(can, neto, cen4, h1, hm, (900., 1900.), p0, ["A", "D", "alpha"],
                var_fn=spectrum.varianza_predicha(
                    r["fondo"][(can >= 900) & (can <= 1900)], r["k_m"], r["k_f"]),
                errores=err)
p = f0["p"].copy(); p[me.I_C], p[me.I_LAM] = 80.0, 350.0
LIB = ["A", "D", "alpha", "w", "C", "lambda"]
f4 = me.ajustar(can, neto, cen4, h1, hm, (300., 1900.), p, LIB,
                var_fn=spectrum.varianza_predicha(
                    r["fondo"][(can >= 300) & (can <= 1900)], r["k_m"], r["k_f"]),
                errores=err)
y4 = me.modelo_ext(x, cen4, h1, hm, f4["p"])
filas.append(("M4  Geant4 + exceso (Etapa 9)", f4["D"],
              2.355*f4["p"][me.I_ALPHA]*np.sqrt(E1)/E1*100, f4["chi2_red"],
              y4, "tab:green"))

print(f"ventana de ajuste = [{V[0]:.0f}, {V[1]:.0f}] mV·ns   (NCC: B = {B0:.0f})\n")
print(f"{'modelo':32s} {'D [mV·ns/keV]':>14s} {'FWHM(963)':>10s} {'χ²/ndf':>8s}")
for et, D, fw, c, *_ in filas:
    print(f"{et:32s} {D:14.4f} {fw:9.1f}% {c:8.3f}")
print(f"\nfit de 2 bordes: A2/A1 = {d['A2']/d['A1']:.3f}  "
      f"(la simulación dice 0.95)  offset = {d['offset_ch']:+.1f} mV·ns")

# ------------------------------- figura -------------------------------
fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))

a = ax[0, 0]
a.errorbar(can[m], neto[m], yerr=err[m], fmt="o", ms=2, color="0.6",
           lw=0.6, label="datos netos", zorder=1)
for et, D, fw, c, y, col in filas:
    a.plot(x, y, color=col, lw=1.8, label=f"{et.split('  ')[1]}  (χ²/ndf={c:.2f})")
a.axvspan(V[0], V[1], color="k", alpha=0.04)
a.set_xlim(0, 2100); a.set_xlabel("Área [mV·ns]"); a.set_ylabel("Cuentas netas")
a.set_title("(a) La misma ventana, cuatro plantillas convolucionadas")
a.legend(fontsize=8); a.grid(alpha=0.3)

a = ax[0, 1]
for et, D, fw, c, y, col in filas:
    a.plot(x, (neto[m]-y)/err[m], color=col, lw=1.0)
a.axhline(0, color="k", lw=0.8)
for v in (-3, 3): a.axhline(v, color="k", ls=":", lw=0.7)
a.set_xlabel("Área [mV·ns]"); a.set_ylabel("(dato − modelo)/σ")
a.set_title("(b) Residuos: dónde falla cada plantilla")
a.grid(alpha=0.3)

a = ax[1, 0]
f1c, f2c, _ = componentes(x, d["B1"], d["sigma1"], d["sigma2"], d["D"], EG1, EG2)
a.errorbar(can[m], neto[m], yerr=err[m], fmt="o", ms=2, color="0.6", lw=0.6)
a.plot(x, d["A1"]*f1c, "b--", lw=1.5, label=f"borde 963.4 keV (A1={d['A1']:.0f})")
a.plot(x, d["A2"]*f2c, "r--", lw=1.5, label=f"borde 1118.1 keV (A2={d['A2']:.0f})")
a.plot(x, y2, color="tab:orange", lw=2, label="suma")
a.axvline(d["B1"], color="b", ls=":"); a.axvline(d["B2"], color="r", ls=":")
a.set_xlim(0, 2100); a.set_xlabel("Área [mV·ns]"); a.set_ylabel("Cuentas netas")
a.set_title(f"(c) Los dos bordes NO se separan: A2/A1 = {d['A2']/d['A1']:.2f}"
            f" (debería ser ≈0.95)")
a.legend(fontsize=8); a.grid(alpha=0.3)

a = ax[1, 1]
et = [f.split("  ")[0] for f, *_ in filas]
Ds = [f[1] for f in filas]; fws = [f[2] for f in filas]
cols = [f[5] for f in filas]
a.bar(np.arange(4)-0.18, Ds, 0.36, color=cols, label="D [mV·ns/keVee]")
a2 = a.twinx()
a2.bar(np.arange(4)+0.18, fws, 0.36, color=cols, alpha=0.45, hatch="//")
a2.set_ylabel("FWHM(963 keV) [%]  (rayado)")
a.axhline(f4["D"], color="tab:green", ls="--", lw=1)
a.set_xticks(range(4)); a.set_xticklabels(et)
a.set_ylabel("D [mV·ns/keVee]  (liso)"); a.set_ylim(0, 1.5)
a2.set_ylim(0, 45)
a.set_title("(d) σ es un basurero: absorbe la física que falta")
a.grid(alpha=0.3, axis="y")

plt.tight_layout(); plt.savefig(rutas.figura("ajuste_verificacion.png"), dpi=130)
print("\najuste_verificacion.png")
