"""
sa.py - Etapa 5: Simulated Annealing (Seccion II-A del paper).

Que es
------
Un optimizador GLOBAL. `least_squares` (Etapas 3 y 6) es LOCAL: baja siempre y
se queda en el primer minimo que encuentra, asi que su resultado depende del
punto de partida (CLAUDE.md 7.4). El SA acepta a veces movimientos que EMPEORAN
el costo, con probabilidad P = exp(-dE/T), y esa es justamente la puerta para
salir de un minimo local. Al enfriar (T -> 0) deja de aceptar subidas y termina
comportandose como un optimizador local.

Receta del paper (Seccion II-A, parrafo de SA)
----------------------------------------------
  T_init = 2000, coeficiente de enfriamiento alpha = 0.95 (rango [0.85, 0.99])
  100 iteraciones de cadena de Markov por temperatura
  x_new = x_current + dx,  dx_i ~ U(-0.1*delta_i, 0.1*delta_i)
  Metropolis: aceptar si dE < 0, si no con probabilidad exp(-dE/T)
  enfriamiento geometrico T_{k+1} = alpha * T_k
  corte: T < 1e-6  o  |dE_mean| < 0.001 en 10 ciclos consecutivos

Dos problemas del paper al implementarlo
----------------------------------------
1. **T_init = 2000 no significa nada por si solo.** La probabilidad de Metropolis
   compara T contra dE, asi que T solo tiene sentido EN LAS UNIDADES DEL COSTO.
   Con la Ec. 14 (suma de cuadrados) sobre estos datos E ~ 1e4 y T=2000 acepta
   casi cualquier cosa; con el error NCC (1-rho ~ 0.01), T=2000 acepta TODO y el
   SA degenera en una caminata aleatoria. Por eso aca T_init="auto" (default)
   calibra T para que la tasa de aceptacion inicial sea ~0.8, midiendola sobre
   movimientos de prueba. Se puede forzar T_init=2000 para reproducir el paper
   al pie de la letra. [fuera del paper]
2. **La ambiguedad Ec.14 vs Metropolis tiene consecuencias.** El paper dice en la
   Ec. 14 que minimiza minimos cuadrados, pero el texto de Metropolis dice que
   E(x) es el "NCC-based matching error". No son intercambiables: Pearson es
   invariante a la amplitud, asi que con el objetivo NCC el parametro A queda
   INDETERMINADO (y en el modelo de dos bordes solo se determina la razon
   A2/A1). Aca estan los dos, y `costo_ncc` proyecta la amplitud por minimos
   cuadrados para que el resultado sea usable.

Corre solo:  python sa.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/sa.py` pone en sys.path la carpeta DEL ARCHIVO, no
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

from etapas.model import forward_model
from etapas.fit_double import modelo_dos_bordes


# ----------------------------------------------------------------------
# 1) el algoritmo
# ----------------------------------------------------------------------
def recocido_simulado(costo, x0, cotas, deltas=None, T_init="auto", alpha=0.95,
                      n_markov=100, T_min=1e-6, tol_energia=1e-3, n_estable=10,
                      frac_paso=0.1, semilla=0, max_evals=200000, verbose=False,
                      paso_adaptativo=True, acept_objetivo=0.35):
    """
    Simulated annealing con la receta del paper.

    costo    : callable(x) -> float, la E(x) a minimizar.
    x0       : punto inicial.
    cotas    : (min, max) por parametro. Define delta_i si no se pasa `deltas`.
    deltas   : rango de busqueda por parametro (el delta_i del paper).
    T_init   : "auto" (calibra para ~80% de aceptacion inicial) o un numero.
    frac_paso: el 0.1 de dx_i ~ U(-0.1*delta_i, 0.1*delta_i).

    paso_adaptativo : [fuera del paper] escala el tamano del paso para mantener
        la tasa de aceptacion cerca de `acept_objetivo`. El paper deja el paso
        FIJO en 0.1*delta_i para siempre, y eso no converge cuando hay muchos
        parametros: a baja temperatura el SA solo sabe proponer saltos grandes,
        que ya no acepta, y no puede refinar. Con 3 parametros da igual; con los
        6 del modelo de dos bordes, el paso fijo termina un factor ~4.5 peor en
        costo que un least_squares bien sembrado (ver CLAUDE.md 7.21).
        Con paso_adaptativo=False se reproduce el paper al pie de la letra.

    Devuelve dict con el mejor x, su costo, y el historial para diagnosticar.
    """
    rng = np.random.default_rng(semilla)
    lo = np.asarray(cotas[0], float)
    hi = np.asarray(cotas[1], float)
    x = np.clip(np.asarray(x0, float), lo, hi)
    d = (hi - lo) if deltas is None else np.asarray(deltas, float)
    d = np.where(np.isfinite(d), d, np.abs(x) + 1.0)      # cotas infinitas -> escala local

    e = float(costo(x))
    n_ev = 1
    mejor_x, mejor_e = x.copy(), e

    # --- calibracion de T_init: T en unidades del costo ---
    if T_init == "auto":
        subidas = []
        for _ in range(200):
            xp = np.clip(x + rng.uniform(-frac_paso * d, frac_paso * d), lo, hi)
            de = float(costo(xp)) - e
            n_ev += 1
            if de > 0:
                subidas.append(de)
        # P = exp(-dE/T) = 0.8  ->  T = dE_medio / ln(1/0.8)
        T = (np.mean(subidas) / np.log(1 / 0.8)) if subidas else 1.0
    else:
        T = float(T_init)
    T_inicial = T

    hist = {"T": [], "mejor": [], "actual": [], "aceptacion": [], "escala": []}
    e_media_prev = None
    estables = 0
    n_ciclos = 0
    escala = 1.0                     # multiplicador del paso (1.0 = el del paper)

    while T > T_min and n_ev < max_evals:
        aceptados = 0
        suma_e = 0.0
        n_it = 0
        for _ in range(n_markov):
            # dx_i ~ U(-0.1 delta_i, 0.1 delta_i)
            xp = np.clip(x + rng.uniform(-frac_paso * escala * d,
                                         frac_paso * escala * d), lo, hi)
            ep = float(costo(xp))
            n_ev += 1
            de = ep - e
            # criterio de Metropolis
            if de < 0 or rng.random() < np.exp(-min(de / T, 700.0)):
                x, e = xp, ep
                aceptados += 1
                if e < mejor_e:
                    mejor_x, mejor_e = x.copy(), e
            suma_e += e
            n_it += 1
            if n_ev >= max_evals:
                break

        tasa = aceptados / max(n_it, 1)
        e_media = suma_e / max(n_it, 1)
        n_ciclos += 1
        hist["T"].append(T)
        hist["mejor"].append(mejor_e)
        hist["actual"].append(e)
        hist["aceptacion"].append(tasa)
        hist["escala"].append(escala)

        # corte por estabilizacion de la energia MEDIA del ciclo.
        # tol relativa: un umbral absoluto de 0.001 (el del paper) depende de las
        # unidades del costo, igual que T_init=2000. Ver el docstring del modulo.
        if e_media_prev is not None and \
                abs(e_media - e_media_prev) < tol_energia * max(abs(e_media_prev), 1e-12):
            estables += 1
            if estables >= n_estable:
                break
        else:
            estables = 0
        e_media_prev = e_media

        if paso_adaptativo:          # [fuera del paper] mantener la aceptacion util
            if tasa > acept_objetivo * 1.5:
                escala *= 1.15
            elif tasa < acept_objetivo * 0.5:
                escala *= 0.85
            escala = min(max(escala, 1e-4), 10.0)

        T *= alpha                                       # enfriamiento geometrico
        if verbose and n_ciclos % 25 == 0:
            print(f"    ciclo {n_ciclos:4d}  T={T:10.4g}  mejor={mejor_e:.6g}  "
                  f"acept={tasa:.2f}  escala={escala:.3g}")

    return {
        "x": mejor_x, "costo": mejor_e, "n_evals": n_ev, "n_ciclos": n_ciclos,
        "T_inicial": T_inicial, "T_final": T, "escala_final": escala, "hist": hist,
    }


def sa_multiarranque(costo, x0, cotas, semillas=(3, 11, 42, 7, 99), **kw):
    """
    Corre el SA con varias semillas y devuelve el de MENOR costo.

    El SA es estocastico: una sola corrida no es un resultado. En la prueba de
    dos bordes, 2 de cada 3 semillas encuentran el optimo global y la tercera se
    queda en un minimo local un factor ~2.7 peor en costo (CLAUDE.md 7.22).
    El costo final es el criterio para elegir: es la misma E(x) para todas.
    """
    mejor = None
    for s in semillas:
        r = recocido_simulado(costo, x0, cotas, semilla=s, **kw)
        if mejor is None or r["costo"] < mejor["costo"]:
            mejor = r
            mejor["semilla"] = s
    return mejor


# ----------------------------------------------------------------------
# 2) las dos funciones de costo que el paper deja ambiguas
# ----------------------------------------------------------------------
def costo_minimos_cuadrados(x_datos, y_datos, modelo, errores=None):
    """E(x) = suma (S1 - S_model)^2  -> la Ec. 14, opcionalmente ponderada."""
    y = np.asarray(y_datos, float)
    w = 1.0 / np.asarray(errores, float) if errores is not None else 1.0

    def E(p):
        r = (modelo(x_datos, p) - y) * w
        return float(r @ r)
    return E


def costo_ncc(x_datos, y_datos, modelo):
    """
    E(x) = 1 - rho(S1, S_model)  -> el "NCC-based matching error" del Metropolis.

    OJO: Pearson es invariante a la amplitud, asi que este costo NO puede
    determinar A (ni A1, A2 por separado: solo su razon). Se usa igual para
    comparar, y despues la amplitud se recupera por proyeccion lineal.
    """
    y = np.asarray(y_datos, float)
    yc = y - y.mean()
    ny = np.sqrt(yc @ yc)

    def E(p):
        m = modelo(x_datos, p)
        mc = m - m.mean()
        nm = np.sqrt(mc @ mc)
        if nm <= 0 or ny <= 0:
            return 2.0
        return 1.0 - float(mc @ yc) / (nm * ny)
    return E


def proyectar_amplitud(x_datos, y_datos, modelo, p, idx_amplitudes):
    """Recupera las amplitudes por minimos cuadrados lineales tras un ajuste NCC."""
    p = np.array(p, float)
    p[idx_amplitudes] = 1.0
    base = modelo(x_datos, p)
    den = float(base @ base)
    if den > 0:
        p[idx_amplitudes] = max(float(base @ np.asarray(y_datos, float)) / den, 0.0)
    return p


# ----------------------------------------------------------------------
# 3) envoltorios para los modelos de este proyecto
# ----------------------------------------------------------------------
def modelo_un_borde(E_gamma):
    """p = [A, B, sigma]"""
    def m(x, p):
        return forward_model(x, E_gamma, p[0], max(p[1], 1e-9), max(p[2], 1e-9))
    return m


def modelo_dos(E_gamma1=1173.2, E_gamma2=1332.5):
    """p = [A1, A2, B1, sigma1, sigma2, D]  (el Theta de la Ec. 11)"""
    def m(x, p):
        return modelo_dos_bordes(x, p[0], p[1], max(p[2], 1e-9),
                                 max(p[3], 1e-9), max(p[4], 1e-9), p[5],
                                 E_gamma1, E_gamma2)
    return m


def costo_dos_bordes_varpro(x_datos, y_datos, E_gamma1=1173.2, E_gamma2=1332.5,
                            errores=None):
    """
    Costo de dos bordes sobre SOLO 4 parametros: p = [B1, sigma1, sigma2, D].

    A1 y A2 son lineales, asi que se resuelven exacto (nnls) en cada evaluacion
    en vez de buscarlas. Le saca dos dimensiones al espacio de busqueda del SA
    y elimina la escala fea de mezclar amplitudes (~5000) con D (~0.6) en el
    mismo vector de pasos. [fuera del paper: el paper busca los 6 a ciegas]
    """
    from etapas.fit_double import componentes, _amplitudes_optimas
    y = np.asarray(y_datos, float)
    w = 1.0 / np.asarray(errores, float) if errores is not None else 1.0

    def E(p):
        f1, f2, _ = componentes(x_datos, max(p[0], 1e-9), max(p[1], 1e-9),
                                max(p[2], 1e-9), p[3], E_gamma1, E_gamma2)
        A1, A2 = _amplitudes_optimas(f1, f2, y)
        r = (A1 * f1 + A2 * f2 - y) * w
        return float(r @ r)
    return E


# ======================================================================
#  DEMO
# ======================================================================
if __name__ == "__main__":
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from etapas.fit_single import _espectro_sintetico, ajustar_borde_unico
    from etapas.fit_double import ajustar_dos_bordes, delta_E
    from etapas.physics import compton_edge_energy

    Eg1, Eg2 = 1173.2, 1332.5
    print("=" * 72)
    print("ETAPA 5 - Simulated Annealing (Seccion II-A)")
    print("=" * 72)

    # =========================================================
    # PRUEBA 1: un borde, arrancando de un punto MALO
    #           (el escenario que rompe least_squares, 7.4)
    # =========================================================
    B_true, s_true, esc = 600.0, 40.0, 5000.0
    canales = np.linspace(0, 1000, 1024)
    y1 = _espectro_sintetico(canales, [(Eg1, esc, B_true, s_true)], semilla=1).astype(float)
    lo, hi = 270.0, 762.0                        # la ventana que eligio el NCC
    m = (canales >= lo) & (canales <= hi)
    xv, yv = canales[m], y1[m]

    # punto inicial deliberadamente malo: B 30% abajo, sigma 3x
    p0_malo = [esc * 0.5, 0.70 * B_true, 3.0 * s_true]
    cotas1 = ([0.0, lo, 1.0], [10 * esc, hi, hi - lo])

    mod1 = modelo_un_borde(Eg1)
    E_ls = costo_minimos_cuadrados(xv, yv, mod1)
    E_ncc = costo_ncc(xv, yv, mod1)

    loc = ajustar_borde_unico(canales, y1, Eg1, (lo, hi), p0=p0_malo)
    t0 = time.time()
    sa_ls = recocido_simulado(E_ls, p0_malo, cotas1, semilla=1)
    t_ls = time.time() - t0
    sa_nc = recocido_simulado(E_ncc, p0_malo, cotas1, semilla=1)
    p_nc = proyectar_amplitud(xv, yv, mod1, sa_nc["x"], 0)

    print(f"\n  PRUEBA 1: un borde (B={B_true:.0f}, sigma={s_true:.0f}), "
          f"arranque malo B0={p0_malo[1]:.0f} sigma0={p0_malo[2]:.0f}")
    print(f"    least_squares desde el mismo p0 : B = {loc['edge_ch']:6.1f}  "
          f"sigma = {loc['sigma']:5.1f}   (error B {loc['edge_ch']-B_true:+.1f})")
    print(f"    SA con costo Ec.14 (cuadrados)  : B = {sa_ls['x'][1]:6.1f}  "
          f"sigma = {sa_ls['x'][2]:5.1f}   (error B {sa_ls['x'][1]-B_true:+.1f})")
    print(f"    SA con costo NCC (1-rho)        : B = {p_nc[1]:6.1f}  "
          f"sigma = {p_nc[2]:5.1f}   (error B {p_nc[1]-B_true:+.1f})")
    print(f"    T_init auto = {sa_ls['T_inicial']:.4g} (con la Ec.14) vs "
          f"{sa_nc['T_inicial']:.4g} (con NCC)  <- por eso 'T=2000' no es universal")
    print(f"    {sa_ls['n_evals']} evaluaciones, {sa_ls['n_ciclos']} ciclos, {t_ls:.1f} s")

    # la amplitud es invisible para el NCC: se comprueba
    p_test = np.array(sa_nc["x"], float)
    e_a = E_ncc(p_test)
    p_test[0] *= 7.3
    print(f"    control: multiplicar A por 7.3 cambia el costo NCC en "
          f"{abs(E_ncc(p_test)-e_a):.2e}  -> A es INVISIBLE para el NCC")

    # =========================================================
    # PRUEBA 2: dos bordes (6 parametros), sintetico
    # =========================================================
    Emax1, Emax2 = compton_edge_energy(Eg1), compton_edge_energy(Eg2)
    B1_true = 600.0
    D_true = B1_true / Emax1
    B2_true = B1_true + D_true * delta_E()
    s1_true = 45.0
    s2_true = s1_true * np.sqrt(Emax2 / Emax1)
    A1_true, A2_true = 5000.0, 4500.0

    can2 = np.linspace(0, 1200, 1200)
    y2 = _espectro_sintetico(can2, [(Eg1, A1_true, B1_true, s1_true),
                                    (Eg2, A2_true, B2_true, s2_true)],
                             semilla=7).astype(float)
    iv2 = (0.45 * B1_true, 1.55 * B1_true)
    m2 = (can2 >= iv2[0]) & (can2 <= iv2[1])
    x2, yv2 = can2[m2], y2[m2]

    mod2 = modelo_dos(Eg1, Eg2)
    E2 = costo_minimos_cuadrados(x2, yv2, mod2)
    E4 = costo_dos_bordes_varpro(x2, yv2, Eg1, Eg2)   # solo [B1, s1, s2, D]

    # arranque malo otra vez: B1 20% abajo, sigmas al doble
    p0_6 = [3000.0, 3000.0, 0.80 * B1_true, 2 * s1_true, 2 * s2_true, 0.8 * D_true]
    cotas6 = ([0, 0, iv2[0], 1.0, 1.0, 0.2 * D_true],
              [5e4, 5e4, iv2[1], iv2[1] - iv2[0], iv2[1] - iv2[0], 3 * D_true])
    p0_4 = [0.80 * B1_true, 2 * s1_true, 2 * s2_true, 0.8 * D_true]
    cotas4 = ([iv2[0], 1.0, 1.0, 0.2 * D_true],
              [iv2[1], iv2[1] - iv2[0], iv2[1] - iv2[0], 3 * D_true])

    loc2_malo = ajustar_dos_bordes(can2, y2, iv2, p0=p0_6)
    loc2_bueno = ajustar_dos_bordes(can2, y2, iv2)          # con barrido de semilla
    e_ref = E2([loc2_bueno["A1"], loc2_bueno["A2"], loc2_bueno["B1"],
                loc2_bueno["sigma1"], loc2_bueno["sigma2"], loc2_bueno["D"]])

    print(f"\n  PRUEBA 2: dos bordes (B1={B1_true:.0f}, B2={B2_true:.0f}, "
          f"D={D_true:.4f}), arranque malo B1_0={p0_6[2]:.0f}")
    print(f"    {'metodo':40s} {'B1':>7s} {'D':>8s} {'err D':>8s} {'costo':>10s}")

    def linea(nom, B1, D, costo):
        print(f"    {nom:40s} {B1:7.1f} {D:8.4f} {100*(D/D_true-1):+7.1f}% {costo:10.4g}")

    linea("least_squares desde p0 malo", loc2_malo["B1"], loc2_malo["D"],
          E2([loc2_malo["A1"], loc2_malo["A2"], loc2_malo["B1"],
              loc2_malo["sigma1"], loc2_malo["sigma2"], loc2_malo["D"]]))
    linea("least_squares con semilla (Etapa 6)", loc2_bueno["B1"], loc2_bueno["D"], e_ref)

    t0 = time.time()
    sa_lit = recocido_simulado(E2, p0_6, cotas6, semilla=3, alpha=0.97,
                               max_evals=50000, paso_adaptativo=False)
    linea("SA 6 params, receta literal", sa_lit["x"][2], sa_lit["x"][5], sa_lit["costo"])

    # con las amplitudes proyectadas fuera, y varias semillas: el SA es estocastico
    for sem in (3, 11, 42):
        r = recocido_simulado(E4, p0_4, cotas4, semilla=sem, alpha=0.97,
                              max_evals=50000, n_estable=20)
        B1s, s1s, s2s, Ds = r["x"]
        from etapas.fit_double import componentes as _cmp, _amplitudes_optimas as _amp
        f1, f2, _ = _cmp(x2, B1s, s1s, s2s, Ds, Eg1, Eg2)
        A1s, A2s = _amp(f1, f2, yv2)
        pul = ajustar_dos_bordes(can2, y2, iv2, p0=[A1s, A2s, B1s, s1s, s2s, Ds])
        linea(f"SA 4p (VARPRO) + pulido, semilla={sem}", pul["B1"], pul["D"],
              E2([pul["A1"], pul["A2"], pul["B1"], pul["sigma1"], pul["sigma2"],
                  pul["D"]]))
    print(f"    (tiempo total prueba 2: {time.time()-t0:.0f} s)")

    print("\n    LECTURA:")
    print("    - Con los 6 parametros crudos (la receta literal) el SA NO converge:")
    print("      mezcla amplitudes (~5000) con D (~0.6) en el mismo vector de pasos.")
    print("    - Proyectando las amplitudes fuera (VARPRO, 4 parametros) + pulido")
    print("      local, el SA llega a un optimo MEJOR que el least_squares sembrado")
    print("      de la Etapa 6, en 2 de cada 3 semillas.")
    print("    - La tercera cae en un minimo local malo: el SA es ESTOCASTICO. Hay")
    print("      que correr varias semillas y quedarse con el de menor costo; una")
    print("      sola corrida no es un resultado.")
    print("    - La dispersion en D entre semillas viene de que los dos bordes se")
    print("      solapan dentro de ~1 sigma (CLAUDE.md 7.16): la superficie de costo")
    print("      es casi plana en esa direccion.")

    # =========================================================
    # PRUEBA 3: los DATOS REALES de Co-60
    # =========================================================
    print("\n  PRUEBA 3: datos reales de Co-60 (multiarranque, 5 semillas)")
    from datos_io import spectrum
    from etapas.fit_double import componentes as cmp3, _amplitudes_optimas as amp3

    rr = spectrum.espectro_neto(rutas.MUESTRA,
                                rutas.FONDO,
                                umbral_ruido=0.05, num_canales=800,
                                rango_max=8000.0, recortar_negativos=False)
    canR, netoR, errR = rr["canales"], rr["neto"], rr["error"]
    ivR = (1960.0, 7886.0)                    # el intervalo que eligio el NCC
    mR = (canR >= ivR[0]) & (canR <= ivR[1])
    xR, yR, eR = canR[mR], netoR[mR], errR[mR]

    ER = costo_dos_bordes_varpro(xR, yR, Eg1, Eg2, errores=eR)
    p0_R = [3000.0, 600.0, 600.0, 3000.0 / Emax1]
    cotasR = ([ivR[0], 50.0, 50.0, 0.2 * 3000.0 / Emax1],
              [ivR[1], 0.5 * (ivR[1] - ivR[0]), 0.5 * (ivR[1] - ivR[0]),
               3.0 * 3000.0 / Emax1])

    t0 = time.time()
    saR = sa_multiarranque(ER, p0_R, cotasR, semillas=(3, 11, 42, 7, 99),
                           alpha=0.97, max_evals=40000, n_estable=20)
    B1R, s1R, s2R, DR = saR["x"]
    f1R, f2R, ce2R = cmp3(xR, B1R, s1R, s2R, DR, Eg1, Eg2)
    A1R, A2R = amp3(f1R, f2R, yR)
    pulR = ajustar_dos_bordes(canR, netoR, ivR, p0=[A1R, A2R, B1R, s1R, s2R, DR])

    print(f"    SA (mejor de 5, semilla {saR['semilla']}): B1 = {B1R:.0f} mV*ns   "
          f"D = {DR:.3f} mV*ns/keV   chi2 = {saR['costo']:.4g}")
    print(f"    + pulido local            : B1 = {pulR['B1']:.0f} mV*ns   "
          f"D = {pulR['D']:.3f} mV*ns/keV")
    print(f"    Etapa 6 (barrido+local)   : B1 = 4179 mV*ns   D = 4.437 mV*ns/keV")
    print(f"    Etapa 6 (origen ligado)   : B1 = 4185 mV*ns   D = 4.344 mV*ns/keV")
    print(f"    forward folding + Geant4  :                   D = 4.23  mV*ns/keVee")
    print(f"    ({time.time()-t0:.0f} s)")

    # ------------------------------ figura ------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

    h = sa_ls["hist"]
    ax[0].plot(h["mejor"], "b", lw=1.8, label="mejor costo")
    ax[0].plot(h["actual"], color="0.6", lw=0.8, alpha=0.8, label="costo actual")
    ax[0].set_yscale("log"); ax[0].set_xlabel("ciclo de temperatura")
    ax[0].set_ylabel("E(x)  (Ec. 14)")
    ax[0].set_title("SA: el costo actual SUBE a veces\n(por eso escapa minimos locales)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    axb = ax[1]
    axb.plot(h["T"], color="crimson", lw=1.8)
    axb.set_yscale("log"); axb.set_xlabel("ciclo de temperatura")
    axb.set_ylabel("T", color="crimson")
    axb2 = axb.twinx()
    axb2.plot(h["aceptacion"], color="steelblue", lw=1.2)
    axb2.set_ylabel("tasa de aceptacion", color="steelblue")
    axb.set_title("Enfriamiento geometrico (alpha=0.95)\ny caida de la aceptacion")
    axb.grid(alpha=0.3)

    xf = np.linspace(lo, hi, 500)
    ax[2].step(canales, y1, where="mid", color="0.75", lw=1, label="datos")
    ax[2].plot(xf, mod1(xf, sa_ls["x"]), "b", lw=2,
               label=f"SA (B={sa_ls['x'][1]:.0f})")
    ax[2].plot(xf, mod1(xf, [loc["A"], loc["edge_ch"], loc["sigma"]]), "r--", lw=2,
               label=f"least_squares (B={loc['edge_ch']:.0f})")
    ax[2].axvline(B_true, color="k", ls="--", alpha=0.7, label=f"B real = {B_true:.0f}")
    ax[2].set_xlim(200, 900); ax[2].set_xlabel("Canal"); ax[2].set_ylabel("Cuentas")
    ax[2].set_title("Mismo punto inicial malo:\nSA llega, el local no")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("sa_demo.png"), dpi=130)
    print("\nFigura guardada: figuras/sa_demo.png")
