"""
ncc.py - Etapa 4: seleccion AUTOMATICA del intervalo del borde (Seccion II-B, Ec. 15).

El problema que resuelve
------------------------
En la Etapa 3 TU elegias el intervalo a mano (`lo,hi` por teclado). Eso es
subjetivo y ademas el ajuste local es sensible a esa eleccion (ver CLAUDE.md
7.4). El paper propone elegir el intervalo por CORRELACION DE FORMA:

  S2 : 100 puntos de la lineshape TEORICA, del plateau al minimo tras el borde.
  S1 : cada intervalo candidato del espectro MEDIDO, reinterpolado a 100 puntos.
  rho: coeficiente de Pearson entre S1 y S2 (Ec. 15). Gana el rho maximo.

Por que Pearson y no una distancia
----------------------------------
Pearson resta la media y divide por la desviacion de CADA senal: es inmune a
offset (fondo residual) y a escala (amplitud A, tiempo de medicion). Solo mira
la FORMA. Y como ambos lados se reinterpolan a 100 puntos, tambien es inmune a
la escala HORIZONTAL: por eso puede encontrar el borde sin saber la ganancia,
que es exactamente lo que queremos calibrar.

Bonus [fuera del paper]
-----------------------
Como el intervalo ganador queda mapeado punto a punto contra la plantilla, el
mapa lineal canal<->energia sale GRATIS: de ahi salen estimaciones de B (canal
del borde) y de sigma para SEMBRAR el ajuste. Eso mata la fragilidad de 7.4
sin necesidad todavia del SA (Etapa 5).

Corre solo:  python ncc.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/ncc.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
from etapas.physics import compton_edge_energy
from etapas.model import forward_model


# ----------------------------------------------------------------------
# 1) S2: la plantilla teorica (100 puntos de referencia)
# ----------------------------------------------------------------------
def plantilla_referencia(E_gamma, sigma_rel=0.06, n_puntos=100,
                         inicio_rel=0.45, fin_rel=None, n_fino=4000):
    """
    Devuelve (E_ref, S2) con n_puntos muestreados uniformemente sobre la
    lineshape teorica CONVOLUCIONADA, desde el plateau hasta despues del borde.

    sigma_rel : ancho gaussiano EN UNIDADES DE E_max (adimensional).
                [fuera del paper] El paper no dice con que sigma se construye
                la plantilla. Se parametriza relativo porque asi la plantilla
                es independiente de la ganancia del detector: la misma S2 sirve
                para cualquier eje de canales.
    inicio_rel: donde empieza la plantilla, en unidades de E_max (0.45 = plateau).
    fin_rel   : donde termina; por defecto 4.5*sigma sobre el borde (el "edge
                minimum" del paper: donde la cola ya cayo a ~0).

    Ojo: se convoluciona sobre una grilla que arranca en E=0 y DESPUES se
    recorta, para que el borde izquierdo de la convolucion no contamine S2.
    """
    E_max = compton_edge_energy(E_gamma)
    if fin_rel is None:
        fin_rel = max(1.15, 1.0 + 4.5 * sigma_rel)

    fino = np.linspace(0.0, fin_rel * E_max, n_fino)
    # forward_model con edge_ch = E_max => el mapa canal->energia es la identidad
    y = forward_model(fino, E_gamma, A=1.0, edge_ch=E_max, sigma=sigma_rel * E_max)

    E_ref = np.linspace(inicio_rel * E_max, fin_rel * E_max, n_puntos)
    return E_ref, np.interp(E_ref, fino, y)


def plantilla_simulada(centros, pesos, sigma_rel=0.20, n_puntos=100,
                       inicio_rel=0.45, fin_rel=None, E_ref=963.4):
    """
    S2 construida con la SIMULACION Geant4 en vez de Klein-Nishina. [fuera del paper]

    Por que: la lineshape analitica es scattering simple en un medio infinito.
    La simulacion ya trae el scattering multiple, el escape, el quenching de
    Birks y las DOS lineas del Co-60 con su intensidad relativa (CLAUDE.md 7.18).
    Si la ventana se va a usar para comparar datos contra simulacion, la
    plantilla que elige la ventana tiene que ser la simulacion.

    Truco de reuso: `modelo_simulado` con D = 1 deja el eje de "canal" en keV,
    asi que la misma funcion del forward folding sirve de plantilla.

    sigma_rel : igual que en plantilla_referencia, sigma(E_ref)/E_ref. Se traduce
                a la parametrizacion de la simulacion como alpha = sigma_rel*sqrt(E_ref).
    E_ref     : energia a la que corresponde el parametro B del barrido. Con la
                plantilla simulada NO hay un unico borde (hay dos), asi que E_ref
                es solo el ancla del mapa por el origen: D = B / E_ref.
    fin_rel   : por defecto llega hasta el SEGUNDO borde mas 3 sigma
                (1118.1/963.4 + 3*sigma_rel): si la ventana corta antes, deja
                fuera justo la parte que distingue una plantilla de otra.
    """
    from ajuste.sim_template import modelo_simulado
    if fin_rel is None:
        fin_rel = fin_rel_simulada(sigma_rel)
    alpha = sigma_rel * np.sqrt(E_ref)
    E_grid = np.linspace(inicio_rel * E_ref, fin_rel * E_ref, n_puntos)
    return E_grid, modelo_simulado(E_grid, centros, pesos, A=1.0, D=1.0, alpha=alpha)


def fin_rel_simulada(sigma_rel, E1=963.4, E2=1118.1):
    """Donde termina la ventana con plantilla simulada: 3 sigma sobre el 2do borde."""
    return E2 / E1 + 3.0 * sigma_rel


def fin_rel_kn(sigma_rel):
    """Donde termina la ventana con Klein-Nishina: 4.5 sigma sobre el unico borde."""
    return max(1.15, 1.0 + 4.5 * sigma_rel)


# ----------------------------------------------------------------------
# 2) Pearson (Ec. 15)
# ----------------------------------------------------------------------
def pearson(s1, s2):
    """Coeficiente de correlacion de Pearson, Ec. 15. rho en [-1, 1]."""
    a = np.asarray(s1, float) - np.mean(s1)
    b = np.asarray(s2, float) - np.mean(s2)
    den = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.dot(a, b) / den) if den > 0 else -1.0


def _z(v):
    """Estandariza un vector (media 0, desv 1). Pearson = promedio de z1*z2."""
    v = np.asarray(v, float)
    sd = v.std()
    return (v - v.mean()) / sd if sd > 0 else np.zeros_like(v)


# ----------------------------------------------------------------------
# 3) Barrido de intervalos: el corazon de la Etapa 4
# ----------------------------------------------------------------------
def _barrido_un_sigma(cuentas, z2, inicios, anchos, n_puntos, acum=None, min_cuentas=0.0):
    """
    Para cada (inicio, ancho) reinterpola la ventana a n_puntos y correlaciona.

    Truco de vectorizacion: para un ancho FIJO, todas las ventanas comparten los
    mismos desplazamientos fraccionarios respecto de su inicio. Se arma entonces
    una matriz (n_inicios x n_puntos) de un tiro con indexado fancy, en vez de
    llamar np.interp una vez por ventana. Requiere canales UNIFORMES (lo son:
    son bins de histograma, ver CLAUDE.md 6).
    """
    N = len(cuentas)
    mejor = (-2.0, None, None)          # (rho, inicio, ancho)
    rho_por_inicio = np.full(N, -2.0)   # diagnostico: mejor rho de cada inicio

    for w in anchos:
        ini = inicios[inicios + w <= N - 1]
        if ini.size == 0:
            continue
        # posiciones (fraccionarias) de los n_puntos dentro de una ventana de ancho w
        off = np.linspace(0.0, w, n_puntos)
        i0 = np.floor(off).astype(int)
        f = off - i0
        i1 = np.minimum(i0 + 1, w)

        idx0 = ini[:, None] + i0[None, :]
        idx1 = ini[:, None] + i1[None, :]
        M = cuentas[idx0] * (1.0 - f) + cuentas[idx1] * f   # S1 de cada ventana

        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, keepdims=True)
        ok = sd[:, 0] > 0
        Z = np.zeros_like(M)
        Z[ok] = (M[ok] - mu[ok]) / sd[ok]

        rho = Z @ z2 / n_puntos            # Pearson vectorizado, Ec. 15
        rho[~ok] = -2.0

        if acum is not None and min_cuentas > 0:
            # Pearson ignora la amplitud -> sin este filtro gana una ventana de
            # la cola vacia, donde el z-score amplifica ruido puro.
            en_ventana = acum[ini + w + 1] - acum[ini]
            rho[en_ventana < min_cuentas] = -2.0

        np.maximum.at(rho_por_inicio, ini, rho)
        k = int(np.argmax(rho))
        if rho[k] > mejor[0]:
            mejor = (float(rho[k]), int(ini[k]), int(w))

    return mejor, rho_por_inicio


def _barrido_origen(canales, cuentas, z2, B_cands, inicio_rel, fin_rel, n_puntos):
    """
    Barrido de UN parametro: el canal del borde B. La ventana no es libre, queda
    amarrada a B por  [inicio_rel*B, fin_rel*B].

    Por que: el forward model (Ec. 8) usa E_eff = E_e * B/E_max, o sea un mapa
    canal<->energia POR EL ORIGEN (area cero <-> energia cero). Si el mapa pasa
    por el origen, fijado B la ventana ya esta determinada: no hay dos extremos
    libres. El barrido 2D del paper esta sobre-parametrizado respecto de su
    propio modelo, y sobre datos ruidosos eso se paga (ver CLAUDE.md 7.10).
    """
    mejor = (-2.0, None)
    rho_de_B = np.full(len(B_cands), -2.0)
    for k, B in enumerate(B_cands):
        lo, hi = inicio_rel * B, fin_rel * B
        if lo < canales[0] or hi > canales[-1]:
            continue
        s1 = np.interp(np.linspace(lo, hi, n_puntos), canales, cuentas)
        sd = s1.std()
        if sd <= 0:
            continue
        rho = float(np.dot((s1 - s1.mean()) / sd, z2) / n_puntos)
        rho_de_B[k] = rho
        if rho > mejor[0]:
            mejor = (rho, float(B))
    return mejor, rho_de_B


def buscar_intervalo(canales, cuentas, E_gamma,
                     sigma_rel=(0.03, 0.045, 0.06, 0.09, 0.12, 0.14, 0.16, 0.18, 0.20),
                     n_puntos=100, inicio_rel=0.45, region=None, modo="origen",
                     ancho_min=None, ancho_max=None, paso=2, umbral_cuentas=0.0,
                     cuentas_crudas=None, min_cuentas=200.0,
                     plantilla=None, E_ref=None, fin_rel_fn=None,
                     modelo_semilla=None):
    """
    Selecciona automaticamente el intervalo del Compton edge (Seccion II-B).

    canales, cuentas : espectro medido (canales uniformes). Conviene SUAVIZAR
                       antes con spectrum.suavizar (Savitzky-Golay).
    E_gamma          : linea gamma [keV] cuyo borde se busca.
    sigma_rel        : escalar o lista. Si es lista se barre y gana el rho mas
                       alto global -> sigma queda estimado por los datos
                       [fuera del paper; el paper no fija este valor].
    n_puntos         : cuantos puntos se reinterpolan de cada ventana. "auto" =
                       UN PUNTO POR BIN de la ventana tipica, o sea sin inventar
                       puntos. Medido (CLAUDE.md 7.27): de 100 a 571 puntos B se
                       mueve 0.4% y se estabiliza; interpolar a POCOS puntos infla
                       rho al promediar ruido. El valor queda FIJO durante todo el
                       barrido: si n cambiara con B, los rho de ventanas distintas
                       dejarian de ser comparables.
    region           : (lo, hi) para limitar la busqueda (util en datos reales
                       para excluir el pico de ruido de bajo canal).
    modo             : "origen" (recomendado) barre solo B con la ventana amarrada
                       por el mapa por el origen del propio forward model.
                       "libre" replica el barrido 2D (inicio, fin) del paper.
    paso             : stride del barrido en bins ("libre") o paso de B en bins
                       ("origen").
    cuentas_crudas   : espectro SIN suavizar; si se da, en modo "libre" se
                       descartan ventanas con menos de min_cuentas cuentas
                       (Pearson ignora la amplitud: sin este filtro gana la cola
                       vacia, donde el z-score amplifica puro ruido).

    plantilla        : callable(sigma_rel, n_puntos, inicio_rel, fin_rel) -> S2.
                       Por defecto la lineshape Klein-Nishina (el paper). Pasarle
                       `plantilla_simulada` amarrada a la simulacion cambia la
                       referencia a Geant4 sin tocar nada mas del barrido.
    E_ref            : energia [keV] a la que corresponde B. Default = E_max del
                       borde. Fija el mapa por el origen: D_seed = B / E_ref.
    fin_rel_fn       : callable(sigma_rel) -> fin_rel. Default `fin_rel_kn`; con
                       plantilla simulada usar `fin_rel_simulada`.
    modelo_semilla   : callable(x, edge_ch, sigma_ch) -> curva sin amplitud, solo
                       para escalar A0. Default el forward model analitico.

    Devuelve dict con: lo, hi, rho, edge_ch (B estimado), sigma (estimado en
    canales), p0 (semilla [A, B, sigma] lista para fit_single), p0_sim (semilla
    [A, D, alpha] lista para sim_template), D_seed, alpha_seed, y diagnosticos.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    N = len(canales)
    paso_ch = float(np.mean(np.diff(canales)))
    sigmas = np.atleast_1d(np.asarray(sigma_rel, float))
    lo_reg = region[0] if region is not None else canales[0]
    hi_reg = region[1] if region is not None else canales[-1]

    if ancho_min is None:
        ancho_min = max(20, int(0.04 * N))
    if ancho_max is None:
        ancho_max = int(0.90 * N)

    if E_ref is None:
        E_ref = compton_edge_energy(E_gamma)
    if fin_rel_fn is None:
        fin_rel_fn = fin_rel_kn
    if plantilla is None:
        def plantilla(s_rel, n_pts, ini, fin):
            return plantilla_referencia(E_gamma, sigma_rel=s_rel, n_puntos=n_pts,
                                        inicio_rel=ini, fin_rel=fin)[1]
    if modelo_semilla is None:
        def modelo_semilla(x, edge_ch, sigma_ch):
            return forward_model(x, E_gamma, A=1.0, edge_ch=edge_ch, sigma=sigma_ch)

    # --- n_puntos = "auto": un punto por bin de la ventana TIPICA -------------
    # Se calcula UNA vez (no por ventana) para que los rho sigan siendo
    # comparables entre si durante el barrido.
    if isinstance(n_puntos, str):
        s_med = float(np.median(sigmas))
        ancho_rel = fin_rel_fn(s_med) - inicio_rel
        if modo == "origen":
            B_tip = 0.5 * (lo_reg / inicio_rel + hi_reg / fin_rel_fn(s_med))
            n_puntos = int(round(ancho_rel * B_tip / paso_ch))
        else:
            n_puntos = int(round(0.5 * (ancho_min + ancho_max)))
        n_puntos = int(min(max(n_puntos, 20), 2000))

    mejor_global = (-2.0, None, None, None, None)   # rho, lo, hi, s_rel, fin_rel
    diag = {}

    if modo == "origen":
        for s_rel in sigmas:
            fin_rel = fin_rel_fn(s_rel)
            z2 = _z(plantilla(s_rel, n_puntos, inicio_rel, fin_rel))
            # B admisibles: la ventana entera debe caber en la region pedida
            B_min = max(lo_reg / inicio_rel,
                        ancho_min * paso_ch / (fin_rel - inicio_rel))
            B_max = hi_reg / fin_rel
            if B_max <= B_min:
                continue
            B_cands = np.arange(B_min, B_max, max(paso, 1) * paso_ch * 0.5)
            (rho, B), rho_de_B = _barrido_origen(canales, cuentas, z2, B_cands,
                                                 inicio_rel, fin_rel, n_puntos)
            if B is not None and rho > mejor_global[0]:
                mejor_global = (rho, inicio_rel * B, fin_rel * B, float(s_rel), fin_rel)
                diag = {"B_cands": B_cands, "rho_de_B": rho_de_B}

    elif modo == "libre":
        valido = cuentas > umbral_cuentas
        valido &= (canales >= lo_reg) & (canales <= hi_reg)
        inicios = np.nonzero(valido)[0][::paso]
        if inicios.size == 0:
            raise ValueError("Ningun canal candidato: revisa region/umbral_cuentas.")
        tope = int(np.searchsorted(canales, hi_reg))
        tope = min(max(tope, 1), N - 1)
        anchos = np.arange(ancho_min, min(ancho_max, tope) + 1, paso)

        # filtro de estadistica: acumulada de cuentas crudas para descartar
        # ventanas vacias en O(1) por ventana
        acum = None
        if cuentas_crudas is not None:
            acum = np.concatenate([[0.0], np.cumsum(np.asarray(cuentas_crudas, float))])

        rho_mapa = np.full(N, -2.0)
        for s_rel in sigmas:
            fin_rel = fin_rel_fn(s_rel)
            S2 = plantilla(s_rel, n_puntos, inicio_rel, fin_rel)
            (rho, i, w), rpi = _barrido_un_sigma(cuentas, _z(S2), inicios, anchos,
                                                 n_puntos, acum, min_cuentas)
            rho_mapa = np.maximum(rho_mapa, rpi)
            if i is not None and rho > mejor_global[0]:
                mejor_global = (rho, canales[i], canales[i + w], float(s_rel), fin_rel)
        diag = {"rho_por_inicio": rho_mapa}
    else:
        raise ValueError("modo debe ser 'origen' o 'libre'")

    rho, lo, hi, s_rel, fin_rel = mejor_global
    if lo is None:
        raise ValueError("El barrido no encontro ningun intervalo valido.")

    # --- el mapa canal<->energia que regala el intervalo ganador --------------
    # el intervalo [lo,hi] quedo alineado con [inicio_rel, fin_rel]*E_max,
    # asi que E_max (el borde) cae en esta fraccion del intervalo:
    frac_borde = (1.0 - inicio_rel) / (fin_rel - inicio_rel)
    edge_ch = lo + frac_borde * (hi - lo)
    sigma_ch = s_rel * edge_ch                    # sigma_rel era relativo a E_max

    # semilla de amplitud: iguala el maximo del modelo al maximo de los datos
    m = (canales >= lo) & (canales <= hi)
    forma = modelo_semilla(canales[m], edge_ch, sigma_ch)
    pico = np.max(forma)
    A0 = (np.max(cuentas[m]) / pico) if pico > 0 else 1.0

    # el mapa por el origen en las coordenadas del forward folding:
    #   D = B / E_ref   y   alpha = sigma_E(E_ref)/sqrt(E_ref) = sigma_rel*sqrt(E_ref)
    D_seed = edge_ch / E_ref
    alpha_seed = s_rel * np.sqrt(E_ref)

    i_ini = int(np.searchsorted(canales, lo))
    salida = {
        "lo": float(lo), "hi": float(hi), "rho": float(rho),
        "edge_ch": float(edge_ch), "sigma": float(sigma_ch),
        "sigma_rel": s_rel, "fin_rel": fin_rel, "inicio_rel": inicio_rel,
        "intervalo": (float(lo), float(hi)),
        "p0": [float(A0), float(edge_ch), float(sigma_ch)],
        "p0_sim": [float(A0), float(D_seed), float(alpha_seed)],
        "D_seed": float(D_seed), "alpha_seed": float(alpha_seed),
        "E_ref": float(E_ref), "n_puntos": int(n_puntos),
        "E_max": compton_edge_energy(E_gamma), "modo": modo,
        "i_ini": i_ini, "ancho_bins": int(np.searchsorted(canales, hi)) - i_ini,
        "paso_ch": paso_ch,
    }
    salida.update(diag)
    return salida


# grilla de sigma_rel para la ruta simulada: la resolucion real de este detector
# es alta (sigma_rel ~ 0.2), asi que la grilla del paper (que llega a 0.20) queda
# corta y el maximo de rho cae en el borde.
GRILLA_SIGMA_SIM = (0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28)


def buscar_intervalo_simulado(canales, cuentas, centros, pesos, E_gamma=1173.2,
                              E_ref=963.4, sigma_rel=GRILLA_SIGMA_SIM,
                              E_ref2=1118.1, **kw):
    """
    `buscar_intervalo` con la plantilla de la SIMULACION. [fuera del paper]

    Es el mismo barrido de un solo parametro (modo origen), cambiando la
    referencia: en vez de preguntar "donde se parece esto a Klein-Nishina",
    pregunta "donde se parece esto a mi Geant4". Como la ventana es la que
    despues se usa para el forward folding, la plantilla que la elige tiene que
    ser la misma que se va a ajustar.

    E_ref2 : borde de la SEGUNDA linea de la fuente. La ventana tiene que
             pasarlo, si no corta justo la parte que distingue una plantilla de
             otra (CLAUDE.md 7.28). Con una sola linea (Cs-137) hay que pasar
             None: ahi no hay segundo borde y se usa el mismo criterio que
             Klein-Nishina, el borde mas 4.5 sigma.

    Devuelve el mismo dict, con `D_seed` y `p0_sim` ya listos para sim_template.
    """
    from ajuste.sim_template import modelo_simulado

    def plantilla(s_rel, n_pts, ini, fin):
        return plantilla_simulada(centros, pesos, sigma_rel=s_rel, n_puntos=n_pts,
                                  inicio_rel=ini, fin_rel=fin, E_ref=E_ref)[1]

    def semilla(x, edge_ch, sigma_ch):
        D = edge_ch / E_ref
        alpha = (sigma_ch / max(edge_ch, 1e-9)) * np.sqrt(E_ref)
        return modelo_simulado(x, centros, pesos, 1.0, D, alpha)

    if E_ref2 and E_ref2 > E_ref:
        def fin_fn(s_rel, E1=E_ref, E2=E_ref2):
            return fin_rel_simulada(s_rel, E1=E1, E2=E2)
    else:
        fin_fn = fin_rel_kn          # una sola linea: no hay segundo borde

    return buscar_intervalo(canales, cuentas, E_gamma, sigma_rel=sigma_rel,
                            plantilla=plantilla, E_ref=E_ref,
                            fin_rel_fn=fin_fn, modelo_semilla=semilla,
                            **kw)


# ======================================================================
#  DEMO / VALIDACION con verdad conocida (asi trabajamos: sintetico primero)
# ======================================================================
if __name__ == "__main__":
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from etapas.fit_single import _espectro_sintetico, ajustar_borde_unico

    Eg = 1173.2                      # linea del Co-60 -> borde 963.4 keV
    B_true, s_true, esc = 600.0, 40.0, 5000.0
    canales = np.linspace(0, 1000, 1024)
    y = _espectro_sintetico(canales, [(Eg, esc, B_true, s_true)], semilla=1).astype(float)

    print("=" * 68)
    print("ETAPA 4 - NCC: seleccion automatica del intervalo (verdad conocida)")
    print("=" * 68)
    print(f"  B verdadero = {B_true:.1f} canales   |   sigma verdadero = {s_true:.1f}"
          f"   (sigma/B = {s_true/B_true:.3f})")

    t0 = time.time()
    r = buscar_intervalo(canales, y, Eg, paso=2)
    dt = time.time() - t0

    print(f"\n  intervalo elegido  = [{r['lo']:.1f}, {r['hi']:.1f}]  canales")
    print(f"  rho maximo         = {r['rho']:.4f}   (sigma_rel ganador = {r['sigma_rel']:.3f})")
    print(f"  B estimado por NCC = {r['edge_ch']:.1f}   (error {r['edge_ch']-B_true:+.1f})")
    print(f"  sigma estimado     = {r['sigma']:.1f}   (error {r['sigma']-s_true:+.1f})")
    print(f"  barrido en {dt:.2f} s")

    # ---- el punto del ejercicio: la semilla del NCC arregla la fragilidad 7.4
    fit_ncc = ajustar_borde_unico(canales, y, Eg, r["intervalo"], p0=r["p0"])
    fit_malo = ajustar_borde_unico(canales, y, Eg, (450.0, 700.0))   # intervalo descentrado

    print("\n  --- ajuste sembrado por NCC vs intervalo malo a mano ---")
    print(f"  fit con NCC   : B = {fit_ncc['edge_ch']:7.1f}  sigma = {fit_ncc['sigma']:6.1f}"
          f"   (error B {fit_ncc['edge_ch']-B_true:+.1f})")
    print(f"  fit a mano    : B = {fit_malo['edge_ch']:7.1f}  sigma = {fit_malo['sigma']:6.1f}"
          f"   (error B {fit_malo['edge_ch']-B_true:+.1f})")

    # ---- sensibilidad a sigma_rel (cuanto importa la eleccion de la plantilla)
    print("\n  --- sensibilidad: forzando cada sigma_rel por separado ---")
    for s_rel in (0.03, 0.045, 0.06, 0.09, 0.12):
        rr = buscar_intervalo(canales, y, Eg, sigma_rel=s_rel, paso=4)
        print(f"    sigma_rel={s_rel:5.3f} -> rho={rr['rho']:.4f}  B={rr['edge_ch']:6.1f}"
              f"  ({rr['edge_ch']-B_true:+.1f})")

    # ------------------------------ figura ------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    E_ref, S2 = plantilla_referencia(Eg, sigma_rel=r["sigma_rel"])
    ax[0, 0].plot(E_ref, S2, "g", lw=2)
    ax[0, 0].axvline(r["E_max"], color="k", ls="--", label=f"E_max = {r['E_max']:.0f} keV")
    ax[0, 0].set_title(f"S2: plantilla de 100 puntos (sigma_rel = {r['sigma_rel']:.3f})")
    ax[0, 0].set_xlabel("Energia del electron [keV]"); ax[0, 0].set_ylabel("u.a.")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)

    ax[0, 1].step(canales, y, where="mid", color="0.6", label="datos sinteticos")
    ax[0, 1].axvspan(r["lo"], r["hi"], color="green", alpha=0.15,
                     label=f"intervalo NCC (rho={r['rho']:.3f})")
    ax[0, 1].axvline(B_true, color="k", ls="--", label=f"B real = {B_true:.0f}")
    ax[0, 1].axvline(r["edge_ch"], color="g", ls=":", lw=2,
                     label=f"B por NCC = {r['edge_ch']:.0f}")
    ax[0, 1].set_xlim(200, 1000); ax[0, 1].set_title("Intervalo elegido AUTOMATICAMENTE")
    ax[0, 1].set_xlabel("Canal"); ax[0, 1].set_ylabel("Cuentas")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)

    # S1 ganador vs S2, ambos estandarizados: esto es LO QUE COMPARA Pearson
    i, w = r["i_ini"], r["ancho_bins"]
    xq = np.linspace(canales[i], canales[i + w], 100)
    S1 = np.interp(xq, canales, y)
    ax[1, 0].plot(np.arange(100), _z(S1), color="0.4", lw=1.2, label="S1 (datos, z)")
    ax[1, 0].plot(np.arange(100), _z(S2), "g", lw=2, label="S2 (teoria, z)")
    ax[1, 0].set_title(f"Lo que compara Pearson: rho = {r['rho']:.4f}")
    ax[1, 0].set_xlabel("punto reinterpolado (1..100)"); ax[1, 0].set_ylabel("z-score")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)

    xf = np.linspace(r["lo"], r["hi"], 400)
    ax[1, 1].step(canales, y, where="mid", color="0.6", label="datos")
    ax[1, 1].plot(xf, forward_model(xf, Eg, fit_ncc["A"], fit_ncc["edge_ch"], fit_ncc["sigma"]),
                  "b", lw=2, label=f"fit sembrado NCC (B={fit_ncc['edge_ch']:.0f})")
    xf2 = np.linspace(450, 700, 400)
    ax[1, 1].plot(xf2, forward_model(xf2, Eg, fit_malo["A"], fit_malo["edge_ch"], fit_malo["sigma"]),
                  "r", lw=2, ls="--", label=f"fit intervalo malo (B={fit_malo['edge_ch']:.0f})")
    ax[1, 1].axvline(B_true, color="k", ls="--", label=f"B real = {B_true:.0f}")
    ax[1, 1].set_xlim(200, 1000); ax[1, 1].set_title("La semilla del NCC arregla el ajuste local")
    ax[1, 1].set_xlabel("Canal"); ax[1, 1].set_ylabel("Cuentas")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("ncc_demo.png"), dpi=130)
    print("\nFigura guardada: figuras/ncc_demo.png")
