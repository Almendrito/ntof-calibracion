"""
modelo_ext.py - Etapa 9: el modelo extendido que convierte los SISTEMATICOS en
parametros ajustados.

De donde sale
-------------
La Etapa 8 entregaba D con un error dominado por dos cosas que NO eran
estadistica y que se trataban a mano:

  1. la PLANTILLA: usar la simulacion completa o solo el scattering simple
     cambiaba D un 5%, y se reportaba esa diferencia como sistematico. Pero
     "cuanto scattering multiple hay" no es una opinion: es un numero que
     depende de la geometria real (tamano del centellador, material de
     alrededor) y que los DATOS pueden medir.
  2. el EXCESO A BAJO CANAL: cuentas que la simulacion no tiene. Se resolvia
     CORTANDO la ventana por arriba del exceso, lo que tira estadistica y deja
     una dependencia con donde se corta.

Aca las dos se meten dentro del modelo:

    S(ch) = A * [ w*h_1sc + (1-w)*h_multi ] (x) G(sigma_E)  +  C*exp(-ch/lambda)
            |______ Geant4 con la mezcla ajustable ____|         |__ exceso __|

y ademas la resolucion se generaliza a los tres terminos habituales

    sigma_E(E) = sqrt( a^2 E^2  +  alpha^2 E  +  beta^2 )
                       ruido de   estadistica   ruido
                       colecta    de fotoelec.  electronico

de los cuales el modelo viejo solo tenia alpha.

Que se gana: el error de D deja de ser "la dispersion entre tres ajustes que
no se pueden comparar" y pasa a ser el ancho MARGINAL del posterior, que ya
incluye la ignorancia sobre w y sobre el exceso. Eso si es propagable.

Que NO se gana: si la simulacion no contiene un proceso, w no lo inventa. w
mide cuanto pesa el scattering multiple SIMULADO, no valida la geometria.

Requisito previo (importante)
-----------------------------
Todo esto corre sobre las areas del integrador con compuerta
(`spectrum._area_pulso_gate`) y la normalizacion por tiempo vivo
(`spectrum.tiempo_vivo`). Sobre las areas viejas no tiene sentido: ahi el
observable mismo estaba roto (CLAUDE.md 7.33-7.34).

Corre solo:  python modelo_ext.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/ajuste/modelo_ext.py` pone en sys.path la carpeta DEL ARCHIVO, no
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

ARCHIVO_SIM = rutas.SIM

# nombres y orden del vector de parametros
NOMBRES = ["A", "D", "alpha", "w", "C", "lambda", "beta", "a_rel"]
I_A, I_D, I_ALPHA, I_W, I_C, I_LAM, I_BETA, I_AREL = range(8)


# ----------------------------------------------------------------------
# 1) las dos componentes de la simulacion, por separado
# ----------------------------------------------------------------------
def cargar_componentes(ruta=ARCHIVO_SIM, rama="Lvis_keVee", bins=175,
                       rango=(0.0, 1400.0), arbol="ComptonSpectrum"):
    """
    Devuelve (centros, h_simple, h_multiple), normalizados al MISMO total.

    h_simple : eventos con NCompton == 1 (un solo scattering; es lo que
               describe Klein-Nishina, o sea la plantilla del paper).
    h_multiple: el resto (>=2 dispersiones, o dispersion + fotoelectrico). Son
               los que rellenan la zona POR ENCIMA del borde.

    Se normalizan por separado a area 1 para que `w` sea directamente la
    FRACCION de eventos de scattering simple, y no quede enredada con A.
    """
    import uproot
    with uproot.open(ruta) as f:
        t = f[arbol]
        d = t.arrays([rama, "NCompton"], library="np")
    x, nc = d[rama], d["NCompton"]

    simple = x[nc == 1]
    multi = x[nc != 1]
    h1, bordes = np.histogram(simple, bins=bins, range=rango)
    hm, _ = np.histogram(multi, bins=bins, range=rango)
    centros = 0.5 * (bordes[:-1] + bordes[1:])

    h1 = h1.astype(float)
    hm = hm.astype(float)
    w_nat = h1.sum() / max(h1.sum() + hm.sum(), 1.0)   # la mezcla de la simulacion
    h1 /= max(h1.sum(), 1.0)
    hm /= max(hm.sum(), 1.0)
    return centros, h1, hm, float(w_nat)


# ----------------------------------------------------------------------
# 2) el modelo
# ----------------------------------------------------------------------
def sigma_energia(E, alpha, beta=0.0, a_rel=0.0):
    """sigma_E(E) = sqrt(a^2 E^2 + alpha^2 E + beta^2)  [keV].

    - a_rel [adimensional] : parte que NO mejora con la estadistica de luz
      (no uniformidad de la colecta, ganancia del SiPM). Da un piso al FWHM
      relativo.
    - alpha [sqrt(keV)]    : estadistica de fotoelectrones, sigma ~ sqrt(N_pe).
    - beta  [keV]          : ruido electronico + linea base, constante en E.
    """
    E = np.maximum(np.asarray(E, float), 0.0)
    return np.sqrt(np.maximum(a_rel, 0.0) ** 2 * E ** 2
                   + np.maximum(alpha, 0.0) ** 2 * E
                   + np.maximum(beta, 0.0) ** 2)


def modelo_ext(canales, centros, h1, hm, p, sigma_min=1e-6):
    """
    S(ch) = A * sum_j [w*h1_j + (1-w)*hm_j] * G(ch; D*E_j, D*sigma_E(E_j))*dch
            + C * exp(-ch/lambda)

    p : array de 8 numeros en el orden de NOMBRES. Los que no se ajusten se
        pasan fijos (C=0 apaga el exceso; beta=a_rel=0 deja la resolucion del
        modelo viejo).
    """
    canales = np.asarray(canales, float)
    p = np.asarray(p, float)
    A, D, alpha, w, C, lam, beta, a_rel = p

    ancho_ch = float(np.mean(np.diff(canales))) if canales.size > 1 else 1.0
    pesos = np.clip(w, 0.0, 1.0) * h1 + (1.0 - np.clip(w, 0.0, 1.0)) * hm

    mu = D * centros
    sd = np.maximum(D * sigma_energia(centros, alpha, beta, a_rel), sigma_min)
    z = (canales[:, None] - mu[None, :]) / sd[None, :]
    G = np.exp(-0.5 * z * z) / (np.sqrt(2.0 * np.pi) * sd[None, :]) * ancho_ch

    y = A * (G @ pesos)
    if C > 0 and lam > 0:
        y = y + C * forma_exceso(canales, lam)
    return y


# La forma del exceso es una eleccion, no fisica medida: cambiarla y ver cuanto
# se mueve D es el sistematico de esa eleccion (ver la demo).
FORMA_EXCESO = "exp"


def forma_exceso(canales, lam):
    """Perfil normalizado a 1 en canal 0 (`exp`) o en canal 1000 (`pot`)."""
    ch = np.maximum(np.asarray(canales, float), 1.0)
    if FORMA_EXCESO == "pot":
        return (ch / 1000.0) ** (-abs(lam) / 100.0)
    return np.exp(-ch / lam)


# ----------------------------------------------------------------------
# 3) ajuste (chi2 de Pearson, iterado)
# ----------------------------------------------------------------------
def _empaquetar(p_ini, libres):
    """Indices libres -> vector reducido, y la funcion que lo re-expande."""
    idx = [NOMBRES.index(n) if isinstance(n, str) else int(n) for n in libres]
    base = np.array(p_ini, float)

    def expandir(q):
        p = base.copy()
        p[idx] = q
        return p
    return idx, base[idx].copy(), expandir


COTAS = {
    "A": (0.0, np.inf), "D": (1e-6, np.inf), "alpha": (0.0, 60.0),
    "w": (0.0, 1.0), "C": (0.0, np.inf), "lambda": (10.0, 5000.0),
    "beta": (0.0, 500.0), "a_rel": (0.0, 1.0),
}


def ajustar(canales, cuentas, centros, h1, hm, intervalo, p_ini, libres,
            var_fn=None, errores=None, n_iter=3, E_ref=963.4):
    """
    Ajusta los parametros de `libres` sobre `intervalo`.

    E_ref : borde de referencia de la fuente [keV], solo para reportar la
            resolucion (963.4 en Co-60, 477.3 en Cs-137). No entra en el ajuste.

    var_fn : callable(modelo) -> varianza esperada. Si se da, el chi2 es el de
             PEARSON (varianza del modelo) y hay que iterar, porque la varianza
             depende de los parametros. Sin el, se usan `errores` (Neyman, que
             sesga D hacia abajo: CLAUDE.md 7.31).
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]

    idx, q0, expandir = _empaquetar(p_ini, libres)
    inf = np.array([COTAS[NOMBRES[i]][0] for i in idx], float)
    sup = np.array([COTAS[NOMBRES[i]][1] for i in idx], float)
    q0 = np.clip(q0, inf + 1e-12, sup - 1e-12)

    err0 = (np.asarray(errores, float)[m] if errores is not None
            else np.ones_like(y))
    w_peso = 1.0 / err0

    sol = None
    for it in range(max(n_iter, 1) if var_fn is not None else 1):
        def residuo(q, w_peso=w_peso):
            return (modelo_ext(x, centros, h1, hm, expandir(q)) - y) * w_peso

        sol = least_squares(residuo, q0, bounds=(inf, sup))
        q0 = sol.x
        if var_fn is not None:
            mod = modelo_ext(x, centros, h1, hm, expandir(sol.x))
            w_peso = 1.0 / np.sqrt(np.asarray(var_fn(mod), float))

    p = expandir(sol.x)
    n_lib = max(len(y) - len(sol.x), 1)
    chi2 = 2.0 * float(sol.cost)

    err_p = np.full(len(sol.x), np.nan)
    cov = None
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
        err_p = np.sqrt(np.abs(np.diag(cov)))
    except Exception:
        pass

    errores_dict = {NOMBRES[i]: float(e) for i, e in zip(idx, err_p)}
    return {
        "p": p, "libres": [NOMBRES[i] for i in idx], "idx": idx,
        "q": sol.x, "cov": cov, "err": errores_dict,
        "chi2": chi2, "ndf": n_lib, "chi2_red": chi2 / n_lib,
        "intervalo": (lo, hi), "exito": bool(sol.success),
        "E_ref": float(E_ref),
        "sigma_E_ref": float(sigma_energia(E_ref, p[I_ALPHA], p[I_BETA],
                                           p[I_AREL])),
        "fwhm_rel_ref": float(2.355 * sigma_energia(E_ref, p[I_ALPHA], p[I_BETA],
                                                    p[I_AREL]) / E_ref),
        "peso_final": w_peso,
        **{n: float(p[k]) for k, n in enumerate(NOMBRES)},
    }


# ----------------------------------------------------------------------
# 4) posterior para el MCMC
# ----------------------------------------------------------------------
def log_posterior(canales, cuentas, centros, h1, hm, intervalo, p_base, libres,
                  var_fn=None, errores=None, cotas_extra=None):
    """
    ln P(theta | datos) = -chi2/2 (+ el termino ln sigma^2 si la varianza
    depende de theta), con priors planos dentro de COTAS.

    Devuelve (lp, expandir, n_bins).
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]
    var_fija = (np.asarray(errores, float)[m] ** 2 if errores is not None
                else np.ones_like(y))

    idx, _, expandir = _empaquetar(p_base, libres)
    cotas = dict(COTAS)
    if cotas_extra:
        cotas.update(cotas_extra)
    inf = np.array([cotas[NOMBRES[i]][0] for i in idx], float)
    sup = np.array([cotas[NOMBRES[i]][1] for i in idx], float)

    def lp(q):
        q = np.asarray(q, float)
        if np.any(q < inf) or np.any(q > sup):
            return -np.inf
        mod = modelo_ext(x, centros, h1, hm, expandir(q))
        if var_fn is None:
            return -0.5 * float(np.sum((mod - y) ** 2 / var_fija))
        var = np.asarray(var_fn(mod), float)
        return -0.5 * float(np.sum((mod - y) ** 2 / var + np.log(var)))

    return lp, expandir, int(m.sum())


# ======================================================================
#  DEMO: la escalera de modelos sobre los datos reales
# ======================================================================
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from datos_io import spectrum

    MUESTRA = rutas.MUESTRA
    FONDO = rutas.FONDO
    N_CAN, RANGO = 500, 2500.0
    E_REF = 963.4

    print("=" * 74)
    print("ETAPA 9 - MODELO EXTENDIDO: mezcla de scattering + exceso a bajo canal")
    print("=" * 74)

    r = spectrum.espectro_neto(MUESTRA, FONDO, num_canales=N_CAN,
                               rango_max=RANGO, recortar_negativos=False)
    can, neto, err = r["canales"], r["neto"], r["error"]
    print(f"  tiempo vivo: muestra {r['t_vivo_muestra']:.0f} s, "
          f"fondo {r['t_vivo_fondo']:.0f} s   (k_fondo = {r['k_f']:.3f})")

    centros, h1, hm, w_nat = cargar_componentes(rama="Lvis_keVee", bins=175)
    print(f"  la simulacion trae w = {w_nat:.3f} de scattering simple")

    def var_de(iv):
        m = (can >= iv[0]) & (can <= iv[1])
        return spectrum.varianza_predicha(r["fondo"][m], r["k_m"], r["k_f"])

    # semilla: ajuste de 3 parametros en la zona limpia (arriba del exceso)
    IV_LIMPIA = (900.0, 1900.0)
    p0 = np.array([1.0, 1.25, 3.5, w_nat, 0.0, 250.0, 0.0, 0.0])
    p0[I_A] = max(neto[(can > 1000) & (can < 1100)].mean(), 1.0) / \
        max(modelo_ext(can, centros, h1, hm, p0)[(can > 1000) & (can < 1100)].mean(), 1e-9)

    escalera = [
        ("M0  3 par, ventana limpia", IV_LIMPIA, ["A", "D", "alpha"]),
        ("M1  + mezcla w", IV_LIMPIA, ["A", "D", "alpha", "w"]),
        ("M2  + exceso, ventana ancha", (300.0, 1900.0),
         ["A", "D", "alpha", "w", "C", "lambda"]),
        ("M3  + resolucion 3 terminos", (300.0, 1900.0),
         ["A", "D", "alpha", "w", "C", "lambda", "a_rel"]),
    ]

    fits = {}
    p_ini = p0.copy()
    for nom, iv, libres in escalera:
        if "C" in libres and p_ini[I_C] == 0.0:
            p_ini[I_C] = 50.0
        f = ajustar(can, neto, centros, h1, hm, iv, p_ini, libres,
                    var_fn=var_de(iv), errores=err)
        fits[nom] = f
        p_ini = f["p"].copy()
        eD = f["err"].get("D", np.nan)
        print(f"\n  {nom}   ventana {iv[0]:.0f}-{iv[1]:.0f}")
        print(f"     D = {f['D']:.4f} +- {eD:.4f} mV*ns/keVee"
              f"   ({1000/f['D']:.1f} keVee por 1000 mV*ns)")
        print(f"     alpha = {f['alpha']:.3f}   a_rel = {f['a_rel']:.4f}"
              f"   -> FWHM({f['E_ref']:.0f}) = {100*f['fwhm_rel_ref']:.1f}%")
        print(f"     w = {f['w']:.3f} +- {f['err'].get('w', float('nan')):.3f}"
              f"   (simulacion: {w_nat:.3f})")
        if f["C"] > 0:
            print(f"     exceso: C = {f['C']:.1f}  lambda = {f['lambda']:.0f} mV*ns")
        print(f"     chi2/ndf = {f['chi2_red']:.3f}  ({f['ndf']} gdl)")

    # ---------------------- figura ----------------------
    mejor = fits[escalera[-1][0]]
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))
    z = (can >= 150) & (can <= 2000)
    ax[0].errorbar(can[z], neto[z], yerr=err[z], fmt=".", ms=2, color="0.6",
                   lw=0.5, label="datos (neto, tiempo vivo)")
    for nom in fits:
        f = fits[nom]
        ax[0].plot(can[z], modelo_ext(can[z], centros, h1, hm, f["p"]),
                   lw=1.6, label=f"{nom.split()[0]}: D={f['D']:.3f}")
    pe = mejor["p"].copy()
    ax[0].plot(can[z], pe[I_C] * np.exp(-can[z] / pe[I_LAM]), "k--", lw=1,
               label="exceso solo")
    ax[0].axvspan(*mejor["intervalo"], color="green", alpha=0.06)
    ax[0].set_xlabel("Area del pulso [mV*ns]")
    ax[0].set_ylabel("Cuentas netas")
    ax[0].set_title("Escalera de modelos")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    p_sin = mejor["p"].copy()
    p_sin[I_C] = 0.0
    mod = modelo_ext(can[z], centros, h1, hm, mejor["p"])
    ax[1].plot(can[z], (neto[z] - mod) / np.maximum(err[z], 1), ".", ms=3)
    ax[1].axhline(0, color="k", lw=1)
    ax[1].axhspan(-1, 1, color="green", alpha=0.1)
    ax[1].axvspan(*mejor["intervalo"], color="green", alpha=0.06)
    ax[1].set_xlabel("Area del pulso [mV*ns]")
    ax[1].set_ylabel("(dato - modelo) / sigma")
    ax[1].set_title(f"Residuos del mejor modelo  (chi2/ndf = {mejor['chi2_red']:.2f})")
    ax[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("modelo_ext.png"), dpi=130)
    print("\n  Figura guardada: figuras/modelo_ext.png")
