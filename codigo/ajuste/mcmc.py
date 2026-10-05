"""
mcmc.py - Etapa 8: el POSTERIOR completo por Markov Chain Monte Carlo.

Que problema resuelve
---------------------
Las Etapas 6 y 7 entregan D con una barra de error que sale de invertir J^T J
(la covarianza LINEALIZADA de `least_squares`). Eso supone que el chi2 es una
parabola perfecta alrededor del minimo. Nadie lo habia verificado, y justo en
este ajuste corr(D, alpha) = -0.70: los dos parametros estan tan ligados que
mirar solo la diagonal de la covarianza esconde la mitad de la historia.

El MCMC no supone nada de eso. Muestrea la distribucion

    P(theta | datos)  proporcional a  exp(-chi2(theta)/2)   (priors planos)

y de las muestras salen directo: la mediana, el intervalo del 68%, la
correlacion real entre parametros, y -- lo mas importante para la pregunta
"?encontramos el optimo?" -- si el posterior tiene UN solo maximo o varios.

Por que Metropolis-Hastings escrito a mano y no una libreria
------------------------------------------------------------
Mismo criterio que `sa.py`: el algoritmo es corto y entenderlo es parte del
trabajo. Ademas aca hay una ventaja concreta: la matriz de propuesta se arma con
la covarianza que YA devuelve `ajustar_forward_folding`, escalada por 2.4^2/d
(el escalado optimo de Roberts & Rosenthal para propuestas gaussianas). Con eso
la cadena mezcla desde el primer paso, y la comparacion entre el posterior y esa
misma elipse es exactamente la pregunta que queremos contestar.

Diagnosticos incluidos (sin ellos una cadena no es un resultado)
----------------------------------------------------------------
  R_hat (Gelman-Rubin) : compara la dispersion ENTRE cadenas con la dispersion
                         DENTRO de cada cadena. Si las cadenas partieron
                         dispersas y terminan describiendo la misma nube,
                         R_hat -> 1. Criterio usual: R_hat < 1.01.
  tau_int              : tiempo de autocorrelacion. Las muestras de una cadena
                         NO son independientes; n_eff = n/tau es cuantas lo son.
  aceptacion           : deberia caer en 0.2-0.5. Muy alta = pasos demasiado
                         cortos (explora lento); muy baja = pasos demasiado
                         largos (se queda quieta).

Corre solo:  python mcmc.py   (test con verdad conocida, sin datos reales)
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/ajuste/mcmc.py` pone en sys.path la carpeta DEL ARCHIVO, no
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


# ----------------------------------------------------------------------
# 1) el muestreador
# ----------------------------------------------------------------------
def metropolis(log_post, x0, cov_prop, n_pasos=12000, semilla=0, escala=None,
               callback=None):
    """
    Metropolis-Hastings con propuesta gaussiana fija.

    log_post : callable(x) -> log de la densidad (sin normalizar). -inf = fuera
               del prior.
    cov_prop : covarianza de la propuesta ANTES de escalar. Lo natural es la
               covarianza del ajuste por minimos cuadrados: dice en que
               direcciones el chi2 es plano y en cuales empinado, asi que la
               propuesta ya viaja a lo largo de la degeneracion D-alpha en vez
               de chocarse contra ella.
    escala   : multiplicador de cov_prop. Default 2.4^2/d, el optimo para un
               objetivo gaussiano en d dimensiones.

    Devuelve dict con la cadena (n_pasos x d), su log_post y la aceptacion.
    """
    rng = np.random.default_rng(semilla)
    x = np.asarray(x0, float).copy()
    d = x.size
    if escala is None:
        escala = 2.4 ** 2 / d
    L = np.linalg.cholesky(np.asarray(cov_prop, float) * escala)

    lp = float(log_post(x))
    if not np.isfinite(lp):
        raise ValueError("El punto inicial tiene log_post = -inf (fuera del prior).")

    cadena = np.empty((n_pasos, d))
    logps = np.empty(n_pasos)
    aceptados = 0
    aviso = max(n_pasos // 20, 1)

    for i in range(n_pasos):
        xp = x + L @ rng.standard_normal(d)
        lpp = float(log_post(xp))
        # log(u) < lpp - lp  es el mismo criterio de Metropolis de sa.py con
        # T = 1, solo que aca NO se enfria: queremos muestrear, no minimizar.
        if lpp - lp > np.log(rng.random()):
            x, lp = xp, lpp
            aceptados += 1
        cadena[i] = x
        logps[i] = lp
        if callback is not None and (i + 1) % aviso == 0:
            callback(i + 1, n_pasos)

    return {"cadena": cadena, "logp": logps,
            "aceptacion": aceptados / n_pasos, "escala": escala}


def correr_cadenas(log_post, x0, cov_prop, n_cadenas=4, n_pasos=12000,
                   dispersion=2.0, semilla=0, callback=None):
    """
    Varias cadenas desde arranques DISPERSOS (x0 + dispersion*sigma al azar).

    Dispersar los arranques no es un detalle: R_hat compara la varianza entre
    cadenas contra la varianza interna. Si todas arrancan en el mismo punto,
    R_hat da ~1 aunque ninguna haya explorado nada.
    """
    rng = np.random.default_rng(semilla)
    x0 = np.asarray(x0, float)
    sig = np.sqrt(np.abs(np.diag(np.asarray(cov_prop, float))))
    salidas = []
    for k in range(n_cadenas):
        xk = x0 if k == 0 else x0 + dispersion * sig * rng.standard_normal(x0.size)
        if not np.isfinite(log_post(xk)):
            xk = x0
        cb = None
        if callback is not None:
            def cb(i, n, k=k):
                callback(k, n_cadenas, i, n)
        salidas.append(metropolis(log_post, xk, cov_prop, n_pasos,
                                  semilla=semilla + 100 * (k + 1), callback=cb))
    return salidas


# ----------------------------------------------------------------------
# 2) diagnosticos de convergencia
# ----------------------------------------------------------------------
def gelman_rubin(cadenas):
    """
    R_hat por parametro. cadenas: array (m, n, d) ya sin quemado.

    W = varianza promedio DENTRO de cada cadena, B/n = varianza ENTRE las medias.
    Si las cadenas exploran la misma distribucion, las dos estiman lo mismo.
    """
    c = np.asarray(cadenas, float)
    m, n, _ = c.shape
    if m < 2 or n < 2:
        return np.full(c.shape[2], np.nan)
    medias = c.mean(axis=1)                       # (m, d)
    W = c.var(axis=1, ddof=1).mean(axis=0)        # (d,)
    B_sobre_n = medias.var(axis=0, ddof=1)        # (d,)
    var_post = (n - 1) / n * W + B_sobre_n
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(np.where(W > 0, var_post / W, np.nan))


def autocorr_integrada(x, c=5.0):
    """
    tau_int de una serie 1D, por FFT y con la ventana automatica de Sokal.

    tau = 1 + 2*sum_t rho_t, cortando la suma en el primer M con M >= c*tau
    (sumar toda la cola solo agrega ruido).
    """
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 10 or np.allclose(x, 0.0):
        return np.nan
    n2 = 1 << (2 * n - 1).bit_length()
    f = np.fft.fft(x, n2)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    if acf[0] <= 0:
        return np.nan
    acf /= acf[0]
    taus = 2.0 * np.cumsum(acf) - 1.0
    m = np.arange(n)
    ok = m >= c * taus
    idx = int(np.argmax(ok)) if ok.any() else n - 1
    return float(taus[idx])


def resumen(cadenas, quemado=None, nombres=None):
    """
    Junta las cadenas y devuelve el resultado listo para reportar.

    quemado : pasos a descartar por cadena. None = 20% (con una propuesta
              sembrada en el optimo alcanza y sobra).
    """
    lista = [np.asarray(r["cadena"], float) for r in cadenas]
    n = min(len(c) for c in lista)
    q = int(0.2 * n) if quemado is None else int(quemado)
    q = min(max(q, 0), n - 2)
    apiladas = np.stack([c[q:n] for c in lista])          # (m, n-q, d)
    plana = apiladas.reshape(-1, apiladas.shape[2])       # (m*(n-q), d)

    med = np.median(plana, axis=0)
    p16, p84 = np.percentile(plana, [16, 84], axis=0)
    R = gelman_rubin(apiladas)
    tau = np.array([np.nanmean([autocorr_integrada(c[q:n, j]) for c in lista])
                    for j in range(plana.shape[1])])
    n_total = apiladas.shape[0] * apiladas.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        n_eff = np.where(np.isfinite(tau) & (tau > 0), n_total / tau, np.nan)
    cov = np.cov(plana.T)
    sd = np.sqrt(np.abs(np.diag(cov)))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(sd, sd)

    return {
        "muestras": plana, "cadenas": apiladas, "quemado": q,
        "mediana": med, "p16": p16, "p84": p84,
        "err_menos": med - p16, "err_mas": p84 - med,
        "sigma": sd, "cov": cov, "corr": corr,
        "R_hat": R, "tau": tau, "n_eff": n_eff,
        "aceptacion": np.array([r["aceptacion"] for r in cadenas]),
        "nombres": list(nombres) if nombres else
                   ["p%d" % j for j in range(plana.shape[1])],
    }


# ----------------------------------------------------------------------
# 3) para dibujar: contornos del posterior y elipse linealizada
# ----------------------------------------------------------------------
def contornos_2d(x, y, niveles=(0.68, 0.95), bins=60, suavizado=1.0):
    """
    Densidad 2D y los niveles que encierran esas fracciones de probabilidad.

    (`corner` no esta instalado y no hace falta: son cuatro lineas con
    histogram2d + el truco de ordenar la densidad y acumular.)
    """
    H, bx, by = np.histogram2d(np.asarray(x, float), np.asarray(y, float), bins=bins)
    if suavizado:
        try:
            from scipy.ndimage import gaussian_filter
            H = gaussian_filter(H, suavizado)
        except Exception:
            pass
    H = H.T
    total = H.sum()
    cx = 0.5 * (bx[:-1] + bx[1:])
    cy = 0.5 * (by[:-1] + by[1:])
    if total <= 0:
        return H, cx, cy, []
    orden = np.sort(H.ravel())[::-1]
    acum = np.cumsum(orden) / total
    cortes = []
    for f in niveles:
        k = int(np.searchsorted(acum, f))
        cortes.append(float(orden[min(k, len(orden) - 1)]))
    return H, cx, cy, sorted(set(cortes))


def elipse_covarianza(mu, cov, n_sigma=1.0, n=200):
    """Elipse de la covarianza LINEALIZADA, para superponerla al posterior."""
    mu = np.asarray(mu, float)
    vals, vecs = np.linalg.eigh(np.asarray(cov, float))
    vals = np.maximum(vals, 0.0)
    t = np.linspace(0, 2 * np.pi, n)
    circ = np.stack([np.cos(t), np.sin(t)])
    pts = (vecs @ (np.sqrt(vals)[:, None] * circ)) * n_sigma
    return mu[0] + pts[0], mu[1] + pts[1]


# ----------------------------------------------------------------------
# 4) el posterior de la calibracion (forward folding con Geant4)
# ----------------------------------------------------------------------
def log_posterior_sim(canales, cuentas, errores, centros, pesos, intervalo,
                      cotas, ajustar_beta=False, var_fn=None):
    """
    ln P(A, D, alpha [, beta] | datos) = -chi2/2, priors planos dentro de cotas.

    Se muestrea A tambien (no se proyecta como en el SA): asi el intervalo de D
    incluye la incertidumbre de la normalizacion, que es parte honesta del error.

    cotas  : (min, max) por parametro, en el orden [A, D, alpha (, beta)].
    var_fn : callable(modelo) -> varianza esperada. Si se da, la varianza depende
             de los parametros (chi2 de Pearson, sin el sesgo de Neyman: ver
             `spectrum.varianza_predicha`) y entonces hay que incluir el termino
             ln(sigma^2) de la gaussiana, que con varianza fija seria constante y
             se podia ignorar. Sin var_fn se usan los `errores` fijos.

    Devuelve (log_post, n_bins) - n_bins para poder calcular chi2/ndf despues.
    """
    from ajuste.sim_template import modelo_simulado

    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]
    var_fija = (np.asarray(errores, float)[m] ** 2 if errores is not None
                else np.ones_like(y))
    inf, sup = np.asarray(cotas[0], float), np.asarray(cotas[1], float)

    def lp(p):
        p = np.asarray(p, float)
        if np.any(p < inf) or np.any(p > sup):
            return -np.inf
        beta = p[3] if ajustar_beta else 0.0
        mod = modelo_simulado(x, centros, pesos, p[0], p[1], p[2], beta)
        if var_fn is None:
            return -0.5 * float(np.sum((mod - y) ** 2 / var_fija))
        var = np.asarray(var_fn(mod), float)
        return -0.5 * float(np.sum((mod - y) ** 2 / var + np.log(var)))

    return lp, int(m.sum())


def posterior_calibracion(espectro, centros, pesos, intervalo, fit,
                          n_cadenas=4, n_pasos=12000, semilla=1, callback=None,
                          usar_pearson=True):
    """
    El posterior de {A, D, alpha} para ESTOS datos y ESTA plantilla.

    espectro : el dict de `spectrum.espectro_neto` (necesita fondo, k_m, k_f
               para la varianza esperada).
    fit      : el resultado de `sim_template.ajustar_pearson` / `_forward_folding`.
               De ahi salen el punto de partida y la covarianza de la propuesta.

    Devuelve el dict de `resumen` con dos extras: 'fit' (el ajuste que lo sembro)
    y 'intervalo'.
    """
    from datos_io.spectrum import varianza_predicha

    canales = np.asarray(espectro["canales"], float)
    m = (canales >= intervalo[0]) & (canales <= intervalo[1])
    var_fn = (varianza_predicha(np.asarray(espectro["fondo"], float)[m],
                                espectro.get("k_m", 1.0), espectro.get("k_f", 1.0))
              if usar_pearson else None)

    x0 = np.array([fit["A"], fit["D"], fit["alpha"]], float)
    cotas = ([0.0, 0.5 * fit["D"], 0.5],
             [10.0 * fit["A"], 2.0 * fit["D"], 30.0])
    lp, n_bins = log_posterior_sim(canales, espectro["neto"], espectro["error"],
                                   centros, pesos, intervalo, cotas, var_fn=var_fn)

    cov = fit.get("cov")
    if cov is None or not np.all(np.isfinite(cov)):
        cov = np.diag((0.02 * np.abs(x0) + 1e-12) ** 2)
    cadenas = correr_cadenas(lp, x0, cov[:3, :3], n_cadenas=n_cadenas,
                             n_pasos=n_pasos, semilla=semilla, callback=callback)
    res = resumen(cadenas, nombres=["A", "D", "alpha"])
    res["fit"] = fit
    res["intervalo"] = tuple(intervalo)
    res["n_bins"] = n_bins
    return res


# ======================================================================
#  DEMO / VALIDACION con verdad conocida (sintetico primero, como siempre)
# ======================================================================
if __name__ == "__main__":
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ajuste.sim_template import (cargar_simulacion, modelo_simulado,
                              ajustar_forward_folding, ajustar_pearson)
    from datos_io.spectrum import varianza_predicha

    # Amplitud elegida para dar ~26 cuentas/bin, la estadistica REAL de este
    # experimento en la ventana de ajuste. Con esa estadistica el sesgo de
    # Neyman es -0.5% en D, mas grande que el error estadistico: por eso el
    # posterior se muestrea con la varianza del MODELO (CLAUDE.md 7.31).
    A_V, D_V, AL_V = 3.0e-2, 4.20, 5.70          # la VERDAD
    INTERVALO = (1850.0, 7850.0)

    print("=" * 72)
    print("ETAPA 8 - MCMC: test con verdad conocida")
    print("=" * 72)

    centros, pesos = cargar_simulacion(rama="Lvis_keVee", bins=175)
    canales = np.arange(5.0, 8000.0, 10.0)
    ideal = modelo_simulado(canales, centros, pesos, A_V, D_V, AL_V)

    rng = np.random.default_rng(7)
    datos = rng.poisson(np.maximum(ideal, 0.0)).astype(float)
    errores = np.sqrt(np.maximum(datos, 1.0))

    print("  verdad       : A = %.4e   D = %.3f   alpha = %.3f" % (A_V, D_V, AL_V))

    m_iv = (canales >= INTERVALO[0]) & (canales <= INTERVALO[1])
    var_fn = varianza_predicha(np.zeros(int(m_iv.sum())), 1.0, 0.0)  # sin fondo

    fit_ney = ajustar_forward_folding(canales, datos, centros, pesos, INTERVALO,
                                      errores=errores)
    fit = ajustar_pearson(canales, datos, centros, pesos, INTERVALO, var_fn,
                          errores=errores)
    print("  Neyman  (sigma del DATO)  : D = %.4f +- %.4f   (%+.2f%%)"
          % (fit_ney["D"], fit_ney["err_D"], 100 * (fit_ney["D"] / D_V - 1)))
    print("  Pearson (sigma del MODELO): D = %.4f +- %.4f   (%+.2f%%)   "
          "alpha = %.3f +- %.3f   chi2/ndf = %.3f"
          % (fit["D"], fit["err_D"], 100 * (fit["D"] / D_V - 1), fit["alpha"],
             fit["err_alpha"], fit["chi2_red"]))

    x0 = np.array([fit["A"], fit["D"], fit["alpha"]])
    cotas = ([0.0, 0.5 * fit["D"], 0.5], [10 * fit["A"], 2.0 * fit["D"], 25.0])
    lp, n_bins = log_posterior_sim(canales, datos, errores, centros, pesos,
                                   INTERVALO, cotas, var_fn=var_fn)

    t0 = time.time()
    cadenas = correr_cadenas(lp, x0, fit["cov"], n_cadenas=4, n_pasos=6000,
                             dispersion=2.0, semilla=1)
    res = resumen(cadenas, nombres=["A", "D", "alpha"])
    dt = time.time() - t0

    print("\n  MCMC: 4 cadenas x 6000 pasos en %.1f s   (aceptacion %.2f)"
          % (dt, res["aceptacion"].mean()))
    verdad = [A_V, D_V, AL_V]
    for j, nom in enumerate(res["nombres"]):
        m, lo68, hi68 = res["mediana"][j], res["p16"][j], res["p84"][j]
        dentro = "SI" if lo68 <= verdad[j] <= hi68 else "NO"
        print("    %-6s = %.4g  -%.3g +%.3g   verdad %.4g  dentro del 68%%: %s"
              "   R_hat = %.4f  n_eff = %.0f"
              % (nom, m, m - lo68, hi68 - m, verdad[j], dentro,
                 res["R_hat"][j], res["n_eff"][j]))

    print("\n  corr(D, alpha) posterior     = %+.3f" % res["corr"][1, 2])
    print("  corr(D, alpha) linealizada   = %+.3f" % fit["corr_D_alpha"])
    print("  sigma_D posterior            = %.4f" % res["sigma"][1])
    print("  sigma_D linealizada          = %.4f" % fit["err_D"])

    # ------------------------------ figura ------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))

    for k in range(res["cadenas"].shape[0]):
        ax[0].plot(res["cadenas"][k, :, 1], lw=0.6, alpha=0.8)
    ax[0].axhline(D_V, color="k", ls="--", label="verdad = %.2f" % D_V)
    ax[0].set_xlabel("paso (post quemado)"); ax[0].set_ylabel("D")
    ax[0].set_title("Trazas de D  (R_hat = %.4f)" % res["R_hat"][1], fontsize=10)
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    D, AL = res["muestras"][:, 1], res["muestras"][:, 2]
    H, cx, cy, niveles = contornos_2d(D, AL)
    ax[1].plot(D[::20], AL[::20], ".", ms=1, color="0.7", alpha=0.5)
    if niveles:
        ax[1].contour(cx, cy, H, levels=niveles, colors="tab:blue",
                      linewidths=[1.0, 2.0][:len(niveles)])
    ex, ey = elipse_covarianza([fit["D"], fit["alpha"]],
                               fit["cov"][1:3, 1:3], n_sigma=1.0)
    ax[1].plot(ex, ey, "crimson", lw=2, ls="--", label="elipse linealizada 1 sigma")
    ax[1].plot(D_V, AL_V, "k*", ms=14, label="verdad")
    ax[1].set_xlabel("D [mV*ns/keV]"); ax[1].set_ylabel("alpha [sqrt(keV)]")
    ax[1].set_title("Posterior 68/95 vs covarianza  (corr = %+.2f)"
                    % res["corr"][1, 2], fontsize=10)
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    ax[2].hist(D, bins=60, color="tab:blue", alpha=0.75, density=True)
    ax[2].axvline(res["mediana"][1], color="tab:blue", lw=2, label="mediana")
    ax[2].axvspan(res["p16"][1], res["p84"][1], color="tab:blue", alpha=0.15,
                  label="68%")
    ax[2].axvline(D_V, color="k", ls="--", label="verdad")
    ax[2].set_xlabel("D [mV*ns/keV]"); ax[2].set_ylabel("densidad")
    ax[2].set_title("Marginal de D", fontsize=10)
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("mcmc_demo.png"), dpi=130)
    print("\n  Figura guardada: figuras/mcmc_demo.png")
