"""
sim_template.py - calibracion por FORWARD FOLDING con la simulacion Geant4.

Que cambia respecto del paper
-----------------------------
El paper usa como plantilla la lineshape ANALITICA de Klein-Nishina (Ec. 10):
un solo scattering Compton en un medio infinito, sin geometria. Eso ignora:

  - scattering Compton MULTIPLE (en esta simulacion: 25% de los eventos, con
    Edep llegando a 1332 keV, muy por encima del borde de 1118) -> rellena la
    zona sobre el borde,
  - escape de los gammas dispersados (depende del tamano del centellador),
  - el quenching de Birks a baja energia (rama Lvis_keVee, QF medio 0.93),
  - la eficiencia relativa de las dos lineas gamma.

Todo eso el ajuste analitico lo tiene que absorber inflando sigma. Por eso daba
42% FWHM (CLAUDE.md 7.13). La simulacion ya contiene esa fisica, asi que aca la
usamos DIRECTO como plantilla y sigma vuelve a significar solo resolucion.

El modelo
---------
    S_model(ch) = A * SUM_j  h_j * G(ch ; D*E_j , D*sigma_E(E_j))

  h_j        : histograma simulado de energia depositada / luz visible
  D          : mV*ns por keV  <- LA CALIBRACION
  sigma_E(E) : resolucion en energia. Por defecto sqrt(alpha^2*E + beta^2), o sea
               un termino de estadistica de fotoelectrones (~sqrt(E), Tabla I del
               paper: sigma_LY y sigma_gain) mas uno constante (ruido electronico,
               sigma_elec). Solo 3-4 parametros: A, D, alpha [, beta].

Es la misma idea del paper (comparar el espectro medido con una prediccion
teorica ensanchada), cambiando la prediccion analitica por la simulada.

Corre solo:  python sim_template.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/ajuste/sim_template.py` pone en sys.path la carpeta DEL ARCHIVO, no
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


# ----------------------------------------------------------------------
# 1) cargar la simulacion
# ----------------------------------------------------------------------
def cargar_simulacion(ruta=ARCHIVO_SIM, rama="Lvis_keVee", bins=700,
                      rango=(0.0, 1400.0), filtro=None, arbol="ComptonSpectrum"):
    """
    Devuelve (E_centros, pesos): histograma de la simulacion en keV.

    rama   : 'Lvis_keVee' (luz visible, con Birks -> lo que ve el SiPM) o
             'Edep_keV' (energia depositada cruda).
    filtro : callable(dict de arrays) -> mascara booleana. Ej. para quedarse
             solo con scattering simple:  lambda d: d['NCompton'] == 1
    bins   : el binning de la plantilla NO cambia el resultado, porque el kernel
             gaussiano es muchisimo mas ancho que un bin: sigma_E(963) = 178 keV
             contra 2 keV de bin. Medido (CLAUDE.md 7.26):

               bins  ancho   D            alpha   chi2/ndf  eval
               700   2 keV   4.2293       5.733   0.963     6.8 ms
               350   4 keV   4.2293       5.733   0.963     4.0 ms
               175   8 keV   4.2294       5.733   0.963     1.1 ms
                70  20 keV   4.2294       5.731   0.963     0.3 ms

             Se deja 700 de default por compatibilidad, pero para el SA y el
             MCMC (miles de evaluaciones) hay que usar 175: mismo numero, 6x
             mas rapido.
    """
    import uproot
    with uproot.open(ruta) as f:
        t = f[arbol]
        ramas = [rama]
        for extra in ("NCompton", "NPhoto", "FullAbs", "E0_keV"):
            if extra in t.keys() and extra not in ramas:
                ramas.append(extra)
        d = t.arrays(ramas, library="np")

    x = d[rama]
    if filtro is not None:
        x = x[filtro(d)]
    pesos, bordes = np.histogram(x, bins=bins, range=rango)
    centros = 0.5 * (bordes[:-1] + bordes[1:])
    return centros, pesos.astype(float)


def lineas_primarias(ruta=ARCHIVO_SIM, arbol="ComptonSpectrum", tol=2.0,
                     min_frac=0.02):
    """
    Que gammas se generaron en la simulacion, leyendo la rama `E0_keV`.

    Devuelve [(E_keV, fraccion), ...] ordenado por energia, o None si el arbol
    no tiene `E0_keV` (simulaciones viejas o de otro generador).

    Para que: la plantilla de forward folding ES la fuente. Ajustar un espectro
    de Cs-137 contra una simulacion de Co-60 no da error en ningun lado -- da un
    D silenciosamente equivocado, porque los bordes de la plantilla estan en
    963 y 1118 keV y los datos tienen uno solo en 477. Este chequeo es barato
    (una rama, sin cortes) y evita justo ese modo de falla.

    tol      : keV con que se agrupan energias vecinas (la simulacion guarda el
               valor exacto, pero conviene tolerar redondeos).
    min_frac : fraccion minima de eventos para contar como linea (descarta
               colas numericas).
    """
    import uproot
    with uproot.open(ruta) as f:
        t = f[arbol]
        if "E0_keV" not in t.keys():
            return None
        E0 = t["E0_keV"].array(library="np")
    E0 = np.asarray(E0, float)
    E0 = E0[np.isfinite(E0) & (E0 > 0)]
    if E0.size == 0:
        return None

    orden = np.sort(E0)
    lineas = []
    inicio = 0
    for i in range(1, orden.size + 1):
        if i == orden.size or orden[i] - orden[inicio] > tol:
            grupo = orden[inicio:i]
            lineas.append((float(np.mean(grupo)), grupo.size / E0.size))
            inicio = i
    return [(E, f) for E, f in lineas if f >= min_frac]


def comparar_con_fuente(lineas_sim, lineas_fuente, tol=5.0):
    """
    (compatible, mensaje) entre lo que simula el .root y la fuente elegida.

    `lineas_sim` es lo que devuelve `lineas_primarias` (o None si no se pudo
    leer: en ese caso no se puede afirmar nada y se devuelve compatible=True
    con un aviso, para no bloquear simulaciones sin la rama E0_keV).
    """
    if lineas_sim is None:
        return True, ("la simulacion no tiene la rama E0_keV: no se puede "
                      "verificar que corresponda a la fuente elegida.")
    Es = np.array([E for E, _ in lineas_sim], float)
    Ef = np.array(sorted(lineas_fuente), float)
    txt_sim = ", ".join(f"{E:.1f}" for E in Es)
    txt_fu = ", ".join(f"{E:.1f}" for E in Ef)
    if Es.size == Ef.size and np.all(np.abs(Es - Ef) <= tol):
        return True, f"la simulacion genera {txt_sim} keV: coincide con la fuente."
    return False, (f"la simulacion genera {txt_sim} keV y la fuente elegida "
                   f"tiene {txt_fu} keV. La plantilla NO corresponde a estos "
                   f"datos: el forward folding daria un D equivocado sin avisar.")


# ----------------------------------------------------------------------
# 2) el forward model: simulacion -> canal, con ensanchamiento
# ----------------------------------------------------------------------
def sigma_energia(E, alpha, beta=0.0):
    """sigma_E(E) = sqrt(alpha^2 * E + beta^2)  [keV].

    El termino alpha^2*E es estadistica de fotoelectrones (sigma ~ sqrt(N_pe));
    beta es el piso constante (ruido electronico). Tabla I del paper.
    """
    return np.sqrt(np.maximum(alpha, 0.0) ** 2 * np.maximum(E, 0.0) + beta ** 2)


def modelo_simulado(canales, centros, pesos, A, D, alpha, beta=0.0, sigma_min=1e-6):
    """
    S_model(ch) = A * sum_j pesos_j * G(ch; D*E_j, D*sigma_E(E_j)) * ancho_bin_ch

    Se arma la matriz de respuesta (n_canales x n_energias) de una vez. Es una
    convolucion de ancho VARIABLE (sigma crece con E), asi que no se puede hacer
    con np.convolve: hay que sumar gaussianas.
    """
    canales = np.asarray(canales, float)
    ancho_ch = float(np.mean(np.diff(canales)))

    mu = D * centros                                   # centro de cada gaussiana
    sd = np.maximum(D * sigma_energia(centros, alpha, beta), sigma_min)

    # matriz de respuesta: filas = canales medidos, columnas = energias simuladas
    z = (canales[:, None] - mu[None, :]) / sd[None, :]
    G = np.exp(-0.5 * z * z) / (np.sqrt(2.0 * np.pi) * sd[None, :]) * ancho_ch
    return A * (G @ pesos)


# ----------------------------------------------------------------------
# 2-bis) costo con la amplitud PROYECTADA (VARPRO), para el SA
# ----------------------------------------------------------------------
def amplitud_optima(base, y, w):
    """
    A que minimiza sum(((A*base - y)*w)^2). Es lineal: se resuelve exacto, no
    se busca. `base` es el modelo con A=1.
    """
    bw = base * w
    den = float(bw @ bw)
    return max(float(bw @ (y * w)) / den, 0.0) if den > 0 else 0.0


def costo_sim_varpro(canales, cuentas, centros, pesos, intervalo, errores=None,
                     ajustar_beta=False):
    """
    Devuelve (E, amplitud): el costo chi2 sobre p = [D, alpha] (+beta) con A
    resuelta EXACTO en cada evaluacion, y una funcion para recuperar esa A.

    Por que sacar A del vector de busqueda: el paso del SA es 0.1*delta_i por
    coordenada, asi que mezclar una amplitud (~1e-3 aca) con D (~4) en el mismo
    vector lo rompe (CLAUDE.md 7.21). A es lineal -> VARPRO, igual que
    `sa.costo_dos_bordes_varpro` hace con A1 y A2.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]
    w = 1.0 / np.asarray(errores, float)[m] if errores is not None else np.ones_like(y)

    def base_de(p):
        beta = p[2] if (ajustar_beta and len(p) > 2) else 0.0
        return modelo_simulado(x, centros, pesos, 1.0, max(p[0], 1e-9),
                               max(p[1], 0.0), beta)

    def E(p):
        base = base_de(p)
        r = (amplitud_optima(base, y, w) * base - y) * w
        return float(r @ r)

    def amplitud(p):
        return amplitud_optima(base_de(p), y, w)

    return E, amplitud


# ----------------------------------------------------------------------
# 3) ajuste
# ----------------------------------------------------------------------
def ajustar_pearson(canales, cuentas, centros, pesos, intervalo, var_fn,
                    errores=None, n_iter=3, **kw):
    """
    Ajuste con la varianza calculada del MODELO en vez de la de los DATOS.

    Es el chi2 de Pearson: mata el sesgo de Neyman (`spectrum.varianza_predicha`)
    a costa de tener que iterar, porque la varianza depende de los parametros.
    Tres iteraciones alcanzan de sobra: el peso cambia en el segundo decimal.

    var_fn : callable(modelo_en_la_ventana) -> varianza esperada, tipicamente
             `spectrum.varianza_predicha(fondo, k_m, k_f)` restringida a la ventana.
    """
    canales = np.asarray(canales, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)

    fit = ajustar_forward_folding(canales, cuentas, centros, pesos, intervalo,
                                  errores=errores, **kw)
    for _ in range(n_iter):
        modelo = modelo_simulado(canales[m], centros, pesos, fit["A"], fit["D"],
                                 fit["alpha"], fit["beta"])
        err = np.array(errores, float) if errores is not None \
            else np.ones_like(canales)
        err = err.copy()
        err[m] = np.sqrt(var_fn(modelo))
        kw2 = dict(kw)
        kw2.pop("p0", None)
        fit = ajustar_forward_folding(canales, cuentas, centros, pesos, intervalo,
                                      errores=err, p0=[fit["A"], fit["D"],
                                                       fit["alpha"]], **kw2)
    fit["errores_pearson"] = err[m]
    return fit


def ajustar_forward_folding(canales, cuentas, centros, pesos, intervalo,
                            p0=None, ajustar_beta=False, D_rango=None,
                            errores=None, E_ref=963.4):
    """
    Ajusta {A, D, alpha [, beta]} por minimos cuadrados sobre el intervalo.
    D es la calibracion: mV*ns por keV.

    E_ref : borde de referencia de la fuente [keV], SOLO para reportar la
            resolucion (`sigma_E_ref`, `fwhm_rel_ref`). No entra en el ajuste.
            963.4 es el borde de la primera linea del Co-60; con Cs-137 hay que
            pasar 477.3, porque un FWHM sin decir a que energia no significa
            nada (sigma va como sqrt(E)).

    errores : barra de error por bin (spectrum.espectro_neto la devuelve como
              'error'). Si se da, el ajuste es PONDERADO (chi2 de verdad) y se
              devuelven chi2/ndf y las incertidumbres de los parametros. Sin
              pesos, los bins del plateau (muchas cuentas) dominan el ajuste y
              no hay forma de poner barras de error al resultado.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    lo, hi = intervalo
    m = (canales >= lo) & (canales <= hi)
    x, y = canales[m], cuentas[m]
    w = 1.0 / np.asarray(errores, float)[m] if errores is not None else np.ones_like(y)

    if p0 is None:
        # semilla por barrido grueso de D (A sale exacto por escala lineal)
        if D_rango is None:
            D_rango = (0.5 * hi / centros.max(), 1.5 * hi / centros.max())
        mejor = (np.inf, None)
        for D in np.linspace(D_rango[0], D_rango[1], 40):
            for al in np.linspace(0.5, 12.0, 24):
                base = modelo_simulado(x, centros, pesos, 1.0, D, al)
                A = amplitud_optima(base, y, w)           # A optimo exacto
                res = (A * base - y) * w
                sse = float(res @ res)
                if sse < mejor[0]:
                    mejor = (sse, [A, D, al])
        p0 = mejor[1] + ([1.0] if ajustar_beta else [])
    elif ajustar_beta and len(p0) == 3:
        p0 = list(p0) + [1.0]

    if ajustar_beta:
        cotas = ([0.0, 1e-6, 0.0, 0.0], [np.inf, np.inf, 100.0, 500.0])

        def residuo(p):
            return (modelo_simulado(x, centros, pesos, p[0], p[1], p[2], p[3]) - y) * w
    else:
        cotas = ([0.0, 1e-6, 0.0], [np.inf, np.inf, 100.0])

        def residuo(p):
            return (modelo_simulado(x, centros, pesos, p[0], p[1], p[2]) - y) * w

    p0 = [min(max(v, a), b) for v, a, b in zip(p0, *cotas)]
    sol = least_squares(residuo, p0, bounds=cotas)

    A, D, alpha = sol.x[0], sol.x[1], sol.x[2]
    beta = sol.x[3] if ajustar_beta else 0.0
    n_lib = max(len(y) - len(sol.x), 1)
    chi2 = 2.0 * float(sol.cost)

    # incertidumbres desde el jacobiano: cov = (J^T J)^-1 (residuos ya pesados)
    err_p = np.full(len(sol.x), np.nan)
    cov = None
    corr_D_alpha = np.nan
    try:
        JTJ = sol.jac.T @ sol.jac
        cov = np.linalg.inv(JTJ)
        if errores is None:                       # sin pesos: escalar por el rms
            cov = cov * chi2 / n_lib
        err_p = np.sqrt(np.abs(np.diag(cov)))
        corr_D_alpha = float(cov[1, 2] / np.sqrt(cov[1, 1] * cov[2, 2]))
    except Exception:
        pass

    return {
        "A": A, "D": D, "alpha": alpha, "beta": beta,
        "err_D": float(err_p[1]), "err_alpha": float(err_p[2]),
        "err_p": err_p, "cov": cov, "x": sol.x,
        "corr_D_alpha": corr_D_alpha,
        "E_ref": float(E_ref),
        "sigma_E_ref": float(sigma_energia(E_ref, alpha, beta)),
        "fwhm_rel_ref": float(2.355 * sigma_energia(E_ref, alpha, beta) / E_ref),
        "coste": float(sol.cost), "chi2": chi2, "ndf": n_lib,
        "chi2_red": chi2 / n_lib,
        "rms": float(np.sqrt(2 * sol.cost / n_lib)),
        "exito": sol.success, "intervalo": (lo, hi), "p0": p0,
    }


# ----------------------------------------------------------------------
# 4) ajuste GLOBAL: recocido simulado + pulido local
# ----------------------------------------------------------------------
def ajustar_sa(canales, cuentas, centros, pesos, intervalo, errores=None,
               D0=None, alpha0=6.0, cotas_D=None, cotas_alpha=(0.5, 20.0),
               semillas=(3, 11, 42), E_ref=963.4, **kw_sa):
    """
    Igual que `ajustar_forward_folding` pero encontrando el optimo con el SA de
    la Etapa 5 (Seccion II-A del paper) en vez de con un barrido de grilla, y
    puliendo despues con `least_squares`.

    Que aporta: `ajustar_forward_folding` necesita una grilla de 40x24 para
    sembrarse. El SA llega al mismo punto SIN esa grilla, partiendo de una sola
    semilla (la del NCC). O sea: verifica que el optimo es global y no un
    artefacto de por donde arrancamos.

    Devuelve dict con 'sa' (el crudo del recocido), 'fit' (el pulido, con
    covarianza y chi2/ndf) y 'hist' para graficar la traza.
    """
    from etapas.sa import sa_multiarranque

    E, amplitud = costo_sim_varpro(canales, cuentas, centros, pesos, intervalo,
                                   errores=errores)
    if D0 is None:
        D0 = intervalo[1] / (1.6 * float(np.max(centros)))
    if cotas_D is None:
        cotas_D = (0.4 * D0, 2.2 * D0)
    cotas = ([cotas_D[0], cotas_alpha[0]], [cotas_D[1], cotas_alpha[1]])

    r = sa_multiarranque(E, [float(D0), float(alpha0)], cotas,
                         semillas=semillas, **kw_sa)
    D_sa, alpha_sa = float(r["x"][0]), float(r["x"][1])
    A_sa = float(amplitud(r["x"]))

    fit = ajustar_forward_folding(canales, cuentas, centros, pesos, intervalo,
                                  p0=[A_sa, D_sa, alpha_sa], errores=errores,
                                  E_ref=E_ref)
    return {"sa": r, "hist": r["hist"], "fit": fit,
            "A_sa": A_sa, "D_sa": D_sa, "alpha_sa": alpha_sa,
            "costo_sa": float(r["costo"]), "n_evals": int(r["n_evals"])}


# ======================================================================
#  DEMO: calibracion de los datos reales con la simulacion
# ======================================================================
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from datos_io import spectrum

    ARCHIVO_MUESTRA = rutas.MUESTRA
    ARCHIVO_FONDO = rutas.FONDO
    INTERVALO = (2000.0, 8000.0)      # sobre el umbral de trigger

    print("=" * 72)
    print("CALIBRACION POR FORWARD FOLDING CON LA SIMULACION")
    print("=" * 72)

    r = spectrum.espectro_neto(ARCHIVO_MUESTRA, ARCHIVO_FONDO, umbral_ruido=0.05,
                               num_canales=800, rango_max=8000.0,
                               recortar_negativos=False)
    can, neto = r["canales"], r["neto"]

    resultados = {}
    for etiqueta, rama, filtro in [
        ("Lvis (luz, con Birks)", "Lvis_keVee", None),
        ("Edep (energia cruda)", "Edep_keV", None),
        ("Edep solo 1 scattering", "Edep_keV", lambda d: d["NCompton"] == 1),
    ]:
        centros, pesos = cargar_simulacion(rama=rama, filtro=filtro)
        fit = ajustar_forward_folding(can, neto, centros, pesos, INTERVALO,
                                      errores=r["error"])
        fit["centros"], fit["pesos"] = centros, pesos
        resultados[etiqueta] = fit
        print(f"\n  {etiqueta}")
        print(f"    D     = {fit['D']:.4f} +- {fit['err_D']:.4f} mV*ns/keV   "
              f"( = {1000/fit['D']:.1f} keV por 1000 mV*ns )")
        print(f"    alpha = {fit['alpha']:.3f} +- {fit['err_alpha']:.3f} sqrt(keV)"
              f"   beta = {fit['beta']:.1f} keV")
        print(f"    sigma({fit['E_ref']:.0f} keV) = {fit['sigma_E_ref']:.1f} keV"
              f"  ->  {100*fit['fwhm_rel_ref']:.1f}% FWHM")
        print(f"    chi2/ndf = {fit['chi2_red']:.2f}  ({fit['ndf']} gdl)"
              f"   corr(D,alpha) = {fit['corr_D_alpha']:+.3f}")

    print("\n" + "=" * 72)
    print("  COMPARACION CON EL METODO ANALITICO (Klein-Nishina, Etapa 6)")
    print("=" * 72)
    D_kn = 4.344       # ligar_origen=True, ver CLAUDE.md 8
    principal = resultados["Lvis (luz, con Birks)"]
    print(f"  D analitico (Klein-Nishina, dos bordes) = {D_kn:.3f} mV*ns/keV")
    print(f"  D con simulacion (Lvis)                 = {principal['D']:.3f} mV*ns/keV")
    print(f"  diferencia                              = "
          f"{100*(principal['D']/D_kn - 1):+.1f}%")
    print(f"\n  FWHM(963 keV): {100*principal['fwhm_rel_ref']:.1f}% con simulacion, "
          f"vs 42% que pedia el modelo analitico.")

    # ------------------------------ figura ------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))

    z = (can >= 1200) & (can <= 7500)
    ax[0].step(can[z], neto[z], where="mid", color="0.7", lw=1, label="datos (neto)")
    colores = {"Lvis (luz, con Birks)": "b", "Edep (energia cruda)": "g",
               "Edep solo 1 scattering": "r"}
    for etiqueta, fit in resultados.items():
        y = modelo_simulado(can[z], fit["centros"], fit["pesos"], fit["A"], fit["D"],
                            fit["alpha"], fit["beta"])
        ax[0].plot(can[z], y, colores[etiqueta], lw=2,
                   label=f"{etiqueta}: D={fit['D']:.3f}")
    ax[0].axvspan(INTERVALO[0], INTERVALO[1], color="green", alpha=0.06)
    ax[0].set_title("Forward folding: simulacion ensanchada vs datos")
    ax[0].set_xlabel("Area del pulso [mV*ns]"); ax[0].set_ylabel("Cuentas norm.")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    # las plantillas crudas, en energia: se ve el efecto del multiple scattering
    for etiqueta, fit in resultados.items():
        p = fit["pesos"] / fit["pesos"].max()
        ax[1].plot(fit["centros"], p, colores[etiqueta], lw=1.5, label=etiqueta)
    for E, txt in [(963.4, "borde 963"), (1118.1, "borde 1118")]:
        ax[1].axvline(E, color="k", ls="--", alpha=0.6)
        ax[1].text(E, 1.02, txt, rotation=90, fontsize=7, ha="right")
    ax[1].set_yscale("log"); ax[1].set_ylim(1e-4, 2)
    ax[1].set_title("Plantillas simuladas (sin ensanchar)")
    ax[1].set_xlabel("Energia [keV]"); ax[1].set_ylabel("normalizado")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("calibracion_simulacion.png"), dpi=130)
    print("\n  Figura guardada: figuras/calibracion_simulacion.png")
