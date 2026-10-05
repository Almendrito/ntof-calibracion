"""
metodo.py — la cadena de la Fig. 2 del paper, de punta a punta.

    Start
      -> Energy Spectrum (Channel, Counts)
      -> Implement NCC Automatic Interval Selection
      -> Calculate Pearson Correlation Coefficient
      -> Maximum Data Interval of rho?
      -> Set Theta = {A1, A2, B1, sigma1, sigma2, D}
      -> Use SA to Solve Local Optimal Problems
      -> Optimal Parameter Solution?
      -> Fit Energy-Band Structure Relationship
    Stop

Dos casos, los dos del paper: una sola linea gamma (Ec. 8-10, Theta = {A,B,sigma})
y dos bordes solapados (Ec. 11-13, Theta = {A1,A2,B1,sigma1,sigma2,D}).

Decisiones que el paper deja abiertas, tomadas aca y marcadas
------------------------------------------------------------
- La plantilla del NCC para el caso de DOS bordes se construye con la linea
  gamma de MENOR energia. El paper describe la plantilla en singular ("the
  predicted lineshape") y su Fig. 2 pone el NCC ANTES de fijar Theta, asi que en
  ese punto todavia no hay D con que armar la composicion de la Ec. 13.
- Las cotas del SA ("d_i defining the parameter search range") no estan en el
  paper. Se derivan de la semilla del NCC y estan en `COTAS_REL`.
- El "Optimal Parameter Solution? No -> modify" del diagrama se implementa como
  un pulido local (least_squares) desde el optimo del SA. Sin ese lazo el
  resultado depende del ultimo salto aceptado. Se puede apagar con pulir=False.
"""

import os
import sys

import numpy as np
from scipy.optimize import least_squares, nnls

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import ecuaciones as ec
from replica_paper import ncc as ncc_mod
from replica_paper import sa as sa_mod

# Grilla de sigma con que se prueba la plantilla S2 del NCC. El paper no dice
# con que sigma se construye S2 (ver ncc.py); gana el de mayor rho, que es el
# unico criterio que el paper da.            [la grilla NO esta en el paper]
SIGMAS_S2 = (0.03, 0.05, 0.07, 0.09, 0.12, 0.16, 0.20, 0.25)

# Cotas del SA relativas a la semilla del NCC.   [NO estan en el paper]
COTAS_REL = dict(A=(0.0, 8.0), B=(0.6, 1.4), sigma=(0.15, 3.0), D=(0.6, 1.4))


def _ventana(canales, cuentas, intervalo):
    ch = np.asarray(canales, float)
    m = (ch >= intervalo[0]) & (ch <= intervalo[1])
    return ch[m], np.asarray(cuentas, float)[m], m


def _cerrar(modelo_fn, xw, yw, r_sa, cotas, pulir):
    """
    El lazo "Optimal Parameter Solution? No -> modify" de la Fig. 2.

    Pule con least_squares CADA salida del SA (no solo la mejor) y se queda con
    la de menor loss. Devuelve tambien que fraccion de los arranques llega al
    mismo optimo, que es la medida de si el SA del paper resuelve el problema o
    depende de la suerte.
    """
    lo = [c[0] for c in cotas]
    hi = [c[1] for c in cotas]
    salidas = r_sa.get("todas", [r_sa])
    finales = []
    for s in salidas:
        p = np.clip(np.asarray(s["x"], float), lo, hi)
        if pulir:
            try:
                # x_scale="jac" no es una eleccion de modelo: A vale ~1e30
                # (absorbe el prefactor de la Ec. 10) y D ~1, o sea 30 ordenes
                # de magnitud de diferencia. Sin reescalar, la region de
                # confianza no puede mover los parametros chicos.
                sol = least_squares(lambda q: modelo_fn(xw, q) - yw, p,
                                    bounds=(lo, hi), x_scale="jac",
                                    max_nfev=40000)
                p = sol.x
            except Exception:
                pass
        r = modelo_fn(xw, p) - yw
        finales.append((float(np.sum(r * r)), p))
    finales.sort(key=lambda t: t[0])
    mejor_loss, mejor_p = finales[0]
    iguales = sum(1 for L, _ in finales if L <= mejor_loss * 1.001)
    return mejor_p, dict(n_arranques=len(finales),
                         n_en_el_optimo=int(iguales),
                         losses=[L for L, _ in finales])


def semilla_por_barrido(xw, yw, eg1, eg2, cotas, delta_e="bordes",
                        n_B=20, n_D=20, n_s=3, k=6, ligar_origen=False):
    """
    Barrido grueso de (B1, D, sigma1) con A1 y A2 exactas por nnls, para SEMBRAR
    el SA.

    Por que hace falta, medido y no supuesto: el paper no dice de donde sale x0.
    El NCC entrega B y sigma pero NO D, que en Theta es libre. Con arranques
    dispersos al azar (24 probados sobre un espectro sintetico de verdad
    conocida) ni el SA ni el pulido local llegan al optimo global: se quedan en
    minimos con A2/A1 absurdo. Sembrando con este barrido si llegan. O sea: el
    SA con los 6 parametros de la Ec. 11 NO es un optimizador global en la
    practica, necesita ayuda. [El barrido no esta en el paper.]

    Devuelve las `k` mejores celdas como vectores Theta.
    """
    E1 = float(ec.e_max(eg1))
    E2 = float(ec.e_max(eg2))
    Bs = np.linspace(cotas[2][0], cotas[2][1], int(n_B))
    Ds = np.linspace(cotas[5][0], cotas[5][1], int(n_D))
    if ligar_origen:
        Ds = None
    ss = np.exp(np.linspace(np.log(cotas[3][0]), np.log(cotas[3][1]), int(n_s) + 2))[1:-1]
    obj = np.maximum(yw, 0.0)
    cache1 = {}
    salidas = []
    for B1 in Bs:
        for s1 in ss:
            key = (float(B1), float(s1))
            if key not in cache1:
                cache1[key] = ec.modelo(xw, eg1, 1.0, B1, s1)
            c1 = cache1[key]
            s2 = s1 * E2 / E1
            for D in ([B1 / E1] if Ds is None else Ds):
                B2 = B1 + D * ec._delta_e(eg1, eg2, delta_e)
                c2 = ec.modelo(xw, eg2, 1.0, B2, s2)
                M = np.column_stack([c1, c2])
                a, _ = nnls(M, obj)
                r = M @ a - yw
                salidas.append((float(np.sum(r * r)),
                                np.array([max(a[0], 1e-8), max(a[1], 1e-8),
                                          float(B1), float(s1), float(s2),
                                          float(D)])))
    salidas.sort(key=lambda t: t[0])
    return [p for _, p in salidas[:int(k)]]


def _amplitud_exacta(base, y):
    """
    A optima por minimos cuadrados con el resto de los parametros fijos.

    No es un agregado al modelo: A entra LINEAL en la Ec. 8, asi que su optimo
    tiene forma cerrada. Se usa solo para SEMBRAR el SA, que despues la ajusta
    como un parametro mas de Theta, como dice el paper.
    """
    d = float(np.dot(base, base))
    return float(np.dot(base, y) / d) if d > 0 else 1.0


# --------------------------------------------------------------- un borde ---
def ajustar_un_borde(canales, cuentas, e_gamma, sigma_rel=None,
                     intervalo=None, objetivo="ec14", pulir=True,
                     semillas_sa=3, ncc_kw=None, sa_kw=None):
    """
    Ec. 8-10 con Theta = {A, B, sigma}.

    `sigma_rel` es el sigma con que se construye la plantilla S2 del NCC. Si es
    None se barre una grilla y gana el de mayor rho, que es el unico criterio
    que el paper da para elegir entre ventanas.
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    ncc_kw = dict(ncc_kw or {})
    sa_kw = dict(sa_kw or {})

    # --- 1. NCC: intervalo y semilla -------------------------------------
    if sigma_rel is None:
        cand = ncc_mod.barrer_sigma_rel(ch, y, e_gamma,
                                        SIGMAS_S2,
                                        **ncc_kw)
        r_ncc = max(cand, key=lambda r: r["rho"])
        barrido = [(c["sigma_rel"], c["rho"], c["B"]) for c in cand]
    else:
        r_ncc = ncc_mod.buscar_intervalo(ch, y, e_gamma, sigma_rel=sigma_rel,
                                         **ncc_kw)
        barrido = [(r_ncc["sigma_rel"], r_ncc["rho"], r_ncc["B"])]

    if intervalo is None:
        intervalo = (r_ncc["lo"], r_ncc["hi"])
    xw, yw, _ = _ventana(ch, y, intervalo)

    B0, s0 = r_ncc["B"], r_ncc["sigma"]
    A0 = _amplitud_exacta(ec.modelo(xw, e_gamma, 1.0, B0, s0), yw)

    # --- 2. SA sobre Theta ------------------------------------------------
    def modelo_fn(x, p):
        return ec.modelo(x, e_gamma, p[0], p[1], p[2])

    cotas = [(COTAS_REL["A"][0] * A0, COTAS_REL["A"][1] * A0),
             (COTAS_REL["B"][0] * B0, COTAS_REL["B"][1] * B0),
             (COTAS_REL["sigma"][0] * s0, COTAS_REL["sigma"][1] * s0)]
    costo = (sa_mod.costo_ec14(modelo_fn) if objetivo == "ec14"
             else sa_mod.costo_ncc(modelo_fn))
    x0s = [np.array([A0, B0, s0])]
    rng = np.random.default_rng(0)
    for _ in range(max(int(semillas_sa) - 1, 0)):
        x0s.append(np.array([A0 * rng.uniform(0.5, 2.0),
                             B0 * rng.uniform(0.85, 1.15),
                             s0 * rng.uniform(0.5, 2.0)]))
    r_sa = sa_mod.multiarranque(costo, x0s, cotas, xw, yw, **sa_kw)

    # --- 3. pulido local (el lazo "Optimal Parameter Solution?") ----------
    p, dispersion = _cerrar(modelo_fn, xw, yw, r_sa, cotas, pulir)

    modelo_y = modelo_fn(xw, p)
    res = yw - modelo_y
    gdl = max(xw.size - p.size, 1)
    var = np.maximum(np.abs(modelo_y), 1.0)      # Poisson del MODELO
    return dict(tipo="un_borde", e_gamma=float(e_gamma),
                A=float(p[0]), B=float(p[1]), sigma=float(p[2]),
                D=float(p[1] / ec.e_max(e_gamma)),
                E_borde=float(ec.e_max(e_gamma)),
                intervalo=(float(intervalo[0]), float(intervalo[1])),
                rho=float(r_ncc["rho"]), sigma_rel_s2=float(r_ncc["sigma_rel"]),
                barrido_sigma_rel=barrido,
                chi2_ndf=float(np.sum(res ** 2 / var) / gdl),
                rms=float(np.sqrt(np.mean(res ** 2))),
                canales=xw, datos=yw, modelo=modelo_y,
                sa=r_sa, ncc=r_ncc, objetivo=objetivo, pulido=bool(pulir),
                dispersion=dispersion)


# ------------------------------------------------------------- dos bordes ---
def ajustar_dos_bordes(canales, cuentas, eg1, eg2, sigma_rel=None,
                       intervalo=None, objetivo="ec14", pulir=True,
                       delta_e="bordes", semillas_sa=8, ligar_origen=False,
                       ncc_kw=None, sa_kw=None):
    """
    Ec. 11-13 con Theta = {A1, A2, B1, sigma1, sigma2, D}.

    `delta_e`: "bordes" (Ec. 12 corregida) o "gammas" (Ec. 12 tal como se
    imprime). Correr las dos y comparar es la forma de MEDIR la errata.

    `ligar_origen`: impone D = B1/E_max,1, o sea Theta pasa a tener 5
    parametros. NO es una restriccion agregada desde afuera: la Ec. 8 del propio
    paper define el mapa canal-energia POR EL ORIGEN (E_eff = E_e * B/E_max), y
    con ese mapa D queda determinado por B1. Dejar D libre ademas, como hace la
    Ec. 11, sobre-parametriza el modelo respecto de su propia Ec. 8. Con
    resolucion buena da casi lo mismo; con la de este detector, no (ver el
    informe de `correr.py`).
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    ncc_kw = dict(ncc_kw or {})
    sa_kw = dict(sa_kw or {})

    # --- 1. NCC con la plantilla de la linea de MENOR energia -------------
    if sigma_rel is None:
        cand = ncc_mod.barrer_sigma_rel(ch, y, eg1,
                                        SIGMAS_S2,
                                        **ncc_kw)
        r_ncc = max(cand, key=lambda r: r["rho"])
        barrido = [(c["sigma_rel"], c["rho"], c["B"]) for c in cand]
    else:
        r_ncc = ncc_mod.buscar_intervalo(ch, y, eg1, sigma_rel=sigma_rel, **ncc_kw)
        barrido = [(r_ncc["sigma_rel"], r_ncc["rho"], r_ncc["B"])]

    if intervalo is None:
        intervalo = (r_ncc["lo"], r_ncc["hi"])
    xw, yw, _ = _ventana(ch, y, intervalo)

    B10, s10 = r_ncc["B"], r_ncc["sigma"]
    D0 = B10 / float(ec.e_max(eg1))            # mapa por el origen (Ec. 8)
    c1, c2, _ = ec.componentes_dos_bordes(xw, eg1, eg2, 1.0, 1.0, B10, s10,
                                          s10 * float(ec.e_max(eg2) / ec.e_max(eg1)),
                                          D0, delta_e)
    M = np.column_stack([c1, c2])
    a, _ = nnls(M, np.maximum(yw, 0.0))
    A10 = float(a[0]) if a[0] > 0 else 1.0
    A20 = float(a[1]) if a[1] > 0 else A10

    E1 = float(ec.e_max(eg1))

    def _D(p):
        return (p[2] / E1) if ligar_origen else p[5]

    def modelo_fn(x, p):
        return ec.modelo_dos_bordes(x, eg1, eg2, p[0], p[1], p[2], p[3], p[4],
                                    _D(p), delta_e=delta_e)

    p0 = np.array([A10, A20, B10, s10,
                   s10 * float(ec.e_max(eg2) / E1), D0])
    cotas = [(0.0, COTAS_REL["A"][1] * max(A10, A20)),
             (0.0, COTAS_REL["A"][1] * max(A10, A20)),
             (COTAS_REL["B"][0] * B10, COTAS_REL["B"][1] * B10),
             (COTAS_REL["sigma"][0] * s10, COTAS_REL["sigma"][1] * s10),
             (COTAS_REL["sigma"][0] * s10, COTAS_REL["sigma"][1] * s10),
             (COTAS_REL["D"][0] * D0, COTAS_REL["D"][1] * D0)]
    if ligar_origen:
        p0 = p0[:5]
        cotas = cotas[:5]
    costo = (sa_mod.costo_ec14(modelo_fn) if objetivo == "ec14"
             else sa_mod.costo_ncc(modelo_fn))

    # Arranques DISPERSOS sobre las cotas de B1, D y sigma, con A1 y A2 puestas
    # por nnls en cada punto. Con 6 parametros el SA no llega desde un solo
    # arranque (ver el demo), y el diagrama de la Fig. 2 dice explicitamente
    # "modify" y volver: eso es esto.
    barrido_cotas = cotas if not ligar_origen else list(cotas) + [
        (COTAS_REL["D"][0] * D0, COTAS_REL["D"][1] * D0)]
    x0s = [p0]
    for q in semilla_por_barrido(xw, yw, eg1, eg2, barrido_cotas, delta_e,
                                 k=max(int(semillas_sa) - 1, 1),
                                 ligar_origen=ligar_origen):
        x0s.append(q[:5] if ligar_origen else q)
    r_sa = sa_mod.multiarranque(costo, x0s, cotas, xw, yw, **sa_kw)

    p, dispersion = _cerrar(modelo_fn, xw, yw, r_sa, cotas, pulir)

    modelo_y = modelo_fn(xw, p)
    res = yw - modelo_y
    gdl = max(xw.size - p.size, 1)
    var = np.maximum(np.abs(modelo_y), 1.0)
    D_fin = _D(p)
    B2 = p[2] + D_fin * ec._delta_e(eg1, eg2, delta_e)
    return dict(tipo="dos_bordes", e_gamma=(float(eg1), float(eg2)),
                A1=float(p[0]), A2=float(p[1]), B1=float(p[2]),
                sigma1=float(p[3]), sigma2=float(p[4]), D=float(D_fin),
                B2=float(B2), ligar_origen=bool(ligar_origen),
                E_borde=(float(ec.e_max(eg1)), float(ec.e_max(eg2))),
                delta_e=delta_e, dE=float(ec._delta_e(eg1, eg2, delta_e)),
                intervalo=(float(intervalo[0]), float(intervalo[1])),
                rho=float(r_ncc["rho"]), sigma_rel_s2=float(r_ncc["sigma_rel"]),
                barrido_sigma_rel=barrido,
                chi2_ndf=float(np.sum(res ** 2 / var) / gdl),
                rms=float(np.sqrt(np.mean(res ** 2))),
                canales=xw, datos=yw, modelo=modelo_y,
                sa=r_sa, ncc=r_ncc, objetivo=objetivo, pulido=bool(pulir),
                dispersion=dispersion)


# ------------------------------------------------------------------ puntos ---
def puntos_calibracion(res):
    """
    Los pares (canal, energia) que este ajuste aporta a la recta del paper
    ("Fit Energy-Band Structure Relationship").
    """
    if res["tipo"] == "un_borde":
        return [(res["B"], res["E_borde"])]
    return [(res["B1"], res["E_borde"][0]), (res["B2"], res["E_borde"][1])]


def ajustar(canales, cuentas, lineas, **kw):
    """Despacha segun cuantas lineas gamma tenga la fuente."""
    lineas = tuple(float(x) for x in lineas)
    if len(lineas) == 1:
        return ajustar_un_borde(canales, cuentas, lineas[0], **kw)
    if len(lineas) == 2:
        return ajustar_dos_bordes(canales, cuentas, lineas[0], lineas[1], **kw)
    raise ValueError("el paper solo trata fuentes de una o dos lineas gamma")


def _demo():
    """Cierre con verdad conocida, para uno y para dos bordes."""
    from replica_paper import rutas
    rutas.salida_utf8()
    rng = np.random.default_rng(19)

    print("=" * 74)
    print("metodo.py - la cadena de la Fig. 2 sobre datos sinteticos")
    print("=" * 74)

    ch = np.linspace(2.5, 1400.0, 280)
    eg = 661.657
    B_real, s_real = 600.0, 45.0
    limpio = ec.modelo(ch, eg, 1.0, B_real, s_real)
    limpio = limpio / limpio.max() * 700.0
    datos = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
    r = ajustar_un_borde(ch, datos, eg)
    print("  UN BORDE   verdad B = %.1f, sigma = %.1f" % (B_real, s_real))
    print("    NCC : [%6.0f, %6.0f]  rho = %.5f  (sigma_rel de S2 = %.2f)"
          % (r["intervalo"][0], r["intervalo"][1], r["rho"], r["sigma_rel_s2"]))
    print("    SA  : B = %.1f (%+.2f%%)   sigma = %.1f (%+.2f%%)   chi2/ndf = %.3f"
          % (r["B"], 100 * (r["B"] / B_real - 1), r["sigma"],
             100 * (r["sigma"] / s_real - 1), r["chi2_ndf"]))

    print()
    ch = np.linspace(5.0, 2500.0, 250)
    eg1, eg2 = 1173.228, 1332.492
    D_real = 1.30
    B1r = D_real * float(ec.e_max(eg1))
    dB = D_real * ec.delta_e_bordes(eg1, eg2)      # separacion de los dos bordes
    print("  DOS BORDES  verdad B1 = %.1f, D = %.4f" % (B1r, D_real))
    print("  La comparacion entre la Ec. 12 impresa ('gammas') y la corregida")
    print("  ('bordes') se hace con la MISMA ventana y el MISMO sigma de")
    print("  plantilla -si no, se estaria midiendo tambien la eleccion de")
    print("  ventana- y sobre 6 realizaciones, porque el SA es estocastico.")
    print()
    print("  sigma/B   separacion   Ec. 12      B1 [%]             D [%]")
    for s_rel in (0.08, 0.04):
        s1r = s_rel * B1r
        s2r = s1r * float(ec.e_max(eg2) / ec.e_max(eg1))
        limpio = ec.modelo_dos_bordes(ch, eg1, eg2, 1.0, 0.95, B1r, s1r, s2r,
                                      D_real)
        limpio = limpio / limpio.max() * 700.0
        d0 = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
        ref = ncc_mod.buscar_intervalo(ch, d0, eg1, sigma_rel=max(s_rel, 0.05))
        iv, sr = (ref["lo"], ref["hi"]), max(s_rel, 0.05)
        for modo in ("bordes", "gammas"):
            eB, eD = [], []
            for k in range(6):
                d = np.random.default_rng(100 + k).poisson(
                    np.maximum(limpio, 0.0)).astype(float)
                r = ajustar_dos_bordes(ch, d, eg1, eg2, intervalo=iv,
                                       sigma_rel=sr, delta_e=modo)
                eB.append(100 * (r["B1"] / B1r - 1))
                eD.append(100 * (r["D"] / D_real - 1))
            print("   %.2f       %.2f sigma   %-8s  %+6.2f +- %4.2f     "
                  "%+6.2f +- %4.2f"
                  % (s_rel, dB / s1r, modo, np.mean(eB), np.std(eB, ddof=1),
                     np.mean(eD), np.std(eD, ddof=1)))
    print()
    print("  Leer asi:")
    print("   - a 2 sigma de separacion el modelo de dos bordes NO sirve: D sale")
    print("     sesgada ~-14% con una dispersion de +-15%, con cualquiera de las")
    print("     dos versiones de la Ec. 12. Ahi la errata es indistinguible")
    print("     porque el parametro no esta determinado. Es exactamente el caso")
    print("     del Co-60 de este proyecto (2.2-2.5 sigma).")
    print("   - a 4 sigma D queda en -0.6% +- 3.0% con la Ec. 12 corregida y en")
    print("     -2.6% +- 2.7% con la impresa: la diferencia entre ambas (~2%) es")
    print("     del orden del 3% que predice (159.3-154.7)/154.7.")
    print("  Donde la errata SI se ve limpia es en los datos publicados del paper")
    print("  (auditoria.ec12): 74 canales medidos entre los dos bordes del Co-60")
    print("  contra 76.0 que da la Ec. 12 impresa y 73.8 la corregida.")


if __name__ == "__main__":
    _demo()
