"""
sa.py — Seccion II-A del paper: recocido simulado, con SUS constantes.

Lo que dice el paper, y esta implementado literal:

    "The SA process initiates with parameter vector x0 at T_init = 2000 and
     cooling coefficient alpha = 0.95. During each temperature stage T_k, the
     algorithm executes 100 Markov chain iterations by generating new states
     x_new = x_current + dx where dx_i ~ U(-0.1 d_i, 0.1 d_i) with d_i defining
     the parameter search range. State transitions follow the Metropolis
     criterion: new configurations are accepted with probability
     P = exp(-dE/T_k) where dE = E(x_new) - E(x_current) and E(x) represents the
     NCC-based matching error objective function. After completing Markov chains
     at current temperature, the system cools geometrically via
     T_k+1 = alpha*T_k (alpha in [0.85, 0.99]). This thermodynamic cycle
     iterates until meeting dual termination criteria: T_k < 1e-6 or energy
     stabilization criteria (|dE_mean| < 0.001 over 10 consecutive cycles)."

Inconsistencia del paper, no resuelta por ellos
-----------------------------------------------
La Ec. 14 define la loss como minimos cuadrados, y el parrafo de Metropolis dice
que E(x) es el "NCC-based matching error". Son objetivos distintos y no se
pueden mezclar, asi que hay que ELEGIR. Aca los dos estan implementados
(`costo_ec14` y `costo_ncc`) y el que se usa es un argumento, para poder mostrar
la diferencia. Ojo con el costo NCC: Pearson ignora la escala, asi que con ese
objetivo las amplitudes A1, A2 son literalmente INDETERMINABLES.

Que significa T_init = 2000
---------------------------
Metropolis compara T con dE, o sea T vive en las UNIDADES DEL COSTO. Con la
Ec. 14 sobre un espectro real la suma de cuadrados es del orden de 1e7-1e9, asi
que T = 2000 acepta practicamente nada desde el primer paso; con el costo NCC
(rango [0, 2]) T = 2000 acepta TODO durante cientos de ciclos. La constante del
paper no es transferible entre objetivos y por eso `t_init` se puede pasar como
"auto": se calibra midiendo la dispersion del costo alrededor de x0. Lo que
CORRE POR DEFECTO es el 2000 del paper; "auto" hay que pedirlo.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import ecuaciones as ec

T_INIT = 2000.0         # del paper
ALFA = 0.95             # del paper
ITERS = 100             # "100 Markov chain iterations" por temperatura
PASO_REL = 0.1          # dx_i ~ U(-0.1 d_i, 0.1 d_i)
T_MIN = 1e-6            # T_k < 1e-6
TOL_E = 1e-3            # |dE_mean| < 0.001
CICLOS_TOL = 10         # "over 10 consecutive cycles"


def costo_ec14(modelo_fn):
    """Objetivo de la Ec. 14: suma de cuadrados de los residuos."""
    def costo(p, x, y):
        return ec.loss_minimos_cuadrados(y, modelo_fn(x, p))
    return costo


def costo_ncc(modelo_fn):
    """
    Objetivo "NCC-based matching error" del parrafo de Metropolis: 1 - rho.

    Con este objetivo las amplitudes no tienen efecto sobre el costo (Pearson es
    invariante a escala), asi que quedan indeterminadas. Es exactamente la
    ambiguedad que el paper deja abierta.
    """
    def costo(p, x, y):
        return 1.0 - ec.pearson(y, modelo_fn(x, p))
    return costo


def _calibrar_t(costo, x0, cotas, x, y, n=60, aceptacion=0.8, rng=None):
    """
    T_init tal que una subida tipica se acepte con probabilidad `aceptacion`.

    NO es del paper: es la unica manera de que el 2000 impreso signifique lo
    mismo con los dos objetivos. Se ofrece como opcion, no como default.
    """
    rng = np.random.default_rng() if rng is None else rng
    d = np.array([b - a for a, b in cotas], float)
    e0 = costo(np.asarray(x0, float), x, y)
    subidas = []
    for _ in range(n):
        p = np.clip(np.asarray(x0, float) + rng.uniform(-PASO_REL, PASO_REL, d.size) * d,
                    [a for a, _ in cotas], [b for _, b in cotas])
        de = costo(p, x, y) - e0
        if de > 0:
            subidas.append(de)
    if not subidas:
        return T_INIT
    return float(np.median(subidas) / -np.log(aceptacion))


def recocido(costo, x0, cotas, x, y, t_init=T_INIT, alfa=ALFA, iters=ITERS,
             paso_rel=PASO_REL, t_min=T_MIN, tol=TOL_E, ciclos_tol=CICLOS_TOL,
             max_ciclos=100000, semilla=None, historia=False):
    """
    Recocido simulado con los parametros del paper.

    `cotas` es la lista de (min, max) por parametro; su ancho es el d_i del
    paper ("the parameter search range").

    Criterio de Metropolis estandar: una bajada se acepta siempre, una subida
    con probabilidad exp(-dE/T). El paper escribe solo P = exp(-dE/T_k), que
    para dE < 0 daria P > 1; se interpreta como min(1, ...), que es lo que
    "Metropolis criterion" significa.

    Criterio de corte por energia: el paper dice |dE_mean| < 0.001 sin definir
    dE_mean. Aca se toma la diferencia entre el costo MEDIO de un ciclo de
    temperatura y el del ciclo anterior, que es la lectura literal de "energy
    stabilization ... over 10 consecutive cycles".
    """
    rng = np.random.default_rng(semilla)
    lo = np.array([a for a, _ in cotas], float)
    hi = np.array([b for _, b in cotas], float)
    d = hi - lo
    xc = np.clip(np.asarray(x0, float), lo, hi)
    ec_actual = costo(xc, x, y)
    mejor_x, mejor_e = xc.copy(), ec_actual

    if isinstance(t_init, str) and t_init == "auto":
        T = _calibrar_t(costo, xc, cotas, x, y, rng=rng)
    else:
        T = float(t_init)

    aceptados = intentos = 0
    media_previa = None
    estables = 0
    hist = []
    ciclos = 0
    while T > t_min and ciclos < max_ciclos:
        suma = 0.0
        for _ in range(int(iters)):
            prop = np.clip(xc + rng.uniform(-paso_rel, paso_rel, d.size) * d, lo, hi)
            e_new = costo(prop, x, y)
            de = e_new - ec_actual
            intentos += 1
            if de <= 0.0 or rng.random() < np.exp(-de / T):
                xc, ec_actual = prop, e_new
                aceptados += 1
                if ec_actual < mejor_e:
                    mejor_x, mejor_e = xc.copy(), ec_actual
            suma += ec_actual
        media = suma / float(iters)
        ciclos += 1
        if historia:
            hist.append((T, media, mejor_e))
        if media_previa is not None and abs(media - media_previa) < tol:
            estables += 1
            if estables >= int(ciclos_tol):
                break
        else:
            estables = 0
        media_previa = media
        T *= float(alfa)

    return dict(x=mejor_x, costo=float(mejor_e), T_final=float(T),
                ciclos=int(ciclos), evaluaciones=int(intentos),
                aceptacion=float(aceptados / max(intentos, 1)),
                t_init=float(T if ciclos == 0 else T / alfa ** ciclos),
                historia=np.array(hist, float) if historia else None,
                corte=("T < %.0e" % t_min) if T <= t_min else
                      ("|dE_mean| < %.0e en %d ciclos" % (tol, ciclos_tol)))


def multiarranque(costo, x0s, cotas, x, y, semilla=0, **kw):
    """
    El SA es estocastico: una sola corrida no es un resultado.

    Corre desde varios arranques y devuelve el mejor, mas la lista completa para
    poder ver la dispersion. [Esto no esta en el paper, que reporta una corrida.]
    """
    salidas = []
    for k, x0 in enumerate(x0s):
        salidas.append(recocido(costo, x0, cotas, x, y, semilla=semilla + k, **kw))
    mejor = min(salidas, key=lambda r: r["costo"])
    mejor = dict(mejor)
    mejor["todas"] = salidas
    return mejor


def _demo():
    """
    Cierre con verdad conocida y, de paso, la medida de que T_init = 2000 no es
    transferible entre los dos objetivos que el paper menciona.
    """
    from replica_paper import rutas
    rutas.salida_utf8()
    rng = np.random.default_rng(11)

    eg = 661.657
    B_real, sigma_real = 600.0, 40.0
    ch = np.linspace(2.5, 1200.0, 240)
    limpio = ec.modelo(ch, eg, 1.0, B_real, sigma_real)
    escala = 800.0 / limpio.max()
    A_real = escala
    datos = rng.poisson(np.maximum(limpio * escala, 0.0)).astype(float)

    def modelo_fn(canales, p):
        return ec.modelo(canales, eg, p[0], p[1], p[2])

    cotas = [(1e-4 * A_real, 10 * A_real), (0.3 * B_real, 2.0 * B_real),
             (0.005 * B_real, 0.5 * B_real)]
    # arranque MALO a proposito: es el caso en que un optimizador local falla
    x0 = [3.0 * A_real, 0.70 * B_real, 3.0 * sigma_real]

    print("=" * 74)
    print("sa.py - Seccion II-A con las constantes del paper")
    print("=" * 74)
    print("  verdad: A = %.4g, B = %.1f, sigma = %.1f" % (A_real, B_real, sigma_real))
    print("  x0    : A = %.4g, B = %.1f, sigma = %.1f (arranque malo a proposito)"
          % tuple(x0))
    print()
    print("  objetivo       T_init   aceptacion  ciclos      B      sigma   corte")
    for nombre, cfn, tin in (("Ec. 14", costo_ec14(modelo_fn), T_INIT),
                             ("Ec. 14", costo_ec14(modelo_fn), "auto"),
                             ("NCC", costo_ncc(modelo_fn), T_INIT),
                             ("NCC", costo_ncc(modelo_fn), "auto")):
        r = recocido(cfn, x0, cotas, ch, datos, t_init=tin, semilla=3)
        print("  %-8s  %8s   %8.3f   %5d   %7.1f  %6.1f   %s"
              % (nombre, tin if isinstance(tin, str) else "%.0f" % tin,
                 r["aceptacion"], r["ciclos"], r["x"][1], r["x"][2], r["corte"]))
    print()
    print("  Leer asi: T = 2000 con la Ec. 14 (costo ~1e6) casi no acepta nada,")
    print("  y con el costo NCC (rango [0,2]) acepta todo. La constante impresa")
    print("  solo tiene sentido en las unidades del costo que se elija.")


if __name__ == "__main__":
    _demo()
