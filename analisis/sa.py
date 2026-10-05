"""
sa.py — Recocido Simulado (Simulated Annealing), Seccion II-A del paper.

Que hace y por que
------------------
Un optimizador local (Nelder-Mead, Levenberg-Marquardt) solo sabe bajar: desde
donde arranca, camina cuesta abajo hasta el fondo del valle en el que cayo. Si
la funcion de costo tiene varios valles, el resultado depende del punto inicial
y no hay manera de saber, desde adentro, si el valle en el que se quedo es el
mas profundo.

El recocido acepta A VECES subir. Con probabilidad `exp(-DeltaE/T)` se queda con
una propuesta que empeora el costo, y esa es la unica manera de salir de un
valle. La temperatura T controla cuanto: alta al principio (el algoritmo pasea y
explora), baja despues (se comporta como un optimizador local y afina). El
nombre viene de la metalurgia: enfriar despacio deja que los atomos encuentren
la red de minima energia, enfriar de golpe los congela en un vidrio.

Lo que dice el paper (Seccion II-A)
-----------------------------------
    x_nuevo = x_actual + Dx,  con  Dx_i ~ U(-0.1*d_i, +0.1*d_i)
    P(aceptar) = exp(-DeltaE / T_k)          criterio de Metropolis
    100 iteraciones de cadena de Markov por cada temperatura
    T_{k+1} = alpha * T_k,  alpha en [0.85, 0.99]  (usan 0.95)
    T_init = 2000
    corte:  T_k < 1e-6   o   |DeltaE_medio| < 0.001 durante 10 ciclos seguidos

Dos de esos numeros NO son transferibles, y el modulo los reemplaza:

  * `T_init = 2000`. Metropolis compara T contra DeltaE, asi que T esta en las
    UNIDADES DEL COSTO. Medido (verificacion 3): sobre nuestro chi2, que vale
    ~470 en la semilla, T = 2000 da una aceptacion inicial de 0.76 y anda
    perfectamente bien; sobre EL MISMO paisaje con el costo dividido por 470
    -la escala que tendria un error de correlacion tipo NCC- el mismo T = 2000
    acepta el 100% y la cadena es una caminata al azar. O sea: el 2000 del paper
    no esta mal, esta atado a la escala de SU costo. `T_init="auto"` lo calibra
    midiendo cuanto empeoran unas propuestas al azar y despejando la T que da la
    aceptacion inicial pedida, asi el numero deja de depender de las unidades.
  * el corte absoluto `|DeltaE| < 0.001`, por lo mismo. Se usa un corte RELATIVO.

Hay ademas un paso ADAPTATIVO disponible (escalar el paso para mantener la
aceptacion cerca de un objetivo), pero NO es el default y conviene decir por que.
La version anterior de este docstring afirmaba que hacia converger donde el paso
fijo se queda pegado. Eso venia del proyecto anterior, sobre un modelo de 4-6
parametros, y NO se sostiene aca: medido con 8 semillas sobre este problema
(`verificar_sa.py`, verificacion 4), el paso fijo del paper da
D = 1.2219 +- 0.0042 y el adaptativo D = 1.2197 +- 0.0046. La dispersion es la
misma. Como no hay evidencia a favor, el default es el del PAPER (paso fijo);
`paso_adaptativo=True` queda disponible para modelos con mas parametros, donde
la escala entre ellos si es un problema (ver verificacion 6).

Correr `python analisis/sa.py` deja `figuras/07_sa.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas


# --------------------------------------------------------------------------
# 1. Calibracion de la temperatura inicial
# --------------------------------------------------------------------------
def temperatura_inicial(costo, x0, delta, aceptacion=0.8, n=60, rng=None):
    """
    T tal que un empeoramiento TIPICO se acepte con probabilidad `aceptacion`.

    Se muestrean `n` propuestas alrededor de x0, se mide cuanto empeoran en
    promedio, y se despeja T de  exp(-<DeltaE+> / T) = aceptacion:

        T = -<DeltaE+> / ln(aceptacion)

    Sin esto hay que adivinar un numero que depende de las unidades del costo, y
    adivinar mal no da error: da un recocido que no explora (T chico) o que
    nunca converge (T grande).
    """
    rng = np.random.default_rng() if rng is None else rng
    x0 = np.asarray(x0, float)
    c0 = costo(x0)
    subidas = []
    for _ in range(n):
        x = x0 + rng.uniform(-1, 1, size=x0.size) * delta
        d = costo(x) - c0
        if d > 0:
            subidas.append(d)
    if not subidas:
        return max(abs(c0) * 1e-3, 1e-9)
    return float(-np.mean(subidas) / np.log(aceptacion))


# --------------------------------------------------------------------------
# 2. El recocido
# --------------------------------------------------------------------------
def recocido(costo, x0, cotas, delta=None, T_init="auto", alfa=0.95,
             n_markov=100, T_min_rel=1e-6, n_estable=10, tol_rel=1e-4,
             paso_adaptativo=False, objetivo_aceptacion=0.35, semilla=None,
             max_ciclos=400):
    """
    Minimiza `costo` con recocido simulado. Devuelve un dict.

    costo  : callable(x) -> float
    x0     : punto inicial (la semilla del NCC)
    cotas  : [(min, max), ...] por parametro. Acota la busqueda a lo fisicamente
             posible; ademas de ahi sale `delta` si no se lo dan.
    delta  : escala de paso por parametro. Por defecto 0.1 * (max - min), que es
             el `0.1*d_i` del paper.
    alfa   : coeficiente de enfriamiento geometrico, T_{k+1} = alfa * T_k.
    n_markov : iteraciones por temperatura (las 100 del paper).
    T_min_rel : corte por temperatura, relativo a T_init.
    n_estable, tol_rel : corte por estabilizacion. El mejor costo tiene que
             mejorar menos de `tol_rel` (relativo) durante `n_estable` ciclos.

    Claves del resultado: x, costo, historia (una fila por ciclo de temperatura
    con T, costo actual, mejor costo, tasa de aceptacion y paso), n_eval.

    Se guarda ademas el CAMINO: `historia_x` es donde estaba la cadena al
    cerrar cada ciclo de temperatura y `historia_mejor_x` el mejor punto
    conocido hasta ese ciclo. Cuesta una fila de floats por ciclo (unos
    cientos en total) y es lo unico que permite ver despues por donde paseo el
    recocido, que es la mitad de lo que hay que mirar para creerle: `historia`
    dice que el costo bajo, `historia_x` dice si para eso hubo que cruzar de
    un valle a otro o alcanzo con bajar en linea recta.
    """
    rng = np.random.default_rng(semilla)
    x = np.asarray(x0, float).copy()
    cotas = np.asarray(cotas, float)
    lo, hi = cotas[:, 0], cotas[:, 1]
    if delta is None:
        delta = 0.1 * (hi - lo)
    delta = np.asarray(delta, float).copy()

    T = (temperatura_inicial(costo, x, delta, rng=rng)
         if T_init == "auto" else float(T_init))
    T0 = T
    escala = 1.0

    c = costo(x)
    c_inicial = c
    mejor_x, mejor_c = x.copy(), c
    historia = []
    historia_x = []
    historia_mejor_x = []
    estables = 0
    n_eval = 1

    for ciclo in range(max_ciclos):
        aceptados = 0
        for _ in range(n_markov):
            prop = x + rng.uniform(-1, 1, size=x.size) * delta * escala
            prop = np.clip(prop, lo, hi)
            cp = costo(prop)
            n_eval += 1
            dE = cp - c
            # Metropolis: si mejora se acepta siempre; si empeora, con exp(-dE/T).
            if dE <= 0 or rng.random() < np.exp(-dE / max(T, 1e-300)):
                x, c = prop, cp
                aceptados += 1
                if c < mejor_c:
                    mejor_x, mejor_c = x.copy(), c
        tasa = aceptados / n_markov
        historia.append((T, c, mejor_c, tasa, escala))
        historia_x.append(x.copy())
        historia_mejor_x.append(mejor_x.copy())

        # paso adaptativo [fuera del paper, apagado por defecto]: si se acepta
        # mucho, el paso es demasiado timido y conviene agrandarlo; si se acepta
        # poco, achicarlo. Sobre ESTE problema no cambia nada (verificacion 4),
        # asi que no se usa; se deja para modelos con mas parametros.
        # El objetivo 0.35 no es un numero canonico: los valores citados en la
        # literatura de MUESTREO por Metropolis son 0.44 (1 dimension) y 0.234
        # (asintotico en muchas). Barrido entre 0.15 y 0.60, D no se mueve mas
        # de 0.011% (verificacion 5).
        if paso_adaptativo:
            if tasa > objetivo_aceptacion + 0.1:
                escala *= 1.25
            elif tasa < objetivo_aceptacion - 0.1:
                escala *= 0.8
            escala = float(np.clip(escala, 1e-4, 20.0))

        # corte por estabilizacion (relativo, no absoluto)
        if len(historia) > 1:
            prev = historia[-2][2]
            if abs(prev - mejor_c) <= tol_rel * max(abs(prev), 1e-12):
                estables += 1
            else:
                estables = 0
            if estables >= n_estable:
                break

        T *= alfa
        if T < T_min_rel * T0:
            break

    return dict(x=mejor_x, costo=mejor_c, historia=np.array(historia),
                historia_x=np.array(historia_x),
                historia_mejor_x=np.array(historia_mejor_x),
                x_inicial=np.asarray(x0, float).copy(),
                costo_inicial=float(c_inicial), delta=delta.copy(),
                n_eval=n_eval, T_init=T0, ciclos=len(historia),
                motivo=("estabilizacion" if estables >= n_estable
                        else "temperatura minima"))


def pulir(costo, x, cotas=None):
    """
    Optimizador local desde donde dejo el recocido.

    El recocido encuentra el valle bueno pero se detiene con la precision del
    ultimo paso, que no es cero. El pulido baja al fondo exacto. Los dos hacen
    falta: el recocido sin pulido da el 3er decimal, el pulido sin recocido da
    el valle equivocado.
    """
    from scipy.optimize import minimize
    r = minimize(costo, np.asarray(x, float), method="Nelder-Mead",
                 options=dict(xatol=1e-7, fatol=1e-4, maxiter=4000))
    return r.x, float(r.fun)


def multiarranque(costo, cotas, n=8, x0=None, dispersion=0.25, semilla=0,
                  callback=None, **kw):
    """
    Corre `n` recocidos desde puntos iniciales DISPERSOS y devuelve todos.

    El recocido es estocastico: una sola corrida no es un resultado, es una
    muestra. Que N arranques independientes y dispersos terminen en el mismo
    punto es la evidencia practica de que ese punto es el minimo global (no una
    demostracion: ver el PDF).
    """
    rng = np.random.default_rng(semilla)
    cotas = np.asarray(cotas, float)
    salidas = []
    for i in range(n):
        if callback is not None:
            # el avance se cuenta por ARRANQUE y no por ciclo de temperatura:
            # cuantos ciclos va a hacer cada recocido no se sabe de antemano
            # (corta por estabilizacion), asi que inventar ese denominador daria
            # una barra que retrocede.
            callback(i / n, "arranque %d de %d" % (i + 1, n))
        if x0 is None:
            xi = rng.uniform(cotas[:, 0], cotas[:, 1])
        else:
            xi = np.asarray(x0, float) * (1 + dispersion *
                                          rng.uniform(-1, 1, size=len(x0)))
            xi = np.clip(xi, cotas[:, 0], cotas[:, 1])
        r = recocido(costo, xi, cotas, semilla=int(rng.integers(1 << 30)), **kw)
        r["x_pulido"], r["costo_pulido"] = pulir(costo, r["x"])
        r["x0"] = xi
        salidas.append(r)
    return salidas


# --------------------------------------------------------------------------
# 3. Como saber si el optimo es bueno
# --------------------------------------------------------------------------
def errores_por_curvatura(costo, x, paso_rel=1e-3, delta_chi2=1.0, paso=None):
    """
    Barras de error a partir de la CURVATURA del costo en el minimo.

    Para un chi2, la matriz de covarianza es `2 * H^-1` con H el hessiano de
    chi2/2, o sea `C = 2 * (d2chi2/dxdx)^-1`. El error de cada parametro es la
    raiz de la diagonal, y corresponde al contorno `chi2_min + 1`.

    Vale solo si el minimo es aproximadamente cuadratico, que es justamente lo
    que hay que mirar en el mapa del paisaje antes de creerle.

    `paso` permite dar el paso de derivacion de cada parametro a mano, en vez
    del `paso_rel * |x|` de siempre. Hace falta cuando la escala natural de un
    parametro no es su propio valor. El caso que lo obligo: la ordenada c0 del
    mapa canal-energia. La plantilla se lleva a canal con `np.interp`, o sea
    que es lineal a trozos sobre la grilla del histograma; con c0 = -163 el
    paso relativo da 0.16 mV*ns, un 3% de un bin de 5.4, y ahi la "curvatura"
    que se mide es el quiebre entre dos trozos y no la del chi2. Daba 21 sigma
    donde el cociente de verosimilitud da 2.1. Para ese parametro el paso tiene
    que ser del orden de UN BIN.

    Control: para D, cuya escala si es su propio valor, el paso no cambia nada
    (err_D = 0.0040 con paso_rel 1e-2 y 0.0044 con 3e-4) y coincide con el
    perfil real de Delta_chi2 = 1, que da 0.0040.
    """
    x = np.asarray(x, float)
    n = x.size
    if paso is not None:
        h = np.abs(np.asarray(paso, float))
    else:
        h = np.maximum(np.abs(x) * paso_rel, 1e-9)
    H = np.empty((n, n))
    c0 = costo(x)
    for i in range(n):
        for j in range(n):
            xp, xm = x.copy(), x.copy()
            if i == j:
                xp[i] += h[i]; xm[i] -= h[i]
                H[i, i] = (costo(xp) - 2 * c0 + costo(xm)) / h[i] ** 2
            else:
                a = x.copy(); a[i] += h[i]; a[j] += h[j]
                b = x.copy(); b[i] += h[i]; b[j] -= h[j]
                cc = x.copy(); cc[i] -= h[i]; cc[j] += h[j]
                d = x.copy(); d[i] -= h[i]; d[j] -= h[j]
                H[i, j] = (costo(a) - costo(b) - costo(cc) + costo(d)) / (4 * h[i] * h[j])
    try:
        C = 2.0 * np.linalg.inv(H)
        err = np.sqrt(np.maximum(np.diag(C), 0.0)) * np.sqrt(delta_chi2)
        return err, C
    except np.linalg.LinAlgError:
        return np.full(n, np.nan), np.full((n, n), np.nan)


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from analisis import modelo

    rutas.salida_utf8()
    clave = "co60"
    print("=" * 74)
    print("Recocido simulado sobre %s" % clave)
    print("=" * 74)
    p, semilla, r, e, c = modelo.armar(clave)
    print("  OJO: esta demo usa el modelo DESNUDO de 2 parametros (D, sigma),")
    print("  que es el paisaje sobre el que se explica el recocido. La")
    print("  calibracion de verdad agrega el exceso de bajo canal y da otro D")
    print("  (1.2512 en vez de 1.2230): ver calibracion.EXCESO.")

    cotas = [(0.6 * semilla[0], 1.6 * semilla[0]), (0.02, 0.30)]

    # --- una corrida, desde la semilla del NCC ---------------------------
    res = recocido(p.costo, semilla, cotas, semilla=1)
    xp, cp = pulir(p.costo, res["x"])
    print("  T_init calibrada = %.3g   (el paper dice 2000, que aca seria "
          "%s)" % (res["T_init"],
                   "absurdo" if res["T_init"] < 20 or res["T_init"] > 2e5 else "razonable"))
    print("  %d ciclos, %d evaluaciones, corte por %s"
          % (res["ciclos"], res["n_eval"], res["motivo"]))
    print("  SA          : D = %.4f  sigma = %.4f  chi2/ndf = %.4f"
          % (res["x"][0], res["x"][1], res["costo"] / p.ndf))
    print("  SA + pulido : D = %.4f  sigma = %.4f  chi2/ndf = %.4f"
          % (xp[0], xp[1], cp / p.ndf))

    err, C = errores_por_curvatura(p.costo, xp)
    corr = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
    print("  errores por curvatura: D = %.4f +- %.4f (%.2f%%), sigma = %.4f "
          "+- %.4f,  corr(D,sigma) = %+.2f"
          % (xp[0], err[0], 100 * err[0] / xp[0], xp[1], err[1], corr))

    # --- multiarranque desde puntos al azar ------------------------------
    print("\n  multiarranque (8 puntos iniciales AL AZAR dentro de las cotas):")
    ms = multiarranque(p.costo, cotas, n=8, semilla=7)
    Ds = np.array([m["x_pulido"][0] for m in ms])
    cs = np.array([m["costo_pulido"] for m in ms])
    for m in ms:
        print("     x0 = (%.3f, %.3f)  ->  D = %.4f  chi2 = %.2f"
              % (m["x0"][0], m["x0"][1], m["x_pulido"][0], m["costo_pulido"]))
    print("     D: mediana %.4f, dispersion %.2g (%.3f%%)"
          % (np.median(Ds), Ds.std(), 100 * Ds.std() / np.median(Ds)))
    print("     %d de %d corridas llegan al mismo chi2 dentro de 0.1"
          % (int(np.sum(cs < cs.min() + 0.1)), len(cs)))

    # --- optimizador local desde un arranque malo ------------------------
    from scipy.optimize import minimize
    malo = (0.70 * xp[0], 0.28)
    loc = minimize(p.costo, malo, method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-4, maxiter=4000))
    sa_malo = recocido(p.costo, malo, cotas, semilla=3)
    xs, csa = pulir(p.costo, sa_malo["x"])
    print("\n  desde un arranque MALO (D = %.3f, sigma = %.2f):" % malo)
    print("     solo local : D = %.4f  chi2 = %.2f" % (loc.x[0], loc.fun))
    print("     SA + pulido: D = %.4f  chi2 = %.2f" % (xs[0], csa))

    # ---------------------------------------------------------- figuras --
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9))

    # (a) paisaje
    Ds_g = np.linspace(cotas[0][0], cotas[0][1], 110)
    Ss_g = np.linspace(cotas[1][0], cotas[1][1], 90)
    Z = p.mapa(Ds_g, Ss_g)
    a = ax[0][0]
    im = a.pcolormesh(Ds_g, Ss_g, np.log10(Z / p.ndf), shading="auto",
                      cmap="viridis")
    a.contour(Ds_g, Ss_g, Z, levels=[cp + 1, cp + 4, cp + 9], colors="w", linewidths=.8)
    a.plot(*xp, "o", color="C3", ms=7, label="optimo")
    a.plot(*semilla, "s", color="C1", ms=6, label="semilla NCC")
    a.set_title("el paisaje: log10(chi2/ndf)")
    a.set_xlabel("D [mV*ns/keVee]"); a.set_ylabel("sigma_rel")
    a.legend(fontsize=8); fig.colorbar(im, ax=a)

    # (b) Metropolis
    a = ax[0][1]
    dE = np.linspace(0, 400, 300)
    for T, col in ((res["T_init"], "C0"), (res["T_init"] * 0.1, "C1"),
                   (res["T_init"] * 0.01, "C2")):
        a.plot(dE, np.exp(-dE / T), color=col, lw=1.4, label="T = %.3g" % T)
    a.set_title("criterio de Metropolis: P = exp(-dE / T)")
    a.set_xlabel("cuanto EMPEORA el costo (dE, unidades de chi2)")
    a.set_ylabel("probabilidad de aceptar igual")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # (c) enfriamiento
    a = ax[0][2]
    k = np.arange(150)
    for al in (0.85, 0.95, 0.99):
        a.semilogy(k, res["T_init"] * al ** k, lw=1.3, label="alfa = %.2f" % al)
    a.set_title("enfriamiento geometrico T_k = alfa^k * T_0")
    a.set_xlabel("ciclo de temperatura"); a.set_ylabel("T")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # (d) la corrida
    h = res["historia"]
    a = ax[1][0]
    a.semilogy(h[:, 1] / p.ndf, lw=1, color="C0", label="costo actual")
    a.semilogy(h[:, 2] / p.ndf, lw=1.6, color="C3", label="mejor hasta aqui")
    a.set_title("la corrida: costo por ciclo de temperatura")
    a.set_xlabel("ciclo"); a.set_ylabel("chi2/ndf")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # (e) aceptacion
    a = ax[1][1]
    a.plot(h[:, 3], lw=1.3, color="C0", label="tasa de aceptacion")
    a.axhline(0.35, color="k", ls="--", lw=1, label="objetivo 0.35")
    a2 = a.twinx()
    a2.semilogy(h[:, 4], lw=1.2, color="C1", label="escala del paso")
    a2.set_ylabel("escala del paso", color="C1")
    a.set_title("aceptacion y paso adaptativo")
    a.set_xlabel("ciclo"); a.set_ylabel("fraccion aceptada")
    a.legend(fontsize=8, loc="upper right"); a.grid(alpha=.3)

    # (f) multiarranque
    a = ax[1][2]
    a.scatter([m["x0"][0] for m in ms], Ds, s=45, c="C0", zorder=3,
              label="8 arranques al azar")
    a.axhline(xp[0], color="C3", ls="--", label="D = %.4f" % xp[0])
    a.set_title("multiarranque: todos caen en el mismo minimo")
    a.set_xlabel("D del punto INICIAL"); a.set_ylabel("D FINAL")
    a.legend(fontsize=8); a.grid(alpha=.3)

    fig.tight_layout()
    salida = rutas.figura("07_sa.png")
    fig.savefig(salida, dpi=130)
    print("\nfigura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
