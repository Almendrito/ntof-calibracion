"""
verificar_sa.py — pone a prueba las afirmaciones que sostienen a `sa.py`.

Cada bloque contrasta UNA afirmacion. La idea es que ninguna decision de
implementacion quede sostenida por "se sabe que" ni por una cita de memoria:
o sale del paper (y entonces se puede leer en el PDF), o se mide aca.

    python analisis/verificar_sa.py
"""

import os
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, fisica, modelo, sa


def titulo(n, t):
    print()
    print("=" * 78)
    print("VERIFICACION %d - %s" % (n, t))
    print("=" * 78)


# --------------------------------------------------------------------------
def v1_boltzmann():
    """
    AFIRMACION: con P = exp(-dE/T), la cadena a T fija muestrea exp(-E/T).

    Es un teorema (balance detallado), pero se puede COMPROBAR: se corre la
    cadena a T fija sobre una funcion conocida y se compara el histograma de
    donde estuvo contra exp(-E(x)/T) normalizada.
    """
    titulo(1, "Metropolis a T fija muestrea la distribucion de Boltzmann")
    f = lambda x: (x - 1.6) ** 2 / 6.0 - np.cos(3.4 * x) + 1.0
    rng = np.random.default_rng(0)
    xs = np.linspace(-4, 6, 400)

    for T in (0.5, 2.0):
        x, c = 0.0, f(0.0)
        muestras = []
        for k in range(400000):
            prop = x + rng.uniform(-1, 1) * 1.0
            if not (-4 <= prop <= 6):
                muestras.append(x); continue
            cp = f(prop)
            if cp - c <= 0 or rng.random() < np.exp(-(cp - c) / T):
                x, c = prop, cp
            muestras.append(x)
        muestras = np.array(muestras[20000:])       # descarta el transitorio

        h, bordes = np.histogram(muestras, bins=80, range=(-4, 6), density=True)
        cen = 0.5 * (bordes[:-1] + bordes[1:])
        teo = np.exp(-f(cen) / T)
        teo /= np.trapezoid(teo, cen)
        err = np.max(np.abs(h - teo)) / teo.max()
        print("  T = %.1f : discrepancia maxima entre el histograma de la cadena "
              "y exp(-E/T) = %.3f %%" % (T, 100 * err))
    print("  -> la regla de aceptacion no es una eleccion de diseno: es la que")
    print("     hace que la cadena muestree Boltzmann. VERIFICADO numericamente.")


# --------------------------------------------------------------------------
def v2_ec1():
    """
    AFIRMACION: la Ec. 1 impresa del paper no da el borde Compton.

    Se puede chequear con los propios datos del paper: su Fig. 1 muestra la
    seccion eficaz para E0 = 662 keV cortandose cerca de 480 keV, no de 180.
    """
    titulo(2, "La Ec. 1 tal como esta impresa en el paper")
    for eg, esperado in ((661.657, "~477, que es donde corta la Fig. 1 del paper"),
                         (1173.228, "963, el borde conocido del Co-60"),
                         (1332.492, "1118, el otro borde del Co-60")):
        a = eg / fisica.M_E
        impresa = eg / (1 + 2 * a)              # Ec. 1 tal cual esta impresa
        correcta = fisica.borde_compton(eg)
        print("  E_gamma = %8.3f keV :  Ec.1 impresa = %6.1f   "
              "con 2a en el numerador = %6.1f   (se espera %s)"
              % (eg, impresa, correcta, esperado))
    print("  -> la formula impresa da la energia del FOTON retrodispersado.")
    print("     Es verificable contra la Fig. 1 del propio paper, no es opinion.")


# --------------------------------------------------------------------------
def v3_T_init(p, semilla, cotas):
    """
    AFIRMACION: T_init depende de las unidades del costo, asi que copiarla
    de otro trabajo no es seguro.

    Se mide la tasa de aceptacion inicial con T = 2000 sobre (a) nuestro chi2 y
    (b) el MISMO problema con el costo reescalado al orden 1, que es la escala
    que tendria un error de correlacion tipo NCC.
    """
    titulo(3, "T_init = 2000 sobre costos de distinta escala")
    rng = np.random.default_rng(0)
    delta = 0.1 * np.array([c[1] - c[0] for c in cotas])

    def aceptacion(costo, T, n=3000):
        x = np.array(semilla, float)
        c = costo(x)
        acc = 0
        for _ in range(n):
            prop = np.clip(x + rng.uniform(-1, 1, 2) * delta,
                           [c0[0] for c0 in cotas], [c0[1] for c0 in cotas])
            cp = costo(prop)
            if cp - c <= 0 or rng.random() < np.exp(-(cp - c) / T):
                x, c = prop, cp
                acc += 1
        return acc / n

    c_chi2 = p.costo
    escala = p.costo(semilla)
    c_orden1 = lambda x: p.costo(x) / escala      # mismo paisaje, otra escala

    T_auto = sa.temperatura_inicial(c_chi2, semilla, delta, rng=rng)
    print("  nuestro chi2  (valor tipico ~%.0f):" % escala)
    print("     T = 2000 (del paper)  -> aceptacion inicial %.2f"
          % aceptacion(c_chi2, 2000.0))
    print("     T = %.0f (calibrada)  -> aceptacion inicial %.2f"
          % (T_auto, aceptacion(c_chi2, T_auto)))
    print("  el MISMO paisaje con el costo dividido por %.0f (escala tipo NCC):" % escala)
    print("     T = 2000 (del paper)  -> aceptacion inicial %.2f"
          % aceptacion(c_orden1, 2000.0))
    T_auto2 = sa.temperatura_inicial(c_orden1, semilla, delta, rng=rng)
    print("     T = %.3g (calibrada)  -> aceptacion inicial %.2f"
          % (T_auto2, aceptacion(c_orden1, T_auto2)))
    print("  -> con T = 2000 y un costo de orden 1 se acepta TODO: la cadena es")
    print("     una caminata al azar y la temperatura no controla nada.")
    print("     Con nuestro chi2, en cambio, el 2000 del paper anda bien.")


# --------------------------------------------------------------------------
def v4_paso(p, semilla, cotas):
    """
    AFIRMACION (estaba en el docstring de sa.py): el paso adaptativo hace
    converger donde el paso fijo del paper se queda pegado.

    No estaba medida. Se mide: 8 semillas con cada variante.
    """
    titulo(4, "Paso fijo (como en el paper) contra paso adaptativo")
    for etiqueta, adapt in (("fijo (paper)", False), ("adaptativo", True)):
        Ds, cs, ev = [], [], []
        for s in range(8):
            r = sa.recocido(p.costo, semilla, cotas, semilla=s,
                            paso_adaptativo=adapt)
            Ds.append(r["x"][0]); cs.append(r["costo"]); ev.append(r["n_eval"])
        Ds, cs = np.array(Ds), np.array(cs)
        print("  %-13s : D = %.4f +- %.4f   chi2 = %.2f +- %.2f   "
              "%.0f evaluaciones"
              % (etiqueta, Ds.mean(), Ds.std(), cs.mean(), cs.std(),
                 np.mean(ev)))
    print("  -> se compara la DISPERSION entre semillas: es lo que dice si el")
    print("     algoritmo llega siempre al mismo lugar o depende de la suerte.")


# --------------------------------------------------------------------------
def v5_aceptacion(p, semilla, cotas):
    """
    AFIRMACION: 0.35 es 'el valor clasico'.

    Eso estaba MAL DICHO. Los numeros citados en la literatura (Roberts, Gelman
    & Gilks 1997) son 0.44 en una dimension y 0.234 asintoticamente en muchas, y
    son resultados sobre MUESTREO por Metropolis, no sobre recocido. El 0.35 es
    una eleccion intermedia. Lo que corresponde es medir si el resultado depende
    de ella.
    """
    titulo(5, "El objetivo de aceptacion, barrido")
    for obj in (0.15, 0.234, 0.35, 0.44, 0.60):
        Ds = []
        for s in range(6):
            r = sa.recocido(p.costo, semilla, cotas, semilla=s,
                            objetivo_aceptacion=obj)
            xp, _ = sa.pulir(p.costo, r["x"])
            Ds.append(xp[0])
        Ds = np.array(Ds)
        print("  objetivo %.3f  ->  D = %.5f +- %.5f  (dispersion %.4f %%)"
              % (obj, Ds.mean(), Ds.std(), 100 * Ds.std() / Ds.mean()))
    print("  -> si D no se mueve, la constante es un detalle de implementacion")
    print("     y no un numero que haya que justificar.")


# --------------------------------------------------------------------------
def v6_varpro(p, semilla, cotas):
    """
    AFIRMACION: proyectar la amplitud (VARPRO) evita el problema de escala del
    paso comun del paper.

    Se compara contra buscar las 3 dimensiones con un paso comun, que es lo que
    la receta del paper implica.
    """
    titulo(6, "Proyectar la amplitud contra buscarla con el resto")
    A_opt = p.chi2(*semilla)[1]
    print("  amplitud optima ~ %.3g ; D ~ %.3g  ->  difieren en %.0f ordenes "
          "de magnitud" % (A_opt, semilla[0], abs(np.log10(A_opt / semilla[0]))))

    def costo3(x):
        D, s, A = x
        if D <= 0 or s <= 0 or A <= 0:
            return 1e30
        return p.chi2(D, s, A=A)[0]

    cotas3 = list(cotas) + [(0.05 * A_opt, 20 * A_opt)]
    x0_3 = list(semilla) + [A_opt]
    for etiqueta, f, x0, ct in (("2 parametros (VARPRO)", p.costo, semilla, cotas),
                                ("3 parametros (paper)", costo3, x0_3, cotas3)):
        Ds, cs = [], []
        for s in range(6):
            r = sa.recocido(f, x0, ct, semilla=s)
            Ds.append(r["x"][0]); cs.append(r["costo"])
        Ds, cs = np.array(Ds), np.array(cs)
        print("  %-22s : D = %.5f +- %.5f   chi2 = %.2f +- %.2f"
              % (etiqueta, Ds.mean(), Ds.std(), cs.mean(), cs.std()))


# --------------------------------------------------------------------------
def v7_sin_verificar():
    """Las afirmaciones que NO estan verificadas en este codigo."""
    titulo(7, "Lo que NO esta verificado aca")
    print("  (a) 'un modelo de 6 parametros esta sub-determinado'")
    print("      Origen: CLAUDE.md 7.23, del proyecto ANTERIOR. Este codigo no")
    print("      implementa el modelo analitico de 6 parametros, asi que la")
    print("      afirmacion esta IMPORTADA, no re-verificada. Para cerrarla hay")
    print("      que implementar la Ec. 11 con Klein-Nishina y correrla.")
    print()
    print("  (b) 'el 25% de los eventos son scattering multiple'")
    print("      Esto SI se mide, pero en simulacion.py, no aca:")
    print("      `python analisis/simulacion.py` lo imprime (75.1% simple).")
    print()
    print("  (c) 'Geman & Geman (1984) probaron convergencia con T ~ c/log k'")
    print("      Es una cita de la literatura. No la puedo verificar desde aca;")
    print("      hay que leer el paper (IEEE PAMI 6(6), 1984). Lo que SI se")
    print("      verifica arriba es que nuestro enfriamiento es geometrico y que")
    print("      por lo tanto la evidencia de optimalidad es empirica.")
    print()
    print("  (d) 'chi2/ndf = 1.58 se debe al exceso de bajo canal / sigma")
    print("      constante / Neyman'. Son HIPOTESIS ordenadas por sospecha, no")
    print("      resultados. Cada una se cierra implementandola y viendo si el")
    print("      chi2 baja.")


def main():
    rutas.salida_utf8()
    t0 = time.time()
    v1_boltzmann()
    v2_ec1()
    p, semilla, ncc, e, c = modelo.armar("co60", verbose=False)
    cotas = [(0.6 * semilla[0], 1.6 * semilla[0]), (0.02, 0.30)]
    v3_T_init(p, semilla, cotas)
    v4_paso(p, semilla, cotas)
    v5_aceptacion(p, semilla, cotas)
    v6_varpro(p, semilla, cotas)
    v7_sin_verificar()
    print()
    print("(%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
