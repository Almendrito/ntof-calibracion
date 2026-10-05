"""
Etapa 4 — el borde Compton y la forma teorica del continuo (paper, Ec. 1, 2, 10).

    python desde_cero/fisica.py

Que resuelve esta etapa
-----------------------
Nada del detector: solo que espectro de energia le deja a un ELECTRON LIBRE Y
QUIETO un foton de energia E_gamma que le hace scattering Compton. Es la
"plantilla" del paper: la forma que despues se lleva al eje medido (Etapa 5),
se compara contra el dato (6) y se ajusta (7).

Por que alcanza con electrones libres en un centellador plastico: el plastico
es carbono e hidrogeno (Z = 6 y 1). La energia de ligadura de sus electrones
es de decenas de eV, contra cientos de keV del foton, y el efecto fotoelectrico
(que va como Z^4-5) es 1e-4 veces el Compton a 662 keV. Por eso no hay
fotopico y el borde Compton es lo unico que se puede usar para calibrar.

Lo que hay aca
--------------
* `e_max`         Ec. 1 CORREGIDA: energia maxima del electron (el borde).
* `e_max_impresa` Ec. 1 como esta IMPRESA en el paper: da otra cosa.
* `dsigma_domega` Ec. 2: Klein-Nishina en angulo solido.
* `dsigma_de`     Ec. 10: la misma seccion eficaz en energia del electron.
* `dsigma_de_jacobiano`  Ec. 10 obtenida de la Ec. 2 a mano, via el
                  jacobiano. Existe para comprobar la Ec. 10.
* `sigma_total_*` la seccion eficaz total, integrando cada forma por su
                  variable y con la formula cerrada. Las tres tienen que dar
                  lo mismo.
* `FUENTES`       las lineas gamma de Cs-137 y Co-60.
"""

import numpy as np

# ------------------------------------------------------------- constantes ----
# CODATA 2018. El paper no da m_e c^2 y usa r0 ~ 2.818e-15 m; replica_paper
# usa 511 y 2.818e-13 porque con esos redibuja la Fig. 1 al digito. Para los
# bordes da igual (477.34 con los dos); para compararse con replica_paper se
# le pasan estas constantes.
M_E = 510.99895          # m_e c^2 [keV]
R0 = 2.8179403262e-13    # radio clasico del electron [cm]
BARN = 1e-24             # [cm^2]

# Lineas gamma: (energia [keV], fotones por desintegracion). Valores de
# ENSDF (NNDC). Del Co-60 importa que las dos lineas salen casi siempre
# JUNTAS y con la misma intensidad: los dos bordes pesan lo mismo, y eso es lo
# que la plantilla de dos bordes tiene que respetar.
FUENTES = {
    "Cs-137": ((661.657, 0.8510),),
    "Co-60": ((1173.228, 0.9985), (1332.492, 0.999826)),
}


# --------------------------------------------------------------------------
# Ec. 1 — el borde
# --------------------------------------------------------------------------
def e_max(e_gamma, mec2=M_E):
    """
    Energia maxima que el foton le puede dar al electron [keV].

    Cinematica (conservar energia y momento en el choque foton-electron):
    el foton sale con

        E' = E_gamma / (1 + alpha (1 - cos theta)),   alpha = E_gamma / m_e c^2

    y el electron se lleva T = E_gamma - E'. T es maxima cuando E' es minima,
    o sea con el foton rebotando hacia atras (theta = 180 grados, cos = -1):

        E'_min = E_gamma / (1 + 2 alpha)
        T_max  = E_gamma - E'_min = E_gamma * 2 alpha / (1 + 2 alpha)

    Esa T_max es el borde Compton.
    """
    e = np.asarray(e_gamma, float)
    a = e / mec2
    return 2.0 * a * e / (1.0 + 2.0 * a)


def e_max_impresa(e_gamma, mec2=M_E):
    """
    La Ec. 1 como esta IMPRESA: E_gamma / (1 + 2 alpha). Es E'_min, la energia
    del FOTON retrodispersado, no la del electron (le falta el 2 alpha del
    numerador). Para Cs-137 da 184.3 keV; la propia Fig. 3 del paper rotula el
    borde en 477. Se conserva para mostrar la errata con numeros.
    """
    e = np.asarray(e_gamma, float)
    return e / (1.0 + 2.0 * e / mec2)


# --------------------------------------------------------------------------
# Ec. 2 — Klein-Nishina en angulo
# --------------------------------------------------------------------------
def dsigma_domega(theta, e_gamma, mec2=M_E, r0=R0):
    """
    Ec. 2: dsigma/dOmega [cm^2/sr] para un foton de `e_gamma` que sale a
    `theta` [rad] de su direccion original, sobre un electron libre en reposo.

    Con P = E'/E_gamma = 1 / (1 + alpha (1 - cos theta)):

        dsigma/dOmega = r0^2/2 * P^2 * (P + 1/P - sin^2 theta)

    que es la forma de la Ec. 2 del paper reescrita. Limite de baja energia
    (alpha -> 0, P -> 1): r0^2 (1 + cos^2 theta)/2, el scattering de Thomson
    de una onda electromagnetica clasica sobre una carga libre.
    """
    c = np.cos(np.asarray(theta, float))
    P = 1.0 / (1.0 + (e_gamma / mec2) * (1.0 - c))
    return 0.5 * r0 ** 2 * P ** 2 * (P + 1.0 / P - (1.0 - c ** 2))


# --------------------------------------------------------------------------
# Ec. 10 — Klein-Nishina en energia del electron
# --------------------------------------------------------------------------
def dsigma_de(e_e, e_gamma, mec2=M_E, r0=R0):
    """
    Ec. 10: dsigma/dT [cm^2/keV], T = energia cinetica del electron.

        dsigma/dT = pi r0^2 / (alpha^2 m_e c^2) *
                    [ 2 + (T/(E-T))^2 ( 1/alpha^2 + (E-T)/E - 2(E-T)/(alpha T) ) ]

    con E = E_gamma. Cero fuera de 0 < T <= e_max: por encima del borde no hay
    angulo que lleve ahi, y ese corte es el escalon que el detector despues
    borronea.

    DE DONDE SALE (el jacobiano). La Ec. 2 reparte la seccion eficaz en angulo;
    el detector mide energia. Cada angulo theta da UNA energia T(theta), asi
    que la probabilidad que cae en dOmega es la que cae en el dT
    correspondiente:

        dsigma/dT = dsigma/dOmega * |dOmega/dT|

    Por simetria alrededor del haz, dOmega = 2 pi d(cos theta). Y de la
    cinematica de `e_max`, con c = cos theta:

        E' = E / (1 + alpha (1 - c))
        dE'/dc = E alpha / (1 + alpha (1 - c))^2 = alpha E'^2 / E = E'^2 / m_e c^2
        T = E - E'   =>   |dT/dc| = |dE'/dc| = E'^2 / m_e c^2

    Entonces

        dsigma/dT = 2 pi (m_e c^2 / E'^2) * dsigma/dOmega(theta(T))

    con E' = E - T y cos theta = 1 - m_e c^2 (1/E' - 1/E). Sustituyendo la
    Ec. 2 y simplificando sale la expresion de arriba. `dsigma_de_jacobiano`
    hace este camino numericamente, sin simplificar, para comprobarlo.

    Que forma tiene: casi plana a T baja (el "plateau Compton"), sube hacia el
    borde porque cerca de 180 grados muchos angulos dan casi la misma T (el
    jacobiano crece), y cae a cero de golpe en e_max.
    """
    T = np.asarray(e_e, float)
    a = e_gamma / mec2
    y = np.zeros_like(T)
    ok = (T > 0.0) & (T <= e_max(e_gamma, mec2))
    t = T[ok]
    s = e_gamma - t                         # energia del foton dispersado
    y[ok] = (np.pi * r0 ** 2 / (a ** 2 * mec2)
             * (2.0 + (t / s) ** 2 * (1.0 / a ** 2 + s / e_gamma
                                      - 2.0 * s / (a * t))))
    return y


def dsigma_de_jacobiano(e_e, e_gamma, mec2=M_E, r0=R0):
    """La Ec. 10 obtenida de la Ec. 2 por el jacobiano, sin simplificar."""
    T = np.asarray(e_e, float)
    y = np.zeros_like(T)
    ok = (T > 0.0) & (T <= e_max(e_gamma, mec2))
    ep = e_gamma - T[ok]
    c = 1.0 - mec2 * (1.0 / ep - 1.0 / e_gamma)
    theta = np.arccos(np.clip(c, -1.0, 1.0))
    y[ok] = 2.0 * np.pi * (mec2 / ep ** 2) * dsigma_domega(theta, e_gamma, mec2, r0)
    return y


# --------------------------------------------------------------------------
# Seccion eficaz total, tres caminos
# --------------------------------------------------------------------------
def sigma_total_angulo(e_gamma, mec2=M_E, r0=R0):
    """Integral de la Ec. 2 sobre la esfera: 2 pi * int dsigma/dOmega d(cos)."""
    from scipy.integrate import quad
    f = lambda c: dsigma_domega(np.arccos(c), e_gamma, mec2, r0)   # noqa: E731
    return 2.0 * np.pi * quad(f, -1.0, 1.0, epsabs=0, epsrel=1e-12)[0]


def sigma_total_energia(e_gamma, mec2=M_E, r0=R0):
    """Integral de la Ec. 10 en energia del electron, de 0 a e_max."""
    from scipy.integrate import quad
    f = lambda t: dsigma_de(np.array([t]), e_gamma, mec2, r0)[0]  # noqa: E731
    return quad(f, 0.0, float(e_max(e_gamma, mec2)), epsabs=0, epsrel=1e-12,
                limit=200)[0]


def sigma_total_cerrada(e_gamma, mec2=M_E, r0=R0):
    """
    Formula cerrada de Klein-Nishina (Evans, The Atomic Nucleus, cap. 23):

      sigma = 2 pi r0^2 { (1+a)/a^2 [ 2(1+a)/(1+2a) - ln(1+2a)/a ]
                          + ln(1+2a)/(2a) - (1+3a)/(1+2a)^2 }

    Para a -> 0 tiende a la de Thomson, 8 pi r0^2 / 3 = 0.6652 b, pero
    evaluada tal cual a `a` chico pierde precision por cancelacion (restas de
    numeros casi iguales, divididas por a^2): medido contra la integral en
    angulo, el error relativo es 8e-8 a 10 eV, 4e-5 a 1 eV y 30 % a 0.01 eV.
    A las energias de este proyecto (cientos de keV) es exacta a 1e-11.
    """
    a = e_gamma / mec2
    L = np.log1p(2.0 * a)
    return 2.0 * np.pi * r0 ** 2 * ((1 + a) / a ** 2 * (2 * (1 + a) / (1 + 2 * a) - L / a)
                                    + L / (2 * a) - (1 + 3 * a) / (1 + 2 * a) ** 2)


def sigma_thomson(r0=R0):
    """8 pi r0^2 / 3: el limite clasico, sin retroceso del electron."""
    return 8.0 * np.pi * r0 ** 2 / 3.0


# --------------------------------------------------------------------------
# Demo: la Fig. 1 del paper y nuestras tres lineas
# --------------------------------------------------------------------------
def _demo():
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("  Ec. 1: el borde Compton")
    print("  %-8s %10s   %12s   %12s" % ("fuente", "E_gamma", "E_max (bien)",
                                         "Ec.1 impresa"))
    for f, lineas in FUENTES.items():
        for eg, _ in lineas:
            print("  %-8s %10.3f   %12.2f   %12.2f"
                  % (f, eg, e_max(eg), e_max_impresa(eg)))
    eg1, eg2 = (e for e, _ in FUENTES["Co-60"])
    print("  Co-60: separacion entre BORDES %.1f keV, entre LINEAS %.1f keV"
          % (e_max(eg2) - e_max(eg1), eg2 - eg1))

    print()
    print("  seccion eficaz total [b]: por angulo / por energia / cerrada")
    for eg in (500.0, 661.657, 1000.0, 1173.228, 1332.492, 2000.0):
        print("  %9.3f keV   %.6f   %.6f   %.6f"
              % (eg, sigma_total_angulo(eg) / BARN, sigma_total_energia(eg) / BARN,
                 sigma_total_cerrada(eg) / BARN))
    print("  Thomson: %.6f b" % (sigma_thomson() / BARN))

    fig, ejes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = ejes[0]
    E = np.linspace(0.0, 5000.0, 40001)
    for e0 in (500.0, 662.0, 1000.0, 2000.0, 5000.0):
        # con las constantes del paper, para comparar al digito con su Fig. 1
        ax.plot(E, dsigma_de(E, e0, mec2=511.0, r0=2.818e-13),
                lw=1.3, label="$E_0$ = %g keV" % e0)
    ax.set(xlabel="E [keV]", ylabel=r"d$\sigma$/dE [cm$^2$/keV]", xlim=(0, 5000),
           title="Fig. 1 del paper (Ec. 10, m$_e$c$^2$ = 511, r$_0$ = 2.818e-13)")
    ax = ejes[1]
    for f, lineas in FUENTES.items():
        for eg, _ in lineas:
            T = np.linspace(0.0, 1.1 * e_max(eg), 4000)
            ax.plot(T, dsigma_de(T, eg) / BARN * 1e3, lw=1.3,
                    label="%s  %.1f keV  (borde %.1f)" % (f, eg, e_max(eg)))
            ax.axvline(e_max(eg), color="0.6", lw=0.6, ls=":")
    ax.set(xlabel="energia del electron T [keV]",
           ylabel=r"d$\sigma$/dT [mb/keV]",
           title="las lineas de este proyecto")
    for ax in ejes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa4_klein_nishina.png")
    fig.savefig(salida, dpi=120)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
