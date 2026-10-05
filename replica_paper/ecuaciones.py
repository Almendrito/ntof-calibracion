"""
ecuaciones.py — las ecuaciones del paper, tal como estan impresas.

Cada funcion cita el numero de ecuacion. Donde la ecuacion impresa esta MAL se
implementan las dos versiones (`*_impresa` y la corregida) para poder mostrar la
diferencia con numeros en vez de afirmarla.

Validacion: `python replica_paper/ecuaciones.py` redibuja la Fig. 1 del paper.
Con A = 1, r0 = 2.818e-13 cm y m0c^2 = 511 keV la Ec. 10 reproduce los valores
del eje de esa figura (1.02e-27 en E->0 para E0 = 500 keV, 0.58e-27 para 662,
0.255e-27 para 1000, 0.064e-27 para 2000, 0.010e-27 para 5000), o sea la Ec. 10
impresa es CORRECTA y sus unidades son cm^2/keV.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"


# --------------------------------------------------------------- constantes --
M_E = 511.0          # m_e c^2 [keV]. El paper la llama "electron rest mass
                     # energy" y no da el valor; 511 es el que reproduce la
                     # Fig. 1 y los bordes rotulados en la Fig. 3.
R0 = 2.818e-13       # radio clasico del electron [cm]. El paper dice
                     # "approximately equal to 2.818e-15 m"; en cm la Ec. 10
                     # da exactamente la escala vertical de la Fig. 1.


# ------------------------------------------------------------------- Ec. 1 ---
def e_max_impresa(e_gamma, mec2=M_E):
    """
    Ec. 1 TAL COMO SE IMPRIME:  E_max = E_gamma / (1 + 2 E_gamma/(m_e c^2)).

    Esta mal: eso es la energia del FOTON retrodispersado, no la del electron.
    Se conserva para poder mostrar el error con numeros (para Cs-137 da 184.3
    keV, y el propio paper rotula 477 en su Fig. 3).
    """
    e = np.asarray(e_gamma, float)
    return e / (1.0 + 2.0 * e / mec2)


def e_max(e_gamma, mec2=M_E):
    """
    Ec. 1 corregida:  E_max = 2*alpha*E_gamma/(1 + 2*alpha), alpha = E_g/m_ec^2.

    Es la energia del electron dispersado a 180 grados. Es la que usa el propio
    paper en sus figuras aunque imprima la otra.
    """
    e = np.asarray(e_gamma, float)
    a = e / mec2
    return 2.0 * a * e / (1.0 + 2.0 * a)


# ------------------------------------------------------------------- Ec. 2 ---
def dsigma_domega(theta, e_gamma, mec2=M_E, r0=R0):
    """
    Ec. 2, Klein-Nishina diferencial en angulo solido. `theta` en radianes.

    No la usa ninguna otra parte del codigo (el paper la cita para llegar a la
    Ec. 10 via el jacobiano), pero se implementa para poder verificar que la
    Ec. 10 es su transformada.
    """
    a = e_gamma / mec2
    c = np.cos(np.asarray(theta, float))
    f = 1.0 / (1.0 + a * (1.0 - c))
    return (r0 ** 2 * f ** 2 * ((1.0 + c ** 2) / 2.0)
            * (1.0 + a ** 2 * (1.0 - c) ** 2
               / ((1.0 + c ** 2) * (1.0 + a * (1.0 - c)))))


# ------------------------------------------------------------------- Ec. 3 ---
# sigma_total = sqrt(sigma_LY^2 + sigma_gain^2 + sigma_elec^2)   (Tabla I)
#   sigma_LY   = sqrt(N_ph) * LY_0
#   sigma_gain = ENF * sqrt(N_pe)
#   sigma_elec = sqrt((4kT + 2qIR) R df)
# El paper NUNCA les asigna valores: son la justificacion fisica de que el
# ensanchamiento sea gaussiano, no un modelo que se ajuste. Se dejan escritas y
# sin evaluar a proposito.
def sigma_total(sigma_ly, sigma_gain, sigma_elec):
    """Ec. 3: suma en cuadratura de los tres terminos de la Tabla I."""
    return np.sqrt(np.asarray(sigma_ly, float) ** 2
                   + np.asarray(sigma_gain, float) ** 2
                   + np.asarray(sigma_elec, float) ** 2)


# ------------------------------------------------------------------ Ec. 10 ---
def dsigma_de(e_e, e_gamma, A=1.0, mec2=M_E, r0=R0):
    """
    Ec. 10: seccion eficaz diferencial en energia del ELECTRON de retroceso.

        dsigma/dE = A*pi*r0^2/(alpha^2 m0c^2) *
                    [ 2 + (E/(hv-E))^2 ( 1/alpha^2 - 2(hv-E)/(alpha E)
                                         + (hv-E)/hv ) ]

    Fuera de 0 < E <= E_max devuelve 0: el corte en E_max es EL borde Compton, y
    es lo que la convolucion de la Ec. 8 despues suaviza.
    """
    E = np.asarray(e_e, float)
    a = e_gamma / mec2
    emax = float(e_max(e_gamma, mec2))
    y = np.zeros_like(E)
    ok = (E > 0.0) & (E <= emax)
    if not np.any(ok):
        return y
    Ee = E[ok]
    s = e_gamma - Ee                      # energia del foton dispersado
    pref = A * np.pi * r0 ** 2 / (a ** 2 * mec2)
    y[ok] = pref * (2.0 + (Ee / s) ** 2
                    * (1.0 / a ** 2 - 2.0 * s / (a * Ee) + s / e_gamma))
    return y


# ------------------------------------------------------------------- Ec. 9 ---
def kernel_gaussiano(sigma, paso, n_sigmas=4.0, largo_senal=None):
    """
    Ec. 9: G(E) = 1/(sigma sqrt(2pi)) exp(-E^2/(2 sigma^2)), muestreada.

    Se NORMALIZA a suma 1 (no a integral 1) porque se usa con np.convolve sobre
    un histograma: asi la convolucion conserva las cuentas y A y sigma quedan
    desacoplados. El radio se capa a (largo-1)//2 porque np.convolve(mode="same")
    devuelve el largo del KERNEL si el kernel es mas largo que la senal.
    """
    r = max(1, int(np.ceil(n_sigmas * abs(sigma) / paso)))
    if largo_senal is not None:
        r = min(r, max(1, (int(largo_senal) - 1) // 2))
    x = np.arange(-r, r + 1) * paso
    k = np.exp(-x ** 2 / (2.0 * sigma ** 2))
    return k / k.sum()


# ------------------------------------------------------------------- Ec. 8 ---
def modelo(canales, e_gamma, A, B, sigma, mec2=M_E, r0=R0, pad=True):
    """
    Ec. 8: la lineshape llevada a canal y convolucionada con la gaussiana.

    El mapa canal<->energia de la Ec. 8 pasa por el ORIGEN: E_eff = E_e * B/E_max,
    o sea canal cero <-> energia cero, sin ordenada. Aca se invierte para poder
    evaluar en canal: E_e = canal * E_max / B.

    `pad` extiende la grilla hasta el canal 0 y 4 sigma por arriba ANTES de
    convolucionar. No es un agregado al modelo: la Ec. 8 es una integral sobre
    todo E_e, y evaluar la convolucion solo dentro de la ventana de ajuste es un
    error NUMERICO (np.convolve rellena con ceros, y el ajuste lo compensa
    achicando sigma).
    """
    ch = np.asarray(canales, float)
    if ch.size < 2:
        raise ValueError("hacen falta al menos 2 canales")
    paso = float(ch[1] - ch[0])
    if not np.allclose(np.diff(ch), paso, rtol=1e-6, atol=1e-9):
        raise ValueError("la grilla de canales tiene que ser uniforme")
    if not np.isfinite(B) or not np.isfinite(sigma) or B <= 0 or sigma <= 0:
        return np.zeros_like(ch)

    if pad:
        # El tope existe porque un optimizador sin cotas puede pedir sigma
        # enorme y la grilla extendida se comeria toda la memoria. 20 veces el
        # largo pedido es mucho mas de lo que cualquier caso sano necesita.
        tope = 10 * ch.size + 100
        n_izq = min(tope,
                    int(np.ceil(max(ch[0], 0.0) / paso))
                    + int(np.ceil(4 * abs(sigma) / paso)))
        n_der = min(tope, int(np.ceil(4 * abs(sigma) / paso)))
        ext = np.arange(-n_izq, ch.size + n_der) * paso + ch[0]
    else:
        ext = ch

    emax = float(e_max(e_gamma, mec2))
    y = dsigma_de(ext * emax / B, e_gamma, A=A, mec2=mec2, r0=r0)
    k = kernel_gaussiano(sigma, paso, largo_senal=ext.size)
    y = np.convolve(y, k, mode="same")
    return y if ext is ch else np.interp(ch, ext, y)


# --------------------------------------------------------------- Ec. 11-13 ---
def delta_e_impresa(eg1, eg2, mec2=M_E):
    """
    Ec. 12 TAL COMO SE IMPRIME: "dE = E2 - E1 energy difference between gamma
    lines". Para Co-60 da 159.3 keV.
    """
    return float(eg2 - eg1)


def delta_e_bordes(eg1, eg2, mec2=M_E):
    """
    Ec. 12 corregida: D convierte energia del ELECTRON a canal, asi que lo que
    se separa en canal son los BORDES. Para Co-60 da 154.7 keV.
    """
    return float(e_max(eg2, mec2) - e_max(eg1, mec2))


def _delta_e(eg1, eg2, modo, mec2=M_E):
    if modo == "gammas":
        return delta_e_impresa(eg1, eg2, mec2)
    if modo == "bordes":
        return delta_e_bordes(eg1, eg2, mec2)
    raise ValueError("delta_e desconocido: %r (usar 'gammas' o 'bordes')" % modo)


def modelo_dos_bordes(canales, eg1, eg2, A1, A2, B1, sigma1, sigma2, D,
                      delta_e="bordes", mec2=M_E, r0=R0, pad=True):
    """
    Ec. 13: S_comp = f1(ch; A1, B1, sigma1) + f2(ch; A2, CE_ch2, sigma2),
    con la Ec. 12 CE_ch2 = B1 + D*dE.  Theta = {A1,A2,B1,sigma1,sigma2,D} (Ec.11).

    `delta_e` elige entre la Ec. 12 impresa ("gammas") y la corregida ("bordes").
    """
    c1, c2, _ = componentes_dos_bordes(canales, eg1, eg2, A1, A2, B1, sigma1,
                                       sigma2, D, delta_e, mec2, r0, pad)
    return c1 + c2


def componentes_dos_bordes(canales, eg1, eg2, A1, A2, B1, sigma1, sigma2, D,
                           delta_e="bordes", mec2=M_E, r0=R0, pad=True):
    """Las dos componentes por separado y el canal del segundo borde."""
    B2 = B1 + D * _delta_e(eg1, eg2, delta_e, mec2)
    return (modelo(canales, eg1, A1, B1, sigma1, mec2, r0, pad),
            modelo(canales, eg2, A2, B2, sigma2, mec2, r0, pad),
            B2)


# ------------------------------------------------------------------ Ec. 14 ---
def loss_minimos_cuadrados(s1, s_modelo):
    """Ec. 14: loss = sum_i [S1(E_i) - S_model(E_i; p)]^2."""
    r = np.asarray(s1, float) - np.asarray(s_modelo, float)
    return float(np.sum(r * r))


# ------------------------------------------------------------------ Ec. 15 ---
def pearson(s1, s2):
    """Ec. 15: coeficiente de correlacion de Pearson entre S1 y S2."""
    a = np.asarray(s1, float)
    b = np.asarray(s2, float)
    a = a - a.mean()
    b = b - b.mean()
    na = np.sqrt(np.sum(a * a))
    nb = np.sqrt(np.sum(b * b))
    if na <= 0 or nb <= 0:
        return -1.0
    return float(np.sum(a * b) / (na * nb))


# ------------------------------------------------------------------- demo ----
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from replica_paper import rutas
    rutas.salida_utf8()

    print("=" * 74)
    print("Ec. 1: borde Compton")
    print("=" * 74)
    print("  linea gamma     Ec.1 impresa    Ec.1 corregida   rotulo del paper")
    for eg, rot in ((511.0, "(no la usan)"), (661.657, "477"), (834.848, "639"),
                    (1173.228, "963"), (1274.537, "1061"), (1332.492, "1117")):
        print("  %10.3f       %8.1f        %8.1f         %s"
              % (eg, e_max_impresa(eg), e_max(eg), rot))
    print("  -> la Ec. 1 impresa da la energia del foton retrodispersado.")
    print("     Los rotulos de la Fig. 3 del paper siguen la CORREGIDA.")

    print()
    print("=" * 74)
    print("Ec. 12: separacion de los dos bordes del Co-60")
    print("=" * 74)
    print("  impresa (lineas gamma) : %.1f keV" % delta_e_impresa(1173.228, 1332.492))
    print("  corregida (bordes)     : %.1f keV" % delta_e_bordes(1173.228, 1332.492))

    print()
    print("=" * 74)
    print("Ec. 10: reproduccion de la Fig. 1 del paper")
    print("=" * 74)
    print("  E0 [keV]   dsigma/dE en E->0    pico en el borde    E_max [keV]")
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for e0 in (500.0, 662.0, 1000.0, 2000.0, 5000.0):
        E = np.linspace(0.0, 5000.0, 40001)
        y = dsigma_de(E, e0)
        y0 = float(dsigma_de(np.array([1e-9]), e0)[0])
        print("  %7.0f     %.4e           %.4e          %7.1f"
              % (e0, y0, y.max(), float(e_max(e0))))
        ax.plot(E, y, lw=1.4, label="$E_0$ = %gkeV" % e0)
    ax.set_xlabel("E/keV")
    ax.set_ylabel(r"d$\sigma$/dE")
    ax.set_xlim(0, 5000)
    ax.grid(alpha=0.4)
    ax.legend()
    fig.tight_layout()
    salida = rutas.figura("fig01_ec10.png")
    fig.savefig(salida, dpi=130)
    print()
    print("  figura: %s" % salida)
    print("  Comparar con la Fig. 1 del paper: los valores en E->0 son")
    print("  1.02 / 0.58 / 0.255 / 0.064 / 0.010 (x1e-27) y los picos")
    print("  1.68 / 1.12 / 0.645 / 0.283 / 0.11.")


if __name__ == "__main__":
    _demo()
