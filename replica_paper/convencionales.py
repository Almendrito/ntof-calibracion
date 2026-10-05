"""
convencionales.py — los dos metodos contra los que el paper compara.

Son la vara de la Fig. 4. El paper los describe en la Introduccion:

  1. "Gaussian FWHM Method": "The first method employs Gaussian function fitting
     to the spectral descent slope, defining the channel position at 50% of the
     maximum amplitude as the Compton edge". (El resumen dice "50% 70%", asi que
     la fraccion queda como parametro.)

  2. "Difference Method": "The second method calculates the first derivative of
     the smoothed measured energy spectrum, identifying the channel
     corresponding to the derivative minimum as the Compton edge".

Lo que el paper NO especifica, y aca queda explicito
----------------------------------------------------
- Con que se suaviza y con que ventana (metodo 2). Se usa Savitzky-Golay, que es
  lo estandar y lo que muestra la Fig. 5 del paper como "Smoothed Data"; la
  ventana es un parametro y `sensibilidad_ventana` mide cuanto cambia el
  resultado con ella.
- Que region se ajusta. Los dos metodos necesitan una, y el paper no da regla.
  En el pipeline se les pasa la MISMA que usa el NCC-SA, que es lo unico que
  hace la comparacion justa; `region_flanco` da una alternativa automatica que
  no depende del NCC.
- Que es "the maximum amplitude" del metodo 1: aca es la amplitud de la
  gaussiana AJUSTADA, que es la unica cantidad bien definida despues de un
  ajuste gaussiano.
"""

import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import ecuaciones as ec

VENTANA = 21            # ventana de Savitzky-Golay [NO esta en el paper]
ORDEN = 3               # orden del polinomio       [NO esta en el paper]
FRACCION = 0.5          # "50% of the maximum amplitude"


def suavizar(y, ventana=VENTANA, orden=ORDEN):
    """Savitzky-Golay, con la ventana capada al largo de la senal."""
    y = np.asarray(y, float)
    v = int(ventana)
    if v % 2 == 0:
        v += 1
    v = min(v, y.size - (1 - y.size % 2))
    if v <= orden + 1:
        return y.copy()
    return savgol_filter(y, v, orden)


def region_flanco(canales, cuentas, ch_min=None, ventana=VENTANA,
                  frac_baja=0.05):
    """
    Una region alrededor del flanco descendente, sin usar el NCC.

    Toma el maximo del espectro suavizado por encima de `ch_min`, baja hacia la
    derecha hasta `frac_baja` de ese maximo, y arma una region simetrica
    alrededor del pico. No esta en el paper: existe para poder correr los
    metodos convencionales sin que hereden la ventana del metodo que se les
    quiere comparar.
    """
    ch = np.asarray(canales, float)
    y = suavizar(cuentas, ventana)
    m = np.ones(ch.size, bool) if ch_min is None else (ch >= float(ch_min))
    if m.sum() < 5:
        raise ValueError("ch_min deja menos de 5 bins")
    idx = np.nonzero(m)[0]
    i_pico = idx[int(np.argmax(y[idx]))]
    umbral = frac_baja * y[i_pico]
    i_fin = i_pico
    while i_fin + 1 < ch.size and y[i_fin] > umbral:
        i_fin += 1
    ancho = max(i_fin - i_pico, 3)
    i_ini = max(idx[0], i_pico - ancho)
    return float(ch[i_ini]), float(ch[i_fin])


def _mascara(canales, region):
    ch = np.asarray(canales, float)
    if region is None:
        return np.ones(ch.size, bool)
    lo, hi = region
    return (ch >= float(lo)) & (ch <= float(hi))


def _gaussiana(x, A, mu, s):
    return A * np.exp(-(x - mu) ** 2 / (2.0 * s ** 2))


def metodo_gauss(canales, cuentas, region, fraccion=FRACCION, p0=None):
    """
    Metodo 1 del paper: gaussiana ajustada al flanco, borde al `fraccion` de la
    amplitud maxima.

    El canal buscado es el del lado DESCENDENTE:
        borde = mu + s*sqrt(2 ln(1/fraccion))
    que para fraccion = 0.5 es mu + 1.17741*s (el semiancho a media altura, de
    ahi el nombre "Gaussian FWHM Method").
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    m = _mascara(ch, region)
    x, yy = ch[m], y[m]
    if x.size < 6:
        raise ValueError("region demasiado corta para ajustar 3 parametros")

    if p0 is None:
        i = int(np.argmax(yy))
        p0 = [max(yy[i], 1.0), x[i], max((x[-1] - x[0]) / 4.0, 1.0)]

    def res(p):
        return _gaussiana(x, p[0], p[1], abs(p[2])) - yy

    sol = least_squares(res, p0, method="lm", max_nfev=20000)
    A, mu, s = float(sol.x[0]), float(sol.x[1]), abs(float(sol.x[2]))
    borde = mu + s * np.sqrt(2.0 * np.log(1.0 / float(fraccion)))
    # Si la region no contiene un maximo (ventanas muy angostas sobre el flanco)
    # la gaussiana ajustada se va al infinito y el "borde" sale fuera de la
    # region, o negativo. Eso no es un borde: es un ajuste que no convergio.
    if not np.isfinite(borde) or A <= 0 or not (x[0] <= borde <= x[-1]):
        borde = float("nan")
    return dict(borde=float(borde), A=A, mu=mu, sigma=s,
                fraccion=float(fraccion), region=(float(x[0]), float(x[-1])),
                exito=bool(sol.success))


def metodo_diferencia(canales, cuentas, region, ventana=VENTANA, orden=ORDEN):
    """
    Metodo 2 del paper: minimo de la primera derivada del espectro suavizado.
    """
    ch = np.asarray(canales, float)
    y = np.asarray(cuentas, float)
    m = _mascara(ch, region)
    x, yy = ch[m], y[m]
    if x.size < max(int(ventana), orden + 2):
        raise ValueError("region demasiado corta para la ventana de suavizado")
    s = suavizar(yy, ventana, orden)
    d = np.gradient(s, x)
    i = int(np.argmin(d))
    return dict(borde=float(x[i]), derivada_min=float(d[i]), suave=s,
                derivada=d, canales=x, ventana=int(ventana),
                region=(float(x[0]), float(x[-1])))


def sensibilidad_ventana(canales, cuentas, region, ventanas=(9, 15, 21, 31, 41)):
    """
    Cuanto mueve el borde la ventana de suavizado, que el paper no especifica.
    [Fuera del paper: es una medida de una ambiguedad del paper, no un metodo.]
    """
    out = []
    for v in ventanas:
        try:
            out.append((int(v), metodo_diferencia(canales, cuentas, region,
                                                  ventana=v)["borde"]))
        except ValueError:
            out.append((int(v), float("nan")))
    return out


def _demo():
    """
    El sesgo de los dos metodos sobre un borde de posicion CONOCIDA.

    Es el punto que el paper hace cualitativamente ("the assumed ideal Gaussian
    distribution significantly deviates from the observed asymmetric
    broadening, leading to systematic errors in edge localization") y que aca
    queda con numero.
    """
    from replica_paper import rutas
    rutas.salida_utf8()
    rng = np.random.default_rng(5)

    eg = 661.657
    B_real = 600.0
    print("=" * 74)
    print("convencionales.py - sesgo con verdad conocida (B = %.0f)" % B_real)
    print("=" * 74)
    print("  sigma/B   gauss 50%   gauss 70%   diferencia   (sin ruido: gauss50/dif)")
    ch = np.linspace(2.5, 1400.0, 280)
    for sr in (0.03, 0.06, 0.09, 0.12, 0.18):
        limpio = ec.modelo(ch, eg, 1.0, B_real, sr * B_real)
        limpio = limpio / limpio.max() * 800.0
        datos = rng.poisson(np.maximum(limpio, 0.0)).astype(float)
        reg = region_flanco(ch, datos)
        g50 = metodo_gauss(ch, datos, reg, 0.5)["borde"]
        g70 = metodo_gauss(ch, datos, reg, 0.7)["borde"]
        dif = metodo_diferencia(ch, datos, reg)["borde"]
        reg0 = region_flanco(ch, limpio)
        g50_0 = metodo_gauss(ch, limpio, reg0, 0.5)["borde"]
        dif_0 = metodo_diferencia(ch, limpio, reg0)["borde"]
        p = lambda v: 100.0 * (v / B_real - 1.0)
        print("   %.2f     %+6.1f%%     %+6.1f%%     %+6.1f%%        %+6.1f%% / %+6.1f%%"
              % (sr, p(g50), p(g70), p(dif), p(g50_0), p(dif_0)))
    print()
    print("  La ultima columna es el MISMO metodo sobre el espectro sin ruido:")
    print("  el sesgo no es estadistico, es del metodo. La lineshape de")
    print("  Klein-Nishina no es un escalon (tiene acumulacion antes del borde)")
    print("  y ninguno de los dos la contiene.")


if __name__ == "__main__":
    _demo()
