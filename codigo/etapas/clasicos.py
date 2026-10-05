"""
clasicos.py - los dos metodos convencionales contra los que compara el paper.

Por que estan aca
-----------------
El paper (Seccion II-C y III) no se limita a proponer NCC-SA: mide cuanto mejor
es que dos metodos clasicos, y ese es la mitad de su resultado. Los dos son:

  1. METODO DE LA DIFERENCIA ("differential minimum detection"): el borde es el
     canal donde la derivada del espectro es MAS NEGATIVA, o sea el punto de
     maxima caida del flanco.
  2. METODO GAUSSIANO FWHM: se ajusta el flanco con un escalon borroneado por
     una gaussiana y se toma el punto de media altura del AJUSTE. De paso el
     mismo ajuste entrega el FWHM, que es lo que le da el nombre.

Hasta ahora el proyecto solo tenia `model.punto_mitad_altura`, que lee el 50%
de altura DIRECTO sobre los datos. Es el mas ingenuo de los tres y se conserva
como referencia, pero no es ninguno de los dos que el paper usa.

Que relacion tienen entre si (esto es lo interesante)
----------------------------------------------------
Si el borde fuera un escalon PERFECTO borroneado por una gaussiana, los tres
metodos darian exactamente lo mismo: la derivada de ese escalon es una gaussiana
centrada justo en el escalon, y su punto de media altura tambien cae ahi. Que en
datos reales NO coincidan mide cuanto se aparta el borde Compton de un escalon:
la lineshape de Klein-Nishina tiene una acumulacion creciente ANTES del borde
(el pico de retrodispersion), y esa pendiente de fondo corre el minimo de la
derivada hacia arriba y el punto de media altura hacia abajo.

Por eso los tres estan implementados y se reportan juntos: la dispersion entre
ellos es una medida honesta del sesgo del metodo clasico, y es exactamente lo
que el paper grafica en su Fig. 4.

Ninguno de los tres es el resultado del proyecto: el numero bueno sale del
ajuste con plantilla (`ajuste/modelo_ext.py`). Estos son la vara de comparacion.

Corre solo:  python codigo/etapas/clasicos.py
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/clasicos.py` pone en sys.path la carpeta DEL ARCHIVO, no
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
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.special import erfc


# ----------------------------------------------------------------------
# utilidades comunes
# ----------------------------------------------------------------------
def _suavizar(cuentas, ventana, orden):
    """Savitzky-Golay, con la ventana saneada (impar y mayor que el orden)."""
    ventana = int(ventana)
    ventana = ventana if ventana % 2 == 1 else ventana + 1
    ventana = max(ventana, orden + 2 + (orden % 2))
    ventana = min(ventana, len(cuentas) - 1 | 1)
    return savgol_filter(np.asarray(cuentas, float), ventana, orden)


def _mascara(canales, region):
    if region is None:
        return np.ones(len(canales), bool)
    return (canales >= region[0]) & (canales <= region[1])


def _ventana_auto(canales, cuentas, region, orden=3):
    """
    Ventana de suavizado deducida del ANCHO del propio flanco.

    Por que no una constante: la derivada del borde es un pozo cuyo ancho es
    proporcional a la resolucion. Con una ventana fija, un borde ancho queda
    sub-suavizado y el pozo se vuelve chato; ahi el minimo deja de estar
    definido y lo elige el ruido. Medido sobre sintetico con B = 600 conocido y
    ventana fija de 21 bins: a sigma/B = 0.09 aparecen TRES minimos locales a
    -29.1, -28.2 y -28.1 repartidos entre los canales 580 y 619, y el metodo
    salta de +1.5% a -3.3% de sesgo entre una resolucion y la siguiente. Sin
    ruido, en cambio, el minimo es unico y se mueve suave: la inestabilidad es
    del suavizado, no del metodo.

    Como se mide el ancho: por la caida del 80% al 20% del flanco, sobre el
    espectro SUAVIZADO. Para un escalon borroneado esa distancia vale 1.68*sigma,
    asi que la ventana termina siendo ~1.26*sigma.

    El factor 0.75 no es a ojo: se eligio barriendolo contra el espectro SIN
    ruido, que da la respuesta que el metodo daria sin necesidad de suavizar.
    Comparando el sesgo con ruido contra ese patron (sigma/B = 0.03 a 0.18):

        f = 0.25  ->  +0.2 +0.9 +1.0 -2.1 -2.6 %   (sub-suaviza, se desestabiliza)
        f = 0.50  ->  +0.7 +1.9 +1.2 +5.6 +4.1 %   (erratico)
        f = 0.75  ->  +0.7 +1.5 +2.7 +4.3 +7.1 %   <- sigue al patron
        f = 1.00  ->  +0.4 +1.6 +3.3 +1.6 +5.8 %
        sin ruido ->  +0.4 +1.4 +2.6 +3.8 +6.2 %   (el metodo puro)

    O sea el criterio es reproducir el metodo, no minimizar el error: bajar el
    sesgo con un suavizado que casualmente cancela seria hacer trampa en la
    comparacion.

    Ojo con la alternativa que NO funciona: medir el ancho a media profundidad
    del pozo de la derivada. Con ruido, a resolucion ancha el minimo cae en una
    espiga y el ancho sale de 3 bins, o sea la ventana se ACHICA justo cuando
    hacia falta agrandarla (medido: 7 bins a sigma/B = 0.18 contra 25 a 0.06).
    El flanco del espectro es mucho mas estable que su derivada, que es
    justamente el motivo por el que hay que suavizar.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    m = _mascara(canales, region)
    idx = np.flatnonzero(m)
    n_reg = max(idx.size, 8)
    v_max = max(int(0.5 * n_reg) | 1, orden + 3)

    v0 = max(int(0.05 * n_reg) | 1, orden + 3)
    suave = _suavizar(cuentas, v0, orden)

    i_top = int(idx[np.argmax(suave[idx])])
    cola = idx[idx > i_top]
    if cola.size < 5:
        return int(np.clip(v0, orden + 3, v_max))
    piso = float(np.percentile(suave[cola], 10))
    alto = float(suave[i_top])
    if alto <= piso:
        return int(np.clip(v0, orden + 3, v_max))

    n80 = piso + 0.80 * (alto - piso)
    n20 = piso + 0.20 * (alto - piso)
    bajo80 = cola[suave[cola] <= n80]
    bajo20 = cola[suave[cola] <= n20]
    if bajo80.size == 0 or bajo20.size == 0:
        return int(np.clip(v0, orden + 3, v_max))
    ancho = max(int(bajo20[0]) - int(bajo80[0]), 3)   # ~1.68 sigma en bins

    return int(np.clip(int(0.75 * ancho) | 1, orden + 3, v_max))


def _vertice_parabola(y0, y1, y2):
    """
    Desplazamiento sub-bin del extremo, ajustando una parabola a 3 puntos.

    Sin esto la resolucion del metodo de la diferencia esta limitada al ancho
    del bin (5 mV*ns con el binning de este proyecto, o sea 0.4% del borde),
    que seria un piso artificial en la comparacion contra un ajuste que si
    interpola.
    """
    den = y0 - 2.0 * y1 + y2
    if abs(den) < 1e-12:
        return 0.0
    return float(np.clip(0.5 * (y0 - y2) / den, -1.0, 1.0))


# ----------------------------------------------------------------------
# 1) metodo de la diferencia
# ----------------------------------------------------------------------
def metodo_diferencia(canales, cuentas, region=None, ventana="auto", orden=3,
                      refinar=True):
    """
    Borde = canal donde la derivada del espectro suavizado es mas negativa.

    canales : eje uniforme [mV*ns]
    cuentas : espectro neto (sin suavizar; se suaviza aca)
    region  : (lo, hi) donde buscar. Sin esto el minimo se lo lleva la subida
              del exceso de bajo canal, que es mucho mas empinada que el borde.
    ventana : bins del Savitzky-Golay, o "auto" para deducirla del ancho del
              flanco (`_ventana_auto`). La derivada de datos crudos es ruido
              puro, asi que suavizar no es cosmetico: es parte del metodo, y
              con ventana fija el metodo se vuelve inestable cuando el borde es
              ancho. Un entero reproduce el comportamiento ingenuo.
    orden   : del polinomio de Savitzky-Golay.
    refinar : interpola el minimo con una parabola de 3 puntos (sub-bin).

    Devuelve dict con `borde`, la derivada y el indice del minimo, para dibujar.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    if ventana == "auto":
        ventana = _ventana_auto(canales, cuentas, region, orden)
    suave = _suavizar(cuentas, ventana, orden)
    paso = float(np.mean(np.diff(canales)))
    deriv = np.gradient(suave, paso)

    m = _mascara(canales, region)
    if m.sum() < 5:
        raise ValueError("region demasiado angosta para el metodo de la diferencia")
    idx = np.flatnonzero(m)
    i = int(idx[np.argmin(deriv[idx])])

    borde = canales[i]
    if refinar and 0 < i < len(canales) - 1:
        borde = canales[i] + paso * _vertice_parabola(deriv[i - 1], deriv[i],
                                                      deriv[i + 1])
    return {"borde": float(borde), "suave": suave, "derivada": deriv,
            "i_min": i, "pendiente": float(deriv[i]), "region": region,
            "ventana": int(ventana), "metodo": "diferencia"}


# ----------------------------------------------------------------------
# 2) metodo gaussiano FWHM
# ----------------------------------------------------------------------
def escalon_gaussiano(ch, B, sigma, C, fondo, pendiente):
    """
    Escalon borroneado por una gaussiana, con un fondo lineal.

        N(ch) = C/2 * erfc((ch - B) / (sqrt(2) sigma)) + fondo + pendiente*(ch-B)

    El `erfc` es la integral de la gaussiana: es EXACTAMENTE lo que se obtiene
    de convolucionar un escalon con la resolucion del detector. Por eso este
    modelo es el que corresponde a "ajustar el borde con una gaussiana": su
    derivada es una gaussiana de ancho sigma, y su punto de media altura cae en
    B por construccion.

    El termino lineal esta porque la lineshape de Klein-Nishina NO es plana
    antes del borde. Sin el, el ajuste inclina el escalon para absorber esa
    pendiente y sigma sale inflado.
    """
    z = (np.asarray(ch, float) - B) / (np.sqrt(2.0) * max(sigma, 1e-9))
    return 0.5 * C * erfc(z) + fondo + pendiente * (np.asarray(ch, float) - B)


def metodo_fwhm_gauss(canales, cuentas, region=None, errores=None,
                      p0=None, ancho_rel=0.45):
    """
    Ajusta el flanco con `escalon_gaussiano` y devuelve su media altura.

    El borde es el parametro B del ajuste, que por construccion es el punto de
    media altura del escalon; el mismo ajuste da `sigma`, y de ahi el FWHM que
    le da el nombre al metodo.

    region : (lo, hi). Si es None se deriva del propio espectro: se localiza el
             flanco con el metodo de la diferencia y se toma +-`ancho_rel`
             alrededor. Asi el metodo queda automatico y comparable con el NCC,
             que tampoco pide intervalo a mano.

    Devuelve dict con `borde`, `sigma`, `fwhm`, `fwhm_rel` y el ajuste completo.
    """
    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)

    if region is None:
        # el flanco primero, la ventana despues: sin esto habria que elegir a
        # mano justo lo que el metodo pretende determinar
        B0 = metodo_diferencia(canales, cuentas,
                               region=(0.15 * canales[-1], canales[-1]))["borde"]
        region = (max((1.0 - ancho_rel) * B0, canales[0]),
                  min((1.0 + ancho_rel) * B0, canales[-1]))
    else:
        B0 = 0.5 * (region[0] + region[1])

    m = _mascara(canales, region)
    x, y = canales[m], cuentas[m]
    if len(x) < 6:
        raise ValueError("region demasiado angosta para el ajuste gaussiano")
    sy = (np.asarray(errores, float)[m] if errores is not None else None)

    if p0 is None:
        alto = float(np.percentile(y, 90) - np.percentile(y, 10))
        p0 = [B0, max(0.08 * B0, 1e-3), max(alto, 1.0),
              float(np.percentile(y, 10)), 0.0]
    cotas = ([x[0], 1e-6, 0.0, -np.inf, -np.inf],
             [x[-1], x[-1] - x[0], np.inf, np.inf, np.inf])

    p, cov = curve_fit(escalon_gaussiano, x, y, p0=p0, sigma=sy,
                       absolute_sigma=sy is not None, bounds=cotas, maxfev=20000)
    err = np.sqrt(np.abs(np.diag(cov))) if cov is not None else np.full(5, np.nan)
    B, sigma = float(p[0]), float(abs(p[1]))
    return {"borde": B, "err_borde": float(err[0]), "sigma": sigma,
            "fwhm": 2.3548 * sigma, "fwhm_rel": 2.3548 * sigma / max(B, 1e-9),
            "C": float(p[2]), "fondo": float(p[3]), "pendiente": float(p[4]),
            "p": p, "cov": cov, "region": (float(x[0]), float(x[-1])),
            "metodo": "fwhm_gauss"}


# ----------------------------------------------------------------------
# 3) los tres juntos, que es como se usan
# ----------------------------------------------------------------------
def comparar(canales, cuentas, E_max, region=None, errores=None,
             ventana="auto", B_ref=None, ancho_rel=0.45):
    """
    Corre los tres metodos clasicos y devuelve el D de cada uno.

    E_max : energia del borde Compton [keV], para pasar canal -> D = B/E_max.
            El mapa por el origen es la convencion de todo el proyecto.
    B_ref : canal aproximado del borde (tipicamente el que estima el NCC). Si se
            da, la region se arma como (1-+ancho_rel)*B_ref.

    ⚠️ NO pasarle la ventana de AJUSTE del modelo extendido. Esa ventana
    arranca en ~0.25*B a proposito, porque ese modelo CONTIENE el exceso de
    bajo canal; los metodos clasicos no lo contienen, y la caida del exceso es
    mucho mas empinada que el borde. Medido sobre el Co-60 real: con la ventana
    de ajuste (317-1963 mV*ns) el minimo de la derivada se va al exceso y da
    D = 0.363 en vez de 1.265, o sea -71%. Con la region anclada al borde da
    1.265. Es el mismo motivo por el que `region` no tiene un default global.

    Devuelve {nombre: {"borde", "D", ...}}. El orden del dict es el orden en que
    conviene leerlos: de mas ingenuo a menos.
    """
    from etapas.model import punto_mitad_altura

    canales = np.asarray(canales, float)
    cuentas = np.asarray(cuentas, float)
    if region is None:
        if B_ref is None:
            B_ref = metodo_diferencia(
                canales, cuentas,
                region=(0.35 * canales[-1], canales[-1]))["borde"]
        region = (max((1.0 - ancho_rel) * B_ref, canales[0]),
                  min((1.0 + ancho_rel) * B_ref, canales[-1]))
    m = _mascara(canales, region)

    v = (_ventana_auto(canales, cuentas, region, 3)
         if ventana == "auto" else int(ventana))
    suave = _suavizar(cuentas, v, 3)
    salida = {}

    mitad = float(punto_mitad_altura(canales[m], suave[m]))
    salida["50% de altura"] = {"borde": mitad, "D": mitad / E_max,
                               "metodo": "mitad_altura"}

    dif = metodo_diferencia(canales, cuentas, region=region, ventana=ventana)
    dif["D"] = dif["borde"] / E_max
    salida["diferencia"] = dif

    try:
        gau = metodo_fwhm_gauss(canales, cuentas, region=region, errores=errores)
        gau["D"] = gau["borde"] / E_max
        salida["FWHM gaussiano"] = gau
    except Exception as e:                    # el flanco puede no ser ajustable
        salida["FWHM gaussiano"] = {"borde": float("nan"), "D": float("nan"),
                                    "error": str(e), "metodo": "fwhm_gauss"}
    return salida


# ======================================================================
#  DEMO: cuanto se equivoca cada metodo, con verdad conocida
# ======================================================================
if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from etapas.fit_single import _espectro_sintetico
    from etapas.model import forward_model
    from etapas.physics import compton_edge_energy

    print("=" * 72)
    print("METODOS CLASICOS - sesgo contra un borde de posicion CONOCIDA")
    print("=" * 72)

    E_GAMMA = 1173.2
    E_MAX = compton_edge_energy(E_GAMMA)
    B_VERDADERO = 600.0
    canales = np.linspace(0, 1000, 1024)

    print(f"\n  E_gamma = {E_GAMMA} keV  ->  E_max = {E_MAX:.1f} keV")
    print(f"  borde sintetico en B = {B_VERDADERO:.1f}  "
          f"(o sea D = {B_VERDADERO/E_MAX:.4f})\n")
    print(f"\n  {'sigma/B':>8}  {'50% altura':>18}  {'diferencia':>18}"
          f"  {'FWHM gauss':>18}  {'dif. SIN ruido':>16}")
    print("  " + "-" * 86)
    for s_rel in (0.03, 0.06, 0.09, 0.12, 0.18):
        sigma = s_rel * B_VERDADERO
        y = _espectro_sintetico(canales, [(E_GAMMA, 5000.0, B_VERDADERO, sigma)],
                                semilla=3).astype(float)
        res = comparar(canales, y, E_MAX)
        texto = ""
        for nombre in ("50% de altura", "diferencia", "FWHM gaussiano"):
            b = res[nombre]["borde"]
            texto += f"  {b:8.1f} ({100.0 * (b / B_VERDADERO - 1.0):+5.1f}%)"
        # el mismo metodo sobre el espectro SIN ruido: el sesgo del metodo solo
        ideal = forward_model(canales, E_GAMMA, 5000.0, B_VERDADERO, sigma)
        puro = metodo_diferencia(canales, ideal,
                                 region=(0.55 * B_VERDADERO, 1.45 * B_VERDADERO),
                                 ventana=7)["borde"]
        texto += f"    {100.0 * (puro / B_VERDADERO - 1.0):+8.1f}%"
        print(f"  {s_rel:8.2f}{texto}")
        if abs(s_rel - 0.09) < 1e-9:
            guardado = (y, res, sigma)

    print("\n  Lecturas:")
    print("  - LOS TRES SOBRESTIMAN, y el sesgo crece con la resolucion. No es")
    print("    ruido: la ultima columna es el metodo de la diferencia sobre el")
    print("    espectro SIN ruido, y va de +0.4% a +6.2%. El sesgo es del")
    print("    METODO, no de la estadistica.")
    print("  - La causa es la misma en los tres: la lineshape de Klein-Nishina")
    print("    NO es un escalon. Tiene una acumulacion creciente antes del")
    print("    borde, y ni el 50% de altura, ni el minimo de la derivada, ni un")
    print("    escalon con fondo lineal la contienen. Solo ajustar la lineshape")
    print("    entera (que es lo que hace el resto del proyecto) la contiene.")
    print("  - El de la diferencia es el menos sesgado de los tres, pero es el")
    print("    mas ruidoso: depende de cuanto se suavice, y suavizar de mas")
    print("    tambien sesga. Los otros dos son estables y peores.")

    # ---- por que la ventana de suavizado no puede ser una constante ----
    print("\n  El metodo de la diferencia con ventana FIJA (lo ingenuo):")
    print(f"  {'sigma/B':>8}  {'ventana 21 fija':>20}  {'ventana auto':>20}"
          f"  {'auto usa':>9}")
    print("  " + "-" * 66)
    for s_rel in (0.03, 0.06, 0.09, 0.12, 0.18):
        y = _espectro_sintetico(canales,
                                [(E_GAMMA, 5000.0, B_VERDADERO,
                                  s_rel * B_VERDADERO)], semilla=3).astype(float)
        reg = (0.55 * B_VERDADERO, 1.45 * B_VERDADERO)
        fija = metodo_diferencia(canales, y, region=reg, ventana=21)
        auto = metodo_diferencia(canales, y, region=reg, ventana="auto")
        print(f"  {s_rel:8.2f}   {fija['borde']:8.1f} "
              f"({100 * (fija['borde'] / B_VERDADERO - 1):+5.1f}%)"
              f"   {auto['borde']:8.1f} "
              f"({100 * (auto['borde'] / B_VERDADERO - 1):+5.1f}%)"
              f"   {auto['ventana']:5d} bins")
    print("\n  Con ventana fija el pozo de la derivada queda chato cuando el")
    print("  borde es ancho, aparecen minimos locales degenerados y el ruido")
    print("  elige cual gana: de ahi los saltos de signo. La ventana tiene que")
    print("  escalar con el ancho del flanco (`_ventana_auto`).")

    # ---------------------- figura ----------------------
    y, res, sigma = guardado
    suave = _suavizar(y, 21, 3)
    dif = res["diferencia"]
    gau = res["FWHM gaussiano"]

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    z = canales < 900
    ax[0].step(canales[z], y[z], where="mid", color="0.75", lw=0.9,
               label="espectro sintetico (Poisson)")
    ax[0].plot(canales[z], suave[z], color="C0", lw=1.4, label="suavizado")
    if np.isfinite(gau["borde"]):
        xf = np.linspace(*gau["region"], 400)
        ax[0].plot(xf, escalon_gaussiano(xf, *gau["p"]), color="C3", lw=1.8,
                   label="ajuste erfc (FWHM gaussiano)")
    ax[0].axvline(B_VERDADERO, color="k", lw=1.6, label=f"B verdadero = {B_VERDADERO:.0f}")
    for nombre, color in (("50% de altura", "C2"), ("diferencia", "C1"),
                          ("FWHM gaussiano", "C3")):
        b = res[nombre]["borde"]
        if np.isfinite(b):
            ax[0].axvline(b, color=color, ls="--", lw=1.2,
                          label=f"{nombre} = {b:.1f} ({100*(b/B_VERDADERO-1):+.1f}%)")
    ax[0].set_xlabel("canal [u.a.]")
    ax[0].set_ylabel("cuentas")
    ax[0].set_title(f"Los tres metodos clasicos  (sigma/B = {sigma/B_VERDADERO:.2f})")
    ax[0].legend(fontsize=7.5)
    ax[0].grid(alpha=0.3)

    ax[1].plot(canales[z], dif["derivada"][z], color="C1", lw=1.3,
               label="derivada del suavizado")
    ax[1].axhline(0, color="0.6", lw=0.7)
    ax[1].axvline(B_VERDADERO, color="k", lw=1.6, label="B verdadero")
    ax[1].axvline(dif["borde"], color="C1", ls="--", lw=1.3,
                  label=f"minimo = {dif['borde']:.1f}")
    ax[1].set_xlabel("canal [u.a.]")
    ax[1].set_ylabel("dN/dcanal")
    ax[1].set_title("Metodo de la diferencia: el minimo de la derivada")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(rutas.figura("clasicos_demo.png"), dpi=130)
    print("\n  Figura guardada: figuras/clasicos_demo.png")
