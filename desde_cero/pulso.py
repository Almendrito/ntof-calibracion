"""
Etapa 1 — mirar UN pulso: linea base, ruido, polaridad, pico, amplitud, onset.

    python desde_cero/pulso.py

Que resuelve esta etapa
-----------------------
De la forma de onda `(t, v)` de un evento (lo que entrega `lector.eventos`)
saca los numeros que describen al pulso, SIN integrar nada todavia:

    base      el nivel de la senal cuando no hay pulso              [mV]
    ruido     cuanto fluctua esa base (desviacion robusta)          [mV]
    signo     polaridad del pulso: -1 o +1, deducida, no supuesta
    pico      donde esta el maximo del pulso                        [ns]
    amplitud  altura del pico sobre la base, siempre positiva       [mV]
    onset     donde arranca el pulso (cruce del 10 % de la amplitud) [ns]

La Etapa 2 ancla la compuerta de integracion al onset de CADA evento y resta
la base de CADA evento, asi que un error aca pasa derecho al area.

La trampa que esta etapa existe para no repetir (CLAUDE.md §7.33)
-----------------------------------------------------------------
El analisis viejo tomaba la linea base como el promedio del primer 20 % del
registro (40 ns). Pero en el Co-60 de junio el pulso arranca a ~25 ns: esa
ventana FIJA se comia la subida del propio pulso, la base salia corrida hacia
el pulso, y al restarla e integrar 200 ns el error se multiplicaba por ~5.
Como el jitter del trigger mueve el pulso entre 18 y 54 ns, cada evento perdia
una fraccion distinta de su carga: el 26 % salia con el signo dado vuelta y el
espectro se ensanchaba al doble. La cura es que la ventana de linea base NO
sea fija sino RELATIVA AL ONSET de cada evento: todo lo que esta antes del
onset menos una guarda. Eso obliga a un orden en dos pasadas (ver `medir`):
primero localizar el pulso con una referencia provisoria, despues medir la
base buena, y recien ahi medir amplitud y onset definitivos.

Tiempos
-------
Todo tiempo se mide como `tau = t - t[0]`, relativo al inicio del registro
(el DRSOSC escribe `t` relativo a la celda de trigger, y `t[0]` cambia entre
corridas: ~0 ns en el Co-60 de junio, 27.6 ns en el Cs-137 del 26-ago).
"""

import os

import numpy as np


# --------------------------------------------------------------- parametros --
FRAC_ONSET = 0.10       # fraccion de la amplitud que define el onset
GUARDA_BASE = 4.0       # [ns] la base termina esta cantidad antes del onset
RECORTE = 4             # muestras que se tiran en cada extremo del registro
MIN_BASE = 10           # muestras minimas para confiar en la base
V_SATURA = 495.0        # [mV] el DRS4 cubre +-500 mV; arriba de esto recorta
K_PICO = 5.0            # un "pulso" menor que K_PICO*ruido es ruido
COLA = 80.0             # [ns] tras el onset: de ahi en adelante, cola sin pulso

# Banderas. Se MARCAN, no se descartan: descartar es una decision que depende
# de para que se use el evento, y esconderla aca haria un corte invisible.
# Mismos valores que analisis/pulso.py para poder comparar bandera a bandera.
SATURADO = 2            # alguna muestra toco el tope del digitalizador
SIN_PRE = 8             # no hay suficiente registro antes del onset para la base
SIN_PULSO = 16          # la "amplitud" no se distingue del ruido

# 1.4826 = 1 / Phi^-1(3/4). Para ruido gaussiano, la mediana de |x - mediana|
# (el MAD) vale sigma * Phi^-1(3/4) = 0.6745*sigma, porque la mitad de una
# gaussiana cae a menos de 0.6745 sigma del centro. Multiplicar por 1.4826
# convierte el MAD en una estimacion de sigma que coincide con la desviacion
# estandar cuando el ruido es gaussiano, pero que NO se entera si una o dos
# muestras se van lejos (un dark count, la cola de un pulso anterior): la
# desviacion estandar se infla con el cuadrado de esas excursiones, el MAD no.
MAD_A_SIGMA = 1.4826


# --------------------------------------------------------------------------
# 0. Registro
# --------------------------------------------------------------------------
def preparar(t, v, recorte=RECORTE):
    """
    `(tau, v)`: tiempo relativo al inicio del registro y bordes recortados.

    Por que recortar: las primeras celdas que lee el DRS4 dan una espiga que
    no es senal. Medida en 2000 eventos de cada corrida de trabajo, el
    percentil 99 de |v - mediana| en las primeras muestras es 60, 44, 29, 15
    y despues 2-4 mV (el ruido); al final, una sola muestra de 15-20 mV.

    Por que 4 y no 2 (el valor de analisis/, y el de esta etapa hasta el
    27-sep): con 2 quedan las muestras de 29 y 15 mV, y en un evento con un
    pulso de menos de ~30 mV la espiga es la excursion mas grande del
    registro. Entonces `buscar_pico` elige la espiga, el onset cae en tau = 0,
    no hay pre-pulso, la compuerta se sale del registro, y el pulso VERDADERO
    se ignora: area ~0 y bandera SIN_PRE|TRUNCADO. Pasaba en ~0.6 % de los
    eventos del Cs-137 del 26-ago y del Co-60 binario, todos de amplitud baja,
    o sea un corte que depende de la amplitud. Lo encontro la Etapa 2.
    """
    t = np.asarray(t, float)
    v = np.asarray(v, float)
    if recorte > 0 and t.size > 4 * recorte:
        t, v = t[recorte:-recorte], v[recorte:-recorte]
    return t - t[0], v


# --------------------------------------------------------------------------
# 1. Las tres piezas
# --------------------------------------------------------------------------
def linea_base(tau, v, hasta):
    """
    `(base, ruido, n)` con las muestras de `tau < hasta`.

    Por que MEDIANA y no media: con ruido gaussiano puro las dos estiman lo
    mismo (la media con un poco menos de varianza). La diferencia aparece
    cuando se cuela algo que no es ruido: el pie de un pulso previo, un dark
    count del fotodetector, una espiga del DRS4. Una sola muestra de 50 mV
    entre 100 corre la media 0.5 mV; a la mediana no la mueve. Y 0.5 mV de
    base mal puesta, integrados sobre 12 ns de compuerta, son 6 mV*ns en CADA
    evento, del orden del ancho de un bin del espectro.
    """
    m = tau < hasta
    x = v[m]
    if x.size == 0:
        return np.nan, np.nan, 0
    base = float(np.median(x))
    ruido = MAD_A_SIGMA * float(np.median(np.abs(x - base)))
    # piso: un tramo de muestras identicas daria ruido 0 y todos los umbrales
    # relativos al ruido se volverian "cualquier cosa es un pulso"
    return base, max(ruido, 1e-3), int(x.size)


def buscar_pico(v, referencia):
    """
    `(i_pico, amplitud, signo)` respecto de `referencia`.

    La polaridad se DEDUCE: es el signo de la excursion mas grande. Los PMT dan
    pulsos negativos, pero un inversor en la cadena, u otro detector, los da
    positivos, y suponerlo convierte el pulso en un "pozo" que el resto del
    analisis integra con el signo equivocado sin avisar.
    """
    d = v - referencia
    i = int(np.argmax(np.abs(d)))
    signo = -1.0 if d[i] < 0 else 1.0
    return i, float(abs(d[i])), signo


def buscar_onset(v, base, i_pico, signo, frac=FRAC_ONSET):
    """
    Indice del onset: la ULTIMA muestra por debajo de `frac`*amplitud antes
    del pico, buscando HACIA ATRAS desde el pico.

    Por que hacia atras: hacia adelante desde el principio del registro, la
    primera muestra que supere el 10 % puede ser una fluctuacion del ruido
    decenas de ns antes del pulso (con un pulso de 10 mV el umbral es 1 mV,
    del orden del ruido). Hacia atras desde el pico se encuentra el pie del
    flanco que efectivamente sube a ESE maximo; lo que haya antes no importa.

    Devuelve un indice de muestra, o sea el onset queda a la izquierda del
    cruce verdadero: sesgo de -dt/2 en promedio (~-0.1 ns). Para anclar una
    compuerta de 12 ns sobra; `onset_fino` interpola si hace falta mas.
    """
    d = signo * (v - base)
    idx = np.flatnonzero(d[:i_pico + 1] < frac * d[i_pico])
    return int(idx[-1]) if idx.size else 0


def onset_fino(tau, v, base, i_on, i_pico, signo, frac=FRAC_ONSET):
    """Cruce del `frac` interpolado linealmente entre `i_on` e `i_on+1`."""
    d = signo * (v - base)
    umbral = frac * d[i_pico]
    if i_on + 1 >= d.size or d[i_on + 1] == d[i_on]:
        return float(tau[i_on])
    x = (umbral - d[i_on]) / (d[i_on + 1] - d[i_on])
    return float(tau[i_on] + x * (tau[i_on + 1] - tau[i_on]))


# --------------------------------------------------------------------------
# 2. El pulso completo, en dos pasadas
# --------------------------------------------------------------------------
def medir(t, v, guarda=GUARDA_BASE, frac=FRAC_ONSET, cola=COLA):
    """
    Todo lo que se sabe de un pulso. Devuelve un dict:

        base, ruido, n_base, signo, i_pico, t_pico, amplitud,
        i_onset, t_onset, t_onset_fino, bandera, tau, v

    (`tau` y `v` son el registro ya preparado, para que quien dibuje o integre
    use exactamente las mismas muestras.)

    Por que dos pasadas: la base se mide ANTES del onset, pero para saber
    donde esta el onset hace falta una base. Se rompe el circulo con una
    referencia provisoria que el pulso no puede contaminar mucho: la mediana
    de TODO el registro. El pulso ocupa ~10 % de las muestras, y la mediana no
    se entera de lo que pase en el 10 % de un extremo. Con esa referencia se
    ubica el pulso; con el onset provisorio se mide la base buena; y con la
    base buena se vuelven a medir pico, amplitud y onset.
    """
    tau, v = preparar(t, v)
    bandera = 0
    if np.max(np.abs(v)) >= V_SATURA:
        bandera |= SATURADO

    # pasada 1: ubicar el pulso con la mediana del registro como referencia
    ref = float(np.median(v))
    i_p, _, signo = buscar_pico(v, ref)
    i_on = buscar_onset(v, ref, i_p, signo, frac)

    # la base buena: todo lo anterior al onset, menos una guarda. La guarda
    # existe porque el onset al 10 % no es donde el pulso EMPIEZA sino donde
    # ya subio un 10 %; las muestras justo antes ya tienen algo de pulso.
    base, ruido, n = linea_base(tau, v, tau[i_on] - guarda)
    if n < MIN_BASE:
        # El pulso llego al principio del registro. No se descarta: el trigger
        # tiene walk (los pulsos grandes disparan antes), asi que descartar por
        # posicion seria un corte que depende de la amplitud, y deformaria
        # justo el espectro que se quiere calibrar. Se usa la COLA, donde el
        # pulso ya se apago, y se marca.
        bandera |= SIN_PRE
        c = tau > tau[i_on] + cola + guarda
        if c.sum() >= MIN_BASE:
            base = float(np.median(v[c]))
            ruido = max(MAD_A_SIGMA * float(np.median(np.abs(v[c] - base))), 1e-3)
        elif not np.isfinite(base):
            base, ruido = ref, np.nan

    # pasada 2: todo de nuevo con la base buena
    i_p, amp, signo = buscar_pico(v, base)
    if np.isfinite(ruido) and amp < K_PICO * ruido:
        bandera |= SIN_PULSO
    i_on = buscar_onset(v, base, i_p, signo, frac)

    return dict(base=base, ruido=ruido, n_base=n, signo=signo,
                i_pico=i_p, t_pico=float(tau[i_p]), amplitud=amp,
                i_onset=i_on, t_onset=float(tau[i_on]),
                t_onset_fino=onset_fino(tau, v, base, i_on, i_p, signo, frac),
                bandera=bandera, tau=tau, v=v)


# --------------------------------------------------------------------------
# Demo: seis pulsos reales con todo marcado
# --------------------------------------------------------------------------
def _elegir(ruta, n=400):
    """De los primeros `n` eventos, uno de amplitud baja y uno alta."""
    from desde_cero.lector import eventos
    ms = [medir(t, v) for t, v in eventos(ruta, limite=n)]
    ms = [m for m in ms if m["bandera"] == 0]
    a = np.array([m["amplitud"] for m in ms])
    orden = np.argsort(a)
    return [ms[orden[len(ms) // 4]], ms[orden[int(0.97 * len(ms))]]]


def _demo():
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero.lector import CORRIDAS_TRABAJO

    fig, ejes = plt.subplots(2, 3, figsize=(15, 7.5))
    for col, (clave, (ruta, _)) in enumerate(CORRIDAS_TRABAJO.items()):
        if not os.path.exists(ruta):
            ejes[0, col].set_title("%s: no esta en disco" % clave)
            continue
        for fila, m in enumerate(_elegir(ruta)):
            ax = ejes[fila, col]
            tau, v = m["tau"], m["v"]
            ax.plot(tau, v, lw=0.8, color="0.25", label="registro")
            fin_base = m["t_onset"] - GUARDA_BASE
            ax.axvspan(tau[0], fin_base, color="tab:blue", alpha=0.10,
                       label="ventana de base")
            ax.axhline(m["base"], color="tab:blue", lw=1, label="base (mediana)")
            ax.axhspan(m["base"] - m["ruido"], m["base"] + m["ruido"],
                       color="tab:blue", alpha=0.2)
            ax.axvline(m["t_onset"], color="tab:red", lw=1,
                       label="onset (10 %)")
            ax.plot(m["t_pico"], m["base"] + m["signo"] * m["amplitud"], "o",
                    color="tab:orange", label="pico")
            # zoom alrededor del pulso, sin perder la base
            ax.set_xlim(max(tau[0], m["t_onset"] - 30), m["t_onset"] + 60)
            ax.set_title("%s   A = %.1f mV   ruido = %.2f mV   signo %+d"
                         % (clave, m["amplitud"], m["ruido"], m["signo"]),
                         fontsize=9)
            ax.set_xlabel("tau [ns]")
            ax.set_ylabel("v [mV]")
            ax.grid(alpha=0.3)
            print("  %-12s A = %7.2f mV   base = %7.2f mV   ruido = %.2f mV   "
                  "onset = %6.2f ns   pico = %6.2f ns   signo %+d"
                  % (clave, m["amplitud"], m["base"], m["ruido"], m["t_onset"],
                     m["t_pico"], m["signo"]))
    ejes[0, 0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras")
    os.makedirs(salida, exist_ok=True)
    salida = os.path.join(salida, "etapa1_pulsos.png")
    fig.savefig(salida, dpi=120)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
