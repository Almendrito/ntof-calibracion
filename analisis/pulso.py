"""
pulso.py — ETAPA 1: de un pulso crudo a un area, y sobre todo: QUE area.

Aca se decide el observable de todo el proyecto. La cadena, en el orden en que
la pediste:

    1. intervalo PREVIO al pulso  ->  promedio  ->  se resta (quita el offset)
    2. integral acumulada Q(t)    ->  donde EMPIEZA a crecer -> t_ini
                                  ->  donde CONVERGE          -> t_fin
    3. area = Q(t_fin) - Q(t_ini)

La idea del punto 2 merece un parrafo, porque es el corazon del modulo.
`Q(t) = integral de (v - base)` tiene una propiedad util: su derivada es la
propia senal. O sea "la integral dejo de crecer" y "la senal volvio a la linea
base" son la MISMA afirmacion... pero no el mismo numero. Mirar la senal es
mirar un punto; mirar la integral es mirar todo lo acumulado hasta ese punto,
que es exactamente la cantidad que despues se histogramea. Por eso el criterio
de convergencia se enuncia sobre Q: cortar donde Q hace plateau garantiza que
no se esta dejando carga afuera, y eso es lo unico que le importa al espectro.

Que se mide y que se supone
---------------------------
NADA esta cableado a este detector: la polaridad se deduce del propio pulso, el
onset sale de la senal y la ventana sale de la curva de convergencia medida
sobre el ensemble. Los unicos numeros fijos son umbrales de calidad
(saturacion del digitalizador) y la fraccion del pico que define el onset.

La trampa que este modulo existe para no repetir
------------------------------------------------
En el analisis anterior la linea base se promediaba sobre el primer 20% del
registro (40 ns). Pero en el Co-60 el pulso arranca a ~25 ns: la ventana de
linea base se comia la SUBIDA del propio pulso. El sesgo se multiplica al
integrar, y como el jitter del trigger mueve el pico entre 20 y 55 ns, cada
evento perdia una fraccion DISTINTA de su carga -> ensanchamiento puro del
espectro (26% de los eventos salian incluso con el signo dado vuelta). De ahi
que aca la linea base se tome siempre RELATIVA al onset de cada evento y nunca
en una ventana absoluta del registro. La demo reproduce el error a proposito
para que se vea la diferencia.

Correr `python analisis/pulso.py` deja `figuras/01_ventana_integracion.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, drs

_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


# --------------------------------------------------------------- parametros --
V_SATURA = 495.0        # |v| a la que satura el DRS4 [mV]
FRAC_ONSET = 0.10       # fraccion de la amplitud que define el onset
RECORTE_BORDES = 2      # muestras a tirar en cada extremo del registro

# Banderas de calidad. Se MARCA, no se descarta: descartar es una decision de
# la Etapa 2, y hacerla aca escondería un corte dentro del calculo del area.
MAL_SATURA = 2          # el digitalizador llego al tope
MAL_TRUNCADO = 4        # la ventana se sale del registro (pulso cortado)
MAL_SIN_PRE = 8         # no hay pre-pulso suficiente para la linea base
MAL_SIN_PICO = 16       # no se distingue un pulso del ruido


# --------------------------------------------------------------------------
# 1. Preparacion del registro
# --------------------------------------------------------------------------
def preparar(t, v, recorte=RECORTE_BORDES):
    """
    Devuelve `(tau, v)` con el tiempo RELATIVO al inicio del registro y los
    bordes recortados.

    Dos razones, las dos medidas en los datos:

    - `tau = t - t[0]`. El DRSOSC escribe el eje relativo a la celda de trigger,
      asi que `t[0]` no es cero ni es el mismo entre corridas: en el Co-60 vale
      ~-0.01 ns y en el Cs-137 del 26-ago vale 27.6 ns. Todo se mide contra
      `tau`, nunca contra `t` absoluto, o las ventanas de dos corridas dejan de
      significar lo mismo.
    - Recorte de bordes. La primera y la ultima celda del anillo del DRS4 dan
      espigas de decenas de mV que no son senal (se ven claritas en
      `figuras/00_pulsos_crudos.png`, en t=0 y t=200 del Cs-137). Si caen en la
      ventana de linea base, envenenan el offset de ese evento.
    """
    if recorte > 0 and t.size > 4 * recorte:
        t = t[recorte:-recorte]
        v = v[recorte:-recorte]
    return t - t[0], v


# --------------------------------------------------------------------------
# 2. Linea base: el intervalo previo al pulso
# --------------------------------------------------------------------------
def linea_base(tau, v, hasta, desde=None, estadistico="mediana", minimo=10):
    """
    Promedia el intervalo PREVIO al pulso y devuelve `(base, ruido, n)`.

    hasta : borde derecho del intervalo [ns]. Tiene que quedar ANTES de que el
            pulso empiece a subir; en la practica se pasa `t_onset - guarda`.
    estadistico : "mediana" (default) o "media".
        La media es la receta literal y es la estimacion optima si el intervalo
        tiene solo ruido gaussiano. La mediana cuesta lo mismo y es inmune a que
        se cuele un dark count o el pie de un pulso anterior, que es algo que en
        un SiPM pasa seguido. Con ruido limpio las dos coinciden dentro de
        0.25*sigma; la demo compara las dos sobre los datos reales.
    minimo : si el intervalo tiene menos de estas muestras, la estimacion no es
        confiable; se devuelve n < minimo y quien llama decide (marca MAL_SIN_PRE).

    `ruido` es la desviacion robusta 1.4826*MAD, que para ruido gaussiano es
    sigma. Se usa despues como ESCALA de todos los umbrales: asi ningun criterio
    del modulo depende de cuantos mV mide este detector en particular.
    """
    m = tau < hasta
    if desde is not None:
        m &= tau >= desde
    muestras = v[m]
    n = muestras.size
    if n == 0:
        return float(np.median(v)), float(1.4826 * np.median(np.abs(v - np.median(v)))), 0
    if estadistico == "media":
        base = float(np.mean(muestras))
    else:
        base = float(np.median(muestras))
    ruido = float(1.4826 * np.median(np.abs(muestras - np.median(muestras))))
    return base, max(ruido, 1e-3), int(n)


# --------------------------------------------------------------------------
# 3. Pico y onset
# --------------------------------------------------------------------------
def detectar_pico(tau, v, base=None):
    """
    Devuelve `(i_pico, t_pico, amplitud, signo)`.

    El signo se DEDUCE de la excursion mas grande respecto de la mediana del
    registro; no se supone "negativo" en ninguna parte. La mediana sirve de
    referencia provisoria porque el pulso ocupa ~10% de las muestras y la
    mediana no lo ve.
    """
    if base is None:
        base = float(np.median(v))
    d = v - base
    i = int(np.argmax(np.abs(d)))
    signo = -1.0 if d[i] < 0 else 1.0
    return i, float(tau[i]), float(abs(d[i])), signo


def detectar_onset(tau, v, base, i_pico, signo, frac=FRAC_ONSET):
    """
    Onset = ULTIMO cruce de `frac`*amplitud RETROCEDIENDO desde el pico.

    Retroceder importa. Buscando hacia adelante desde el principio del registro,
    cualquier fluctuacion de la linea base que pase el umbral dispara un onset
    falso decenas de ns antes del pulso, y con eso la ventana entera se corre.
    Retrocediendo desde el pico se encuentra el pie del flanco que efectivamente
    lleva a ESE maximo.
    """
    d = signo * (v - base)
    amp = d[i_pico]
    idx = np.flatnonzero(d[:i_pico + 1] < frac * amp)
    i_on = int(idx[-1]) if idx.size else 0
    return i_on, float(tau[i_on])


# --------------------------------------------------------------------------
# 4. La integral acumulada y su convergencia
# --------------------------------------------------------------------------
def integral_acumulada(tau, v, base=0.0, signo=1.0):
    """
    Q(t) = integral de `signo*(v - base)` desde el inicio del registro, con el
    trapecio y sobre el eje `tau` REAL.

    Nunca `np.cumsum(v)*dt`: el dt del DRS4 va de 0.10 a 0.32 ns entre celdas en
    la corrida de Co-60 (154% de variacion, medido en la Etapa 0), asi que un dt
    constante mete un error de forma que se ve directo en el espectro.

    Devuelve un array del mismo largo que `tau`, con Q[0] = 0.
    """
    d = signo * (np.asarray(v, float) - base)
    paso = 0.5 * (d[1:] + d[:-1]) * np.diff(tau)
    Q = np.empty(tau.size)
    Q[0] = 0.0
    np.cumsum(paso, out=Q[1:])
    return Q


def _suavizar(tau, v, ancho):
    """Promedio movil de `ancho` ns. Solo para DECIDIR umbrales, jamas para integrar."""
    if ancho <= 0:
        return v
    n = max(int(round(ancho / max(float(np.median(np.diff(tau))), 1e-9))), 1)
    if n <= 1:
        return v
    nucleo = np.ones(n) / n
    return np.convolve(v, nucleo, mode="same")


def ventana_por_convergencia(tau, v, base, signo, ruido, k=3.0, tramo=4.0,
                             ancho_suave=1.0, i_pico=None):
    """
    Encuentra `(t_ini, t_fin)` mirando donde la integral acumulada EMPIEZA a
    crecer y donde CONVERGE. Es el criterio que define el area.

    Como se enuncia
    ---------------
    Q crece exactamente donde `v - base` se aparta de cero, asi que las dos
    preguntas se contestan sobre la senal suavizada, con el ruido de la propia
    linea base como vara:

      t_ini : ULTIMO instante antes del pico en que la senal estuvo quieta
              (|v - base| < k*ruido) durante un tramo continuo de `tramo` ns.
              Pedir un TRAMO y no un solo punto es lo que evita que una
              fluctuacion aislada en mitad del flanco se haga pasar por "todavia
              no empezo".
      t_fin : PRIMER instante despues del pico en que la senal vuelve a estar
              quieta durante `tramo` ns seguidos. Un solo punto plano no alcanza:
              los pulsos de centellador tienen un hombro en la cola y ahi la
              derivada pasa cerca de cero sin que la integral haya convergido.

    Parametros
    ----------
    k : umbral en unidades del ruido de la linea base. Con k=3 se corta cuando
        la senal es indistinguible de la linea base a 3 sigma.
    tramo : cuantos ns seguidos hay que estar quieto para declarar convergencia.
    ancho_suave : promedio movil [ns] aplicado SOLO para evaluar el criterio.
        Sin el, el ruido muestra a muestra cruza el umbral todo el tiempo y el
        `tramo` casi nunca se cumple.

    Devuelve `(t_ini, t_fin, diag)`. `diag` trae la carga total del registro, la
    carga dentro de la ventana y la fraccion entre ambas, que es lo que permite
    GRAFICAR la convergencia en vez de creerle al criterio.
    """
    d = signo * (v - base)
    s = _suavizar(tau, d, ancho_suave)
    if i_pico is None:
        i_pico = int(np.argmax(s))

    dt = max(float(np.median(np.diff(tau))), 1e-9)
    n_tramo = max(int(round(tramo / dt)), 2)
    quieto = np.abs(s) < k * ruido

    # Rachas de `n_tramo` puntos quietos consecutivos, en una sola pasada:
    # `racha[i]` es True si desde i hay n_tramo puntos quietos seguidos.
    acum = np.concatenate([[0], np.cumsum(quieto)])
    racha = (acum[n_tramo:] - acum[:-n_tramo]) == n_tramo   # largo N - n_tramo + 1

    # --- t_ini: la ultima racha que TERMINA antes del pico -------------------
    ini_ok = np.flatnonzero(racha[:max(i_pico - n_tramo + 1, 1)])
    i_ini = int(ini_ok[-1] + n_tramo - 1) if ini_ok.size else 0

    # --- t_fin: la primera racha que EMPIEZA despues del pico ----------------
    fin_ok = np.flatnonzero(racha[i_pico:])
    i_fin = int(i_pico + fin_ok[0]) if fin_ok.size else tau.size - 1

    Q = integral_acumulada(tau, v, base, signo)
    diag = {
        "i_ini": i_ini, "i_fin": i_fin,
        "Q_total": float(Q[-1]),
        "Q_ventana": float(Q[i_fin] - Q[i_ini]),
        "ancho": float(tau[i_fin] - tau[i_ini]),
        "convergio": bool(fin_ok.size),
        "hubo_pre": bool(ini_ok.size),
    }
    diag["fraccion"] = (diag["Q_ventana"] / diag["Q_total"]
                        if abs(diag["Q_total"]) > 1e-12 else np.nan)
    return float(tau[i_ini]), float(tau[i_fin]), diag


# --------------------------------------------------------------------------
# 5. El pulso promedio: la convergencia medida con estadistica
# --------------------------------------------------------------------------
def pulso_promedio(T, V, tau_min=-25.0, tau_max=140.0, paso=0.05,
                   amp_minima=20.0, guarda_base=4.0, normalizar=True):
    """
    Alinea los pulsos por su onset y devuelve `(tau_rel, forma, n_usados)`.

    Por que hace falta alinear: el jitter del trigger mueve el pico entre ~20 y
    ~55 ns en la corrida de Co-60. Promediar los registros crudos no da el pulso
    promedio sino el pulso convolucionado con el jitter, que es mas ancho y
    converge mas tarde. La ventana que saliera de ahi seria mas larga de lo
    necesario y meteria ruido de mas en cada area.

    Cada pulso se lleva a una grilla comun `tau - t_onset` por interpolacion
    lineal (obligatorio: el eje del DRS4 no es uniforme, asi que no se puede
    alinear corriendo indices) y se normaliza por su amplitud, de modo que el
    promedio sea la FORMA y no este dominado por los pulsos mas grandes.

    `amp_minima` deja afuera los pulsos chicos, cuyo onset es ruidoso. Es solo
    para MEDIR la forma; el area despues se calcula sobre todos los eventos.
    """
    grilla = np.arange(tau_min, tau_max + 0.5 * paso, paso)
    suma = np.zeros_like(grilla)
    n = 0
    for t, v in zip(T, V):
        tau, vv = preparar(t, v)
        i_p, _, amp, signo = detectar_pico(tau, vv)
        if amp < amp_minima:
            continue
        base0 = float(np.median(vv))
        i_on, t_on = detectar_onset(tau, vv, base0, i_p, signo)
        base, ruido, n_base = linea_base(tau, vv, hasta=t_on - guarda_base)
        if n_base < 10:
            continue
        i_p, _, amp, signo = detectar_pico(tau, vv, base)
        if amp < amp_minima:
            continue
        d = signo * (vv - base)
        suma += np.interp(grilla, tau - t_on, d / (amp if normalizar else 1.0),
                          left=0.0, right=0.0)
        n += 1
    if n == 0:
        raise ValueError("ningun pulso paso el corte de amplitud (%.0f mV)"
                         % amp_minima)
    return grilla, suma / n, n


def curva_convergencia(tau_rel, forma, t_ini=None, referencia="plateau"):
    """
    Fraccion de la carga acumulada en funcion de donde se corte.

    Devuelve `(tau_rel, fraccion, info)`.

    `referencia` decide contra que se normaliza, y NO es un detalle cosmetico:

      "plateau" (default) : contra el MAXIMO de Q. Es lo correcto, porque el
          pulso de este SiPM tiene un rebote de polaridad opuesta despues de la
          cola: Q sube, llega a un techo, y despues DECRECE lentamente. Medido
          en el Co-60, ese rebote se lleva ~1.7% del area. Normalizando contra
          el final de la grilla, la curva pasa por 101.7% y "converger al 100%"
          deja de significar nada.
      "total" : contra Q al final de la grilla. Sirve para VER el rebote.

    `info` trae `t_plateau` (donde Q es maxima), `sobrepaso` (cuanto sube la
    curva por encima de 1 antes de bajar) y `caida_final`.
    """
    Q = integral_acumulada(tau_rel, forma)
    i0 = 0 if t_ini is None else int(np.searchsorted(tau_rel, t_ini))
    Qv = Q - Q[i0]
    i_max = int(np.argmax(Qv))
    ref = Qv[i_max] if referencia == "plateau" else Qv[-1]
    frac = Qv / ref if abs(ref) > 1e-12 else Qv
    info = {"t_plateau": float(tau_rel[i_max]),
            "sobrepaso": float(Qv[i_max] / Qv[-1] - 1.0) if abs(Qv[-1]) > 1e-12 else np.nan,
            "caida_final": float(1.0 - Qv[-1] / Qv[i_max]) if abs(Qv[i_max]) > 1e-12 else np.nan,
            # Si el maximo cae en el ULTIMO tramo de la grilla, la integral no
            # convergio ADENTRO de lo que se miro: el "100%" es una convencion,
            # no una medicion, y hay que decirlo en vez de reportar el numero.
            "en_borde": bool(i_max >= tau_rel.size - 3)}
    return tau_rel, frac, info


# --------------------------------------------------------------------------
# 6. El area de UN evento, con las dos ventanas
# --------------------------------------------------------------------------
def area_evento(t, v, pre=None, post=None, guarda_base=4.0,
                estadistico="mediana", k=3.0, tramo=4.0, ancho_suave=1.0,
                frac_onset=FRAC_ONSET, v_satura=V_SATURA, recorte=RECORTE_BORDES):
    """
    Calcula el area de un evento y devuelve un dict con todo lo que se midio.

    Dos modos, segun se pase o no `(pre, post)`:

    * `pre`/`post` dados  -> METODO GLOBAL. La ventana es `[onset-pre,
      onset+post]` con el ancho fijo que la Etapa 1 midio sobre el ensemble,
      anclada al onset DE ESTE evento. El jitter del trigger no lo toca.
    * `pre`/`post` = None -> METODO POR CONVERGENCIA. Este evento decide su
      propia ventana con `ventana_por_convergencia`. Es la receta literal;
      sobre pulsos chicos el criterio se vuelve ruidoso, y medir cuanto se
      vuelve es justamente el punto de compararlos.

    Claves del dict: area [mV*ns], amplitud [mV], t_onset, t_ini, t_fin, base,
    ruido, bandera, y `fraccion` (que parte de la carga del registro quedo
    adentro; solo en el modo por convergencia).
    """
    tau, v = preparar(t, v, recorte)
    if tau.size < 16:
        return dict(area=0.0, amplitud=0.0, t_onset=np.nan, t_ini=np.nan,
                    t_fin=np.nan, base=np.nan, ruido=np.nan,
                    bandera=MAL_TRUNCADO, fraccion=np.nan)

    # --- 1. localizar el pulso con una referencia provisoria -----------------
    base0 = float(np.median(v))          # el pulso ocupa ~10% del registro
    i_p, t_pico, amp0, signo = detectar_pico(tau, v, base0)
    i_on, t_on = detectar_onset(tau, v, base0, i_p, signo, frac_onset)

    bandera = 0
    if np.max(np.abs(v)) >= v_satura:
        bandera |= MAL_SATURA

    # --- 2. linea base en el intervalo PREVIO al pulso ----------------------
    base, ruido, n_base = linea_base(tau, v, hasta=t_on - guarda_base,
                                     estadistico=estadistico)
    if n_base < 10:
        # Sin pre-pulso se cae al post-pulso (la cola despues de la ventana).
        # No se DESCARTA el evento: rechazar por posicion del onset seria un
        # corte que depende de la amplitud, porque el trigger tiene walk (los
        # pulsos grandes cruzan el umbral antes), y un corte dependiente de la
        # amplitud deforma justo lo que se quiere calibrar.
        bandera |= MAL_SIN_PRE
        cola = tau > (t_on + (post if post is not None else 80.0) + guarda_base)
        if cola.sum() >= 10:
            base = float(np.median(v[cola]))
            ruido = max(float(1.4826 * np.median(np.abs(v[cola] - base))), 1e-3)

    # --- 3. re-medir el pulso con la linea base buena -----------------------
    i_p, t_pico, amp, signo = detectar_pico(tau, v, base)
    if amp < 5.0 * ruido:
        bandera |= MAL_SIN_PICO
    i_on, t_on = detectar_onset(tau, v, base, i_p, signo, frac_onset)

    # --- 4. la ventana -------------------------------------------------------
    fraccion = np.nan
    if pre is not None and post is not None:
        t_ini, t_fin = t_on - pre, t_on + post
    else:
        t_ini, t_fin, diag = ventana_por_convergencia(
            tau, v, base, signo, ruido, k=k, tramo=tramo,
            ancho_suave=ancho_suave, i_pico=i_p)
        fraccion = diag["fraccion"]
        if not diag["convergio"]:
            bandera |= MAL_TRUNCADO

    if t_ini < tau[0] or t_fin > tau[-1]:
        bandera |= MAL_TRUNCADO
    i_ini = int(np.clip(np.searchsorted(tau, t_ini), 0, tau.size - 2))
    i_fin = int(np.clip(np.searchsorted(tau, t_fin), i_ini + 2, tau.size - 1))

    # --- 5. integrar ---------------------------------------------------------
    area = signo * _TRAPZ(v[i_ini:i_fin + 1] - base, x=tau[i_ini:i_fin + 1])

    return dict(area=float(area), amplitud=float(amp), t_onset=float(t_on),
                t_ini=float(tau[i_ini]), t_fin=float(tau[i_fin]),
                base=float(base), ruido=float(ruido), bandera=int(bandera),
                fraccion=float(fraccion))


# --------------------------------------------------------------------------
# 7. TERCER METODO: ajuste de una plantilla del pulso, evento a evento
# --------------------------------------------------------------------------
# Por que un tercer metodo (24-sep-2026)
# -------------------------------------
# La compuerta comun de 12 ns es lineal y comparable entre corridas, pero mide
# una FRACCION de la carga (~81-83% de lo acumulado hasta 40 ns, CLAUDE.md
# 7.66), no la carga. Y medir la carga entera integrando no se puede: el pulso
# necesita ~100 ns tras el onset para converger y el registro deja 47-158 ns
# segun la sesion. El metodo por convergencia lo intentaba evento a evento y
# fallaba porque cada evento corta donde su senal se hunde en SU ruido (13 ns
# en los chicos, 42 en los grandes: CLAUDE.md 7.47).
#
# La salida es no integrar la carga sino AJUSTARLA: si todos los pulsos tienen
# la misma forma, un evento queda descrito por tres numeros,
#
#     v(t) = b + A * f(t - t0)
#
# con f la forma normalizada a amplitud 1. La carga es A * integral(f), y esa
# integral se mide UNA vez sobre la forma promedio de las corridas de registro
# largo. Lo que el registro de un evento no contiene (la cola, en las sesiones
# de septiembre) se recupera por la forma, no por la integral: el ajuste usa
# el tramo que hay, que es donde esta casi toda la informacion (el pico).
#
# Lo que se SUPONE, y hay que verificar antes de usarlo: que f no depende del
# detector, del voltaje ni de la amplitud. Si dependiera de la amplitud, A dejaria
# de ser proporcional a la carga y el metodo romperia la linealidad igual que el
# de convergencia; por eso la validacion mide area/amplitud y el cociente contra
# la compuerta evento a evento.

T_PLANTILLA = (-20.0, 150.0)   # grilla de la forma, relativa al onset [ns]
PASO_PLANTILLA = 0.05          # [ns]; el DRS4 muestrea a ~0.2 ns
AMP_PLANTILLA = 20.0           # [mV] pulsos que entran a MEDIR la forma
# Tramo que decide t0 [ns desde el onset]. El t0 lo fijan la subida y el pico;
# la cola casi no cambia con un corrimiento de fracciones de ns y solo agrega
# muestras (costo) sin agregar informacion sobre t0.
TRAMO_T0 = (-10.0, 40.0)
# Hasta donde se cuenta la carga de la plantilla [ns desde el onset]. NO es "el
# pulso entero": en la configuracion de ago/sep el pulso no converge dentro de
# ningun registro (a 130 ns la integral acumulada todavia sube ~0.1%/ns por una
# cola casi plana de ~1.5% de la amplitud), asi que "hasta el plateau" no esta
# definido. Un T fijo, el MISMO para todas las plantillas, es la unica
# definicion que significa lo mismo en todas las corridas. 120 ns cabe en los
# registros largos de las dos configuraciones (junio llega a ~150, el 26-ago a
# ~130).
T_CARGA = 120.0
# Busqueda de t0 alrededor del onset detectado: gruesa y despues fina. El onset
# de un pulso de 5 mV con 1 mV de ruido puede estar corrido varios ns, por eso
# la gruesa cubre +-4 ns.
T0_GRUESO = np.arange(-4.0, 4.0001, 0.4)
T0_FINO = np.arange(-0.4, 0.4001, 0.04)


def cruce_onset(tau, d, i_pico, frac=FRAC_ONSET):
    """
    Onset SUB-MUESTRA: el cruce de `frac`*amplitud interpolado linealmente
    entre las dos muestras que lo rodean.

    `detectar_onset` devuelve la muestra, y eso alcanza para anclar una
    compuerta de 12 ns, pero no para PROMEDIAR pulsos: alinear por la muestra
    mete un jitter uniforme de un paso (~0.2 ns) que ensancha la subida de la
    forma promedio. `d` ya tiene la linea base restada y el signo corregido.
    """
    amp = d[i_pico]
    idx = np.flatnonzero(d[:i_pico + 1] < frac * amp)
    if idx.size == 0:
        return float(tau[0])
    i = int(idx[-1])
    if i + 1 >= d.size or d[i + 1] == d[i]:
        return float(tau[i])
    x = (frac * amp - d[i]) / (d[i + 1] - d[i])
    return float(tau[i] + x * (tau[i + 1] - tau[i]))


def _ajuste_lineal(tau, d, grilla, forma, t0s):
    """
    Para cada t0 de `t0s`, ajusta `d = b + A*f(tau - t0)` por minimos cuadrados
    LINEALES y devuelve `(A, b, chi2_sin_normalizar)` como arrays del largo de
    `t0s`.

    A y b tienen solucion cerrada (dos parametros lineales), asi que la unica
    busqueda no lineal es la de t0, y se vectoriza: una sola llamada a
    `np.interp` sobre la matriz (t0, muestra). La plantilla se evalua sobre el
    eje de tiempo REAL del evento, porque el DRS4 no muestrea uniforme (dt de
    0.02 a 0.34 ns): interpolar el evento a una grilla uniforme meteria un
    error de forma en la subida, que es justo donde se decide t0.

    Pesos uniformes por MUESTRA, no por dt: el ruido es por celda del DRS4, no
    por ns, y el chi2 correcto cuenta muestras.
    """
    P = np.interp(tau[None, :] - t0s[:, None], grilla, forma, left=0.0, right=0.0)
    n = float(tau.size)
    Sp = P.sum(axis=1)
    Spp = np.einsum("ij,ij->i", P, P)
    Spv = P @ d
    Sv = float(d.sum())
    Svv = float(d @ d)
    det = n * Spp - Sp * Sp
    det = np.where(np.abs(det) < 1e-12, np.nan, det)
    A = (n * Spv - Sp * Sv) / det
    b = (Sv - A * Sp) / n
    # residuo de la solucion de minimos cuadrados: sum(d^2) - A*sum(p d) - b*sum(d)
    chi2 = Svv - A * Spv - b * Sv
    return A, b, chi2


def ajustar_plantilla(tau, d, plantilla, t_on, tramo_t0=TRAMO_T0):
    """
    Ajuste de un evento: `d = b + A * f(tau - t0)`. `d` es la senal con el
    signo ya corregido (el pulso positivo), SIN restar la linea base: b es un
    parametro del ajuste, y restarla antes haria depender el resultado de la
    ventana de linea base, que es lo que la compuerta ya mide por su lado.

    Tres pasos:
      1. t0 en una grilla GRUESA (+-4 ns) y despues FINA (+-0.4 ns) alrededor
         del onset detectado, ajustando solo `tramo_t0` (subida y pico).
      2. Refinamiento parabolico sobre los tres puntos finos del minimo:
         lleva la resolucion de 0.04 ns a ~0.005 sin mas evaluaciones.
      3. A y b FINALES con t0 fijo, sobre TODO el tramo que el registro
         contiene y la plantilla cubre. Es aca donde entra la cola: donde el
         registro la tiene, pesa; donde no, la pone la plantilla.

    Devuelve un dict con A, b, t0, chi2 (suma de residuos al cuadrado),
    n (muestras del ajuste final), `t_ult` (ultimo instante del registro
    relativo a t0) y `en_borde` (el minimo en t0 cayo en el extremo de la
    grilla gruesa: el onset detectado estaba lejos, tipico de pulsos al nivel
    del ruido).
    """
    grilla, forma = plantilla["tau"], plantilla["forma"]
    m = (tau >= t_on + tramo_t0[0]) & (tau <= t_on + tramo_t0[1])
    if m.sum() < 20:
        m = np.ones(tau.size, bool)
    ts, ds = tau[m], d[m]

    t0g = t_on + T0_GRUESO
    _, _, c = _ajuste_lineal(ts, ds, grilla, forma, t0g)
    c = np.where(np.isfinite(c), c, np.inf)
    ig = int(np.argmin(c))
    en_borde = ig in (0, t0g.size - 1)

    t0f = t0g[ig] + T0_FINO
    _, _, c = _ajuste_lineal(ts, ds, grilla, forma, t0f)
    c = np.where(np.isfinite(c), c, np.inf)
    i = int(np.argmin(c))
    t0 = float(t0f[i])
    if 0 < i < t0f.size - 1 and np.all(np.isfinite(c[i - 1:i + 2])):
        # vertice de la parabola por tres puntos equiespaciados
        den = c[i - 1] - 2.0 * c[i] + c[i + 1]
        if den > 0:
            h = float(t0f[1] - t0f[0])
            t0 += 0.5 * h * (c[i - 1] - c[i + 1]) / den

    mf = (tau >= t0 + grilla[0]) & (tau <= t0 + grilla[-1])
    A, b, chi2 = _ajuste_lineal(tau[mf], d[mf], grilla, forma, np.array([t0]))
    return dict(A=float(A[0]), b=float(b[0]), t0=t0, chi2=float(chi2[0]),
                n=int(mf.sum()), t_ult=float(tau[-1] - t0), en_borde=bool(en_borde))


def area_evento_plantilla(t, v, plantilla, pre, post, recorte=RECORTE_BORDES,
                          **kw):
    """
    Area por ajuste de plantilla, y en la MISMA pasada el area por compuerta.

    Las dos cosas juntas por dos razones:
      - La BANDERA es la del metodo global, calculada igual. Asi muestra y
        fondo, y el metodo global y este, descartan exactamente los mismos
        eventos: la eficiencia no cambia y el cociente area_plantilla /
        area_compuerta se puede mirar evento a evento, que es el control mas
        directo de que el ajuste no tiene un sesgo con la amplitud.
      - El onset, la linea base y el signo de `area_evento` son la semilla del
        ajuste: no hay que recalcularlos.

    Claves: todas las de `area_evento` (con `area` REEMPLAZADA por la del
    ajuste), mas `area_compuerta`, `A`, `t0`, `b_ajuste`, `chi2_red` (chi2 por
    grado de libertad en unidades del ruido de la linea base: ~1 para un pulso
    limpio, grande para apilamiento o una forma ajena), y `extrapolada` (el
    registro termina antes de que la plantilla converja, o sea que parte de la
    carga la pone la forma y no el dato).
    """
    r = area_evento(t, v, pre=pre, post=post, recorte=recorte, **kw)
    extra = dict(area_compuerta=r["area"], A=np.nan, t0=np.nan, b_ajuste=np.nan,
                 chi2_red=np.nan, extrapolada=False, t0_en_borde=False)
    if not np.isfinite(r["base"]):          # registro inservible (< 16 muestras)
        r.update(extra)
        return r
    tau, vv = preparar(t, v, recorte)
    # el signo con la linea base buena, igual que el paso 3 de `area_evento`
    _, _, _, signo = detectar_pico(tau, vv, r["base"])
    d = signo * vv                      # b queda libre: NO se resta la base
    dd = d - signo * r["base"]
    # onset SUB-muestra para sembrar t0: la plantilla tiene su cruce del 10%
    # en cero, asi que t0 ~ onset interpolado.
    t_on = cruce_onset(tau, dd, int(np.argmax(dd)))
    aj = ajustar_plantilla(tau, d, plantilla, t_on)
    r["area_compuerta"] = r["area"]
    r["area"] = aj["A"] * plantilla["integral"]
    r["A"] = aj["A"]
    r["t0"] = aj["t0"]
    r["b_ajuste"] = aj["b"] * signo
    ndf = max(aj["n"] - 3, 1)
    ruido = r["ruido"] if np.isfinite(r["ruido"]) and r["ruido"] > 0 else 1.0
    r["chi2_red"] = aj["chi2"] / ndf / ruido ** 2
    r["extrapolada"] = bool(aj["t_ult"] < plantilla["t_carga"])
    r["t0_en_borde"] = aj["en_borde"]
    return r


def plantilla_pulso(rutas_archivos, n=1500, amp_minima=AMP_PLANTILLA,
                    t_rango=T_PLANTILLA, paso=PASO_PLANTILLA, iteraciones=1,
                    guarda_base=4.0, canal=None, t_carga=None):
    """
    La FORMA del pulso: promedio alineado por onset, normalizado a amplitud 1.

    rutas_archivos : archivos de donde se mide. Tienen que ser de registro
        LARGO (los XML del detector A: ~130-158 ns tras el onset), porque la
        plantilla tiene que llegar hasta donde el pulso converge. De un binario
        de septiembre saldria una forma cortada a 47-81 ns, y extrapolar con
        ella seria extrapolar con nada.

    Dos diferencias con `pulso_promedio`, las dos por lo mismo (una plantilla
    se usa para AJUSTAR, asi que su ancho se hereda en cada evento):
      - Alineacion SUB-muestra (`cruce_onset`) en la primera pasada, y en las
        `iteraciones` siguientes por el t0 AJUSTADO contra la plantilla previa,
        que usa todo el flanco y no un solo cruce. Cada pulso entra
        normalizado por su A ajustada, que es menos ruidosa que el maximo de
        una muestra.
      - Promedio con CONTEO por punto de la grilla: un registro que no cubre
        un tramo no aporta ahi, en vez de aportar ceros (`pulso_promedio`
        rellena con ceros, y eso hace "converger" la forma donde termina el
        registro mas corto, no donde termina el pulso).

    Devuelve un dict: tau, forma (maximo = 1), n_por_punto, integral (carga
    hasta el plateau de la integral acumulada, que es la CARGA TOTAL del pulso;
    el rebote de polaridad opuesta que viene despues no se cuenta, igual que
    en `curva_convergencia`), integral_total (hasta el final de la grilla),
    t_plateau, frac (curva de convergencia), n_pulsos, archivos y hash.
    """
    import hashlib
    from analisis import drs

    # t_carga="fin": hasta donde llegue la grilla medida. Solo para DIAGNOSTICO
    # (medir la forma de una corrida de registro corto y compararla); una
    # plantilla para ajustar tiene que llegar a T_CARGA.
    hasta_el_fin = isinstance(t_carga, str) and t_carga == "fin"
    t_carga = T_CARGA if (t_carga is None or hasta_el_fin) else float(t_carga)
    grilla = np.arange(t_rango[0], t_rango[1] + 0.5 * paso, paso)
    eventos = []                   # (tau relativo al onset, d normalizada)
    for ruta in rutas_archivos:
        T, V, _, _ = drs.leer_pulsos(ruta, n=n, canal=canal)
        for t, v in zip(T, V):
            tau, vv = preparar(t, v)
            ip, _, amp, signo = detectar_pico(tau, vv)
            if amp < amp_minima or np.max(np.abs(vv)) >= V_SATURA:
                continue
            _, t_on0 = detectar_onset(tau, vv, float(np.median(vv)), ip, signo)
            base, ruido, nb = linea_base(tau, vv, hasta=t_on0 - guarda_base)
            if nb < 10:
                continue
            d = signo * (vv - base)
            ip = int(np.argmax(d))
            if d[ip] < amp_minima:
                continue
            t_on = cruce_onset(tau, d, ip)
            eventos.append((tau, d, t_on, float(d[ip]), signo * vv))
    if not eventos:
        raise ValueError("ningun pulso de mas de %.0f mV en %s"
                         % (amp_minima, rutas_archivos))

    def promediar(desplazamientos, escalas):
        # No es un promedio de pulsos normalizados sino una REGRESION punto a
        # punto de la grilla: d(t) = alfa(t)*A + beta(t) sobre todos los
        # eventos. Medido el 24-sep-2026: la forma normalizada NO es la misma
        # en todos los tamanos (la cola 40-75 ns pesa 7.9% de Q40 en el cuartil
        # mas chico del Cs-137 del 26-ago y 4.6% en el mas grande), y ese
        # patron es el de una componente ADITIVA -algo que esta en cada
        # registro con el mismo tamano, tipicamente interferencia sincronizada
        # con el trigger- que en la forma normalizada pesa como 1/A. Promediar
        # d/A la mete en la plantilla con el peso de la amplitud tipica; la
        # regresion la separa: alfa es la parte que escala con el pulso (la
        # que se ajusta) y beta la que no (se reporta).
        S1 = np.zeros_like(grilla); SA = np.zeros_like(grilla)
        SAA = np.zeros_like(grilla); Sd = np.zeros_like(grilla)
        SAd = np.zeros_like(grilla)
        for (tau, d, _, _, _), t0, a in zip(eventos, desplazamientos, escalas):
            x = tau - t0
            cubre = (grilla >= x[0]) & (grilla <= x[-1])
            y = np.interp(grilla[cubre], x, d)
            S1[cubre] += 1; SA[cubre] += a; SAA[cubre] += a * a
            Sd[cubre] += y; SAd[cubre] += a * y
        det = S1 * SAA - SA * SA
        ok = (S1 > 2) & (det > 0)
        alfa = np.where(ok, (S1 * SAd - SA * Sd) / np.where(ok, det, 1), 0.0)
        beta = np.where(ok, (Sd - alfa * SA) / np.maximum(S1, 1), 0.0)
        promediar.beta = beta
        return alfa, S1

    t0s = [e[2] for e in eventos]
    amps = [e[3] for e in eventos]
    forma, cuenta = promediar(t0s, amps)
    beta = promediar.beta
    for _ in range(int(iteraciones)):
        # la forma, re-anclada: maximo = 1 y cruce del 10% en cero
        forma, beta = _reanclar(grilla, forma / forma.max(), beta)
        pl = dict(tau=grilla, forma=forma)
        nuevos, t0s, amps = [], [], []
        for tau, _, t_on, amp_pico, cruda in eventos:
            aj = ajustar_plantilla(tau, cruda, pl, t_on)
            # Solo entran los ajustes SANOS. La regresion es por minimos
            # cuadrados, y un solo evento con A absurda (un apilamiento que el
            # ajuste lee como un pulso enorme) tiene palanca para arrastrar la
            # forma entera: medido, con 800 eventos del C a 1600 V la forma
            # iterada terminaba con el maximo en 0.35 ns y un pozo de -0.7.
            if aj["en_borde"] or not (0.7 < aj["A"] / amp_pico < 1.4):
                continue
            a = aj["A"]
            # la linea base del AJUSTE reemplaza a la mediana pre-onset: usa
            # todo el registro y no depende de donde se corte el pre-pulso
            nuevos.append((tau, cruda - aj["b"], t_on, a, cruda))
            t0s.append(aj["t0"])
            amps.append(a)
        eventos = nuevos
        forma, cuenta = promediar(t0s, amps)
        beta = promediar.beta
    escala = forma.max()
    forma, beta = _reanclar(grilla, forma / escala, beta)

    # Donde hay pocos registros que lleguen, el promedio es ruido: se corta la
    # grilla donde la cobertura cae por debajo de un cuarto de los pulsos.
    ok = cuenta >= 0.25 * len(eventos)
    ultimo = int(np.flatnonzero(ok)[-1]) if ok.any() else grilla.size - 1
    grilla, forma, cuenta, beta = (grilla[:ultimo + 1], forma[:ultimo + 1],
                                   cuenta[:ultimo + 1], beta[:ultimo + 1])

    _, frac, info = curva_convergencia(grilla, forma)
    Q = integral_acumulada(grilla, forma)
    Qb = integral_acumulada(grilla, beta)
    i_pl = int(np.argmax(Q))
    if hasta_el_fin:
        t_carga = float(grilla[-1])
    if t_carga > grilla[-1] + 1e-9:
        raise ValueError("la plantilla llega a %.1f ns y se pidio la carga hasta "
                         "%.1f ns: hacen falta registros mas largos"
                         % (grilla[-1], t_carga))
    i_c = int(np.searchsorted(grilla, t_carga))
    h = hashlib.sha1(np.round(forma, 6).tobytes() + grilla.tobytes()
                     + ("%.3f" % t_carga).encode()).hexdigest()[:12]
    return dict(tau=grilla, forma=forma, n_por_punto=cuenta,
                integral=float(Q[i_c]), t_carga=float(t_carga),
                integral_plateau=float(Q[i_pl]), integral_total=float(Q[-1]),
                t_plateau=float(grilla[i_pl]), frac=frac,
                caida_final=info["caida_final"],
                convergio=not info["en_borde"],
                aditiva=beta, integral_aditiva=float(Qb[i_c]),
                n_pulsos=len(eventos),
                archivos=[os.path.basename(r) for r in rutas_archivos], hash=h)


def _reanclar(grilla, forma, beta, frac=FRAC_ONSET):
    """Corre la forma (y la componente aditiva, con ella) para que el cruce del
    10% quede exactamente en t = 0: es la convencion con la que
    `ajustar_plantilla` siembra t0 con el onset."""
    ip = int(np.argmax(forma))
    t_c = cruce_onset(grilla, forma, ip, frac)
    return (np.interp(grilla + t_c, grilla, forma, left=0.0, right=forma[-1]),
            np.interp(grilla + t_c, grilla, beta, left=0.0, right=beta[-1]))


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _tabla_convergencia(tau_rel, frac, cortes):
    filas = []
    for c in cortes:
        i = int(np.searchsorted(tau_rel, c))
        i = min(i, frac.size - 1)
        filas.append((tau_rel[i], frac[i]))
    return filas


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES if rutas.corrida(k).disponible()]
    if not claves:
        print("No hay corridas accesibles.")
        return

    n_leer = 600
    fig, axes = plt.subplots(len(claves), 4, figsize=(20, 4.4 * len(claves)),
                             squeeze=False)
    resumen = {}

    for fila, clave in enumerate(claves):
        c = rutas.corrida(clave)
        print("=" * 76)
        print("%s  —  %s" % (clave, c.descripcion))
        print("=" * 76)
        T, V, _, _ = drs.leer_pulsos(c.muestra, n=n_leer)

        # --- forma promedio alineada por onset ------------------------------
        tau_rel, forma, n_us = pulso_promedio(T, V)
        _, frac, info = curva_convergencia(tau_rel, forma)
        print("  pulso promedio: %d de %d eventos (los de amplitud > 20 mV)"
              % (n_us, len(T)))
        print("  la integral hace techo en t - onset = %.1f ns; despues cae "
              "%.2f%% (rebote de polaridad opuesta en la cola)"
              % (info["t_plateau"], 100 * info["caida_final"]))

        # ruido del promedio, medido en su propia zona pre-onset
        pre_zona = tau_rel < -8.0
        ruido_f = float(1.4826 * np.median(np.abs(forma[pre_zona]
                                                  - np.median(forma[pre_zona]))))
        t_ini_f, t_fin_f, diag = ventana_por_convergencia(
            tau_rel, forma, 0.0, 1.0, ruido_f, k=3.0, tramo=4.0, ancho_suave=1.0)
        print("  ventana sobre la forma promedio: [%.1f, %.1f] ns respecto del "
              "onset  (ancho %.1f ns)" % (t_ini_f, t_fin_f, t_fin_f - t_ini_f))
        i_fin_f = int(np.searchsorted(tau_rel, t_fin_f))
        adentro = float(frac[i_fin_f])
        print("  carga dentro de esa ventana: %.2f%% de la que llega a acumularse"
              % (100 * adentro))
        if adentro < 0.99:
            print("  ATENCION: queda %.1f%% de la carga FUERA de la ventana. La "
                  "cola de este pulso sigue aportando mas alla de t_fin, o sea "
                  "que la integral no llego a converger de verdad" % (100 * (1 - adentro)))
            print("     dentro de lo que se puede mirar (el registro dura %.0f ns "
                  "y el pulso arranca a mitad de camino). Eso NO invalida la "
                  "calibracion mientras la fraccion perdida sea la misma para "
                  "todos los pulsos" % tau_rel[-1])
            print("     — que es exactamente lo que verifica el test de "
                  "area/amplitud de la Etapa 2 —, pero hay que tenerlo anotado.")

        print("  curva de convergencia (carga acumulada desde el onset):")
        for tt, ff in _tabla_convergencia(tau_rel, frac,
                                          [-10, 0, 5, 10, 20, 30, 40, 50, 60,
                                           80, 100, 120]):
            print("      t - onset = %6.1f ns   ->  %7.2f %%" % (tt, 100 * ff))

        pre = max(-t_ini_f, 2.0)
        post = t_fin_f
        resumen[clave] = dict(pre=pre, post=post, tau_rel=tau_rel, forma=forma,
                              frac=frac, t_ini=t_ini_f, t_fin=t_fin_f)

        # --- comparacion de lineas base: la correcta y la del error viejo ----
        # Sobre TODO el ensemble, no sobre un evento: en un solo pulso las dos
        # ventanas pueden coincidir por casualidad, y lo que importa es la
        # DISPERSION que la ventana fija le agrega al espectro.
        d_base, d_media = [], []
        t_on0 = np.nan
        for t, v in zip(T, V):
            tt, vv = preparar(t, v)
            i_p0, _, a0, sg0 = detectar_pico(tt, vv)
            if a0 < 20.0:
                continue
            _, t_on0 = detectar_onset(tt, vv, float(np.median(vv)), i_p0, sg0)
            b_ok, _, n_ok = linea_base(tt, vv, hasta=t_on0 - 4.0)
            if n_ok < 10:
                continue
            b_mal, _, _ = linea_base(tt, vv, hasta=0.20 * (tt[-1] - tt[0]))
            b_med, _, _ = linea_base(tt, vv, hasta=t_on0 - 4.0, estadistico="media")
            d_base.append(sg0 * (b_mal - b_ok))
            d_media.append(b_med - b_ok)
        d_base = np.array(d_base); d_media = np.array(d_media)
        ancho_g = pre + post
        print("  LINEA BASE (%d eventos):" % d_base.size)
        print("     media vs mediana en el intervalo pre-onset: desvio %.4f +- "
              "%.4f mV  ->  irrelevante" % (d_media.mean(), d_media.std()))
        print("     ventana FIJA del 20% del registro (el error viejo) menos la "
              "pre-onset:")
        print("        %.3f +- %.3f mV  ->  sobre la compuerta de %.0f ns son "
              "%.0f +- %.0f mV*ns de sesgo"
              % (d_base.mean(), d_base.std(), ancho_g,
                 d_base.mean() * ancho_g, d_base.std() * ancho_g))
        print("        (lo grave no es el promedio sino la DISPERSION: es ancho "
              "que se le suma al espectro evento por evento)")

        # --- que tan estable es el criterio evento a evento -----------------
        anchos, fracs, amps = [], [], []
        for t, v in zip(T[:400], V[:400]):
            r = area_evento(t, v)                    # modo por convergencia
            if r["bandera"] & (MAL_SIN_PICO | MAL_TRUNCADO):
                continue
            anchos.append(r["t_fin"] - r["t_ini"])
            fracs.append(r["fraccion"])
            amps.append(r["amplitud"])
        anchos = np.array(anchos); amps = np.array(amps)
        print("  VENTANA QUE ELIGE CADA EVENTO POR SI SOLO (%d eventos):" % anchos.size)
        print("     mediana %.1f ns, rango [%.1f, %.1f], desviacion %.1f ns"
              % (np.median(anchos), anchos.min(), anchos.max(), anchos.std()))
        cortes = np.percentile(amps, [0, 25, 50, 75, 100])
        print("     por cuartil de amplitud:")
        for j in range(4):
            m = (amps >= cortes[j]) & (amps <= cortes[j + 1])
            print("        amp %5.0f-%5.0f mV  ->  ancho mediano %5.1f ns "
                  "(+- %.1f)" % (cortes[j], cortes[j + 1],
                                 np.median(anchos[m]), anchos[m].std()))
        print("     ^ el ancho CRECE con la amplitud: el criterio corta donde la "
              "senal se hunde en el ruido DE ESE evento, y en un pulso chico eso "
              "pasa antes. Es un corte que depende de la amplitud, o sea que")
        print("       deforma el espectro que se quiere calibrar. La Etapa 2 mide "
              "cuanto.")
        print()

        # ---------------------------------------------------------- figuras --
        ax = axes[fila][0]
        onsets = []
        for k in range(min(25, len(V))):
            tt, vv = preparar(T[k], V[k])
            ax.plot(tt, vv, lw=0.6, alpha=0.45)
            i_pk, _, a_k, sg_k = detectar_pico(tt, vv)
            if a_k >= 20.0:
                onsets.append(detectar_onset(tt, vv, float(np.median(vv)),
                                             i_pk, sg_k)[1])
        t_on_med = float(np.median(onsets)) if onsets else 20.0
        largo = float(tt[-1] - tt[0])
        ax.axvspan(0, t_on_med - 4.0, color="C2", alpha=0.18,
                   label="linea base pre-onset (onset mediano %.0f ns)" % t_on_med)
        ax.axvspan(0, 0.20 * largo, fc="none", ec="C3", hatch="//",
                   lw=1.2, label="ventana fija del 20% (MAL: se come el pulso)")
        ax.set_title("%s: pulsos crudos y de donde sale el offset" % clave)
        ax.set_xlabel("t desde el inicio del registro [ns]")
        ax.set_ylabel("v [mV]"); ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[fila][1]
        ax.plot(tau_rel, forma, lw=1.4, color="C0")
        ax.axhline(0, color="k", lw=0.6)
        ax.axhspan(-3 * ruido_f, 3 * ruido_f, color="0.7", alpha=0.5,
                   label="banda de +-3 sigma del ruido")
        ax.axvline(t_ini_f, color="C2", ls="--", label="t_ini = %.1f ns" % t_ini_f)
        ax.axvline(t_fin_f, color="C3", ls="--", label="t_fin = %.1f ns" % t_fin_f)
        ax.set_title("pulso promedio alineado por onset (%d eventos)" % n_us)
        ax.set_xlabel("t - onset [ns]"); ax.set_ylabel("forma normalizada")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[fila][2]
        ax.plot(tau_rel, 100 * frac, lw=1.6, color="C0")
        ax.axhline(100, color="k", ls=":", lw=1)
        ax.axhline(99.5, color="0.5", ls=":", lw=1)
        ax.axvline(t_ini_f, color="C2", ls="--")
        ax.axvline(t_fin_f, color="C3", ls="--")
        ax.annotate("t_fin = %.1f ns\n%.2f%% de la carga"
                    % (t_fin_f, 100 * frac[np.searchsorted(tau_rel, t_fin_f)]),
                    xy=(t_fin_f, 100 * frac[np.searchsorted(tau_rel, t_fin_f)]),
                    xytext=(t_fin_f + 12, 60), fontsize=8,
                    arrowprops=dict(arrowstyle="->", lw=0.8))
        ax.set_title("INTEGRAL ACUMULADA: donde converge")
        ax.set_xlabel("t - onset [ns]"); ax.set_ylabel("carga acumulada [%]")
        ax.set_ylim(-8, 118); ax.grid(alpha=0.3)

        ax = axes[fila][3]
        ax.scatter(amps, anchos, s=5, alpha=0.35)
        ax.axhline(post + pre, color="C1", lw=1.5,
                   label="ancho global = %.1f ns" % (post + pre))
        ax.set_title("ancho que elige CADA evento por si solo")
        ax.set_xlabel("amplitud [mV]"); ax.set_ylabel("t_fin - t_ini [ns]")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

    fig.tight_layout()
    salida = rutas.figura("01_ventana_integracion.png")
    fig.savefig(salida, dpi=125)
    print("figura -> %s" % rutas.relativa(salida))
    print()
    print("VENTANA GLOBAL que sale de esta etapa (respecto del onset de cada evento):")
    for clave, r in resumen.items():
        print("   %-12s pre = %.1f ns   post = %.1f ns   (ancho %.1f ns)"
              % (clave, r["pre"], r["post"], r["pre"] + r["post"]))


if __name__ == "__main__":
    _demo()
