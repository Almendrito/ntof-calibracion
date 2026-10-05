"""
Etapa 2 — el observable: el area de cada pulso en una compuerta anclada al onset.

    python desde_cero/area.py

Que resuelve esta etapa
-----------------------
A cada evento se le asigna UN numero, el area del pulso en mV*ns, integrada en

    [onset - PRE, onset + POST]        con el onset de ESE evento (Etapa 1)

Ese numero es el eje x del espectro, y todo lo que viene despues (borde
Compton, D, ley de ganancia) es una afirmacion sobre el. Si el area no es
proporcional a la luz del centellador, ningun ajuste lo arregla.

Las decisiones, y por que
-------------------------
1. **Anclada al onset de cada evento, no a un tiempo fijo del registro.** El
   jitter del trigger mueve el pulso ~35 ns entre eventos (CLAUDE.md §7.33).
   Una compuerta fija cortaria a cada pulso en un lugar distinto de su forma,
   y la fraccion de carga recogida variaria evento a evento: ensanchamiento
   puro del espectro.

2. **Trapecio sobre el eje de tiempo REAL.** El DRS4 no muestrea uniforme
   (Etapa 0: el ancho de celda varia ~5 % y hasta 0.1-0.3 ns entre celdas).
   `sum(v)*dt` meteria un error que cambia con la celda donde cayo el pulso.

3. **Compuerta COMUN, POST = 12 ns, igual en todas las corridas** (§7.53,
   §7.66). El pulso completo necesita ~100 ns para converger, y el registro no
   siempre los tiene: la sesion del 1-sep deja ~81 ns tras el onset y el D del
   8-sep 47 ns. Si cada archivo usara la compuerta que su pulso pide, dos
   corridas recogerian fracciones distintas de la carga y esa diferencia
   entraria derecha en D (la ley de ganancia pasaba de chi2 = 744 a 9.2 con 1
   gdl al igualar la compuerta). 12 ns recoge ~72-82 % de la carga, pero la
   MISMA fraccion en todas, y despues del pico el pulso ya bajo a ~20 % de su
   amplitud: seguir integrando agrega ruido de linea base mas rapido que
   senal. Precio: D es "la ganancia medida en esta compuerta", no la de pulso
   completo. Para comparar corridas da igual; para una ganancia absoluta, no.

4. **PRE = 2 ns.** Antes del cruce del 10 % el pulso casi no tiene carga:
   con la forma media, [-2, 0] ns es el 0.8 % de la carga de la compuerta, y
   estirarlo a 8.8 ns solo sube eso a 0.9 %. `analisis/area.py` MIDE pre en
   cada archivo (2.0-2.3 ns en las corridas de trabajo); aca se fija, porque
   lo que se agrega mas alla es linea base: ruido, no senal.

5. **Banderas que se marcan, no se descartan.** Que hacer con un evento
   saturado o truncado depende de para que se use (un espectro de forma no es
   un conteo absoluto), y descartarlo aca esconderia un corte dentro del
   calculo del area. Descartar es una decision de la Etapa 3, a la vista.

6. **Tramas repetidas SI se tiran.** El DRSOSC a veces vuelca dos veces el
   mismo registro. No es un evento con un defecto: es el mismo evento contado
   dos veces, y dejarlo meteria un pico falso en el espectro.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import lector, pulso                        # noqa: E402

_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

PRE = 2.0               # [ns] antes del onset
POST = 12.0             # [ns] despues del onset: la compuerta comun
UMBRAL_PP = 3.0         # [mV] pico a pico por debajo de esto no hay pulso
T_COMUN = 40.0          # [ns] tras el onset: referencia de "carga" (ver abajo)

# Banderas: las tres de la Etapa 1 mas una. Mismos valores que analisis/.
SATURADO = pulso.SATURADO
TRUNCADO = 4            # la compuerta se sale del registro
SIN_PRE = pulso.SIN_PRE
SIN_PULSO = pulso.SIN_PULSO


# --------------------------------------------------------------------------
# 1. El area de un evento
# --------------------------------------------------------------------------
def area_evento(t, v, pre=PRE, post=POST):
    """
    Area de un pulso en `[onset - pre, onset + post]`, positiva, en mV*ns.

    Devuelve el dict de `pulso.medir` con tres claves mas: `area`, `t_ini` y
    `t_fin` (los bordes de la integral que SE USARON).

    Los bordes caen sobre muestras, no exactamente en onset-pre y onset+post:
    la integral va de la primera muestra >= onset-pre a la primera muestra
    >= onset+post. Es lo que hace el oraculo, y el error que mete es a lo sumo
    un dt (~0.2 ns) en cada extremo. En el extremo izquierdo no pesa (ahi no
    hay pulso); en el derecho el pulso vale ~22 % de su amplitud, asi que el
    area fluctua +-0.2 ns * 0.22 / 7.1 ns ~ 0.6 % pico a pico de evento a
    evento. Como no depende de la amplitud, ensancha un poco pero no deforma.
    `verificar_area.py` lo mide sobre pulsos fabricados.
    """
    m = pulso.medir(t, v)
    tau, vv = m["tau"], m["v"]
    t_ini, t_fin = m["t_onset"] - pre, m["t_onset"] + post
    bandera = m["bandera"]
    if t_ini < tau[0] or t_fin > tau[-1]:
        # El registro termina antes de que termine la compuerta: el area que
        # sale es menor que la verdadera, en una fraccion que depende de
        # CUANDO llego el pulso. Se marca; la Etapa 3 decide.
        bandera |= TRUNCADO
    i_ini = int(np.clip(np.searchsorted(tau, t_ini), 0, tau.size - 2))
    i_fin = int(np.clip(np.searchsorted(tau, t_fin), i_ini + 2, tau.size - 1))
    # el signo deja el area positiva sea cual sea la polaridad del pulso
    a = m["signo"] * _TRAPZ(vv[i_ini:i_fin + 1] - m["base"],
                            x=tau[i_ini:i_fin + 1])
    m.update(area=float(a), t_ini=float(tau[i_ini]), t_fin=float(tau[i_fin]),
             bandera=bandera)
    return m


# --------------------------------------------------------------------------
# 2. Un archivo entero
# --------------------------------------------------------------------------
CAMPOS = ("area", "amplitud", "t_onset", "base", "ruido", "pico_a_pico")


def areas_de_archivo(ruta, pre=PRE, post=POST, limite=None):
    """
    Recorre el archivo y devuelve un dict de arrays, un elemento por evento.

    Claves: las de CAMPOS, `bandera` (int), `indice` (numero de evento en el
    archivo, para poder volver a mirar uno), y los escalares `n_leidos` y
    `repetidos`.

    Ademas `reloj` (str) y `scaler` [Hz] de TODOS los eventos leidos, tramas
    repetidas incluidas: son para el tiempo vivo de la Etapa 3, que cuenta
    cuantos eventos GUARDO el DRSOSC contra cuantos triggers disparo la placa,
    y eso no depende de si despues el evento sirvio o no.

    Los eventos con pico a pico <= UMBRAL_PP se CONSERVAN con la bandera
    SIN_PULSO: no hay pulso y su area es ruido alrededor de cero, pero tirarlos
    aca haria que la Etapa 3 no pudiera ni contarlos. (El oraculo si los tira;
    `verificar_area.py` alinea las dos listas por `pico_a_pico`.)
    """
    acum = {c: [] for c in CAMPOS}
    banderas, indices = [], []
    relojes, scalers = [], []
    anterior, n, rep = None, 0, 0
    for k, (t, v, reloj, sc) in enumerate(lector.eventos(ruta, limite=limite,
                                                          con_reloj=True)):
        n += 1
        relojes.append(reloj)
        scalers.append(sc)
        if anterior is not None and np.array_equal(v, anterior):
            rep += 1
            continue
        anterior = v
        r = area_evento(t, v, pre, post)
        # sobre el registro CRUDO, como el oraculo: es un filtro previo que
        # pregunta "hay algo aca", antes de cualquier recorte
        r["pico_a_pico"] = float(np.max(v) - np.min(v))
        b = r["bandera"]
        if r["pico_a_pico"] <= UMBRAL_PP:
            b |= SIN_PULSO
        for c in CAMPOS:
            acum[c].append(r[c])
        banderas.append(b)
        indices.append(k)
    out = {c: np.asarray(acum[c], float) for c in CAMPOS}
    out.update(bandera=np.asarray(banderas, int), indice=np.asarray(indices),
               reloj=relojes, scaler=np.asarray(scalers, float),
               n_leidos=n, repetidos=rep, pre=pre, post=post)
    if out["area"].size == 0:
        # un vacio silencioso aparece cinco etapas despues como un histograma
        # de ceros (§7.36): mejor explotar aca, donde esta la causa
        raise ValueError("%s: cero eventos leidos" % os.path.basename(ruta))
    return out


# --------------------------------------------------------------------------
# 3. El pulso promedio y su carga acumulada
# --------------------------------------------------------------------------
def pulso_promedio(ruta, n=1500, amp_min=20.0, tau_min=-20.0, tau_max=140.0,
                   paso=0.05):
    """
    Forma media del pulso, alineada por onset y normalizada a amplitud 1.

    Devuelve `(grilla, forma, cuantos)`: `cuantos[i]` es el numero de pulsos
    que tenian registro en `grilla[i]`.

    * **Alinear por onset** antes de promediar: promediar los registros tal
      cual daria el pulso convolucionado con el jitter del trigger, mas ancho
      de lo que es.
    * **Onset interpolado** (`t_onset_fino`), no el de muestra: alinear por
      muestra mete un jitter uniforme de un dt (~0.2 ns) que ensancha la
      subida del promedio.
    * **Normalizado por amplitud**, para que sea la FORMA y no la dominen los
      pulsos grandes. `amp_min` deja afuera los chicos, cuyo onset es ruidoso
      (Etapa 1): es solo para medir la forma.
    * **Se promedia solo donde hay registro.** Rellenar con ceros fuera del
      registro haria que en un archivo de registro corto la forma "llegue a
      cero" donde en realidad se acabo el dato, y la carga acumulada parezca
      converger antes (la trampa de §7.66 con el detector D).
    """
    grilla = np.arange(tau_min, tau_max + 0.5 * paso, paso)
    suma = np.zeros_like(grilla)
    cuantos = np.zeros_like(grilla)
    usados = 0
    for t, v in lector.eventos(ruta, limite=n):
        m = pulso.medir(t, v)
        if m["amplitud"] < amp_min or m["bandera"]:
            continue
        rel = m["tau"] - m["t_onset_fino"]
        d = m["signo"] * (m["v"] - m["base"]) / m["amplitud"]
        dentro = (grilla >= rel[0]) & (grilla <= rel[-1])
        suma[dentro] += np.interp(grilla[dentro], rel, d)
        cuantos[dentro] += 1
        usados += 1
    if usados == 0:
        raise ValueError("ningun pulso de mas de %.0f mV" % amp_min)
    with np.errstate(invalid="ignore"):
        forma = np.where(cuantos > 0, suma / np.maximum(cuantos, 1), np.nan)
    return grilla, forma, cuantos


def carga_acumulada(grilla, forma, desde=-PRE):
    """
    Q(tau) = integral de la forma desde `desde` hasta tau. NaN donde no hay
    registro (ahi no se sabe cuanta carga hay).
    """
    ok = np.isfinite(forma)
    Q = np.full_like(grilla, np.nan)
    m = ok & (grilla >= desde)
    g, f = grilla[m], forma[m]
    Q[m] = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(g))))
    return Q


def fraccion_de_carga(grilla, forma, pre=PRE, post=POST, hasta=T_COMUN):
    """
    Que fraccion de la carga recoge la compuerta, relativa a la carga hasta
    `hasta` ns tras el onset.

    Por que relativa a 40 ns y no "a la carga total": el total no existe en
    todos los registros (el D del 8-sep termina a 47 ns del onset, y en junio
    el pulso sigue vivo a 130 ns). Una fraccion normalizada al final de CADA
    registro compararia registros, no pulsos (§7.66). 40 ns cabe en todos.
    """
    Q = carga_acumulada(grilla, forma, -pre)
    q_gate = np.interp(post, grilla, Q)
    q_ref = np.interp(hasta, grilla, Q)
    return float(q_gate / q_ref)


# --------------------------------------------------------------------------
# 4. El control que decide: area/amplitud contra amplitud
# --------------------------------------------------------------------------
def perfil_area_amplitud(d, bins=12, minimo=30):
    """
    Mediana de area/amplitud en bines de amplitud (cuantiles), eventos sin
    bandera. Devuelve `(amp_centro, mediana, error)`.

    Si todos los pulsos tienen la misma forma, area/amplitud es la misma
    constante para todos: la integral de la forma sobre la compuerta. Si se
    inclina con la amplitud, la compuerta recoge una fraccion distinta segun
    el tamano del pulso, y entonces el area deja de ser proporcional a la luz:
    el eje del espectro se curva y ninguna calibracion lo endereza. No
    necesita conocer ni la fuente ni la energia.
    """
    ok = d["bandera"] == 0
    a, r = d["amplitud"][ok], d["area"][ok] / d["amplitud"][ok]
    bordes = np.quantile(a, np.linspace(0, 1, bins + 1))
    c, med, err = [], [], []
    for lo, hi in zip(bordes[:-1], bordes[1:]):
        m = (a >= lo) & (a < hi)
        if m.sum() < minimo:
            continue
        c.append(np.median(a[m]))
        med.append(np.median(r[m]))
        # error de la mediana: 1.2533 * sigma / sqrt(n) para una gaussiana
        err.append(1.2533 * pulso.MAD_A_SIGMA
                   * np.median(np.abs(r[m] - np.median(r[m]))) / np.sqrt(m.sum()))
    return np.array(c), np.array(med), np.array(err)


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def _demo(n=5000):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ejes = plt.subplots(1, 3, figsize=(16, 4.6))
    for clave, (ruta, _) in lector.CORRIDAS_TRABAJO.items():
        if not os.path.exists(ruta):
            print("  %s: no esta en disco" % clave)
            continue
        g, f, cu = pulso_promedio(ruta)
        Q = carga_acumulada(g, f)
        q40 = np.interp(T_COMUN, g, Q)
        fr = fraccion_de_carga(g, f)
        fin = g[np.flatnonzero(cu > 0.5 * cu.max())[-1]]
        ejes[0].plot(g, f, lw=1, label="%s (registro hasta %.0f ns)" % (clave, fin))
        ejes[1].plot(g, Q / q40, lw=1, label="%s: compuerta = %.1f %%"
                     % (clave, 100 * fr))

        d = areas_de_archivo(ruta, limite=n)
        c, med, err = perfil_area_amplitud(d)
        sobre = c >= np.median(d["amplitud"][d["bandera"] == 0])
        plano = (med[sobre].max() - med[sobre].min()) / np.median(med[sobre])
        ejes[2].errorbar(c, med, err, fmt="o-", ms=3, lw=1, capsize=2,
                         label="%s: variacion %.1f %% sobre la mediana"
                         % (clave, 100 * plano))
        b = d["bandera"]
        print("  %-12s %5d ev (%d repetidos)  compuerta [-%.0f, +%.0f] ns = "
              "%.1f %% de la carga a 40 ns   area/amp sobre la mediana: "
              "variacion %.1f %%   banderas: sat %d, trunc %d, sin pre %d, "
              "sin pulso %d"
              % (clave, d["area"].size, d["repetidos"], PRE, POST, 100 * fr,
                 100 * plano, (b & SATURADO > 0).sum(), (b & TRUNCADO > 0).sum(),
                 (b & SIN_PRE > 0).sum(), (b & SIN_PULSO > 0).sum()))

    ax = ejes[0]
    ax.axvspan(-PRE, POST, color="0.85", label="compuerta")
    ax.set(xlabel="tiempo desde el onset [ns]", ylabel="forma / amplitud",
           title="pulso promedio alineado por onset", xlim=(-15, 140))
    ax = ejes[1]
    ax.axvspan(-PRE, POST, color="0.85")
    ax.axvline(T_COMUN, color="0.5", ls=":")
    ax.set(xlabel="tiempo desde el onset [ns]",
           ylabel="carga acumulada / carga a 40 ns",
           title="carga acumulada", xlim=(-15, 140), ylim=(-0.05, 1.25))
    ax = ejes[2]
    ax.set_xscale("log")
    ax.set(xlabel="amplitud [mV]", ylabel="area / amplitud [ns]",
           title="control de linealidad (eventos sin bandera)")
    for ax in ejes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa2_area.png")
    fig.savefig(salida, dpi=120)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
