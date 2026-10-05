"""
Etapa 3 — de las areas al espectro neto: histograma, tiempo vivo, resta, error.

    python desde_cero/espectro.py

Que resuelve esta etapa
-----------------------
Con el area de cada evento (Etapa 2) de la MUESTRA (fuente puesta) y del
FONDO (sin fuente), arma el espectro de lo que hizo la fuente sola:

    neto(bin) = k_m * N_m(bin) - k_f * N_f(bin)
    error(bin) = sqrt(k_m^2 N_m + k_f^2 N_f)        (piso 1)

con `k_m`, `k_f` que llevan los dos conteos al MISMO tiempo vivo. Eso es lo
que la Etapa 6 (NCC) y la 7 (ajuste) comparan contra el modelo.

Las decisiones, y por que
-------------------------
1. **Los mismos bordes de bin para los dos.** Si muestra y fondo se
   histogramean cada uno con su rango, el bin 37 de uno no es el bin 37 del
   otro y la resta no significa nada. Un solo `rango` para los dos, siempre.

2. **Tope por percentil, no por maximo.** Unos pocos eventos de area enorme
   (cosmicos, apilamiento) estiran el maximo a 20-40 veces donde vive la
   fisica; con 500 bines el borde Compton quedaria en 20 bines (§7.7). El tope
   es 1.6 x el percentil 99.5 de la muestra: el borde queda con aire arriba
   (hace falta ver la caida entera, §7.6) y la cola rara no manda.

3. **Tiempo VIVO desde el scaler, no el reloj de pared.** Ver `tiempo_vivo`.

4. **Poisson propagado y NO recortar negativos.** Donde el fondo le gana a la
   muestra, el neto puede ser negativo: es una fluctuacion legitima. Recortar
   a cero sesgaria el ajuste hacia arriba en cada bin de pocas cuentas. El
   piso 1 en la varianza es para que un bin con cero cuentas no tenga error
   cero (y peso infinito en el chi2).

5. **El mismo corte de calidad para muestra y fondo.** Se usan los eventos sin
   bandera, salvo SIN_PRE (base tomada de la cola: el evento es bueno, y
   tirarlo seria un corte que depende de la amplitud por el walk del
   trigger, Etapa 1). Si cada histograma tuviera su propio corte, tendrian
   eficiencias distintas y la resta dejaria de ser una resta.
"""

import os
import sys
from datetime import datetime

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import area                                 # noqa: E402

BINES = 500
TOPE_FACTOR = 1.6
TOPE_PERCENTIL = 99.5
RECORTE_SCALER = 10.0      # veces la mediana
PERMITIDAS = area.SIN_PRE  # banderas que NO descartan
SUAVE_UMBRAL = 9           # bines del promedio movil de `fin_del_umbral`


# --------------------------------------------------------------------------
# 1. Tiempo vivo
# --------------------------------------------------------------------------
def a_segundos(relojes):
    """Marcas "AAAA/MM/DD hh:mm:ss.mmm" -> segundos desde la primera."""
    fmt = "%Y/%m/%d %H:%M:%S.%f"
    t0 = datetime.strptime(relojes[0], fmt)
    return np.array([(datetime.strptime(r, fmt) - t0).total_seconds()
                     for r in relojes])


def tiempo_vivo(relojes, scaler, n_eventos, recorte=RECORTE_SCALER):
    """
    Devuelve `(t_vivo, fraccion_viva, t_pared, scaler_mediano)` en segundos.

    Por que el reloj de pared no sirve: el DRS4 dispara un trigger por cada
    pulso que pasa el umbral, pero el programa que guarda (DRSOSC) tarda en
    escribir cada evento (~35 kB de XML) y mientras escribe la placa no esta
    escuchando. Lo guardado es una muestra de lo disparado. Si la tasa real es
    R y se guardaron N eventos en T de pared, la placa estuvo "viva" una
    fraccion N / (R T) del tiempo, y SOLO en ese tiempo vivo se contaron
    los eventos del histograma. Comparar dos corridas por tiempo de pared es
    suponer que las dos estuvieron vivas igual. No lo estuvieron: el Co-60 de
    junio guardo 71.5 ev/s con 170 Hz reales (40 % vivo) y su fondo 10 ev/s con
    10 Hz (100 %). Por pared, el fondo entraba con k = 0.709 en vez de 0.287:
    2.5 veces sobre-restado (§7.34).

    Como se cuenta R*T: el `scaler` de cada evento es la tasa de trigger que
    la placa midio en ese momento [Hz]. Los disparos verdaderos son la
    integral de esa tasa en el tiempo, `sum(scaler_i * dt_i)`, con `dt_i` el
    tiempo desde el evento anterior. Las marcas de tiempo vienen en rafagas
    (miles de eventos con el mismo milisegundo), asi que dt es cero casi
    siempre; da igual, porque la SUMA de los dt es la duracion entera.

    `recorte`: el scaler tiene rafagas de ruido electronico (hasta 90 kHz en
    un fondo que promedia 10 Hz). Se recorta a `recorte` veces su mediana
    antes de integrar; si no, una rafaga sola se come el resultado.

    Sin scaler (nan o cero) se devuelve la duracion de pared con 100 % vivo:
    sin dato, sub-restar el fondo es preferible a sobre-restarlo, porque lo
    segundo puede fabricar un neto negativo que el ojo no detecta.
    """
    t = a_segundos(relojes)
    pared = float(t[-1] - t[0])
    s = np.asarray(scaler, float)
    ok = np.isfinite(s) & (s > 0)
    if t.size < 2 or ok.sum() < 10 or pared <= 0:
        return pared, 1.0, pared, np.nan
    med = float(np.median(s[ok]))
    s = np.clip(np.nan_to_num(s, nan=med), 0.0, recorte * max(med, 1.0))
    dt = np.clip(np.diff(t, prepend=t[0]), 0.0, None)
    disparos = float(np.sum(s * dt))
    if disparos <= 0:
        return pared, 1.0, pared, med
    # >1 solo puede salir por la cuantizacion del scaler (Hz enteros a tasa
    # baja); fisicamente no se puede guardar mas de lo que se disparo
    frac = min(n_eventos / disparos, 1.0)
    return frac * pared, frac, pared, med


# --------------------------------------------------------------------------
# 2. El espectro neto
# --------------------------------------------------------------------------
def buenos(d, permitidas=PERMITIDAS):
    """Mascara de eventos sin banderas, salvo las `permitidas`."""
    return (d["bandera"] & ~permitidas) == 0


def espectro_neto(d_m, d_f=None, bines=BINES, tope=None):
    """
    Espectro neto a partir de los dicts de `area.areas_de_archivo`.

    Devuelve un dict: `canales` (centros), `bordes`, `muestra` y `fondo` ya
    escalados, `neto`, `error`, `k_m`, `k_f`, y los tiempos.

    Sin fondo NO se inventa una resta: neto = muestra cruda, k = 1, y
    `hay_fondo = False`. Prestarle el fondo de otra corrida seria peor que no
    restar: el fondo depende del punto de operacion (a otro HV el mismo
    cosmico deja otra area), asi que meteria una forma ajena (§7.57).
    """
    a_m = d_m["area"][buenos(d_m)]
    if tope is None:
        tope = TOPE_FACTOR * float(np.percentile(a_m, TOPE_PERCENTIL))
    h_m, bordes = np.histogram(a_m, bins=bines, range=(0.0, tope))
    canales = 0.5 * (bordes[:-1] + bordes[1:])
    tv_m, fr_m, pa_m, sc_m = tiempo_vivo(d_m["reloj"], d_m["scaler"],
                                         d_m["n_leidos"])
    r = dict(canales=canales, bordes=bordes, tope=tope, hay_fondo=d_f is not None,
             t_vivo_m=tv_m, frac_viva_m=fr_m, t_pared_m=pa_m, scaler_m=sc_m,
             n_m=int(a_m.size))
    if d_f is None:
        r.update(muestra=h_m.astype(float), fondo=np.zeros(bines), k_m=1.0,
                 k_f=0.0, neto=h_m.astype(float),
                 error=np.sqrt(np.maximum(h_m, 1.0)))
        return r

    a_f = d_f["area"][buenos(d_f)]
    h_f, _ = np.histogram(a_f, bins=bordes)                # MISMOS bordes
    tv_f, fr_f, pa_f, sc_f = tiempo_vivo(d_f["reloj"], d_f["scaler"],
                                         d_f["n_leidos"])
    # Se lleva todo al MENOR de los dos tiempos vivos: asi ningun conteo se
    # AMPLIFICA (k <= 1). Un histograma multiplicado por 3 tendria "cuentas"
    # que nadie conto y un error que ya no es el de Poisson.
    t_norm = min(tv_m, tv_f)
    k_m, k_f = t_norm / tv_m, t_norm / tv_f
    r.update(muestra=k_m * h_m, fondo=k_f * h_f, k_m=k_m, k_f=k_f,
             neto=k_m * h_m - k_f * h_f,
             # var(k N) = k^2 var(N) = k^2 N (Poisson); en una resta de dos
             # conteos independientes las varianzas se SUMAN
             error=np.sqrt(np.maximum(k_m ** 2 * h_m + k_f ** 2 * h_f, 1.0)),
             t_vivo_f=tv_f, frac_viva_f=fr_f, t_pared_f=pa_f, scaler_f=sc_f,
             n_f=int(a_f.size))
    return r


# --------------------------------------------------------------------------
# 3. El piso de cualquier ventana: donde termina la subida del umbral
# --------------------------------------------------------------------------
def fin_del_umbral(canales, neto, suave=SUAVE_UMBRAL):
    """
    Canal donde termina la subida de la eficiencia del trigger (§7.48).

    Debajo del umbral del discriminador el espectro no es fisica sino la
    curva de eficiencia del trigger, que sube de 0 a 1. Ninguna plantilla la
    contiene, asi que una ventana de ajuste que entre ahi no sesga, rompe
    (en el Cs-137 a 1600 V: chi2/ndf 4.3 arrancando en 229 mV*ns, 0.95 en 320).

    Criterio: el primer maximo local del espectro suavizado, empezando donde
    deja de ser cero. La eficiencia CRECE con el area y todo lo que hay debajo
    del borde Compton DECRECE, asi que donde el espectro deja de subir una
    cosa le pasa el mando a la otra. Si el espectro ya baja desde el primer
    bin (el Co-60 de junio, con el discriminador en 0.2 mV) no hay subida y
    se devuelve el primer canal: no restringe nada, que es lo correcto.
    """
    x, y = np.asarray(canales, float), np.asarray(neto, float)
    if y.size < 2 * suave:
        return float(x[0])
    ys = np.convolve(y, np.ones(suave) / suave, mode="same")
    if ys.max() <= 0:
        return float(x[0])
    # por debajo del umbral hay CERO eventos y la derivada no dice nada: se
    # arranca donde el espectro supera el 2 % de su maximo
    i0 = int(np.flatnonzero(ys > 0.02 * ys.max())[0])
    sube = np.diff(ys) > 0
    if i0 >= sube.size or not sube[i0]:
        return float(x[0])
    deja = np.flatnonzero(~sube[i0:])
    return float(x[i0 + deja[0]]) if deja.size else float(x[0])


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def de_corrida(clave, limite=None):
    """El espectro neto de una corrida de trabajo (muestra y, si hay, fondo)."""
    from desde_cero.lector import CORRIDAS_TRABAJO
    ruta_m, ruta_f = CORRIDAS_TRABAJO[clave]
    d_m = area.areas_de_archivo(ruta_m, limite=limite)
    d_f = area.areas_de_archivo(ruta_f, limite=limite) if ruta_f else None
    return espectro_neto(d_m, d_f), d_m, d_f


def resumen(clave, e):
    t = ("  %-12s vivo %.1f de %.1f s de pared (%.2f %%, scaler mediano %.0f Hz)"
         % (clave, e["t_vivo_m"], e["t_pared_m"], 100 * e["frac_viva_m"],
            e["scaler_m"]))
    if e["hay_fondo"]:
        t += ("\n  %-12s fondo: vivo %.1f de %.1f s (%.2f %%, %.0f Hz)   "
              "k_m = %.4f   k_f = %.5f"
              % ("", e["t_vivo_f"], e["t_pared_f"], 100 * e["frac_viva_f"],
                 e["scaler_f"], e["k_m"], e["k_f"]))
    else:
        t += "\n  %-12s SIN fondo: el neto es la muestra cruda" % ""
    return t


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero.lector import CORRIDAS_TRABAJO

    claves = [k for k, (m, _) in CORRIDAS_TRABAJO.items() if os.path.exists(m)]
    fig, ejes = plt.subplots(1, len(claves), figsize=(6 * len(claves), 4.6),
                             squeeze=False)
    for ax, clave in zip(ejes[0], claves):
        e, _, _ = de_corrida(clave)
        piso = fin_del_umbral(e["canales"], e["neto"])
        print(resumen(clave, e))
        print("  %-12s %d bines de %.2f mV*ns hasta %.0f; fin del umbral en "
              "%.0f mV*ns" % ("", BINES, e["tope"] / BINES, e["tope"], piso))
        x = e["canales"]
        ax.step(x, e["muestra"], where="mid", lw=0.8, color="0.5",
                label="muestra x k_m")
        if e["hay_fondo"]:
            ax.step(x, e["fondo"], where="mid", lw=0.8, color="tab:red",
                    label="fondo x k_f")
        ax.errorbar(x, e["neto"], e["error"], fmt="none", lw=0.5,
                    color="tab:blue", alpha=0.6)
        ax.step(x, e["neto"], where="mid", lw=1, color="tab:blue", label="neto")
        ax.axvline(piso, color="k", ls=":", label="fin del umbral")
        ax.set(xlabel="area [mV*ns]", ylabel="cuentas / bin (al tiempo vivo comun)",
               title=clave)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa3_espectro.png")
    fig.savefig(salida, dpi=120)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
