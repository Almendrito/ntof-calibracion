"""
espectro.py — ETAPA 3: de las areas al espectro neto.

Tres operaciones, y las tres tienen una trampa adentro:

1. HISTOGRAMEAR muestra y fondo con EL MISMO eje. Si cada uno se histogramea
   con su propio `range`, los bins de uno no son los del otro y restarlos es
   restar peras de manzanas. Un solo `rango` fijo para los dos, siempre.

2. NORMALIZAR POR TIEMPO VIVO, no por tiempo de pared. El DRSOSC escribe ~35 kB
   de XML por evento y a tasa alta NO alcanza a guardar todos los triggers que
   la placa dispara. El `<ScalerN>` de cada evento trae la tasa REAL, asi que

       triggers verdaderos  ~  integral de Scaler dt
       fraccion viva        =  eventos guardados / triggers verdaderos
       T_vivo               =  fraccion viva * duracion de pared

   Sin esto la resta esta mal en un factor que puede ser enorme: la corrida de
   Co-60 estuvo viva el 40% del tiempo y su fondo el 100%, o sea que
   normalizando por pared el fondo entra 2.5 veces mas fuerte de lo que debe.

3. PROPAGAR EL ERROR. Cada histograma crudo tiene var = N (Poisson); al
   escalarlo por k la varianza va como k^2*N. El neto es una RESTA, asi que las
   dos varianzas se suman. Sin barras de error no hay chi2, no hay bondad de
   ajuste y no hay barra de error sobre la calibracion: es la diferencia entre
   un numero y un numero con incertidumbre.

Correr `python analisis/espectro.py` deja `figuras/03_espectro_neto.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, area, pulso

# Tope del histograma cuando no se fija a mano: `RANGO_FACTOR` veces el
# percentil `RANGO_PERCENTIL` de las areas de la muestra. Estan a nivel de
# modulo para que el informe los cite en vez de tipearlos.
RANGO_FACTOR = 1.6
RANGO_PERCENTIL = 99.5


# --------------------------------------------------------------------------
# 1. Tiempo vivo
# --------------------------------------------------------------------------
def tiempo_vivo(reloj_s, scaler, duracion, n_eventos, recorte=10.0):
    """
    Tiempo VIVO de la adquisicion [s] a partir del scaler de la placa.

    `recorte` : el scaler tiene rafagas absurdas (se han visto 90 kHz en una
        corrida de fondo que promedia 10 Hz, ruido electronico del disparo). Se
        lo recorta a `recorte` veces su mediana antes de integrar; si no, una
        sola rafaga se come el resultado.

    Si no hay scaler utilizable se devuelve la duracion de pared y se asume
    100% vivo, que es lo conservador: sub-restar fondo es preferible a
    sobre-restarlo, porque lo segundo puede dejar el neto negativo y no hay como
    darse cuenta mirando el espectro.

    Las marcas `<Time>` NO sirven para esto: el DRSOSC las escribe en rafagas
    (miles de eventos comparten el mismo milisegundo), asi que dan la duracion
    total pero no la tasa instantanea.
    """
    reloj_s = np.asarray(reloj_s, float)
    scaler = np.asarray(scaler, float)
    ok = np.isfinite(scaler) & (scaler > 0)
    if reloj_s.size < 2 or ok.sum() < 10 or duracion <= 0:
        return float(duracion), 1.0, np.nan

    med = float(np.median(scaler[ok]))
    s_lim = np.clip(np.nan_to_num(scaler, nan=med), 0.0, recorte * max(med, 1.0))
    dt = np.clip(np.diff(reloj_s, prepend=reloj_s[0]), 0.0, None)
    n_true = float(np.sum(s_lim * dt))
    if n_true <= 0:
        return float(duracion), 1.0, med
    frac = min(n_eventos / n_true, 1.0)      # >1 solo por cuantizacion del scaler
    return float(frac * duracion), float(frac), med


# --------------------------------------------------------------------------
# 2. El espectro neto
# --------------------------------------------------------------------------
def espectro_neto(d_muestra, d_fondo=None, bins=500, rango_max=None,
                  usar_tiempo_vivo=True, permitir=pulso.MAL_SIN_PRE,
                  recortar_negativos=False, verbose=True):
    """
    Arma el espectro a partir de los dicts que devuelve `area.areas_de_archivo`.

    rango_max : tope del histograma [mV*ns]. Si es None se toma `RANGO_FACTOR`
        veces el percentil `RANGO_PERCENTIL` de la muestra. Fijarlo importa: unos pocos eventos de
        area gigante (pile-up, cosmicos) estiran el rango automatico y dejan
        toda la fisica apretada en los primeros bins.
    permitir : banderas que NO descartan. El MISMO corte se aplica a muestra y
        fondo; si no, la resta deja de tener sentido porque cada histograma
        tendria una eficiencia distinta.
    recortar_negativos : True SOLO para graficar. Para ajustar hay que conservar
        el neto real, negativos incluidos: son fluctuaciones estadisticas
        legitimas y recortarlas sesga el ajuste hacia arriba.

    Devuelve un dict con canales, muestra, fondo (ya normalizados), neto, error,
    los factores k y los tiempos.
    """
    a_m = d_muestra["area"][area.buenos(d_muestra, permitir)]
    if rango_max is None:
        rango_max = float(np.percentile(a_m, RANGO_PERCENTIL) * RANGO_FACTOR)
    rango = (0.0, float(rango_max))

    h_m, bordes = np.histogram(a_m, bins=bins, range=rango)
    canales = 0.5 * (bordes[:-1] + bordes[1:])

    tv_m, frac_m, sc_m = tiempo_vivo(d_muestra["reloj_s"], d_muestra["scaler"],
                                     d_muestra["duracion"],
                                     d_muestra["n_eventos"])
    if not usar_tiempo_vivo:
        tv_m, frac_m = float(d_muestra["duracion"]), 1.0

    if d_fondo is None:
        # Sin fondo NO se inventa una resta. Se devuelve el crudo y se dice.
        h_f = np.zeros_like(h_m)
        tv_f, frac_f, sc_f = 0.0, np.nan, np.nan
        k_m = k_f = 1.0
        neto = h_m.astype(float)
        error = np.sqrt(np.maximum(h_m, 1.0))
        m_norm, f_norm = neto, h_f.astype(float)
    else:
        a_f = d_fondo["area"][area.buenos(d_fondo, permitir)]
        h_f, _ = np.histogram(a_f, bins=bins, range=rango)   # MISMOS bordes
        tv_f, frac_f, sc_f = tiempo_vivo(d_fondo["reloj_s"], d_fondo["scaler"],
                                         d_fondo["duracion"],
                                         d_fondo["n_eventos"])
        if not usar_tiempo_vivo:
            tv_f, frac_f = float(d_fondo["duracion"]), 1.0
        tv_m = tv_m if tv_m > 0 else float(d_muestra["duracion"])
        tv_f = tv_f if tv_f > 0 else float(d_fondo["duracion"])

        # Se normaliza al MENOR de los dos tiempos vivos: asi ninguno de los dos
        # histogramas se AMPLIFICA, y las barras de error siguen siendo las de
        # un conteo real escalado hacia abajo.
        t_norm = min(tv_m, tv_f)
        k_m = t_norm / tv_m if tv_m > 0 else 1.0
        k_f = t_norm / tv_f if tv_f > 0 else 1.0
        m_norm = h_m * k_m
        f_norm = h_f * k_f
        neto = m_norm - f_norm
        error = np.sqrt(np.maximum(h_m * k_m ** 2 + h_f * k_f ** 2, 1.0))

    if recortar_negativos:
        neto = np.maximum(neto, 0.0)

    res = dict(canales=canales, bordes=bordes, muestra=m_norm, fondo=f_norm,
               neto=neto, error=error, k_m=float(k_m), k_f=float(k_f),
               t_pared_muestra=float(d_muestra["duracion"]), t_vivo_muestra=tv_m,
               frac_viva_muestra=frac_m, scaler_muestra=sc_m,
               t_pared_fondo=float(d_fondo["duracion"]) if d_fondo else 0.0,
               t_vivo_fondo=tv_f, frac_viva_fondo=frac_f, scaler_fondo=sc_f,
               rango_max=float(rango_max), bins=int(bins),
               n_muestra=int(a_m.size),
               n_fondo=int(d_fondo["area"][area.buenos(d_fondo, permitir)].size)
               if d_fondo else 0,
               metodo=str(d_muestra["metodo"]), pre=float(d_muestra["pre"]),
               post=float(d_muestra["post"]), hay_fondo=d_fondo is not None)
    if verbose:
        informe(res)
    return res


def fin_del_umbral(canales, neto, suave=9):
    """
    Canal donde termina la subida del trigger. Piso duro para cualquier ventana.

    Por debajo del umbral del discriminador el espectro no es fisica: es la
    curva de eficiencia del trigger, que sube de cero hasta uno. Ni la plantilla
    de Geant4 ni la componente de exceso la contienen, asi que ajustar ahi
    adentro no sesga: rompe. Medido sobre el Cs-137 a 1600 V, que tiene cero
    cuentas debajo de 121 mV*ns y sube hasta ~350:

        arranque de la ventana    229      320
        chi2/ndf                 4.255    0.953

    El criterio es el primer maximo local del espectro suavizado: la eficiencia
    del trigger CRECE y todo lo que hay debajo del borde (continuo Compton mas
    el exceso de bajo canal) DECRECE, asi que el punto donde deja de crecer es
    donde una cosa le pasa el mando a la otra. Cuando no hay umbral visible
    -el Co-60, tomado con el discriminador en 0.2 mV- el maximo cae en los
    primeros bines y la funcion no restringe nada, que es lo que corresponde.
    """
    x = np.asarray(canales, float)
    y = np.asarray(neto, float)
    if y.size < 2 * suave:
        return float(x[0])
    k = np.ones(int(suave)) / float(suave)
    ys = np.convolve(y, k, mode="same")
    if not np.isfinite(ys).any() or ys.max() <= 0:
        return float(x[0])
    # arrancar donde el espectro deja de ser CERO: por debajo del umbral no hay
    # ni un evento, y ahi la derivada es 0 y no dice nada. Sin esto la busqueda
    # termina en el primer bin y la funcion no protege de nada.
    vivos = np.flatnonzero(ys > 0.02 * ys.max())
    if vivos.size == 0:
        return float(x[0])
    i0 = int(vivos[0])
    sube = np.diff(ys) > 0
    if i0 >= sube.size or not sube[i0]:
        return float(x[0])            # ya viene bajando: no hay subida de trigger
    resto = ~sube[i0:]
    if not resto.any():
        return float(x[0])
    return float(x[i0 + int(np.argmax(resto))])


def informe(e):
    """Diagnostico del espectro en texto: lo que hay que mirar antes de ajustar."""
    print("  espectro: %d bins de %.2f mV*ns, rango [0, %.0f], metodo %s "
          "(ventana %.1f + %.1f ns)"
          % (e["bins"], e["rango_max"] / e["bins"], e["rango_max"], e["metodo"],
             e["pre"], e["post"]))
    print("  muestra : %d pulsos, %.0f s de pared -> %.0f s VIVOS (%.1f%%), "
          "scaler %.0f Hz"
          % (e["n_muestra"], e["t_pared_muestra"], e["t_vivo_muestra"],
             100 * e["frac_viva_muestra"], e["scaler_muestra"]))
    if e["hay_fondo"]:
        print("  fondo   : %d pulsos, %.0f s de pared -> %.0f s VIVOS (%.1f%%), "
              "scaler %.0f Hz"
              % (e["n_fondo"], e["t_pared_fondo"], e["t_vivo_fondo"],
                 100 * e["frac_viva_fondo"], e["scaler_fondo"]))
        print("  factores: k_muestra = %.4f, k_fondo = %.4f   (con tiempo de "
              "PARED habrian sido %.4f y %.4f)"
              % (e["k_m"], e["k_f"],
                 min(e["t_pared_muestra"], e["t_pared_fondo"]) / e["t_pared_muestra"],
                 min(e["t_pared_muestra"], e["t_pared_fondo"]) / e["t_pared_fondo"]))
        peso = e["fondo"].sum() / max(e["muestra"].sum(), 1e-9)
        print("  el fondo pesa el %.2f%% de la muestra" % (100 * peso))
        neg = int((e["neto"] < -2 * e["error"]).sum())
        print("  bins con neto negativo a mas de 2 sigma: %d de %d %s"
              % (neg, e["bins"],
                 "(normal: son fluctuaciones)" if neg < 0.05 * e["bins"]
                 else "<- MUCHOS: revisar la normalizacion temporal"))
    else:
        print("  fondo   : NO HAY. El espectro es la muestra CRUDA, sin restar.")
    x, y = e["canales"], e["neto"]
    m = y > 0
    if m.any():
        acum = np.cumsum(y[m]) / y[m].sum()
        for q in (0.5, 0.9, 0.99):
            print("     el %2d%% del neto esta por debajo de %.0f mV*ns"
                  % (100 * q, x[m][np.searchsorted(acum, q)]))


# --------------------------------------------------------------------------
# 3. Guardar / cargar
# --------------------------------------------------------------------------
def guardar(e, nombre):
    """Escribe el espectro en `.cache/<nombre>.npz`. Es la entrada del ajuste."""
    ruta = rutas.cache("%s.npz" % nombre)
    np.savez_compressed(ruta, **{k: v for k, v in e.items() if v is not None})
    return ruta


def cargar(nombre):
    ruta = rutas.cache("%s.npz" % nombre)
    d = np.load(ruta, allow_pickle=False)
    return {k: (d[k] if d[k].ndim else d[k].item()) for k in d.files}


# --------------------------------------------------------------------------
# 4. La cadena completa para una corrida
# --------------------------------------------------------------------------
def de_corrida(clave, metodo="global", limite=None, bins=500, rango_max=None,
               pre=None, post=None, verbose=True, **kw):
    """
    XML -> areas -> espectro, para una clave del catalogo de `rutas`.

    Si no se dan `pre`/`post` la ventana se MIDE sobre la propia muestra y se
    usa la MISMA para el fondo. Esto ultimo no es un atajo: muestra y fondo se
    tomaron dias distintos y con retardos de trigger distintos, asi que si cada
    uno eligiera su ventana estarian midiendo observables diferentes y la resta
    no significaria nada. La ventana viaja anclada al onset de cada evento, asi
    que usar la misma no impone ninguna suposicion sobre CUANDO llego el pulso.
    """
    c = rutas.corrida(clave)
    if verbose:
        print("=" * 78)
        print("%s  -  %s" % (clave, c.descripcion))
        print("=" * 78)

    if metodo == "global" and (pre is None or post is None):
        v = area.calibrar_ventana(c.muestra, verbose=verbose)
        pre, post = v["pre"], v["post"]

    d_m = area.areas_de_archivo(c.muestra, metodo=metodo, pre=pre, post=post,
                                limite=limite, verbose=verbose, **kw)
    d_f = None
    if c.fondo and os.path.exists(c.fondo):
        d_f = area.areas_de_archivo(c.fondo, metodo=metodo, pre=pre, post=post,
                                    limite=limite, verbose=verbose, **kw)
    elif verbose:
        print("  (esta corrida no tiene fondo propio: no se resta nada)")

    e = espectro_neto(d_m, d_f, bins=bins, rango_max=rango_max, verbose=verbose)
    e["corrida"] = clave
    e["fuente"] = c.fuente or ""
    return e


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo(limite=8000, bins=400):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES if rutas.corrida(k).disponible()]
    fig, axes = plt.subplots(len(claves), 2, figsize=(14, 4.6 * len(claves)),
                             squeeze=False)

    for fila, clave in enumerate(claves):
        e = de_corrida(clave, metodo="global", limite=limite, bins=bins)
        print()

        ax = axes[fila][0]
        ax.step(e["canales"], e["muestra"], where="mid", lw=1.1, label="muestra")
        if e["hay_fondo"]:
            ax.step(e["canales"], e["fondo"], where="mid", lw=1.1,
                    label="fondo (normalizado por tiempo VIVO)")
        ax.set_yscale("log")
        ax.set_title("%s: muestra y fondo sobre el MISMO eje" % clave)
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[fila][1]
        ax.errorbar(e["canales"], e["neto"], yerr=e["error"], fmt="none",
                    ecolor="0.75", lw=0.6)
        ax.step(e["canales"], e["neto"], where="mid", lw=1.2, color="C2",
                label="neto = muestra - fondo")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title("espectro NETO (%s, %d eventos leidos)"
                     % (e["fuente"], limite))
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        # Ojo con el nombre: esta demo corre con `limite`, o sea produce un
        # espectro PARCIAL. Si se guardara como `espectro_<corrida>_<metodo>`
        # pisaria el resultado completo de `hacer_espectro.py` sin dejar rastro
        # (paso). El prefijo `demo_` lo mantiene aparte.
        guardar(e, "demo_espectro_%s_%s" % (clave, e["metodo"]))

    fig.tight_layout()
    salida = rutas.figura("03_espectro_neto.png")
    fig.savefig(salida, dpi=130)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    _demo(limite=n)
