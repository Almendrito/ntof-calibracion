"""
area.py — ETAPA 2: aplicar la Etapa 1 a millones de eventos.

Lo que agrega sobre `pulso.py`:

  * `calibrar_ventana`   la ventana GLOBAL se MIDE sobre una submuestra del
                         propio archivo (no se tipea, no se hereda de otra corrida)
  * `areas_de_archivo`   recorre el XML entero, con cache
  * la comparacion de los dos metodos, que es el entregable de esta etapa

Los dos metodos, y por que se comparan
--------------------------------------
A) GLOBAL   : ancho fijo `[onset-pre, onset+post]`, con `pre`/`post` medidos una
              vez sobre el pulso promedio alineado, y anclado al onset DE CADA
              EVENTO. El jitter del trigger no lo toca porque la ventana viaja
              con el pulso.
B) CONVERGENCIA : cada evento decide su ventana con el criterio de la Etapa 1.
              Es la receta literal.

La Etapa 1 ya mostro que en B el ancho crece con la amplitud (13 -> 42 ns entre
el primer y el ultimo cuartil del Co-60). Eso NO es un detalle de implementacion:
un corte que depende de la amplitud le cambia la forma al espectro, que es lo
unico que la calibracion mide. Aca se cuantifica en mV*ns.

El test que decide
------------------
`area / amplitud` en funcion de la amplitud tiene que ser PLANO. Si la ventana
recorta una fraccion distinta de la carga segun el tamano del pulso, esa razon
se inclina, y con ella se va la proporcionalidad entre area y energia: el eje
del espectro deja de ser lineal y ninguna calibracion posterior lo arregla. Es
un test barato que no depende de conocer la fuente ni la energia.

Correr `python analisis/area.py` deja `figuras/02_comparacion_metodos.png`.
"""

import hashlib
import os
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, drs, pulso

# Umbral de ruido pico a pico [mV]: por debajo de esto el registro es linea base
# y no hay pulso que integrar. Es deliberadamente flojo (el ruido rms mide entre
# 0.3 y 1.0 mV segun la corrida): la Etapa 3 no puede recuperar un evento que
# esta etapa tiro, asi que se prefiere marcar antes que descartar.
UMBRAL_RUIDO = 3.0

VERSION = 3          # sube esto para invalidar todos los caches a proposito

# --------------------------------------------------------------------------
# 0. Las plantillas del pulso (metodo="plantilla", 24-sep-2026)
# --------------------------------------------------------------------------
# De que corridas se MIDE cada forma. Tienen que ser de registro largo (los XML:
# ~130-158 ns tras el onset) porque la plantilla tiene que llegar hasta
# `pulso.T_CARGA`.
#
# Son DOS y no una, y no por detector sino por MONTAJE. Medido sobre las 28
# corridas (verificaciones/_verif_area_plantilla.py): el Co-60 de junio decae a
# cero hacia los 50 ns con un rebote chico a ~25 ns, mientras que TODAS las de
# ago/sep -detectores A, B, C y D, XML y binario, a los tres voltajes- tienen
# un segundo pulso de ~4.5% a ~78 ns del onset (una reflexion: el retardo es
# fijo) y una cola casi plana de ~1.5% de la amplitud que a 130 ns sigue viva.
# Entre las de ago/sep el residuo contra la plantilla es 0.3-0.8% de la
# amplitud; el de junio contra la de ago/sep es 1.4%. Cambio el cableado entre
# junio y agosto, no el detector.
#
# Los detectores B, C y D no tienen corridas de registro largo: su forma se
# verifico hasta donde llega su registro (81 ns los del 1/2-sep, 47 ns el D del
# 8-sep). La cola que se les EXTRAPOLA es la del detector A, y con estos datos
# eso es una hipotesis, no una medicion.
PLANTILLAS = {
    "junio": ("co60",),
    "ago-sep": ("cs137_1600", "cs137_1650", "cs137_1700"),
}
CAMBIO_DE_MONTAJE = "2026-08-01"    # antes de esta fecha -> "junio"
VERSION_PLANTILLA = 1                # sube esto si cambia `pulso.plantilla_pulso`


def plantilla_para(ruta):
    """
    Que plantilla le corresponde a un archivo, por su FECHA de medicion (la
    forma cambio con el montaje, CLAUDE.md 7.74). Muestra y fondo tienen que
    usar la MISMA: `pipeline` le pasa al fondo la de la muestra, igual que con
    la compuerta, porque un fondo con otra plantilla mediria otro observable y
    la resta no significaria nada.
    """
    f = rutas.fecha(ruta)
    return "junio" if f and f < CAMBIO_DE_MONTAJE else "ago-sep"


_PLANTILLAS_MEM = {}


def cargar_plantilla(nombre="ago-sep", n=1500, rehacer=False, verbose=False):
    """
    La plantilla `nombre`, de `.cache/plantilla_pulso_<nombre>.npz` si esta y
    sus archivos fuente no cambiaron; si no, se mide (~7 s) y se guarda.

    La clave es la misma idea que la de las areas: tamano + fecha de cada
    archivo fuente, mas los parametros que cambian la forma. El HASH de la
    forma resultante es lo que entra despues en la clave de las areas, asi que
    una plantilla nueva invalida sola todas las areas que se calcularon con la
    vieja.
    """
    claves = PLANTILLAS[nombre]
    archivos = [rutas.corrida(k).muestra for k in claves]
    firma = "|".join(["%d" % VERSION_PLANTILLA, "%d" % n,
                      "%.3f" % pulso.T_CARGA, repr(pulso.T_PLANTILLA),
                      "%.1f" % pulso.AMP_PLANTILLA]
                     + ["%d:%d" % (os.stat(a).st_size, int(os.stat(a).st_mtime))
                        for a in archivos])
    if not rehacer and nombre in _PLANTILLAS_MEM \
            and _PLANTILLAS_MEM[nombre][1] == firma:
        return _PLANTILLAS_MEM[nombre][0]
    f = rutas.cache("plantilla_pulso_%s.npz" % nombre)
    if not rehacer and os.path.exists(f):
        d = np.load(f, allow_pickle=False)
        if str(d["firma"]) == firma:
            p = {k: (d[k] if d[k].ndim else d[k].item()) for k in d.files
                 if k != "firma"}
            p["archivos"] = [str(a) for a in p["archivos"]]
            _PLANTILLAS_MEM[nombre] = (p, firma)
            return p
    if verbose:
        print("  [plantilla] midiendo la forma '%s' sobre %s ..."
              % (nombre, ", ".join(claves)))
    p = pulso.plantilla_pulso(archivos, n=n)
    p["nombre"] = nombre
    np.savez_compressed(f, firma=firma, **{k: np.asarray(v) for k, v in p.items()})
    _PLANTILLAS_MEM[nombre] = (p, firma)
    return p


# --------------------------------------------------------------------------
# 1. La ventana global se MIDE
# --------------------------------------------------------------------------
def espacio_tras_el_onset(T, V, n_base=40):
    """
    Cuanto registro queda despues del onset, evento a evento [ns], y con que
    amplitud. Devuelve `(espacio, amplitud, t_fin)`.

    Las dos cosas van juntas porque el trigger tiene WALK: los pulsos grandes
    cruzan el umbral antes, asi que su onset es mas temprano y les sobra mas
    registro. La consecuencia practica es que un percentil crudo del espacio
    esta dominado por los pulsos chicos -los que menos importan- y recorta la
    compuerta de mas.
    """
    fin = float(np.median(T[:, -1]))
    espacio = np.empty(T.shape[0])
    amplitud = np.empty(T.shape[0])
    for k in range(T.shape[0]):
        t, v = T[k], V[k]
        base = float(np.median(v[:n_base]))
        ip, _, amp, signo = pulso.detectar_pico(t, v - base)
        _, t_on = pulso.detectar_onset(t, v, base, ip, signo)
        espacio[k] = fin - t_on
        amplitud[k] = amp
    return espacio, amplitud, fin


def calibrar_ventana(ruta, n=600, canal=None, k=3.0, tramo=4.0,
                     ancho_suave=1.0, margen=0.0, frac_eventos=0.95,
                     margen_registro=5.0, verbose=True):
    """
    Mide `(pre, post)` sobre el pulso promedio alineado del propio archivo.

    Se hace con el PROMEDIO y no evento por evento porque el criterio de
    convergencia compara la senal contra el ruido, y el ruido del promedio de
    ~400 pulsos es ~20 veces menor que el de uno solo. O sea: el promedio ve la
    cola hasta mucho mas abajo, y por eso la ventana que da es la que de verdad
    contiene el pulso. La misma ventana aplicada evento por evento no depende ya
    de la amplitud de ese evento — que es todo el punto.

    `margen` [ns] se le suma a `post` si uno quiere ser conservador.

    La compuerta tiene que CABER en el registro
    -------------------------------------------
    Una ventana que el registro no contiene no es una ventana: el evento sale
    truncado y pierde una fraccion de carga que depende de cuando llego el
    pulso, o sea justo el error que la compuerta anclada al onset existe para
    evitar. Asi que `post` se recorta al espacio disponible.

    El espacio se mide en los eventos de amplitud por encima de la mediana (ver
    `espacio_tras_el_onset`) y se le descuenta `margen_registro`. Medido en el
    Co-60 del 01-sep-2026, donde el trigger quedo en t ~ 123 ns de un registro
    de 200 ns:

        post pedido por el pulso : 83.0 ns  ->  99.9% de la carga, 16% de
                                                eventos utiles
        post = 80 ns (sin margen): 99.5% de la carga, 71% de eventos utiles
        post = 75 ns (con margen): 98.0% de la carga, 94% de eventos utiles

    y la D corregida por la fraccion de carga es la MISMA en los dos ultimos
    (1.948 contra 1.949), asi que los 5 ns de margen no cuestan nada y
    recuperan un cuarto de la estadistica.

    En los archivos de junio y del 26-ago no recorta nada: ahi el pulso cae a
    t ~ 48 ns y sobran ~120-155 ns de registro.

    Devuelve un dict con `pre`, `post`, la fraccion de carga que la ventana
    recoge, la forma promedio (para poder graficarla) y, si hubo recorte,
    `post_pedido` y `recortado`.
    """
    T, V, _, _ = drs.leer_pulsos(ruta, n=n, canal=canal)
    tau_rel, forma, n_us = pulso.pulso_promedio(T, V)
    pre_zona = tau_rel < -8.0
    ruido = float(1.4826 * np.median(np.abs(forma[pre_zona]
                                            - np.median(forma[pre_zona]))))
    t_ini, t_fin, _ = pulso.ventana_por_convergencia(
        tau_rel, forma, 0.0, 1.0, ruido, k=k, tramo=tramo,
        ancho_suave=ancho_suave)
    _, frac, info = pulso.curva_convergencia(tau_rel, forma)

    post_pedido = float(t_fin + margen)
    espacio, amplitud, t_reg = espacio_tras_el_onset(T, V)
    grandes = espacio[amplitud >= np.median(amplitud)]
    disponible = float(np.percentile(grandes if grandes.size > 20 else espacio,
                                     100.0 * (1.0 - frac_eventos))
                       - margen_registro)
    post = min(post_pedido, max(disponible, 20.0))
    recortado = post < post_pedido - 0.05

    i_fin = min(int(np.searchsorted(tau_rel, post)), frac.size - 1)
    res = dict(pre=float(max(-t_ini, 2.0)), post=float(post),
               post_pedido=post_pedido, disponible=disponible,
               recortado=bool(recortado),
               frac_eventos_enteros=float((espacio >= post).mean()),
               espacio=espacio, amplitud=amplitud, t_registro=t_reg,
               fraccion=float(frac[i_fin]), t_plateau=info["t_plateau"],
               caida_final=info["caida_final"], n_usados=n_us,
               tau_rel=tau_rel, forma=forma, ruido_forma=ruido)
    if verbose:
        print("  ventana medida sobre %d pulsos: pre = %.1f ns, post = %.1f ns "
              "(ancho %.1f ns) -> recoge el %.2f%% de la carga"
              % (n_us, res["pre"], res["post"], res["pre"] + res["post"],
                 100 * res["fraccion"]))
        if recortado:
            print("  OJO: el pulso pedia post = %.1f ns pero el registro deja "
                  "%.1f ns" % (post_pedido, disponible + margen_registro))
            print("       tras el onset (eventos de amplitud > mediana). Con "
                  "post = %.1f ns" % post)
            print("       entra entero el %.1f%% de los eventos."
                  % (100 * res["frac_eventos_enteros"]))
    return res


def fraccion_de_carga(v, post):
    """
    Que fraccion de la carga del pulso promedio recoge una compuerta de `post`.

    Sale de la curva de convergencia que `calibrar_ventana` ya calculo. Es el
    numero que hace comparables (o no) las D de dos corridas con ventanas
    distintas: dos corridas con `post` distinto miden fracciones distintas del
    mismo pulso, y eso entra derecho en D.
    """
    _, frac, _ = pulso.curva_convergencia(v["tau_rel"], v["forma"])
    i = int(np.searchsorted(v["tau_rel"], float(post)))
    return float(frac[min(i, frac.size - 1)])


# --------------------------------------------------------------------------
# 2. Recorrer el archivo entero
# --------------------------------------------------------------------------
def _clave(ruta, canal, metodo, params):
    """
    Huella del CONTENIDO + de todo lo que cambia el numero.

    El canal entra en la clave, y no es paranoia: un XML leido del canal
    equivocado da cero pulsos sin decir nada, y sin el canal en la clave ese
    vacio se cachearia para siempre, incluso despues de arreglar el codigo.

    La RUTA no entra, y eso es un arreglo del 03-sep-2026. Estaba, y por lo
    tanto mover un archivo de carpeta -o copiarlo, que es lo que se hizo al
    traer los 12 GB del disco externo a `datos/`- invalidaba el cache y
    obligaba a re-parsear 80 s por XML, aunque el archivo fuera byte a byte el
    mismo. Lo que identifica a estos datos es tamano + fecha de modificacion,
    que es justo lo que `copy2` y `mv` conservan (es la misma clave que usa
    `codigo/datos_io/spectrum.py`, y por eso ahi la reorganizacion de carpetas
    de §7.43 no costo nada). Dos archivos DISTINTOS con el mismo tamano y la
    misma fecha al segundo no existen en este proyecto.
    """
    st = os.stat(ruta)
    firma = "|".join([
        str(st.st_size), str(int(st.st_mtime)),
        str(canal), str(metodo), str(VERSION),
        "|".join("%s=%s" % (k, params[k]) for k in sorted(params)),
    ])
    return hashlib.sha1(firma.encode("utf-8")).hexdigest()[:16], firma


def _params_cache(pre, post, umbral_ruido, limite, kw, plantilla=None):
    """Los parametros que entran en la clave del cache, en un solo lugar.

    Con `metodo="plantilla"` entra ademas el HASH de la forma (no su nombre):
    si la plantilla se re-mide y cambia, las areas viejas dejan de valer solas.
    Sin plantilla la clave queda identica a la de antes, asi que los caches del
    metodo global siguen valiendo.
    """
    params = dict(pre=("%.3f" % pre) if pre is not None else "auto",
                  post=("%.3f" % post) if post is not None else "auto",
                  umbral=umbral_ruido, limite=limite or 0,
                  **{k: kw[k] for k in sorted(kw)})
    if plantilla is not None:
        params["plantilla"] = plantilla["hash"]
    return params


def _resolver_plantilla(ruta, metodo, plantilla):
    """None si el metodo no la usa; si no, el dict de la plantilla (acepta el
    dict, un nombre de `PLANTILLAS`, o None = la que le toca por fecha)."""
    if metodo != "plantilla":
        return None
    if plantilla is None:
        plantilla = plantilla_para(ruta)
    if isinstance(plantilla, str):
        plantilla = cargar_plantilla(plantilla)
    return plantilla


def hay_cache(ruta, metodo="global", pre=None, post=None, canal=None,
              umbral_ruido=UMBRAL_RUIDO, limite=None, plantilla=None, **kw):
    """
    True si las areas de este archivo YA estan cacheadas con estos parametros.

    Sirve para pesar la barra de progreso ANTES de empezar: parsear un XML de
    1.76 GB son ~40 s y leer su cache son milisegundos, asi que una etapa que va
    a salir del cache no tiene que reservarse la mitad de la barra. Sin esto el
    tiempo restante estimado sale absurdo desde el primer segundo.
    """
    try:
        if canal is None:
            canal, _ = drs.detectar_canal(ruta)
        params = _params_cache(pre, post, umbral_ruido, limite, kw,
                               _resolver_plantilla(ruta, metodo, plantilla))
        clave, firma = _clave(ruta, canal, metodo, params)
        f = rutas.cache("areas_%s.npz" % clave)
        if not os.path.exists(f):
            return False
        d = np.load(f, allow_pickle=False)
        return str(d["firma"]) == firma
    except Exception:
        return False


def areas_de_archivo(ruta, metodo="global", pre=None, post=None, canal=None,
                     limite=None, umbral_ruido=UMBRAL_RUIDO, usar_cache=True,
                     verbose=True, callback=None, plantilla=None, **kw):
    """
    Recorre el XML y devuelve un dict de arrays, uno por evento util.

    metodo : "global" (usa `pre`/`post`), "convergencia" (cada evento el suyo)
             o "plantilla" (ajuste de la forma del pulso, ver
             `pulso.area_evento_plantilla`). Si `metodo` es "global" o
             "plantilla" y no se pasa `pre`/`post`, se miden aca mismo con
             `calibrar_ventana` sobre el propio archivo.
    plantilla : solo con `metodo="plantilla"`. El dict de `cargar_plantilla`,
             un nombre de `PLANTILLAS`, o None para la que le toca al archivo
             por fecha (`plantilla_para`).

    Con "plantilla" el area es la del AJUSTE, y en la misma pasada se calcula
    la de la compuerta global: la bandera -y por lo tanto los descartes- es
    exactamente la del metodo global, y se guarda ademas `area_compuerta`, `A`,
    `t0`, `b_ajuste`, `chi2_red`, `extrapolada` y `t0_en_borde` por evento.

    Claves: area, amplitud, t_onset, t_ini, t_fin, base, ruido, bandera,
    fraccion, reloj, scaler; mas los escalares n_eventos, duracion, canal,
    metodo, pre, post y `descartados` (por umbral de ruido / trama repetida).

    Dos filtros baratos antes de integrar, los dos heredados de la practica:
      - trama REPETIDA consecutiva: el DRSOSC a veces vuelca dos veces el mismo
        registro; contarlo dos veces mete un pico falso en el espectro.
      - ruido pico a pico por debajo de `umbral_ruido`: no hay pulso.

    Si el archivo termina con CERO pulsos utiles se levanta ValueError en vez de
    devolver arrays vacios. Un vacio silencioso se cachea y despues aparece
    cinco etapas mas abajo como un histograma de ceros, sin ninguna pista de
    donde vino.
    """
    if canal is None:
        canal, _ = drs.detectar_canal(ruta)

    if metodo not in ("global", "convergencia", "plantilla"):
        raise ValueError("metodo desconocido: %r" % metodo)
    if metodo in ("global", "plantilla") and (pre is None or post is None):
        if verbose:
            print("  (no me dieron ventana: la mido sobre este mismo archivo)")
        v = calibrar_ventana(ruta, canal=canal, verbose=verbose)
        pre, post = v["pre"], v["post"]
    plantilla = _resolver_plantilla(ruta, metodo, plantilla)

    params = _params_cache(pre, post, umbral_ruido, limite, kw, plantilla)
    clave, firma = _clave(ruta, canal, metodo, params)
    archivo_cache = rutas.cache("areas_%s.npz" % clave)

    if usar_cache and os.path.exists(archivo_cache):
        d = np.load(archivo_cache, allow_pickle=False)
        if str(d["firma"]) == firma:
            if verbose:
                print("  [cache] %s: %d pulsos (%s)"
                      % (os.path.basename(ruta), d["area"].size, metodo))
            return {k: (d[k] if d[k].ndim else d[k].item()) for k in d.files
                    if k != "firma"}

    if verbose:
        print("  [parseo] %s (%.2f GB, canal %s, metodo %s)..."
              % (os.path.basename(ruta), os.path.getsize(ruta) / 1e9, canal,
                 metodo))
    t0 = time.time()

    campos = ("area", "amplitud", "t_onset", "t_ini", "t_fin", "base", "ruido",
              "fraccion")
    if plantilla is not None:
        campos = campos + ("area_compuerta", "A", "t0", "b_ajuste", "chi2_red",
                           "extrapolada", "t0_en_borde")
    acum = {c: [] for c in campos}
    banderas, relojes, scalers = [], [], []
    n_ev = n_rep = n_ruido = 0
    v_anterior = None
    t_primero = t_ultimo = None

    for t, v, reloj, scaler in drs.iter_eventos(ruta, canal=canal,
                                                    limite=limite,
                                                    callback=callback):
        n_ev += 1
        if reloj:
            if t_primero is None:
                t_primero = reloj
            t_ultimo = reloj
            relojes.append(reloj)
            scalers.append(scaler)
        if t.size < 16:
            continue
        if v_anterior is not None and v.size == v_anterior.size and \
                np.array_equal(v, v_anterior):
            n_rep += 1
            continue
        v_anterior = v
        if v.max() - v.min() <= umbral_ruido:
            n_ruido += 1
            continue

        if plantilla is not None:
            r = pulso.area_evento_plantilla(t, v, plantilla, pre, post, **kw)
        else:
            r = pulso.area_evento(t, v, pre=pre if metodo == "global" else None,
                                  post=post if metodo == "global" else None,
                                  **kw)
        for c in campos:
            acum[c].append(r[c])
        banderas.append(r["bandera"])

    out = {c: np.asarray(acum[c], float) for c in campos}
    out["bandera"] = np.asarray(banderas, np.int16)
    out["scaler"] = np.asarray(scalers, float)
    out["reloj_s"] = _relojes_a_segundos(relojes)
    out["n_eventos"] = n_ev
    out["duracion"] = (out["reloj_s"][-1] - out["reloj_s"][0]
                       if out["reloj_s"].size > 1 else 0.0)
    out["canal"] = canal
    out["metodo"] = metodo
    out["pre"] = float(pre) if pre is not None else np.nan
    out["post"] = float(post) if post is not None else np.nan
    if plantilla is not None:
        out["plantilla"] = str(plantilla.get("nombre", ""))
        out["plantilla_hash"] = str(plantilla["hash"])
        out["plantilla_integral"] = float(plantilla["integral"])
        out["plantilla_t_carga"] = float(plantilla["t_carga"])
    out["descartados_repetidos"] = n_rep
    out["descartados_ruido"] = n_ruido

    if out["area"].size == 0:
        raise ValueError(
            "%s: %d eventos leidos y CERO pulsos utiles en el canal %s. "
            "Revisa el canal (los <Data> tienen que estar dentro de <%s>) y el "
            "umbral de ruido (umbral_ruido=%.1f mV pico a pico)."
            % (os.path.basename(ruta), n_ev, canal, canal, umbral_ruido))

    if verbose:
        marcados = int((out["bandera"] != 0).sum())
        print("  [parseo] %d eventos, %d pulsos utiles (%d repetidos, %d bajo "
              "umbral), %.0f s de pared, en %.0f s"
              % (n_ev, out["area"].size, n_rep, n_ruido, out["duracion"],
                 time.time() - t0))
        print("           %d marcados (%.2f%%): %d saturados, %d truncados, "
              "%d sin pre-pulso, %d sin pico"
              % (marcados, 100 * marcados / max(out["area"].size, 1),
                 int((out["bandera"] & pulso.MAL_SATURA).astype(bool).sum()),
                 int((out["bandera"] & pulso.MAL_TRUNCADO).astype(bool).sum()),
                 int((out["bandera"] & pulso.MAL_SIN_PRE).astype(bool).sum()),
                 int((out["bandera"] & pulso.MAL_SIN_PICO).astype(bool).sum())))

    if usar_cache:
        try:
            np.savez_compressed(archivo_cache, firma=firma, **out)
        except Exception as e:
            print("  [aviso] no pude escribir el cache: %s" % e)
    return out


def _relojes_a_segundos(relojes):
    """
    Convierte las marcas `<Time>` a segundos desde el primer evento.

    Se hace con `datetime` una sola vez por evento (no por muestra), asi que el
    costo es despreciable frente al parseo. Ojo: el DRSOSC escribe estas marcas
    en RAFAGAS —miles de eventos comparten el mismo milisegundo—, asi que sirven
    para medir la DURACION de la corrida pero no para medir tiempo entre
    eventos. El tiempo vivo sale del scaler, en la Etapa 3.
    """
    from datetime import datetime
    if not relojes:
        return np.zeros(0)
    fmt = "%Y/%m/%d %H:%M:%S.%f"
    try:
        t0 = datetime.strptime(relojes[0], fmt)
        return np.array([(datetime.strptime(r, fmt) - t0).total_seconds()
                         for r in relojes], float)
    except (ValueError, TypeError):
        return np.zeros(len(relojes))


# --------------------------------------------------------------------------
# 3. Utilidades de diagnostico
# --------------------------------------------------------------------------
def buenos(d, permitir=pulso.MAL_SIN_PRE):
    """
    Mascara de eventos utilizables.

    `permitir` es el OR de las banderas que NO descartan. Por defecto se acepta
    MAL_SIN_PRE (linea base tomada del post-pulso): el evento sirve, solo tiene
    un pedestal chico. Descartar por ahi seria un corte que depende de la
    posicion del onset, y el onset depende de la amplitud porque el trigger
    tiene walk (los pulsos grandes cruzan el umbral antes) — o sea, otra vez, un
    corte dependiente de la amplitud.
    """
    return (d["bandera"] & ~permitir) == 0


def perfil(x, y, bins=14, minimo=15):
    """
    Mediana de `y` en tramos de `x` (percentiles iguales). Devuelve
    `(centros, medianas, dispersion)`. Es como se mira si una razon es plana
    sin que unos pocos puntos extremos manden.
    """
    bordes = np.percentile(x, np.linspace(0, 100, bins + 1))
    bordes = np.unique(bordes)
    cx, cy, ce = [], [], []
    for j in range(bordes.size - 1):
        m = (x >= bordes[j]) & (x < bordes[j + 1])
        if m.sum() < minimo:
            continue
        cx.append(np.median(x[m])); cy.append(np.median(y[m]))
        ce.append(1.2533 * y[m].std() / np.sqrt(m.sum()))   # error de la mediana
    return np.array(cx), np.array(cy), np.array(ce)


# --------------------------------------------------------------------------
# 4. La forma del pulso en funcion de la amplitud
# --------------------------------------------------------------------------
def amplitud_rapida(t, v):
    """Amplitud de un registro, para poder separar los eventos por tamano."""
    tau, vv = pulso.preparar(t, v)
    return pulso.detectar_pico(tau, vv)[2]


def formas_por_amplitud(ruta, n=1500, cuantiles=(0, 25, 50, 75, 100),
                        canal=None):
    """
    Pulso promedio alineado por onset, separado por CUARTIL DE AMPLITUD.

    Es la pregunta de fondo detras del test de linealidad: si la forma del pulso
    fuera la misma para todos los tamanos, cualquier ventana fija recogeria la
    misma FRACCION de la carga y el area seria proporcional a la energia sin
    mas. Si la forma cambia con la amplitud, ninguna eleccion de ventana lo
    arregla: es electronica del detector, y hay que saberlo ANTES de calibrar,
    no descubrirlo cuando la recta no cierra.

    Devuelve `(tau_rel, formas, etiquetas)`, todas normalizadas a amplitud 1.
    """
    T, V, _, _ = drs.leer_pulsos(ruta, n=n, canal=canal)
    amps = np.array([amplitud_rapida(t, v) for t, v in zip(T, V)])
    cortes = np.percentile(amps[amps > 5], list(cuantiles))
    formas, etiquetas, tau_rel = [], [], None
    for j in range(len(cortes) - 1):
        m = (amps >= cortes[j]) & (amps <= cortes[j + 1])
        if m.sum() < 20:
            continue
        tau_rel, forma, n_us = pulso.pulso_promedio(T[m], V[m], amp_minima=0.0)
        formas.append(forma)
        etiquetas.append("%.0f-%.0f mV (n=%d)" % (cortes[j], cortes[j + 1], n_us))
    return tau_rel, formas, etiquetas


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo(limite=6000):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES if rutas.corrida(k).disponible()]
    fig, axes = plt.subplots(len(claves), 5, figsize=(25, 4.4 * len(claves)),
                             squeeze=False)

    for fila, clave in enumerate(claves):
        c = rutas.corrida(clave)
        print("=" * 78)
        print("%s  -  %s   (primeros %d eventos)" % (clave, c.descripcion, limite))
        print("=" * 78)

        vent = calibrar_ventana(c.muestra)
        dA = areas_de_archivo(c.muestra, metodo="global", pre=vent["pre"],
                              post=vent["post"], limite=limite)
        dB = areas_de_archivo(c.muestra, metodo="convergencia", limite=limite)
        # Ventana ANCHA de referencia. No es un metodo candidato (mete ruido y
        # el rebote de la cola): es la vara contra la cual se mide QUE FRACCION
        # de la carga recoge cada ventana.
        post_ancho = 130.0
        dR = areas_de_archivo(c.muestra, metodo="global", pre=vent["pre"],
                              post=post_ancho, limite=limite)

        mA, mB, mR = buenos(dA), buenos(dB), buenos(dR)
        m = mA & mB & (dA["amplitud"] > 5.0)
        a, b = dA["area"][m], dB["area"][m]
        print("  %d eventos comparables" % m.sum())
        rel = (b - a) / np.maximum(a, 1e-9)
        print("  area(convergencia)/area(global): mediana %.3f  ->  el metodo por "
              "convergencia pierde %.1f%% del area"
              % (np.median(b / a), -100 * np.median(rel)))

        # --- TEST 1: la fraccion de carga recogida, contra una ventana ancha --
        # El eje x es la AMPLITUD y no un area, y eso no es casual: si se binea
        # por el area que aparece en el cociente, el bin mas bajo se llena
        # preferentemente con eventos cuya area de referencia fluctuo hacia
        # abajo, y el cociente sale inflado por pura seleccion (probado: daba
        # una "deriva" de -41% que era enteramente artefacto). La amplitud es
        # independiente de la eleccion de ventana.
        #
        # El corte en amplitud tambien hace falta: la amplitud es el MAXIMO de
        # un registro ruidoso, asi que en pulsos comparables al ruido esta
        # sesgada hacia arriba. Por debajo de ~20*ruido el eje x deja de ser
        # confiable, y ese es el unico lugar donde el test no dice nada.
        amp_min = 20.0 * float(np.median(dA["ruido"]))
        mr = mA & mB & mR & (dR["area"] > 0) & (dA["amplitud"] > amp_min)
        ar_ref = dR["area"][mr]
        amp_r = dA["amplitud"][mr]
        print("  TEST DE LINEALIDAD (fraccion recogida contra una ventana ancha "
              "de %.0f ns; tiene que ser PLANA)" % post_ancho)
        print("     %d eventos con la ventana ancha completa dentro del registro "
              "y amplitud > %.0f mV (20 sigma de ruido)"
              % (int(mr.sum()), amp_min))
        derivas = {}
        for nombre, dd in (("global", dA), ("convergencia", dB)):
            razon = dd["area"][mr] / ar_ref
            x, y, e = perfil(amp_r, razon)
            if x.size < 4:
                continue
            derivas[nombre] = 100 * (y[-1] - y[0]) / y[0]
            print("     %-12s : recoge %.1f%% de la carga; entre amplitud %.0f y "
                  "%.0f mV va de %.3f a %.3f  ->  deriva %+6.1f%%"
                  % (nombre, 100 * np.median(razon), x[0], x[-1], y[0], y[-1],
                     derivas[nombre]))

        # Un cociente que se aparta de 1 no dice POR QUE. La diferencia si:
        #   constante con la amplitud  -> es un pedestal (error de linea base
        #       multiplicado por los ns de mas que integra la ventana ancha), y
        #       el problema es de la REFERENCIA, no de la ventana que se juzga.
        #   proporcional a la amplitud -> es forma: la ventana recorta una
        #       fraccion distinta segun el tamano del pulso, y eso si rompe la
        #       proporcionalidad area-energia.
        dif = ar_ref - dA["area"][mr]
        xd, yd, ed = perfil(amp_r, dif)
        c1, c0 = np.polyfit(xd, yd, 1)          # dif = c0 + c1 * amplitud
        amp_tip = float(np.median(amp_r))
        area_tip = float(np.median(ar_ref))
        print("     de donde viene: area(ancha) - area(global) contra amplitud")
        print("        amp %5.0f mV -> %+7.1f mV*ns        "
              "amp %5.0f mV -> %+7.1f mV*ns" % (xd[0], yd[0], xd[-1], yd[-1]))
        print("        ajuste  dif = %+.1f mV*ns  %+.4f * amplitud" % (c0, c1))
        print("           termino CONSTANTE : %+.1f mV*ns  -> pedestal de la "
              "ventana ancha (error de linea base x los ns de mas que integra)."
              % c0)
        print("                               son %+.2f mV de offset sobre los "
                "%.0f ns extra, y %+.1f%% del area tipica"
              % (c0 / max(post_ancho - vent["post"], 1e-9),
                 post_ancho - vent["post"], 100 * c0 / area_tip))
        print("           termino PROPORCIONAL: %+.1f%% del area tipica -> esto "
              "si es forma: cola que la ventana no recoge (o rebote que la ancha "
              "si integra)." % (100 * c1 * amp_tip / area_tip))

        # --- TEST 2: cambia la FORMA del pulso con la amplitud? ---------------
        tau_rel, formas, etiquetas = formas_por_amplitud(c.muestra)
        fr = []
        for f in formas:
            _, ff, _ = pulso.curva_convergencia(tau_rel, f)
            i_v = min(int(np.searchsorted(tau_rel, vent["post"])), ff.size - 1)
            fr.append(float(ff[i_v]))
        if len(fr) >= 2:
            print("  FORMA DEL PULSO por cuartil de amplitud: fraccion de la carga "
                  "que cae dentro de la ventana de %.0f ns"
                  % (vent["pre"] + vent["post"]))
            for et, f in zip(etiquetas, fr):
                print("     %-22s -> %6.2f %%" % (et, 100 * f))
            print("     variacion extremo a extremo: %+.2f%%" % (100 * (fr[-1] - fr[0])))

        print("  percentiles del area [mV*ns]  (p50 / p90 / p99 / p99.9):")
        for nombre, ar in (("global      ", a), ("convergencia", b)):
            print("     %s : %8.1f %8.1f %8.1f %8.1f"
                  % (nombre, *np.percentile(ar, [50, 90, 99, 99.9])))
        print()

        # ------------------------------------------------------------ figura --
        ax = axes[fila][0]
        ax.scatter(a, b, s=3, alpha=0.2)
        lim = float(np.percentile(a, 99.5))
        ax.plot([0, lim], [0, lim], "k--", lw=1, label="identidad")
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_title("%s: area global vs area por convergencia" % clave)
        ax.set_xlabel("area GLOBAL [mV*ns]")
        ax.set_ylabel("area CONVERGENCIA [mV*ns]")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[fila][1]
        ax.hist(100 * rel, bins=120, range=(-60, 20), color="C1")
        ax.axvline(100 * np.median(rel), color="k", ls="--",
                   label="mediana %.1f%%" % (100 * np.median(rel)))
        ax.set_title("cuanto se pierde por evento")
        ax.set_xlabel("(convergencia - global) / global  [%]")
        ax.set_ylabel("eventos"); ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[fila][2]
        for nombre, dd, col in (("global", dA, "C0"), ("convergencia", dB, "C3")):
            x, y, e = perfil(amp_r, dd["area"][mr] / ar_ref)
            ax.errorbar(x, y, yerr=e, marker="o", ms=3, lw=1, color=col,
                        label="%s (deriva %+.1f%%)" % (nombre, derivas.get(nombre, np.nan)))
        ax.axhline(1.0, color="k", ls=":", lw=1)
        ax.set_title("TEST: fraccion de carga recogida (tiene que ser PLANA)")
        ax.set_xlabel("amplitud [mV]")
        ax.set_ylabel("area(ventana) / area(ancha de %.0f ns)" % post_ancho)
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        ax = axes[fila][3]
        for f, et in zip(formas, etiquetas):
            ax.plot(tau_rel, f, lw=1.2, label=et)
        ax.axvline(vent["post"], color="k", ls="--", lw=1,
                   label="t_fin = %.0f ns" % vent["post"])
        ax.set_title("forma normalizada por cuartil de amplitud")
        ax.set_xlabel("t - onset [ns]"); ax.set_ylabel("v / amplitud")
        ax.set_ylim(-0.05, 0.30)
        ax.legend(fontsize=6); ax.grid(alpha=0.3)

        ax = axes[fila][4]
        rango = (0.0, float(np.percentile(a, 99.5) * 1.6))
        hA, bordes = np.histogram(a, bins=400, range=rango)
        hB, _ = np.histogram(b, bins=400, range=rango)
        ch = 0.5 * (bordes[:-1] + bordes[1:])
        ax.step(ch, hA, where="mid", lw=1.1, label="global")
        ax.step(ch, hB, where="mid", lw=1.1, label="convergencia")
        ax.set_yscale("log")
        ax.set_title("el espectro que sale de cada uno")
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

    fig.tight_layout()
    salida = rutas.figura("02_comparacion_metodos.png")
    fig.savefig(salida, dpi=125)
    print("figura -> %s" % rutas.relativa(salida))


def listar_cache():
    """
    Qué hay en `.cache/`. Los nombres son hashes, asi que sin esto no hay manera
    de saber que guardo cada uno ni cual borrar.
    """
    import glob
    archivos = sorted(glob.glob(os.path.join(rutas.CACHE, "areas_*.npz")))
    if not archivos:
        print("no hay areas cacheadas todavia")
        return
    print("%d archivos de areas en %s" % (len(archivos), rutas.relativa(rutas.CACHE)))
    for f in archivos:
        try:
            d = np.load(f, allow_pickle=False)
            partes = str(d["firma"]).split("|")
            print("  %-22s %6d pulsos   %s"
                  % (os.path.basename(f), d["area"].size,
                     " ".join(p for p in partes[3:] if "=" in p or p.isalpha())))
            print("  %-22s %s" % ("", os.path.basename(partes[0])))
        except Exception as e:
            print("  %-22s ilegible (%s)" % (os.path.basename(f), e))


if __name__ == "__main__":
    if "--cache" in sys.argv:
        rutas.salida_utf8()
        listar_cache()
    else:
        n = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
        _demo(limite=n)
