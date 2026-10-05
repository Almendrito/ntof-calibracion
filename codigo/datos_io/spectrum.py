"""
spectrum.py — Etapa 0: del archivo crudo a un espectro (canal, cuentas) limpio.

Lee los DOS formatos del DRSOSC sin que haya que decirle cual es: el XML (junio
y 26-ago-2026) y el binario del DRS4 (1-sep-2026, sin extension en el nombre).
El formato se decide por los primeros 4 bytes (`formato`/`es_binario`) y el
binario se delega en `analisis/drs_bin.py`, que es donde vive el parser byte a
byte; aca no hay una segunda copia de ese formato.

Refactor de lecturaxml.py con dos arreglos clave:
  1. BINNING COMUN: muestra y fondo se histogramean con EL MISMO range fijo,
     asi comparten bordes de bin y la resta de fondo tiene sentido fisico.
  2. LINEA BASE: se resta el baseline de cada pulso antes de integrar el area,
     para que el offset no meta ancho artificial en el espectro.

Conserva de tu codigo: dedup de tramas consecutivas, umbral de ruido
peak-to-peak, y normalizacion por tiempo de adquisicion.
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/datos_io/spectrum.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------


import os
import time
import numpy as np
from datetime import datetime

_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy nuevo o viejo
_FMT_TIEMPO = "%Y/%m/%d %H:%M:%S.%f"


# ==========================================================================
# De que FORMATO es el archivo
# ==========================================================================
# El DRSOSC guarda de dos maneras y en este proyecto hay de las dos:
#
#   - XML     (junio 2026 y la sesion del 26-ago-2026)
#   - binario ("Save as binary", la sesion del 1-sep-2026)
#
# Los binarios de la sesion del 1-sep NO TIENEN EXTENSION en el nombre (se
# llaman `Co60`, `Cs137`, `B1700V`), asi que decidir por el nombre habria sido
# decidir por nada: el formato se mira en los primeros 4 BYTES del archivo.
#
# El parseo byte a byte NO se duplica aca. Vive en `analisis/drs_bin.py`, que es
# donde se dedujo el formato y se verifico contra los archivos de la sesion. Un
# segundo parser del mismo formato es exactamente la clase de cosa que se
# desincroniza en silencio, que es la moraleja de CLAUDE.md 7.36.
MAGIA_BIN = (b"DRS2", b"DRS8")


def es_binario(ruta_archivo):
    """True si el archivo es el binario del DRS4 en vez del XML del DRSOSC."""
    try:
        with open(ruta_archivo, "rb") as f:
            return f.read(4) in MAGIA_BIN
    except OSError:
        return False


def formato(ruta_archivo):
    """'binario' o 'xml'. Mira el CONTENIDO, no el nombre."""
    return "binario" if es_binario(ruta_archivo) else "xml"


def _drs_bin():
    """
    El lector binario, importado tarde.

    Tarde a proposito: asi `codigo/` sigue corriendo entero sobre los XML aunque
    `analisis/` no este, y solo pide el modulo cuando de verdad hay un binario
    que leer.
    """
    _RAIZ = _os.path.dirname(_CODIGO)                # .../Calibracion CE
    if _RAIZ not in _sys.path:
        _sys.path.insert(0, _RAIZ)
    from analisis import drs_bin                     # noqa: PLC0415
    return drs_bin


def pares_a_pulso(pares):
    """
    `(t, v)` a partir de lo que entrega `_iter_eventos`, venga del formato que venga.

    Del XML llega una lista de cadenas "t,v" (una por muestra) y del binario un
    array (n, 2) ya numerico. Ningun consumidor tiene por que saber cual de los
    dos le toco, asi que la conversion se hace en un solo lugar.

    Devuelve `(None, None)` si el evento no trae un pulso utilizable.
    """
    if pares is None:
        return None, None
    if isinstance(pares, np.ndarray):
        arr = pares
    else:
        if not pares:
            return None, None
        arr = np.fromstring(",".join(pares), sep=",")
        if arr.size < 4:
            return None, None
        arr = arr[:(arr.size // 2) * 2].reshape(-1, 2)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return None, None
    return arr[:, 0], arr[:, 1]


def _area_pulso(t, v, baseline_frac=0.2):
    """
    LEGACY - se conserva solo para poder REPRODUCIR el analisis viejo.

    Area del pulso con la linea base estimada sobre el primer `baseline_frac`
    del registro (20% = 40 ns) y la integral sobre TODO el registro.

    ⚠️ ESTA DEFINICION ESTA ROTA para estos datos, ver `_area_pulso_gate` y
    CLAUDE.md 7.33: el pulso arranca a ~25 ns y el pico esta a ~48 ns, asi que
    la ventana de linea base SE COME LA SUBIDA del pulso. Medido sobre 3639
    eventos limpios, integral/amplitud [ns]:

        base 0-12 ns : mediana -8.48   std  0.84   eventos con signo + :   0
        base 0-40 ns : mediana -7.71   std 11.33   eventos con signo + : 954 (26%)

    El `abs()` de abajo recicla esos 954 como si fueran buenos.
    """
    n_base = max(1, int(len(v) * baseline_frac))
    baseline = np.mean(v[:n_base])
    return abs(_TRAPZ(v - baseline, x=t))


# ---- parametros del integrador con compuerta (en ns; ver CLAUDE.md 7.33) ----
GATE_PRE = 8.0        # cuanto integrar ANTES del onset
GATE_POST = 55.0      # cuanto integrar DESPUES del onset (el pulso muere a ~85 ns)
GUARDA_BASE = 4.0     # margen que se le saca a la linea base alrededor del pulso
V_SATURA = 495.0      # |v| de saturacion del DRS4 [mV]
FRAC_ONSET = 0.10     # fraccion de la amplitud que define el onset


def _area_pulso_gate(t, v, pre=None, post=None, guarda=None,
                     frac_onset=None, v_satura=None):
    """
    Area del pulso con linea base FUERA del pulso y compuerta de integracion.

    Los arreglos respecto de `_area_pulso`, en orden de importancia:

    1. LINEA BASE FUERA DEL PULSO. Se detecta el onset (primer cruce del 10% de
       la amplitud retrocediendo desde el pico) y la linea base es la MEDIANA de
       las muestras PRE-onset (hasta `guarda` ns antes de que abra la
       compuerta); si no hay al menos 10, se cae al complemento de la compuerta
       y se marca el evento. Con la ventana fija del 20% del registro el pulso
       caia dentro de su propia linea base: el sesgo `b -> b + q_b/T_b` se
       multiplica por `T/T_b ~ 5` al integrar todo el registro, y como el jitter
       del trigger mueve el pico entre 18 y 54 ns, cada evento perdia una
       fraccion DISTINTA de su carga. Eso es un ensanchamiento multiplicativo
       puro, que es lo que inflaba la resolucion al 43% FWHM.
    2. COMPUERTA. Se integra `[onset-pre, onset+post]` en vez de los 200 ns. El
       pulso medio devuelve el 100% de su area en 80 ns, asi que los 120 ns
       restantes solo aportan ruido y deriva.
    3. TIEMPO RELATIVO AL REGISTRO. El DRSOSC escribe el eje de tiempo relativo
       a la celda de trigger, y en ~2% de los eventos arranca corrido (t[0] entre
       -78 y +121 ns). Todo se mide contra `t - t[0]`, nunca contra t absoluto.
    4. MEDIANA en vez de media: inmune a un dark count o a un pile-up aislado.

    Que NO se descarta [importante]: la version anterior de esta funcion exigia
    que el onset cayera despues de 10 ns para tener linea base, y eso rechazaba
    entre 4.7% y 12.2% de los eventos SEGUN LA AMPLITUD (el trigger tiene walk:
    los pulsos grandes cruzan el umbral antes). Un corte que depende de la
    amplitud deforma el espectro, que es exactamente lo que no se puede hacer en
    una calibracion. Cayendo al complemento de la compuerta cuando no hay
    pre-pulso, el corte deja de hacer falta.

    Devuelve (area, amplitud, t_onset, bandera). bandera = 0 si el evento es
    bueno; si no, un bit por problema:
        2 = saturacion del digitalizador
        4 = el registro arranca DENTRO del pulso (queda truncado)
        8 = no habia pre-pulso suficiente y la linea base salio del complemento
            de la compuerta. Medido sobre 7379 eventos con las dos disponibles,
            la mediana del post queda 0 o 1 LSB (0.1 mV) por encima de la del
            pre: sobre la compuerta de 63 ns eso es un pedestal de ~6 mV*ns, o
            sea 0.5% en el borde. Como el mapa canal-energia pasa por el origen,
            ese pedestal se traduce casi 1:1 en D. Se marca pero NO se descarta.
    """
    # Los defaults se leen de las globales EN CADA LLAMADA, no se congelan en la
    # firma: si no, cambiar `spectrum.GATE_POST` desde afuera no tendria efecto y
    # los barridos de sistematicos darian todos el mismo numero (pasa de verdad).
    pre = GATE_PRE if pre is None else pre
    post = GATE_POST if post is None else post
    guarda = GUARDA_BASE if guarda is None else guarda
    frac_onset = FRAC_ONSET if frac_onset is None else frac_onset
    v_satura = V_SATURA if v_satura is None else v_satura

    n = len(v)
    if n < 16:
        return 0.0, 0.0, 0.0, 4
    tau = t - t[0]                                # tiempo relativo al registro

    # linea base provisoria: mediana de TODO el registro. El pulso ocupa ~30% de
    # las muestras, asi que la mediana no lo ve.
    vc = v - np.median(v)

    i_pico = int(np.argmax(np.abs(vc)))
    signo = -1.0 if vc[i_pico] < 0 else 1.0
    amp = abs(float(vc[i_pico]))

    # onset: ultimo cruce del 10% de la amplitud ANTES del pico
    idx = np.flatnonzero(np.abs(vc[:i_pico + 1]) < frac_onset * amp)
    i_on = int(idx[-1]) if idx.size else 0
    t_on = float(tau[i_on])

    bandera = 0
    if np.max(np.abs(v)) >= v_satura:
        bandera |= 2

    ini = int(np.searchsorted(tau, t_on - pre))
    fin = int(np.searchsorted(tau, t_on + post))
    ini_c, fin_c = max(ini, 0), min(fin, n - 1)
    if fin_c - ini_c < 4:
        return 0.0, amp, t_on, bandera | 4

    # linea base definitiva: mediana del COMPLEMENTO de la compuerta (con guarda)
    g = int(round(guarda / max(float(np.median(np.diff(tau))), 1e-9)))
    fuera = np.ones(n, bool)
    fuera[max(ini_c - g, 0):min(fin_c + g, n)] = False
    # La linea base sale del PRE-pulso, no de todo el complemento: la zona
    # posterior al pulso queda ~1 LSB por encima (rebote de polaridad opuesta),
    # y sobre la compuerta eso es un pedestal de ~6 mV*ns que el mapa por el
    # origen convierte en 0.5% de D. Si no hay pre-pulso suficiente se usa el
    # complemento y se marca, pero NO se descarta: descartar por posicion del
    # onset seria un corte dependiente de la amplitud (el trigger tiene walk).
    pre_ok = np.zeros(n, bool)
    pre_ok[:max(ini_c - g, 0)] = True
    # 10 muestras alcanzan: la mediana de 10 tiene sigma = 0.42/sqrt(10) = 0.13 mV,
    # que sobre la compuerta de 63 ns son 8 mV*ns contra un borde de 1229.
    if pre_ok.sum() >= 10:
        muestras = v[pre_ok]
    else:
        bandera |= 8
        muestras = v[fuera] if fuera.sum() >= 20 else v
    base = float(np.median(muestras))
    # ruido robusto (1.4826*MAD = sigma de una gaussiana)
    ruido = max(1.4826 * float(np.median(np.abs(muestras - base))), 1e-3)

    # Truncamiento: el registro arranca DENTRO del pulso. El umbral lleva un piso
    # de ruido, si no los pulsos chicos se auto-marcan (con amp = 10 mV, el 10%
    # son 1 mV y el rms es 0.42: se marcaria el 40% por fluctuacion).
    if fin >= n or abs(float(v[0] - base)) > max(frac_onset * amp, 5.0 * ruido):
        bandera |= 4
    ini, fin = ini_c, fin_c

    area = signo * _TRAPZ(v[ini:fin + 1] - base, x=t[ini:fin + 1])
    return float(max(area, 0.0)), amp, t_on, bandera


def detectar_canal(ruta_archivo, max_lineas=200000):
    """
    Averigua en QUE canal esta el pulso, mirando el primer evento del XML.

    Por que hace falta: el canal no es una propiedad del analisis sino de como
    estaba enchufado el detector ese dia. Las corridas de junio (Co-60) salieron
    por CHN2 con <Scaler1>; las del 26-ago-2026 (Cs-137, Detector A a 1600 V)
    salieron por CHN1 con <Scaler0>. Con el canal cableado a CHN2, esas corridas
    daban CERO pulsos sin decir por que: el escaner recorria el archivo entero
    sin encontrar un solo <Data> adentro del tag que buscaba.

    El scaler del DRSOSC va indexado desde 0, asi que a CHN{n} le corresponde
    <Scaler{n-1}>; si ese tag no esta, se usa el primer <ScalerN> que aparezca.

    En el BINARIO no hay tags que buscar: los canales estan declarados en la
    cabecera del archivo, asi que la pregunta se le hace a `drs_bin` y la
    respuesta se devuelve con el mismo vocabulario ("CHN1", "Scaler0") para que
    la clave del cache y los informes no tengan que distinguir formatos.

    Devuelve (canal, tag_scaler), p.ej. ("CHN1", "Scaler0").
    """
    if es_binario(ruta_archivo):
        return _drs_bin().detectar_canal(ruta_archivo)

    import re
    canales, scalers, canal_actual = [], [], None
    with open(ruta_archivo, "r", encoding="ISO-8859-1", errors="replace") as f:
        for k, linea in enumerate(f):
            s = linea.strip()
            if s.startswith("<Data>"):
                if canal_actual and canal_actual not in canales:
                    canales.append(canal_actual)    # ese canal SI trae datos
                continue
            if s.startswith("</Event>") or k > max_lineas:
                break
            m = re.match(r"<(CHN\d+)>", s)
            if m:
                canal_actual = m.group(1)
                continue
            if s.startswith("</CHN"):
                canal_actual = None
                continue
            m = re.match(r"<(Scaler\d+)>", s)
            if m:
                scalers.append(m.group(1))
    if not canales:
        raise ValueError(
            f"{os.path.basename(ruta_archivo)}: no encontre ningun <CHNn> con "
            f"<Data> en el primer evento. Es un XML del DRSOSC?")
    canal = canales[0]
    esperado = "Scaler%d" % (int(canal[3:]) - 1)
    tag = esperado if esperado in scalers else (scalers[0] if scalers else None)
    return canal, tag


def _iter_eventos(ruta_archivo, canal=None, tag_scaler=None, callback=None,
                  cada=400):
    """
    Recorre el XML como TEXTO y entrega (time_str, pares_data, scaler) por evento.

    Por que no ElementTree: el archivo de Co-60 pesa 1.76 GB y tiene ~5e7 tags
    <Data>. ET.parse construye un objeto Element por tag -> decenas de GB de RAM
    y varios minutos. El formato del DRSOSC es rigido (un tag por linea), asi que
    un escaneo lineal da lo mismo con RAM plana y mucho mas rapido.

    Estructura:  DRSOSC > Event > Board_XXXX > CHN{n} > Data ("t,v")
    Solo se leen los <Data> DENTRO del canal del detector. Si `canal` es None se
    detecta solo con `detectar_canal`: CHN2 en las corridas de Co-60 de junio,
    CHN1 en las del 26-ago-2026.

    callback : callable(fraccion, n_eventos) llamado cada `cada` eventos. Este
        bucle es el unico paso del proyecto que tarda ~40 s sin decir nada, asi
        que la GUI necesita poder mostrar por donde va. La fraccion se mide en
        BYTES leidos y no en eventos, porque el numero total de eventos no se
        sabe hasta terminar mientras que el tamano del archivo si.

        Ojo con el conteo: el archivo se abre en modo texto, asi que Python
        traduce CRLF a LF y `len(linea)` queda un byte corto por linea. Se
        detecta el fin de linea real con `f.newlines` en el primer evento y se
        compensa; sin eso la barra llegaria al 78% y ahi terminaria.

    Con un archivo BINARIO nada de esto aplica y el trabajo se delega en
    `drs_bin`, que devuelve el pulso ya convertido a numeros. Para que el
    consumidor no tenga que ramificar, lo que se entrega en el lugar de `pares`
    es un array (n, 2) en vez de una lista de cadenas: `pares_a_pulso` sabe
    leer los dos.
    """
    if es_binario(ruta_archivo):
        yield from _iter_eventos_bin(ruta_archivo, canal=canal,
                                     callback=callback, cada=cada)
        return

    if canal is None:
        canal, tag_auto = detectar_canal(ruta_archivo)
        if tag_scaler is None:
            tag_scaler = tag_auto
    if tag_scaler is None:
        tag_scaler = "Scaler%d" % (int(canal[3:]) - 1)
    ab, ce = "<%s>" % canal, "</%s>" % canal
    ab_sc, n_sc = "<%s>" % tag_scaler, len(tag_scaler) + 2

    pares, time_str, en_canal, scaler = [], None, False, np.nan
    total = max(os.path.getsize(ruta_archivo), 1)
    leidos, extra_nl, n_ev = 0, 0, 0
    with open(ruta_archivo, "r", encoding="ISO-8859-1", errors="replace") as f:
        for linea in f:
            leidos += len(linea) + extra_nl
            s = linea.strip()
            if s.startswith("<Data>"):
                if en_canal:
                    pares.append(s[6:-7])           # quita <Data> y </Data>
            elif s.startswith(ab):
                en_canal = True
            elif s.startswith(ce):
                en_canal = False
            elif s.startswith("<Time>"):
                time_str = s[6:-7]
            elif s.startswith(ab_sc):
                # tasa de trigger REAL de la placa [Hz]. Es la clave del tiempo
                # vivo: el DRSOSC escribe XML por evento y no alcanza a guardar
                # todos los triggers (ver `tiempo_vivo`).
                try:
                    scaler = float(s[n_sc:-(n_sc + 1)])
                except ValueError:
                    scaler = np.nan
            elif s.startswith("</Event>"):
                n_ev += 1
                if n_ev == 1:
                    nl = f.newlines
                    nl = nl if isinstance(nl, str) else (nl[0] if nl else "\n")
                    extra_nl = len(nl) - 1
                if callback is not None and n_ev % cada == 0:
                    callback(min(leidos / total, 1.0), n_ev)
                yield time_str, pares, scaler
                pares, time_str, scaler = [], None, np.nan


def _iter_eventos_bin(ruta_archivo, canal=None, callback=None, cada=400):
    """
    Lo mismo que `_iter_eventos` pero sobre el binario del DRS4.

    Entrega `(reloj, arr, scaler)` con `arr` de forma (n, 2) — columnas t [ns] y
    v [mV] — que es lo que `pares_a_pulso` espera del lado binario.

    Dos cosas que el binario da gratis y el XML no:
      * el numero de eventos se sabe de la cabecera, asi que la fraccion del
        callback es exacta desde el primer evento (en el XML hay que medirla en
        bytes leidos porque el total de eventos no se conoce hasta terminar);
      * el eje temporal viene del ancho real de cada celda del DRS4 y de la
        celda de trigger del evento, o sea NO es uniforme, igual que en el XML.

    Ojo con el scaler: las versiones viejas del DRSOSC no lo escriben. Cuando
    falta, `drs_bin` devuelve nan y `tiempo_vivo` se cae solo al tiempo de
    pared (supone 100% vivo), que es lo unico honesto sin ese dato.
    """
    drs_bin = _drs_bin()
    for t, v, reloj, scaler in drs_bin.iter_eventos(ruta_archivo, canal=canal,
                                                    callback=callback,
                                                    cada=cada):
        yield reloj, np.column_stack((t, v)), scaler


_CODIGO_METODO = {"gate": 1.0, "legacy": 0.0}
_VERSION_GATE = 6.0        # subir esto invalida el cache de areas a proposito


def tiempo_vivo(relojes, scalers, duracion, n_eventos):
    """
    Tiempo VIVO de la adquisicion [s], a partir del scaler de la placa.

    Por que hace falta [fuera del paper]
    ------------------------------------
    El DRSOSC escribe un XML de ~35 kB por evento, y a tasa alta NO alcanza a
    guardar todos los triggers. El `<Scaler1>` de cada evento es la tasa de
    trigger REAL de la placa [Hz], asi que

        triggers verdaderos ~ integral de Scaler1 dt      (todos los que hubo)
        fraccion viva       = eventos guardados / triggers verdaderos
        T_vivo              = fraccion viva * duracion de pared

    Medido en estos datos:

        Co-60 : 71.5 ev/s guardados, scaler ~170 Hz -> fraccion viva 0.40
        fondo : 10.1 ev/s guardados, scaler ~ 10 Hz -> fraccion viva 1.00

    O sea la corrida de Co-60 estuvo VIVA solo 283 s de los 699 s de pared, y la
    de fondo practicamente los 986 s. Normalizar la resta de fondo por el tiempo
    de PARED (lo que se hacia antes) sobre-resta el fondo por un factor 2.5.

    El scaler tiene rafagas absurdas (hasta 90 kHz en la corrida de fondo, ruido
    electronico), asi que se recorta a 10x la mediana antes de integrar.
    """
    relojes = np.asarray(relojes, float)
    scalers = np.asarray(scalers, float)
    ok = np.isfinite(scalers) & (scalers > 0)
    if relojes.size < 2 or ok.sum() < 10 or duracion <= 0:
        return float(duracion)                    # sin scaler: se asume 100% vivo

    med = float(np.median(scalers[ok]))
    s_lim = np.clip(np.nan_to_num(scalers, nan=med), 0.0, 10.0 * max(med, 1.0))
    dt = np.clip(np.diff(relojes, prepend=relojes[0]), 0.0, None)
    n_true = float(np.sum(s_lim * dt))
    if n_true <= 0:
        return float(duracion)
    frac = min(n_eventos / n_true, 1.0)            # >1 es cuantizacion del scaler
    return float(frac * duracion)


def _clave_cache(ruta_archivo, umbral_ruido, baseline_frac, metodo="legacy",
                 canal=None):
    """
    Huella del archivo + parametros: si algo cambia, el cache se invalida.

    El CANAL entra en la clave: un XML leido del canal equivocado deja un cache
    de cero pulsos que, sin esto, se reusaria para siempre (paso con las
    corridas del 26-ago-2026, que estan en CHN1). Agregarlo alarga el vector, o
    sea que la primera corrida despues de este cambio re-parsea todo una vez.
    """
    st = os.stat(ruta_archivo)
    n_canal = float(canal[3:]) if canal else 0.0
    return np.array([st.st_size, int(st.st_mtime), umbral_ruido, baseline_frac,
                     _CODIGO_METODO.get(metodo, -1.0), GATE_PRE, GATE_POST,
                     GUARDA_BASE, FRAC_ONSET, _VERSION_GATE, n_canal], float)


def leer_areas_y_tiempo(ruta_archivo, umbral_ruido, baseline_frac=0.2,
                        usar_cache=True, verbose=True, metodo="gate",
                        devolver_extra=False, canal=None, callback=None):
    """
    Lee un XML y devuelve (areas_pulsos, duracion_segundos).
    NO histograma todavia: eso se hace despues con un range COMUN.

    usar_cache: guarda las areas en '<archivo>.areas.npz'. El parseo del XML es
    el paso caro (minutos); el analisis posterior se itera muchas veces, asi que
    conviene pagarlo UNA vez. El cache se invalida solo si cambia el tamano o la
    fecha del XML, o el umbral/baseline con que se calcularon las areas.

    callback : callable(fraccion, n_eventos) para reportar avance del parseo.
        Solo se llama cuando de verdad hay que parsear; si el cache sirve, esta
        funcion devuelve al instante y nunca avisa nada.
    """
    if canal is None:
        canal, tag_scaler = detectar_canal(ruta_archivo)
    else:
        tag_scaler = None
    sufijo = ".areas.npz" if metodo == "legacy" else ".areas_gate.npz"
    cache = ruta_archivo + sufijo
    clave = _clave_cache(ruta_archivo, umbral_ruido, baseline_frac, metodo,
                         canal)

    if usar_cache and os.path.exists(cache):
        try:
            d = np.load(cache)
            if d["clave"].shape == clave.shape and np.allclose(d["clave"], clave):
                if verbose:
                    print(f"  [cache] {os.path.basename(ruta_archivo)}: "
                          f"{len(d['areas'])} pulsos ({metodo})")
                if devolver_extra:
                    return (d["areas"], float(d["duracion"]),
                            {k: d[k] for k in ("amp", "t_on", "bandera",
                                               "t_vivo", "scaler")
                             if k in d.files})
                return d["areas"], float(d["duracion"])
        except Exception:
            pass                              # cache corrupto -> se reconstruye

    if verbose:
        mb = os.path.getsize(ruta_archivo) / 1e6
        print(f"  [parseo] {os.path.basename(ruta_archivo)} ({mb:.0f} MB, "
              f"{formato(ruta_archivo)}, canal {canal})...")
    t0 = time.time()

    areas, amps, t_ons, banderas = [], [], [], []
    relojes, scalers = [], []
    pulso_anterior = None
    t_ini = t_fin = None
    n_ev = 0

    for time_str, pares, scaler in _iter_eventos(ruta_archivo, canal,
                                                 tag_scaler,
                                                 callback=callback):
        n_ev += 1
        if time_str:                          # duracion = reloj del 1er al ultimo trigger
            ts = datetime.strptime(time_str, _FMT_TIEMPO)
            if t_ini is None:
                t_ini = ts
            t_fin = ts
            relojes.append((ts - t_ini).total_seconds())
            scalers.append(scaler)
        # `pares` es una lista de cadenas "t,v" si el archivo es XML y un array
        # (n, 2) si es binario. `pares_a_pulso` unifica los dos y devuelve
        # (None, None) cuando el evento no trae pulso.
        t, v = pares_a_pulso(pares)
        if t is None:
            continue

        if pulso_anterior is not None and np.array_equal(v, pulso_anterior):
            continue                          # descarta trama repetida consecutiva
        pulso_anterior = v

        if v.max() - v.min() <= umbral_ruido:  # umbral peak-to-peak (ruido plano)
            continue
        if metodo == "legacy":
            areas.append(_area_pulso(t, v, baseline_frac))
            amps.append(np.nan); t_ons.append(np.nan); banderas.append(0)
        else:
            a, amp, t_on, bandera = _area_pulso_gate(t, v)
            areas.append(a); amps.append(amp)
            t_ons.append(t_on); banderas.append(bandera)

    areas = np.array(areas)
    amps = np.array(amps, float)
    t_ons = np.array(t_ons, float)
    banderas = np.array(banderas, np.int16)
    relojes = np.array(relojes, float)
    scalers = np.array(scalers, float)
    duracion = (t_fin - t_ini).total_seconds() if (t_ini and t_fin) else 0.0
    t_vivo = tiempo_vivo(relojes, scalers, duracion, n_ev)
    if verbose:
        malos = int((banderas != 0).sum())
        print(f"  [parseo] {n_ev} eventos, {len(areas)} pulsos utiles, "
              f"{duracion:.0f} s de pared / {t_vivo:.0f} s vivos, "
              f"en {time.time()-t0:.0f} s"
              + (f"   ({malos} marcados, {100*malos/max(len(areas),1):.1f}%)"
                 if metodo != "legacy" else ""))

    if len(areas) == 0:
        # Sin esto se guardaba un cache de cero pulsos y el error se propagaba
        # callado hasta un histograma de ceros varias etapas mas abajo.
        donde = ("los <Data> tienen que estar dentro de <%s>" % canal
                 if not es_binario(ruta_archivo)
                 else "la cabecera del binario declara %s" % canal)
        raise ValueError(
            f"{os.path.basename(ruta_archivo)}: {n_ev} eventos leidos pero 0 "
            f"pulsos utiles en el canal {canal} ({formato(ruta_archivo)}). "
            f"Revisa el canal ({donde}) y el umbral de ruido "
            f"(umbral_ruido={umbral_ruido} mV pico a pico).")

    if usar_cache:
        try:
            np.savez(cache, areas=areas, duracion=duracion, clave=clave,
                     amp=amps, t_on=t_ons, bandera=banderas,
                     t_vivo=t_vivo, reloj=relojes, scaler=scalers)
        except Exception as e:
            print(f"  [aviso] no se pudo escribir el cache: {e}")

    if devolver_extra:
        return areas, duracion, {"amp": amps, "t_on": t_ons, "bandera": banderas,
                                 "t_vivo": t_vivo, "scaler": scalers}
    return areas, duracion


def espectro_neto(ruta_muestra, ruta_fondo=None, umbral_ruido=0.05,
                  num_canales=1024, rango_max=None, recortar_negativos=False,
                  metodo="gate", descartar_marcados=True,
                  usar_tiempo_vivo=True, canal=None, callback=None):
    """
    Devuelve un dict con canales, cuentas de muestra/fondo NORMALIZADAS y el neto,
    todo sobre EL MISMO eje de canales.

    ruta_fondo : puede ser None. Las corridas binarias de la sesion del
                1-sep-2026 se tomaron SIN fondo propio, y forzar uno prestado de
                otra sesion seria peor que no restar nada: el fondo depende del
                punto de operacion (a 1700 V el mismo cosmico deja el triple de
                area) y restar el de otro voltaje mete una forma ajena en el
                espectro. Sin fondo se devuelve neto = muestra, k_f = 0, y el
                informe lo dice para que no se lea como un neto.

    rango_max : tope del histograma [mV*ns]. Si None, usa el maximo global de ambos
                archivos -> garantiza bins identicos para muestra y fondo.
    recortar_negativos : True solo para graficar; False para ajustar (conserva el neto real).
    callback : callable(fraccion, texto) para la barra de progreso de la GUI.
               Los dos archivos se reparten la fraccion en proporcion a lo que
               PESAN, no mitad y mitad: el XML de la muestra es 5 veces el del
               fondo, y con el reparto ingenuo la barra se quedaba clavada en la
               primera mitad y el tiempo restante estimado salia el doble de lo
               que era. Si un archivo ya tiene cache su tramo se saltea de un
               golpe, pero el avance sigue siendo monotono.
    """
    hay_fondo = bool(ruta_fondo) and os.path.exists(ruta_fondo)
    try:
        tam = [float(os.path.getsize(ruta_muestra)),
               float(os.path.getsize(ruta_fondo)) if hay_fondo else 0.0]
    except OSError:
        tam = [1.0, 1.0]
    corte = 1.0 if not hay_fondo else tam[0] / max(tam[0] + tam[1], 1.0)
    tramos = ((0.0, corte), (corte, 1.0))

    def _sub(indice, etiqueta):
        if callback is None:
            return None
        a, b = tramos[indice]

        def cb(frac, n_ev):
            callback(a + frac * (b - a),
                     "%s: %s eventos" % (etiqueta, format(n_ev, ",d")))
        return cb

    # canal=None -> se detecta en CADA archivo por separado: muestra y fondo
    # pueden venir de sesiones distintas (y de hecho vienen).
    areas_m, t_m, ex_m = leer_areas_y_tiempo(ruta_muestra, umbral_ruido,
                                             metodo=metodo, devolver_extra=True,
                                             canal=canal,
                                             callback=_sub(0, "muestra"))
    if hay_fondo:
        areas_f, t_f, ex_f = leer_areas_y_tiempo(ruta_fondo, umbral_ruido,
                                                 metodo=metodo,
                                                 devolver_extra=True,
                                                 canal=canal,
                                                 callback=_sub(1, "fondo"))
    else:
        areas_f = np.empty(0)
        t_f, ex_f = 0.0, {"bandera": np.empty(0, np.int16),
                          "amp": np.empty(0), "t_on": np.empty(0),
                          "t_vivo": 0.0, "scaler": np.empty(0)}

    # Descarte por calidad del pulso (solo con metodo="gate"). Se aplica IGUAL a
    # muestra y fondo: si no, la resta deja de tener sentido. Es una eficiencia
    # comun, no una correccion.
    n_desc = (0, 0)
    if descartar_marcados and metodo != "legacy":
        bm = ex_m.get("bandera"); bf = ex_f.get("bandera")
        if bf is None:
            bf = np.empty(0, np.int16)
        if bm is not None and len(bm) == len(areas_m) and len(bf) == len(areas_f):
            # el bit 8 (linea base del post-pulso) NO descarta: el evento es
            # utilizable, solo tiene el area ~1.6% alta.
            mal_m, mal_f = (bm & ~8) != 0, (bf & ~8) != 0
            n_desc = (int(mal_m.sum()), int(mal_f.sum()))
            areas_m = areas_m[~mal_m]
            areas_f = areas_f[~mal_f]

    # ---- RANGE COMUN: la clave para que la resta tenga sentido ----
    if rango_max is None:
        todo = np.concatenate([areas_m, areas_f]) if len(areas_m)+len(areas_f) else np.array([1.0])
        rango_max = todo.max()
    rango = (0.0, rango_max)

    cuentas_m, bordes = np.histogram(areas_m, bins=num_canales, range=rango)
    cuentas_f, _      = np.histogram(areas_f, bins=num_canales, range=rango)  # MISMOS bordes
    canales = 0.5 * (bordes[:-1] + bordes[1:])

    # ---- normalizacion por TIEMPO VIVO ----
    # Con tiempo de pared la resta esta mal: la corrida de Co-60 perdio el 60%
    # de los triggers por el volcado a XML y la de fondo casi ninguno, asi que
    # el fondo entraba 2.5 veces mas fuerte de lo que corresponde (`tiempo_vivo`).
    tv_m = float(ex_m.get("t_vivo", t_m)) if usar_tiempo_vivo else t_m
    tv_f = float(ex_f.get("t_vivo", t_f)) if usar_tiempo_vivo else t_f
    tv_m = tv_m if tv_m > 0 else t_m
    tv_f = tv_f if tv_f > 0 else t_f

    # Sin fondo no hay nada que igualar: la muestra entra tal cual (k_m = 1) y
    # el fondo pesa cero. Poner t_norm = min(...) con tv_f = 0 dejaria la
    # muestra multiplicada por 0, o sea un espectro vacio sin decir por que.
    t_norm = min(tv_m, tv_f) if hay_fondo else tv_m
    m_norm = cuentas_m * (t_norm / tv_m) if tv_m > 0 else cuentas_m.astype(float)
    f_norm = (cuentas_f * (t_norm / tv_f) if (hay_fondo and tv_f > 0)
              else np.zeros_like(cuentas_f, dtype=float))

    neto = m_norm - f_norm
    if recortar_negativos:
        neto = np.maximum(neto, 0)

    # ---- barra de error del neto (Poisson propagado) ----
    # cada histograma crudo tiene var = N; al escalar por k la var va como k^2*N.
    # El piso de 1 evita peso infinito en los bins vacios.
    k_m = (t_norm / tv_m) if tv_m > 0 else 1.0
    k_f = (t_norm / tv_f) if (hay_fondo and tv_f > 0) else 0.0
    var_m = cuentas_m * k_m ** 2
    var_f = cuentas_f * k_f ** 2
    error = np.sqrt(np.maximum(var_m + var_f, 1.0))

    return {
        "canales": canales, "muestra": m_norm, "fondo": f_norm, "neto": neto,
        "error": error, "k_m": float(k_m), "k_f": float(k_f),
        "t_muestra": t_m, "t_fondo": t_f, "rango_max": rango_max,
        "metodo": metodo, "n_descartados": n_desc, "hay_fondo": hay_fondo,
        "t_vivo_muestra": tv_m, "t_vivo_fondo": tv_f,
        "extra_muestra": ex_m, "extra_fondo": ex_f,
    }


def varianza_predicha(fondo_norm, k_m, k_f):
    """
    Devuelve var(modelo) -> varianza ESPERADA del neto si el modelo es cierto.

    Por que hace falta [fuera del paper]: `error` de arriba se calcula con las
    cuentas OBSERVADAS. Eso es el chi2 de Neyman, y tiene un sesgo conocido: un
    bin que fluctua bajo recibe una barra mas chica, o sea MAS peso, asi que el
    ajuste se deja arrastrar hacia abajo. El sesgo va como 1/N por bin. Medido
    sobre sintetico con verdad conocida (CLAUDE.md 7.31):

        2.6 cuentas/bin -> D sesgado -1.9%
         26 cuentas/bin -> D sesgado -0.7%
        258 cuentas/bin -> D sesgado -0.05%

    y este experimento tiene ~15-25 cuentas/bin en la ventana de ajuste, o sea
    el sesgo es del orden del error estadistico que reportamos (0.6%).

    La varianza esperada se arma con el MODELO en la parte de muestra y con el
    fondo medido en la del fondo:
        var = k_m*(modelo + fondo_norm) + k_f*fondo_norm
    (usando var(k*N) = k^2*N = k*E[k*N], y E[muestra_norm] = neto + fondo_norm).
    """
    fondo_norm = np.asarray(fondo_norm, float)

    def var(modelo):
        return np.maximum(k_m * (np.asarray(modelo, float) + fondo_norm)
                          + k_f * fondo_norm, 1.0)
    return var


def suavizar(cuentas, ventana=15, orden=3):
    """Savitzky-Golay para el NCC (Etapa 4). Preserva la forma del borde mejor que
    un promedio movil. Requiere scipy; ventana impar y > orden."""
    from scipy.signal import savgol_filter
    ventana = ventana if ventana % 2 == 1 else ventana + 1
    return savgol_filter(cuentas, ventana, orden)
