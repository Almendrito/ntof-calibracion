"""
drs_xml.py — ETAPA 0: abrir el XML del DRSOSC y reconstruir los pulsos.

Este modulo es SOLO lectura: entrega `(t, v)` en (ns, mV) tal como salieron del
digitalizador y no decide nada de fisica. La linea base, el onset y el area son
de `pulso.py`; separarlo asi permite mirar los pulsos crudos sin que ninguna
decision de analisis se haya colado antes.

Estructura del archivo
----------------------
    DRSOSC
      Event
        Serial, Time, HUnit, VUnit
        Board_2859
          Trigger_Cell
          Scaler0 / Scaler1 / ...     <- tasa de trigger REAL de la placa [Hz]
          CHN1 / CHN2 / ...
            Data   "t,v"              <- un tag por linea, t en ns y v en mV

Tres cosas que hay que saber antes de escribir una sola linea de codigo:

1. **`ET.parse` es inviable.** El XML del Co-60 pesa 1.76 GB y tiene ~5e7 tags
   `<Data>`; ElementTree construye un objeto por tag, o sea decenas de GB de RAM
   y varios minutos. Como el DRSOSC escribe un tag por linea, un escaneo lineal
   de TEXTO da exactamente lo mismo con RAM plana. Los pares "t,v" de un evento
   se juntan en un solo string y se convierten de una sola pasada en C con
   `np.fromstring(..., sep=",")`.

2. **El canal no es una constante del analisis.** Es como estaba enchufado el
   detector ese dia: las corridas de Co-60 de junio salieron por `CHN2` con
   `<Scaler1>`, y las del 26-ago-2026 por `CHN1` con `<Scaler0>`. Buscar el
   canal equivocado no da error: recorre el archivo entero y devuelve CERO
   pulsos, en silencio. Por eso `detectar_canal` lo saca del primer evento y
   `iter_eventos` lo usa por defecto.

3. **El eje temporal NO es uniforme.** El DRS4 muestrea con un anillo de celdas
   de capacitores y cada celda tiene su propio retardo, asi que el XML trae el
   tiempo de CADA muestra. Los `dt` van de ~0.13 a ~0.28 ns alrededor de una
   media de 0.195 ns (5 GS/s). Consecuencia: toda integral se hace con `x=t`,
   nunca multiplicando por un `dt` constante.

Correr `python analisis/drs_xml.py` imprime los metadatos de las corridas
disponibles y deja `figuras/00_pulsos_crudos.png`.
"""

import os
import re
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    # `python analisis/drs_xml.py` pone en sys.path la carpeta DEL ARCHIVO, no
    # la de arriba, asi que `from analisis...` no resolveria. Correr cada modulo
    # solo es una caracteristica del proyecto, no un accidente.
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas


_RE_CHN = re.compile(r"<(CHN\d+)>")
_RE_SCALER = re.compile(r"<(Scaler\d+)>")
_RE_BOARD = re.compile(r"<(Board_\w+)>")


# --------------------------------------------------------------------------
# 1. Que canal trae el pulso
# --------------------------------------------------------------------------
def detectar_canal(ruta, max_lineas=200000):
    """
    Devuelve `(canal, tag_scaler)` mirando el PRIMER evento del archivo.

    Un `<CHNn>` cuenta como canal del detector solo si adentro hay al menos un
    `<Data>`: el DRSOSC a veces abre canales vacios.

    El scaler del DRSOSC va indexado desde 0, asi que a `CHN{n}` le corresponde
    `<Scaler{n-1}>`. Si ese tag no aparece se usa el primero que haya, y si no
    hay ninguno se devuelve None (sin scaler no hay tiempo vivo, se asume 100%).
    """
    canales, scalers, actual = [], [], None
    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for k, linea in enumerate(f):
            s = linea.strip()
            if s.startswith("<Data>"):
                if actual and actual not in canales:
                    canales.append(actual)          # ese canal SI trae datos
                continue
            if s.startswith("</Event>") or k > max_lineas:
                break
            m = _RE_CHN.match(s)
            if m:
                actual = m.group(1)
                continue
            if s.startswith("</CHN"):
                actual = None
                continue
            m = _RE_SCALER.match(s)
            if m:
                scalers.append(m.group(1))
    if not canales:
        raise ValueError(
            "%s: no encontre ningun <CHNn> con <Data> en el primer evento. "
            "Es un XML del DRSOSC?" % os.path.basename(ruta))
    canal = canales[0]
    esperado = "Scaler%d" % (int(canal[3:]) - 1)
    tag = esperado if esperado in scalers else (scalers[0] if scalers else None)
    return canal, tag


# --------------------------------------------------------------------------
# 2. Metadatos baratos (solo el primer evento)
# --------------------------------------------------------------------------
def metadatos(ruta):
    """
    Cabecera del archivo sin recorrerlo entero: unidades, placa, canal, muestras
    por evento y el paso temporal medido. Cuesta milisegundos.
    """
    info = {"ruta": ruta, "archivo": os.path.basename(ruta),
            "bytes": os.path.getsize(ruta), "creado": None, "board": None,
            "h_unit": None, "v_unit": None, "reloj_primero": None}
    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for k, linea in enumerate(f):
            s = linea.strip()
            if s.startswith("<!--") and "created by" in s:
                info["creado"] = s.strip("<!->").strip()
            elif s.startswith("<Time>"):
                info["reloj_primero"] = s[6:-7]
            elif s.startswith("<HUnit>"):
                info["h_unit"] = s[7:-8]
            elif s.startswith("<VUnit>"):
                info["v_unit"] = s[7:-8]
            else:
                m = _RE_BOARD.match(s)
                if m:
                    info["board"] = m.group(1)
            if s.startswith("</Event>") or k > 200000:
                break

    info["canal"], info["scaler"] = detectar_canal(ruta)
    t, v, _, _ = next(iter_eventos(ruta, canal=info["canal"], limite=1))
    dt = np.diff(t)
    info.update(muestras=len(t), t_inicial=float(t[0]), t_final=float(t[-1]),
                duracion_registro=float(t[-1] - t[0]),
                dt_medio=float(dt.mean()), dt_min=float(dt.min()),
                dt_max=float(dt.max()),
                muestreo_GSs=1e-3 / float(dt.mean()) * 1e3 / 1e3)
    # 1/dt[ns] = GHz = GS/s. Se escribe explicito para que no haya duda:
    info["muestreo_GSs"] = 1.0 / info["dt_medio"]
    return info


def imprimir_metadatos(info):
    """Vuelca el dict de `metadatos` en un bloque legible."""
    print("  archivo   : %s  (%.2f GB)" % (info["archivo"], info["bytes"] / 1e9))
    if info["creado"]:
        print("  creado    : %s" % info["creado"])
    print("  placa     : %s     canal: %s     scaler: %s"
          % (info["board"], info["canal"], info["scaler"]))
    print("  primer ev.: %s" % info["reloj_primero"])
    print("  registro  : %d muestras, %.1f %s (de %.2f a %.2f)"
          % (info["muestras"], info["duracion_registro"], info["h_unit"] or "ns",
             info["t_inicial"], info["t_final"]))
    print("  paso dt   : medio %.4f ns  (min %.4f, max %.4f)  ->  %.2f GS/s"
          % (info["dt_medio"], info["dt_min"], info["dt_max"],
             info["muestreo_GSs"]))
    print("  ojo       : dt varia %.0f%% entre celdas -> integrar SIEMPRE con x=t"
          % (100 * (info["dt_max"] - info["dt_min"]) / info["dt_medio"]))


# --------------------------------------------------------------------------
# 3. El iterador: el unico lugar que toca el archivo grande
# --------------------------------------------------------------------------
def iter_eventos(ruta, canal=None, tag_scaler=None, limite=None,
                 callback=None, cada=400):
    """
    Generador que entrega `(t, v, reloj, scaler)` por evento.

    t : ndarray [ns], eje temporal REAL del DRS4 (no uniforme)
    v : ndarray [mV], la senal cruda, sin linea base restada
    reloj  : str con la marca `<Time>` del evento, o None
    scaler : float, tasa de trigger de la placa [Hz], o nan

    `limite=N` corta despues de N eventos. Es lo que hace posible iterar el
    codigo en segundos sobre un archivo de 2.9 GB.

    `callback(fraccion, n_eventos)` se llama cada `cada` eventos. La fraccion se
    mide en BYTES leidos y no en eventos, porque cuantos eventos tiene el
    archivo no se sabe hasta terminarlo mientras que su tamano se sabe desde el
    principio. Detalle que parece menor y no lo es: el archivo se abre en modo
    TEXTO, asi que Python traduce CRLF a LF y `len(linea)` queda un byte corto
    POR LINEA; con ~5e7 lineas eso es un 20% de error acumulado. Se mide el fin
    de linea real con `f.newlines` en el primer evento y se compensa.
    """
    if canal is None:
        canal, tag_auto = detectar_canal(ruta)
        if tag_scaler is None:
            tag_scaler = tag_auto
    if tag_scaler is None:
        tag_scaler = "Scaler%d" % (int(canal[3:]) - 1)

    ab, ce = "<%s>" % canal, "</%s>" % canal
    ab_sc, n_sc = "<%s>" % tag_scaler, len(tag_scaler) + 2

    pares, reloj, en_canal, scaler = [], None, False, np.nan
    total = max(os.path.getsize(ruta), 1)
    leidos = extra_nl = n_ev = 0

    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for linea in f:
            leidos += len(linea) + extra_nl
            s = linea.strip()
            if s.startswith("<Data>"):
                if en_canal:
                    pares.append(s[6:-7])           # saca <Data> y </Data>
            elif s.startswith(ab):
                en_canal = True
            elif s.startswith(ce):
                en_canal = False
            elif s.startswith("<Time>"):
                reloj = s[6:-7]
            elif s.startswith(ab_sc):
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
                t, v = _a_arrays(pares)
                yield t, v, reloj, scaler
                pares, reloj, scaler = [], None, np.nan
                if limite is not None and n_ev >= limite:
                    return


def _a_arrays(pares):
    """
    Convierte la lista de strings "t,v" del evento en dos ndarray.

    Se hace con UNA sola llamada a `np.fromstring` sobre todo el evento junto:
    parsear 1024 strings de a uno en Python cuesta ~30 veces mas, y eso, por 50
    mil eventos, es la diferencia entre 40 segundos y 20 minutos.
    """
    if not pares:
        return np.empty(0), np.empty(0)
    arr = np.fromstring(",".join(pares), sep=",")
    if arr.size < 4:
        return np.empty(0), np.empty(0)
    arr = arr[:(arr.size // 2) * 2].reshape(-1, 2)
    return arr[:, 0], arr[:, 1]


def leer_pulsos(ruta, n=200, canal=None, saltar=0):
    """
    Devuelve `(T, V, relojes, scalers)` con los primeros `n` eventos de igual
    largo, apilados en matrices `(n, muestras)`.

    Para inspeccionar, no para analizar: carga todo en memoria. Los eventos de
    largo distinto al del primero se descartan (los hay, muy pocos, cuando el
    volcado a XML se corta).
    """
    T, V, relojes, scalers = [], [], [], []
    largo = None
    for k, (t, v, reloj, scaler) in enumerate(
            iter_eventos(ruta, canal=canal, limite=None if saltar else n + saltar)):
        if k < saltar:
            continue
        if t.size == 0:
            continue
        if largo is None:
            largo = t.size
        if t.size != largo:
            continue
        T.append(t); V.append(v); relojes.append(reloj); scalers.append(scaler)
        if len(T) >= n:
            break
    if not T:
        raise ValueError("%s: no pude leer ni un pulso (canal correcto?)"
                         % os.path.basename(ruta))
    return (np.array(T), np.array(V), relojes,
            np.array(scalers, dtype=float))


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES if rutas.corrida(k).disponible()]
    if not claves:
        print("No hay ninguna corrida accesible. Corre `python analisis/rutas.py`.")
        return

    fig, axes = plt.subplots(len(claves), 3, figsize=(15, 4.2 * len(claves)),
                             squeeze=False)

    for fila, clave in enumerate(claves):
        c = rutas.corrida(clave)
        print("=" * 74)
        print("%s  —  %s" % (clave, c.descripcion))
        print("=" * 74)
        info = metadatos(c.muestra)
        imprimir_metadatos(info)

        T, V, _, sc = leer_pulsos(c.muestra, n=40, canal=info["canal"])
        t0 = T[0] - T[0][0]
        print("  scaler    : mediana %.0f Hz sobre los primeros %d eventos"
              % (np.nanmedian(sc), len(sc)))
        # polaridad: se DEDUCE, no se supone
        centro = np.median(V, axis=1, keepdims=True)
        desv = V - centro
        signo = -1 if abs(desv.min()) > abs(desv.max()) else +1
        print("  polaridad : %s (excursion max %.1f mV, min %.1f mV)"
              % ("NEGATIVA" if signo < 0 else "POSITIVA",
                 desv.max(), desv.min()))

        ax = axes[fila][0]
        for k in range(len(V)):
            ax.plot(t0, V[k], lw=0.6, alpha=0.5)
        ax.set_title("%s: 40 pulsos crudos (canal %s)" % (clave, info["canal"]))
        ax.set_xlabel("t desde el inicio del registro [ns]")
        ax.set_ylabel("v [mV]")
        ax.grid(alpha=0.3)

        ax = axes[fila][1]
        ax.plot(t0[:-1], np.diff(T[0]), lw=0.8)
        ax.axhline(info["dt_medio"], color="k", ls="--", lw=1,
                   label="medio %.4f ns" % info["dt_medio"])
        ax.set_title("paso temporal muestra a muestra (NO es uniforme)")
        ax.set_xlabel("t [ns]"); ax.set_ylabel("dt [ns]")
        ax.legend(); ax.grid(alpha=0.3)

        ax = axes[fila][2]
        ax.hist(np.diff(T, axis=1).ravel(), bins=80, color="C1")
        ax.axvline(info["dt_medio"], color="k", ls="--", lw=1)
        ax.set_title("distribucion de dt (%d eventos)" % len(T))
        ax.set_xlabel("dt [ns]"); ax.set_ylabel("cuentas")
        ax.grid(alpha=0.3)
        print()

    fig.tight_layout()
    salida = rutas.figura("00_pulsos_crudos.png")
    fig.savefig(salida, dpi=130)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
