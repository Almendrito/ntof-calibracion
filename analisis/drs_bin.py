"""
drs_bin.py — ETAPA 0 bis: leer el formato BINARIO del DRS4.

Por que existe
--------------
La sesion del 1-sep-2026 (los cuatro detectores, y el Co-60 del detector A) se
guardo con "Save as binary" del DRSOSC en vez de XML. Es el mismo instrumento y
el mismo pulso; cambia solo el envase. Este modulo expone EXACTAMENTE la misma
interfaz que `drs_xml` (`detectar_canal`, `metadatos`, `iter_eventos`,
`leer_pulsos`), asi que aguas abajo no hay que tocar nada: `drs.py` elige el
lector mirando los primeros bytes del archivo.

El formato, verificado byte a byte contra los archivos de la sesion
------------------------------------------------------------------
Cabecera:
    "DRS2"                                   4 bytes   (version del archivo)
    "TIME"                                   4 bytes
    por placa:
        "B#" + numero de serie (uint16)      4 bytes
        por canal:
            "C00n"                           4 bytes
            ancho de cada celda (1024 f32)   4096 bytes   [ns]
Evento (se repite hasta el final):
    "EHDR"                                   4 bytes
    numero de serie (int32)                  4 bytes
    ano, mes, dia, hora, min, seg, ms, rango 8 x uint16 = 16 bytes
    por placa:
        "B#" + serie (uint16)                4 bytes
        "T#" + celda de trigger (uint16)     4 bytes
        por canal:
            "C00n"                           4 bytes
            scaler (int32) [Hz]              4 bytes    <- puede no estar
            muestras (1024 x uint16)         2048 bytes

Dos cosas que NO se dan por supuestas y se deducen del archivo:

  - **Si el scaler esta o no.** Las versiones viejas del DRSOSC no lo escriben.
    Se prueban los dos tamanos de evento y se toma el que divide EXACTO el
    tamano del archivo; despues se verifica que "EHDR" caiga en el lugar
    predicho en varios eventos repartidos por todo el archivo. Si ninguno de los
    dos cierra, se levanta ValueError en vez de devolver basura.
  - **Cuantos canales y cuantas placas.** Se cuentan en la cabecera.

Conversion de unidades (la misma que el `read_binary.cpp` del DRS4):

    v[mV] = (u16/65536 + rango/1000 - 0.5) * 1000
    t[i]  = suma de los anchos de celda (j + celda_de_trigger) % 1024, j < i

O sea el eje temporal NO es uniforme y depende de la celda de trigger de cada
evento, igual que en el XML. La polaridad no se toca: `pulso.py` la deduce del
propio pulso.

`python analisis/drs_bin.py` imprime los metadatos de los archivos binarios que
encuentre y dibuja unos pulsos.
"""

import os
import struct
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "analisis"

from analisis import rutas

MAGIA = (b"DRS2", b"DRS8")      # cabeceras de archivo conocidas
N_MUESTRAS = 1024               # el DRS4 siempre vuelca 1024 celdas por canal


def es_binario(ruta):
    """True si el archivo empieza con la cabecera binaria del DRS4."""
    try:
        with open(ruta, "rb") as f:
            return f.read(4) in MAGIA
    except OSError:
        return False


# --------------------------------------------------------------------------
# 1. La cabecera: placas, canales, anchos de celda, tamano de evento
# --------------------------------------------------------------------------
def estructura(ruta):
    """
    Lee la cabecera y deduce el tamano de evento. Cuesta milisegundos.

    Devuelve un dict con: version, placas (lista de dicts con serie, canales y
    anchos), n_canales, offset_datos, con_scaler, tam_evento, n_eventos.
    """
    tam = os.path.getsize(ruta)
    with open(ruta, "rb") as f:
        cab = f.read(4 + 4 + 64 * (4 + 4 + 4 * (4 + 4 * N_MUESTRAS)))
    if cab[0:4] not in MAGIA:
        raise ValueError("%s: no empieza con %s; no es un binario del DRS4"
                         % (os.path.basename(ruta),
                            " ni ".join(m.decode() for m in MAGIA)))
    if cab[4:8] != b"TIME":
        raise ValueError("%s: falta la cabecera 'TIME'" % os.path.basename(ruta))

    p = 8
    placas = []
    while cab[p:p + 2] == b"B#":
        serie = struct.unpack_from("<H", cab, p + 2)[0]
        p += 4
        canales = []
        while cab[p:p + 3] == b"C00" and cab[p + 3:p + 4].isdigit():
            n = int(cab[p + 3:p + 4])
            ancho = np.frombuffer(cab, dtype="<f4", count=N_MUESTRAS,
                                  offset=p + 4).copy()
            canales.append(dict(n=n, etiqueta="CHN%d" % n, ancho=ancho))
            p += 4 + 4 * N_MUESTRAS
        if not canales:
            raise ValueError("%s: la placa B#%d no declara ningun canal"
                             % (os.path.basename(ruta), serie))
        placas.append(dict(serie=serie, canales=canales))
    if not placas:
        raise ValueError("%s: no encontre ninguna placa 'B#' en la cabecera"
                         % os.path.basename(ruta))

    offset = p
    n_ch = sum(len(b["canales"]) for b in placas)
    resto = tam - offset

    # El scaler puede estar o no segun la version del DRSOSC. Se elige el
    # tamano de evento que divide exacto Y que deja "EHDR" en su lugar.
    for con_scaler in (True, False):
        por_canal = 4 + (4 if con_scaler else 0) + 2 * N_MUESTRAS
        tam_ev = 4 + 4 + 16 + sum(4 + 4 + len(b["canales"]) * por_canal
                                  for b in placas)
        if resto <= 0 or resto % tam_ev:
            continue
        n_ev = resto // tam_ev
        if _verificar(ruta, offset, tam_ev, n_ev):
            return dict(version=cab[0:4].decode(), placas=placas,
                        n_canales=n_ch, offset_datos=offset,
                        con_scaler=con_scaler, tam_evento=tam_ev,
                        n_eventos=int(n_ev), bytes=tam, ruta=ruta)

    raise ValueError(
        "%s: la cabecera cierra pero el cuerpo no. Quedan %d bytes de datos, "
        "que no son multiplo de ningun tamano de evento esperado (%d canales, "
        "con o sin scaler). Archivo truncado o formato distinto."
        % (os.path.basename(ruta), resto, n_ch))


def _verificar(ruta, offset, tam_ev, n_ev, n_pruebas=8):
    """Comprueba que 'EHDR' caiga donde el tamano de evento predice."""
    if n_ev < 1:
        return False
    idx = np.unique(np.linspace(0, n_ev - 1, min(n_pruebas, n_ev)).astype(int))
    with open(ruta, "rb") as f:
        for k in idx:
            f.seek(offset + int(k) * tam_ev)
            if f.read(4) != b"EHDR":
                return False
    return True


# --------------------------------------------------------------------------
# 2. La misma interfaz que drs_xml
# --------------------------------------------------------------------------
def detectar_canal(ruta, max_lineas=None):
    """
    `(canal, tag_scaler)` como `drs_xml.detectar_canal`, para que la clave del
    cache y los informes hablen el mismo idioma ("CHN1", "Scaler0").

    `max_lineas` se acepta y se ignora: existe solo para que las dos firmas sean
    intercambiables.
    """
    e = estructura(ruta)
    ch = e["placas"][0]["canales"][0]
    return ch["etiqueta"], ("Scaler%d" % (ch["n"] - 1) if e["con_scaler"] else None)


def _localizar(estr, canal):
    """
    Indice de placa, indice de canal y offset del canal pedido dentro del
    evento. `canal` puede ser "CHN1", "C001" o None (el primero que haya).
    """
    if canal is None:
        ib, ic = 0, 0
    else:
        s = str(canal).upper()
        n = int("".join(c for c in s if c.isdigit()) or 0)
        ib = ic = None
        for i, b in enumerate(estr["placas"]):
            for j, c in enumerate(b["canales"]):
                if c["n"] == n:
                    ib, ic = i, j
                    break
            if ib is not None:
                break
        if ib is None:
            hay = ", ".join(c["etiqueta"] for b in estr["placas"]
                            for c in b["canales"])
            raise ValueError("%s: no tiene %s. Canales del archivo: %s"
                             % (os.path.basename(estr["ruta"]), canal, hay))

    por_canal = 4 + (4 if estr["con_scaler"] else 0) + 2 * N_MUESTRAS
    off = 4 + 4 + 16                                   # EHDR + serie + fecha
    for i, b in enumerate(estr["placas"]):
        off_placa = off
        off += 4 + 4                                   # "B#"+serie, "T#"+celda
        if i == ib:
            off_celda = off_placa + 4 + 2              # los 2 bytes de la celda
            return ib, ic, off + ic * por_canal, off_celda
        off += len(b["canales"]) * por_canal
    raise AssertionError("placa no encontrada")        # no deberia pasar


def _dtype(estr, off_canal, off_celda):
    """Vista estructurada de un evento: solo los campos que se usan."""
    sc = [("scaler", "<i4", (), off_canal + 4)] if estr["con_scaler"] else []
    campos = [("serie", "<i4", (), 4),
              ("fecha", "<u2", (8,), 8),
              ("celda", "<u2", (), off_celda)] + sc + \
             [("v", "<u2", (N_MUESTRAS,),
               off_canal + 4 + (4 if estr["con_scaler"] else 0))]
    return np.dtype(dict(names=[c[0] for c in campos],
                         formats=[(c[1], c[2]) if c[2] else c[1] for c in campos],
                         offsets=[c[3] for c in campos],
                         itemsize=estr["tam_evento"]))


def _reloj(fecha):
    """La marca de tiempo con el MISMO formato que el `<Time>` del XML."""
    return "%04d/%02d/%02d %02d:%02d:%02d.%03d" % (
        fecha[0], fecha[1], fecha[2], fecha[3], fecha[4], fecha[5], fecha[6])


def _eje_temporal(ancho, celda):
    """
    t[i] = suma de los anchos de las celdas (j + celda) % 1024 para j < i.

    Es la misma cuenta del `read_binary.cpp` del DRS4. El eje NO es uniforme y
    cambia evento a evento con la celda de trigger: por eso no se puede cachear
    uno solo.
    """
    rodado = np.roll(ancho, -int(celda))
    t = np.empty(N_MUESTRAS, float)
    t[0] = 0.0
    np.cumsum(rodado[:-1], out=t[1:])
    return t


def iter_eventos(ruta, canal=None, tag_scaler=None, limite=None,
                 callback=None, cada=400, bloque=2000):
    """
    Generador `(t, v, reloj, scaler)` por evento, igual que `drs_xml`.

    t : ndarray [ns], eje temporal real del DRS4 (no uniforme)
    v : ndarray [mV], senal cruda sin linea base restada
    reloj  : str "YYYY/MM/DD HH:MM:SS.mmm"
    scaler : float [Hz], o nan si el archivo no lo trae

    Lee por bloques con `np.memmap`, asi que la RAM queda plana aunque el
    archivo tenga cientos de MB.
    """
    estr = estructura(ruta)
    ib, ic, off_canal, off_celda = _localizar(estr, canal)
    ancho = estr["placas"][ib]["canales"][ic]["ancho"]
    dt = _dtype(estr, off_canal, off_celda)

    n_total = estr["n_eventos"] if limite is None else min(limite, estr["n_eventos"])
    mm = np.memmap(ruta, dtype=dt, mode="r", offset=estr["offset_datos"],
                   shape=(estr["n_eventos"],))
    try:
        hecho = 0
        while hecho < n_total:
            n = min(int(bloque), n_total - hecho)
            trozo = np.array(mm[hecho:hecho + n])          # copia el bloque
            # rango: uint16 en el archivo, pero representa mV con signo
            rango = trozo["fecha"][:, 7].astype(np.int32)
            rango = np.where(rango > 32767, rango - 65536, rango)
            volt = (trozo["v"].astype(np.float64) / 65536.0
                    + rango[:, None] / 1000.0 - 0.5) * 1000.0
            sc = (trozo["scaler"].astype(float) if estr["con_scaler"]
                  else np.full(n, np.nan))
            for k in range(n):
                yield (_eje_temporal(ancho, trozo["celda"][k]), volt[k],
                       _reloj(trozo["fecha"][k]), float(sc[k]))
            hecho += n
            if callback is not None and hecho % cada < n:
                callback(min(hecho / max(n_total, 1), 1.0), hecho)
    finally:
        del mm


def metadatos(ruta):
    """Cabecera legible, sin recorrer el archivo. Mismas claves que `drs_xml`."""
    estr = estructura(ruta)
    canal, tag = detectar_canal(ruta)
    t, v, reloj, scaler = next(iter_eventos(ruta, canal=canal, limite=1))
    d = np.diff(t)
    return dict(ruta=ruta, archivo=os.path.basename(ruta), bytes=estr["bytes"],
                formato="binario DRS4 (%s)" % estr["version"], creado=None,
                board=str(estr["placas"][0]["serie"]), canal=canal, scaler=tag,
                reloj_primero=reloj, h_unit="ns", v_unit="mV",
                muestras=int(t.size), t_inicial=float(t[0]), t_final=float(t[-1]),
                duracion_registro=float(t[-1] - t[0]),
                dt_medio=float(d.mean()), dt_min=float(d.min()),
                dt_max=float(d.max()), muestreo_GSs=1.0 / float(d.mean()),
                n_eventos=estr["n_eventos"], con_scaler=estr["con_scaler"],
                canales=[c["etiqueta"] for b in estr["placas"]
                         for c in b["canales"]],
                scaler_primero=scaler)


def imprimir_metadatos(info):
    """Vuelca el dict de `metadatos` en un bloque legible."""
    print("  archivo   : %s  (%.2f GB)" % (info["archivo"], info["bytes"] / 1e9))
    print("  formato   : %s" % info["formato"])
    print("  placa     : %s     canales: %s     scaler: %s"
          % (info["board"], ", ".join(info["canales"]),
             "si" if info["con_scaler"] else "NO"))
    print("  eventos   : %d" % info["n_eventos"])
    print("  primer ev.: %s   (scaler %.0f Hz)"
          % (info["reloj_primero"], info["scaler_primero"]))
    print("  registro  : %d muestras, %.1f ns (de %.2f a %.2f)"
          % (info["muestras"], info["duracion_registro"], info["t_inicial"],
             info["t_final"]))
    print("  muestreo  : %.2f GS/s   (dt medio %.4f ns, min %.4f, max %.4f)"
          % (info["muestreo_GSs"], info["dt_medio"], info["dt_min"],
             info["dt_max"]))


def leer_pulsos(ruta, n=200, canal=None, saltar=0):
    """
    `(T, V, relojes, scalers)` con `n` eventos apilados. Para inspeccionar.
    """
    T, V, relojes, scalers = [], [], [], []
    for k, (t, v, reloj, scaler) in enumerate(
            iter_eventos(ruta, canal=canal, limite=n + saltar)):
        if k < saltar:
            continue
        T.append(t); V.append(v); relojes.append(reloj); scalers.append(scaler)
        if len(T) >= n:
            break
    if not T:
        raise ValueError("%s: no pude leer ni un pulso" % os.path.basename(ruta))
    return np.array(T), np.array(V), relojes, np.array(scalers, float)


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.CORRIDAS
              if rutas.corrida(k).disponible() and es_binario(rutas.corrida(k).muestra)]
    print("=" * 74)
    print("drs_bin.py - archivos binarios del DRS4 en el catalogo")
    print("=" * 74)
    if not claves:
        print("  Ninguna corrida del catalogo apunta a un binario.")
        print("  Corre `python analisis/rutas.py` para ver que hay.")
        return

    fig, ax = plt.subplots(1, min(len(claves), 3), squeeze=False,
                           figsize=(5.0 * min(len(claves), 3), 3.8))
    for i, clave in enumerate(claves[:3]):
        c = rutas.corrida(clave)
        print()
        print("  %s  -  %s" % (clave, c.descripcion))
        imprimir_metadatos(metadatos(c.muestra))
        T, V, _, _ = leer_pulsos(c.muestra, n=40)
        a = ax[0][i]
        for k in range(0, 40, 4):
            a.plot(T[k], V[k], lw=0.7, alpha=0.8)
        a.set_title(clave, fontsize=10)
        a.set_xlabel("t [ns]")
        a.set_ylabel("v [mV]")
        a.grid(alpha=0.3)
    fig.tight_layout()
    s = rutas.figura("drs_bin_pulsos.png")
    fig.savefig(s, dpi=130)
    print()
    print("  figura: %s" % rutas.relativa(s))


if __name__ == "__main__":
    _demo()
