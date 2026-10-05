"""
Etapa 0 — leer los pulsos crudos que entrega el DRS4.

    python desde_cero/lector.py

Que resuelve esta etapa
-----------------------
El DRS4 guarda cada disparo (evento) como un registro de 1024 muestras. Todo el
analisis posterior —area del pulso, espectro, calibracion— se hace sobre esas
muestras, asi que si esto esta mal, esta mal todo lo de aguas abajo. Aca no se
calcula nada de fisica: solo se convierte "archivo en disco" en

    (t, v)   con t en ns y v en mV,  un par de arrays por evento.

El DRSOSC escribe en DOS formatos y en este proyecto hay de los dos, asi que hay
un lector para cada uno y UNA sola puerta de entrada (`eventos`) que elige.

Las cuatro decisiones de esta etapa, y por que
----------------------------------------------
1. **El formato se decide por el CONTENIDO, no por el nombre.** Los archivos
   binarios de la sesion del 1-sep no tienen extension (se llaman `Co60`,
   `Cs137`), asi que mirar el nombre no serviria de nada. Se miran los primeros
   4 bytes: el binario empieza con "DRS2".

2. **Se lee en STREAMING, un evento a la vez.** El XML de Co-60 pesa 1.76 GB y
   tiene ~5e7 tags. `xml.etree.ElementTree` construiria un objeto por tag:
   decenas de GB de RAM. Por eso `eventos` es un GENERADOR (`yield`): el que lo
   usa procesa un evento, lo tira, y pide el siguiente. La RAM queda plana sin
   importar cuanto pese el archivo.

3. **El canal NO esta cableado.** En que `CHN{n}` sale el pulso depende de como
   estaba enchufado el detector ese dia, no del analisis. Si uno lo fija a mano
   y se equivoca, el lector recorre el archivo entero, no encuentra nada, y
   devuelve cero pulsos SIN ERROR: el problema aparece cinco pasos despues como
   un histograma vacio. Aca el canal se DEDUCE del archivo.

4. **El eje temporal no es uniforme y hay que reconstruirlo por evento.** Ver
   la explicacion larga en `_eventos_bin`.

Lo que esta etapa NO hace a proposito (viene despues)
-----------------------------------------------------
No integra el pulso, no resta linea base, no cachea, no calcula tiempo vivo, no
histogramea. Cada una de esas cosas es una etapa aparte, y mezclarlas aca haria
imposible probar esta sola.
"""

import os
import struct

import numpy as np


N_MUESTRAS = 1024          # el DRS4 siempre vuelca sus 1024 celdas
MAGIA_BIN = (b"DRS2", b"DRS8")


# ==========================================================================
# 1. De que formato es el archivo
# ==========================================================================
def formato(ruta):
    """'binario' o 'xml', mirando los primeros 4 bytes del archivo."""
    with open(ruta, "rb") as f:
        return "binario" if f.read(4) in MAGIA_BIN else "xml"


# ==========================================================================
# 2. El XML del DRSOSC
# ==========================================================================
# Como se ve por dentro (un tag por linea, y eso es lo que lo hace facil):
#
#     <DRSOSC>
#       <Event>
#         <Time>2026/06/24 11:56:18.456</Time>
#         <Board_2859>
#           <Trigger_Cell>793</Trigger_Cell>
#           <CHN2>
#             <Data>0.081,2.3</Data>      <- "t,v" : ns y mV
#             ... 1024 veces ...
#           </CHN2>
#         </Board_2859>
#       </Event>
#       ... un <Event> por disparo ...
#
# El eje temporal ya viene calculado adentro del archivo, asi que del XML se
# lee y listo. En el binario hay que reconstruirlo (punto 3).

def _canal_xml(ruta, max_lineas=200_000):
    """
    En que `CHN{n}` estan los datos, mirando SOLO el primer evento.

    Un `<CHNn>` cuenta como el canal del detector si adentro hay al menos un
    `<Data>`: el DRSOSC a veces abre canales vacios, y quedarse con el primero
    que aparece sin mirar si trae datos es como no deducir nada.
    """
    canal_abierto, con_datos = None, []
    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for k, linea in enumerate(f):
            s = linea.strip()
            if s.startswith("<Data>"):
                if canal_abierto and canal_abierto not in con_datos:
                    con_datos.append(canal_abierto)
            elif s.startswith("<CHN"):
                canal_abierto = s[1:s.index(">")]
            elif s.startswith("</CHN"):
                canal_abierto = None
            elif s.startswith("</Event>") or k > max_lineas:
                break
    if not con_datos:
        raise ValueError("%s: ningun <CHNn> con <Data> en el primer evento. "
                         "Es un XML del DRSOSC?" % os.path.basename(ruta))
    return con_datos[0]


def _scaler_xml(ruta, canal, max_lineas=200_000):
    """
    El tag del scaler que corresponde a `canal`, mirando el primer evento.

    El DRSOSC numera los scalers desde 0 y los canales desde 1: a `CHN{n}` le
    toca `<Scaler{n-1}>` (junio: CHN2 + Scaler1; agosto: CHN1 + Scaler0). Si
    ese no esta se usa el primero que aparezca, y si no hay ninguno, None: sin
    scaler no hay tiempo vivo y la Etapa 3 lo tiene que decir.
    """
    vistos = []
    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for k, linea in enumerate(f):
            s = linea.strip()
            if s.startswith("<Scaler"):
                vistos.append(s[1:s.index(">")])
            elif s.startswith("</Event>") or k > max_lineas:
                break
    esperado = "Scaler%d" % (int(canal[3:]) - 1)
    return esperado if esperado in vistos else (vistos[0] if vistos else None)


def _eventos_xml(ruta, canal=None, limite=None, con_reloj=False):
    """Generador de `(t, v)` leyendo el XML como TEXTO, linea por linea."""
    canal = canal or _canal_xml(ruta)
    abre, cierra = "<%s>" % canal, "</%s>" % canal
    tag_sc = _scaler_xml(ruta, canal) if con_reloj else None
    abre_sc = "<%s>" % tag_sc if tag_sc else None

    pares, dentro, n = [], False, 0
    reloj, scaler = None, np.nan
    with open(ruta, "r", encoding="ISO-8859-1", errors="replace") as f:
        for linea in f:
            s = linea.strip()
            if s.startswith("<Data>"):
                if dentro:
                    pares.append(s[6:-7])          # saca <Data> y </Data>
            elif s.startswith(abre):
                dentro = True
            elif s.startswith(cierra):
                dentro = False
            elif con_reloj and s.startswith("<Time>"):
                reloj = s[6:-7]
            elif abre_sc and s.startswith(abre_sc):
                try:
                    scaler = float(s[len(abre_sc):-(len(abre_sc) + 1)])
                except ValueError:
                    scaler = np.nan
            elif s.startswith("</Event>"):
                if pares:
                    # UNA sola conversion en C para todo el pulso: pegar los
                    # 1024 pares con comas y dejar que numpy parsee de una.
                    # Convertir par por par en Python es ~20x mas lento.
                    a = np.fromstring(",".join(pares), sep=",")
                    a = a[:(a.size // 2) * 2].reshape(-1, 2)
                    if con_reloj:
                        yield a[:, 0], a[:, 1], reloj, scaler
                    else:
                        yield a[:, 0], a[:, 1]
                    n += 1
                    if limite is not None and n >= limite:
                        return
                pares, reloj, scaler = [], None, np.nan


# ==========================================================================
# 3. El binario del DRS4
# ==========================================================================
# El formato, byte a byte:
#
#   CABECERA (una vez)
#     "DRS2"                                      4 B
#     "TIME"                                      4 B
#     "B#" + numero de serie de la placa (u16)    4 B
#     por cada canal:
#        "C00n"                                   4 B
#        ancho de las 1024 celdas (float32)    4096 B   [ns]
#
#   EVENTO (se repite hasta el final del archivo)
#     "EHDR"                                      4 B
#     numero de serie del evento (int32)          4 B
#     ano, mes, dia, hora, min, seg, ms, rango    8 x u16 = 16 B
#     "B#" + serie                                4 B
#     "T#" + celda de trigger (u16)               4 B
#     por cada canal:
#        "C00n"                                   4 B
#        scaler (int32) [Hz]                      4 B   <- puede NO estar
#        1024 muestras (uint16)                2048 B

def _cabecera_bin(ruta):
    """
    Lee la cabecera y deduce el tamano de un evento. Cuesta milisegundos.

    Dos cosas se DEDUCEN del archivo en vez de darlas por sabidas:

    * **cuantos canales hay** — se cuentan los bloques "C00n" de la cabecera;
    * **si el scaler esta o no** — las versiones viejas del DRSOSC no lo
      escriben, y eso cambia el tamano del evento en 4 bytes por canal. Se
      prueban los dos tamanos y se elige el que (a) divide exacto lo que queda
      del archivo y (b) deja "EHDR" justo donde predice. Si ninguno cierra se
      levanta un error, en vez de devolver basura con cara de dato.
    """
    tam = os.path.getsize(ruta)
    with open(ruta, "rb") as f:
        cab = f.read(8 + 4 + 8 * (4 + 4 * N_MUESTRAS))    # de sobra
    if cab[0:4] not in MAGIA_BIN or cab[4:8] != b"TIME":
        raise ValueError("%s: no es un binario del DRS4" % os.path.basename(ruta))

    p = 8
    serie = struct.unpack_from("<H", cab, p + 2)[0]
    p += 4                                                # "B#" + serie
    canales = []
    while cab[p:p + 3] == b"C00":
        canales.append(dict(
            etiqueta="CHN%d" % int(cab[p + 3:p + 4]),
            ancho=np.frombuffer(cab, "<f4", N_MUESTRAS, p + 4).astype(float)))
        p += 4 + 4 * N_MUESTRAS
    if not canales:
        raise ValueError("%s: la cabecera no declara ningun canal"
                         % os.path.basename(ruta))

    resto = tam - p
    for con_scaler in (True, False):
        por_canal = 4 + (4 if con_scaler else 0) + 2 * N_MUESTRAS
        tam_ev = 4 + 4 + 16 + 4 + 4 + len(canales) * por_canal
        if resto > 0 and resto % tam_ev == 0 and _cae_el_EHDR(ruta, p, tam_ev,
                                                              resto // tam_ev):
            return dict(serie=serie, canales=canales, inicio=p,
                        con_scaler=con_scaler, tam_evento=tam_ev,
                        n_eventos=resto // tam_ev)
    raise ValueError("%s: la cabecera cierra pero el cuerpo no (%d bytes). "
                     "Archivo truncado o formato distinto."
                     % (os.path.basename(ruta), resto))


def _cae_el_EHDR(ruta, inicio, tam_ev, n_ev, pruebas=8):
    """Comprueba que 'EHDR' aparezca donde el tamano de evento predice."""
    if n_ev < 1:
        return False
    with open(ruta, "rb") as f:
        for k in np.unique(np.linspace(0, n_ev - 1, min(pruebas, n_ev)).astype(int)):
            f.seek(inicio + int(k) * tam_ev)
            if f.read(4) != b"EHDR":
                return False
    return True


def _eje_temporal(ancho, celda):
    """
    El eje de tiempos de UN evento, en ns.

    Por que no es simplemente `np.arange(1024) * 0.2`:

    El DRS4 no es un digitalizador comun: es un anillo de 1024 condensadores que
    se van cargando en cadena, y se congela cuando llega el trigger. De ahi salen
    dos cosas que no se pueden ignorar:

      * cada celda tiene su propio ancho temporal (variacion de fabricacion, ~5%),
        y por eso la cabecera trae los 1024 anchos medidos;
      * el anillo se congela en un punto distinto cada vez. La celda donde quedo
        el trigger (`T#`) dice DONDE empieza a leerse el anillo, asi que el orden
        en que hay que sumar los anchos cambia evento a evento.

    Entonces: t[i] = suma de los anchos de las celdas (celda + j) % 1024, j < i.
    Usar un eje uniforme mete un error de tiempo de hasta ~5% por muestra, que es
    justo del tamano del ancho de subida del pulso.
    """
    rodado = np.roll(ancho, -int(celda))
    t = np.empty(N_MUESTRAS)
    t[0] = 0.0
    np.cumsum(rodado[:-1], out=t[1:])
    return t


def _eventos_bin(ruta, canal=None, limite=None, con_reloj=False):
    """Generador de `(t, v)` sobre el binario."""
    cab = _cabecera_bin(ruta)
    i = 0 if canal is None else [c["etiqueta"] for c in cab["canales"]].index(canal)
    ancho = cab["canales"][i]["ancho"]

    # posicion del bloque de muestras del canal pedido, dentro del evento
    por_canal = 4 + (4 if cab["con_scaler"] else 0) + 2 * N_MUESTRAS
    off_celda = 4 + 4 + 16 + 4 + 2                     # los 2 bytes tras "T#"
    off_datos = (4 + 4 + 16 + 4 + 4 + i * por_canal
                 + 4 + (4 if cab["con_scaler"] else 0))

    n = cab["n_eventos"] if limite is None else min(limite, cab["n_eventos"])
    with open(ruta, "rb") as f:
        for k in range(n):
            base = cab["inicio"] + k * cab["tam_evento"]
            # El "rango" es el OCTAVO de los 8 uint16 de la fecha, o sea
            # 4 (EHDR) + 4 (serie) + 7*2 = byte 22. Contarlos mal es un error
            # mudo: en el byte 14 esta la HORA, que tambien es un numero
            # chico y plausible, y el pulso sale desplazado ~16 mV sin que
            # nada falle. Asi se encontro: comparando contra `analisis/`.
            f.seek(base + 4 + 4)
            fecha = struct.unpack("<7H", f.read(14))  # ano ... milisegundo
            rango = struct.unpack("<h", f.read(2))[0]  # con signo: son mV
            f.seek(base + off_celda)
            celda = struct.unpack("<H", f.read(2))[0]
            f.seek(base + off_datos)
            crudo = np.frombuffer(f.read(2 * N_MUESTRAS), "<u2")
            # La conversion del propio DRS4 (su read_binary.cpp): el u16 cubre
            # 1 V centrado en `rango`, asi que 0 -> rango-500 mV y 65535 -> +500.
            v = (crudo / 65536.0 + rango / 1000.0 - 0.5) * 1000.0
            if not con_reloj:
                yield _eje_temporal(ancho, celda), v
                continue
            scaler = np.nan
            if cab["con_scaler"]:
                # el scaler son los 4 bytes justo antes de las muestras
                f.seek(base + off_datos - 4)
                scaler = float(struct.unpack("<i", f.read(4))[0])
            # mismo formato que el <Time> del XML, para que lo de aguas abajo
            # no se entere de que archivo vino
            reloj = "%04d/%02d/%02d %02d:%02d:%02d.%03d" % fecha
            yield _eje_temporal(ancho, celda), v, reloj, scaler


# ==========================================================================
# 4. La puerta unica
# ==========================================================================
def eventos(ruta, canal=None, limite=None, con_reloj=False):
    """
    Generador de `(t, v)` por evento — t en ns, v en mV — sea cual sea el formato.

    Con `con_reloj=True` entrega `(t, v, reloj, scaler)`:
      reloj  : str "AAAA/MM/DD hh:mm:ss.mmm", la marca del evento
      scaler : float, tasa de trigger REAL de la placa [Hz], o nan si el
               archivo no la trae
    Son lo que hace falta para el tiempo vivo (Etapa 3). Van aparte, y no por
    defecto, para que todo lo que ya pedia `(t, v)` siga igual.

    Todo lo que venga despues (area, espectro, ajuste) pide los datos por aca y
    no se entera de si el archivo era XML o binario. Ese es el punto: el envase
    cambia, el pulso no.
    """
    lector = _eventos_bin if formato(ruta) == "binario" else _eventos_xml
    return lector(ruta, canal=canal, limite=limite, con_reloj=con_reloj)


def info(ruta):
    """Lo que se sabe del archivo sin recorrerlo entero. Para mirarlo antes de usarlo."""
    d = dict(archivo=os.path.basename(ruta), bytes=os.path.getsize(ruta),
             formato=formato(ruta))
    if d["formato"] == "binario":
        cab = _cabecera_bin(ruta)
        d.update(placa=cab["serie"], n_eventos=cab["n_eventos"],
                 canales=[c["etiqueta"] for c in cab["canales"]],
                 scaler=cab["con_scaler"])
    else:
        d.update(placa=None, n_eventos=None, canales=[_canal_xml(ruta)],
                 scaler=None)
    t, v = next(eventos(ruta, limite=1))
    dt = np.diff(t)
    d.update(muestras=t.size, duracion=float(t[-1] - t[0]),
             dt_medio=float(dt.mean()), dt_min=float(dt.min()),
             dt_max=float(dt.max()), v_min=float(v.min()), v_max=float(v.max()))
    return d


# ==========================================================================
# Demo: mirar lo que se acaba de leer
# ==========================================================================
def _revisar(d):
    """
    Comprobaciones baratas sobre el primer evento.

    No prueban que el analisis este bien; prueban que lo LEIDO tenga sentido, que
    es lo primero que hay que descartar cuando un numero sale raro.
    """
    avisos = []
    if d["muestras"] != N_MUESTRAS:
        avisos.append("el registro no tiene %d muestras sino %d"
                      % (N_MUESTRAS, d["muestras"]))
    if d["dt_min"] <= -0.10 * d["dt_medio"]:
        # Un dt levemente negativo NO es un error de lectura: la calibracion
        # temporal del DRS4 puede dejar alguna celda con ancho ligeramente
        # negativo. En el Co-60 a 1600 V la celda 498 mide -0.0022 ns: 1 celda
        # de 1024, el 1% de un dt y el 0.001% del registro. Eso no es un error
        # de lectura. Se avisa solo si el retroceso es de al menos un 10% de un
        # dt, que ya seria un eje mal reconstruido.
        avisos.append("el eje temporal retrocede de verdad (dt_min = %.3f ns "
                      "contra dt medio %.3f)" % (d["dt_min"], d["dt_medio"]))
    if not 0.05 < d["dt_medio"] < 1.0:
        avisos.append("dt medio = %.3f ns: fuera de lo esperable para el DRS4"
                      % d["dt_medio"])
    if max(abs(d["v_min"]), abs(d["v_max"])) > 600:
        avisos.append("tensiones fuera del rango del DRS4 (+-500 mV)")
    return avisos


RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DET_A = os.path.join(RAIZ, "datos", "Detector A")

# Las tres corridas de trabajo de PLAN.md, con la MISMA clave que en
# analisis/rutas.py. Las rutas se escriben aca y no se importan de alla porque
# desde_cero/ no importa analisis/ salvo en los verificar_*.py.
# (muestra, fondo); fondo None = la corrida no tiene fondo propio.
CORRIDAS_TRABAJO = {
    "cs137_1600": (os.path.join(_DET_A, "1600V", "Cs137.xml"),
                   os.path.join(_DET_A, "1600V", "BG.xml")),
    "co60": (os.path.join(_DET_A, "A Co60 5GaS 50k 02mV trigg 1550 att.xml"),
             os.path.join(_DET_A, "A BG 5GaS 10k 02mV trigg 1550 att.xml")),
    "a_co60_1600": (os.path.join(_DET_A, "1600V", "Co60"), None),
}


def _archivos_de_ejemplo():
    """Un XML y un binario, si estan a mano."""
    # el XML de la MUESTRA, no el de fondo: en el fondo casi no hay pulsos y la
    # figura sale de puro ruido. Se leen 200 eventos, o sea unos MB, aunque el
    # archivo pese 1.76 GB: para eso el lector es un generador.
    # (Hasta el 27-sep apuntaba a datos/<xml> y a D:\..., que ya no existen:
    # la campana se copio a datos/Detector X/ en §7.61.)
    candidatos = [CORRIDAS_TRABAJO["co60"][0],
                  CORRIDAS_TRABAJO["a_co60_1600"][0]]
    return [c for c in candidatos if os.path.exists(c)]


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas = _archivos_de_ejemplo()
    if not rutas:
        print("No encontre archivos de ejemplo. Pasale una ruta a `eventos()`.")
        return

    fig, ejes = plt.subplots(1, len(rutas), figsize=(6.0 * len(rutas), 4.0),
                             squeeze=False)
    for eje, ruta in zip(ejes[0], rutas):
        d = info(ruta)
        print()
        print("=" * 74)
        print("  %s   [%s]" % (d["archivo"], d["formato"]))
        print("=" * 74)
        print("  tamano    : %.2f GB" % (d["bytes"] / 1e9))
        print("  canales   : %s   (se usa %s)"
              % (", ".join(d["canales"]), d["canales"][0]))
        if d["n_eventos"] is not None:
            print("  eventos   : %d   placa %s   scaler: %s"
                  % (d["n_eventos"], d["placa"], "si" if d["scaler"] else "no"))
        print("  registro  : %d muestras, %.1f ns"
              % (d["muestras"], d["duracion"]))
        print("  muestreo  : %.2f GS/s   dt medio %.4f ns (min %.4f, max %.4f)"
              % (1.0 / d["dt_medio"], d["dt_medio"], d["dt_min"], d["dt_max"]))
        print("              ^ dt_min != dt_max: el DRS4 NO muestrea uniforme")
        avisos = _revisar(d)
        print("  revision  : %s" % ("todo coherente" if not avisos
                                    else "; ".join(avisos)))

        for k, (t, v) in enumerate(eventos(ruta, limite=200)):
            if k % 25 == 0:
                eje.plot(t, v, lw=0.8, alpha=0.8)
        eje.set_title("%s  [%s]" % (d["archivo"][:34], d["formato"]), fontsize=10)
        eje.set_xlabel("t [ns]")
        eje.set_ylabel("v [mV]")
        eje.grid(alpha=0.3)

    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "figuras")
    os.makedirs(salida, exist_ok=True)
    salida = os.path.join(salida, "etapa0_pulsos.png")
    fig.savefig(salida, dpi=130)
    print()
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
