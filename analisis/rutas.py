"""
rutas.py — el unico modulo que sabe DONDE esta cada cosa.

Por que existe
--------------
Si cada modulo nombra sus archivos con rutas relativas al directorio de trabajo,
el proyecto solo funciona si uno esta parado en la carpeta correcta. Aca todo se
resuelve desde `__file__`, asi que `python analisis/pulso.py` da lo mismo desde
donde se lo llame.

Ademas centraliza el CATALOGO DE CORRIDAS: que archivo es la muestra, cual el
fondo, y de que fuente es cada una. Un dato del catalogo puede no estar montado
(las corridas de Cs-137 viven en un disco externo), asi que nada de esto explota
al importarse: se pregunta con `disponible()`.

Correr `python analisis/rutas.py` lista lo que ve y lo que falta.
"""

import os
import sys
import time

# ---------------------------------------------------------------- carpetas ---
AQUI = os.path.dirname(os.path.abspath(__file__))           # .../analisis
RAIZ = os.path.dirname(AQUI)                                # .../Calibracion CE
DATOS = os.path.join(RAIZ, "datos")
FIGURAS = os.path.join(AQUI, "figuras")
CACHE = os.path.join(AQUI, ".cache")

# El disco externo. Puede no estar conectado.
#
# Ojo con el nombre de la carpeta: dice 26_08_2026 pero adentro conviven DOS
# sesiones. Las marcas de tiempo que el propio DRS4 escribe en cada evento lo
# dicen sin ambiguedad:
#   - 26-ago-2026: Cs-137 del detector A a 1600/1650/1700 V, en XML, mas el
#     unico fondo que hay (BG.xml, solo a 1600 V).
#   - 01-sep-2026: los cuatro detectores (A, B, C, D) a 1600/1650/1700 V, en
#     formato BINARIO y sin extension en el nombre. Cs-137 en los cuatro y
#     Co-60 solo en el A. Sin fondo.
# La bitacora esta en `Bitacora.txt` dentro de esa misma carpeta.
NOMBRE_SESION = "Sesión Calibración 26_08_2026"
LETRA_HABITUAL = "D:\\"


def _buscar_sesion(nombre=NOMBRE_SESION):
    """
    Donde esta montada la carpeta de la sesion, ahora.

    La letra de unidad de un disco externo NO es una propiedad del disco: se la
    asigna Windows segun que haya conectado antes, asi que el mismo disco puede
    ser D: hoy y E: manana. Cablear "D:" hacia que el catalogo entero apareciera
    [NO MONTADO] sin decir por que, que es la misma clase de falla silenciosa
    que el canal cableado a CHN2 (CLAUDE.md 7.36): no falla nada, simplemente no
    hay datos.

    Orden de busqueda:
      1. la variable de entorno CALIBRACION_SESION, si apunta a algo que existe
         (escotilla para poner la carpeta donde uno quiera, p.ej. una copia
         local);
      2. la letra de siempre, D:, que es lo que ahorra el barrido casi siempre;
      3. el resto de las unidades, de la A a la Z.

    Si no esta montada devuelve la ruta con la letra habitual, para que los
    mensajes sigan diciendo donde deberia estar en vez de quedar vacios.
    """
    de_entorno = os.environ.get("CALIBRACION_SESION", "").strip()
    if de_entorno and os.path.isdir(de_entorno):
        return de_entorno
    for raiz in _unidades():
        ruta = os.path.join(raiz, nombre)
        try:
            if os.path.isdir(ruta):
                return ruta
        except OSError:
            continue
    return os.path.join(LETRA_HABITUAL, nombre)


def _unidades():
    """
    Las unidades que EXISTEN, empezando por la habitual.

    Ojo con la tentacion de probar de la A a la Z a ciegas: en Windows tocar una
    unidad extraible que existe pero esta vacia (una lectora de tarjetas sin
    tarjeta, un DVD vacio) tarda segundos y puede abrir el dialogo modal "No hay
    disco en la unidad", que desde un programa sin consola parece un cuelgue.
    `GetLogicalDrives` devuelve el mapa de bits de las que hay montadas y cuesta
    microsegundos, asi que solo se sondean esas.
    """
    letras = []
    if sys.platform == "win32":
        try:
            import ctypes
            mapa = ctypes.windll.kernel32.GetLogicalDrives()
            letras = ["%s:\\" % chr(ord("A") + i)
                      for i in range(26) if mapa >> i & 1]
        except Exception:                                      # noqa: BLE001
            letras = [LETRA_HABITUAL]
    if LETRA_HABITUAL in letras:                # la de siempre, primero
        letras.remove(LETRA_HABITUAL)
        letras.insert(0, LETRA_HABITUAL)
    return letras


SESION_EXTERNA = _buscar_sesion()
EXTERNO = os.path.join(SESION_EXTERNA, "Detector A")     # compatibilidad


# Los nombres que en el disco no dicen de que fuente son. En la copia local se
# renombraron para que se entiendan sin abrir el catalogo.
_RENOMBRE_LOCAL = {"B1700V": "Cs137"}


def _local(detector, carpeta, nombre):
    """
    La misma medicion, en la copia que vive dentro de `datos/`.

    La copia esta ordenada por detector y voltaje igual que el disco, con dos
    arreglos: la carpeta del detector A a 1700 V se llama "1700" en el disco y
    "1700V" aca (todas las demas ya tenian la V), y `B1700V` -el unico archivo
    cuyo nombre no dice que es- pasa a `Cs137`.
    """
    if not carpeta.upper().endswith("V"):
        carpeta += "V"
    raiz, ext = os.path.splitext(nombre)
    return os.path.join(DATOS, "Detector %s" % detector, carpeta,
                        _RENOMBRE_LOCAL.get(raiz, raiz) + ext)


def _ext(detector, carpeta, nombre):
    """
    Donde esta esa medicion: la copia local si esta, si no el disco externo.

    La copia local gana a proposito. El disco externo no siempre esta
    conectado, y cuando lo esta puede tener otra letra (`_buscar_sesion`);
    leyendo de `datos/` el proyecto funciona igual sin el disco. Las dos rutas
    dan el mismo archivo: la copia conserva la fecha de modificacion, y desde
    que la clave del cache de areas dejo de llevar la ruta (`area._clave`) las
    dos comparten tambien las areas ya calculadas.
    """
    local = _local(detector, carpeta, nombre)
    if os.path.exists(local):
        return local
    return os.path.join(SESION_EXTERNA, "Detector %s" % detector, carpeta, nombre)


def _junio(nombre):
    """
    Los dos XML de junio, que se movieron a `datos/Detector A/`.

    Estaban sueltos en la raiz de `datos/` cuando se armo el catalogo y despues
    pasaron a la carpeta de su detector, donde corresponde. Se prueban las dos
    ubicaciones en vez de cambiar la ruta y listo, porque el .npz de areas vive
    al lado del XML y el que ya existe -en la ubicacion vieja o en la nueva-
    sigue valiendo: la clave del cache es tamano + fecha, no la ruta.
    """
    nuevo = os.path.join(DATOS, "Detector A", nombre)
    return nuevo if os.path.exists(nuevo) else os.path.join(DATOS, nombre)


def _crear(carpeta):
    """Crea la carpeta si falta. Se llama tarde (al guardar), no al importar."""
    os.makedirs(carpeta, exist_ok=True)
    return carpeta


def figura(nombre):
    """Ruta absoluta de una figura de salida, creando `figuras/` si hace falta."""
    return os.path.join(_crear(FIGURAS), nombre)


def cache(nombre):
    """Ruta absoluta de un archivo de cache, creando `.cache/` si hace falta."""
    return os.path.join(_crear(CACHE), nombre)


def relativa(ruta):
    """Ruta corta para imprimir: relativa a la raiz si esta adentro, si no absoluta."""
    ruta = os.path.abspath(ruta)
    try:
        rel = os.path.relpath(ruta, RAIZ)
    except ValueError:                    # otra unidad (D:) -> no existe relativa
        return ruta
    return ruta if rel.startswith("..") else rel.replace("\\", "/")


# --------------------------------------------------------------- corridas ----
class Corrida:
    """
    Una medicion: el XML de la muestra, el del fondo, y de que fuente es.

    `fondo` puede ser None (corridas sin fondo propio); en ese caso el espectro
    se arma sin restar nada y hay que decirlo en el informe, no dejarlo pasar.

    El PUNTO DE OPERACION va aparte del detector
    --------------------------------------------
    `detector` es el detector fisico ("A"); `hv`, `atenuador`, `umbral` y
    `muestreo` son como estaba configurado ese dia. La separacion no es
    cosmetica: D no es una propiedad del detector sino del detector EN UN PUNTO
    DE OPERACION, y mezclar corridas de puntos distintos en una misma recta de
    calibracion da un numero sin sentido.

    Que pasa si no se separan: hasta el 28-ago-2026 el Co-60 se llamaba
    "A / SiPM MicroFC-30020" y el Cs-137 "A @ 1600 V", asi que el registro los
    mandaba a grupos distintos y el test de linealidad NUNCA se disparaba. Y al
    reves, si los dos se hubieran llamado "A" se habrian combinado en una recta
    aunque estan a 1550 y 1600 V, que es peor todavia.

    `atenuador` importa por lo mismo: el area se mide a la entrada del
    digitalizador, asi que un atenuador de 20 dB divide D por 10. Ningun modulo
    lo corrige (y esta bien: el observable es lo que entra al DRS4), pero dos
    corridas con atenuacion distinta NO son comparables.

    CONFIRMADO (22-sep-2026, por el operador): las 25 corridas del catalogo se
    tomaron con 20 dB en linea, puestos para que el pulso del PMT no pase el
    rango de entrada del DRS4. Hasta esa fecha el campo estaba en None en 24 de
    las 25 y `misma_cadena` avisaba por "sin registrar" en toda combinacion. El
    valor se escribe aca y NO se corrige en ningun lado: el observable sigue
    siendo lo que entra al digitalizador, asi que D esta medida ahi y la carga a
    la salida del PMT es 10 veces mayor.
    """

    def __init__(self, clave, descripcion, muestra, fondo=None, fuente=None,
                 detector=None, nota="", sim=None, hv=None, atenuador=None,
                 umbral=None, muestreo=None, distancia=None):
        self.clave = clave
        self.descripcion = descripcion
        self.muestra = muestra
        self.fondo = fondo
        self.fuente = fuente
        self.detector = detector
        self.nota = nota
        self.sim = sim          # .root de Geant4 de ESTA fuente, o None
        self.hv = hv                    # alto voltaje [V]
        self.atenuador = atenuador      # atenuacion en linea [dB], 0 = sin
        self.umbral = umbral            # nivel del discriminador [mV]
        self.muestreo = muestreo        # velocidad de muestreo [GS/s]
        self.distancia = distancia      # fuente-detector [mm], si se anoto

    def punto(self):
        """
        Etiqueta del punto de operacion. Es la clave por la que se agrupa.

        Compartir punto es NECESARIO para poner dos D en la misma recta, pero
        no alcanza: el resto de la cadena (atenuador, muestreo) tambien tiene
        que coincidir, y eso lo revisa `misma_cadena`. Se deja afuera de la
        etiqueta a proposito: un dato NO REGISTRADO no debe partir el grupo en
        silencio, tiene que salir como aviso.
        """
        s = self.detector or "sin nombre"
        if self.hv is not None:
            s += " @ %g V" % self.hv
        return s

    def cadena(self):
        """La cadena de senal en una linea, para el informe y la bitacora."""
        campos = [("HV", self.hv, "V"), ("atenuador", self.atenuador, "dB"),
                  ("umbral", self.umbral, "mV"), ("muestreo", self.muestreo, "GS/s"),
                  ("distancia", self.distancia, "mm")]
        vistos = ["%s %g %s" % (n, v, u) for n, v, u in campos if v is not None]
        return ", ".join(vistos) if vistos else "sin registrar"

    def disponible(self):
        """True si al menos el XML de la muestra esta accesible."""
        return os.path.exists(self.muestra)

    def faltantes(self):
        """Rutas del catalogo que no estan en disco."""
        return [r for r in (self.muestra, self.fondo)
                if r is not None and not os.path.exists(r)]

    def tamanos(self):
        """(bytes_muestra, bytes_fondo); None donde el archivo no esta."""
        def t(r):
            try:
                return os.path.getsize(r)
            except (OSError, TypeError):
                return None
        return t(self.muestra), t(self.fondo)

    def __repr__(self):
        return "<Corrida %s: %s>" % (self.clave, self.descripcion)


def duplicados(corridas=None):
    """
    Grupos de corridas cuya MUESTRA es el mismo archivo, byte a byte.

    Existe por un caso real (10-sep-2026): aparecieron tres `Co60.dat` en las
    carpetas del detector B que eran copias exactas de los del detector C
    -mismo tamano, misma fecha y mismo contenido-, y calibrarlos como si fueran
    del B habria dado tres numeros perfectamente formateados del detector
    equivocado. El catalogo no tiene como saber de que detector salio un
    archivo; lo unico que puede hacer es avisar cuando dos entradas apuntan al
    mismo contenido, que es una contradiccion: una medicion es de un detector.

    Barato a proposito, para poder llamarlo siempre: primero agrupa por tamano
    y fecha, que no cuesta nada, y solo a los que coinciden en las dos cosas
    les lee 4 MB del principio y 4 MB del final. Con archivos del DRS4 eso
    alcanza; no pretende ser una prueba criptografica sino un aviso.
    """
    import hashlib
    cs = list(corridas if corridas is not None else CORRIDAS.values())
    por_firma = {}
    for c in cs:
        if not c.muestra or not os.path.exists(c.muestra):
            continue
        st = os.stat(c.muestra)
        por_firma.setdefault((st.st_size, int(st.st_mtime)), []).append(c)

    def huella(ruta, n=4 << 20):
        m = hashlib.sha1()
        with open(ruta, "rb") as f:
            m.update(f.read(n))
            if os.path.getsize(ruta) > 2 * n:
                f.seek(-n, os.SEEK_END)
                m.update(f.read(n))
        return m.hexdigest()

    grupos = []
    for candidatos in por_firma.values():
        if len(candidatos) < 2:
            continue
        por_hash = {}
        for c in candidatos:
            por_hash.setdefault(huella(c.muestra), []).append(c)
        grupos += [g for g in por_hash.values() if len(g) > 1]
    return grupos


def misma_cadena(corridas):
    """
    Avisos por diferencias de cadena de senal entre corridas que se van a
    combinar. Lista vacia = se pueden combinar.

    El caso que esto tiene que atrapar: el area se integra a la ENTRADA del
    digitalizador, asi que un atenuador de 20 dB divide D por 10. Combinar una
    corrida atenuada con una que no lo esta da una "no linealidad" de un factor
    10 que no tiene nada que ver con el detector. Un campo sin registrar cuenta
    como aviso y no como coincidencia: no saber no es lo mismo que ser igual.

    Desde el 22-sep-2026 esta funcion no deberia avisar por el atenuador: las 25
    corridas estan registradas con 20 dB. Sigue existiendo para la corrida que
    venga con otro valor, que es cuando el aviso serviria de algo.
    """
    avisos = []
    for campo, unidad in (("atenuador", "dB"), ("muestreo", "GS/s")):
        vals = [getattr(c, campo) for c in corridas]
        if any(v is None for v in vals):
            faltan = [c.clave for c, v in zip(corridas, vals) if v is None]
            avisos.append("%s sin registrar en: %s" % (campo, ", ".join(faltan)))
        elif len(set(vals)) > 1:
            avisos.append("%s DISTINTO entre corridas: %s"
                          % (campo, ", ".join("%s = %g %s" % (c.clave, v, unidad)
                                              for c, v in zip(corridas, vals))))
    return avisos


SIM_CS137 = os.path.join(DATOS, "ComptonSpectrum_Cs137_10M.root")

CORRIDAS = {
    "co60": Corrida(
        "co60", "Co-60, Detector A a 1550 V (junio 2026)",
        muestra=_junio("A Co60 5GaS 50k 02mV trigg 1550 att.xml"),
        fondo=_junio("A BG 5GaS 10k 02mV trigg 1550 att.xml"),
        fuente="Co-60", detector="A", hv=1550, atenuador=20, umbral=0.2,
        muestreo=5,
        sim=os.path.join(DATOS, "ComptonSpectrum_Co60_3M.root"),
        nota="dos lineas gamma (1173.2 y 1332.5 keV) -> DOS bordes Compton "
             "solapados. El '1550' del nombre es el HV, no el trigger; el "
             "fotodetector NO esta confirmado (a 1550 V no puede ser el SiPM "
             "MicroFC-30020 que decia la etiqueta vieja: un MicroFC se "
             "polariza a ~25 V)"),
    "cs137_1600": Corrida(
        "cs137_1600", "Cs-137, Detector A a 1600 V (26-ago-2026)",
        muestra=_ext("A", "1600V", "Cs137.xml"),
        fondo=_ext("A", "1600V", "BG.xml"),
        fuente="Cs-137", detector="A", hv=1600, atenuador=20, muestreo=5,
        sim=SIM_CS137,
        nota="una sola linea (661.657 keV) -> un solo borde. Vivo el 0.13% "
             "(46 kHz de trigger real): el umbral estaba muy abajo. El "
             "espectro no tiene ni una cuenta debajo de 121 mV*ns y sube "
             "hasta ~334: eso es la eficiencia del trigger, no fisica"),
    "cs137_1650": Corrida(
        "cs137_1650", "Cs-137, Detector A a 1650 V (26-ago-2026)",
        muestra=_ext("A", "1650V", "Cs137.xml"),
        fondo=None, fuente="Cs-137", detector="A", hv=1650, atenuador=20,
        muestreo=5, sim=SIM_CS137,
        nota="sin fondo propio; el de 1600 V NO sirve, es otra ganancia"),
    "cs137_1700": Corrida(
        "cs137_1700", "Cs-137, Detector A a 1700 V (26-ago-2026)",
        muestra=_ext("A", "1700", "Cs137.xml"),
        fondo=None, fuente="Cs-137", detector="A", hv=1700, atenuador=20,
        muestreo=5, sim=SIM_CS137,
        nota="sin fondo propio; la carpeta se llama '1700', sin la V"),
}


# --------------------------------------------------------------------------
# Sesion del 01-sep-2026: cuatro detectores, tres voltajes, formato BINARIO.
#
# Se agregan por programa porque son 12 entradas con la misma forma. Ninguna
# tiene fondo propio: el unico BG que existe es el XML del detector A a 1600 V
# del 26-ago, que NO sirve para estas (otro dia, otro formato y, para B/C/D,
# otro detector). Sin fondo el espectro es el CRUDO y hay que decirlo.
#
# El Co-60 del detector A a 1600/1650/1700 V es lo que faltaba para la recta
# canal-energia: con el Cs-137 del 26-ago al mismo voltaje nominal quedan TRES
# bordes por punto de operacion (477.3, 963.4 y 1118.1 keV). Ojo con darlo por
# hecho: son dias distintos y formatos distintos, asi que hay que verificar que
# la cadena de senal sea la misma antes de combinarlos (`misma_cadena`).
# --------------------------------------------------------------------------
_NOTA_BIN = ("formato binario del DRS4 (sin extension en el nombre). Sin fondo "
             "propio. Fuente en un recipiente de plomo segun la bitacora. "
             "Atenuador de 20 dB, como todas las corridas del catalogo")

# La misma sin el recipiente de plomo, para las corridas de otras sesiones: de
# esas no consta, y "no consta" no es lo mismo que "no estaba" (rutas.misma_cadena
# hace la misma distincion con el atenuador).
_NOTA_BIN_SIN_BITACORA = (
    "formato binario del DRS4 (sin extension en el nombre). Sin fondo propio. "
    "Atenuador de 20 dB, como todas las corridas del catalogo. Del montaje de "
    "esta sesion no consta nada mas en la bitacora")


def fecha(ruta):
    """
    Cuando se midio ese archivo, leido de su fecha de modificacion.

    No se tabula a mano porque ya se equivoco una vez: el catalogo decia
    "01-sep-2026" para las 12 corridas binarias, y despues aparecieron el Co-60
    del detector C (02-sep) y el Co-60 y el Cs-137 del D (08-sep). La fecha
    importa para el test de linealidad: dos fuentes medidas en sesiones
    distintas pueden no estar en el mismo punto de operacion real aunque el
    voltaje nominal sea el mismo (CLAUDE.md 7.55), y eso hay que poder verlo.
    `shutil.copy2` conserva la fecha, asi que la copia local dice lo mismo que
    el disco.
    """
    try:
        return time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(ruta)))
    except OSError:
        return ""

SIM_CO60 = os.path.join(DATOS, "ComptonSpectrum_Co60_3M.root")

# (detector, fuente, archivo, sufijo de la clave, nota extra, nota base)
#
# El sufijo existe por el detector D, que tiene DOS corridas de Cs-137 por
# voltaje: `Cs137` y `Cs137 old`. Sin sufijo la segunda pisaria a la primera en
# el catalogo, que es la clase de error que no falla en ningun lado: quedaria
# una sola corrida y nadie se enteraria de cual.
#
# La nota base es None para las de la sesion del 1-sep (usan `_NOTA_BIN`); las
# de otras sesiones traen la suya, porque la bitacora que `_NOTA_BIN` cita es
# la de aquel dia.
_SESION_BIN = (
    ("A", "Co-60", "Co60", "", "", None),
    ("B", "Cs-137", "Cs137", "", "", None),
    ("C", "Cs-137", "Cs137", "", "", None),
    ("C", "Co-60", "Co60.dat", "", "", None),
    ("D", "Cs-137", "Cs137", "", "", None),
    ("D", "Co-60", "Co60", "", "", None),
    ("D", "Cs-137", "Cs137 old", "old",
     "DESCARTADA: en la sesion del 1-sep el detector D tenia la tierra mal "
     "conectada. Las corridas buenas del D son las del 8-sep. Esta se calibra "
     "igual, y se guarda, porque mide el tamano de esa falla (la ganancia "
     "cambio ~25%); no entra en ningun ajuste ni en la recta de calibracion",
     None),
    # 22-sep-2026: el Co-60 del detector B, que era la medicion que faltaba
    # para poder testear su linealidad (hasta hoy el B tenia una sola fuente).
    # Verificado antes de agregarlo, por el incidente del 10-sep: los tres
    # archivos son unicos byte a byte, no son copias de los del detector C ni
    # del Cs-137 del propio B (`rutas.duplicados`). Ojo al combinarlas con el
    # Cs-137 del B: son de SESIONES DISTINTAS (1-sep contra 22-sep), que es
    # exactamente la ambiguedad que dejo abierto el test A del detector A
    # (CLAUDE.md 7.55).
    ("B", "Co-60", "Co60", "", "", _NOTA_BIN_SIN_BITACORA),
)

for _det, _fuente, _archivo, _suf, _nota_extra, _base in _SESION_BIN:
    for _hv, _carpeta in ((1600, "1600V"), (1650, "1650V"), (1700, "1700V")):
        # el detector A guardo el de 1700 V en una carpeta llamada "1700",
        # y el B guardo su CS-137 de 1700 V en un archivo llamado "B1700V".
        # La condicion lleva la fuente a proposito: sin ella, el Co-60 del B a
        # 1700 V (22-sep) se resolveria al archivo del Cs-137 y se calibraria
        # la corrida equivocada sin que nada fallara.
        _c = "1700" if (_det == "A" and _hv == 1700) else _carpeta
        _a = ("B1700V" if (_det == "B" and _hv == 1700 and _fuente == "Cs-137")
              else _archivo)
        _clave = "%s_%s%s_%d" % (_det.lower(),
                                 _fuente.lower().replace("-", ""), _suf, _hv)
        _ruta = _ext(_det, _c, _a)
        CORRIDAS[_clave] = Corrida(
            _clave,
            "%s%s, Detector %s a %d V (%s, binario)"
            % (_fuente, " DESCARTADA" if _suf else "", _det, _hv,
               fecha(_ruta) or "fecha desconocida"),
            muestra=_ruta, fondo=None, fuente=_fuente,
            detector=_det, hv=_hv, atenuador=20, muestreo=5,
            sim=(SIM_CO60 if _fuente == "Co-60" else SIM_CS137),
            nota=((_base or _NOTA_BIN)
                  + ((". " + _nota_extra) if _nota_extra else "")))
del _det, _fuente, _archivo, _suf, _nota_extra, _base, _hv, _carpeta
del _c, _a, _clave, _ruta


# Na-22 (Detector C a 1600 V) queda FUERA del catalogo a proposito. No es que
# el archivo falte: esta en `datos/Detector C/1600V/Na22.dat`. Lo que falta son
# las dos cosas sin las que esta cadena no puede calibrarlo, y conviene que se
# vea aca y no como un error tres etapas mas adelante:
#   - `fisica.FUENTES` no tiene Na-22 (511.0 y 1274.5 keV), y
#   - no hay simulacion de Geant4 de esa fuente con esta geometria, asi que la
#     unica ruta posible seria la plantilla analitica, que sesga D un +4%
#     respecto de las demas (CLAUDE.md 7.18). Un punto con ese sesgo metido en
#     la recta de linealidad mediria la eleccion de plantilla, no el detector.
NA22_C_1600 = _local("C", "1600V", "Na22.dat")


# Las dos con las que se desarrolla el pipeline.
PRINCIPALES = ("co60", "cs137_1600")


def corrida(clave):
    """Devuelve la Corrida del catalogo, con un error claro si la clave no existe."""
    if clave not in CORRIDAS:
        raise KeyError("corrida desconocida: %r. Conocidas: %s"
                       % (clave, ", ".join(sorted(CORRIDAS))))
    return CORRIDAS[clave]


def descartada(clave):
    """
    True si la corrida del catalogo esta DESCARTADA como medicion valida.

    El unico caso hoy: los tres Cs-137 del detector D del 1-sep (`Cs137 old`),
    tomados con la tierra mal conectada; los buenos son los del 8-sep. Vive
    aca y no en `calibrar_todo` porque el registro tambien tiene que saberlo:
    ahi entraban al test A y se promediaban en la curva D(V) del detector D
    como si fueran una repeticion de la misma medicion.
    """
    return str(clave or "").endswith(("old_1600", "old_1650", "old_1700"))


def disponibles():
    """Claves de las corridas cuyo XML de muestra esta accesible ahora."""
    return [k for k, c in CORRIDAS.items() if c.disponible()]


def resolver(x):
    """
    Devuelve una `Corrida`, venga como clave del catalogo o ya como objeto.

    Existe para que el analisis no dependa del catalogo. La interfaz tiene que
    poder abrir los XML de un detector nuevo, que por definicion no esta en
    CORRIDAS; sin esto, cada funcion de arriba tendria que aprender a recibir
    rutas sueltas por separado.
    """
    return x if isinstance(x, Corrida) else corrida(x)


# ------------------------------------------------------------------ salida ---
def salida_utf8():
    """
    Hace que print() pueda escribir acentos.

    En este equipo `sys.stdout.encoding` es cp1252, asi que un print con 'o'
    acentuada o con el punto medio de 'mV*ns' sale mal segun la terminal. Se
    llama al principio de cada demo, nunca al importar el modulo: cambiarle el
    stdout a un programa que nos importa seria de mala educacion.
    """
    for flujo in (sys.stdout, sys.stderr):
        try:
            flujo.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


def _demo():
    salida_utf8()
    print("=" * 74)
    print("rutas.py - que ve el proyecto")
    print("=" * 74)
    print("  raiz    :", RAIZ)
    print("  datos   :", relativa(DATOS), "" if os.path.isdir(DATOS) else " [FALTA]")
    print("  figuras :", relativa(FIGURAS))
    print("  cache   :", relativa(CACHE))
    print("  sesion  :", SESION_EXTERNA,
          "" if os.path.isdir(SESION_EXTERNA) else " [NO MONTADA]")
    if not os.path.isdir(SESION_EXTERNA):
        print("            (se buscó la carpeta %r en todas las unidades; si el"
              % NOMBRE_SESION)
        print("             disco está en otro lado: set CALIBRACION_SESION=<ruta>)")
    print("  externo :", EXTERNO,
          "" if os.path.isdir(EXTERNO) else " [NO MONTADO]")
    print()
    print("Corridas del catalogo:")
    for clave in sorted(CORRIDAS):
        c = CORRIDAS[clave]
        tm, tf = c.tamanos()
        marca = "OK" if c.disponible() else "--"
        print()
        print("  [%s] %-12s %s" % (marca, clave, c.descripcion))
        print("       punto  : %s     cadena: %s" % (c.punto(), c.cadena()))
        print("       muestra: %s%s" % (
            relativa(c.muestra),
            "  (%.2f GB)" % (tm / 1e9) if tm else "   [FALTA]"))
        if c.fondo is None:
            print("       fondo  : (no hay fondo propio para esta corrida)")
        else:
            print("       fondo  : %s%s" % (
                relativa(c.fondo),
                "  (%.2f GB)" % (tf / 1e9) if tf else "   [FALTA]"))
        if c.sim:
            print("       sim    : %s%s" % (
                relativa(c.sim),
                "" if os.path.exists(c.sim) else "   [FALTA]"))
        if c.nota:
            print("       nota   : %s" % c.nota)
    print()
    hay = disponibles()
    print("Disponibles ahora: %s" % (", ".join(sorted(hay)) if hay else "NINGUNA"))
    faltan = sorted(k for k in CORRIDAS if k not in hay)
    if faltan:
        print("No accesibles    : %s" % ", ".join(faltan))
        print("  (si son las de Cs-137, conecta el disco externo D:)")


if __name__ == "__main__":
    _demo()
