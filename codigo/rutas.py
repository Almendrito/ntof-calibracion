"""
rutas.py - donde vive cada cosa del proyecto.

Por que existe
--------------
Hasta la reorganizacion del 27-ago-2026 los archivos de datos se nombraban
relativos al directorio de trabajo ("A Co60 ....xml"), asi que el codigo solo
funcionaba si uno estaba parado exactamente en la carpeta del proyecto. Con las
carpetas separadas eso deja de alcanzar, y ademas conviene que no importe desde
donde se ejecute: todo se resuelve desde `__file__`.

La estructura que asume:

    Calibracion CE/            <- RAIZ
    |- codigo/                 <- AQUI (este archivo)
    |   |- rutas.py
    |   |- gui.py  calibracion.py  documentacion.py
    |   |- etapas/     physics, model, fit_single, fit_double, ncc, sa
    |   |- datos_io/   spectrum
    |   |- ajuste/     sim_template, modelo_ext, mcmc
    |   `- scripts/    calibrar_v2, calibrar_ncc, calibrar
    |- datos/          los .xml del DRSOSC, sus caches .npz y los .root
    |- figuras/        todo lo que los modulos dibujan al correr solos
    |- docs/           el paper de referencia y notas historicas
    |- verificaciones/ scripts sueltos de comprobacion
    `- .cache/         caches de la GUI (preferencias, MCMC, registro)

Los nombres de los archivos de datos que estan aca son solo los DEFAULTS con
que arranca la interfaz. Cualquier otro archivo se elige en la pestana
"Archivos", y de ahi en mas manda esa eleccion.
"""

import os

AQUI = os.path.dirname(os.path.abspath(__file__))       # .../codigo
RAIZ = os.path.dirname(AQUI)                            # .../Calibracion CE

DATOS = os.path.join(RAIZ, "datos")
FIGURAS = os.path.join(RAIZ, "figuras")
DOCS = os.path.join(RAIZ, "docs")
CACHE = os.path.join(RAIZ, ".cache")


def dato(nombre):
    """Ruta a un archivo de `datos/`. No verifica que exista."""
    return os.path.join(DATOS, nombre)


def figura(nombre):
    """Ruta donde guardar una figura, creando `figuras/` si hace falta."""
    os.makedirs(FIGURAS, exist_ok=True)
    return os.path.join(FIGURAS, nombre)


def cache(*partes):
    """Ruta dentro de `.cache/`, creando la carpeta si hace falta."""
    ruta = os.path.join(CACHE, *partes)
    os.makedirs(os.path.dirname(ruta) or CACHE, exist_ok=True)
    return ruta


def relativa(ruta):
    """
    La ruta escrita relativa a la raiz del proyecto, si esta adentro.

    Es para MOSTRAR (etiquetas de la GUI, preferencias guardadas): una ruta
    absoluta de OneDrive con espacios y acentos no se lee, y guardar la
    relativa hace que el proyecto siga funcionando si se mueve de carpeta.
    """
    try:
        rel = os.path.relpath(ruta, RAIZ)
        return ruta if rel.startswith("..") else rel
    except ValueError:            # otra unidad de disco
        return ruta


def absoluta(ruta):
    """Inversa de `relativa`: resuelve contra la raiz lo que sea relativo."""
    if not ruta:
        return ruta
    return ruta if os.path.isabs(ruta) else os.path.normpath(
        os.path.join(RAIZ, ruta))


# --- los defaults con que arranca la interfaz ---
MUESTRA = dato("A Co60 5GaS 50k 02mV trigg 1550 att.xml")
FONDO = dato("A BG 5GaS 10k 02mV trigg 1550 att.xml")
SIM = dato("ComptonSpectrum_Co60_3M.root")
SIM_CS137 = dato("ComptonSpectrum_Cs137_10M.root")
ARBOL_SIM = "ComptonSpectrum"


if __name__ == "__main__":
    print("RAIZ    :", RAIZ)
    for nombre in ("DATOS", "FIGURAS", "DOCS", "CACHE"):
        ruta = globals()[nombre]
        print("%-8s: %s   %s" % (nombre, relativa(ruta),
                                 "" if os.path.isdir(ruta) else "(no existe)"))
    print()
    for nombre in ("MUESTRA", "FONDO", "SIM", "SIM_CS137"):
        ruta = globals()[nombre]
        print("%-10s %-58s %s"
              % (nombre, os.path.basename(ruta),
                 "ok" if os.path.exists(ruta) else "FALTA"))
