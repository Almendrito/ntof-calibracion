"""
rutas.py — donde escribe esta carpeta.

Resuelve todo desde `__file__` para que los modulos corran desde cualquier
directorio de trabajo, igual que `analisis/rutas.py`. No repite el catalogo de
corridas: para eso `datos.py` importa el catalogo de `analisis`.
"""

import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))       # .../replica_paper
RAIZ = os.path.dirname(AQUI)                            # .../Calibracion CE
FIGURAS = os.path.join(AQUI, "figuras")
CACHE = os.path.join(AQUI, ".cache")


def _crear(carpeta):
    os.makedirs(carpeta, exist_ok=True)
    return carpeta


def figura(nombre):
    """Ruta de una figura de salida, creando la carpeta si falta."""
    return os.path.join(_crear(FIGURAS), nombre)


def cache(nombre):
    """Ruta de un archivo de cache de esta carpeta."""
    return os.path.join(_crear(CACHE), nombre)


def relativa_o(ruta):
    """Ruta relativa a la raiz del proyecto si se puede, absoluta si no."""
    try:
        r = os.path.relpath(ruta, RAIZ)
        return r if not r.startswith("..") else ruta
    except ValueError:
        return ruta


def salida_utf8():
    """Permite acentos en print(); el stdout de este equipo es cp1252."""
    for flujo in (sys.stdout, sys.stderr):
        try:
            flujo.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


def bootstrap():
    """
    Pone la RAIZ del proyecto en sys.path.

    Hace falta porque `python replica_paper/correr.py` pone en sys.path la
    carpeta del archivo, no la raiz, y entonces `from replica_paper import ...`
    no resolveria. Correr cada modulo solo es una caracteristica del proyecto.
    """
    if RAIZ not in sys.path:
        sys.path.insert(0, RAIZ)
