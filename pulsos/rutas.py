"""
rutas.py — donde vive cada cosa DENTRO de `pulsos/`.

El catalogo de corridas (que XML es cada fuente, a que HV, etc.) es
responsabilidad de `analisis/`, no se duplica: se re-exporta aca para que
`pulsos/*.py` haga `from pulsos import rutas` y tenga ahi mismo tanto sus
propias carpetas de salida como el catalogo, igual que el resto del proyecto.

Correr `python pulsos/rutas.py` lista lo que ve.
"""

import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))          # .../pulsos
RAIZ = os.path.dirname(AQUI)                                # .../Calibracion CE
FIGURAS = os.path.join(AQUI, "figuras")
CACHE = os.path.join(AQUI, ".cache")

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, RAIZ)

from analisis.rutas import (disponibles, corrida, resolver, CORRIDAS,
                            Corrida, salida_utf8)


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


def _demo():
    salida_utf8()
    print("=" * 74)
    print("pulsos/rutas.py - que ve el paquete")
    print("=" * 74)
    print("  raiz    :", RAIZ)
    print("  figuras :", relativa(FIGURAS))
    print("  cache   :", relativa(CACHE))
    print()
    hay = disponibles()
    print("Corridas disponibles (catalogo de analisis/rutas.py): %s"
          % (", ".join(sorted(hay)) if hay else "NINGUNA"))
    faltan = sorted(k for k in CORRIDAS if k not in hay)
    if faltan:
        print("No accesibles                                   : %s"
              % ", ".join(faltan))


if __name__ == "__main__":
    _demo()
