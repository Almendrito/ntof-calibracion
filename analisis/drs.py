"""
drs.py — elige el lector segun el archivo, y expone una sola interfaz.

El DRSOSC guarda en dos formatos y en este proyecto hay de los dos:

  - XML   (junio 2026 y sesion del 26-ago-2026)  -> `drs_xml`
  - binario "Save as binary" (sesion del 1-sep)  -> `drs_bin`

Es el mismo instrumento y el mismo pulso; cambia el envase. Todo lo que hay
aguas abajo (`area`, `pulso`, `espectro`, `pipeline`, la GUI) pide los eventos
por aca y no se entera.

La eleccion se hace mirando los primeros 4 BYTES del archivo, no la extension:
los archivos binarios de la sesion del 1-sep no tienen extension (se llaman
`Co60`, `Cs137`, `B1700V`), asi que decidir por el nombre habria sido decidir
por nada. Un archivo que no empieza con "DRS2"/"DRS8" se trata como XML, que es
lo que era todo hasta ahora.
"""

import os
import sys

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "analisis"

from analisis import drs_bin, drs_xml


def formato(ruta):
    """'binario' o 'xml'. Mira el contenido, no el nombre."""
    return "binario" if drs_bin.es_binario(ruta) else "xml"


def lector(ruta):
    """El modulo que sabe leer este archivo."""
    return drs_bin if drs_bin.es_binario(ruta) else drs_xml


# --- la interfaz, delegada -------------------------------------------------
def detectar_canal(ruta, **kw):
    return lector(ruta).detectar_canal(ruta, **kw)


def iter_eventos(ruta, **kw):
    return lector(ruta).iter_eventos(ruta, **kw)


def leer_pulsos(ruta, **kw):
    return lector(ruta).leer_pulsos(ruta, **kw)


def metadatos(ruta):
    d = lector(ruta).metadatos(ruta)
    d.setdefault("formato", "XML del DRSOSC")
    return d


def imprimir_metadatos(info):
    mod = drs_bin if str(info.get("formato", "")).startswith("binario") else drs_xml
    return mod.imprimir_metadatos(info)


def _demo():
    from analisis import rutas
    rutas.salida_utf8()
    print("=" * 74)
    print("drs.py - que formato tiene cada corrida del catalogo")
    print("=" * 74)
    for clave in sorted(rutas.CORRIDAS):
        c = rutas.corrida(clave)
        if not c.disponible():
            print("  %-16s [no montado]" % clave)
            continue
        f = formato(c.muestra)
        try:
            canal, sc = detectar_canal(c.muestra)
        except Exception as e:                            # noqa: BLE001
            canal, sc = "?", str(e)[:40]
        print("  %-16s %-8s  %-6s  scaler %s   %s"
              % (clave, f, canal, sc, os.path.basename(c.muestra)))


if __name__ == "__main__":
    _demo()
