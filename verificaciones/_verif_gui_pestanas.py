"""
Que las pestanas de `analisis/gui.py` DIBUJEN, y que dibujen el resultado.

Por que existe
--------------
Las figuras de la interfaz se dibujan dentro de `Lienzo.dibujar`, que atrapa
cualquier excepcion y escribe "no se pudo dibujar: ..." en el propio grafico.
Eso es lo correcto para no matar la ventana por una figura, pero convierte un
error de programacion en un cartel gris que hay que ir a mirar pestana por
pestana. Asi estuvo roto el paisaje del chi2: `Problema.mapa` no sabia en que
valor congelar `lam`, tiraba TypeError, y la pestana mostraba el cartel.

Este archivo corre la cadena sobre una corrida del catalogo y llama a CADA
metodo de figura directamente, sin `Lienzo` de por medio, asi que una excepcion
se ve. Ademas guarda los PNG en `figuras/gui_*.png` para poder mirar lo que la
interfaz mostraria sin abrirla.

No necesita Tk vivo: se instancia un objeto minimo con los atributos de clase
que las figuras usan, y el lienzo es un canvas Agg. O sea que prueba las
figuras, no los widgets.

Si el archivo crudo de la corrida no esta a mano -las corridas de agosto y
septiembre viven en un disco externo que no siempre esta conectado- se arma lo
mismo desde el espectro cacheado. Lo unico que no se puede dibujar asi es la
ventana de integracion, que se mide sobre los pulsos, y esa pestana queda
salteada en vez de fallar.

    python verificaciones/_verif_gui_pestanas.py               # co60
    python verificaciones/_verif_gui_pestanas.py cs137_1600
"""

import os as _os
import sys as _sys
_RAIZ = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _RAIZ not in _sys.path:
    _sys.path.insert(0, _RAIZ)

import time

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from analisis import (rutas, espectro, simulacion, calibracion, pipeline)
from analisis import gui as g


# Las figuras solo usan atributos de CLASE (los diccionarios de nombres), asi
# que no hace falta una App de verdad -que pediria un root de Tk.
class _Stub:
    _NOMBRE_LIBRE = g.App._NOMBRE_LIBRE
    _CORTO = g.App._CORTO


# nombre del PNG, metodo, y que claves de `salida` necesita para poder dibujar
PESTANAS = [("resultado", g.App._fig_resultado, ()),
            ("ventana", g.App._fig_ventana, ("ventana",)),
            ("1_datos_y_simulacion", g.App._fig_espectro, ()),
            ("2_ncc", g.App._fig_ncc, ()),
            ("3_recocido", g.App._fig_recocido, ("sa",)),
            ("4_inicial_vs_final", g.App._fig_comparacion, ("sa",)),
            ("parametros", g.App._fig_parametros, ()),
            ("paisaje", g.App._fig_paisaje, ())]


def salida_de(clave, metodo="global"):
    """
    Todo lo que las pestanas necesitan, por la cadena completa o por el cache.

    Devuelve `(salida, como)`. `como` dice por que camino se armo, que hay que
    imprimir: un resultado que salio del espectro cacheado no prueba nada sobre
    el parseo del archivo.
    """
    c = rutas.corrida(clave)
    if c.muestra and _os.path.exists(c.muestra):
        return pipeline.calibrar_completo(c, {"metodo": metodo}), "cadena completa"

    nombre = "espectro_%s_%s" % (clave, metodo)
    if not _os.path.exists(rutas.cache(nombre + ".npz")):
        raise SystemExit(
            "no esta el archivo crudo de %s (%s) ni el espectro cacheado (%s).\n"
            "Con el disco externo conectado corre primero:\n"
            "   python analisis/hacer_espectro.py --corrida %s"
            % (clave, c.muestra, nombre, clave))
    e = espectro.cargar(nombre)
    plantilla = simulacion.cargar(c.sim, "Lvis_keVee")
    det = {}
    r, p, _, ncc = calibracion.calibrar(c, esp=e, plantilla=plantilla,
                                        metodo=metodo, verbose=False,
                                        detalles=det)
    salida = dict(corrida=c, opciones={}, espectro=e, plantilla=plantilla,
                  resultado=r, problema=p, ncc=ncc, sa=det,
                  lineas_sim=simulacion.lineas_primarias(c.sim))
    return salida, "espectro cacheado (el archivo crudo no esta montado)"


def main():
    rutas.salida_utf8()
    clave = _sys.argv[1] if len(_sys.argv) > 1 else "co60"
    c = rutas.corrida(clave)
    print("=" * 74)
    print("pestanas de la interfaz sobre %s  (%s)" % (clave, c.descripcion))
    print("=" * 74)

    t0 = time.time()
    salida, como = salida_de(clave)
    r = salida["resultado"]
    print("  %s en %.1f s   ->   D = %.4f +- %.4f, chi2/ndf = %.3f"
          % (como, time.time() - t0, r["D"], r["err_D"], r["chi2_ndf"]))

    # lo que las pestanas 3 y 4 necesitan y antes no salia de `calibrar`
    det = salida.get("sa") or {}
    assert det.get("corridas"), "no se devolvio el detalle del recocido"
    print("  detalle del recocido: %d arranques, semilla %s, chi2 semilla %.1f"
          % (len(det["corridas"]),
             ", ".join("%.4g" % v for v in det["semilla"]),
             det["chi2_semilla"]))
    for i, m in enumerate(det["corridas"]):
        assert m["historia_x"].shape[1] == len(det["libres"]), \
            "la historia del arranque %d no tiene los parametros del modelo" % i
    print("  camino guardado: %s ciclos por arranque"
          % ", ".join(str(m["ciclos"]) for m in det["corridas"]))

    stub = _Stub()
    fallas, hechas = [], 0
    for nombre, metodo, necesita in PESTANAS:
        faltan = [k for k in necesita if k not in salida]
        if faltan:
            print("  %-22s salteada (falta %s)" % (nombre, ", ".join(faltan)))
            continue
        fig = Figure(figsize=(15, 8.6), dpi=110)
        FigureCanvasAgg(fig)
        t = time.time()
        try:
            metodo(stub, fig, salida)
        except Exception as e:                                  # noqa: BLE001
            fallas.append((nombre, e))
            print("  %-22s FALLA: %s: %s" % (nombre, type(e).__name__, e))
            continue
        # el nombre lleva la corrida: si no, verificar una segunda
        # corrida pisa las figuras de la primera
        destino = rutas.figura("gui_%s_%s.png" % (clave, nombre))
        fig.savefig(destino)
        hechas += 1
        print("  %-22s ok  (%.1f s)  ->  %s"
              % (nombre, time.time() - t, rutas.relativa(destino)))

    print()
    if fallas:
        print("  %d de %d pestanas NO dibujan" % (len(fallas), len(PESTANAS)))
        for nombre, e in fallas:
            print("     %s: %r" % (nombre, e))
        return 1
    print("  las %d pestanas dibujan" % hechas)
    return 0


if __name__ == "__main__":
    _sys.exit(main())
