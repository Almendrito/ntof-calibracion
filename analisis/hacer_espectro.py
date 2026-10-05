"""
hacer_espectro.py — el script de cabecera: XML -> espectro de energia.

Encadena las cuatro etapas sobre una corrida del catalogo y deja el espectro
guardado y dibujado. Es lo unico que hay que correr para producir el resultado;
los modulos sueltos son para entender cada paso.

    python analisis/hacer_espectro.py                      # co60, ventana medida
    python analisis/hacer_espectro.py --corrida cs137_1600
    python analisis/hacer_espectro.py --metodo ambos       # global vs convergencia
    python analisis/hacer_espectro.py --limite 20000       # rapido, para iterar
    python analisis/hacer_espectro.py --todas --limite 5000

La primera corrida sobre un archivo nuevo parsea el XML (~40 s por GB) y cachea
las areas en `analisis/.cache/`; las siguientes arrancan al instante.
"""

import argparse
import os
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, area, espectro


def _barra(frac, n_eventos=0):
    """
    Barra de progreso de una linea.

    La fraccion viene medida en BYTES leidos y no en eventos: cuantos eventos
    tiene el XML no se sabe hasta terminarlo, mientras que su tamano se sabe
    desde el principio (ver `drs_xml.iter_eventos`).
    """
    n = int(40 * frac)
    sys.stdout.write("\r    [%s%s] %5.1f%%   %s eventos" %
                     ("#" * n, "." * (40 - n), 100 * frac,
                      format(int(n_eventos), ",d")))
    sys.stdout.flush()


def procesar(clave, metodo="global", limite=None, bins=500, rango_max=None,
             progreso=True, dibujar=True):
    """Corre la cadena entera para una corrida y devuelve el/los espectros."""
    c = rutas.corrida(clave)
    if not c.disponible():
        print("  [!] %s no esta accesible: %s" % (clave, c.faltantes()))
        return {}

    print("=" * 78)
    print("%s  -  %s" % (clave, c.descripcion))
    print("     fuente: %s   detector: %s" % (c.fuente, c.detector))
    if c.nota:
        print("     nota  : %s" % c.nota)
    print("=" * 78)

    cb = _barra if progreso else None
    t0 = time.time()

    # La ventana se MIDE sobre la muestra y se usa la misma para el fondo.
    print("  [1/3] midiendo la ventana de integracion sobre el pulso promedio...")
    v = area.calibrar_ventana(c.muestra)

    metodos = ["global", "convergencia"] if metodo == "ambos" else [metodo]
    salida = {}
    for met in metodos:
        print("  [2/3] integrando los pulsos (metodo %s)..." % met)
        d_m = area.areas_de_archivo(c.muestra, metodo=met, pre=v["pre"],
                                    post=v["post"], limite=limite, callback=cb)
        if progreso:
            sys.stdout.write("\r" + " " * 70 + "\r")
        d_f = None
        if c.fondo and os.path.exists(c.fondo):
            d_f = area.areas_de_archivo(c.fondo, metodo=met, pre=v["pre"],
                                        post=v["post"], limite=limite,
                                        callback=cb)
            if progreso:
                sys.stdout.write("\r" + " " * 70 + "\r")

        print("  [3/3] histograma con binning comun y resta de fondo...")
        e = espectro.espectro_neto(d_m, d_f, bins=bins, rango_max=rango_max)
        # El binning comun no vale solo entre muestra y fondo: si cada METODO
        # elige su propio tope (el default sale de un percentil de sus areas),
        # los dos espectros quedan en ejes distintos y superponerlos es
        # mentira. El primero fija el eje y el resto lo hereda.
        if rango_max is None:
            rango_max = e["rango_max"]
        e["corrida"] = clave
        e["fuente"] = c.fuente or ""
        e["limite"] = int(limite) if limite else 0
        # El limite va en el NOMBRE: si no, una corrida rapida de prueba
        # sobrescribe el espectro completo con uno de 2000 eventos y no queda
        # ninguna senal de que paso.
        sufijo = "" if not limite else "_lim%d" % limite
        ruta = espectro.guardar(e, "espectro_%s_%s%s" % (clave, met, sufijo))
        print("  guardado -> %s" % rutas.relativa(ruta))
        salida[met] = e

    print("  (%.0f s en total)" % (time.time() - t0))

    if dibujar and salida:
        _dibujar(clave, salida, v)
    return salida


def _dibujar(clave, salida, vent):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    ax = axes[0]
    ax.plot(vent["tau_rel"], vent["forma"], lw=1.3)
    ax.axvline(-vent["pre"], color="C2", ls="--",
               label="t_ini = %.1f ns" % -vent["pre"])
    ax.axvline(vent["post"], color="C3", ls="--",
               label="t_fin = %+.1f ns" % vent["post"])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_title("la ventana, medida sobre %d pulsos" % vent["n_usados"])
    ax.set_xlabel("t - onset [ns]"); ax.set_ylabel("forma normalizada")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    e0 = next(iter(salida.values()))
    ax = axes[1]
    ax.step(e0["canales"], e0["muestra"], where="mid", lw=1.1, label="muestra")
    if e0["hay_fondo"]:
        ax.step(e0["canales"], e0["fondo"], where="mid", lw=1.1,
                label="fondo (norm. por tiempo VIVO, k=%.4f)" % e0["k_f"])
    ax.set_yscale("log")
    ax.set_title("%s: muestra y fondo, mismo eje" % clave)
    ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[2]
    for met, e in salida.items():
        ax.errorbar(e["canales"], e["neto"], yerr=e["error"], fmt="none",
                    ecolor="0.8", lw=0.5)
        ax.step(e["canales"], e["neto"], where="mid", lw=1.2, label=met)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("ESPECTRO NETO  (%s)" % (e0["fuente"] or clave))
    ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.tight_layout()
    ruta = rutas.figura("espectro_%s.png" % clave)
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    print("  figura   -> %s" % rutas.relativa(ruta))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corrida", default="co60",
                   help="clave del catalogo (ver `python analisis/rutas.py`)")
    p.add_argument("--todas", action="store_true",
                   help="procesa todas las corridas accesibles")
    p.add_argument("--metodo", default="global",
                   choices=["global", "convergencia", "ambos"])
    p.add_argument("--limite", type=int, default=None,
                   help="procesa solo los primeros N eventos")
    p.add_argument("--bins", type=int, default=500)
    p.add_argument("--rango-max", type=float, default=None,
                   help="tope del histograma [mV*ns]; por defecto 1.6*p99.5")
    p.add_argument("--sin-figura", action="store_true")
    a = p.parse_args(argv)

    rutas.salida_utf8()
    claves = sorted(rutas.disponibles()) if a.todas else [a.corrida]
    for clave in claves:
        procesar(clave, metodo=a.metodo, limite=a.limite, bins=a.bins,
                 rango_max=a.rango_max, dibujar=not a.sin_figura)
        print()


if __name__ == "__main__":
    main()
