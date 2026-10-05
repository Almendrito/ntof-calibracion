"""
datos.py — el UNICO puente con el resto del proyecto.

El diagrama de flujo del paper (Fig. 2) arranca en "Energy Spectrum (Channel,
Counts)": como se llega desde el XML del DRS4 hasta ese histograma no es parte
del paper. Asi que ese tramo se toma prestado de `analisis/` (area con
compuerta, tiempo vivo, resta de fondo con binning comun) en vez de reescribirlo
peor, y queda aislado aca para que se vea exactamente donde termina lo prestado
y donde empieza la replica.

Lo que SI es de este archivo, porque el paper lo pide y el proyecto no lo tenia:

  "Ten independent pulse-height spectrum measurements per scintillator-source
   combination yielded 80 datasets ... error bars representing standard
   deviations across measurements."

Nosotros tenemos UNA corrida por fuente, no diez. `bloques()` parte los eventos
de una corrida en 10 tramos consecutivos y los trata como las diez mediciones.
OJO CON QUE PRUEBA ESO: la dispersion entre bloques es solo la ESTADISTICA DE
CONTEO. Las diez mediciones del paper se tomaron por separado, asi que sus
barras incluyen ademas la variacion de montaje, de posicion de la fuente y de
deriva de ganancia. La barra de aca es una COTA INFERIOR de la de ellos.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import rutas as _rutas

_rutas.bootstrap()

from analisis import area as _area          # noqa: E402
from analisis import espectro as _espectro  # noqa: E402
from analisis import rutas as _cat          # noqa: E402

# Lineas gamma de las fuentes medidas [keV]. Las mismas que usa `analisis`.
LINEAS = {
    "Co-60": (1173.228, 1332.492),
    "Cs-137": (661.657,),
}

# Binning por corrida. Se fija a mano y NO se deja en automatico: unos pocos
# eventos de area enorme (pile-up, cosmicos) estiran el rango automatico y
# dejan toda la fisica en los primeros bins.
# El numero de bins importa para replicar la Fig. 5: los intervalos mas angostos
# del paper cubren el 6% del canal del borde, y sobre sus espectros eso son ~100
# bins. Con 800 bins nuestros espectros dan ~30 en esa ventana. Menos que eso y
# el ajuste no tiene con que; el informe lo dice intervalo por intervalo.
BINNING = {
    "co60": dict(bins=800, rango_max=2500.0),
    "cs137_1600": dict(bins=800, rango_max=1600.0),
    "cs137_1650": dict(bins=800, rango_max=2000.0),
    "cs137_1700": dict(bins=800, rango_max=2600.0),
}


def _areas(corrida, verbose=False):
    """Areas por evento de la muestra y del fondo, con la MISMA ventana."""
    v = _area.calibrar_ventana(corrida.muestra, verbose=verbose)
    d_m = _area.areas_de_archivo(corrida.muestra, pre=v["pre"], post=v["post"],
                                 verbose=verbose)
    d_f = None
    if corrida.fondo and os.path.exists(corrida.fondo):
        d_f = _area.areas_de_archivo(corrida.fondo, pre=v["pre"], post=v["post"],
                                     verbose=verbose)
    return d_m, d_f


def espectro(clave, bins=None, rango_max=None, verbose=False):
    """
    (canal, cuentas) de una corrida del catalogo, mas su metadata.

    `cuentas` es el NETO (muestra - fondo normalizado). Se conservan los
    negativos: son fluctuaciones legitimas y recortarlas sesga el ajuste.
    """
    c = _cat.corrida(clave)
    cfg = dict(BINNING.get(clave, dict(bins=400, rango_max=None)))
    if bins is not None:
        cfg["bins"] = bins
    if rango_max is not None:
        cfg["rango_max"] = rango_max
    d_m, d_f = _areas(c, verbose=verbose)
    e = _espectro.espectro_neto(d_m, d_f, verbose=verbose, **cfg)
    return dict(clave=clave, corrida=c, fuente=c.fuente,
                lineas=LINEAS[c.fuente], punto=c.punto(),
                canales=e["canales"], cuentas=e["neto"], error=e["error"],
                muestra=e["muestra"], fondo=e["fondo"],
                bins=cfg["bins"], rango_max=e["rango_max"],
                n_pulsos=e["n_muestra"], hay_fondo=e["hay_fondo"],
                k_m=e["k_m"], k_f=e["k_f"])


def _trozo(d, i0, i1, j0, j1):
    """Sub-dict de areas con los pulsos [i0:i1) y los eventos [j0:j1)."""
    out = {}
    for k, v in d.items():
        if isinstance(v, np.ndarray) and v.ndim == 1:
            out[k] = v[i0:i1] if v.size == d["area"].size else v[j0:j1]
        else:
            out[k] = v
    rel = out.get("reloj_s")
    if rel is not None and rel.size >= 2:
        out["duracion"] = float(rel[-1] - rel[0])
    else:
        out["duracion"] = float(d["duracion"]) * (j1 - j0) / max(d["n_eventos"], 1)
    out["n_eventos"] = int(j1 - j0)
    return out


def bloques(clave, n=10, bins=None, rango_max=None, verbose=False):
    """
    Las "diez mediciones independientes" del paper, hechas con lo que hay.

    Parte los eventos de la MUESTRA en `n` tramos consecutivos (el fondo entero
    se resta de cada uno, escalado por el tiempo vivo del tramo) y devuelve una
    lista de espectros con el MISMO binning que `espectro()`.
    """
    c = _cat.corrida(clave)
    cfg = dict(BINNING.get(clave, dict(bins=400, rango_max=None)))
    if bins is not None:
        cfg["bins"] = bins
    if rango_max is not None:
        cfg["rango_max"] = rango_max
    d_m, d_f = _areas(c, verbose=verbose)

    if cfg.get("rango_max") is None:
        cfg["rango_max"] = float(np.percentile(d_m["area"], 99.5) * 1.6)

    npul = int(d_m["area"].size)
    neve = int(d_m["reloj_s"].size)
    cortes_p = np.linspace(0, npul, int(n) + 1).astype(int)
    cortes_e = np.linspace(0, neve, int(n) + 1).astype(int)

    out = []
    for b in range(int(n)):
        sub = _trozo(d_m, cortes_p[b], cortes_p[b + 1],
                     cortes_e[b], cortes_e[b + 1])
        e = _espectro.espectro_neto(sub, d_f, verbose=False, **cfg)
        out.append(dict(clave=clave, bloque=b, corrida=c, fuente=c.fuente,
                        lineas=LINEAS[c.fuente], punto=c.punto(),
                        canales=e["canales"], cuentas=e["neto"],
                        error=e["error"], bins=cfg["bins"],
                        rango_max=e["rango_max"],
                        n_pulsos=int(cortes_p[b + 1] - cortes_p[b])))
    return out


def disponibles():
    """Claves del catalogo cuyo XML esta accesible ahora."""
    return _cat.disponibles()


def con_fondo():
    """
    Claves accesibles que ademas tienen fondo propio.

    Sin fondo no hay espectro NETO, y el paper trabaja sobre el espectro medido
    sin decir nada de la resta: usar una corrida sin su fondo mezclaria cosmicos
    y ruido con la fisica sin que se note.
    """
    return [k for k in disponibles() if _cat.corrida(k).fondo]


def corrida(clave):
    """La `Corrida` del catalogo (para no importar `analisis.rutas` afuera)."""
    return _cat.corrida(clave)


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _rutas.salida_utf8()
    print("=" * 74)
    print("datos.py - de que espectros dispone la replica")
    print("=" * 74)
    hay = disponibles()
    fig, ax = plt.subplots(1, len(hay) or 1, figsize=(5.2 * max(len(hay), 1), 4.0),
                           squeeze=False)
    for k, clave in enumerate(hay):
        e = espectro(clave)
        print("  %-12s %-8s  %6d pulsos  %3d bins hasta %.0f mV*ns  fondo: %s"
              % (clave, e["fuente"], e["n_pulsos"], e["bins"], e["rango_max"],
                 "si" if e["hay_fondo"] else "NO"))
        print("               lineas gamma %s -> bordes %s keV"
              % (", ".join("%.1f" % x for x in e["lineas"]),
                 ", ".join("%.1f" % (2 * (x / 511.0) * x / (1 + 2 * x / 511.0))
                           for x in e["lineas"])))
        a = ax[0][k]
        a.step(e["canales"], e["cuentas"], where="mid", lw=0.9)
        a.set_title("%s (%s)" % (clave, e["punto"]))
        a.set_xlabel("canal [mV*ns]")
        a.set_ylabel("cuentas netas")
        a.grid(alpha=0.3)
    fig.tight_layout()
    s = _rutas.figura("datos_espectros.png")
    fig.savefig(s, dpi=130)
    print("  figura: %s" % s)


if __name__ == "__main__":
    _demo()
