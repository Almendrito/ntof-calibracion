"""
revisar.py — primera mirada a una tanda de corridas, antes de calibrar nada.

Contesta las preguntas que hay que contestar ANTES de ajustar un borde:

  1. Se puede leer? Que formato, que placa, que canal, cuantos eventos.
  2. Cuanto tiempo estuvo viva la adquisicion, y a que tasa real disparo.
  3. La ventana de integracion ENTRA en el registro? Esta es la que se paso por
     alto en la sesion del 01-sep: el trigger quedo en t ~ 123 ns de un registro
     de 200 ns, asi que despues del onset sobran ~81 ns y el pulso necesita 83
     para converger al 99.9%. Con la ventana que pide `calibrar_ventana` el 84%
     de los eventos sale TRUNCADO.
  4. Que pinta tiene el espectro, donde termina la subida del umbral y donde
     cae el borde (estimacion cruda del 50% de altura, solo indicativa).

`python analisis/revisar.py`              revisa las corridas binarias
`python analisis/revisar.py --todas`      revisa todo el catalogo accesible
`python analisis/revisar.py a_co60_1600`  revisa las que se le nombren
"""

import argparse
import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "analisis"

from analisis import area, drs, espectro, fisica, pulso, rutas

# La compuerta (incluido el recorte al espacio que deja el registro) la mide
# `area.calibrar_ventana`: es donde vive esa regla y no se duplica aca.


def borde_crudo(canales, neto, ch_min=0.0, ventana=15, prominencia=0.05,
                piso=0.02):
    """
    Estimacion CRUDA del borde: 50% de la altura del hombro Compton.

    No es una calibracion -el 50% de altura sesga alto, medido en
    `codigo/etapas/clasicos.py`- sino un numero para saber si el espectro esta
    donde uno espera y si el binning tiene sentido.

    Dos trampas que hay que esquivar, las dos vistas en estos datos:

      - **El hombro no es el maximo del espectro.** El maximo esta pegado al
        umbral del trigger (subida de eficiencia + exceso de bajo canal).
        Quedarse con el da un "borde" al 30% del verdadero.
      - **El valle no es el minimo global.** Si se busca el minimo sobre todo el
        rango, cae en la cola vacia de arriba y el hombro se queda sin region
        donde vivir.

    Receta: recortar la cola vacia, buscar los MAXIMOS LOCALES con prominencia
    apreciable, y quedarse con el ULTIMO. Ese es el hombro Compton; el valle es
    el minimo entre el pico del umbral y el hombro.

    Devuelve `(canal del 50%, canal del hombro, canal del valle)`.
    """
    from scipy.signal import find_peaks, savgol_filter
    ch = np.asarray(canales, float)
    y = np.asarray(neto, float)
    v = int(ventana) | 1
    if y.size > v + 2:
        y = savgol_filter(y, v, 3)
    idx = np.nonzero(ch >= float(ch_min))[0]
    if idx.size < 10:
        return float("nan"), float("nan"), float("nan")

    i0 = idx[int(np.argmax(y[idx]))]                    # pico pegado al umbral
    encima = idx[(idx > i0) & (y[idx] > piso * y[i0])]  # recorte de la cola
    if encima.size:
        idx = idx[idx <= encima[-1]]
    if idx.size < 10:
        return float("nan"), float(ch[i0]), float("nan")

    tramo = y[idx]
    picos, _ = find_peaks(tramo, prominence=prominencia * float(tramo.max()))
    i_pico = idx[picos[-1]] if picos.size else i0
    if i_pico <= i0:                                    # no hay hombro aparte
        i_pico, i_valle = i0, np.nan
    else:
        entre = idx[(idx > i0) & (idx < i_pico)]
        i_valle = entre[int(np.argmin(y[entre]))] if entre.size else np.nan

    mitad = 0.5 * y[i_pico]
    i = int(i_pico)
    while i + 1 < ch.size and y[i] > mitad:
        i += 1
    return (float(ch[i]), float(ch[i_pico]),
            float(ch[i_valle]) if np.isfinite(i_valle) else float("nan"))


def revisar(clave, bins=600, rango_max=None, verbose=True):
    """Todo lo de arriba para una corrida del catalogo."""
    c = rutas.corrida(clave)
    if not c.disponible():
        return dict(clave=clave, error="no montado")

    info = drs.metadatos(c.muestra)
    v = area.calibrar_ventana(c.muestra, canal=info["canal"], verbose=False)
    post, post_pedido = v["post"], v["post_pedido"]
    esp, t_fin = v["espacio"], v["t_registro"]
    frac_entra, frac_carga = v["frac_eventos_enteros"], v["fraccion"]

    d_m = area.areas_de_archivo(c.muestra, pre=v["pre"], post=post,
                                canal=info["canal"], verbose=False)
    d_f = None
    if c.fondo and os.path.exists(c.fondo):
        d_f = area.areas_de_archivo(c.fondo, pre=v["pre"], post=post,
                                    verbose=False)
    e = espectro.espectro_neto(d_m, d_f, bins=bins, rango_max=rango_max,
                               verbose=False)
    ok = area.buenos(d_m)
    umbral = espectro.fin_del_umbral(e["canales"], e["neto"])
    b50, pico, valle = borde_crudo(e["canales"], e["neto"], ch_min=umbral)

    res = dict(clave=clave, corrida=c, info=info, ventana=v, post=post,
               post_pedido=post_pedido, frac_entra=frac_entra,
               espacio=esp, t_fin=t_fin, areas=d_m, espectro=e,
               frac_carga=frac_carga,
               n_utiles=int(ok.sum()), n_total=int(ok.size),
               umbral=umbral, borde50=b50, pico=pico, valle=valle,
               banderas={n: int(((d_m["bandera"] & b) > 0).sum())
                         for n, b in (("satura", pulso.MAL_SATURA),
                                      ("truncado", pulso.MAL_TRUNCADO),
                                      ("sin_pre", pulso.MAL_SIN_PRE),
                                      ("sin_pico", pulso.MAL_SIN_PICO))})
    if verbose:
        informe(res)
    return res


def informe(r):
    """Vuelca una revision en un bloque legible."""
    c, i, e, d = r["corrida"], r["info"], r["espectro"], r["areas"]
    print("=" * 78)
    print("%s  -  %s" % (r["clave"], c.descripcion))
    print("=" * 78)
    print("  formato   : %s   placa %s   canal %s   %d eventos"
          % (i.get("formato", "?"), i["board"], i["canal"],
             i.get("n_eventos", r["n_total"])))
    print("  registro  : %d muestras, %.1f ns, %.2f GS/s"
          % (i["muestras"], i["duracion_registro"], i["muestreo_GSs"]))
    print("  tiempo    : %.0f s de pared, %.0f s VIVOS (%.1f%%), scaler %.0f Hz"
          % (e["t_pared_muestra"], e["t_vivo_muestra"],
             100 * e["frac_viva_muestra"], e["scaler_muestra"]))
    print("             %.1f eventos/s guardados" % (d["n_eventos"]
                                                     / max(d["duracion"], 1e-9)))
    print("  COMPUERTA : pre = %.1f ns, post = %.1f ns -> recoge el %.2f%% de "
          "la carga" % (r["ventana"]["pre"], r["post"], 100 * r["frac_carga"]))
    if r["post"] < r["post_pedido"] - 0.05:
        p = np.percentile(r["espacio"], [5, 50, 95])
        print("             la ventana que pide el pulso es %.1f ns, pero el"
              % r["post_pedido"])
        print("             registro deja %.0f ns despues del onset (p5 %.0f, p95"
              % (p[1], p[0]))
        print("             %.0f), asi que se recorta a %.1f ns. Con eso entra el"
              % (p[2], r["post"]))
        print("             %.1f%% de los eventos." % (100 * r["frac_entra"]))
    print("  calidad   : %d/%d utiles (%.1f%%)   satura %d, truncado %d, "
          "sin_pre %d" % (r["n_utiles"], r["n_total"],
                          100 * r["n_utiles"] / max(r["n_total"], 1),
                          r["banderas"]["satura"], r["banderas"]["truncado"],
                          r["banderas"]["sin_pre"]))
    a = d["area"][area.buenos(d)]
    print("  area      : p1 %.0f   p50 %.0f   p90 %.0f   p99 %.0f   max %.0f mV*ns"
          % tuple(np.percentile(a, [1, 50, 90, 99, 100])))
    print("  espectro  : %d bins hasta %.0f mV*ns   fondo: %s"
          % (e["bins"], e["rango_max"], "restado" if e["hay_fondo"] else "NO HAY"))
    print("  umbral    : la subida del trigger termina en %.0f mV*ns" % r["umbral"])
    bordes = fisica.bordes(c.fuente) if c.fuente in fisica.FUENTES else ()
    print("  borde 50%% : %.0f mV*ns   (valle en %.0f, hombro en %.0f)"
          % (r["borde50"], r["valle"], r["pico"]))
    if bordes:
        print("              -> D ~ %.3f mV*ns/keV para el borde de %.1f keV"
              % (r["borde50"] / bordes[0], bordes[0]))
    print("              [crudo y sesgado alto; sirve para ubicarse, no para")
    print("               calibrar]")


def figura(revisiones, nombre="revision.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rs = [r for r in revisiones if "error" not in r]
    if not rs:
        return None
    n = len(rs)
    cols = min(3, n)
    filas = int(np.ceil(n / cols))
    fig, ax = plt.subplots(filas, cols, figsize=(4.8 * cols, 3.4 * filas),
                           squeeze=False)
    for k, r in enumerate(rs):
        a = ax[k // cols][k % cols]
        e = r["espectro"]
        a.step(e["canales"], e["neto"], where="mid", lw=0.8, color="0.35")
        a.axvline(r["umbral"], color="tab:blue", ls=":", lw=1,
                  label="fin del umbral")
        if np.isfinite(r["pico"]):
            a.axvline(r["pico"], color="tab:green", ls="--", lw=1, label="hombro")
        a.axvline(r["borde50"], color="tab:red", lw=1.2, label="borde 50% (crudo)")
        a.set_title("%s   %s" % (r["clave"], r["corrida"].punto()), fontsize=9)
        a.set_xlabel("area [mV*ns]")
        a.set_ylabel("cuentas")
        a.grid(alpha=0.3)
        a.legend(fontsize=6)
    for k in range(n, filas * cols):
        ax[k // cols][k % cols].axis("off")
    fig.tight_layout()
    s = rutas.figura(nombre)
    fig.savefig(s, dpi=125)
    return s


def _demo(claves=None, todas=False):
    rutas.salida_utf8()
    if claves:
        pass
    elif todas:
        claves = sorted(rutas.disponibles())
    else:
        claves = [k for k in sorted(rutas.disponibles())
                  if drs.formato(rutas.corrida(k).muestra) == "binario"]
    if not claves:
        print("No hay corridas accesibles. Corre `python analisis/rutas.py`.")
        return

    rs = []
    for k in claves:
        rs.append(revisar(k))
        print()
    s = figura(rs, "revision.png")

    print("=" * 78)
    print("Resumen")
    print("=" * 78)
    print("  corrida        punto           post  carga  utiles   viva    "
          "borde50    D ~     D/carga")
    for r in rs:
        if "error" in r:
            continue
        b = fisica.bordes(r["corrida"].fuente)[0]
        d = r["borde50"] / b
        print("  %-14s %-13s %5.1f %5.1f%% %5.1f%%  %5.1f%%   %7.0f   %6.3f   %6.3f"
              % (r["clave"], r["corrida"].punto(), r["post"],
                 100 * r["frac_carga"],
                 100 * r["n_utiles"] / max(r["n_total"], 1),
                 100 * r["espectro"]["frac_viva_muestra"],
                 r["borde50"], d, d / max(r["frac_carga"], 1e-9)))
    print()
    print("  'carga' es la fraccion del pulso que recoge la compuerta de esa")
    print("  corrida. La columna D/carga extrapola a carga completa, que es lo")
    print("  unico que hace comparables dos corridas con `post` distinto.")
    print("  'D ~' sale del 50% de altura, que SESGA ALTO: no es la calibracion.")
    if s:
        print()
        print("  figura: %s" % rutas.relativa(s))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("claves", nargs="*", help="claves del catalogo a revisar")
    ap.add_argument("--todas", action="store_true",
                    help="revisar todo el catalogo accesible")
    a = ap.parse_args()
    _demo(a.claves or None, a.todas)
