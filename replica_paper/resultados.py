"""
resultados.py — el analisis de resultados del paper (Seccion II-C y III).

Lo que el paper hace con los bordes una vez encontrados, y nada mas:

  1. Recta canal-energia ("Fit Energy-Band Structure Relationship"), reportada
     con pendiente, ordenada, SD y R^2. La definicion de SD no esta escrita:
     se dedujo reproduciendo la Fig. 3 desde sus propios puntos publicados y
     resulta ser el RMS de los residuos (ddof = 0). Ver `auditoria.py`.
  2. Similitud coseno entre relaciones energia-canal. El paper reporta 99.9993%
     (LaBr3) y 99.9999% (NaI) y NUNCA dice entre que vectores. Aca estan las
     lecturas razonables, todas, con el aviso de que ninguna reproduce los dos
     numeros.
  3. Energias reconstruidas desde la recta, comparadas con las teoricas ->
     error relativo por punto, con barra = desviacion estandar entre mediciones.
  4. Reduccion del error medio de NCC-SA contra cada metodo convencional (los
     80%/85% y 44%/66% del resumen).
  5. Robustez frente al intervalo de ajuste (Fig. 5).

La recta se ajusta SIEMPRE dentro de un mismo punto de operacion (detector +
alto voltaje). Mezclar voltajes distintos no da una recta canal-energia: da la
ley de ganancia del fotodetector, que es otra medicion.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"


# ----------------------------------------------------------------- 1. recta --
def ajuste_lineal(canal, energia):
    """
    Recta energia = m*canal + b, con SD y R^2 como los reporta la Fig. 3.

    SD = sqrt(mean(residuos^2)). No es una eleccion: es lo que reproduce los
    cinco decimales impresos en los dos paneles de la Fig. 3 del paper.
    R^2 = 1 - SS_res/SS_tot.

    Con 2 puntos la recta es exacta (0 grados de libertad): SD = 0 y R^2 = 1 no
    significan nada, y por eso salen marcados con `gdl`.
    """
    x = np.asarray(canal, float)
    y = np.asarray(energia, float)
    if x.size != y.size:
        raise ValueError("canal y energia tienen que tener el mismo largo")
    if x.size < 2:
        raise ValueError("hacen falta al menos 2 puntos para una recta")
    m, b = np.polyfit(x, y, 1)
    res = y - (m * x + b)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    return dict(m=float(m), b=float(b), sd=float(np.sqrt(np.mean(res ** 2))),
                r2=float(r2), residuos=res, n=int(x.size),
                gdl=int(x.size - 2), canal=x, energia=y)


def energia_reconstruida(recta, canal):
    """E = m*canal + b."""
    return recta["m"] * np.asarray(canal, float) + recta["b"]


# ------------------------------------------------------- 2. similitud coseno --
def similitud_coseno(a, b):
    """cos(a, b) = a.b / (|a| |b|)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def similitudes(recta1, recta2, canales=None):
    """
    Las lecturas posibles de "cosine similarity of the energy-channel
    relationship", porque el paper no dice cual usa.

    - "coeficientes": entre los vectores (m, b) de las dos rectas.
    - "en los canales medidos": entre las energias que las dos rectas
      reconstruyen en los canales donde efectivamente hay puntos.
    - "en toda la escala": lo mismo sobre una grilla de canales de 0 al maximo.
    """
    ch = (np.asarray(canales, float) if canales is not None
          else np.asarray(recta1["canal"], float))
    grilla = np.linspace(0.0, float(max(ch.max(), 1.0)), 2001)
    return dict(
        coeficientes=similitud_coseno([recta1["m"], recta1["b"]],
                                      [recta2["m"], recta2["b"]]),
        en_los_canales=similitud_coseno(energia_reconstruida(recta1, ch),
                                        energia_reconstruida(recta2, ch)),
        en_toda_la_escala=similitud_coseno(energia_reconstruida(recta1, grilla),
                                           energia_reconstruida(recta2, grilla)))


# ----------------------------------------------------- 3. errores relativos --
def error_relativo(e_reconstruida, e_teorica):
    """Error relativo en %, con el signo conservado."""
    er = np.asarray(e_reconstruida, float)
    et = np.asarray(e_teorica, float)
    return 100.0 * (er - et) / et


def errores_de_metodo(canales, sd_canales, energias_teoricas, recta=None):
    """
    La tabla de la Fig. 4 para UN metodo: error relativo por punto y su barra.

    La barra sale propagando la dispersion del canal entre mediciones a traves
    de la recta: sigma_E = |m| * sd_canal, y en % sobre la energia teorica.
    Es la lectura literal de "error bars representing standard deviations across
    measurements" cuando lo que se repite es la medicion del canal.
    """
    ch = np.asarray(canales, float)
    sd = np.asarray(sd_canales, float)
    et = np.asarray(energias_teoricas, float)
    r = ajuste_lineal(ch, et) if recta is None else recta
    er = energia_reconstruida(r, ch)
    return dict(recta=r, canal=ch, e_teorica=et, e_reconstruida=er,
                error_rel=error_relativo(er, et),
                barra_rel=100.0 * np.abs(r["m"]) * sd / et)


def reduccion_media(error_referencia, error_nuevo):
    """
    "reduced mean errors by X% compared to ..." — sobre el error absoluto medio.

        reduccion = (media|e_ref| - media|e_nuevo|) / media|e_ref| * 100
    """
    a = float(np.mean(np.abs(np.asarray(error_referencia, float))))
    b = float(np.mean(np.abs(np.asarray(error_nuevo, float))))
    if a <= 0:
        return float("nan")
    return 100.0 * (a - b) / a


# ---------------------------------------------------------- 4. agrupamiento --
def agrupar_por_punto(puntos):
    """
    Agrupa los puntos (canal, energia) por punto de operacion.

    `puntos` es una lista de dicts con al menos: punto, metodo, canal, energia.
    La recta canal-energia SOLO tiene sentido dentro de un punto de operacion.
    """
    grupos = {}
    for p in puntos:
        grupos.setdefault((p["punto"], p["metodo"]), []).append(p)
    return grupos


def rectas_por_punto(puntos):
    """
    Una recta por (punto de operacion, metodo). Los grupos con menos de 2 puntos
    no se ajustan: se devuelven con `recta = None` y su motivo.
    """
    salida = {}
    for clave, lista in agrupar_por_punto(puntos).items():
        lista = sorted(lista, key=lambda p: p["canal"])
        ch = [p["canal"] for p in lista]
        en = [p["energia"] for p in lista]
        if len(ch) < 2:
            salida[clave] = dict(puntos=lista, recta=None,
                                 motivo="un solo borde medido en este punto de "
                                        "operacion: no hay recta que ajustar")
            continue
        r = ajuste_lineal(ch, en)
        motivo = None
        if r["gdl"] == 0:
            motivo = ("2 puntos = 0 grados de libertad: la recta pasa exacta, "
                      "asi que SD = 0 y R^2 = 1 no miden nada")
        salida[clave] = dict(puntos=lista, recta=r, motivo=motivo)
    return salida


# ------------------------------------------------------------ 5. robustez ----
def robustez_intervalo(funcion_borde, intervalos):
    """
    Fig. 5: el mismo espectro, distintos intervalos de ajuste.

    `funcion_borde(intervalo) -> canal del borde`. Devuelve la lista de
    (intervalo, borde) y las estadisticas con que el paper la resume
    cualitativamente ("the calibration results remained highly consistent").
    """
    filas = []
    for iv in intervalos:
        try:
            filas.append((tuple(float(v) for v in iv), float(funcion_borde(iv))))
        except Exception as e:                       # noqa: BLE001
            filas.append((tuple(float(v) for v in iv), float("nan")))
            del e
    v = np.array([b for _, b in filas], float)
    ok = np.isfinite(v)
    return dict(filas=filas, valores=v,
                media=float(np.mean(v[ok])) if ok.any() else float("nan"),
                sd=float(np.std(v[ok], ddof=0)) if ok.any() else float("nan"),
                rango=float(v[ok].max() - v[ok].min()) if ok.any() else float("nan"),
                dispersion_rel=(float((v[ok].max() - v[ok].min()) / np.mean(v[ok]))
                                if ok.any() and np.mean(v[ok]) != 0 else float("nan")))


# ------------------------------------------------------------------- demo ----
def _demo():
    from replica_paper import rutas
    rutas.salida_utf8()
    print("=" * 74)
    print("resultados.py - las herramientas del analisis del paper")
    print("=" * 74)
    print("  Se validan contra los datos publicados del propio paper.")
    print("  Correr `python replica_paper/auditoria.py` para ver eso.")
    print()
    ch = [239, 314, 469, 518, 543]
    en = [477, 639, 963, 1061, 1117]
    r = ajuste_lineal(ch, en)
    print("  LaBr3, bordes Compton (puntos leidos de la Fig. 3a):")
    print("    y = %.4f x %+.3f    SD = %.6f   R2 = %.6f   (n = %d, gdl = %d)"
          % (r["m"], r["b"], r["sd"], r["r2"], r["n"], r["gdl"]))
    print("    el paper imprime  y = 2.10x - 21.73")


if __name__ == "__main__":
    _demo()
