"""
lote.py — cargar N pulsos de un XML y calcular area + FWHM de todos.

Conecta `analisis/drs_xml.py` (lectura), `analisis/pulso.py` (area, metodo
GLOBAL) y `pulsos/ancho.py` (FWHM) sobre un LOTE que vive en memoria, no sobre
el archivo entero: esta herramienta es para inspeccionar cientos o miles de
pulsos, no para recorrer los millones de eventos de una corrida completa (eso
es lo que hace `analisis/area.py`).

Correr `python pulsos/lote.py [n]` deja `figuras/01_lote.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import drs_xml, pulso
from pulsos import rutas, ancho


def cargar_pulsos(ruta, n=500, canal=None, callback=None, cada=200):
    """
    Apila los primeros `n` eventos UTILES en memoria: `(T, V, relojes, scalers, canal)`.

    Se llama directo a `drs_xml.iter_eventos` (el generador primitivo del que
    esta hecho `drs_xml.leer_pulsos`) y no a `leer_pulsos`, porque esa funcion
    no acepta `callback` y aca la barra de progreso importa: el XML puede
    pesar gigabytes aunque solo se lean los primeros `n` eventos.

    Eventos de largo distinto al del primero se descartan (el volcado a XML
    a veces se corta), igual criterio que `leer_pulsos`.
    """
    if canal is None:
        canal, _ = drs_xml.detectar_canal(ruta)
    T, V, relojes, scalers, largo = [], [], [], [], None
    for t, v, reloj, scaler in drs_xml.iter_eventos(
            ruta, canal=canal, limite=n, callback=callback, cada=cada):
        if t.size == 0:
            continue
        if largo is None:
            largo = t.size
        if t.size != largo:
            continue
        T.append(t); V.append(v); relojes.append(reloj); scalers.append(scaler)
    if not T:
        raise ValueError("%s: no pude leer ni un pulso util (canal %s)"
                         % (os.path.basename(ruta), canal))
    return (np.array(T), np.array(V), relojes,
            np.array(scalers, dtype=float), canal)


def sugerir_ventana(T, V):
    """
    Sugiere `(pre, post)` sobre pulsos YA cargados en memoria, sin releer el
    archivo. Misma receta que `analisis.area.calibrar_ventana` (`pulso_promedio`
    -> ruido de la forma -> `ventana_por_convergencia`), pero esa funcion
    vuelve a abrir el XML desde disco; aca se reusan directo las piezas de
    `pulso.py` sobre los datos que ya estan en memoria.

    Devuelve dict: `pre`, `post`, `tau_rel`, `forma`, `n_usados`.
    """
    tau_rel, forma, n_us = pulso.pulso_promedio(T, V)
    pre_zona = tau_rel < -8.0
    ruido = float(1.4826 * np.median(np.abs(
        forma[pre_zona] - np.median(forma[pre_zona]))))
    t_ini, t_fin, _ = pulso.ventana_por_convergencia(
        tau_rel, forma, 0.0, 1.0, ruido, k=3.0, tramo=4.0, ancho_suave=1.0)
    return dict(pre=float(max(-t_ini, 2.0)), post=float(t_fin),
               tau_rel=tau_rel, forma=forma, n_usados=n_us)


CAMPOS_AREA = ("area", "amplitud", "t_onset", "t_ini", "t_fin", "base", "ruido")

# de CAMPOS_AREA + lo que agrega ancho.py, cuales son VOLTAJE (mV o mV*ns) y
# por lo tanto escalan con el atenuador; t_onset/t_ini/t_fin/fwhm/t_cruce_* son
# tiempo y el atenuador no los toca.
_CAMPOS_VOLTAJE = ("area", "amplitud", "base", "ruido", "altura_media")


def factor_atenuacion(atenuador_db):
    """
    Factor por el que multiplicar un voltaje MEDIDO (a la entrada del
    digitalizador) para recuperar el voltaje a la SALIDA del detector, antes
    del atenuador en linea.

    Un atenuador en dB actua sobre voltaje (no potencia): 20 dB = factor 10,
    no 100. `analisis/rutas.py` documenta el mismo numero (20 dB divide D por
    10) pero para el pipeline de calibracion; aca se aplica la correccion en
    vez de dejarla, porque el proposito de `pulsos/` es mirar la intensidad
    REAL de cada pulso, no un observable de calibracion.

    `atenuador_db=None` o `0` -> factor 1.0 (sin correccion).
    """
    if not atenuador_db:
        return 1.0
    return 10.0 ** (float(atenuador_db) / 20.0)


def calcular_lote(T, V, pre, post, callback=None, cada=50, atenuador_db=None):
    """
    Area (metodo GLOBAL, `pulso.area_evento(t, v, pre, post)`, el vigente en
    todo el proyecto) y FWHM (`ancho.ancho_media_altura`) de cada pulso
    cargado.

    `atenuador_db`: si se da, todos los campos en VOLTAJE (`area`, `amplitud`,
    `base`, `ruido`, `altura_media`) se multiplican por `factor_atenuacion`
    ANTES de guardarse, para reportar el voltaje real a la salida del
    detector. El FWHM y los tiempos (`t_onset`, `t_ini`, `t_fin`,
    `t_cruce_sube`, `t_cruce_baja`) no cambian: un atenuador no distorsiona el
    tiempo. `pre`/`post` tampoco: son el intervalo de integracion, en ns.

    Devuelve un dict de arrays NumPy, uno por campo de `area_evento` mas
    `fwhm`, `t_cruce_sube`, `t_cruce_baja`, `altura_media`, `bandera_area`,
    `bandera_fwhm`, `pre`, `post`, `atenuador_db`, `factor_atenuacion`.
    """
    factor = factor_atenuacion(atenuador_db)
    acum = {c: [] for c in CAMPOS_AREA}
    acum.update(fwhm=[], t_cruce_sube=[], t_cruce_baja=[], altura_media=[])
    bandera_area, bandera_fwhm = [], []
    n = len(T)
    for i, (t, v) in enumerate(zip(T, V)):
        ra = pulso.area_evento(t, v, pre=pre, post=post)
        rf = ancho.ancho_media_altura(t, v)
        for c in CAMPOS_AREA:
            acum[c].append(ra[c])
        acum["fwhm"].append(rf["fwhm"])
        acum["t_cruce_sube"].append(rf["t_cruce_sube"])
        acum["t_cruce_baja"].append(rf["t_cruce_baja"])
        acum["altura_media"].append(rf["altura_media"])
        bandera_area.append(ra["bandera"])
        bandera_fwhm.append(rf["bandera"])
        if callback is not None and (i + 1) % cada == 0:
            callback((i + 1) / n, i + 1)
    out = {c: np.asarray(acum[c], float) for c in acum}
    if factor != 1.0:
        for c in _CAMPOS_VOLTAJE:
            out[c] = out[c] * factor
    out["bandera_area"] = np.asarray(bandera_area, np.int16)
    out["bandera_fwhm"] = np.asarray(bandera_fwhm, np.int16)
    out["pre"], out["post"] = float(pre), float(post)
    out["atenuador_db"] = float(atenuador_db) if atenuador_db else 0.0
    out["factor_atenuacion"] = factor
    if callback is not None:
        callback(1.0, n)
    return out


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    claves = [k for k in rutas.CORRIDAS if rutas.corrida(k).disponible()]
    if not claves:
        print("No hay corridas accesibles. Corre `python analisis/rutas.py`.")
        return
    clave = claves[0]
    c = rutas.corrida(clave)
    print("=" * 74)
    print("%s  —  %s  (%d pulsos)" % (clave, c.descripcion, n))
    print("=" * 74)

    T, V, _, _, canal = cargar_pulsos(c.muestra, n=n)
    print("  leidos %d pulsos por el canal %s" % (len(T), canal))

    vent = sugerir_ventana(T, V)
    print("  ventana sugerida: pre = %.1f ns, post = %.1f ns (%d pulsos usados "
          "para la forma promedio)" % (vent["pre"], vent["post"], vent["n_usados"]))

    d = calcular_lote(T, V, vent["pre"], vent["post"])
    ok = (d["bandera_area"] == 0) & (d["bandera_fwhm"] == 0)
    print("  %d/%d pulsos sin ninguna bandera" % (int(ok.sum()), len(T)))
    print("  area  : mediana %.1f mV*ns (std %.1f)" % (np.median(d["area"][ok]), d["area"][ok].std()))
    print("  FWHM  : mediana %.2f ns (std %.2f)" % (np.median(d["fwhm"][ok]), d["fwhm"][ok].std()))
    print("  amp   : mediana %.1f mV (std %.1f)" % (np.median(d["amplitud"][ok]), d["amplitud"][ok].std()))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    ax = axes[0]
    for i in range(min(60, len(T))):
        tau, vv = pulso.preparar(T[i], V[i])
        ax.plot(tau - d["t_onset"][i], vv - d["base"][i], lw=0.5, alpha=0.35)
    ax.axvspan(-vent["pre"], vent["post"], color="C2", alpha=0.15,
              label="ventana aplicada")
    ax.set_title("%s: %d pulsos superpuestos (alineados por onset)"
                % (clave, min(60, len(T))))
    ax.set_xlabel("t - t_onset [ns]"); ax.set_ylabel("v - base [mV]")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.hist(d["fwhm"][ok], bins=40, color="C0")
    ax.set_title("distribucion del FWHM")
    ax.set_xlabel("FWHM [ns]"); ax.set_ylabel("cuentas"); ax.grid(alpha=0.3)

    ax = axes[2]
    ax.hist(d["area"][ok], bins=40, color="C1")
    ax.set_title("distribucion del area")
    ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas"); ax.grid(alpha=0.3)

    fig.tight_layout()
    salida = rutas.figura("01_lote.png")
    fig.savefig(salida, dpi=125)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
