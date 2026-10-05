"""
ancho.py — el ancho a media altura (FWHM) TEMPORAL de un pulso individual.

No existe en ningun otro lugar del repo. `codigo/etapas/model.py::punto_mitad_altura`
es para el borde Compton del ESPECTRO (canal vs cuentas, sin interpolacion); esto
es sobre UN pulso crudo (t vs v, en ns y mV).

Reusa literal de `analisis/pulso.py` la localizacion de pico/onset/linea base
(`preparar`, `detectar_pico`, `detectar_onset`, `linea_base`): es la MISMA
receta que usa `area_evento`, para que "donde esta el pulso" signifique lo
mismo en toda la herramienta. Lo unico nuevo es la busqueda de los dos cruces
a media altura, por interpolacion LINEAL sub-muestra entre las dos muestras
que encierran el cruce.

Los cruces se buscan igual que el onset: RETROCEDIENDO desde el pico para la
subida, AVANZANDO desde el pico para la bajada. Nunca desde los bordes del
registro, o una fluctuacion de linea base lejana contaria como "medio pulso".

Correr `python pulsos/ancho.py` deja `figuras/00_fwhm.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import drs_xml, pulso
from pulsos import rutas

# Banderas propias: continuan la numeracion de bits de pulso.py (llega a 16).
# Las de pulso.py (MAL_SATURA, MAL_SIN_PRE, MAL_SIN_PICO, MAL_TRUNCADO) se
# reusan por import, no se redefinen.
MAL_SIN_CRUCE_SUBIDA = 32   # el flanco de subida no baja de la mitad en el registro
MAL_SIN_CRUCE_BAJADA = 64   # el flanco de bajada no baja de la mitad (pulso truncado)


def _cruce_lineal(tau, d, i0, i1, nivel):
    """Cruce de `d` por `nivel` entre las muestras i0 < i1, por interpolacion lineal."""
    d0, d1 = d[i0], d[i1]
    if d1 == d0:
        return float(tau[i0])
    return float(tau[i0] + (nivel - d0) / (d1 - d0) * (tau[i1] - tau[i0]))


def ancho_media_altura(t, v, guarda_base=4.0, estadistico="mediana",
                       frac_onset=pulso.FRAC_ONSET, v_satura=pulso.V_SATURA,
                       recorte=pulso.RECORTE_BORDES):
    """
    FWHM de un pulso: tiempo entre el cruce de SUBIDA y el de BAJADA por la
    mitad de la amplitud, con interpolacion lineal entre muestras.

    Devuelve un dict: `fwhm` [ns], `t_cruce_sube`, `t_cruce_baja` (en `tau`,
    la MISMA convencion que `area_evento`: relativas al inicio del registro,
    no al onset), `t_onset`, `t_pico`, `amplitud`, `altura_media`, `base`,
    `ruido`, `bandera`.

    Si un cruce no aparece dentro del registro (pulso truncado, o el pico cae
    casi sobre el arranque) se MARCA la bandera correspondiente y se usa el
    extremo del registro como reemplazo — igual filosofia que `area_evento`:
    marcar, no descartar; descartar es decision de quien llama.
    """
    tau, v = pulso.preparar(t, v, recorte)
    if tau.size < 16:
        return dict(fwhm=np.nan, t_cruce_sube=np.nan, t_cruce_baja=np.nan,
                    t_onset=np.nan, t_pico=np.nan, amplitud=0.0,
                    altura_media=np.nan, base=np.nan, ruido=np.nan,
                    bandera=pulso.MAL_TRUNCADO)

    # --- localizar el pulso con una referencia provisoria, igual que area_evento ---
    base0 = float(np.median(v))
    i_p, t_pico, amp0, signo = pulso.detectar_pico(tau, v, base0)
    i_on, t_on = pulso.detectar_onset(tau, v, base0, i_p, signo, frac_onset)

    bandera = 0
    if np.max(np.abs(v)) >= v_satura:
        bandera |= pulso.MAL_SATURA

    base, ruido, n_base = pulso.linea_base(tau, v, hasta=t_on - guarda_base,
                                           estadistico=estadistico)
    if n_base < 10:
        bandera |= pulso.MAL_SIN_PRE
        base = base0            # sin pre-pulso util, se cae a la mediana global

    i_p, t_pico, amp, signo = pulso.detectar_pico(tau, v, base)
    if amp < 5.0 * ruido:
        bandera |= pulso.MAL_SIN_PICO
    i_on, t_on = pulso.detectar_onset(tau, v, base, i_p, signo, frac_onset)

    d = signo * (v - base)
    nivel = amp / 2.0

    idx_sub = np.flatnonzero(d[:i_p + 1] < nivel)
    if idx_sub.size:
        i0 = int(idx_sub[-1])
        i1 = min(i0 + 1, tau.size - 1)
        t_sube = _cruce_lineal(tau, d, i0, i1, nivel)
    else:
        bandera |= MAL_SIN_CRUCE_SUBIDA
        t_sube = float(tau[0])

    idx_baj = np.flatnonzero(d[i_p:] < nivel)
    if idx_baj.size:
        i1 = int(i_p + idx_baj[0])
        i0 = max(i1 - 1, 0)
        t_baja = _cruce_lineal(tau, d, i0, i1, nivel)
    else:
        bandera |= MAL_SIN_CRUCE_BAJADA
        t_baja = float(tau[-1])

    return dict(fwhm=float(t_baja - t_sube), t_cruce_sube=t_sube,
               t_cruce_baja=t_baja, t_onset=float(t_on), t_pico=float(t_pico),
               amplitud=float(amp), altura_media=float(nivel),
               base=float(base), ruido=float(ruido), bandera=int(bandera))


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.CORRIDAS if rutas.corrida(k).disponible()]
    if not claves:
        print("No hay corridas accesibles. Corre `python analisis/rutas.py`.")
        return
    clave = claves[0]
    c = rutas.corrida(clave)
    print("=" * 74)
    print("%s  —  %s" % (clave, c.descripcion))
    print("=" * 74)

    n_leer = 400
    T, V, _, _ = drs_xml.leer_pulsos(c.muestra, n=n_leer)

    fwhm, amp, banderas = [], [], []
    resultados = []
    for t, v in zip(T, V):
        r = ancho_media_altura(t, v)
        resultados.append(r)
        fwhm.append(r["fwhm"])
        amp.append(r["amplitud"])
        banderas.append(r["bandera"])
    fwhm = np.array(fwhm)
    amp = np.array(amp)
    banderas = np.array(banderas, dtype=np.int16)

    ok = banderas == 0
    print("  %d pulsos, %d sin ninguna bandera" % (len(fwhm), int(ok.sum())))
    print("  FWHM (solo los sin bandera): mediana %.2f ns, desviacion %.2f ns"
          % (np.median(fwhm[ok]), fwhm[ok].std()))
    for nombre, bit in (("SATURA", pulso.MAL_SATURA),
                       ("TRUNCADO", pulso.MAL_TRUNCADO),
                       ("SIN_PRE", pulso.MAL_SIN_PRE),
                       ("SIN_PICO", pulso.MAL_SIN_PICO),
                       ("SIN_CRUCE_SUBIDA", MAL_SIN_CRUCE_SUBIDA),
                       ("SIN_CRUCE_BAJADA", MAL_SIN_CRUCE_BAJADA)):
        n = int((banderas & bit).astype(bool).sum())
        if n:
            print("     %-18s %d" % (nombre, n))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    ax = axes[0]
    ax.hist(fwhm[ok], bins=40, color="C0")
    ax.set_title("%s: distribucion del FWHM (%d pulsos limpios)" % (clave, int(ok.sum())))
    ax.set_xlabel("FWHM [ns]"); ax.set_ylabel("cuentas"); ax.grid(alpha=0.3)

    ax = axes[1]
    idx_buenos = np.flatnonzero(ok)[:12]
    for i in idx_buenos:
        tau, vv = pulso.preparar(T[i], V[i])
        r = resultados[i]
        signo = 1.0 if (vv[np.argmax(np.abs(vv - r["base"]))] - r["base"]) >= 0 else -1.0
        ax.plot(tau, signo * (vv - r["base"]), lw=0.8, alpha=0.6)
        ax.plot([r["t_cruce_sube"], r["t_cruce_baja"]],
               [r["altura_media"], r["altura_media"]], "o", color="k",
               ms=3, alpha=0.6)
    ax.set_title("cruces de media altura (12 pulsos)")
    ax.set_xlabel("t [ns]"); ax.set_ylabel("v - base, orientado [mV]")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    salida = rutas.figura("00_fwhm.png")
    fig.savefig(salida, dpi=125)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
