"""
Comprobar la Etapa 2: contra el oraculo y contra pulsos de integral conocida.

    python desde_cero/verificar_area.py

1. ORACULO. `analisis/area.py:areas_de_archivo(metodo="global")` con el MISMO
   pre/post (2 y 12 ns), sin cache, sobre 5000 eventos de cada corrida de
   trabajo. El oraculo TIRA los eventos de pico a pico <= 3 mV; nosotros los
   marcamos, asi que se alinean quitando esos de nuestra lista.

2. VERDAD CONOCIDA. Pulsos fabricados con la forma medida del pulso (la de la
   Etapa 1), sobre un eje temporal real del DRS4. El area verdadera es
   A * integral de la forma en [t0 - pre, t0 + post], con t0 el onset
   VERDADERO: se conoce sin medir nada. Tres pruebas:
     a. sin ruido: lo unico que queda es que los bordes caen sobre muestras;
     b. con ruido de 1 mV y amplitudes de 5 a 300 mV: sesgo del area contra
        amplitud (el onset tarde de la Etapa 1 a bajo SNR) y, aparte, sesgo
        del COCIENTE area/amplitud, que es el control que se usa en datos;
     c. registros truncados: la bandera tiene que prender siempre y el area
        salir por debajo de la verdadera en lo que predice la forma.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import area as oraculo                        # noqa: E402
from analisis import pulso as oraculo_pulso                 # noqa: E402
from desde_cero import area, lector, pulso                  # noqa: E402
from desde_cero.verificar_pulso import PLANTILLA, fabricar  # noqa: E402

N_EVENTOS = 5000


def contra_oraculo(ruta, n=N_EVENTOS):
    nos = area.areas_de_archivo(ruta, limite=n)
    ora = oraculo.areas_de_archivo(ruta, metodo="global", pre=area.PRE,
                                   post=area.POST, limite=n, usar_cache=False,
                                   verbose=False, recorte=pulso.RECORTE)
    m = nos["pico_a_pico"] > area.UMBRAL_PP
    a1, a2 = nos["area"][m], ora["area"]
    if a1.size != a2.size:
        return dict(ok=False, motivo="%d eventos contra %d" % (a1.size, a2.size))
    rel = np.abs(a1 - a2) / np.maximum(np.abs(a2), 1e-9)
    b1, b2 = nos["bandera"][m], ora["bandera"]
    # Los eventos SIN_PRE se comparan aparte: ahi la linea base sale de la
    # COLA, y la cola empieza en lugares distintos a proposito (prueba d).
    sp = ((b1 | b2) & area.SIN_PRE) > 0
    return dict(ok=True, n=a1.size, tirados=int((~m).sum()),
                rel_mediana=float(np.median(rel[~sp])),
                rel_max=float(rel[~sp].max()),
                banderas=int((b1[~sp] != b2[~sp]).sum()),
                n_sin_pre=int(sp.sum()),
                sin_pre_mediana=float(np.median(((a1 - a2) / np.abs(a2))[sp]))
                if sp.any() else 0.0)


# --------------------------------------------------------------------------
def _forma():
    z = np.load(PLANTILLA)
    return z["tau"], z["forma"]


def integral_verdadera(ftau, f, amp, pre=area.PRE, post=area.POST):
    """A * integral de la forma entre -pre y +post (onset verdadero en 0)."""
    g = np.linspace(-pre, post, 20001)
    return amp * np.trapezoid(np.interp(g, ftau, f), g)


def fabricados(eje_t, sigma, amps, n, rng, t0_rango=(40.0, 70.0)):
    """
    Para cada amplitud: media y dispersion de area/verdad - 1, y media de
    (area/amplitud_medida) / (verdad/amplitud_verdadera) - 1.
    """
    ftau, f = _forma()
    t = np.asarray(eje_t, float)
    filas = []
    for amp in amps:
        verdad = integral_verdadera(ftau, f, amp)
        e_area, e_ctrl, band = [], [], 0
        for _ in range(n):
            t0 = rng.uniform(*t0_rango)
            s = rng.choice((-1.0, 1.0))
            v = fabricar(t, ftau, f, amp, t0, rng.uniform(-3, 3), sigma, s, rng) \
                if sigma > 0 else \
                rng.uniform(-3, 3) + s * amp * np.interp(t - t0, ftau, f,
                                                         left=0.0, right=0.0)
            r = area.area_evento(t, v)
            e_area.append(r["area"] / verdad - 1.0)
            e_ctrl.append((r["area"] / r["amplitud"]) / (verdad / amp) - 1.0)
            band += r["bandera"] != 0
        e_area, e_ctrl = np.asarray(e_area), np.asarray(e_ctrl)
        filas.append(dict(amp=amp, sesgo=e_area.mean(),
                          err=e_area.std(ddof=1) / np.sqrt(n),
                          disp=e_area.std(ddof=1), ctrl=e_ctrl.mean(),
                          banderas=band))
    return filas


def truncados(eje_t, n, rng, amp=100.0, sigma=1.0):
    """
    Pulsos con el onset tan tarde que el registro termina entre 3 y 11 ns
    despues del onset. Devuelve (fraccion con bandera, error medio del area,
    error medio PREDICHO por la forma: la carga que queda despues del fin).
    """
    ftau, f = _forma()
    t = np.asarray(eje_t, float)
    fin = t[-1 - pulso.RECORTE]
    band, err, pred = 0, [], []
    for _ in range(n):
        queda = rng.uniform(3.0, 11.0)
        t0 = fin - queda
        v = fabricar(t, ftau, f, amp, t0, 0.0, sigma, -1.0, rng)
        r = area.area_evento(t, v)
        band += bool(r["bandera"] & area.TRUNCADO)
        verdad = integral_verdadera(ftau, f, amp)
        err.append(r["area"] / verdad - 1.0)
        pred.append(integral_verdadera(ftau, f, amp, post=queda) / verdad - 1.0)
    return band / n, float(np.mean(err)), float(np.mean(pred))


def cola_de_respaldo(eje_t, n, rng, sigma=1.0):
    """
    Pulsos con el onset 1-5 ns despues del inicio del registro: no hay
    pre-pulso y la base sale de la cola. Compara nuestra cola (desde onset +
    pulso.COLA + guarda = +84 ns) con la del oraculo (desde onset + post +
    guarda = +16 ns con la compuerta de 12 ns).

    La del oraculo depende de la COMPUERTA, y no deberia: donde termina el
    pulso es una propiedad del pulso, no de cuanto se integra. A +16 ns el
    pulso todavia vale ~5 % de su amplitud y tiene una cola de ~1.5 % que
    sigue viva a 130 ns (§7.66), asi que la mediana de la cola sale corrida
    hacia el pulso y el area, por debajo.
    """
    ftau, f = _forma()
    t = np.asarray(eje_t, float)
    filas = []
    for amp in (30.0, 100.0):
        verdad = integral_verdadera(ftau, f, amp)
        en, eo = [], []
        for _ in range(n):
            t0 = t[pulso.RECORTE] + rng.uniform(1.0, 5.0)
            v = fabricar(t, ftau, f, amp, t0, rng.uniform(-3, 3), sigma, -1.0, rng)
            en.append(area.area_evento(t, v)["area"] / verdad - 1.0)
            eo.append(oraculo_pulso.area_evento(t, v, pre=area.PRE,
                                                post=area.POST,
                                                recorte=pulso.RECORTE)["area"]
                      / verdad - 1.0)
        filas.append((amp, float(np.mean(en)), float(np.mean(eo))))
    return filas


def main():
    ok = True
    print("=" * 78)
    print("Etapa 2 contra analisis/area.py  (pre = %.0f ns, post = %.0f ns, "
          "%d eventos)" % (area.PRE, area.POST, N_EVENTOS))
    print("=" * 78)
    eje = None
    for clave, (ruta, _) in lector.CORRIDAS_TRABAJO.items():
        if not os.path.exists(ruta):
            print("  %-12s no esta en disco" % clave)
            continue
        r = contra_oraculo(ruta)
        if not r["ok"]:
            ok = False
            print("  %-12s NO ALINEAN: %s" % (clave, r["motivo"]))
            continue
        bien = r["rel_max"] < 1e-6 and r["banderas"] == 0
        ok &= bien
        print("  %-12s %5d ev (%d sin pulso, que el oraculo tira)  "
              "|d area|/area: mediana %.1e, max %.1e   banderas distintas %d   %s"
              % (clave, r["n"] - r["n_sin_pre"], r["tirados"], r["rel_mediana"],
                 r["rel_max"], r["banderas"], "ok" if bien else "DIFIEREN"))
        print("  %-12s %5d ev SIN_PRE (base de la cola): nuestra area / la del "
              "oraculo - 1 = %+.2f %% (mediana)"
              % ("", r["n_sin_pre"], 100 * r["sin_pre_mediana"]))
        if lector.formato(ruta) == "binario":
            eje = next(lector.eventos(ruta, limite=1))[0]
    if eje is None:
        eje = np.arange(1024) * 0.2

    rng = np.random.default_rng(2)
    print()
    print("=" * 78)
    print("Pulsos fabricados: area contra la integral verdadera en la compuerta")
    print("=" * 78)
    print("  a. SIN ruido (solo el efecto de que los bordes caen en muestras)")
    for f in fabricados(eje, 0.0, (10.0, 100.0), 2000, rng):
        bien = abs(f["sesgo"]) < 2e-3 and f["disp"] < 5e-3
        ok &= bien
        print("     A = %5.0f mV   sesgo %+.4f %%   dispersion %.3f %%   %s"
              % (f["amp"], 100 * f["sesgo"], 100 * f["disp"],
                 "ok" if bien else "FALLA"))

    print("  b. ruido 1 mV")
    print("     %7s  %22s  %12s  %22s"
          % ("A [mV]", "sesgo del area", "dispersion",
             "sesgo de area/amplitud"))
    for f in fabricados(eje, 1.0, (5.0, 10.0, 30.0, 100.0, 300.0), 3000, rng):
        exigible = f["amp"] >= 30.0          # SNR >= 30, como en la Etapa 1
        bien = abs(f["sesgo"]) < 0.005
        if exigible:
            ok &= bien
        print("     %7.0f  %+7.3f %% +- %.3f %%  %10.2f %%  %+10.2f %%   %s"
              % (f["amp"], 100 * f["sesgo"], 100 * f["err"], 100 * f["disp"],
                 100 * f["ctrl"],
                 ("ok" if bien else "FALLA") if exigible else "(SNR < 30)"))

    print("  c. registro que termina 3-11 ns despues del onset")
    fr, e, p = truncados(eje, 2000, rng)
    bien = fr == 1.0 and abs(e - p) < 0.01
    ok &= bien
    print("     bandera TRUNCADO en el %.1f %%   area %+.1f %% contra %+.1f %% "
          "predicho por la forma   %s"
          % (100 * fr, 100 * e, 100 * p, "ok" if bien else "FALLA"))
    print("  d. sin pre-pulso: la linea base sale de la cola")
    for amp, en, eo in cola_de_respaldo(eje, 2000, rng):
        bien = abs(en) < 0.005
        ok &= bien
        print("     A = %4.0f mV   nuestra cola (+84 ns) %+.2f %%   la del "
              "oraculo (+16 ns) %+.2f %%   %s"
              % (amp, 100 * en, 100 * eo, "ok" if bien else "FALLA"))
    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
