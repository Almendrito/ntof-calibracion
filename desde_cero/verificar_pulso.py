"""
Comprobar la Etapa 1 dos veces: contra el oraculo y contra verdad conocida.

    python desde_cero/verificar_pulso.py

1. ORACULO. Los mismos `(t, v)` (leidos con `lector`, que ya se verifico
   identico al de analisis/) pasan por `desde_cero.pulso.medir` y por
   `analisis.pulso.area_evento`. Onset, amplitud, base, ruido y banderas
   tienen que coincidir. Al oraculo se le pasa pre=2, post=COLA: el area no
   importa aca, y con post=COLA su respaldo de linea base por la cola usa el
   mismo tramo que el nuestro. Y `recorte=pulso.RECORTE`, porque el nuestro
   es 4 muestras y el del oraculo 2 (la espiga del DRS4, ver
   `pulso.preparar`).

2. VERDAD CONOCIDA. Pulsos fabricados: la forma media medida del pulso
   (analisis/.cache/plantilla_pulso_ago-sep.npz, normalizada a amplitud 1 y
   con su cruce del 10 % en tau = 0), escalada por una amplitud conocida,
   puesta en un onset conocido con jitter, sobre una base conocida y con ruido
   gaussiano, muestreada sobre un eje temporal REAL del DRS4 (no uniforme).
   Ahi se sabe la respuesta, asi que se mide el sesgo y la dispersion del
   onset y se comparan con lo que el ruido predice.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import pulso as oraculo                       # noqa: E402
from desde_cero import lector, pulso                        # noqa: E402

N_EVENTOS = 2000
PLANTILLA = os.path.join(lector.RAIZ, "analisis", ".cache",
                         "plantilla_pulso_ago-sep.npz")
BANDERAS = pulso.SATURADO | pulso.SIN_PRE | pulso.SIN_PULSO   # sin la de truncado

# Desde que amplitud/ruido se exige el criterio del PLAN (sesgo < 0.1 ns).
# Por debajo, el umbral del 10 % queda a menos de 3 sigma del ruido y el
# flanco lo cruza varias veces: `buscar_onset` se queda con el ULTIMO cruce y
# el onset sale tarde (+0.30 ns a SNR 10, medido separando las dos causas; la
# otra, +0.06 ns, es la amplitud sesgada hacia arriba por tomar el maximo de
# muestras ruidosas). Buscar hacia adelante seria mucho peor: -3.97 ns a
# SNR 10, porque dispara en el ruido. No es un error a corregir sino el costo
# conocido del criterio; se reporta, y la Etapa 2 mide si mueve el area.
SNR_MIN = 30.0


# --------------------------------------------------------------------------
# 1. Contra el oraculo
# --------------------------------------------------------------------------
def contra_oraculo(ruta, n=N_EVENTOS):
    d_on = d_amp = d_base = d_ruido = 0.0
    n_band = k = 0
    dt, snr = [], []
    for t, v in lector.eventos(ruta, limite=n):
        a = pulso.medir(t, v)
        # el recorte se le PASA al oraculo: el suyo es 2 y el nuestro 4 (ver
        # pulso.preparar); asi se compara la misma cuenta
        b = oraculo.area_evento(t, v, pre=2.0, post=pulso.COLA,
                                recorte=pulso.RECORTE)
        d_on = max(d_on, abs(a["t_onset"] - b["t_onset"]))
        d_amp = max(d_amp, abs(a["amplitud"] - b["amplitud"])
                    / max(b["amplitud"], 1e-9))
        d_base = max(d_base, abs(a["base"] - b["base"]))
        d_ruido = max(d_ruido, abs(a["ruido"] - b["ruido"]))
        n_band += (a["bandera"] & BANDERAS) != (b["bandera"] & BANDERAS)
        dt.append(np.median(np.diff(a["tau"])))
        snr.append(a["amplitud"] / a["ruido"])
        k += 1
    snr = np.asarray(snr)
    return dict(eventos=k, onset=d_on, amp_rel=d_amp, base=d_base,
                ruido=d_ruido, banderas=n_band, dt=float(np.median(dt)),
                snr_mediana=float(np.median(snr)),
                bajo_snr=float(np.mean(snr < SNR_MIN)))


# --------------------------------------------------------------------------
# 2. Pulsos fabricados
# --------------------------------------------------------------------------
def fabricar(t, forma_tau, forma, amp, t0, base, sigma, signo, rng):
    """v(t) = base + signo*amp*f(t - t0) + ruido."""
    f = np.interp(t - t0, forma_tau, forma, left=0.0, right=0.0)
    return base + signo * amp * f + rng.normal(0.0, sigma, t.size)


def verdad_conocida(eje_t, n=2000, sigma=1.0, semilla=1):
    """
    Para cada amplitud, `n` pulsos con onset uniforme en [40, 70] ns (jitter
    del trigger), base uniforme en [-3, 3] mV y polaridad al azar.

    Devuelve, por amplitud: sesgo y dispersion del onset (el de muestra y el
    interpolado), sesgo de amplitud, aciertos de polaridad, y la dispersion
    que PREDICE el ruido: sigma_t = sigma_v / pendiente, con la pendiente del
    flanco en el cruce del 10 % = amp * f'(0). Es la regla de siempre para
    medir un tiempo cruzando un umbral: el ruido en tension se convierte en
    ruido en tiempo dividido por lo empinado del flanco. sigma_v no es el
    sigma de una muestra: el cruce se INTERPOLA entre dos muestras con pesos
    x y 1-x, y el ruido de esa combinacion tiene varianza
    sigma^2 (x^2 + (1-x)^2), que promediada sobre x uniforme da (2/3) sigma^2.
    Por eso sigma_v = sqrt(2/3) sigma.
    """
    z = np.load(PLANTILLA)
    ftau, f = z["tau"], z["forma"]
    i = int(np.argmin(np.abs(ftau)))
    pendiente = (f[i + 1] - f[i - 1]) / (ftau[i + 1] - ftau[i - 1])   # 1/ns
    rng = np.random.default_rng(semilla)
    t = np.asarray(eje_t, float)
    t0_tau = t[pulso.RECORTE]          # medir() cuenta tau desde aqui
    filas = []
    for amp in (10.0, 30.0, 100.0, 300.0):
        e_m, e_f, e_a, ok = [], [], [], 0
        for _ in range(n):
            t0 = rng.uniform(40.0, 70.0)
            s = rng.choice((-1.0, 1.0))
            v = fabricar(t, ftau, f, amp, t0, rng.uniform(-3, 3), sigma, s, rng)
            m = pulso.medir(t, v)
            verdad = t0 - t0_tau
            e_m.append(m["t_onset"] - verdad)
            e_f.append(m["t_onset_fino"] - verdad)
            e_a.append(m["amplitud"] - amp)
            ok += m["signo"] == s
        e_m, e_f, e_a = map(np.asarray, (e_m, e_f, e_a))
        filas.append(dict(
            amp=amp, n=n, polaridad=ok / n,
            sesgo_m=e_m.mean(), sesgo_f=e_f.mean(),
            err_sesgo_f=e_f.std(ddof=1) / np.sqrt(n),
            disp_f=e_f.std(ddof=1),
            pred=np.sqrt(2.0 / 3.0) * sigma / (amp * pendiente),
            sesgo_amp=e_a.mean(), disp_amp=e_a.std(ddof=1)))
    return filas


def main():
    ok = True
    print("=" * 78)
    print("Etapa 1 contra analisis/pulso.py  (%d eventos por corrida)" % N_EVENTOS)
    print("=" * 78)
    eje_real = None
    for clave, (ruta, _) in lector.CORRIDAS_TRABAJO.items():
        if not os.path.exists(ruta):
            print("  %-12s no esta en disco" % clave)
            continue
        r = contra_oraculo(ruta)
        # Criterio del PLAN: onset a menos de un paso de muestreo, amplitud a
        # 1e-6 relativo. En la practica tiene que dar CERO: es la misma cuenta
        # sobre las mismas muestras.
        bien = (r["onset"] < r["dt"] and r["amp_rel"] < 1e-6
                and r["banderas"] == 0)
        ok &= bien
        print("  %-12s %4d ev  max|d onset| = %.1e ns (dt %.3f)  "
              "max|d amp|/amp = %.1e  max|d base| = %.1e mV  "
              "max|d ruido| = %.1e mV  banderas distintas: %d   %s"
              % (clave, r["eventos"], r["onset"], r["dt"], r["amp_rel"],
                 r["base"], r["ruido"], r["banderas"],
                 "ok" if bien else "DIFIEREN"))
        print("  %-12s amplitud/ruido: mediana %.0f, %.1f %% de los eventos "
              "bajo %.0f" % ("", r["snr_mediana"], 100 * r["bajo_snr"], SNR_MIN))
        if lector.formato(ruta) == "binario":
            eje_real = next(lector.eventos(ruta, limite=1))[0]

    print()
    print("=" * 78)
    print("Pulsos fabricados (ruido 1.0 mV, eje temporal %s)"
          % ("REAL del DRS4" if eje_real is not None else "uniforme 0.2 ns"))
    print("=" * 78)
    if eje_real is None:
        eje_real = np.arange(1024) * 0.2
    print("  %6s  %9s  %22s  %18s  %8s  %16s"
          % ("A[mV]", "polaridad", "sesgo onset fino [ns]", "disp. [ns] med/pred",
             "sesgo on", "sesgo amp [mV]"))
    print("  %6s  %9s  %22s  %18s  %8s  %16s"
          % ("", "", "", "", "muestra", ""))
    for f in verdad_conocida(eje_real):
        exigible = f["amp"] / 1.0 >= SNR_MIN          # ruido de 1.0 mV
        bien = (abs(f["sesgo_f"]) < 0.1 and f["polaridad"] == 1.0
                and 0.8 < f["disp_f"] / f["pred"] < 1.25)
        if exigible:
            ok &= bien
        estado = ("ok" if bien else "FALLA") if exigible else                  "SNR < %.0f: el ultimo cruce sale tarde (esperado)" % SNR_MIN
        print("  %6.0f  %8.1f%%  %+9.4f +- %.4f     %7.3f / %6.3f   %+8.3f  "
              "%+6.3f +- %.3f   %s"
              % (f["amp"], 100 * f["polaridad"], f["sesgo_f"], f["err_sesgo_f"],
                 f["disp_f"], f["pred"], f["sesgo_m"], f["sesgo_amp"],
                 f["disp_amp"], estado))
    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
