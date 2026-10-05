"""
Que los dos lectores del formato BINARIO del DRS4 den lo MISMO.

Por que existe
--------------
Desde la sesion del 1-sep-2026 hay archivos en el formato binario del DRS4
("Save as binary" del DRSOSC) y sin extension en el nombre. Los lee
`analisis/drs_bin.py`, y `codigo/datos_io/spectrum.py` los lee DELEGANDO en ese
mismo modulo en vez de tener su propio parser (un segundo parser del mismo
formato es lo que se desincroniza en silencio, CLAUDE.md 7.36).

"Delega" es una afirmacion que hay que poder verificar, y ademas el disco
externo no siempre esta conectado. Por eso este archivo FABRICA un binario del
DRS4 con verdad conocida —pulsos de amplitud y posicion elegidas por nosotros— y
comprueba cuatro cosas:

  1. la deteccion de formato NO se guia por el nombre: el archivo de prueba no
     tiene extension, igual que los de la sesion;
  2. `drs_bin` reconstruye la tension y el eje temporal que se escribieron;
  3. `spectrum` (la cadena de `codigo/`) lee el mismo archivo y saca las mismas
     areas que `area` (la cadena de `analisis/`);
  4. un archivo SIN scaler (versiones viejas del DRSOSC) no rompe nada.

    python verificaciones/_verif_binario.py
"""

import os as _os
import sys as _sys
_RAIZ = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_CODIGO = _os.path.join(_RAIZ, "codigo")
for _p in (_RAIZ, _CODIGO):
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import struct
import tempfile

import numpy as np

from analisis import drs, drs_bin, area as a_area
from datos_io import spectrum


N_MUESTRAS = 1024
RANGO_MV = 0                     # centro de la ventana de tension del DRS4
TAU_SUBIDA = 1.5
TAU_BAJADA = 12.0


# --------------------------------------------------------------------------
# 1. Fabricar un binario del DRS4 con verdad conocida
# --------------------------------------------------------------------------
def _u16_de_mV(v, rango=RANGO_MV):
    """La inversa exacta de la conversion del `read_binary.cpp` del DRS4."""
    u = (v / 1000.0 + 0.5 - rango / 1000.0) * 65536.0
    return np.clip(np.rint(u), 0, 65535).astype("<u2")


def pulso_modelo(t, t0, amplitud, tau_subida=TAU_SUBIDA, tau_bajada=TAU_BAJADA):
    """
    Pulso NEGATIVO tipo centellador: subida rapida, cola exponencial.

    Es la forma que tienen los datos reales (el pulso muere a ~85 ns). Lo que
    importa aca no es que sea fiel sino que su area sea CALCULABLE, para poder
    comparar contra lo que reconstruye la cadena.
    """
    x = np.maximum(t - t0, 0.0)
    f = np.exp(-x / tau_bajada) - np.exp(-x / tau_subida)
    return -amplitud * f / f.max()


def escribir_binario(ruta, n_eventos=400, semilla=0, con_scaler=True,
                     canal=1, serie_placa=2859, scaler_hz=170.0,
                     dt_nominal=0.2, t0=40.0):
    """
    Escribe un binario del DRS4 valido y devuelve la verdad de cada evento.

    El eje temporal se arma con anchos de celda LEVEMENTE distintos entre si (el
    DRS4 real no muestrea uniforme) y con una celda de trigger distinta por
    evento, que es justo lo que obliga a reconstruir el tiempo evento a evento.
    """
    rng = np.random.default_rng(semilla)
    ancho = (dt_nominal * (1.0 + 0.05 * rng.standard_normal(N_MUESTRAS))
             ).astype("<f4")
    amplitudes = rng.uniform(30.0, 300.0, n_eventos)
    celdas = rng.integers(0, N_MUESTRAS, n_eventos)

    with open(ruta, "wb") as f:
        f.write(b"DRS2")
        f.write(b"TIME")
        f.write(b"B#" + struct.pack("<H", serie_placa))
        f.write(("C%03d" % canal).encode())
        f.write(ancho.tobytes())

        for k in range(n_eventos):
            f.write(b"EHDR")
            f.write(struct.pack("<i", k + 1))
            seg, ms = divmod(int(k * 1000 / 70.0), 1000)      # ~70 eventos/s
            mm, ss = divmod(seg, 60)
            f.write(struct.pack("<8H", 2026, 9, 1, 12, 30 + mm, ss, ms,
                                RANGO_MV & 0xFFFF))
            f.write(b"B#" + struct.pack("<H", serie_placa))
            f.write(b"T#" + struct.pack("<H", int(celdas[k])))
            f.write(("C%03d" % canal).encode())
            if con_scaler:
                f.write(struct.pack("<i", int(scaler_hz)))
            t = drs_bin._eje_temporal(ancho.astype(float), celdas[k])
            f.write(_u16_de_mV(pulso_modelo(t, t0, amplitudes[k])).tobytes())

    return dict(ancho=ancho.astype(float), amplitudes=amplitudes,
                celdas=celdas, n_eventos=n_eventos, t0=t0, canal=canal)


def area_analitica(amplitud, pre, post, tau_subida=TAU_SUBIDA,
                   tau_bajada=TAU_BAJADA):
    """
    Area exacta del pulso modelo dentro de la compuerta, en mV*ns.

    Sirve para que la comparacion no sea solo "los dos lectores coinciden": dos
    lectores pueden coincidir y estar los dos equivocados.
    """
    t = np.linspace(0.0, float(post), 200001)
    f = np.exp(-t / tau_bajada) - np.exp(-t / tau_subida)
    return float(amplitud * abs(np.trapezoid(f / f.max(), t)))


# --------------------------------------------------------------------------
# 2. Las comprobaciones
# --------------------------------------------------------------------------
def _titulo(t):
    print()
    print("-" * 74)
    print(t)
    print("-" * 74)


def main():
    fallos = []
    tmp = tempfile.mkdtemp(prefix="drs_verif_")
    # SIN EXTENSION a proposito: es como se llaman los archivos de la sesion
    ruta = _os.path.join(tmp, "Cs137")
    verdad = escribir_binario(ruta, n_eventos=400)
    print("archivo de prueba: %s  (%.1f kB, sin extension en el nombre)"
          % (ruta, _os.path.getsize(ruta) / 1e3))

    # ---- 1. el formato se decide por el CONTENIDO -------------------------
    _titulo("1. deteccion de formato (por los primeros 4 bytes, no por el nombre)")
    for etiqueta, fn in (("analisis/drs", drs.formato),
                         ("codigo/spectrum", spectrum.formato)):
        got = fn(ruta)
        ok = got == "binario"
        print("  %-18s -> %-8s  %s" % (etiqueta, got, "ok" if ok else "MAL"))
        if not ok:
            fallos.append("%s no reconocio el binario" % etiqueta)

    canal_a = drs.detectar_canal(ruta)
    canal_c = spectrum.detectar_canal(ruta)
    ok = canal_a == canal_c == ("CHN1", "Scaler0")
    print("  canal detectado    -> %s / %s  %s"
          % (canal_a, canal_c, "ok" if ok else "MAL"))
    if not ok:
        fallos.append("el canal detectado no coincide: %s vs %s"
                      % (canal_a, canal_c))

    # ---- 2. la tension y el eje temporal vuelven como se escribieron ------
    _titulo("2. reconstruccion del pulso (verdad conocida)")
    T, V, _, _ = drs_bin.leer_pulsos(ruta, n=5)
    err_v, err_t = [], []
    for k in range(5):
        t_esp = drs_bin._eje_temporal(verdad["ancho"], verdad["celdas"][k])
        v_esp = pulso_modelo(t_esp, verdad["t0"], verdad["amplitudes"][k])
        err_t.append(float(np.abs(T[k] - t_esp).max()))
        err_v.append(float(np.abs(V[k] - v_esp).max()))
    cuanto = 1000.0 / 65536.0        # el cuanto del ADC del DRS4, en mV
    print("  error maximo en t : %.3e ns   (esperado 0)" % max(err_t))
    print("  error maximo en v : %.5f mV  (cuanto del ADC %.5f, tolerancia %.5f)"
          % (max(err_v), cuanto, 0.51 * cuanto))
    if max(err_t) > 1e-9:
        fallos.append("el eje temporal no se reconstruye")
    if max(err_v) > 0.51 * cuanto:
        fallos.append("la tension no se reconstruye dentro del cuanto del ADC")

    # ---- 3. las dos cadenas sacan las MISMAS areas ------------------------
    _titulo("3. mismas areas por las dos cadenas (analisis/ y codigo/)")
    pre, post = spectrum.GATE_PRE, spectrum.GATE_POST

    d = a_area.areas_de_archivo(ruta, metodo="global", pre=pre, post=post,
                                canal="CHN1", usar_cache=False, verbose=False)
    areas_a = np.asarray(d["area"], float)

    areas_c, dur_c, extra_c = spectrum.leer_areas_y_tiempo(
        ruta, umbral_ruido=0.05, metodo="gate", usar_cache=False,
        verbose=False, devolver_extra=True)

    print("  analisis/area   : %d pulsos, mediana %.2f mV*ns"
          % (areas_a.size, np.median(areas_a)))
    print("  codigo/spectrum : %d pulsos, mediana %.2f mV*ns"
          % (areas_c.size, np.median(areas_c)))
    if areas_a.size != areas_c.size:
        fallos.append("distinto numero de pulsos: %d vs %d"
                      % (areas_a.size, areas_c.size))
    else:
        rel = np.abs(areas_c - areas_a) / np.maximum(np.abs(areas_a), 1e-9)
        print("  diferencia relativa: maxima %.2e, mediana %.2e"
              % (rel.max(), np.median(rel)))
        if rel.max() > 1e-3:
            fallos.append("las areas no coinciden (max rel %.2e)" % rel.max())

    # dos lectores pueden coincidir y estar los dos mal: contra la integral
    esperado = area_analitica(float(np.median(verdad["amplitudes"])), pre, post)
    coc = float(np.median(areas_c)) / esperado
    print("  area analitica del pulso mediano dentro de la compuerta: %.1f mV*ns"
          % esperado)
    print("  cociente medido/analitico: %.4f   (1.00 +- unos pocos %% es lo "
          "esperable: la compuerta arranca `pre` antes del onset)" % coc)
    if not (0.90 < coc < 1.10):
        fallos.append("el area absoluta no cierra con la integral analitica "
                      "(cociente %.3f)" % coc)

    # ---- 4. scaler y tiempo vivo -----------------------------------------
    _titulo("4. scaler y tiempo vivo")
    print("  duracion de pared : %.1f s" % dur_c)
    print("  tiempo vivo       : %.1f s" % float(extra_c["t_vivo"]))
    print("  scaler mediano    : %.0f Hz  (se escribieron 170)"
          % float(np.median(extra_c["scaler"])))

    # ---- 5. sin scaler (versiones viejas del DRSOSC) ----------------------
    ruta2 = _os.path.join(tmp, "Co60_sin_scaler")
    escribir_binario(ruta2, n_eventos=120, semilla=7, con_scaler=False, canal=2)
    canal2, tag2 = spectrum.detectar_canal(ruta2)
    a2, dur2, ex2 = spectrum.leer_areas_y_tiempo(
        ruta2, umbral_ruido=0.05, metodo="gate", usar_cache=False,
        verbose=False, devolver_extra=True)
    print("  archivo sin scaler: canal %s, tag %s, %d pulsos, t_vivo %.1f s "
          "(= pared, que es lo unico honesto sin scaler)"
          % (canal2, tag2, a2.size, float(ex2["t_vivo"])))
    if tag2 is not None or a2.size == 0:
        fallos.append("el archivo sin scaler no se manejo bien")

    # ---- 6. espectro sin fondo -------------------------------------------
    _titulo("6. espectro neto SIN fondo (las corridas binarias no tienen)")
    esp = spectrum.espectro_neto(ruta, None, umbral_ruido=0.05,
                                 num_canales=200, rango_max=8000.0)
    print("  hay_fondo = %s, k_m = %.3f, k_f = %.3f, cuentas netas = %d"
          % (esp["hay_fondo"], esp["k_m"], esp["k_f"], int(esp["neto"].sum())))
    if esp["hay_fondo"] or esp["k_f"] != 0.0 or esp["neto"].sum() <= 0:
        fallos.append("el espectro sin fondo no salio bien")

    # ---- resumen ---------------------------------------------------------
    print()
    print("=" * 74)
    if fallos:
        print("FALLARON %d comprobaciones:" % len(fallos))
        for f_ in fallos:
            print("  - %s" % f_)
    else:
        print("Todo ok: los dos lectores dan lo mismo sobre el mismo binario,")
        print("y el area absoluta cierra con la integral analitica del pulso.")
    print("=" * 74)

    for r in (ruta, ruta2):
        for extra in ("", ".areas.npz", ".areas_gate.npz"):
            try:
                _os.remove(r + extra)
            except OSError:
                pass
    try:
        _os.rmdir(tmp)
    except OSError:
        pass
    return 1 if fallos else 0


if __name__ == "__main__":
    _sys.exit(main())
