"""
Comprobar la Etapa 3: contra el oraculo y contra muestra/fondo fabricados.

    python desde_cero/verificar_espectro.py            (~5 min: archivos enteros)
    python desde_cero/verificar_espectro.py --rapido   (20000 eventos por archivo)

1. ORACULO. Las areas de la Etapa 2 (archivos enteros de `cs137_1600` y
   `co60`, muestra y fondo) se pasan a `desde_cero.espectro.espectro_neto` y
   a `analisis.espectro.espectro_neto`. Asi se compara SOLO la logica de esta
   etapa (tope, bines, tiempo vivo, k, resta, error, fin del umbral), con la
   entrada identica. Las marcas de tiempo se convierten a segundos con la
   funcion del oraculo para el oraculo y con la nuestra para nosotros.

2. VERDAD CONOCIDA. Una "corrida de muestra" con una fuente de tasa conocida
   mas fondo, viva solo una fraccion conocida del tiempo, y una "corrida de
   fondo" viva el 100 %, con marcas de tiempo en rafagas y un scaler con
   espigas de ruido, como los reales. Se sabe cuantas cuentas de FUENTE tiene
   que haber en cada bin; el neto tiene que caer ahi dentro de su barra
   (pulls de media 0 y ancho 1), y el tiempo vivo tiene que ser el que se
   puso.
"""

import os
import sys
from datetime import datetime, timedelta

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import area as oraculo_area                   # noqa: E402
from analisis import espectro as oraculo                    # noqa: E402
from desde_cero import area, espectro, lector               # noqa: E402


# --------------------------------------------------------------------------
# 1. Contra el oraculo
# --------------------------------------------------------------------------
def al_formato_del_oraculo(d):
    """Nuestro dict de areas con las claves que espera analisis/espectro."""
    t = oraculo_area._relojes_a_segundos(d["reloj"])
    return dict(area=d["area"], bandera=d["bandera"], reloj_s=t,
                scaler=d["scaler"], duracion=float(t[-1] - t[0]),
                n_eventos=d["n_leidos"], metodo="global", pre=d["pre"],
                post=d["post"])


def contra_oraculo(clave, limite):
    ruta_m, ruta_f = lector.CORRIDAS_TRABAJO[clave]
    d_m = area.areas_de_archivo(ruta_m, limite=limite)
    d_f = area.areas_de_archivo(ruta_f, limite=limite)
    nos = espectro.espectro_neto(d_m, d_f)
    ora = oraculo.espectro_neto(al_formato_del_oraculo(d_m),
                                al_formato_del_oraculo(d_f), bins=espectro.BINES,
                                verbose=False)
    rel = lambda a, b: abs(a - b) / abs(b)                  # noqa: E731
    return dict(
        canales=float(np.abs(nos["canales"] - ora["canales"]).max()),
        neto=float(np.abs(nos["neto"] - ora["neto"]).max()),
        error=float(np.abs(nos["error"] - ora["error"]).max()),
        vivo_m=rel(nos["t_vivo_m"], ora["t_vivo_muestra"]),
        vivo_f=rel(nos["t_vivo_f"], ora["t_vivo_fondo"]),
        k_f=rel(nos["k_f"], ora["k_f"]),
        umbral=(espectro.fin_del_umbral(nos["canales"], nos["neto"]),
                oraculo.fin_del_umbral(ora["canales"], ora["neto"])),
        nos=nos)


# --------------------------------------------------------------------------
# 2. Fabricados
# --------------------------------------------------------------------------
T0 = datetime(2026, 9, 27, 12, 0, 0)


def _relojes(n, pared, rng, rafaga=200):
    """
    n marcas ordenadas en [0, pared] s, escritas en RAFAGAS: cada `rafaga`
    eventos seguidos comparten la marca, al milisegundo, como en el DRSOSC.
    """
    t = np.sort(rng.uniform(0.0, pared, n))
    t = np.repeat(t[::rafaga], rafaga)[:n]
    t[0], t[-1] = 0.0, pared
    # Se redondea UNA vez, a milisegundos enteros, y despues se formatea.
    # Redondear los milisegundos aparte de los segundos (lo que hacia esto
    # antes) convierte 12.9996 s en "12.000": la marca retrocede un segundo,
    # el dt negativo se recorta a 0 y el siguiente vale 1 s de mas. Salia el
    # tiempo vivo 0.5 % bajo sin espigas: error del FABRICADO, no del modulo.
    out = []
    for s in t:
        d = T0 + timedelta(milliseconds=int(round(1000.0 * float(s))))
        out.append(d.strftime("%Y/%m/%d %H:%M:%S.") + "%03d" % (d.microsecond // 1000))
    return out


def _corrida(tasa_fuente, tasa_fondo, pared, viva, fuente, rng, espigas=0.0,
             rafaga=200):
    """
    Una corrida fabricada. Devuelve (dict como el de area.areas_de_archivo,
    cuentas de fuente guardadas por area).

    La placa dispara a tasa_fuente + tasa_fondo durante `pared` segundos, y
    se guarda una fraccion `viva` de los disparos. El scaler de cada evento
    es la tasa real con ruido de Poisson (Hz enteros) y, con probabilidad
    `espigas`, una rafaga de 90 kHz como las del fondo real.
    """
    tasa = tasa_fuente + tasa_fondo
    n = rng.poisson(viva * tasa * pared)
    de_fuente = rng.random(n) < tasa_fuente / tasa
    a = np.where(de_fuente, fuente(n, rng), rng.exponential(120.0, n))
    sc = rng.poisson(tasa, n).astype(float)
    sc[rng.random(n) < espigas] = 9.0e4
    d = dict(area=a, bandera=np.zeros(n, int),
             reloj=_relojes(n, pared, rng, rafaga),
             scaler=sc, n_leidos=n)
    return d, a[de_fuente]


def _fuente(n, rng):
    """Un "Compton" de juguete: plateau hasta 600 mV*ns con borde borroneado."""
    return np.abs(rng.uniform(0, 600, n) + rng.normal(0, 25, n))


def fabricado(rng, tasa_fuente=150.0, tasa_fondo=10.0, pared_m=600.0,
              viva_m=0.40, pared_f=900.0, espigas=0.0):
    """
    La muestra escribe sus marcas en rafagas de 200 (tasa alta, como el
    Co-60 real); el fondo NO (a 10 ev/s cada evento tiene su marca: medido,
    89-100 % de los eventos de los fondos reales tienen dt > 0).
    """
    d_m, fuente_m = _corrida(tasa_fuente, tasa_fondo, pared_m, viva_m, _fuente, rng)
    d_f, _ = _corrida(0.0, tasa_fondo, pared_f, 1.0, _fuente, rng,
                      espigas=espigas, rafaga=1)
    e = espectro.espectro_neto(d_m, d_f)
    # Lo que tendria que dar: las cuentas de FUENTE que la muestra guardo,
    # llevadas a t_norm. Se usa la esperanza, no las cuentas fabricadas, para
    # que el pull mida solo la estadistica que el error dice cubrir.
    vivo_verdadero = viva_m * pared_m
    t_norm = min(vivo_verdadero, pared_f)
    p, _ = np.histogram(_fuente(2_000_000, rng), bins=e["bordes"])
    esperado = tasa_fuente * t_norm * p / 2_000_000
    m = esperado > 20
    pull = (e["neto"][m] - esperado[m]) / e["error"][m]
    return dict(vivo=e["t_vivo_m"], vivo_verdadero=vivo_verdadero,
                vivo_f=e["t_vivo_f"], pared_f=pared_f, k_f=e["k_f"],
                k_f_verdadero=t_norm / pared_f, pull_media=pull.mean(),
                pull_ancho=pull.std(ddof=1), n_bines=int(m.sum()))


def main():
    rapido = "--rapido" in sys.argv
    limite = 20000 if rapido else None
    ok = True
    print("=" * 78)
    print("Etapa 3 contra analisis/espectro.py  (mismas areas, %s)"
          % ("20000 eventos por archivo" if rapido else "archivos enteros"))
    print("=" * 78)
    for clave in ("cs137_1600", "co60"):
        ruta_m, ruta_f = lector.CORRIDAS_TRABAJO[clave]
        if not (os.path.exists(ruta_m) and os.path.exists(ruta_f)):
            print("  %-12s no esta en disco" % clave)
            continue
        r = contra_oraculo(clave, limite)
        # misma cuenta con la misma entrada: tiene que dar lo mismo salvo
        # redondeo de punto flotante; el criterio del PLAN para el tiempo
        # vivo es 0.1 %
        bien = (r["canales"] < 1e-9 and r["neto"] < 1e-9 and r["error"] < 1e-9
                and r["vivo_m"] < 1e-3 and r["vivo_f"] < 1e-3
                and r["umbral"][0] == r["umbral"][1])
        ok &= bien
        print("  %-12s max|d canal| %.1e  max|d neto| %.1e  max|d error| %.1e  "
              "tiempo vivo %.1e / %.1e (rel)  fin del umbral %.1f / %.1f   %s"
              % (clave, r["canales"], r["neto"], r["error"], r["vivo_m"],
                 r["vivo_f"], r["umbral"][0], r["umbral"][1],
                 "ok" if bien else "DIFIEREN"))
        print(espectro.resumen(clave, r["nos"]))

    print()
    print("=" * 78)
    print("Muestra y fondo fabricados (fuente 150 Hz + fondo 10 Hz, muestra viva "
          "el 40 %)")
    print("=" * 78)
    rng = np.random.default_rng(4)
    medias, anchos = [], []
    print("  a. scaler limpio")
    for k in range(5):
        f = fabricado(rng)
        medias.append(f["pull_media"])
        anchos.append(f["pull_ancho"])
        print("  realizacion %d: vivo %.2f s (verdad %.1f)  fondo vivo %.1f "
              "(verdad %.1f)  k_f %.4f (verdad %.4f)  pulls en %d bines: "
              "media %+.3f, ancho %.3f"
              % (k, f["vivo"], f["vivo_verdadero"], f["vivo_f"], f["pared_f"],
                 f["k_f"], f["k_f_verdadero"], f["n_bines"], f["pull_media"],
                 f["pull_ancho"]))
        # el tiempo vivo sale de n / sum(scaler*dt): su error relativo es el
        # de Poisson del conteo de eventos, ~1/sqrt(n) ~ 0.2 %
        ok &= abs(f["vivo"] / f["vivo_verdadero"] - 1) < 0.01
        ok &= abs(f["vivo_f"] / f["pared_f"] - 1) < 0.01
    m, a = float(np.mean(medias)), float(np.mean(anchos))
    # media de los pulls sobre ~5 x 100 bines: su error es ~1/sqrt(500) ~ 0.05
    bien = abs(m) < 0.15 and 0.9 < a < 1.1
    ok &= bien
    print("  promedio: pulls de media %+.3f y ancho %.3f   %s"
          % (m, a, "ok" if bien else "FALLA"))

    # b. Las espigas del scaler. El recorte las LIMITA a 10 veces la mediana,
    # no las elimina: cada espiga todavia cuenta ~9 veces la tasa de mas, asi
    # que una fraccion p de espigas infla los disparos en ~9p y baja el
    # tiempo vivo en lo mismo. No se le pide nada: se mide, porque en los
    # datos reales NO pasa (ver PLAN.md, bitacora de la Etapa 3: a tasa baja
    # el scaler real lee 20-30 % de menos, la fraccion viva sale > 1 y el
    # tope la deja en 1).
    p = 0.002
    print("  b. scaler con %.1f %% de espigas de 90 kHz en el fondo" % (100 * p))
    # fondo 10 veces mas largo: con ~9000 eventos la fluctuacion de Poisson
    # (1 %) es del tamano del sesgo que se quiere ver
    # (y promedio de 4: cada una fluctua ~0.4 %)
    fs = [fabricado(rng, espigas=p, pared_f=9000.0) for _ in range(4)]
    vivo = float(np.mean([f["vivo_f"] for f in fs]))
    pred = 9000.0 / (1 + p * (espectro.RECORTE_SCALER - 1))
    bien = abs(vivo / pred - 1) < 0.004
    ok &= bien
    print("     fondo vivo %.1f s en promedio (verdad 9000.0, predicho con el "
          "sesgo del recorte %.1f)   %s" % (vivo, pred, "ok" if bien else "FALLA"))
    print()
    print("Todo en orden." if ok else "Hay algo que revisar.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
