"""
Comprobar la Etapa 0 contra un lector que ya sabemos que funciona.

    python desde_cero/verificar.py

Por que existe
--------------
`analisis/` ya lee estos archivos y sus resultados estan publicados y revisados.
Mientras construimos de cero, eso es un ORACULO: no hace falta creerle al codigo
nuevo, se lo puede contrastar contra uno que da la respuesta correcta.

Sirve mientras exista `analisis/`. Cuando la version nueva tenga sus propias
pruebas con verdad conocida (datos fabricados por nosotros), este archivo pasa a
ser opcional. Pero para la primera etapa es lo mas barato que hay.

Ya encontro un error real: `_eventos_bin` leia el "rango" del evento en el byte
14 en vez del 22. El byte 14 es la HORA, que tambien es un numero chico y
plausible, asi que nada fallaba: los pulsos salian corridos ~16 mV. Sin comparar
contra algo, ese error habria viajado hasta el espectro.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import drs_bin, drs_xml            # el oraculo  # noqa: E402
from desde_cero import lector                    # lo nuestro  # noqa: E402


def comparar(ruta, n=200):
    """
    Lee los primeros `n` eventos con los dos lectores y mide la diferencia.

    Dos pasadas por el nuestro: la de siempre, `(t, v)`, y la de
    `con_reloj=True`, que tiene que dar el MISMO pulso mas reloj y scaler
    iguales a los del oraculo (el scaler puede ser nan en los dos).
    """
    fmt = lector.formato(ruta)
    oraculo = drs_bin if fmt == "binario" else drs_xml
    dt = dv = 0.0
    k = relojes = scalers = 0
    for (t1, v1), (t1b, v1b, r1, s1), (t2, v2, r2, s2) in zip(
            lector.eventos(ruta, limite=n),
            lector.eventos(ruta, limite=n, con_reloj=True),
            oraculo.iter_eventos(ruta, limite=n)):
        dt = max(dt, float(np.abs(t1 - t2).max()), float(np.abs(t1b - t2).max()))
        dv = max(dv, float(np.abs(v1 - v2).max()), float(np.abs(v1b - v2).max()))
        relojes += r1 != r2
        scalers += not (s1 == s2 or (np.isnan(s1) and np.isnan(s2)))
        k += 1
    return dict(formato=fmt, eventos=k, dt=dt, dv=dv, relojes=relojes,
                scalers=scalers, ultimo=(r1, s1))


def main():
    # las muestras y los fondos de las corridas de trabajo: los dos formatos
    # y los dos canales/scalers (CHN2+Scaler1 en junio, CHN1+Scaler0 despues)
    rutas = [r for par in lector.CORRIDAS_TRABAJO.values() for r in par
             if r is not None and os.path.exists(r)]
    if not rutas:
        print("No hay archivos de ejemplo a mano.")
        return 0
    print("=" * 74)
    print("Etapa 0 contra el lector de analisis/  (los dos deben dar lo MISMO)")
    print("=" * 74)
    malos = 0
    for ruta in rutas:
        r = comparar(ruta)
        # Tolerancia CERO a proposito: los dos lectores hacen la misma cuenta
        # sobre los mismos bytes, asi que cualquier diferencia es un error de
        # alguno de los dos, no ruido numerico.
        ok = (r["dt"] == 0.0 and r["dv"] == 0.0 and r["relojes"] == 0
              and r["scalers"] == 0)
        malos += not ok
        print("  %-38s [%-7s] %3d ev   max|dt| = %.1e ns   max|dv| = %.1e mV   "
              "relojes distintos %d   scalers distintos %d   %s"
              % (os.path.basename(ruta)[:38], r["formato"], r["eventos"],
                 r["dt"], r["dv"], r["relojes"], r["scalers"],
                 "ok" if ok else "DIFIEREN"))
        print("  %-38s ultimo evento: reloj %s, scaler %s Hz"
              % ("", r["ultimo"][0], r["ultimo"][1]))
    print()
    print("Todo coincide." if not malos else "Hay %d archivo(s) que difieren." % malos)
    return 1 if malos else 0


if __name__ == "__main__":
    sys.exit(main())
