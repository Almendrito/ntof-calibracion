"""
verificar_gui.py — que `gui.py` corra entero y que no invente ningun numero.

    python replica_paper/verificar_gui.py            las dos cosas
    python replica_paper/verificar_gui.py --rapido   sin las pestanas lentas
    python replica_paper/verificar_gui.py --humo     solo que corra
    python replica_paper/verificar_gui.py --numeros  solo los numeros

Son dos chequeos distintos y los dos hacen falta.

1. HUMO. Instancia la `App` y, pestana por pestana, llama `leer_controles()`,
   `trabajo()` y `pintar()` EN EL HILO PRINCIPAL, sin `mainloop`. Es la forma
   barata de ejercitar una interfaz entera sin tocarla: no simula clicks, llama
   a lo mismo que el boton «Recalcular» termina llamando. Encontro, por
   ejemplo, que `"50% de altura %8.1f"` es un format string de TRES
   conversiones y no de una, porque `% d` es un `%d` con flag de espacio.

2. NUMEROS. La regla de `gui.py` es que no calcula nada por su cuenta: todo
   sale de llamar a los modulos de esta carpeta. La forma de verificarlo es
   comparar lo que la interfaz devuelve contra `informe.txt`, que lo genero
   `correr.py`. Si algun default de la interfaz se desviara -y hay uno que se
   presta, `ligar_origen`, que vale True en `correr.tres_metodos` y False en
   `metodo.ajustar_dos_bordes`- el numero saldria distinto sin que nada falle.
"""

import argparse
import os
import sys
import time
import traceback

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "replica_paper"

from replica_paper import datos as dat
from replica_paper import gui
from replica_paper import rutas

# Las pestanas mas caras. Con --rapido no se corren; el resto tarda ~1 min.
LENTAS = {"8 · Fig. 5 · Robustez", "9 · II-C · Las 10 mediciones",
          "11 · Sec. III · La recta", "12 · Cierre con verdad conocida"}

# Lo que dice `replica_paper/informe.txt`, que genero `correr.py`. Es la
# referencia: la interfaz tiene que devolver esto mismo.
INFORME = {
    "cs137_1600": dict(lo=579.0, hi=1389.0, rho=0.99618, sigma_rel=0.16,
                       B=915.2, sigma=115.5, chi2=2.02,
                       gauss=951.7, dif=829.0),
    "co60": dict(lo=952.0, hi=1967.0, rho=0.99095, sigma_rel=0.25,
                 B=1124.2, B2=1304.7, D=1.1669, chi2=1.18,
                 gauss=1341.3, dif=1110.9),
}

TOL_REL = 0.005          # medio por ciento: los del informe estan redondeados


def humo(clave, rapido=False):
    """Corre cada pestana de punta a punta. Devuelve cuantas fallaron."""
    print("=" * 74)
    print("HUMO - las %d pestanas, sobre %s" % (len(gui.PESTANAS), clave))
    print("=" * 74)
    app = gui.App()
    app.withdraw()
    fallo = 0
    for p in app.pestanas:
        if rapido and p.titulo in LENTAS:
            print("  saltada  %-34s (--rapido)" % p.titulo)
            continue
        # La pestana 0 es la que manda sobre que corrida se usa: su combo
        # sobrescribe `app.estado` al leer los controles.
        if isinstance(p, gui.TabDatos):
            p.v_cor.set(gui.etiqueta_corrida(clave))
        t0 = time.time()
        try:
            p.cfg = p.leer_controles()
            r = p.trabajo()
            p.fig.clear()
            p.pintar(r)
            print("  ok       %-34s %6.1f s" % (p.titulo, time.time() - t0))
        except Exception:                                    # noqa: BLE001
            print("  FALLO    %-34s %6.1f s" % (p.titulo, time.time() - t0))
            print(traceback.format_exc())
            fallo += 1
    app.destroy()
    return fallo


def numeros(claves=None):
    """Compara la interfaz contra informe.txt, valor por valor."""
    print()
    print("=" * 74)
    print("NUMEROS - la interfaz contra replica_paper/informe.txt")
    print("=" * 74)
    mal = 0
    for clave in (claves or list(INFORME)):
        esp = INFORME.get(clave)
        if esp is None:
            print("  %s no esta en el informe: no hay contra que comparar" % clave)
            continue
        if clave not in dat.disponibles():
            print("  %s no esta accesible ahora (disco externo?)" % clave)
            continue
        e = dat.espectro(clave)
        bordes, r = gui.tres_metodos(e["canales"], e["cuentas"], e["lineas"])
        filas = [("NCC lo", r["intervalo"][0], esp["lo"]),
                 ("NCC hi", r["intervalo"][1], esp["hi"]),
                 ("rho", r["rho"], esp["rho"]),
                 ("sigma_rel de S2", r["sigma_rel_s2"], esp["sigma_rel"]),
                 ("chi2/ndf", r["chi2_ndf"], esp["chi2"]),
                 ("Gauss FWHM", bordes["Gauss FWHM"], esp["gauss"]),
                 ("Diferencia", bordes["Diferencia"], esp["dif"])]
        if r["tipo"] == "un_borde":
            filas += [("B", r["B"], esp["B"]),
                      ("sigma", r["sigma"], esp["sigma"])]
        else:
            filas += [("B1", r["B1"], esp["B"]), ("B2", r["B2"], esp["B2"]),
                      ("D", r["D"], esp["D"])]
        print("  %s" % clave)
        for nombre, tengo, quiero in filas:
            tol = max(abs(quiero) * TOL_REL, 0.005)
            ok = abs(tengo - quiero) <= tol
            mal += 0 if ok else 1
            print("     %-16s %12.4f   informe %10.4f   %s"
                  % (nombre, tengo, quiero, "ok" if ok else "<-- DISTINTO"))
    return mal


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rapido", action="store_true",
                    help="sin las pestanas que tardan minutos")
    ap.add_argument("--humo", action="store_true", help="solo el humo")
    ap.add_argument("--numeros", action="store_true", help="solo los numeros")
    ap.add_argument("--corrida", default=None,
                    help="clave del catalogo (por defecto, la primera con fondo)")
    args = ap.parse_args()
    rutas.salida_utf8()

    hay = dat.con_fondo() or dat.disponibles()
    if not hay:
        print("No hay ninguna corrida accesible: nada que verificar.")
        return 2
    clave = args.corrida or hay[0]

    fallo = 0
    if not args.numeros:
        fallo += humo(clave, rapido=args.rapido)
    if not args.humo:
        fallo += numeros()
    print()
    print("=" * 74)
    print("%s" % ("TODO OK" if fallo == 0 else "%d PROBLEMAS" % fallo))
    print("=" * 74)
    return 1 if fallo else 0


if __name__ == "__main__":
    sys.exit(main())
