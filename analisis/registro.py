"""
registro.py — las calibraciones que se van acumulando, por detector.

Para que sirve
--------------
Una calibracion suelta da D para UN detector con UNA fuente. Lo que hoy falta en
este proyecto es el **test de linealidad**, y no se puede hacer con una sola
fuente: el mapa canal-energia pasa por el origen, asi que un ajuste entrega un
solo parametro y todos los bordes quedan atados a el. Los dos bordes del Co-60
NO son dos puntos independientes (es la Ec. 12 del paper: el segundo borde se
amarra al primero con la misma D).

Con DOS fuentes distintas medidas en el MISMO detector si lo son. Este modulo
guarda cada calibracion con su detector y, cuando hay dos fuentes del mismo
detector, arma los puntos y ajusta la recta de las dos maneras que importan:

  * por el ORIGEN, que es lo que el metodo supone;
  * con ORDENADA LIBRE, que es lo que la testea. Si la ordenada sale compatible
    con cero, la suposicion se sostiene; si no, hay un pedestal y la calibracion
    de una sola fuente esta sesgada.

Archivo: `.cache/registro_calibraciones.json`.
"""

import json
import os
import sys
import time

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas

ARCHIVO = "registro_calibraciones.json"


def _ruta():
    return rutas.cache(ARCHIVO)


def huella(corrida):
    """
    Identifica los DATOS de una calibracion: nombre, tamano y fecha de cada
    archivo. Dos calibraciones del mismo detector y la misma fuente pero con
    archivos distintos son entradas distintas, no una que pisa a la otra.
    """
    partes = []
    for r in (corrida.muestra, corrida.fondo, corrida.sim):
        if not r or not os.path.exists(r):
            partes.append("-")
            continue
        st = os.stat(r)
        partes.append("%s:%d:%d" % (os.path.basename(r), st.st_size,
                                    int(st.st_mtime)))
    return "|".join(partes)


def punto(r):
    """
    El punto de operacion de una entrada: detector + HV.

    Se lee del campo `punto` que guarda `calibracion.calibrar`; las entradas
    viejas no lo tienen y se reconstruye del detector, que es lo que habia.
    """
    p = (r.get("punto") or "").strip()
    if p:
        return p
    det = (r.get("detector") or "sin nombre").strip()
    return "%s @ %g V" % (det, r["hv"]) if r.get("hv") else det


def escala(r):
    return "keVee" if r["rama"].startswith("Lvis") else "keV"


def clave(resultado, corrida=None):
    h = resultado.get("huella") or (huella(corrida) if corrida else "")
    return "%s | %s | %s | %s" % (punto(resultado), resultado["fuente"],
                                  escala(resultado), h[:40])


def cargar():
    try:
        with open(_ruta(), encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def guardar(resultado, corrida=None, detector=None):
    """Agrega o actualiza una entrada. Devuelve (clave, registro completo)."""
    r = dict(resultado)
    if detector:
        r["detector"] = detector
    if corrida is not None:
        r["huella"] = huella(corrida)
        r["archivos"] = dict(
            muestra=rutas.relativa(corrida.muestra) if corrida.muestra else "",
            fondo=rutas.relativa(corrida.fondo) if corrida.fondo else "",
            sim=rutas.relativa(corrida.sim) if corrida.sim else "")
    r["guardado"] = time.strftime("%Y-%m-%d %H:%M")
    reg = cargar()
    k = clave(r, corrida)
    reg[k] = r
    with open(_ruta(), "w", encoding="utf-8") as f:
        json.dump(reg, f, indent=2, ensure_ascii=False)
    return k, reg


def borrar(k):
    reg = cargar()
    reg.pop(k, None)
    with open(_ruta(), "w", encoding="utf-8") as f:
        json.dump(reg, f, indent=2, ensure_ascii=False)
    return reg


def por_punto(reg=None, escala_pedida="keVee", metodo="global"):
    """
    Agrupa por PUNTO DE OPERACION (detector + HV): es el test de linealidad.

    Filtra tambien por metodo de integracion, no solo por escala. Sobre el
    Cs-137 a 1600 V, "global" y "convergencia" dan D = 1.945 y 1.611: un 17% de
    diferencia, muchisimo mas que cualquier no linealidad que uno quiera medir.
    Mezclarlos en una recta seria medir la eleccion del metodo.
    """
    reg = cargar() if reg is None else reg
    d = {}
    for k, r in reg.items():
        if escala(r) != escala_pedida:
            continue
        if metodo is not None and r.get("metodo", "global") != metodo:
            continue
        if rutas.descartada(r.get("corrida")):
            continue
        d.setdefault(punto(r), []).append((k, r))
    return d


def por_fuente(reg=None, escala_pedida="keVee", metodo="global"):
    """
    Agrupa por (detector, fuente): es la ganancia contra voltaje.

    Al reves que `por_punto`: aca lo que se mantiene fijo es la FUENTE y lo que
    varia es el HV, asi que cada grupo son los puntos de la curva D(V) de un
    detector medidos con la misma linea gamma.
    """
    reg = cargar() if reg is None else reg
    d = {}
    for k, r in reg.items():
        if escala(r) != escala_pedida:
            continue
        if metodo is not None and r.get("metodo", "global") != metodo:
            continue
        if rutas.descartada(r.get("corrida")):
            continue
        if r.get("hv") is None:
            continue
        det = (r.get("detector") or "sin nombre").strip()
        d.setdefault((det, r["fuente"]), []).append((k, r))
    return d


# --------------------------------------------------------------------------
# El test de linealidad
# --------------------------------------------------------------------------
def puntos(entradas):
    """
    De una lista de calibraciones saca los puntos (E_borde, canal, error).

    Se toma UN punto por fuente: el borde de referencia. Los demas bordes de una
    misma fuente estan atados al mismo D y no aportan informacion nueva; meterlos
    como puntos independientes seria inflar artificialmente la evidencia.
    """
    E, ch, err, etiq = [], [], [], []
    for k, r in entradas:
        E.append(r["E_ref"])
        ch.append(r["D"] * r["E_ref"])
        err.append(r["err_D"] * r["E_ref"])
        etiq.append(r["fuente"])
    o = np.argsort(E)
    return (np.array(E)[o], np.array(ch)[o], np.array(err)[o],
            [etiq[i] for i in o])


def linealidad(entradas):
    """
    Ajusta la recta canal(E) por el origen y con ordenada libre.

    Devuelve None si no hay al menos dos FUENTES distintas: con una sola, la
    recta con ordenada libre pasa exacto por el punto y por el origen y el test
    no dice nada.
    """
    E, ch, err, etiq = puntos(entradas)
    if len(set(etiq)) < 2:
        return None
    w = 1.0 / np.maximum(err, 1e-12) ** 2

    # por el origen: ch = D*E  ->  D = sum(w E ch) / sum(w E^2)
    D0 = float(np.sum(w * E * ch) / np.sum(w * E * E))
    res0 = ch - D0 * E

    # con ordenada libre: ch = a + D*E, minimos cuadrados ponderados
    S, Sx = np.sum(w), np.sum(w * E)
    Sxx, Sy, Sxy = np.sum(w * E * E), np.sum(w * ch), np.sum(w * E * ch)
    det = S * Sxx - Sx * Sx
    a = float((Sxx * Sy - Sx * Sxy) / det)
    D1 = float((S * Sxy - Sx * Sy) / det)
    err_a = float(np.sqrt(Sxx / det))
    err_D1 = float(np.sqrt(S / det))
    res1 = ch - (a + D1 * E)

    return dict(E=E, canal=ch, err=err, fuentes=etiq,
                D_origen=D0, residuos_origen=res0,
                ordenada=a, err_ordenada=err_a, D_libre=D1, err_D_libre=err_D1,
                residuos_libre=res1,
                compatible_con_origen=bool(abs(a) < 2 * err_a),
                n_fuentes=len(set(etiq)))


def informe(escala_pedida="keVee", metodo="global"):
    reg = cargar()
    if not reg:
        print("  el registro esta vacio")
        return
    grupos = por_punto(reg, escala_pedida, metodo)
    print("  %d calibraciones guardadas, %d en escala %s con metodo '%s'"
          % (len(reg), sum(len(v) for v in grupos.values()), escala_pedida,
             metodo))
    ordenadas = []
    for det, entradas in sorted(grupos.items()):
        print()
        print("  PUNTO DE OPERACION: %s   (%d calibracion(es))"
              % (det, len(entradas)))
        for k, r in entradas:
            print("     %-8s D = %.4f +- %.4f   borde %7.1f mV*ns   "
                  "chi2/ndf %.2f   %s"
                  % (r["fuente"], r["D"], r["err_D"], r["D"] * r["E_ref"],
                     r["chi2_ndf"], r.get("guardado", "")))
        lin = linealidad(entradas)
        if lin is None:
            print("     -> con una sola fuente NO se puede testear linealidad.")
            print("        Los bordes de una misma fuente estan atados a la "
                  "misma D por el mapa por el origen.")
            continue
        print("     TEST DE LINEALIDAD (%d fuentes):" % lin["n_fuentes"])
        print("        recta por el origen : D = %.4f" % lin["D_origen"])
        for E, f, r_ in zip(lin["E"], lin["fuentes"], lin["residuos_origen"]):
            print("           %-8s E = %7.1f  residuo = %+7.1f mV*ns"
                  % (f, E, r_))
        print("        ordenada libre      : D = %.4f +- %.4f,  "
              "ordenada = %+.1f +- %.1f mV*ns"
              % (lin["D_libre"], lin["err_D_libre"], lin["ordenada"],
                 lin["err_ordenada"]))
        print("        -> la ordenada %s con cero."
              % ("ES compatible" if lin["compatible_con_origen"]
                 else "NO es compatible"))
        if lin["compatible_con_origen"]:
            print("           El mapa por el origen se sostiene en este punto.")
        else:
            print("           Con UN punto de operacion no se puede decir si es")
            print("           un pedestal o una diferencia de ganancia entre las")
            print("           corridas: ver el chequeo del final.")
        ordenadas.append((det, entradas, lin))

    _pedestal_o_ganancia(ordenadas)


def _pedestal_o_ganancia(ordenadas):
    """
    Una ordenada distinta de cero en UN punto de operacion no prueba que haya
    pedestal.

    Un pedestal del digitalizador es una carga fija por evento: no depende del
    alto voltaje, asi que tiene que valer LO MISMO en canal a todos los puntos.
    Si en cambio las corridas que se combinan estan a ganancias distintas -por
    ejemplo tomadas en sesiones distintas con fuentes de poder distintas- la
    ordenada aparente CRECE con el voltaje, porque lo que es constante es el
    COCIENTE de las D.

    Con un solo punto de operacion las dos hipotesis son indistinguibles; con
    dos o mas se separan solas, y este chequeo es lo que las separa.
    """
    utiles = [(d, e, l) for d, e, l in ordenadas if l is not None]
    if len(utiles) < 2:
        return
    a = np.array([l["ordenada"] for _, _, l in utiles])
    ea = np.array([l["err_ordenada"] for _, _, l in utiles])
    w = 1.0 / np.maximum(ea, 1e-12) ** 2
    m = float(np.sum(w * a) / np.sum(w))
    chi2 = float(np.sum(w * (a - m) ** 2))
    gdl = len(utiles) - 1
    print()
    print("  PEDESTAL O DIFERENCIA DE GANANCIA?")
    print("     Un pedestal no depende del alto voltaje: tiene que dar la misma")
    print("     ordenada en canal a todos los puntos de operacion.")
    for d, _, l in utiles:
        print("        %-14s ordenada = %+8.1f +- %.1f mV*ns   (D libre %.4f)"
              % (d, l["ordenada"], l["err_ordenada"], l["D_libre"]))
    print("     media ponderada %+.1f mV*ns, chi2 = %.0f con %d gdl"
          % (m, chi2, gdl))
    if chi2 < 2.0 * gdl + 5.0:
        print("     -> compatible con un PEDESTAL comun. Conviene ajustarlo y")
        print("        restarlo, o dejar la ordenada libre en la calibracion.")
    else:
        print("     -> NO es un pedestal: la ordenada cambia con el voltaje.")
        print("        Lo que se ve es que las corridas que se estan combinando")
        print("        no comparten la ganancia. Los cocientes de D son:")
        for d, e, _ in utiles:
            por_fuente_ = {r["fuente"]: r["D"] for _, r in e}
            if len(por_fuente_) == 2:
                (f1, D1), (f2, D2) = sorted(por_fuente_.items(),
                                            key=lambda kv: kv[1])
                print("           %-14s %s/%s = %.3f" % (d, f2, f1, D2 / D1))
        print("        Hasta que las dos fuentes se midan en la MISMA sesion,")
        print("        este test de linealidad no se puede cerrar.")


# --------------------------------------------------------------------------
# Rehacer el registro
# --------------------------------------------------------------------------
def rehacer(escalas=("Lvis_keVee", "Edep_keV"), metodo="global", n_arranques=8):
    """
    Recalcula y regraba todas las corridas del catalogo que tengan espectro.

    Hace falta cada vez que cambia el MODELO o la clave del registro: una
    entrada vieja sigue ahi con su numero de antes y nada avisa que ya no es
    comparable con las nuevas. Paso el 28-ago-2026 dos veces en el mismo dia:
    una entrada de Cs-137 guardada con el metodo `convergencia` (17% por debajo
    del `global`, ver `por_punto`) y las claves viejas, que llamaban al mismo
    detector de dos maneras distintas segun la corrida.

    Las entradas superadas -mismo corrida + escala, otra clave- se borran, y se
    dice cual y por que. No se vacia el registro entero: puede haber
    calibraciones de archivos sueltos que el catalogo no conoce y que nadie
    puede recalcular desde aca.
    """
    from analisis import calibracion, espectro

    hechas, borradas = [], []
    for clave_c in sorted(rutas.CORRIDAS):
        c = rutas.corrida(clave_c)
        nombre = "espectro_%s_%s" % (c.clave, metodo)
        if not os.path.exists(rutas.cache("%s.npz" % nombre)):
            print("  %-14s sin espectro '%s': se saltea" % (c.clave, metodo))
            continue
        if not (c.sim and os.path.exists(c.sim)):
            print("  %-14s sin simulacion: se saltea" % c.clave)
            continue
        e = espectro.cargar(nombre)
        for rama in escalas:
            r, _, _, _ = calibracion.calibrar(c, rama=rama, metodo=metodo,
                                              n_arranques=n_arranques,
                                              verbose=False, esp=e)
            k, _ = guardar(r, c, detector=c.detector)
            hechas.append(k)
            print("  %-14s %-11s D = %.4f +- %.4f   -> %s"
                  % (c.clave, escala(r), r["D"], r["err_D"], k))

    reg = cargar()
    for k, r in list(reg.items()):
        if k in hechas:
            continue
        gemela = [h for h in hechas
                  if reg[h]["corrida"] == r.get("corrida")
                  and escala(reg[h]) == escala(r)]
        if gemela:
            borrar(k)
            borradas.append(k)
            print("  superada y borrada: %s" % k)
    print("  %d entradas escritas, %d borradas" % (len(hechas), len(borradas)))
    return hechas, borradas


# --------------------------------------------------------------------------
# Cierre con verdad conocida
# --------------------------------------------------------------------------
def _entrada(fuente, E_ref, D, err, hv=1550):
    """Una entrada de registro fabricada, para probar el test sin medir."""
    return ("sintetica|%s" % fuente,
            dict(corrida="sint_%s" % fuente, fuente=fuente, detector="X",
                 punto="X @ %g V" % hv, hv=hv, rama="Lvis_keVee",
                 metodo="global", D=D, err_D=err, E_ref=E_ref,
                 chi2_ndf=1.0, guardado="sintetico"))


def cierre(D_real=1.25, c0_real=0.0, err_rel=0.004):
    """
    Verifica el test de linealidad con verdad conocida, sin tocar el registro.

    Fabrica los puntos de cuatro fuentes con `canal = D*E + c0` y comprueba que
    el ajuste recupera D y detecta (o no) la ordenada. Con `c0_real = 0` la
    ordenada tiene que salir compatible con cero; con un pedestal metido a mano,
    NO. Un test que solo se probo con datos que pasan no es un test.
    """
    E = dict([("Na-22", 340.7), ("Cs-137", 477.3), ("Mn-54", 639.5),
              ("Co-60", 963.4)])
    entradas = []
    for f, e in E.items():
        ch = D_real * e + c0_real
        entradas.append(_entrada(f, e, ch / e, err_rel * ch / e))
    lin = linealidad(entradas)
    print("  verdad: D = %.4f, c0 = %+.1f mV*ns" % (D_real, c0_real))
    print("     recuperado por el origen : D = %.4f" % lin["D_origen"])
    print("     ordenada libre           : D = %.4f, c0 = %+.1f +- %.1f"
          % (lin["D_libre"], lin["ordenada"], lin["err_ordenada"]))
    print("     residuos por el origen   : %s"
          % "  ".join("%+.1f" % r for r in lin["residuos_origen"]))
    print("     -> ordenada %s con cero  (esperado: %s)"
          % ("compatible" if lin["compatible_con_origen"] else "NO compatible",
             "compatible" if abs(c0_real) < 1e-9 else "NO compatible"))
    return lin


if __name__ == "__main__":
    rutas.salida_utf8()
    print("registro:", rutas.relativa(_ruta()))
    if "--cierre" in sys.argv:
        print("cierre con verdad conocida:")
        cierre(c0_real=0.0)
        print()
        cierre(c0_real=-60.0)
        print()
    if "--rehacer" in sys.argv:
        print("recalculando todo el catalogo...")
        rehacer()
        print()
    informe()
