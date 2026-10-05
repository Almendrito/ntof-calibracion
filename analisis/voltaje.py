"""
voltaje.py — la ganancia contra el alto voltaje, D(V).

Que pregunta contesta, y cual NO
--------------------------------
Hay dos linealidades en este proyecto y se confunden todo el tiempo:

  TEST A — el EJE:  `canal = D * E`.  Es el que valida la calibracion. Necesita
      varias FUENTES al mismo voltaje. Vive en `registro.linealidad`.
  TEST B — la GANANCIA:  `D(V)`.  Caracteriza el fotodetector, no la
      calibracion. Necesita la MISMA fuente a varios voltajes. Es este modulo.

Y la respuesta esperada al test B no es "si": para un tubo fotomultiplicador la
ganancia es una LEY DE POTENCIA,

    G(V) = a * V^b        con  b = k*n  (k por dinodo, n dinodos), tipico 6 a 10

asi que preguntar "¿D es lineal en V?" esta mal planteado. La pregunta util es
cuanto vale b, y si los datos alcanzan para RECHAZAR la recta.

Cuando los datos NO alcanzan (y hay que decirlo)
------------------------------------------------
En un rango angosto toda funcion suave parece una recta. El desvio maximo entre
la ley de potencia y la cuerda que une sus extremos vale, a segundo orden,

    desvio_relativo ~ (1/8) * b*(b-1) * (DV / V_medio)^2

o sea que con b = 8:

    rango de HV      span    desvio recta vs potencia
    1600-1700         6.1%        2.6%   <- del orden del sistematico: no decide
    1500-1700        12.5%         11%
    1450-1700        15.9%         18%   <- decide con holgura

Por eso `ajustar` devuelve `discriminacion`, y el informe dice explicitamente
cuando el span es insuficiente. Dibujar una recta por tres puntos juntos y
declarar linealidad seria el error clasico.

El test sin modelo: el COLAPSO de los espectros
-----------------------------------------------
Ajustar D en cada voltaje ya supone el modelo. Hay un test mas honesto y que no
depende de ninguna plantilla: dividir el eje de canal de cada corrida por su
D(V) y superponer. Si toda la cadena es lineal en amplitud, los espectros
COLAPSAN sobre una sola curva, porque lo unico que cambio es un factor de
escala. Si el colapso falla por arriba, hay saturacion; si falla abajo, es el
umbral del trigger, que no escala con la ganancia salvo que se lo haya movido a
proposito.

Correr `python analisis/voltaje.py` hace el cierre con verdad conocida (genera
D(V) de una ley de potencia elegida, con ruido, y verifica que se recupera b) y
despues, si el registro tiene dos o mas voltajes de la misma fuente, ajusta los
datos reales. Deja `figuras/09_voltaje.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, registro

# Rango de b fisicamente razonable para un PMT: k ~ 0.7-0.8 por dinodo y 8-12
# dinodos. Fuera de esto, o no es un PMT, o entre las corridas cambio algo mas
# que el voltaje.
B_PMT = (6.0, 10.0)


# --------------------------------------------------------------------------
# 1. Los dos ajustes
# --------------------------------------------------------------------------
def _recta_ponderada(x, y, err):
    """y = a + m*x por minimos cuadrados ponderados. Devuelve (a, m, ea, em)."""
    w = 1.0 / np.maximum(err, 1e-12) ** 2
    S, Sx = np.sum(w), np.sum(w * x)
    Sxx, Sy, Sxy = np.sum(w * x * x), np.sum(w * y), np.sum(w * x * y)
    det = S * Sxx - Sx * Sx
    if det <= 0:
        return np.nan, np.nan, np.nan, np.nan
    a = (Sxx * Sy - Sx * Sxy) / det
    m = (S * Sxy - Sx * Sy) / det
    return float(a), float(m), float(np.sqrt(Sxx / det)), float(np.sqrt(S / det))


def discriminacion(V, b=8.0):
    """
    Cuanto se separan una recta y la ley de potencia en el rango medido.

    Es el desvio relativo de la cuerda respecto de la curva en el punto medio,
    `(1/8)*b*(b-1)*(DV/V_medio)^2`. Compararlo con el error de D dice si el
    experimento puede distinguir las dos hipotesis ANTES de mirar el chi2.
    """
    V = np.asarray(V, float)
    if V.size < 2:
        return 0.0
    dV = V.max() - V.min()
    return float(0.125 * b * (b - 1.0) * (dV / np.mean(V)) ** 2)


def ajustar(V, D, err=None):
    """
    Ajusta D(V) con recta y con ley de potencia, y compara.

    La ley de potencia se ajusta en log-log, donde es una recta: log D =
    log D0 + b*log(V/V0). El pivote V0 es la media geometrica de los voltajes,
    no un numero redondo: centrar el pivote descorrelaciona D0 de b, asi que el
    error de b sale limpio en vez de estar dominado por la palanca del pivote.
    """
    V = np.asarray(V, float)
    D = np.asarray(D, float)
    err = np.full_like(D, 0.01 * np.mean(D)) if err is None else np.asarray(err, float)
    n = V.size
    if n < 2:
        raise ValueError("hacen falta al menos dos voltajes; hay %d" % n)

    # --- recta ---
    a, m, ea, em = _recta_ponderada(V, D, err)
    res_recta = D - (a + m * V)
    chi2_recta = float(np.sum((res_recta / err) ** 2))

    # --- ley de potencia, en log-log ---
    V0 = float(np.exp(np.mean(np.log(V))))
    x, y = np.log(V / V0), np.log(D)
    ey = err / D                                  # propagacion a log
    lD0, b, elD0, eb = _recta_ponderada(x, y, ey)
    D0 = float(np.exp(lD0))
    pred_pot = D0 * (V / V0) ** b
    res_pot = D - pred_pot
    chi2_pot = float(np.sum((res_pot / err) ** 2))

    gdl = max(n - 2, 0)
    sep = discriminacion(V, b if np.isfinite(b) and b > 1 else 8.0)
    err_rel = float(np.mean(err / D))
    return dict(
        V=V, D=D, err=err, n=n, gdl=gdl,
        span=float((V.max() - V.min()) / np.mean(V)),
        # recta
        ordenada=a, pendiente=m, err_ordenada=ea, err_pendiente=em,
        residuos_recta=res_recta, chi2_recta=chi2_recta,
        # potencia
        V0=V0, D0=D0, err_D0=float(D0 * elD0), b=float(b), err_b=float(eb),
        residuos_potencia=res_pot, chi2_potencia=chi2_pot,
        # veredicto
        separacion=sep, err_rel=err_rel,
        # con gdl = 0 los dos modelos pasan exacto por los puntos y comparar
        # chi2 no significa nada; hace falta un voltaje mas que parametros
        decide=bool(gdl > 0 and sep > 3.0 * err_rel),
        b_plausible=bool(np.isfinite(b) and B_PMT[0] <= b <= B_PMT[1]),
    )


def evaluar(res, V):
    """D predicha por la ley de potencia ajustada, en voltajes cualesquiera."""
    return res["D0"] * (np.asarray(V, float) / res["V0"]) ** res["b"]


# --------------------------------------------------------------------------
# 2. El test sin modelo: colapso de los espectros
# --------------------------------------------------------------------------
def colapso(espectros, Ds, e_lo=0.25, e_hi=1.15, n=200, E_ref=1.0):
    """
    Superpone los espectros con el eje dividido por su D, y mide la dispersion.

    espectros : lista de dicts con `canales` y `neto` (lo que deja `espectro.py`)
    Ds        : la D ajustada de cada uno
    e_lo, e_hi: la ventana de comparacion, en unidades de E_ref (el borde)

    Cada espectro se lleva a energia (E = ch/D), se interpola a una grilla comun
    y se normaliza por su integral EN LA VENTANA -no por el total, que depende
    del umbral de trigger y de cuantos eventos se guardaron-. Devuelve la
    dispersion relativa punto a punto: si la cadena es lineal en amplitud tiene
    que ser compatible con el ruido estadistico.
    """
    if len(espectros) != len(Ds):
        raise ValueError("un D por espectro")
    Eg = np.linspace(e_lo * E_ref, e_hi * E_ref, int(n))
    curvas = []
    for e, D in zip(espectros, Ds):
        E = np.asarray(e["canales"], float) / D
        y = np.asarray(e["neto"], float) * D          # jacobiano: cuentas por keV
        c = np.interp(Eg, E, y, left=np.nan, right=np.nan)
        area = np.nansum(c)
        curvas.append(c / area if area > 0 else c * np.nan)
    C = np.array(curvas)
    med = np.nanmean(C, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        disp = np.nanstd(C, axis=0) / med
    return dict(E=Eg, curvas=C, media=med, dispersion=disp,
                dispersion_media=float(np.nanmean(disp)),
                dispersion_max=float(np.nanmax(disp)) if np.isfinite(disp).any()
                else np.nan)


# --------------------------------------------------------------------------
# 3. Informe
# --------------------------------------------------------------------------
def informe(res, etiquetas=None):
    etiquetas = etiquetas or ["%.0f V" % v for v in res["V"]]
    print("  D(V) con %d voltajes, span del %.1f%% (%.0f a %.0f V)"
          % (res["n"], 100 * res["span"], res["V"].min(), res["V"].max()))
    print("     %-12s %10s %10s %12s %12s"
          % ("punto", "D", "error", "res. recta", "res. potencia"))
    for i, et in enumerate(etiquetas):
        print("     %-12s %10.4f %10.4f %12.4f %12.4f"
              % (et, res["D"][i], res["err"][i], res["residuos_recta"][i],
                 res["residuos_potencia"][i]))
    print("  RECTA        D = %.4g %+.5g*V" % (res["ordenada"], res["pendiente"]))
    print("               chi2 = %.2f con %d gdl" % (res["chi2_recta"], res["gdl"]))
    print("  POTENCIA     D = %.4f * (V/%.0f)^b   con  b = %.2f +- %.2f"
          % (res["D0"], res["V0"], res["b"], res["err_b"]))
    print("               chi2 = %.2f con %d gdl" % (res["chi2_potencia"], res["gdl"]))
    print("  ¿DECIDE?")
    print("     separacion esperada entre las dos hipotesis: %.2f%%"
          % (100 * res["separacion"]))
    print("     error relativo de D                        : %.2f%%"
          % (100 * res["err_rel"]))
    if res["gdl"] == 0:
        print("     -> NO. Con %d puntos y 2 parametros los dos modelos pasan "
              "exacto por los datos." % res["n"])
    elif not res["decide"]:
        print("     -> NO. El span es muy angosto: en este rango la potencia y "
              "la recta son indistinguibles.")
        print("        Hacen falta voltajes mas separados, no mas estadistica.")
    else:
        mejor = "potencia" if res["chi2_potencia"] < res["chi2_recta"] else "recta"
        print("     -> SI, y gana la %s (chi2 %.2f contra %.2f)."
              % (mejor, min(res["chi2_potencia"], res["chi2_recta"]),
                 max(res["chi2_potencia"], res["chi2_recta"])))
    if np.isfinite(res["b"]):
        print("  b = %.2f  ->  %s" % (
            res["b"],
            "dentro del rango tipico de un PMT (%.0f a %.0f)" % B_PMT
            if res["b_plausible"] else
            "FUERA del rango tipico de un PMT (%.0f a %.0f): entre las corridas "
            "cambio algo mas que el voltaje" % B_PMT))


# --------------------------------------------------------------------------
# 4. Los puntos que hay en el registro
# --------------------------------------------------------------------------
def puntos_del_registro(escala_pedida="keVee", metodo="global"):
    """
    Devuelve {(detector, fuente): dict} con los puntos D(V) ya calculados.

    Se toma UN punto por voltaje: si hay repeticiones del mismo punto se
    promedian, porque son la misma medicion hecha dos veces y meterlas como
    puntos independientes solo achica las barras sin agregar informacion.
    """
    salida = {}
    for llave, entradas in registro.por_fuente(None, escala_pedida, metodo).items():
        por_v = {}
        for _, r in entradas:
            por_v.setdefault(float(r["hv"]), []).append(r)
        V = np.array(sorted(por_v))
        D = np.array([np.mean([r["D"] for r in por_v[v]]) for v in V])
        err = np.array([np.mean([r["err_D"] for r in por_v[v]])
                        / np.sqrt(len(por_v[v])) for v in V])
        salida[llave] = dict(
            V=V, D=D, err=err,
            etiquetas=["%.0f V (n=%d)" % (v, len(por_v[v])) for v in V],
            claves=[r.get("corrida") for v in V for r in por_v[v][:1]])
    return salida


def avisos_de_cadena(claves):
    """
    Avisos de `rutas.misma_cadena` para las corridas que esten en el catalogo.

    Se filtra por las que el catalogo conoce porque una calibracion puede venir
    de archivos sueltos (la interfaz lo permite) y ahi no hay cadena registrada
    que comparar. Que falte no es motivo para no avisar de las que si estan.
    """
    corridas = [rutas.CORRIDAS[k] for k in claves if k in rutas.CORRIDAS]
    sueltas = [k for k in claves if k not in rutas.CORRIDAS]
    avisos = rutas.misma_cadena(corridas) if len(corridas) > 1 else []
    if sueltas:
        avisos.append("fuera del catalogo (cadena desconocida): %s"
                      % ", ".join(str(k) for k in sueltas))
    return avisos


# --------------------------------------------------------------------------
# demo: cierre con verdad conocida + lo que haya en el registro
# --------------------------------------------------------------------------
def _cierre(b_real=8.0, voltajes=(1450, 1500, 1550, 1600, 1650, 1700),
            err_rel=0.01, semilla=3):
    """Genera D(V) de una potencia conocida y verifica que se la recupera."""
    rng = np.random.default_rng(semilla)
    V = np.array(voltajes, float)
    D_real = 1.25 * (V / 1550.0) ** b_real
    err = err_rel * D_real
    D = D_real + rng.normal(0, err)
    res = ajustar(V, D, err)
    print("  verdad: b = %.2f    recuperado: b = %.2f +- %.2f  (%.1f sigma)"
          % (b_real, res["b"], res["err_b"],
             abs(res["b"] - b_real) / max(res["err_b"], 1e-9)))
    return res


def _cierre_colapso():
    """
    Verifica que el test sin modelo distingue lo que tiene que distinguir.

    Toma un espectro real, lo reescala a mano (la MISMA fisica con otra
    ganancia) y comprueba dos cosas: que con las D correctas el colapso es
    exacto, y que con una D corrida un 5% se nota. Sin la segunda parte, un
    colapso perfecto no probaria nada: una funcion que devuelve cero siempre
    tambien lo daria.
    """
    from analisis import espectro as _esp
    try:
        e = _esp.cargar("espectro_co60_global")
    except OSError:
        print("  (sin espectro de co60 cacheado, se saltea)")
        return
    D0 = 1.2512
    e2 = dict(canales=e["canales"] * 1.6, neto=e["neto"] / 1.6)
    bien = colapso([e, e2], [D0, D0 * 1.6], E_ref=963.4)
    mal = colapso([e, e2], [D0, D0 * 1.6 * 1.05], E_ref=963.4)
    print("  espectro contra su copia reescalada 1.6x:")
    print("     con las D correctas      : dispersion media %.1e, maxima %.1e"
          % (bien["dispersion_media"], bien["dispersion_max"]))
    print("     con la 2a D corrida un 5%%: dispersion media %.3f, maxima %.3f"
          % (mal["dispersion_media"], mal["dispersion_max"]))
    print("     -> el colapso es exacto cuando debe serlo y se rompe cuando debe")


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    print("=" * 74)
    print("1. CIERRE CON VERDAD CONOCIDA")
    print("=" * 74)
    res_ok = _cierre()
    informe(res_ok)
    print()
    _cierre_colapso()
    print()
    print("  el mismo experimento con solo los tres voltajes que hay hoy:")
    res_mal = _cierre(voltajes=(1600, 1650, 1700))
    informe(res_mal)

    print()
    print("=" * 74)
    print("2. LOS DATOS REALES DEL REGISTRO")
    print("=" * 74)
    reales = puntos_del_registro()
    if not reales:
        print("  no hay ninguna entrada del registro con HV anotado.")
        print("  (las entradas viejas no lo traen: hay que recalcularlas)")
    hechos = []
    for (det, fuente), d in sorted(reales.items()):
        print()
        print("  detector %s, fuente %s: %d voltaje(s)"
              % (det, fuente, d["V"].size))
        for a in avisos_de_cadena(d["claves"]):
            print("     AVISO: %s" % a)
        if d["V"].size < 2:
            print("     -> con un solo voltaje no hay curva D(V). Falta medir "
                  "la misma fuente a otros HV.")
            continue
        r = ajustar(d["V"], d["D"], d["err"])
        informe(r, d["etiquetas"])
        hechos.append(((det, fuente), r))

    fig, axes = plt.subplots(1, 2 + len(hechos), figsize=(6 * (2 + len(hechos)), 4.6),
                             squeeze=False)
    for ax, res, tit in ((axes[0][0], res_ok, "cierre: 1450-1700 V (span 16%)"),
                         (axes[0][1], res_mal, "cierre: 1600-1700 V (span 6%)")):
        _dibujar(ax, res, tit)
    for k, ((det, fuente), r) in enumerate(hechos):
        _dibujar(axes[0][2 + k], r, "%s con %s" % (det, fuente))
    fig.tight_layout()
    salida = rutas.figura("09_voltaje.png")
    fig.savefig(salida, dpi=130)
    print()
    print("figura -> %s" % rutas.relativa(salida))


def _dibujar(ax, res, titulo):
    V = res["V"]
    Vf = np.linspace(V.min() * 0.99, V.max() * 1.01, 200)
    ax.errorbar(V, res["D"], yerr=res["err"], fmt="o", color="k", ms=5,
                capsize=3, label="medido", zorder=4)
    ax.plot(Vf, res["ordenada"] + res["pendiente"] * Vf, ls="--", color="C0",
            label="recta")
    ax.plot(Vf, evaluar(res, Vf), color="C3",
            label="potencia, b = %.2f ± %.2f" % (res["b"], res["err_b"]))
    ax.set_title("%s\n%s" % (titulo,
                             "el span DECIDE" if res["decide"]
                             else "el span NO alcanza para decidir"), fontsize=10)
    ax.set_xlabel("alto voltaje [V]"); ax.set_ylabel("D [mV·ns/keVee]")
    ax.legend(fontsize=8); ax.grid(alpha=.3)


if __name__ == "__main__":
    _demo()
