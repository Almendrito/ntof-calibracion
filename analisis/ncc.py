"""
ncc.py — seleccion AUTOMATICA del intervalo del borde por correlacion cruzada
normalizada (Seccion II-B del paper, Ec. 15).

El problema que resuelve
------------------------
Para ajustar el borde Compton hay que decir en que tramo del espectro ajustar. A
mano eso es subjetivo y el resultado depende de la eleccion. El NCC lo decide
correlacionando la FORMA de los datos contra la plantilla teorica: se prueba un
tramo, se estandarizan las dos curvas y se calcula el coeficiente de Pearson.
Como Pearson resta la media y divide por la desviacion, es inmune a offset y a
escala: compara forma y nada mas.

Un solo parametro libre, no dos
-------------------------------
El paper barre (inicio, fin) de la ventana por separado. Eso sobra respecto de su
propio modelo: la Ec. 8 usa un mapa canal-energia POR EL ORIGEN (area cero <->
energia cero), asi que fijado el canal del borde B la ventana ya esta
determinada, `[inicio_rel*B, fin_rel*B]`. Barrer el ancho aparte agrega un grado
de libertad que el modelo no tiene, y sobre datos reales se lo gasta en elegir
ventanas degeneradas en la cola vacia (donde el z-score amplifica ruido puro
hasta parecerse a cualquier cosa). Aca se barre solo B: 1 dimension.

De B sale la calibracion: `D = B / E_ref`, con `E_ref` el borde de la primera
linea gamma. Esa D es la semilla del ajuste que viene despues.

Correr `python analisis/ncc.py` deja `figuras/05_ncc.png`.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, fisica, simulacion, espectro


# --------------------------------------------------------------------------
# 1. Pearson (Ec. 15)
# --------------------------------------------------------------------------
def pearson(a, b):
    """Coeficiente de correlacion de Pearson entre dos vectores del mismo largo."""
    a = np.asarray(a, float) - np.mean(a)
    b = np.asarray(b, float) - np.mean(b)
    d = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / d) if d > 0 else 0.0


def suavizar(y, ventana=21, orden=3):
    """Savitzky-Golay: quita ruido sin achatar el flanco como haria un promedio movil."""
    from scipy.signal import savgol_filter
    ventana = min(ventana if ventana % 2 else ventana + 1, len(y) - (1 - len(y) % 2))
    if ventana <= orden:
        return np.asarray(y, float)
    return savgol_filter(np.asarray(y, float), ventana, orden)


# --------------------------------------------------------------------------
# 2. La plantilla llevada a canal
# --------------------------------------------------------------------------
def _gaussiana(n_sigma_bins, radio_max, n_sigmas=6.0):
    """
    Kernel gaussiano normalizado, de radio `n_sigmas` sigmas.

    El radio es un numero ENTERO de bines, asi que crece a saltos cuando sigma
    crece: el modelo es discontinuo en sigma. Lo que decide si eso importa es
    cuanta masa del kernel se agrega en el salto. Con radio 3 sigma la cola que
    entra y sale pesa ~0.4% del pico, y al renormalizar mueve el modelo entero:
    medido sobre el Co-60, el chi2 pega un salto de 0.067 justo donde el radio
    pasa de 83 a 84 bines, o sea 29 veces el paso suave. Eso rompe el hessiano
    numerico -de ahi salen las barras de error y la correlacion D-sigma-, que
    con paso 1e-4 puede caer de un lado o del otro del salto.

    Con 6 sigmas la cola pesa exp(-18) ~ 1e-8 y el salto desaparece. El costo es
    el doble de largo el kernel, que es despreciable.
    """
    r = int(min(max(n_sigmas * n_sigma_bins, 1), radio_max))
    x = np.arange(-r, r + 1, dtype=float)
    k = np.exp(-0.5 * (x / max(n_sigma_bins, 1e-9)) ** 2)
    return k / k.sum()


def plantilla_en_canal(canales, cen_E, pesos, D, sigma_rel, E_ref, c0=0.0):
    """
    Lleva la plantilla de energia a espacio de canal y le aplica la resolucion.

    canal = D * energia (mapa por el origen), asi que la plantilla evaluada en
    `canales/D` es la misma curva reexpresada. El kernel gaussiano esta
    NORMALIZADO a suma 1: conserva el area, o sea que ensanchar no cambia las
    cuentas totales y la amplitud queda desacoplada de la resolucion.

    `sigma_rel` es el ancho relativo en `E_ref`, asi que en canal vale
    `sigma_rel * D * E_ref = sigma_rel * B`.

    `c0` es una ORDENADA en el mapa: `canal = D*E + c0`. El default 0 es la
    hipotesis del paper (area cero <-> energia cero) y es lo que se usa para
    calibrar. Liberarlo sirve para TESTEARLA con una sola fuente: el continuo
    Compton va de 0 a E_max, asi que su forma tiene informacion sobre c0
    aunque no haya un segundo borde. Es un test debil -c0 se confunde con el
    exceso de bajo canal y con sigma- y da una cota, no una medicion.
    """
    canales = np.asarray(canales, float)
    y = np.interp((canales - c0) / D, cen_E, pesos, left=0.0, right=0.0)
    paso = float(np.median(np.diff(canales)))
    sigma_ch = sigma_rel * D * E_ref
    n_sig = sigma_ch / max(paso, 1e-9)
    if n_sig < 0.5:
        return y
    k = _gaussiana(n_sig, radio_max=(len(y) - 1) // 2)
    return np.convolve(y, k, mode="same")


# --------------------------------------------------------------------------
# 3. El barrido
# --------------------------------------------------------------------------
def buscar_intervalo(canales, cuentas, cen_E, pesos, bordes,
                     sigma_rel=0.09, inicio_rel=0.45, fin_rel=None,
                     B_rango=None, n_B=240, suave=21, callback=None,
                     n_pts=None):
    """
    Barre el canal del borde B y devuelve la ventana que mejor correlaciona.

    bordes : bordes Compton de la fuente [keV], en orden. `E_ref = bordes[0]`
        es el ancla del mapa por el origen.
    inicio_rel : borde izquierdo de la ventana en unidades de B (0.45 es el del
        paper). Sube esto si el espectro tiene contaminacion a bajo canal.
    fin_rel : borde derecho en unidades de B. Por defecto se deriva de la fuente:
        tiene que pasar el ULTIMO borde, o la ventana corta justo la parte que
        distingue una plantilla de otra.
    B_rango : (min, max) del barrido. Por defecto, de 0.4 a 1.1 veces el canal
        donde el espectro suavizado cae a la mitad de su plateau.
    n_pts : a cuantos puntos se remuestrean las dos curvas. Por defecto, un
        punto por bin de la ventana media (no se inventan puntos). Se puede
        forzar para MEDIR el efecto: remuestrear a menos puntos promedia ruido
        e infla rho, y a mas no agrega informacion.

    Devuelve un dict con lo, hi, rho, B, D, sigma_rel y las curvas del barrido.
    """
    canales = np.asarray(canales, float)
    E_ref = float(bordes[0])
    if fin_rel is None:
        fin_rel = bordes[-1] / E_ref + 3.0 * sigma_rel

    y = suavizar(cuentas, suave)
    if B_rango is None:
        B0 = _borde_aproximado(canales, y)
        B_rango = (0.40 * B0, 1.10 * B0)

    Bs = np.linspace(B_rango[0], B_rango[1], n_B)
    # n FIJO durante todo el barrido: si cambiara con B, los rho de ventanas
    # distintas dejarian de ser comparables (la interpolacion a pocos puntos
    # infla rho promediando ruido).
    paso = float(np.median(np.diff(canales)))
    if n_pts is None:
        n_pts = max(int((fin_rel - inicio_rel) * np.mean(Bs) / paso), 24)
    n_pts = int(n_pts)

    rhos = np.full(n_B, -1.0)
    for i, B in enumerate(Bs):
        if callback is not None and i % 8 == 0:
            callback(i / n_B, "B = %.0f mV*ns" % B)
        lo, hi = inicio_rel * B, fin_rel * B
        m = (canales >= lo) & (canales <= hi)
        if m.sum() < 12:
            continue
        D = B / E_ref
        p = plantilla_en_canal(canales, cen_E, pesos, D, sigma_rel, E_ref)
        u = np.linspace(lo, hi, n_pts)
        s1 = np.interp(u, canales[m], y[m])
        s2 = np.interp(u, canales[m], p[m])
        rhos[i] = pearson(s1, s2)

    j = int(np.argmax(rhos))
    B = float(Bs[j])
    return dict(B=B, D=B / E_ref, rho=float(rhos[j]),
                lo=inicio_rel * B, hi=fin_rel * B,
                sigma_rel=sigma_rel, inicio_rel=inicio_rel, fin_rel=fin_rel,
                E_ref=E_ref, Bs=Bs, rhos=rhos, n_pts=n_pts)


def _borde_aproximado(canales, y):
    """Canal donde el espectro cae a la mitad del plateau. Solo para acotar el barrido."""
    pos = np.maximum(y, 0.0)
    cum = np.cumsum(pos) / max(pos.sum(), 1e-9)
    reg = (cum > 0.5) & (cum < 0.9)
    if reg.sum() < 5:
        return float(canales[len(canales) // 2])
    alto = float(np.median(y[reg]))
    i0 = int(np.flatnonzero(reg)[-1])
    j = np.flatnonzero((y < 0.5 * alto) & (np.arange(y.size) > i0))
    return float(canales[j[0]]) if j.size else float(canales[i0])


def barrer_sigma(canales, cuentas, cen_E, pesos, bordes, sigmas, **kw):
    """
    Repite el barrido para varios `sigma_rel`. El paper no dice con que ancho
    construir la plantilla (sigma es justo lo que el ajuste estima despues), asi
    que conviene ver cuanto se mueve B: si se mueve poco, la eleccion no importa.
    """
    return [buscar_intervalo(canales, cuentas, cen_E, pesos, bordes,
                             sigma_rel=s, **kw) for s in sigmas]


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES
              if os.path.exists(rutas.cache("espectro_%s_global.npz" % k))]
    if not claves:
        print("Faltan espectros. Corre primero:")
        print("   python analisis/hacer_espectro.py --corrida co60")
        return

    fig, axes = plt.subplots(len(claves), 3, figsize=(16.5, 4.4 * len(claves)),
                             squeeze=False)
    sigmas = [0.05, 0.07, 0.09, 0.11, 0.13]

    for fila, clave in enumerate(claves):
        c = rutas.corrida(clave)
        e = espectro.cargar("espectro_%s_global" % clave)
        bordes = fisica.bordes(c.fuente)
        print("=" * 74)
        print("%s  -  %s   (bordes teoricos: %s keV)"
              % (clave, c.fuente, ", ".join("%.1f" % b for b in bordes)))
        print("=" * 74)

        if not simulacion.coincide_fuente(c.sim, c.fuente):
            print("  la simulacion no es de esta fuente; salteo")
            continue
        cen, pesos = simulacion.cargar(c.sim, "Lvis_keVee")

        r = buscar_intervalo(e["canales"], e["neto"], cen, pesos, bordes)
        print("  ventana elegida: [%.0f, %.0f] mV*ns    rho = %.5f"
              % (r["lo"], r["hi"], r["rho"]))
        print("  B = %.1f mV*ns  <->  E_ref = %.1f keVee   ->   D = %.4f mV*ns/keVee"
              % (r["B"], r["E_ref"], r["D"]))
        print("  (%d puntos de interpolacion, inicio_rel=%.2f, fin_rel=%.2f)"
              % (r["n_pts"], r["inicio_rel"], r["fin_rel"]))

        print("  sensibilidad al sigma con que se construye la plantilla:")
        rs = barrer_sigma(e["canales"], e["neto"], cen, pesos, bordes, sigmas)
        for s, rr in zip(sigmas, rs):
            print("     sigma_rel = %.2f  ->  B = %7.1f   D = %.4f   rho = %.5f"
                  % (s, rr["B"], rr["D"], rr["rho"]))
        Bs_s = np.array([rr["B"] for rr in rs])
        print("     B varia %.2f%% en todo el rango; el rho maximo esta en "
              "sigma_rel = %.2f"
              % (100 * (Bs_s.max() - Bs_s.min()) / Bs_s.mean(),
                 sigmas[int(np.argmax([rr["rho"] for rr in rs]))]))
        print()

        # ---- panel 1: rho contra B
        ax = axes[fila][0]
        for s, rr in zip(sigmas, rs):
            ax.plot(rr["Bs"], rr["rhos"], lw=1, label="sigma_rel %.2f" % s)
        ax.axvline(r["B"], color="k", ls="--", lw=1, label="B = %.0f" % r["B"])
        ax.set_title("%s: barrido de 1 parametro" % clave)
        ax.set_xlabel("canal del borde B [mV*ns]"); ax.set_ylabel("rho (Pearson)")
        ax.set_ylim(max(0.0, np.nanmin([rr["rhos"].max() - 0.35 for rr in rs])), 1.005)
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        # ---- panel 2: la ventana sobre el espectro
        ax = axes[fila][1]
        ax.step(e["canales"], e["neto"], where="mid", lw=1, color="0.55",
                label="neto")
        p = plantilla_en_canal(e["canales"], cen, pesos, r["D"], r["sigma_rel"],
                               r["E_ref"])
        m = (e["canales"] >= r["lo"]) & (e["canales"] <= r["hi"])
        esc = np.sum(e["neto"][m]) / max(np.sum(p[m]), 1e-9)
        ax.plot(e["canales"], p * esc, lw=1.3, color="C1", label="plantilla Geant4")
        ax.axvspan(r["lo"], r["hi"], color="C0", alpha=0.13,
                   label="ventana automatica")
        ax.axvline(r["B"], color="C3", ls="--", lw=1.2,
                   label="B = %.0f mV*ns" % r["B"])
        ax.set_xlim(0, min(e["canales"][-1], 1.9 * r["B"]))
        ax.set_title("intervalo elegido (rho = %.4f)" % r["rho"])
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

        # ---- panel 3: lo que el NCC compara de verdad
        ax = axes[fila][2]
        u = np.linspace(r["lo"], r["hi"], r["n_pts"])
        s1 = np.interp(u, e["canales"], suavizar(e["neto"]))
        s2 = np.interp(u, e["canales"], p)
        z = lambda v: (v - v.mean()) / v.std()
        ax.plot(u, z(s1), lw=1.4, label="datos (S1)")
        ax.plot(u, z(s2), lw=1.4, label="plantilla (S2)")
        ax.set_title("lo que Pearson compara: forma estandarizada")
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("z-score")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

    fig.tight_layout()
    salida = rutas.figura("05_ncc.png")
    fig.savefig(salida, dpi=130)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
