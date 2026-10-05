"""
_verif_area_plantilla.py — valida el TERCER metodo de area (ajuste de una
plantilla del pulso, `metodo="plantilla"`) contra los dos que ya existian.

Cinco partes, en el orden en que conviene creerlas (CLAUDE.md 7.74):

  1. FORMA. La plantilla supone que todos los pulsos tienen la misma forma.
     Antes de usarla se mide la forma de las 28 corridas y se la compara
     contra la plantilla (residuo rms, carga hasta 12/40/75 ns) y por cuartil
     de amplitud.                                   -> plantilla_pulso.png
  2. VERDAD CONOCIDA. Pulsos fabricados con la plantilla x amplitud conocida,
     ruido con el ESPECTRO medido en la linea base real, jitter del onset y
     registros truncados a 47, 81 y 158 ns. Sesgo de los tres metodos contra
     amplitud y contra truncamiento. Dos variantes que no son el caso ideal:
     con la componente aditiva medida (interferencia sincronizada con el
     trigger) y con una forma AJENA (pulsos de junio ajustados con la
     plantilla de ago-sep).                          -> plantilla_sintetico.png
  3. LINEALIDAD EN DATOS REALES. area/amplitud contra amplitud para los tres
     metodos y el cociente area_plantilla/area_compuerta evento a evento, en
     cs137_1600 (XML), a_co60_1600 (binario 1-sep) y d_cs137_1600 (8-sep, el
     registro mas corto).                            -> plantilla_linealidad.png
  4. CALIBRACION. La cadena completa (NCC -> SA -> curvatura) con los tres
     metodos en esas tres corridas.                  -> plantilla_espectros.png
  5. TESTS. b del test B y cocientes Cs/Co del test A con las 28 corridas
     recalibradas por plantilla (`calibrar_todo.py --area plantilla`) contra
     los vigentes.

La tabla resumen queda en `analisis/figuras/validacion_area_plantilla.txt`.

    python verificaciones/_verif_area_plantilla.py            # todo
    python verificaciones/_verif_area_plantilla.py forma sintetico
"""

import os as _os
import sys as _sys
_RAIZ = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _RAIZ not in _sys.path:
    _sys.path.insert(0, _RAIZ)

import json
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analisis import rutas, drs, pulso, area, pipeline

# La misma configuracion que el lote (`calibrar_todo.py`), para que los numeros
# de la parte 4 sean comparables con las tablas de los tests.
from analisis import calibrar_todo as ct

CORRIDAS_REALES = ("cs137_1600", "a_co60_1600", "d_cs137_1600")
METODOS = ("global", "convergencia", "plantilla")
COLOR = {"global": "C0", "convergencia": "C3", "plantilla": "C2"}
NOMBRE = {"global": "compuerta 12 ns", "convergencia": "convergencia",
          "plantilla": "plantilla"}

_TXT = []


def P(s=""):
    print(s)
    _TXT.append(s)


def _Q(tau, f, x):
    """Carga acumulada de la forma hasta `x` ns."""
    Q = pulso.integral_acumulada(tau, f)
    return float(Q[min(int(np.searchsorted(tau, x)), Q.size - 1)])


# ==========================================================================
# 1. FORMA
# ==========================================================================
def parte_forma():
    P("=" * 78)
    P("1. LA FORMA DEL PULSO: una plantilla por montaje, verificada corrida a corrida")
    P("=" * 78)
    pls = {n: area.cargar_plantilla(n, verbose=True) for n in area.PLANTILLAS}
    for n, p in pls.items():
        P("  plantilla '%s' (%s): %d pulsos de > %.0f mV, grilla hasta %.0f ns"
          % (n, ", ".join(area.PLANTILLAS[n]), p["n_pulsos"], pulso.AMP_PLANTILLA,
             p["tau"][-1]))
        P("     carga hasta T_CARGA = %.0f ns: %.3f ns x amplitud   (hasta el "
          "maximo de Q: %.3f, en t = %.1f ns; %s)"
          % (p["t_carga"], p["integral"], p["integral_plateau"], p["t_plateau"],
             "converge dentro de la grilla" if p["convergio"]
             else "NO converge dentro de la grilla"))
        P("     fraccion de la carga (hasta %.0f ns) recogida hasta 12 / 40 / 47 "
          "/ 81 ns: %s" % (p["t_carga"], " / ".join(
              "%.3f" % (_Q(p["tau"], p["forma"], x) / p["integral"])
              for x in (12, 40, 47, 81))))
        P("     componente ADITIVA (regresion d = alfa*A + beta): integral de "
          "beta hasta 12 / 40 / %.0f ns = %s mV*ns"
          % (p["t_carga"], " / ".join("%.2f" % _Q(p["tau"], p["aditiva"], x)
                                      for x in (12, 40, p["t_carga"]))))

    # --- las 28 corridas contra su plantilla ------------------------------
    P()
    P("  cada corrida contra la plantilla que le toca (forma normalizada; rms "
      "del residuo en unidades de la amplitud,")
    P("  sobre [-10 ns, fin del registro o 120 ns]; Q_X = carga hasta X ns del onset)")
    P("  %-16s %-8s %5s %8s %8s %8s %8s %9s" % ("corrida", "plantil.", "n",
                                                "Q12/Q40", "QX/Q40", "X[ns]",
                                                "rms", "rms otra"))
    filas = {}
    for k in sorted(rutas.CORRIDAS):
        c = rutas.corrida(k)
        if not c.disponible():
            continue
        nom = area.plantilla_para(c.muestra)
        otra = "junio" if nom == "ago-sep" else "ago-sep"
        p = pulso.plantilla_pulso([c.muestra], n=800, t_carga="fin")
        X = min(75.0, p["tau"][-1] - 1.0)
        m = (p["tau"] >= -10) & (p["tau"] <= min(p["tau"][-1], 120.0))
        rms = {}
        for n_ in (nom, otra):
            r = np.interp(p["tau"][m], pls[n_]["tau"], pls[n_]["forma"])
            rms[n_] = float(np.sqrt(np.mean((p["forma"][m] - r) ** 2)))
        q12 = _Q(p["tau"], p["forma"], 12) / _Q(p["tau"], p["forma"], 40)
        qx = _Q(p["tau"], p["forma"], X) / _Q(p["tau"], p["forma"], 40)
        ref = pls[nom]
        qx_ref = _Q(ref["tau"], ref["forma"], X) / _Q(ref["tau"], ref["forma"], 40)
        filas[k] = dict(p=p, nom=nom, q12=q12, qx=qx, X=X, qx_ref=qx_ref,
                        rms=rms[nom], rms_otra=rms[otra])
        P("  %-16s %-8s %5d %8.4f %8.4f %8.0f %8.4f %9.4f%s"
          % (k, nom, p["n_pulsos"], q12, qx, X, rms[nom], rms[otra],
             "   (descartada)" if rutas.descartada(k) else ""))
    buenas = [f for k, f in filas.items() if f["nom"] == "ago-sep"
              and not rutas.descartada(k)]
    rq = np.array([f["q12"] for f in buenas])
    P("  ago-sep (sin descartadas): Q12/Q40 = %.4f +- %.4f (rango %.4f-%.4f); "
      "rms mediano %.4f"
      % (rq.mean(), rq.std(), rq.min(), rq.max(),
         np.median([f["rms"] for f in buenas])))

    # --- por cuartil de amplitud -------------------------------------------
    P()
    P("  la forma POR CUARTIL DE AMPLITUD (pulso promedio alineado; lo que "
      "importa es si Q_X/Q40 depende de la amplitud)")
    cuart = {}
    for k in ("cs137_1600", "a_co60_1700"):
        T, V, _, _ = drs.leer_pulsos(rutas.corrida(k).muestra, n=4000)
        amps = np.array([pulso.detectar_pico(*pulso.preparar(t, v))[2]
                         for t, v in zip(T, V)])
        q = np.percentile(amps[amps > 8], [0, 25, 50, 75, 100])
        cuart[k] = []
        for j in range(4):
            mm = (amps >= q[j]) & (amps <= q[j + 1])
            g, f, n = pulso.pulso_promedio(T[mm], V[mm], amp_minima=0.0)
            cuart[k].append((q[j], q[j + 1], g, f))
            P("     %-12s %4.0f-%4.0f mV  Q12/Q40 %.4f   Q75/Q40 %.4f"
              % (k, q[j], q[j + 1], _Q(g, f, 12) / _Q(g, f, 40),
                 _Q(g, f, 75) / _Q(g, f, 40)))

    # --- figura -------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    a = ax[0][0]
    for n, p in pls.items():
        a.plot(p["tau"], p["forma"], lw=1.3, label="plantilla '%s'" % n)
    a.set_xlabel("t - onset [ns]"); a.set_ylabel("forma (maximo = 1)")
    a.set_title("(a) las dos plantillas")
    a.legend(fontsize=8); a.grid(alpha=0.3)
    ins = a.inset_axes([0.42, 0.35, 0.55, 0.55])
    for n, p in pls.items():
        ins.plot(p["tau"], p["forma"], lw=1.0)
    ins.set_xlim(20, 150); ins.set_ylim(-0.03, 0.10); ins.grid(alpha=0.3)
    ins.set_title("cola (x10)", fontsize=7); ins.tick_params(labelsize=7)

    a = ax[0][1]
    for n, p in pls.items():
        Q = pulso.integral_acumulada(p["tau"], p["forma"]) / p["integral"]
        a.plot(p["tau"], 100 * Q, lw=1.4, label="'%s'" % n)
    for x, et in ((12, "compuerta"), (47, "fin reg. D 8-sep"),
                  (81, "fin reg. 1-sep"), (pulso.T_CARGA, "T_CARGA")):
        a.axvline(x, color="0.5", ls=":", lw=1)
        a.text(x + 1, 52, et, rotation=90, fontsize=7, color="0.35")
    a.set_ylim(50, 104); a.set_xlabel("t - onset [ns]")
    a.set_ylabel("carga acumulada [% de la carga hasta T_CARGA]")
    a.set_title("(b) convergencia: lo que el registro contiene y lo que se extrapola")
    a.legend(fontsize=8); a.grid(alpha=0.3)

    a = ax[1][0]
    cmap = {"A": "C0", "B": "C1", "C": "C2", "D": "C4"}
    for k, f in filas.items():
        c = rutas.corrida(k)
        ref = pls[f["nom"]]
        r = f["p"]["forma"] - np.interp(f["p"]["tau"], ref["tau"], ref["forma"])
        if rutas.descartada(k) or f["nom"] == "junio":
            continue
        a.plot(f["p"]["tau"], r, lw=0.6, alpha=0.7, color=cmap.get(c.detector, "k"))
    for d_, col in cmap.items():
        a.plot([], [], color=col, label="detector %s" % d_)
    a.axhline(0, color="k", lw=0.6)
    a.set_xlim(-15, 125); a.set_ylim(-0.06, 0.06)
    a.set_xlabel("t - onset [ns]"); a.set_ylabel("forma - plantilla 'ago-sep'")
    a.set_title("(c) residuo de las corridas de ago/sep (24, sin descartadas)")
    a.legend(fontsize=8, ncol=4); a.grid(alpha=0.3)

    a = ax[1][1]
    for k, st in (("cs137_1600", "-"), ("a_co60_1700", "--")):
        for j, (q0, q1, g, f) in enumerate(cuart[k]):
            Q = pulso.integral_acumulada(g, f)
            a.plot(g, Q / _Q(g, f, 40), ls=st, lw=1.1, color="C%d" % j,
                   label="%s %0.f-%0.f mV" % (k, q0, q1))
    a.set_xlim(0, 125); a.set_ylim(0.6, 1.12)
    a.set_xlabel("t - onset [ns]"); a.set_ylabel("Q(t) / Q(40 ns)")
    a.set_title("(d) forma por cuartil de amplitud")
    a.legend(fontsize=6.5, ncol=2); a.grid(alpha=0.3)
    fig.tight_layout()
    ruta = rutas.figura("plantilla_pulso.png")
    fig.savefig(ruta, dpi=120)
    plt.close(fig)
    P("  figura -> %s" % rutas.relativa(ruta))
    return pls, filas


# ==========================================================================
# 2. VERDAD CONOCIDA
# ==========================================================================
class Sintetico:
    """
    Fabrica eventos con verdad conocida.

    El RUIDO no es blanco: se genera con el espectro de potencia medido en la
    linea base real del detector D del 8-sep (sus registros tienen ~150 ns de
    pre-onset limpio). El ruido de un DRS4 tiene correlacion entre celdas, y un
    ruido blanco del mismo rms sobreestimaria cuanto promedia la compuerta.
    El eje de TIEMPO es el de eventos reales (no uniforme, como el del DRS4).
    """

    def __init__(self, semilla=1):
        self.rng = np.random.default_rng(semilla)
        T, V, _, _ = drs.leer_pulsos(rutas.corrida("d_co60_1600").muestra, n=400)
        seg = np.array([v[5:605] - np.median(v[5:605]) for v in V])
        self.psd = np.mean(np.abs(np.fft.rfft(seg, axis=1)) ** 2, axis=0)
        self.sigma = float(seg.std())
        self.ejes = [t - t[0] for t in T[:60]]

    def ruido(self, n):
        f = np.fft.rfftfreq(n)
        amp = np.sqrt(np.interp(f, np.fft.rfftfreq(600), self.psd))
        x = np.fft.irfft(amp * np.exp(2j * np.pi * self.rng.random(f.size)), n)
        return x / x.std() * self.sigma

    def evento(self, forma, amp, L, aditiva=None, jitter=1.0):
        """Un registro con el onset a L ns del final (L = registro tras el
        onset), polaridad negativa como la real, linea base al azar."""
        t = self.ejes[self.rng.integers(len(self.ejes))]
        ton = t[-1] - L + self.rng.normal(0.0, jitter)
        s = amp * np.interp(t - ton, forma["tau"], forma["forma"],
                            left=0.0, right=0.0)
        if aditiva is not None:
            s = s + np.interp(t - ton, aditiva["tau"], aditiva["aditiva"],
                              left=0.0, right=0.0)
        v = -(s + self.ruido(t.size)) + self.rng.normal(0.0, 3.0)
        return t, v


def _tres_metodos(t, v, pl, pre, post):
    r = pulso.area_evento_plantilla(t, v, pl, pre, post)
    c = pulso.area_evento(t, v)
    return dict(global_=r["area_compuerta"], convergencia=c["area"],
                plantilla=r["area"], amplitud=r["amplitud"])


def parte_sintetico(pls, n_ev=250):
    P()
    P("=" * 78)
    P("2. VERDAD CONOCIDA")
    P("=" * 78)
    s = Sintetico()
    pl = pls["ago-sep"]
    pre, post = 2.2, 12.0
    P("  ruido con el espectro de la linea base real: rms %.2f mV; jitter del "
      "onset 1 ns; compuerta [onset-%.1f, onset+%.0f] ns" % (s.sigma, pre, post))
    P("  verdad = amplitud x (carga de la plantilla hasta %.0f ns); se reporta "
      "area/verdad (mediana +- error de la mediana)" % pl["t_carga"])
    P("  la compuerta mide una FRACCION fija por construccion (~0.78): lo que "
      "importa es que sea PLANA, no que valga 1.")
    amps = [5, 10, 20, 50, 100, 200, 400]
    Ls = [47, 81, 158]
    casos = [("ideal", pl, None), ("con aditiva", pl, pl),
             ("forma ajena (junio)", pls["junio"], None)]
    res = {}
    for nom, forma, adit in casos:
        P()
        P("  --- %s ---" % nom)
        P("  %4s %5s   %-17s %-17s %-17s" % ("L", "amp", "compuerta", "convergencia",
                                           "plantilla"))
        for L in Ls if nom == "ideal" else [81]:
            for a in amps:
                acc = {m: [] for m in METODOS}
                verdad = a * forma["integral"]
                for _ in range(n_ev):
                    t, v = s.evento(forma, a, L, aditiva=adit)
                    r = _tres_metodos(t, v, pl, pre, post)
                    acc["global"].append(r["global_"] / verdad)
                    acc["convergencia"].append(r["convergencia"] / verdad)
                    acc["plantilla"].append(r["plantilla"] / verdad)
                fila = {}
                for m in METODOS:
                    x = np.array(acc[m])
                    fila[m] = (float(np.median(x)),
                               float(1.2533 * x.std() / np.sqrt(x.size)),
                               float(x.std()))
                res[(nom, L, a)] = fila
                P("  %4d %5d   %s" % (L, a, "  ".join(
                    "%.4f+-%.4f" % fila[m][:2] for m in METODOS)))

    # resumen: deriva entre 20 y 400 mV y entre L = 47 y 158
    P()
    for m in METODOS:
        d_amp = res[("ideal", 81, 400)][m][0] - res[("ideal", 81, 20)][m][0]
        d_L = res[("ideal", 158, 100)][m][0] - res[("ideal", 47, 100)][m][0]
        ref = res[("ideal", 81, 400)][m][0]
        P("  %-13s deriva 20->400 mV: %+.2f%%   deriva L 47->158 ns (100 mV): "
          "%+.2f%%   resolucion a 100 mV: %.2f%%"
          % (NOMBRE[m], 100 * d_amp / ref, 100 * d_L / ref,
             100 * res[("ideal", 81, 100)][m][2] / ref))

    # figura
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    for m in METODOS:
        for L, ls in zip(Ls, (":", "--", "-")):
            y = [res[("ideal", L, a)][m][0] for a in amps]
            e = [res[("ideal", L, a)][m][1] for a in amps]
            ax[0].errorbar(amps, y, yerr=e, ls=ls, marker="o", ms=3,
                           color=COLOR[m],
                           label="%s, L=%d" % (NOMBRE[m], L) if L == 81 or m == "plantilla"
                           else None)
    ax[0].set_xscale("log"); ax[0].set_ylim(0.45, 1.08)
    ax[0].set_xlabel("amplitud [mV]"); ax[0].set_ylabel("area / verdad")
    ax[0].set_title("(a) sesgo contra amplitud (:47  --81  -158 ns)")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    for m in METODOS:
        for a, mk in ((20, "s"), (100, "o")):
            y = [res[("ideal", L, a)][m][0] / res[("ideal", 158, a)][m][0]
                 for L in Ls]
            ax[1].plot(Ls, y, marker=mk, color=COLOR[m],
                       label="%s, %d mV" % (NOMBRE[m], a))
    ax[1].set_xlabel("registro tras el onset L [ns]")
    ax[1].set_ylabel("area(L) / area(158 ns)")
    ax[1].set_title("(b) sesgo contra truncamiento")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    for nom, ls in (("ideal", "-"), ("con aditiva", "--"),
                    ("forma ajena (junio)", ":")):
        for m in ("global", "plantilla"):
            y = np.array([res[(nom, 81, a)][m][0] for a in amps])
            ax[2].plot(amps, y / y[-1], ls=ls, marker="o", ms=3, color=COLOR[m],
                       label="%s, %s" % (NOMBRE[m], nom))
    ax[2].set_xscale("log"); ax[2].set_ylim(0.9, 1.1)
    ax[2].set_xlabel("amplitud [mV]"); ax[2].set_ylabel("area / area(400 mV)")
    ax[2].set_title("(c) variantes: componente aditiva y forma ajena (L = 81)")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)
    fig.tight_layout()
    ruta = rutas.figura("plantilla_sintetico.png")
    fig.savefig(ruta, dpi=120)
    plt.close(fig)
    P("  figura -> %s" % rutas.relativa(ruta))
    return res


# ==========================================================================
# 3 y 4. DATOS REALES: linealidad y calibracion
# ==========================================================================
def _opciones(metodo):
    return dict(ct.OPCIONES, metodo=metodo, guardar=False)


def calibrar(clave, metodo):
    return pipeline.calibrar_completo(rutas.corrida(clave), _opciones(metodo))


def _pedestal_predicho(pl, pre, post=12.0):
    """
    Cuanto le suma la componente aditiva `beta` a cada metodo [mV*ns] y el `k`
    de la ley `cociente = R_inf - k/amplitud` que eso predice.

    Compuerta: integral de beta sobre la compuerta, con la base tomada antes
    del onset. Plantilla: la proyeccion de beta sobre la forma con b libre
    (lo que el ajuste lineal le atribuye a A).
    """
    g, f, b = pl["tau"], pl["forma"], pl["aditiva"]
    base = float(np.median(b[(g > -20) & (g < -4)]))
    m = (g >= -pre) & (g <= post)
    cg = float(pulso._TRAPZ(b[m] - base, g[m]))
    gA = float(pulso._TRAPZ(f[m], g[m]))
    X = np.vstack([np.ones(g.size), f]).T
    alfa = float(np.linalg.lstsq(X, b, rcond=None)[0][1])
    R = pl["integral"] / gA
    return dict(cg=cg, cp=alfa * pl["integral"], R=R,
                k=R * (cg / gA - alfa))


def parte_real(pls):
    P()
    P("=" * 78)
    P("3. LINEALIDAD EN DATOS REALES y 4. CALIBRACION (misma configuracion que "
      "el lote:")
    P("   escala %s, compuerta %g ns, %s)"
      % (ct.RAMA, ct.POST_FIJO, "con exceso" if ct.EXCESO else "sin exceso"))
    P("=" * 78)
    sal = {}
    for k in CORRIDAS_REALES:
        sal[k] = {}
        for m in METODOS:
            t0 = time.time()
            print("  [%s / %s] ..." % (k, m), flush=True)
            sal[k][m] = calibrar(k, m)
            print("     %.0f s" % (time.time() - t0), flush=True)

    # --- 3. linealidad ------------------------------------------------------
    fig, ax = plt.subplots(2, len(CORRIDAS_REALES), figsize=(17, 8.5),
                           squeeze=False)
    P()
    P("  3. area/amplitud contra amplitud (mediana por tramo, normalizada al "
      "tramo de mayor amplitud) y cociente plantilla/compuerta")
    for j, k in enumerate(CORRIDAS_REALES):
        d_p = sal[k]["plantilla"]["areas_muestra"]
        d_c = sal[k]["convergencia"]["areas_muestra"]
        ok = area.buenos(d_p) & area.buenos(d_c) & (d_p["amplitud"] > 5)
        amp = d_p["amplitud"][ok]
        P("  %s (%d eventos buenos, plantilla '%s', %.0f%% con la cola "
          "extrapolada, chi2_red mediano %.2f)"
          % (k, ok.sum(), d_p["plantilla"], 100 * d_p["extrapolada"][ok].mean(),
             np.median(d_p["chi2_red"][ok])))
        for m, areas in (("global", d_p["area_compuerta"][ok]),
                         ("convergencia", d_c["area"][ok]),
                         ("plantilla", d_p["area"][ok])):
            x, y, e = area.perfil(amp, areas / amp, bins=12)
            ax[0][j].errorbar(x, y / y[-1], yerr=e / y[-1], marker="o", ms=3,
                              color=COLOR[m], label=NOMBRE[m])
            sobre = x >= np.median(amp)
            P("     %-15s area/amp: de %.3f (%.0f mV) a %.3f (%.0f mV) -> "
              "%+.1f%%;  sobre la amplitud mediana %+.1f%%"
              % (NOMBRE[m], y[0], x[0], y[-1], x[-1],
                 100 * (y[-1] - y[0]) / y[-1],
                 100 * (y[sobre][-1] - y[sobre][0]) / y[-1]))
        ax[0][j].set_xscale("log"); ax[0][j].set_ylim(0.6, 1.15)
        ax[0][j].set_xlabel("amplitud [mV]")
        ax[0][j].set_ylabel("area/amplitud (norm. al ultimo tramo)")
        ax[0][j].set_title("%s: linealidad" % k)
        ax[0][j].legend(fontsize=7); ax[0][j].grid(alpha=0.3)

        # cociente evento a evento, y la ley R_inf - k/A
        r = d_p["area"][ok] / d_p["area_compuerta"][ok]
        m20 = amp > 15
        x, y, e = area.perfil(amp[m20], r[m20], bins=14)
        A = np.vstack([np.ones(x.size), 1.0 / x]).T
        (R_inf, mk), *_ = np.linalg.lstsq(A, y, rcond=None)
        pred = _pedestal_predicho(pls[d_p["plantilla"]],
                                  sal[k]["plantilla"]["ventana"]["pre"])
        xx = np.geomspace(x[0], x[-1], 50)
        ax[1][j].errorbar(x, y, yerr=e, marker="o", ms=3, color="k",
                          label="medido")
        ax[1][j].plot(xx, R_inf + mk / xx, color="C1",
                      label="ajuste R - k/A: k = %.2f mV" % -mk)
        ax[1][j].plot(xx, pred["R"] - pred["k"] / xx, color="C2", ls="--",
                      label="predicho por la aditiva: k = %.2f mV" % pred["k"])
        ax[1][j].set_xscale("log")
        ax[1][j].set_xlabel("amplitud [mV]")
        ax[1][j].set_ylabel("area plantilla / area compuerta")
        ax[1][j].set_title("%s: cociente evento a evento" % k)
        ax[1][j].legend(fontsize=7); ax[1][j].grid(alpha=0.3)
        P("     cociente plantilla/compuerta: %.4f (%.0f mV) -> %.4f (%.0f mV), "
          "deriva %+.2f%%; ajuste R - k/A: R = %.4f, k = %.2f mV"
          % (y[0], x[0], y[-1], x[-1], 100 * (y[-1] - y[0]) / y[-1], R_inf, -mk))
        P("        predicho por la componente aditiva de la plantilla: R = %.4f, "
          "k = %.2f mV (la compuerta recoge %+.2f mV*ns de la aditiva, la "
          "plantilla %+.2f)" % (pred["R"], pred["k"], pred["cg"], pred["cp"]))
        sobre = x >= np.median(amp)
        P("        sobre la amplitud mediana: deriva %+.2f%%"
          % (100 * (y[sobre][-1] - y[sobre][0]) / y[-1]))
    fig.tight_layout()
    ruta = rutas.figura("plantilla_linealidad.png")
    fig.savefig(ruta, dpi=120)
    plt.close(fig)
    P("  figura -> %s" % rutas.relativa(ruta))

    # --- 4. calibracion -----------------------------------------------------
    P()
    P("  4. calibracion con los tres metodos")
    P("  %-14s %-15s %10s %8s %7s %8s %8s %11s"
      % ("corrida", "metodo", "D", "err", "FWHM", "chi2/ndf", "pulsos",
         "D_12/D_pl"))
    fracciones = []
    for k in CORRIDAS_REALES:
        Dp = sal[k]["plantilla"]["resultado"]["D"]
        for m in METODOS:
            r = sal[k][m]["resultado"]
            P("  %-14s %-15s %10.4f %8.4f %6.1f%% %8.3f %8d %11s"
              % (k, NOMBRE[m], r["D"], r["err_D"], 100 * r["fwhm_rel"],
                 r["chi2_ndf"], sal[k][m]["espectro"]["n_muestra"],
                 "%.4f" % (r["D"] / Dp) if m == "global" else ""))
        rg = sal[k]["global"]["resultado"]
        f = rg["D"] / Dp
        ef = f * np.hypot(rg["err_D"] / rg["D"],
                          sal[k]["plantilla"]["resultado"]["err_D"] / Dp)
        fracciones.append((k, f, ef))
    fr = np.array([f for _, f, _ in fracciones])
    P("  D_12ns/D_plantilla: %s  ->  dispersion %.2f%%"
      % (", ".join("%.4f" % f for f in fr), 100 * fr.std() / fr.mean()))
    P("  (ojo: las dos D salen de ajustes INDEPENDIENTES del mismo espectro; su "
      "error no es independiente, asi que la")
    P("   barra de cada cociente es una cota superior)")

    # figura de espectros
    fig, ax = plt.subplots(1, len(CORRIDAS_REALES), figsize=(17, 4.8),
                           squeeze=False)
    for j, k in enumerate(CORRIDAS_REALES):
        a = ax[0][j]
        for m in METODOS:
            e = sal[k][m]["espectro"]
            r = sal[k][m]["resultado"]
            E = e["canales"] / r["D"]
            a.step(E, e["neto"] / e["neto"][(E > 0.3 * r["E_ref"])
                                            & (E < 0.6 * r["E_ref"])].mean(),
                   where="mid", lw=1.1, color=COLOR[m],
                   label="%s: D = %.3f, FWHM %.1f%%"
                   % (NOMBRE[m], r["D"], 100 * r["fwhm_rel"]))
        a.axvline(sal[k]["global"]["resultado"]["E_ref"], color="0.4", ls=":")
        a.set_xlim(0, 1.6 * sal[k]["global"]["resultado"]["E_ref"])
        a.set_ylim(0, None)
        a.set_xlabel("area / D  [keV]")
        a.set_ylabel("cuentas / plateau")
        a.set_title("%s: los tres espectros en su propia escala" % k)
        a.legend(fontsize=7); a.grid(alpha=0.3)
    fig.tight_layout()
    ruta = rutas.figura("plantilla_espectros.png")
    fig.savefig(ruta, dpi=120)
    plt.close(fig)
    P("  figura -> %s" % rutas.relativa(ruta))
    return sal


# ==========================================================================
# 5. TESTS
# ==========================================================================
def parte_tests():
    P()
    P("=" * 78)
    P("5. LOS TESTS con las 28 corridas recalibradas por plantilla")
    P("=" * 78)
    ruta_g = os.path.join(rutas.FIGURAS, "calibraciones", "resultados.json")
    ruta_p = os.path.join(rutas.FIGURAS, "calibraciones_plantilla",
                          "resultados.json")
    if not os.path.exists(ruta_p):
        P("  falta %s: correr `python analisis/calibrar_todo.py --area "
          "plantilla`" % rutas.relativa(ruta_p))
        return
    with open(ruta_g, encoding="utf-8") as f:
        rg = json.load(f)
    with open(ruta_p, encoding="utf-8") as f:
        rp = json.load(f)
    comunes = sorted(k for k in rg if not k.startswith("_") and k in rp)
    P("  %d corridas en las dos tablas" % len(comunes))
    P()
    P("  %-16s %9s %9s %9s %7s %7s %6s %6s %6s %6s"
      % ("corrida", "D_12ns", "D_pl", "D12/Dpl", "FWHM12", "FWHMpl", "chi2_12",
         "chi2pl", "extrap", "chi2ev"))
    razones, dfw = [], []
    for k in comunes:
        a, b = rg[k], rp[k]
        q = a["D"] / b["D"]
        if not a["descartada"]:
            razones.append((k, q, rutas.corrida(k).detector,
                            area.plantilla_para(rutas.corrida(k).muestra)))
            dfw.append(b["fwhm_rel"] / a["fwhm_rel"] - 1)
        P("  %-16s %9.4f %9.4f %9.4f %6.1f%% %6.1f%% %6.2f %6.2f %5.0f%% %6.2f%s"
          % (k, a["D"], b["D"], q, 100 * a["fwhm_rel"], 100 * b["fwhm_rel"],
             a["chi2_ndf"], b["chi2_ndf"], 100 * b.get("frac_extrapolada", np.nan),
             b.get("chi2_red_mediana", np.nan),
             "  (descartada)" if a["descartada"] else ""))
    q = np.array([r for _, r, _, n in razones if n == "ago-sep"])
    P()
    P("  D_12ns/D_plantilla, ago-sep (sin descartadas, %d corridas): %.4f +- %.4f "
      "(dispersion %.2f%%, rango %.4f-%.4f)"
      % (q.size, q.mean(), q.std(), 100 * q.std() / q.mean(), q.min(), q.max()))
    for det in "ABCD":
        qd = np.array([r for _, r, d_, n in razones if d_ == det and n == "ago-sep"])
        if qd.size:
            P("     detector %s: %.4f +- %.4f (%d)" % (det, qd.mean(), qd.std(), qd.size))
    qj = [r for k, r, _, n in razones if n == "junio"]
    if qj:
        P("  junio (co60, otra plantilla): %.4f" % qj[0])
    dfw = np.array(dfw)
    P("  FWHM plantilla / FWHM compuerta - 1: mediana %+.2f%%, rango %+.2f%% a "
      "%+.2f%%" % (100 * np.median(dfw), 100 * dfw.min(), 100 * dfw.max()))

    # test B y test A con las MISMAS funciones del lote
    tb_g, tb_p = ct.test_b(rg), ct.test_b(rp)
    P()
    P("  TEST B: exponente b de D ~ V^b")
    P("  %-18s %16s %16s %8s" % ("grupo", "compuerta 12 ns", "plantilla", "dif"))
    for g_ in sorted(tb_g):
        ag, ap = tb_g[g_]["ajuste"], tb_p.get(g_, {}).get("ajuste")
        if not ag or not ap:
            continue
        bg, bp = ag["b"], ap["b"]
        eg, ep = ag["err_b"], ap["err_b"]
        P("  %-18s %7.3f +- %5.3f %7.3f +- %5.3f %+8.3f"
          % ("%s %s" % g_, bg, eg, bp, ep, bp - bg))
    ta_g, ta_p = ct.test_a(rg), ct.test_a(rp)
    P()
    P("  TEST A: cociente D(Cs-137)/D(Co-60) en cada punto (esperado por Birks "
      "%.4f)" % ct.cociente_esperado())
    P("  %-12s %16s %16s %8s" % ("punto", "compuerta 12 ns", "plantilla", "dif"))
    for g_ in sorted(ta_g):
        cg_, cp_ = ta_g[g_]["cociente"], ta_p.get(g_, {}).get("cociente")
        if not cg_ or not cp_:
            continue
        P("  %-12s %8.4f+-%.4f %8.4f+-%.4f %+7.2f%%"
          % ("%s %g V" % g_, cg_["valor"], cg_["err"], cp_["valor"], cp_["err"],
             100 * (cp_["valor"] / cg_["valor"] - 1)))


# ==========================================================================
def main():
    rutas.salida_utf8()
    partes = [a for a in _sys.argv[1:]] or ["forma", "sintetico", "real",
                                            "tests"]
    t0 = time.time()
    pls = {n: area.cargar_plantilla(n) for n in area.PLANTILLAS}
    if "forma" in partes:
        parte_forma()
    if "sintetico" in partes:
        parte_sintetico(pls)
    if "real" in partes:
        parte_real(pls)
    if "tests" in partes:
        parte_tests()
    P()
    P("(%.1f min)" % ((time.time() - t0) / 60))
    nombre = "validacion_area_plantilla%s.txt" % (
        "" if len(partes) == 4 else "_" + "_".join(partes))
    with open(rutas.figura(nombre), "w", encoding="utf-8") as f:
        f.write("\n".join(_TXT) + "\n")
    print("tabla -> %s" % rutas.relativa(rutas.figura(nombre)))


os = _os
if __name__ == "__main__":
    main()
