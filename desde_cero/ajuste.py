"""
Etapa 7 — el ajuste: chi2, recocido simulado y la barra de error -> D.

    python desde_cero/ajuste.py

Que resuelve esta etapa
-----------------------
Con la ventana y la semilla del NCC (Etapa 6), encuentra los (D, sigma_rel)
que mejor explican el espectro neto (Etapa 3) con el modelo directo (Etapa
5), y dice cuanto se puede confiar en D. Es el ultimo paso del paper: con
esto la replica sobre nuestros datos esta completa, con plantilla
Klein-Nishina.

Las decisiones, y por que
-------------------------
1. **chi2 con pesos, no la Ec. 14.** El paper minimiza la suma de cuadrados
   SIN pesos. Eso le da a cada bin el mismo voto, y un bin del plateau con
   400 cuentas (error 20) vota igual que uno de la cola con 4 (error 2):
   el ajuste se deja llevar por donde hay mas cuentas y no hay manera de
   decir que tan bueno es. Pesando por 1/error^2 cada residuo se mide en
   unidades de SU incertidumbre, chi2/ndf ~ 1 dice "el modelo describe el
   dato dentro del ruido", y la curvatura del chi2 da la barra de error.
   El error es el de la Etapa 3, calculado con las cuentas OBSERVADAS: es el
   chi2 de NEYMAN. Tiene un sesgo conocido hacia abajo (un bin que fluctuo
   bajo recibe error chico y pesa de mas, §7.31) que va como 1/N por bin;
   el cierre de `verificar_ajuste.py` lo mide a la estadistica real.

2. **A se despeja, no se busca.** El modelo es lineal en A: m = A * f(D, s).
   Para (D, s) fijos, el A que minimiza el chi2 tiene forma cerrada,
       A = sum(y f / e^2) / sum(f^2 / e^2).
   Asi el recocido busca en 2 dimensiones y no en 3, y la mas dificil de las
   tres (A abarca ordenes de magnitud) desaparece. Es la idea de "VARPRO".
   El paper mete A en el recocido; §7.21 mostro que mezclar escalas tan
   distintas en un mismo paso rompe la busqueda.

3. **Recocido simulado (Metropolis).** A temperatura T, un paso al azar se
   acepta siempre si baja el chi2, y con probabilidad exp(-dchi2/T) si lo
   sube. A T alta la busqueda pasea por todo el espacio y puede salir de
   valles falsos; T baja 5 % por ciclo (alfa = 0.95, 100 propuestas por
   ciclo, como el paper) y la busqueda se va quedando en el valle mas hondo.
   El paso es U(-0.1, 0.1) x el ancho de la cota, como el paper.

4. **T_init se MIDE.** El paper dice 2000. Pero T se compara con dchi2, asi
   que 2000 significa cosas distintas segun las unidades del costo (§7.22:
   su calibracion automatica pedia 2.1e9 con la Ec. 14 y 3.2 con el error
   NCC). Se prueban 60 pasos desde el arranque y se elige T para que una
   subida tipica se acepte el 80 % de las veces: T = mediana(subidas)/ln(1/0.8).
   El corte tambien es relativo: el chi2 medio del ciclo cambia menos de 1e-6
   relativo durante 10 ciclos seguidos.

5. **8 arranques al azar dentro de las cotas, no desde la semilla.** El
   recocido es estocastico: una corrida no es un resultado. Y si los 8
   arrancaran de D0, que terminen juntos no probaria nada. Arrancando de
   puntos uniformes en D en [0.6, 1.6] D0 y sigma_rel en [0.02, 0.30], que
   terminen en el mismo minimo SI prueba que el minimo es el global.

6. **Pulido con Nelder-Mead.** El recocido encuentra el valle; no lo baja
   hasta el fondo con precision (su ultimo paso es de ~0.1 % de la cota).
   Un simplex local desde el mejor punto termina el trabajo.

7. **Barra de error por la curvatura.** Cerca del minimo el chi2 es una
   parabola: chi2(x) ~ chi2_min + (x - x0)^T H (x - x0) / 2, con H el
   hessiano (derivadas segundas, por diferencias finitas centradas). Para
   un chi2 gaussiano la covarianza de los parametros es C = 2 H^-1 (el 2
   porque chi2 = -2 ln L). Como A esta despejada en cada punto, esto es la
   curvatura del chi2 PERFILADO en A: la barra de D ya incluye no saber A.
   Si chi2/ndf > 1 el modelo no describe el dato al nivel del ruido (en los
   datos reales, seguro: Klein-Nishina no tiene scattering multiple), y la
   barra se multiplica por sqrt(chi2/ndf): es suponer que el desajuste es
   ruido de mas, repartido parejo. Es una convencion, no una ley, y se
   reporta la barra con y sin escalar.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import espectro, modelo, ncc                # noqa: E402

ALFA = 0.95
PROPUESTAS = 100
PASO = 0.1
ACEPTACION_INICIAL = 0.8
TOL_REL = 1e-6
CICLOS_ESTABLES = 10
MAX_CICLOS = 400
N_ARRANQUES = 8
SEMILLA = 7
PASO_HESSIANO = 1e-3     # relativo a cada parametro


class Problema:
    """El chi2 de UNA ventana de un espectro, con A despejada."""

    def __init__(self, canales, cuentas, errores, fuente, lo, hi, plantilla=None):
        m = (canales >= lo) & (canales <= hi)
        self.x, self.y, self.e = canales[m], cuentas[m], errores[m]
        self.w = 1.0 / self.e ** 2
        self.fuente = fuente
        self.plantilla = plantilla
        self.evaluaciones = 0

    def forma(self, D, s):
        return modelo.modelo(self.x, 1.0, D, s, self.fuente,
                             plantilla=self.plantilla)

    def amplitud(self, f):
        den = np.sum(self.w * f * f)
        return max(float(np.sum(self.w * self.y * f) / den), 0.0) if den > 0 else 0.0

    def chi2(self, p):
        D, s = p
        self.evaluaciones += 1
        if D <= 0 or s <= 0:
            return 1e30
        f = self.forma(D, s)
        A = self.amplitud(f)
        return float(np.sum(self.w * (self.y - A * f) ** 2))

    def prediccion(self, p):
        f = self.forma(*p)
        return self.amplitud(f) * f


# --------------------------------------------------------------------------
# Recocido
# --------------------------------------------------------------------------
def temperatura_inicial(costo, x0, lo, hi, rng, n=60):
    """T para aceptar el 80 % de las subidas tipicas desde x0."""
    d = hi - lo
    c0 = costo(x0)
    subidas = []
    for _ in range(n):
        c = costo(np.clip(x0 + rng.uniform(-PASO, PASO, d.size) * d, lo, hi)) - c0
        if c > 0:
            subidas.append(c)
    return float(np.median(subidas) / np.log(1.0 / ACEPTACION_INICIAL)) if subidas else 1.0


def recocido(costo, x0, lo, hi, semilla):
    """
    Un recocido. Devuelve el mejor punto visto y la historia por ciclo
    (T, chi2 medio, mejor chi2, aceptacion, punto actual).
    """
    rng = np.random.default_rng(semilla)
    d = hi - lo
    x = np.clip(np.asarray(x0, float), lo, hi)
    c = costo(x)
    mejor_x, mejor_c = x.copy(), c
    T = temperatura_inicial(costo, x, lo, hi, rng)
    T0 = T
    hist, previa, estables = [], None, 0
    for ciclo in range(MAX_CICLOS):
        suma, acept = 0.0, 0
        for _ in range(PROPUESTAS):
            p = np.clip(x + rng.uniform(-PASO, PASO, d.size) * d, lo, hi)
            cp = costo(p)
            if cp <= c or rng.random() < np.exp(-(cp - c) / T):
                x, c = p, cp
                acept += 1
                if c < mejor_c:
                    mejor_x, mejor_c = x.copy(), c
            suma += c
        media = suma / PROPUESTAS
        hist.append((T, media, mejor_c, acept / PROPUESTAS, x[0], x[1]))
        if previa is not None and abs(media - previa) <= TOL_REL * abs(previa):
            estables += 1
            if estables >= CICLOS_ESTABLES:
                break
        else:
            estables = 0
        previa = media
        T *= ALFA
    return dict(x=mejor_x, chi2=mejor_c, T0=T0, ciclos=len(hist),
                historia=np.array(hist))


def pulir(costo, x):
    from scipy.optimize import minimize
    r = minimize(costo, x, method="Nelder-Mead",
                 options=dict(xatol=1e-7, fatol=1e-7, maxiter=4000))
    return r.x if r.fun <= costo(x) else x


def hessiano(costo, x, paso_rel=PASO_HESSIANO):
    """Derivadas segundas por diferencias finitas centradas."""
    x = np.asarray(x, float)
    h = paso_rel * np.abs(x)
    n = x.size
    H = np.empty((n, n))
    f0 = costo(x)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h[i]
        H[i, i] = (costo(x + ei) - 2 * f0 + costo(x - ei)) / h[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = h[j]
            H[i, j] = H[j, i] = (costo(x + ei + ej) - costo(x + ei - ej)
                                 - costo(x - ei + ej) + costo(x - ei - ej)) \
                / (4 * h[i] * h[j])
    return H


# --------------------------------------------------------------------------
# La cadena
# --------------------------------------------------------------------------
def calibrar(canales, cuentas, errores, fuente, ventana=None, semilla=SEMILLA,
             n_arranques=N_ARRANQUES, plantilla=None):
    """
    NCC -> recocido (8 arranques) -> pulido -> barra de error.

    `ventana`: el dict de `ncc.barrido_B`, o None para calcularlo aca (con
    la MISMA plantilla que el ajuste).
    `plantilla`: None = Klein-Nishina; o la de Geant4 (`plantilla.py`).

    La ventana nunca arranca por debajo del fin de la subida del umbral
    (`espectro.fin_del_umbral`, §7.48): ahi el espectro es la eficiencia del
    trigger y ninguna plantilla la contiene. En las corridas de trabajo el
    NCC ya la deja arriba, asi que hoy no restringe nada.
    Devuelve un dict con D, err_D (escalado), err_D_crudo, sigma_rel, A,
    chi2, ndf, la covarianza, la ventana, y los recocidos (para dibujarlos).
    """
    x, y, e = (np.asarray(v, float) for v in (canales, cuentas, errores))
    w = ncc.barrido_B(x, y, fuente, plantilla=plantilla) if ventana is None         else ventana
    piso = espectro.fin_del_umbral(x, y)
    w = dict(w, lo=max(w["lo"], piso), piso=piso)
    P = Problema(x, y, e, fuente, w["lo"], w["hi"], plantilla)
    lo = np.array([w["cota_D"][0], w["cota_sigma"][0]])
    hi = np.array([w["cota_D"][1], w["cota_sigma"][1]])

    rng = np.random.default_rng(semilla)
    arranques = lo + rng.random((n_arranques, 2)) * (hi - lo)
    corridas = [recocido(P.chi2, a, lo, hi, semilla + k)
                for k, a in enumerate(arranques)]
    mejor = min(corridas, key=lambda r: r["chi2"])
    p = pulir(P.chi2, mejor["x"])
    chi2 = P.chi2(p)
    ndf = P.x.size - 3                       # A, D y sigma_rel
    H = hessiano(P.chi2, p)
    C = 2.0 * np.linalg.inv(H)
    escala = np.sqrt(max(chi2 / ndf, 1.0))
    err = np.sqrt(np.diag(C))
    return dict(D=float(p[0]), sigma_rel=float(p[1]),
                A=P.amplitud(P.forma(*p)), chi2=chi2, ndf=int(ndf),
                err_D_crudo=float(err[0]), err_D=float(err[0] * escala),
                err_sigma=float(err[1] * escala), escala=float(escala),
                corr=float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])), cov=C,
                ventana=w, problema=P, arranques=arranques, corridas=corridas,
                # cada arranque PULIDO: si todos caen en el mismo punto, el
                # minimo es unico en las cotas. Sin pulir, la dispersion es la
                # del ultimo paso del recocido (~0.1 % de la cota), no la de
                # minimos distintos.
                D_arranques=np.array([pulir(P.chi2, r["x"])[0] for r in corridas]),
                D_recocido=np.array([r["x"][0] for r in corridas]),
                evaluaciones=P.evaluaciones, p=p)


def perfil_D(res, n=21, ancho=3.0):
    """
    chi2 perfilado en D (sigma y A re-optimizados en cada D) alrededor del
    minimo. Donde sube 1 esta la barra SIN escalar: la comprobacion de que
    la curvatura no miente (§7.49).
    """
    from scipy.optimize import minimize_scalar
    P, p = res["problema"], res["p"]
    Ds = p[0] + np.linspace(-ancho, ancho, n) * res["err_D_crudo"]
    ch = []
    for D in Ds:
        r = minimize_scalar(lambda s: P.chi2((D, s)),
                            bounds=(0.5 * p[1], 1.5 * p[1]), method="bounded",
                            options=dict(xatol=1e-9))
        ch.append(r.fun)
    return Ds, np.array(ch) - res["chi2"]


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def _demo():
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from desde_cero import espectro
    from desde_cero.lector import CORRIDAS_TRABAJO

    fuentes = {"cs137_1600": "Cs-137", "co60": "Co-60", "a_co60_1600": "Co-60"}
    claves = [k for k in fuentes if os.path.exists(CORRIDAS_TRABAJO[k][0])]
    fig, ejes = plt.subplots(4, len(claves), figsize=(6 * len(claves), 15),
                             squeeze=False)
    for col, clave in enumerate(claves):
        e, _, _ = espectro.de_corrida(clave)
        t0 = time.time()
        r = calibrar(e["canales"], e["neto"], e["error"], fuentes[clave])
        P, w = r["problema"], r["ventana"]
        print("  %-12s D = %.4f +- %.4f (sin escalar %.4f) mV*ns/keV   "
              "sigma_rel = %.4f   chi2/ndf = %.1f/%d = %.2f   corr(D, s) = %+.2f"
              % (clave, r["D"], r["err_D"], r["err_D_crudo"], r["sigma_rel"],
                 r["chi2"], r["ndf"], r["chi2"] / r["ndf"], r["corr"]))
        print("  %-12s semilla NCC D0 = %.4f (%+.2f %%); 8 arranques terminan en "
              "D = %.4f a %.4f; %d evaluaciones, %.0f s"
              % ("", w["D0"], 100 * (w["D0"] / r["D"] - 1), r["D_arranques"].min(),
                 r["D_arranques"].max(), r["evaluaciones"], time.time() - t0))

        ax = ejes[0, col]
        x = e["canales"]
        m = (x > 0.2 * w["lo"]) & (x < 1.4 * w["hi"])
        ax.errorbar(x[m], e["neto"][m], e["error"][m], fmt=".", ms=2, lw=0.5,
                    color="0.5", label="neto")
        ax.plot(P.x, P.prediccion(r["p"]), color="tab:red", lw=1.5,
                label="ajuste Klein-Nishina")
        ax.axvline(r["D"] * ncc.bordes_de(fuentes[clave])[0], color="tab:red",
                   ls=":", lw=1)
        ax.set(title="%s: D = %.4f +- %.4f, chi2/ndf = %.2f"
               % (clave, r["D"], r["err_D"], r["chi2"] / r["ndf"]),
               xlabel="area [mV*ns]", ylabel="cuentas")
        ax.legend(fontsize=7)

        ax = ejes[1, col]
        ax.plot(P.x, (P.y - P.prediccion(r["p"])) / P.e, ".", ms=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set(xlabel="area [mV*ns]", ylabel="(dato - modelo) / error",
               title="residuos normalizados")

        ax = ejes[2, col]
        for c in r["corridas"]:
            ax.loglog(np.arange(1, len(c["historia"]) + 1),
                      c["historia"][:, 1] / r["ndf"], lw=1)
        ax.set(xlabel="ciclo de temperatura", ylabel="chi2 medio del ciclo / ndf",
               title="8 recocidos desde arranques al azar")

        ax = ejes[3, col]
        for c, a in zip(r["corridas"], r["arranques"]):
            ax.plot(c["historia"][:, 4], c["historia"][:, 5], lw=0.6, alpha=0.8)
            ax.plot(*a, "o", ms=3, color="k")
        ax.plot(*r["p"], "*", ms=12, color="tab:red")
        ax.plot(w["D0"], w["sigma_rel"], "s", ms=6, color="tab:blue",
                label="semilla NCC")
        ax.set(xlabel="D [mV*ns/keV]", ylabel="sigma_rel",
               title="camino en el plano (D, sigma_rel)")
        ax.legend(fontsize=7)
    for ax in ejes.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa7_ajuste.png")
    fig.savefig(salida, dpi=100)
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
