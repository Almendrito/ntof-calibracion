"""
modelo.py — el modelo directo: de la plantilla teorica al espectro medido.

    S_modelo(ch) = A * [ plantilla(ch / D) (*) G(sigma) ]

con tres parametros:

    A        amplitud (cuentas). Solo escala: no cambia la forma.
    D        ganancia [mV*ns / keV]. Es LO QUE SE BUSCA. El mapa canal-energia
             pasa por el origen, ch = D * E, asi que D es el unico parametro que
             relaciona el eje medido con el eje fisico.
    sigma_rel  resolucion, como fraccion de B = D * E_ref. En canal el ancho
             gaussiano vale sigma_ch = sigma_rel * B.

Es la Ec. 8-9 del paper, con una diferencia: donde el paper pone la lineshape
analitica de Klein-Nishina (su Ec. 10), aca va el histograma de Geant4. La
plantilla simulada ya trae adentro el scattering multiple, el escape y el
quenching de Birks, y ademas trae las DOS lineas del Co-60 con su intensidad
relativa correcta. Por eso nuestro vector de parametros tiene 3 componentes y no
las 6 de la Ec. 11 del paper: A2/A1, B2/B1 y la mezcla de scattering no son
parametros libres, los fija la simulacion.

El costo es chi2 con las barras de error del espectro neto:

    chi2(theta) = suma_i [ (dato_i - modelo_i(theta)) / error_i ]^2

El exceso de bajo canal (opcional, `exceso=True`)
-------------------------------------------------
El espectro medido tiene, por debajo del borde, MAS cuentas de las que la
simulacion predice: 41% de mas a 0.35*B y creciendo hacia abajo. No es el umbral
de trigger (eso seria un DEFICIT) ni el fondo (ya esta restado); los candidatos
son dispersion en el material de alrededor y en el portafuente, que la geometria
simulada no tiene.

Ignorarlo no es gratis. Si la ventana de ajuste entra en esa zona, el ajuste lo
absorbe ensanchando sigma y BAJANDO D, y el sesgo es una deriva monotona con el
limite inferior de la ventana:

    inicio_rel   0.45     0.55     0.65     0.75     0.80
    D (Co-60)   1.2230   1.2302   1.2391   1.2422   1.2430
    chi2/ndf     1.576    1.372    1.165    1.090    1.067

O sea que "elegir bien la ventana" no es una eleccion neutral: es elegir cuanto
sesgo se acepta. Por eso el exceso se MODELA en vez de esquivarse,

    S(ch) = A * [plantilla(ch/D) (*) G(sigma)]  +  C * exp(-(ch - lo) / lambda)

con lambda como tercer parametro no lineal y (A, C) proyectadas por nnls, asi
que el recocido sigue buscando en pocas dimensiones. Con el exceso adentro se
puede ajustar MUCHO mas abajo (0.25*B en vez de 0.45*B), que es donde el borde
tiene su flanco completo.

La FORMA del exceso es una eleccion, no fisica medida: cambiarla (exponencial
contra ley de potencia) y ver cuanto se mueve D es el sistematico de esa
eleccion, y hoy es el mas grande del presupuesto. Cerrarlo pide una medicion:
la misma fuente alejada o apantallada.

Correr `python analisis/modelo.py` deja `figuras/06_modelo.png`.
"""

import os
import sys

import numpy as np
from scipy.optimize import nnls

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, fisica, simulacion, espectro, ncc


def prediccion(canales, cen_E, pesos, A, D, sigma_rel, E_ref):
    """Espectro predicho en las unidades de los datos (cuentas por bin)."""
    p = ncc.plantilla_en_canal(canales, cen_E, pesos, D, sigma_rel, E_ref)
    return A * p


def amplitud_optima(dato, error, forma):
    """
    A que minimiza el chi2, en forma cerrada.

    A entra LINEALMENTE en el modelo, asi que no hace falta buscarla: dado
    (D, sigma) la mejor A se despeja. Esto se llama proyeccion (VARPRO) y tiene
    dos consecuencias practicas grandes:
      - el optimizador busca en 2 dimensiones en vez de 3;
      - desaparece el problema de escala entre parametros (A vale ~1e3 y D ~1,
        y un paso comun para los dos no sirve para ninguno).
    """
    w = 1.0 / np.maximum(error, 1e-9) ** 2
    den = float(np.sum(w * forma * forma))
    if den <= 0:
        return 0.0
    return float(np.sum(w * forma * dato) / den)


def forma_exceso(canales, lam, lo=0.0, tipo="exp"):
    """
    La forma del exceso de bajo canal, con amplitud 1 en `lo`.

    Se normaliza en el borde izquierdo de la ventana y no en el canal 0 para que
    la columna valga O(1) donde hay datos: con exp(-ch/lambda) crudo y lambda de
    ~130 mV*ns, la columna vale 1e-2 en el arranque de la ventana y su amplitud
    sale 1e2 veces mas grande que la de la plantilla, lo que descondiciona el
    sistema lineal sin ninguna necesidad.

    `tipo` existe para poder MEDIR el sistematico de la eleccion: "exp" y "pot"
    describen los datos casi igual de bien y dan D distintas.
    """
    ch = np.asarray(canales, float)
    lam = max(float(lam), 1e-9)
    if tipo == "pot":
        # ley de potencia con el mismo parametro de escala: (lo/ch)^(lo/lambda)
        base = np.maximum(ch, 1e-9) / max(lo, 1e-9)
        return base ** (-max(lo, 1e-9) / lam)
    return np.exp(-(ch - lo) / lam)


def amplitudes_optimas(dato, error, columnas):
    """
    Las amplitudes lineales optimas de varias componentes, a la vez y positivas.

    Mismo argumento que `amplitud_optima` (VARPRO), extendido a mas de una
    columna: fijados los parametros no lineales, el problema es lineal y se
    resuelve exacto. Se usa nnls y no minimos cuadrados libres porque una
    amplitud negativa no significa nada fisico -ni la plantilla ni el exceso
    pueden restar cuentas- y dejarla libre le da al ajuste una via de escape
    para compensar errores de forma.
    """
    w = 1.0 / np.maximum(np.asarray(error, float), 1e-9)
    M = np.column_stack([np.asarray(c, float) * w for c in columnas])
    a, _ = nnls(M, np.asarray(dato, float) * w)
    return a


# Paso de derivacion del hessiano, relativo al valor de cada parametro (D,
# sigma_rel, lam). `c0` usa un bin (ver `Problema.pasos_curvatura`). Esta a
# nivel de modulo para que el informe lo cite sin tipearlo.
PASO_CURV_REL = 1e-3


class Problema:
    """
    Junta datos + plantilla + ventana, y expone el costo que el optimizador ve.

    Se arma una vez y se llama muchas veces: guardar aca la mascara de la
    ventana y los arrays recortados evita rehacer ese trabajo en cada una de las
    ~20 000 evaluaciones que pide un recocido.
    """

    def __init__(self, canales, cuentas, error, cen_E, pesos, E_ref, lo, hi,
                 exceso=False, tipo_exceso="exp", ordenada=False,
                 varianza="dato"):
        m = (canales >= lo) & (canales <= hi)
        self.canales_full = np.asarray(canales, float)
        self.m = m
        self.canales = self.canales_full[m]
        self.dato = np.asarray(cuentas, float)[m]
        self.error = np.maximum(np.asarray(error, float)[m], 1.0)
        self.cen_E, self.pesos, self.E_ref = cen_E, pesos, E_ref
        self.lo, self.hi = float(lo), float(hi)
        self.exceso = bool(exceso)
        self.tipo_exceso = tipo_exceso
        self.ordenada = bool(ordenada)
        # Con que varianza se pesa cada bin. "dato" es el chi2 de NEYMAN, el
        # de siempre: sigma_i^2 = la barra del espectro neto, que viene de las
        # cuentas OBSERVADAS. "modelo" es el chi2 de PEARSON: sigma_i^2 = la
        # prediccion, recalculada en cada evaluacion. La diferencia no es
        # cosmetica: con Neyman un bin que fluctuo hacia abajo recibe una barra
        # mas chica y por lo tanto MAS peso, y el ajuste se deja arrastrar hacia
        # abajo (Baker y Cousins, NIM 221, 1984). El sesgo va como 1/N por bin.
        # La campana se ajusta con "dato"; "modelo" existe para MEDIR ese sesgo
        # con verdad conocida (informe.cierre_conocido).
        if varianza not in ("dato", "modelo"):
            raise ValueError("varianza tiene que ser 'dato' o 'modelo'")
        self.varianza = varianza
        # Los parametros NO LINEALES que busca el optimizador, en orden. Las
        # amplitudes (A y, con exceso, C) no estan: se proyectan exacto en cada
        # evaluacion (VARPRO), asi que el recocido busca en 2, 3 o 4
        # dimensiones y no en 4, 5 o 6.
        self.libres = ["D", "sigma_rel"]
        if self.exceso:
            self.libres.append("lam")
        if self.ordenada:
            self.libres.append("c0")
        self.n_par = len(self.libres) + (2 if self.exceso else 1)
        self.ndf = int(m.sum()) - self.n_par
        self.C_exceso = 0.0          # ultima amplitud del exceso proyectada
        self.n_eval = 0
        self.paso_canal = float(np.median(np.diff(self.canales_full)))             if self.canales_full.size > 1 else 1.0

    def pasos_curvatura(self, x):
        """
        Paso de derivacion de cada parametro libre, para el hessiano.

        Casi todos escalan con su propio valor; `c0` no, porque su escala
        natural es el BIN del histograma (la plantilla se interpola sobre esa
        grilla). Ver `sa.errores_por_curvatura`.
        """
        x = np.asarray(x, float)
        h = np.maximum(np.abs(x) * PASO_CURV_REL, 1e-9)
        if self.ordenada:
            h[self.libres.index("c0")] = self.paso_canal
        return h

    def valores(self, x):
        """El vector del optimizador -> dict con TODOS los parametros.

        Los que no son libres quedan en su valor fijo: `lam` no se usa y `c0`
        vale 0, que es la hipotesis del mapa por el origen."""
        v = dict(D=np.nan, sigma_rel=np.nan, lam=None, c0=0.0)
        for nombre, valor in zip(self.libres, x):
            v[nombre] = float(valor)
        return v

    def forma(self, D, sigma_rel, c0=0.0):
        """La plantilla en canal, sin amplitud. Se evalua en TODA la grilla y
        despues se recorta: la convolucion necesita los canales de afuera de la
        ventana, porque el broadening mete cuentas hacia adentro desde abajo.
        Evaluar solo dentro de la ventana hace que el ajuste compense
        achicando sigma."""
        p = ncc.plantilla_en_canal(self.canales_full, self.cen_E, self.pesos,
                                   D, sigma_rel, self.E_ref, c0)
        return p[self.m]

    def exceso_col(self, lam):
        """La columna del exceso dentro de la ventana, con amplitud 1 en `lo`."""
        return forma_exceso(self.canales, lam, self.lo, self.tipo_exceso)

    def chi2(self, D, sigma_rel, lam=None, A=None, c0=0.0):
        """
        chi2 con las amplitudes proyectadas (si no se las dan).

        `lam` solo se usa cuando el problema se armo con `exceso=True`. La
        amplitud del exceso queda en `self.C_exceso`, no se devuelve, para que
        la firma `(chi2, A)` siga valiendo para todo lo que ya la usa.
        """
        self.n_eval += 1
        f = self.forma(D, sigma_rel, c0)
        if not self.exceso:
            if self.varianza == "modelo":
                return self._chi2_pearson(f, A)
            if A is None:
                A = amplitud_optima(self.dato, self.error, f)
            self.C_exceso = 0.0
            r = (self.dato - A * f) / self.error
            return float(np.sum(r * r)), float(A)

        g = self.exceso_col(lam)
        if A is None:
            a = amplitudes_optimas(self.dato, self.error, (f, g))
            A, C = float(a[0]), float(a[1])
        else:
            C = max(amplitud_optima(self.dato - A * f, self.error, g), 0.0)
        self.C_exceso = float(C)
        r = (self.dato - A * f - C * g) / self.error
        return float(np.sum(r * r)), float(A)

    def _chi2_pearson(self, f, A=None):
        """
        chi2 de Pearson: la varianza de cada bin es la PREDICCION, no el dato.

        Como la varianza depende de A y A se proyecta con esa varianza, se
        itera: se arranca con las barras del dato, se proyecta A, se
        recalcula la varianza con A*f y se vuelve a proyectar. Tres vueltas
        alcanzan (la tercera mueve A en menos de 1e-6 relativo). El piso de 1
        cuenta es el mismo que el de las barras del dato.
        """
        err = self.error
        for _ in range(3 if A is None else 1):
            A_ = amplitud_optima(self.dato, err, f) if A is None else A
            err = np.sqrt(np.maximum(A_ * f, 1.0))
        self.C_exceso = 0.0
        r = (self.dato - A_ * f) / err
        return float(np.sum(r * r)), float(A_)

    def chi2_de(self, x, A=None):
        """chi2 y amplitud para un vector del optimizador."""
        v = self.valores(x)
        return self.chi2(v["D"], v["sigma_rel"], v["lam"], A=A, c0=v["c0"])

    def prediccion(self, x, A=None):
        """El modelo completo dentro de la ventana, exceso incluido."""
        v = self.valores(x)
        _, A = self.chi2(v["D"], v["sigma_rel"], v["lam"], A=A, c0=v["c0"])
        f = self.forma(v["D"], v["sigma_rel"], v["c0"])
        if not self.exceso:
            return A * f
        return A * f + self.C_exceso * self.exceso_col(v["lam"])

    def costo(self, x):
        """Costo para el optimizador: x sigue el orden de `self.libres`."""
        v = self.valores(x)
        if not (v["D"] > 0 and v["sigma_rel"] > 0):
            return 1e30
        if self.exceso and not v["lam"] > 0:
            return 1e30
        return self.chi2(v["D"], v["sigma_rel"], v["lam"], c0=v["c0"])[0]

    def mapa(self, Ds, sigmas, x_base=None):
        """
        Matriz de chi2 sobre la grilla (D, sigma), para dibujar el paisaje.

        `x_base` es el vector COMPLETO de parametros libres del que se toman
        los que no son D ni sigma (`lam`, `c0`). Hace falta en cuanto el
        modelo tiene mas de dos: sin el, `valores` deja `lam` en None y el
        costo revienta con un TypeError. Eso es exactamente lo que pasaba con
        el modelo por defecto (que lleva el exceso adentro): la pestana del
        paisaje mostraba "no se pudo dibujar" y el error no llegaba a ninguna
        parte.

        Lo que se dibuja entonces es un CORTE del chi2: D y sigma barridos con
        el resto de los parametros congelados en el valor que se pase -el
        optimo, si se quiere ver el valle donde el ajuste se quedo-. No es el
        chi2 perfilado (eso pediria reoptimizar el resto en cada punto de la
        grilla, unas 10 000 optimizaciones), y por eso el contorno que se
        dibuje encima hay que leerlo como el del corte.
        """
        resto = ()
        if len(self.libres) > 2:
            if x_base is None:
                raise ValueError(
                    "este problema tiene %d parametros libres (%s): hay que "
                    "pasar x_base para decir en que valor se congelan los que "
                    "no se barren." % (len(self.libres), ", ".join(self.libres)))
            resto = tuple(float(v) for v in np.asarray(x_base, float)[2:])
        Z = np.empty((len(sigmas), len(Ds)))
        for i, s in enumerate(sigmas):
            for j, D in enumerate(Ds):
                Z[i, j] = self.costo((D, s) + resto)
        return Z


def armar(corrida, metodo="global", rama="Lvis_keVee", sigma_rel=0.09,
          inicio_rel=0.45, verbose=True, esp=None, plantilla=None,
          cb_ncc=None, exceso=False, tipo_exceso="exp",
          inicio_ajuste_rel=None, ordenada=False, varianza="dato"):
    """
    Carga espectro + simulacion, corre el NCC y devuelve `(problema, semilla)`.

    corrida : clave del catalogo o una `rutas.Corrida` armada al vuelo. Lo
        segundo es lo que permite calibrar un detector nuevo, que por definicion
        no esta en el catalogo.
    esp : espectro ya calculado. Si no se da, se lee del `.npz` por nombre. La
        interfaz lo tiene en memoria y volver a leerlo del disco no tendria
        sentido; ademas, para archivos sueltos ese `.npz` ni existe.
    plantilla : `(centros, pesos)` ya cargados de la simulacion, por lo mismo.
    exceso : agrega la componente de bajo canal al modelo (ver la cabecera).
    inicio_ajuste_rel : baja el arranque de la ventana de AJUSTE por debajo del
        que eligio el NCC, en unidades de B. Existe porque las dos ventanas
        responden a criterios distintos: el NCC compara FORMA contra una
        plantilla que no tiene el exceso, asi que conviene que mire donde el
        exceso ya murio; el ajuste, en cambio, si el exceso esta modelado, gana
        con el flanco entero. Sin `exceso=True` no tiene sentido usarlo.

    La semilla sale del NCC: D0 = B/E_ref y sigma0 = el sigma_rel del barrido.
    Es exactamente el eslabon que el paper describe: NCC elige el intervalo Y da
    el punto de partida del optimizador.
    """
    c = rutas.resolver(corrida)
    e = esp if esp is not None else espectro.cargar(
        "espectro_%s_%s" % (c.clave, metodo))
    bordes = fisica.bordes(c.fuente)
    cen, pesos = plantilla if plantilla is not None \
        else simulacion.cargar(c.sim, rama)
    r = ncc.buscar_intervalo(e["canales"], e["neto"], cen, pesos, bordes,
                             sigma_rel=sigma_rel, inicio_rel=inicio_rel,
                             callback=cb_ncc)
    lo = r["lo"]
    piso = espectro.fin_del_umbral(e["canales"], e["neto"])
    if inicio_ajuste_rel is not None:
        # B = D*E_ref es el canal del borde de referencia que estimo el NCC
        lo = min(lo, float(inicio_ajuste_rel) * r["D"] * bordes[0])
    # el umbral del trigger es un piso duro: ninguna ventana puede entrar ahi
    lo = max(lo, piso)
    p = Problema(e["canales"], e["neto"], e["error"], cen, pesos,
                 bordes[0], lo, r["hi"],
                 exceso=exceso, tipo_exceso=tipo_exceso, ordenada=ordenada,
                 varianza=varianza)
    if verbose:
        print("  NCC -> ventana [%.0f, %.0f] mV*ns, rho = %.5f, "
              "semilla D = %.4f, sigma_rel = %.3f"
              % (r["lo"], r["hi"], r["rho"], r["D"], r["sigma_rel"]))
        if lo < r["lo"]:
            print("  ventana de AJUSTE extendida hacia abajo: [%.0f, %.0f] "
                  "(el exceso esta modelado)" % (lo, r["hi"]))
        if lo <= piso:
            print("  el umbral de trigger MANDA sobre lo pedido: la ventana no "
                  "puede arrancar antes de %.0f mV*ns" % piso)
        print("  problema: %d bins, %d grados de libertad"
              % (p.m.sum(), p.ndf))
    return p, (r["D"], r["sigma_rel"]), r, e, c


def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES
              if os.path.exists(rutas.cache("espectro_%s_global.npz" % k))]
    fig, axes = plt.subplots(len(claves), 3, figsize=(16.5, 4.4 * len(claves)),
                             squeeze=False)

    for fila, clave in enumerate(claves):
        print("=" * 74)
        print(clave)
        print("=" * 74)
        # la demo muestra el modelo DESNUDO, que es el que se explica arriba;
        # la calibracion usa ademas el exceso de bajo canal (calibracion.EXCESO)
        p, semilla, r, e, c = armar(clave)

        # pulido local desde la semilla del NCC (Nelder-Mead: sin derivadas)
        res = minimize(p.costo, semilla, method="Nelder-Mead",
                       options=dict(xatol=1e-5, fatol=1e-3, maxiter=2000))
        D, s = res.x
        chi2, A = p.chi2(D, s)
        print("  semilla NCC : D = %.4f  sigma_rel = %.4f  ->  chi2/ndf = %.3f"
              % (semilla[0], semilla[1], p.chi2(*semilla)[0] / p.ndf))
        print("  pulido local: D = %.4f  sigma_rel = %.4f  ->  chi2/ndf = %.3f "
              "(A = %.3g, %d evaluaciones)"
              % (D, s, chi2 / p.ndf, A, p.n_eval))

        ax = axes[fila][0]
        ax.errorbar(p.canales, p.dato, yerr=p.error, fmt="none", ecolor="0.8", lw=.5)
        ax.step(p.canales, p.dato, where="mid", lw=1, color="0.4", label="medido")
        ax.plot(p.canales, A * p.forma(D, s), lw=1.6, color="C1", label="modelo")
        ax.plot(p.canales, p.chi2(*semilla)[1] * p.forma(*semilla), lw=1,
                ls="--", color="C0", label="semilla del NCC")
        ax.set_title("%s: modelo vs medido" % clave)
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

        ax = axes[fila][1]
        ax.axhline(0, color="k", lw=.8)
        ax.step(p.canales, (p.dato - A * p.forma(D, s)) / p.error, where="mid", lw=1)
        ax.set_title("residuos normalizados (chi2/ndf = %.3f)" % (chi2 / p.ndf))
        ax.set_xlabel("area [mV*ns]"); ax.set_ylabel("(dato - modelo) / error")
        ax.grid(alpha=.3)

        ax = axes[fila][2]
        Ds = np.linspace(0.80 * D, 1.20 * D, 90)
        chis = np.array([p.costo((d, s)) for d in Ds]) / p.ndf
        ax.plot(Ds, chis, lw=1.4)
        ax.axvline(D, color="C1", ls="--", label="optimo D = %.4f" % D)
        ax.axvline(semilla[0], color="C0", ls=":", label="semilla NCC")
        ax.set_title("costo contra D (sigma fijo en el optimo)")
        ax.set_xlabel("D [mV*ns/keVee]"); ax.set_ylabel("chi2/ndf")
        ax.legend(fontsize=8); ax.grid(alpha=.3)
        print()

    fig.tight_layout()
    salida = rutas.figura("06_modelo.png")
    fig.savefig(salida, dpi=130)
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
