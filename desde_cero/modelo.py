"""
Etapa 5 — el modelo directo: la forma teorica en el eje medido, con resolucion
(paper, Ec. 8, 9 y 11-13).

    python desde_cero/modelo.py

Que resuelve esta etapa
-----------------------
La Etapa 4 da una forma en keV con un escalon perfecto en el borde. El dato
(Etapa 3) es un histograma en mV*ns con un borde borroneado. El modelo directo
hace el camino de la teoria al dato:

    1. mapa energia -> canal:  ch = D * E          (Ec. 8, por el origen)
    2. resolucion:             convolucion con una gaussiana de ancho sigma
                               (Ec. 9)
    3. cuentas:                A = numero total de cuentas de la fuente

Tres parametros: A, D [mV*ns/keV] y sigma_rel. D es lo que se quiere medir.

Las decisiones, y por que
-------------------------
1. **El mapa pasa por el origen.** Area cero <-> energia cero, sin ordenada.
   Es lo que dice la Ec. 8 (E_eff = E * B/E_max) y lo que tiene que pasar si
   el area es proporcional a la luz. Con eso, fijar D fija el modelo entero.

2. **A son CUENTAS.** La forma se normaliza a integral 1 en energia y cada
   canal recibe A * forma(E) * (ancho del canal en keV = paso/D). Sumando
   todos los canales da A. En el paper A absorbe el prefactor de la Ec. 10
   y no significa nada; aca es el numero de eventos de la fuente, y la
   conservacion de cuentas es una prueba (ver `verificar_modelo.py`).

3. **El kernel suma 1** (no integra 1). Se aplica con `np.convolve` sobre un
   histograma, y un kernel de suma 1 reparte las cuentas de cada canal entre
   sus vecinos sin crear ni destruir ninguna. Por eso A y sigma quedan
   DESACOPLADOS: ensanchar baja el pico pero no cambia el total. Con otra
   normalizacion, subir sigma cambiaria la altura y el ajuste podria
   compensarlo con A: dos parametros haciendo lo mismo.

4. **La convolucion se hace en una grilla mas ancha que la ventana** y
   despues se recorta. `np.convolve` rellena con CEROS fuera de lo que se le
   da: si se le da solo la ventana de ajuste, el borde izquierdo de la
   ventana no recibe las cuentas que la resolucion trae desde abajo. El
   modelo sale bajo ahi, y el ajuste lo compensa ACHICANDO sigma (-18 % en
   sigma y +3.9 % en D en la prueba de §7.14). La grilla va desde energia
   cero menos 4 sigma hasta el ultimo borde (o la ventana) mas 4 sigma.

5. **Co-60: dos lineas, UNA sola D.** Los dos bordes estan en D*E_max1 y
   D*E_max2, asi que su separacion en canal es D*(E_max2 - E_max1) = D*154.7
   keV. La Ec. 12 impresa usa la diferencia de las LINEAS (159.3 keV), que
   esta mal: D convierte energia del ELECTRON, y lo que se separa en canal
   son los bordes. Con una sola D eso no hay que escribirlo: sale solo.
   Las dos lineas pesan su intensidad (ENSDF) y cada una ya trae su seccion
   eficaz de Klein-Nishina. Eso es la aproximacion de detector DELGADO (cada
   gamma interactua a lo sumo una vez); el scattering multiple es la
   Etapa 8.

6. **Un solo sigma, en canales: sigma_ch = sigma_rel * D * E_ref**, con E_ref
   el borde mas bajo. El paper usa un sigma por borde (Ec. 11). La
   resolucion de un centellador crece como sqrt(E), asi que el del segundo
   borde del Co-60 tendria que ser sqrt(1118/963) = 1.08 veces el del
   primero; con un solo sigma eso se ignora. Se elige asi por la leccion de
   §7.23: con sigma1 y sigma2 libres el ajuste los usa para cualquier cosa
   (salian invertidos). Si el 8 % importa, lo va a decir el chi2.
"""

import os
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from desde_cero import fisica                               # noqa: E402

N_SIGMAS = 4.0     # radio del kernel y del relleno de la grilla


def lineas_de(fuente):
    """(energias, pesos) de una fuente de `fisica.FUENTES`."""
    l = fisica.FUENTES[fuente]
    return np.array([e for e, _ in l]), np.array([w for _, w in l])


def forma_energia(E, energias, pesos):
    """
    Densidad de probabilidad en energia del electron [1/keV], integral 1.

    Suma de las Klein-Nishina de cada linea pesadas por su intensidad. Las
    secciones eficaces NO se normalizan por separado: una linea con mas
    sigma interactua mas, y eso ya esta en la Ec. 10 en cm^2/keV.
    """
    E = np.asarray(E, float)
    tot = sum(w * fisica.sigma_total_cerrada(e) for e, w in zip(energias, pesos))
    return sum(w * fisica.dsigma_de(E, e) for e, w in zip(energias, pesos)) / tot


_ACUMULADAS = {}      # energia de la linea -> (grilla T, F(T)); se arma una vez


def _acumulada(e_gamma, n=20001):
    """
    F(T) = integral de 0 a T de dsigma/dT [cm^2], en una grilla fina de 0 a
    E_max. F(E_max) es la seccion eficaz total de la linea. Es un calculo de
    la TEORIA (no depende de ningun dato), asi que se guarda en memoria.
    """
    if e_gamma not in _ACUMULADAS:
        T = np.linspace(0.0, float(fisica.e_max(e_gamma)), n)
        f = fisica.dsigma_de(T, e_gamma)
        f[-1] = fisica.dsigma_de(np.array([T[-1] * (1 - 1e-12)]), e_gamma)[0]
        F = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(T))))
        _ACUMULADAS[e_gamma] = (T, F)
    return _ACUMULADAS[e_gamma]


def cuentas_por_bin(E_inf, E_sup, energias, pesos):
    """
    Fraccion de las cuentas que cae entre E_inf y E_sup [keV], sumando las
    lineas por intensidad x seccion eficaz: F(E_sup) - F(E_inf), normalizado
    a que el total sea 1.

    Por que integrar el bin y no evaluar la forma en su centro: la forma de
    Klein-Nishina tiene un ESCALON en el borde. Evaluada en el centro de cada
    bin, cuando D cambia y el borde cruza el centro de un bin, ese bin pasa de
    0 al valor del pico DE GOLPE: el modelo, y con el el chi2, es discontinuo
    en D. Medido en el Cs-137 real: el chi2 perfilado en D saltaba 22 unidades
    a una fraccion de barra de error del minimo, y el hessiano -la barra de
    error- no significaba nada. Integrado, el bin del borde se llena de a poco
    a medida que el borde lo cruza: continuo. (El paper, y replica_paper,
    evaluan en el centro.)
    """
    tot, acum = 0.0, 0.0
    for e, w in zip(energias, pesos):
        T, F = _acumulada(float(e))
        acum = acum + w * (np.interp(E_sup, T, F, left=0.0, right=F[-1])
                           - np.interp(E_inf, T, F, left=0.0, right=F[-1]))
        tot += w * F[-1]
    return acum / tot


def kernel(sigma_ch, paso):
    """Ec. 9 muestreada en la grilla y normalizada a SUMA 1."""
    r = max(1, int(np.ceil(N_SIGMAS * sigma_ch / paso)))
    x = np.arange(-r, r + 1) * paso
    k = np.exp(-0.5 * (x / sigma_ch) ** 2)
    return k / k.sum()


def grilla_extendida(canales, D, sigma_ch, e_tope):
    """
    Grilla uniforme ALINEADA con `canales` (mismo paso, mismos puntos) que va
    de energia cero - 4 sigma hasta max(ventana, borde) + 4 sigma.

    Devuelve (ext, i0): `ext[i0:i0 + len(canales)]` son exactamente los
    canales pedidos, asi que recortar es indexar y no interpolar.
    """
    ch = np.asarray(canales, float)
    paso = float(ch[1] - ch[0])
    if not np.allclose(np.diff(ch), paso, rtol=1e-9, atol=1e-12):
        raise ValueError("la grilla de canales tiene que ser uniforme")
    margen = N_SIGMAS * sigma_ch
    abajo = int(np.ceil((ch[0] - (0.0 - margen)) / paso))
    arriba = int(np.ceil((max(ch[-1], D * e_tope) + margen - ch[-1]) / paso))
    abajo, arriba = max(abajo, 0), max(arriba, 0)
    ext = ch[0] + paso * np.arange(-abajo, ch.size + arriba)
    return ext, abajo


def modelo(canales, A, D, sigma_rel, fuente=None, energias=None, pesos=None,
           extender=True, por_bin=True, plantilla=None):
    """
    Cuentas esperadas en cada canal (centros de bin, grilla uniforme).

    A         : cuentas totales de la fuente
    D         : ganancia [mV*ns/keV]; el borde de la linea i esta en D*E_max_i
    sigma_rel : resolucion, sigma_ch / (D * E_ref), E_ref = el borde mas bajo
    fuente    : "Cs-137" / "Co-60", o pasar `energias` y `pesos` a mano
    extender  : False reproduce el error de §7.14 (convolucion solo en la
                ventana). Esta para poder mostrarlo, no para usarlo.
    por_bin   : True (default) integra la forma sobre cada bin
                (`cuentas_por_bin`); False la evalua en el centro del bin, como
                el paper. False solo sirve para compararse con replica_paper:
                el modelo sale discontinuo en D.
    plantilla : None = Klein-Nishina (el paper). Si no, cualquier objeto con
                `fraccion(E_inf, E_sup)` (cuentas entre dos energias, total 1),
                `densidad(E)` [1/keV] y `e_tope` [keV] (Etapa 8: la
                simulacion de Geant4). Lo demas no cambia: el mapa, la
                resolucion y sigma_rel se definen igual, con E_ref = el borde
                Compton mas bajo de la `fuente`, asi que D y sigma_rel de las
                dos plantillas se comparan directo.

    La suma sobre canales >= 0 es MENOR que A: la resolucion manda a canales
    negativos parte de las cuentas de energia ~0 (6 % con sigma_rel = 0.15 en
    el Cs-137). Es lo correcto -esas cuentas existen en el modelo, solo que
    fuera de lo que se mira- y no afecta al ajuste, cuya ventana arranca
    muy por encima de cero.
    """
    if fuente is not None:
        energias, pesos = lineas_de(fuente)
    energias = np.atleast_1d(np.asarray(energias, float))
    pesos = np.ones_like(energias) if pesos is None else np.asarray(pesos, float)
    ch = np.asarray(canales, float)
    paso = float(ch[1] - ch[0])
    bordes = fisica.e_max(energias)
    sigma_ch = sigma_rel * D * float(bordes.min())
    if not (A >= 0 and D > 0 and sigma_ch > 0):
        return np.zeros_like(ch)

    # la simulacion llega mas arriba que el borde (scattering multiple y
    # absorcion total), y la grilla tiene que cubrir hasta donde haya cuentas
    e_tope = float(bordes.max()) if plantilla is None else         max(float(bordes.max()), float(plantilla.e_tope))
    if extender:
        ext, i0 = grilla_extendida(ch, D, sigma_ch, e_tope)
    else:
        ext, i0 = ch, 0
    # cuentas por canal SIN resolucion: densidad en E por el ancho del canal
    # en keV. E = ch/D es la Ec. 8 invertida.
    E_inf, E_sup = (ext - 0.5 * paso) / D, (ext + 0.5 * paso) / D
    if plantilla is not None:
        crudo = A * (plantilla.fraccion(E_inf, E_sup) if por_bin
                     else plantilla.densidad(ext / D) * (paso / D))
    elif por_bin:
        crudo = A * cuentas_por_bin(E_inf, E_sup, energias, pesos)
    else:
        crudo = A * forma_energia(ext / D, energias, pesos) * (paso / D)
    k = kernel(sigma_ch, paso)
    if k.size > ext.size:
        # np.convolve(mode="same") devuelve el largo del MAS LARGO de los dos:
        # con un kernel mas ancho que la senal el recorte de abajo saldria
        # corrido (§7.2). Solo puede pasar con extender=False y sigma enorme.
        r = (ext.size - 1) // 2
        k = k[k.size // 2 - r:k.size // 2 + r + 1]
        k = k / k.sum()
    y = np.convolve(crudo, k, mode="same")
    return y[i0:i0 + ch.size]


def mitad_de_altura(canales, y, desde):
    """
    El metodo clasico: el canal donde el espectro cae al 50 % de su maximo,
    buscando por encima de `desde`. Interpolado linealmente.

    El maximo es la acumulacion de Klein-Nishina justo antes del borde,
    borroneada. Esta aca para MOSTRAR su sesgo, no para usarlo.
    """
    ch, y = np.asarray(canales, float), np.asarray(y, float)
    m = ch >= desde
    c, v = ch[m], y[m]
    i = int(np.argmax(v))
    mitad = 0.5 * v[i]
    j = i + int(np.flatnonzero(v[i:] < mitad)[0])
    return float(c[j - 1] + (mitad - v[j - 1]) / (v[j] - v[j - 1]) * (c[j] - c[j - 1]))


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    D, paso = 1.9, 2.5
    ch = np.arange(0.5 * paso, 3000.0, paso)
    sigmas = (0.03, 0.08, 0.15)
    fig, ejes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, fuente in zip(ejes[:2], ("Cs-137", "Co-60")):
        E, _ = lineas_de(fuente)
        bordes = D * fisica.e_max(E)
        for s in sigmas:
            y = modelo(ch, 1e5, D, s, fuente)
            ax.plot(ch, y, lw=1.2, label="sigma_rel = %.2f  (suma %.1f)" % (s, y.sum()))
            if fuente == "Cs-137":
                h = mitad_de_altura(ch, y, 0.5 * bordes[0])
                ax.plot(h, np.interp(h, ch, y), "k|", ms=14)
        for b in bordes:
            ax.axvline(b, color="0.4", ls=":", lw=1)
        ax.set(xlabel="canal [mV*ns]", ylabel="cuentas / canal",
               title="%s, D = %.2f mV*ns/keV, A = 1e5 cuentas\n(punteado: borde "
                     "verdadero D*E_max; | : 50 %% de altura)" % (fuente, D),
               xlim=(0, 1.35 * bordes.max()))
        ax.legend(fontsize=8)

    ax = ejes[2]
    B = D * float(fisica.e_max(661.657))
    sr = np.linspace(0.01, 0.20, 20)
    sesgo = [mitad_de_altura(ch, modelo(ch, 1e5, D, s, "Cs-137"), 0.5 * B) / B - 1
             for s in sr]
    ax.plot(sr, 100 * np.array(sesgo), "o-", ms=4)
    ax.set(xlabel="sigma_rel = sigma / borde", ylabel="(50 % de altura) / borde - 1 [%]",
           title="sesgo del metodo del 50 % (Cs-137, sin ruido)")
    for ax in ejes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    salida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras",
                          "etapa5_modelo.png")
    fig.savefig(salida, dpi=120)
    print("  sesgo del 50 %% de altura: %s"
          % ", ".join("%.2f -> %+.1f %%" % (s, 100 * b)
                      for s, b in zip(sr[::3], sesgo[::3])))
    print("  figura: %s" % salida)


if __name__ == "__main__":
    _demo()
