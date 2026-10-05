"""
physics.py — Etapa 1 del metodo de calibracion por Compton edge.

Codifica la fisica de la Seccion II del paper:
  - compton_edge_energy : E_max, la posicion del Compton edge   [Ec. 1, corregida]
  - dsigma_dE           : forma teorica del espectro de electrones
                          de retroceso (Klein-Nishina en energia) [Ec. 10]

Esta es la "plantilla" (S2) contra la que despues se correlaciona (NCC)
y se ajusta (SA) el espectro medido. No necesita datos: es pura teoria.
Todo el modulo trabaja en keV.
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/etapas/physics.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas


import numpy as np

MEC2 = 510.999  # energia en reposo del electron, m_e c^2 [keV]


def compton_edge_energy(E_gamma):
    """
    Energia maxima del electron de retroceso = posicion del Compton edge.
    Ocurre en retrodispersion del foton (theta = 180 grados).

        E_max = 2*alpha*E_gamma / (1 + 2*alpha),   alpha = E_gamma / m_e c^2

    OJO (erratum del paper): la Ec.1 impresa, E_gamma/(1+2 alpha), es la
    energia del FOTON retrodispersado, no la del electron. Aqui va la forma
    correcta; reproduce los valores de los datos del propio paper (Fig. 3).
    """
    alpha = E_gamma / MEC2
    return 2.0 * alpha * E_gamma / (1.0 + 2.0 * alpha)


def dsigma_dE(E_e, E_gamma, A=1.0):
    """
    dSigma/dE_e : Klein-Nishina angular (Ec.2) pasada a energia del electron
    via el Jacobiano que derivamos en clase -> Ec. 10 del paper, termino a termino.

    E_e     : array de energias del electron [keV]
    E_gamma : energia del gamma incidente [keV]
    A       : factor de amplitud (el 'A' del paper). El constante pi*r0^2/m_e c^2
              queda absorbido en A: para correlacionar/ajustar solo importa la FORMA.

    Vale 0 fuera de (0, E_max): arriba de E_max no hay cinematica -> el escalon.
    """
    E_e = np.asarray(E_e, dtype=float)
    alpha = E_gamma / MEC2
    E_max = compton_edge_energy(E_gamma)

    valido = (E_e > 0.0) & (E_e < E_max)
    Ee = np.where(valido, E_e, np.nan)     # evita 0/0 fuera del rango valido
    resto = E_gamma - Ee                    # (h*nu - E_e)

    bracket = 2.0 + (Ee / resto) ** 2 * (
        1.0 / alpha ** 2
        - 2.0 * resto / (alpha * Ee)
        + resto / E_gamma
    )

    out = (A / alpha ** 2) * bracket
    return np.where(valido, out, 0.0)


# ----------------------------------------------------------------------
# 3) que fuente se esta midiendo
# ----------------------------------------------------------------------
# Catalogo de fuentes gamma de calibracion, ordenado por energia del borde.
# Cada entrada son las lineas gamma DOMINANTES, en keV. A lo sumo dos: el modelo
# de bordes solapados (Ec. 11-13) solo sabe combinar dos componentes, asi que
# meter una tercera linea la dejaria silenciosamente afuera. Para fuentes con
# muchas lineas (Ba-133, Eu-152) va la dominante y las demas quedan como parte
# del fondo del ajuste; en esos casos el numero sale peor y hay que mirarlo.
FUENTES = {
    "Am-241": (59.5409,),
    "Co-57": (122.0613,),
    "Ba-133": (356.0129,),
    "Cs-137": (661.657,),
    "Mn-54": (834.848,),
    "Na-22": (511.0, 1274.537),
    "Co-60": (1173.2, 1332.5),
    "K-40": (1460.822,),
    "Tl-208": (2614.511,),
}

# Lo que hay que saber antes de usar cada una para calibrar por Compton edge.
NOTAS_FUENTE = {
    "Am-241": "borde en 11 keV: por debajo del umbral de casi cualquier SiPM. "
              "Sirve para ver el fotopico en NaI, no para esto.",
    "Co-57": "borde en 39 keV: muy bajo, cae dentro del exceso de bajo canal.",
    "Ba-133": "356 keV es la dominante (62%), pero 303 y 384 keV aportan y no "
              "estan en el modelo: el borde sale ensanchado.",
    "Cs-137": "el caso limpio: UNA linea, un borde en 477.3 keV, D = B/E_max "
              "directo.",
    "Mn-54": "una sola linea, borde en 639.2 keV. Media vida corta (312 d).",
    "Na-22": "511 (aniquilacion) y 1274.5 keV: los bordes quedan MUY separados "
             "(341 y 1062 keV), asi que da dos puntos casi independientes. La "
             "de 511 es de aniquilacion del positron, no del nucleo.",
    "Co-60": "dos lineas casi iguales en intensidad; los bordes (963.4 y "
             "1118.1 keV) se solapan. Es la fuente de este proyecto.",
    "K-40": "1460.8 keV, borde en 1243.4. Actividad bajisima: sirve como linea "
            "de FONDO natural, no como fuente de calibracion.",
    "Tl-208": "2614.5 keV (cadena del Th-232), borde en 2381.8. El punto mas "
              "alto disponible sin acelerador.",
}


class Fuente:
    """
    Las lineas gamma de lo que se esta midiendo, y sus bordes Compton.

    Por que existe: hasta la Etapa 9 el proyecto era Co-60 y sus dos lineas
    estaban cableadas en el codigo (EG1, EG2 en gui.py). Con la sesion del
    26-ago-2026 entro Cs-137, que tiene UNA sola linea: el modelo de dos bordes
    (Ec. 11-13) directamente no aplica y hay que caer al de un borde. Esta clase
    es lo unico que sabe cuantas lineas hay; el resto del codigo le pregunta.

    `E_ref` es el ancla del mapa por el origen (D = B / E_ref): el borde de la
    PRIMERA linea, que es contra el que se mide B en todo el proyecto (963.4 keV
    para Co-60, 477.3 para Cs-137).

    Las lineas se guardan ORDENADAS, asi que E_ref siempre es el borde mas bajo.
    """

    def __init__(self, nombre, lineas=None):
        if lineas is None:
            if nombre not in FUENTES:
                raise KeyError(f"fuente desconocida: {nombre}. Conocidas: "
                               + ", ".join(FUENTES))
            lineas = FUENTES[nombre]
        lineas = tuple(sorted(float(E) for E in lineas
                              if E is not None and float(E) > 0))
        if not lineas:
            raise ValueError("una fuente necesita al menos una linea gamma")
        self.nombre = str(nombre)
        self.lineas = lineas
        self.bordes = tuple(compton_edge_energy(E) for E in lineas)

    @property
    def E_ref(self):
        """Borde de la primera linea [keV]: el ancla del mapa por el origen."""
        return self.bordes[0]

    @property
    def E_ref2(self):
        """Borde de la segunda linea, o None si la fuente tiene una sola."""
        return self.bordes[1] if len(self.bordes) > 1 else None

    @property
    def n_lineas(self):
        return len(self.lineas)

    @property
    def doble(self):
        """True si el modelo de dos bordes solapados (Ec. 11-13) tiene sentido."""
        return len(self.lineas) >= 2

    def nota(self):
        """Que hay que saber de esta fuente para calibrar por Compton edge."""
        return NOTAS_FUENTE.get(self.nombre, "")

    def avisos(self):
        """
        Lista de advertencias DERIVADAS de las lineas, no de una tabla.

        Se calculan aca y no en la GUI para que valgan tambien cuando la fuente
        se tipea a mano ("manual"): el usuario puede poner cualquier energia y
        conviene que el programa diga en voz alta cuando esa eleccion rompe
        alguna hipotesis del metodo.
        """
        avisos = []
        if self.E_ref < 100.0:
            avisos.append(
                f"el borde queda en {self.E_ref:.0f} keV, dentro de la zona "
                "donde vive el exceso de bajo canal y cerca del umbral de "
                "trigger: el ajuste va a estar dominado por lo que NO es "
                "Compton (CLAUDE.md 7.32, 7.35).")
        if len(self.lineas) > 2:
            extra = ", ".join(f"{E:.1f}" for E in self.lineas[2:])
            avisos.append(
                f"hay {len(self.lineas)} lineas y el modelo solo combina dos: "
                f"{extra} keV quedan sin modelar.")
        if self.doble:
            sep = (self.bordes[1] - self.bordes[0]) / self.bordes[0]
            if sep < 0.25:
                avisos.append(
                    f"los dos bordes estan a {100*sep:.0f}% uno de otro: se "
                    "solapan y la descomposicion interna aporta poco. Conviene "
                    "ligar el origen (D = B1/E_max,1), CLAUDE.md 7.16.")
        return avisos

    def es_preset(self):
        """True si las lineas son exactamente las del preset de ese nombre."""
        base = FUENTES.get(self.nombre)
        return base is not None and self.lineas == tuple(sorted(base))

    def descripcion(self):
        return (f"{self.nombre}: "
                + " y ".join(f"{E:.1f}" for E in self.lineas)
                + " keV  ->  borde"
                + ("s " if self.doble else " ")
                + " y ".join(f"{E:.1f}" for E in self.bordes) + " keV")

    def __repr__(self):
        return f"Fuente({self.nombre!r}, {self.lineas})"


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for nombre in FUENTES:
        fu = Fuente(nombre)
        print("  " + fu.descripcion())
        if fu.nota():
            print("      " + fu.nota())
        for a in fu.avisos():
            print("      ! " + a)
    print()
    lineas = {"1173 keV": 1173.2, "1333 keV": 1332.5}   # lineas gamma del Co-60
    E = np.linspace(0, 1400, 4000)                       # eje E_e [keV]

    plt.figure(figsize=(9, 5.5))
    for etiqueta, Eg in lineas.items():
        y = dsigma_dE(E, Eg)
        Emax = compton_edge_energy(Eg)
        ln, = plt.plot(E, y, lw=1.8, label=f"{etiqueta}  ->  edge {Emax:.0f} keV")
        plt.axvline(Emax, color=ln.get_color(), ls="--", alpha=0.6)
        print(f"E_gamma = {Eg:7.1f} keV  ->  Compton edge E_max = {Emax:7.1f} keV")

    plt.title("Etapa 1 - Lineshape teorica de electrones de retroceso (Co-60)")
    plt.xlabel("Energia del electron de retroceso  E_e  [keV]")
    plt.ylabel("dSigma/dE_e   (forma, u.a.)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(rutas.figura("lineshape_co60.png"), dpi=130)
    print("Figura guardada: figuras/lineshape_co60.png")
