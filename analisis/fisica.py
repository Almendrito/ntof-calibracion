"""
fisica.py — lo poco de teoria que hace falta: el borde Compton y las fuentes.

El borde Compton es la maxima energia que un foton de `E_gamma` puede darle a un
electron, que es la retrodispersion a 180 grados:

    E_max = 2*alpha*E_gamma / (1 + 2*alpha)        con  alpha = E_gamma / m_e c^2

Ojo con la Ec. 1 impresa del paper: lleva `alpha` y no `2*alpha` en el numerador,
y asi da la energia del FOTON retrodispersado, no la del electron. Verificado
contra los propios datos del paper (Fig. 3): para Co-60 tiene que dar 963 y 1118 keV.
"""

M_E = 510.999          # masa del electron [keV]

# Solo las fuentes que efectivamente se midieron. Energias en keV.
FUENTES = {
    "Co-60": (1173.228, 1332.492),
    "Cs-137": (661.657,),
}


def borde_compton(e_gamma):
    """Energia del borde Compton [keV] de una linea gamma."""
    a = e_gamma / M_E
    return 2.0 * a * e_gamma / (1.0 + 2.0 * a)


def bordes(fuente):
    """Bordes Compton de una fuente del catalogo, en keV."""
    return tuple(borde_compton(e) for e in FUENTES[fuente])


def e_ref(fuente):
    """
    Energia de referencia del mapa canal-energia: el PRIMER borde.

    El mapa pasa por el origen (area cero <-> energia cero), asi que basta un
    ancla: D = B / E_ref. Con dos lineas se toma la de menor energia y la otra
    queda determinada por la misma D.
    """
    return bordes(fuente)[0]


if __name__ == "__main__":
    import sys
    sys.path.insert(0, __file__.rsplit("analisis", 1)[0])
    print("fuente     linea gamma [keV]   borde Compton [keV]")
    for f, lineas in FUENTES.items():
        for eg in lineas:
            print("%-9s  %10.3f        %10.1f" % (f, eg, borde_compton(eg)))
