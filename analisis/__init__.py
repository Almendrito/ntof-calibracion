"""
analisis — procesamiento de datos del DRS4, reescrito de cero (27-ago-2026).

Cadena, en el orden en que se ejecuta:

    rutas      donde esta cada cosa
    drs_xml    E0 . abrir el XML y reconstruir los pulsos (t, v)
    pulso      E1 . linea base, onset, integral acumulada -> ventana de integracion
    area       E2 . integrar todos los eventos (dos metodos) + cache
    espectro   E3 . histograma con binning comun, tiempo vivo, resta de fondo

Cada modulo corre solo (`python analisis/<modulo>.py`) y deja una figura en
`analisis/figuras/` que demuestra lo que hace.

Este paquete NO importa nada de `codigo/`: es una implementacion independiente.
"""
