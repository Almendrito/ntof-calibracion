"""
pulsos — caracterizacion de PULSOS INDIVIDUALES: ancho a media altura,
intensidad y forma.

No es una etapa del pipeline de calibracion (eso vive en `analisis/`): es una
herramienta exploratoria para mirar los pulsos uno a uno. Se construye ENCIMA
de `analisis/drs_xml.py` (lectura del XML) y `analisis/pulso.py` (linea base,
onset, pico, integral acumulada), que se importan y nunca se reimplementan.

    rutas   donde esta cada cosa (mas el catalogo de corridas de `analisis`)
    ancho   FWHM temporal de un pulso individual
    lote    cargar N pulsos + area y FWHM del lote + sugerir ventana
    gui     la interfaz (`python pulsos/gui.py`)
"""
