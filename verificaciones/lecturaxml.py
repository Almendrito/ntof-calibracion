# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python verificaciones/lecturaxml.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.join(_os.path.dirname(_os.path.dirname(
    _os.path.abspath(__file__))), "codigo")
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------

import rutas

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def generar_espectro_y_tiempo(ruta_archivo, umbral_ruido, num_canales=1024):
    """
    Lee el archivo XML, calcula la duración total de la adquisición y genera el espectro.
    Retorna: canales, cuentas, duracion_segundos
    """
    tree = ET.parse(ruta_archivo)
    root = tree.getroot()
    
    # Extraer todos los eventos
    eventos = root.findall('.//Event')
    if not eventos:
        return np.array([]), np.array([]), 0
    
    # 1. Calcular el tiempo total de adquisición
    formato_tiempo = "%Y/%m/%d %H:%M:%S.%f"
    tiempo_inicial_str = eventos[0].find('Time').text
    tiempo_final_str = eventos[-1].find('Time').text
    
    t_inicio = datetime.strptime(tiempo_inicial_str, formato_tiempo)
    t_fin = datetime.strptime(tiempo_final_str, formato_tiempo)
    duracion_segundos = (t_fin - t_inicio).total_seconds()
    
    # 2. Procesamiento de pulsos (con filtros de ruido y tramas repetidas)
    areas_pulsos = []
    pulso_anterior = None
    
    for evento in eventos:
        canal = evento.find('.//CHN2')
        if canal is None:
            continue
            
        tiempos_pulso, voltajes = [], []
        for data in canal.findall('Data'):
            t_str, v_str = data.text.split(',')
            tiempos_pulso.append(float(t_str))
            voltajes.append(float(v_str))
            
        tiempos_pulso = np.array(tiempos_pulso)
        voltajes = np.array(voltajes)
        
        if len(voltajes) == 0:
            continue

        if pulso_anterior is not None and np.array_equal(voltajes, pulso_anterior):
            continue
        pulso_anterior = np.copy(voltajes)
        
        if np.max(voltajes) - np.min(voltajes) > umbral_ruido:
            # Usando np.trapezoid para versiones nuevas de NumPy
            area = abs(np.trapezoid(voltajes, x=tiempos_pulso))
            areas_pulsos.append(area)

    # 3. Generar el histograma
    cuentas, bordes_bins = np.histogram(areas_pulsos, bins=num_canales, range=(0, max(areas_pulsos) if areas_pulsos else 1))
    canales = (bordes_bins[:-1] + bordes_bins[1:]) / 2
    
    return canales, cuentas, duracion_segundos

def obtener_espectro_neto(cuentas_muestra, tiempo_muestra, cuentas_fondo, tiempo_fondo):
    """
    Normaliza los espectros al tiempo de adquisición más corto y resta el fondo.
    """
    tiempo_norm = min(tiempo_muestra, tiempo_fondo)
    print(f"Normalizando a tiempo base: {tiempo_norm:.2f} segundos")
    
    # Escalar cuentas proporcionalmente al tiempo normalizado
    cuentas_muestra_norm = cuentas_muestra * (tiempo_norm / tiempo_muestra)
    cuentas_fondo_norm = cuentas_fondo * (tiempo_norm / tiempo_fondo)
    
    # Resta canal por canal
    espectro_neto = cuentas_muestra_norm - cuentas_fondo_norm
    
    # Evitar cuentas negativas (físicamente imposibles, ocurren por fluctuaciones estadísticas)
    espectro_neto = np.maximum(espectro_neto, 0)

    print("--- Diagnóstico de Datos ---")
    print(f"Duración Muestra:    {tiempo_muestra:.2f} s | Total cuentas: {np.sum(cuentas_muestra)}")
    print(f"Duración Fondo:      {tiempo_fondo:.2f} s | Total cuentas: {np.sum(cuentas_fondo)}")
    print(f"Tasa Muestra:        {np.sum(cuentas_muestra)/tiempo_muestra:.2f} cps")
    print(f"Tasa Fondo:          {np.sum(cuentas_fondo)/tiempo_fondo:.2f} cps")
    print(f"Total cuentas netas: {np.sum(espectro_neto):.0f}")
    print("----------------------------")

    return espectro_neto


canales_muestra,cuentas_muestra, tiempo_muestra = generar_espectro_y_tiempo(rutas.MUESTRA, umbral_ruido=0.05)
canales_fondo,cuentas_fondo, tiempo_fondo = generar_espectro_y_tiempo(rutas.FONDO, umbral_ruido=0.05)

espectro_neto = obtener_espectro_neto(cuentas_muestra, tiempo_muestra, cuentas_fondo, tiempo_fondo)

plt.figure(figsize=(10, 6))

# Factor de escala temporal
factor_muestra = min(tiempo_muestra, tiempo_fondo) / tiempo_muestra
factor_fondo = min(tiempo_muestra, tiempo_fondo) / tiempo_fondo

# Graficar las tres señales
plt.step(canales_muestra, cuentas_muestra * factor_muestra, where="mid", color='blue', alpha=0.6, label='Muestra (normalizada)')
plt.step(canales_fondo, cuentas_fondo * factor_fondo, where="mid", color='gray', alpha=0.6, label='Fondo (normalizado)')
plt.step(canales_muestra, espectro_neto, where="mid", color='red', linewidth=1.8, label='Espectro Neto')

plt.title('Validación: Muestra vs Fondo vs Espectro Neto', fontsize=13)
plt.xlabel('Canal (Área del pulso [mV*ns])', fontsize=11)
plt.ylabel('Cuentas normalizadas', fontsize=11)
#plt.yscale('log')
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.show()