import numpy as np
from scipy.optimize import curve_fit

def extraer_crudos(ruta):
    """
    Extrae la matriz de datos desde un archivo .npz.
    Retorna: (datos_crudos, is_raw_bytes)
    """
    try:
        datos = np.load(ruta, allow_pickle=True)
        if "datos_bytes" in datos:
            return datos["datos_bytes"], True
        elif "datos_voltios" in datos:
            return datos["datos_voltios"], False
        else:
            raise ValueError("Formato NPZ no reconocido. Faltan claves 'datos_bytes' o 'datos_voltios'.")
    except Exception as e:
        raise RuntimeError(f"Error al cargar el archivo NPZ: {str(e)}")


def filtrar_duplicados_ventana(matriz, ventana=100):
    """
    Filtra pulsos duplicados en una ventana de tiempo (Post-análisis).
    """
    n = matriz.shape[0]
    if n <= 1:
        return matriz, 0
    
    mascara = np.ones(n, dtype=bool)
    for i in range(1, ventana + 1):
        if n > i:
            distintos = np.any(matriz[i:] != matriz[:-i], axis=1)
            mascara[i:] &= distintos
            
    filtrada = matriz[mascara]
    return filtrada, n - filtrada.shape[0]


# ── Herramientas matemáticas para ajuste de decaimiento ──
def modelo_exponencial(t, A, tau, C):
    """Modelo matemático para la cola del pulso."""
    return A * np.exp(-t / tau) + C


def ajustar_cola_pulso(t_ns, v_corr, idx_pico):
    """
    Ajusta una exponencial a la cola del pulso con tolerancia al ruido 
    de alta frecuencia y forzado de ventana para after-pulses tempranos.
    """
    max_muestras_busqueda = 60 
    idx_seguro_max = min(len(t_ns) - 1, idx_pico + max_muestras_busqueda)
    
    v_abs = v_corr 

    # Empezamos el ajuste 1 muestra DESPUÉS del pico para saltarnos 
    # el redondeo que produce el ancho de banda analógico del osciloscopio
    idx_inicio_ajuste = idx_pico + 1 
    idx_fin_ajuste = idx_inicio_ajuste + 1 
    
    consecutivos_subiendo = 0
    
    while idx_fin_ajuste < idx_seguro_max:
        # 1. Si cruza cero (sobreimpulso muy fuerte): detener.
        if v_abs[idx_fin_ajuste] <= 0:
            break
            
        # 2. Tolerancia al ruido: Paramos SOLO si sube 2 veces seguidas
        if v_abs[idx_fin_ajuste] > v_abs[idx_fin_ajuste - 1]:
            consecutivos_subiendo += 1
            if consecutivos_subiendo >= 2: # Detectado after-pulse real
                idx_fin_ajuste -= 1 # Retrocedemos antes de la subida
                break
        else:
            consecutivos_subiendo = 0 # Reset si vuelve a bajar
            
        idx_fin_ajuste += 1
        
    t_fit = t_ns[idx_inicio_ajuste:idx_fin_ajuste]
    v_fit = v_abs[idx_inicio_ajuste:idx_fin_ajuste]
    
    # 3. MECANISMO DE RESCATE: 
    # Si la ventana es menor a 4 puntos (porque el after-pulse estaba pegadísimo),
    # forzamos una ventana de 8 muestras para obligarlo a sacar una curva.
    if len(t_fit) < 4: 
        idx_fin_forzado = min(len(t_ns) - 1, idx_inicio_ajuste + 8)
        t_fit = t_ns[idx_inicio_ajuste:idx_fin_forzado]
        v_fit = v_abs[idx_inicio_ajuste:idx_fin_forzado]
        if len(t_fit) < 4: # Si definitivamente no hay datos (final del vector)
            return None, None
            
    t_shifted = t_fit - t_fit[0]
    
    # Parámetros iniciales: A = voltaje inicial, Tau = 2.5ns, Offset = 0
    p0 = [v_fit[0], 2.5, 0.0] 
    
    # Límites más holgados para permitir que curve_fit pelee con los after-pulses
    bounds = (
        (0.001, 0.2, -0.05),  
        (0.500, 20.0, 0.05)   
    )
    
    try:
        popt, _ = curve_fit(modelo_exponencial, t_shifted, v_fit, p0=p0, bounds=bounds, maxfev=3000)
        
        # Generamos la curva teórica proyectándola desde el pico
        t_completo_shifted = t_ns[idx_pico:] - t_ns[idx_inicio_ajuste]
        curva_teorica = modelo_exponencial(t_completo_shifted, *popt)
        return popt, curva_teorica
        
    except (RuntimeError, ValueError):
        return None, None


# ── Cálculo Principal ──
def calcular_area(t_ns, v, baseline_muestras, polaridad, umbral_amplitud,
                  metodo_integracion, pre_pico, post_pico, umbral_inicio, umbral_fin):
    """
    Calcula el área del pulso (integración).
    Retorna: (area, offset, idx_inicio, idx_fin, v_corregido, v_teorico_full)
    """
    nb = baseline_muestras
    offset = np.mean(v[:nb]) if nb > 0 else v[0]
    v_corr = v - offset
    v_abs  = v_corr if polaridad == "Positivo" else -v_corr
    v_teorico_full = None 

    if np.max(v_abs) < umbral_amplitud:
        return 0.0, offset, None, None, v_corr, None

    idx_pico = int(np.argmax(v_abs))

    if metodo_integracion == "ventana_fija":
        idx_i = max(0, idx_pico - pre_pico)
        idx_f = min(len(v) - 1, idx_pico + post_pico)

    elif metodo_integracion == "umbral":
        hits = np.where(v_abs > umbral_inicio)[0]
        if len(hits) == 0:
            return 0.0, offset, None, None, v_corr, None
        idx_i    = hits[0]
        fin      = np.where(v_abs[idx_pico:] < umbral_fin)[0]
        idx_f    = idx_pico + fin[0] if len(fin) > 0 else len(v) - 1

    elif metodo_integracion == "ajuste_exponencial":
        idx_i = max(0, idx_pico - pre_pico)
        
        popt, curva_teorica = ajustar_cola_pulso(t_ns, v_abs, idx_pico)
        
        if popt is not None:
            indices_debajo = np.where(curva_teorica < umbral_fin)[0]
            if len(indices_debajo) > 0:
                idx_f = idx_pico + indices_debajo[0]
            else:
                idx_f = min(len(v) - 1, idx_pico + post_pico)
                
            v_teorico_full = np.zeros_like(v_corr)
            v_teorico_full[idx_pico:] = curva_teorica if polaridad == "Positivo" else -curva_teorica
        else:
            idx_f = min(len(v) - 1, idx_pico + post_pico)

    area = float(np.abs(np.trapezoid(v_corr[idx_i:idx_f + 1], t_ns[idx_i:idx_f + 1])))
    
    return area, offset, idx_i, idx_f, v_corr, v_teorico_full