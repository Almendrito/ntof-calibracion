import numpy as np
import pandas as pd
from core import processing

class DataManager:
    """
    Administra el estado de los datos en memoria (RAM) y orquesta el procesamiento.
    No interactúa con la interfaz gráfica (Tkinter) directamente.
    """
    def __init__(self):
        # Datos crudos
        self.pulsos_raw_muestra = None
        self.is_raw_muestra = False
        self.pulsos_raw_fondo = None
        self.is_raw_fondo = False
        
        # Datos calibrados
        self.pulsos_muestra = None
        self.pulsos_fondo = None
        
        # Vectores temporales
        self.t_vector = None
        self.t_ns = None
        
        # Resultados de integración
        self.areas_acumuladas = []
        self.areas_fondo = []

    def cargar_muestra(self, ruta):
        raw_data, is_raw = processing.extraer_crudos(ruta)
        self.pulsos_raw_muestra, elim = processing.filtrar_duplicados_ventana(raw_data, ventana=3)
        self.is_raw_muestra = is_raw
        
        num_p, num_m = self.pulsos_raw_muestra.shape
        self.t_vector = np.arange(num_m)
        
        return num_p, elim

    def cargar_fondo(self, ruta):
        raw_data, is_raw = processing.extraer_crudos(ruta)
        self.pulsos_raw_fondo, elim = processing.filtrar_duplicados_ventana(raw_data, ventana=3)
        self.is_raw_fondo = is_raw
        
        return self.pulsos_raw_fondo.shape[0], elim

    def aplicar_conversion(self, v_scale, v_offset, y_ref, f_muestreo_ghz, factor_escala_hw):
        y_inc = (v_scale * 8) / 256
        
        # Conversión de la muestra
        if self.pulsos_raw_muestra is not None:
            if self.is_raw_muestra:
                self.pulsos_muestra = (self.pulsos_raw_muestra.astype(float) - y_ref) * y_inc + v_offset
            else:
                self.pulsos_muestra = self.pulsos_raw_muestra.astype(float)
                
        # Conversión del background
        if self.pulsos_raw_fondo is not None:
            if self.is_raw_fondo:
                self.pulsos_fondo = (self.pulsos_raw_fondo.astype(float) - y_ref) * y_inc + v_offset
            else:
                self.pulsos_fondo = self.pulsos_raw_fondo.astype(float)
                
        # Conversión temporal
        if self.t_vector is not None:
            dt = (1.0 / f_muestreo_ghz) / factor_escala_hw
            self.t_ns = self.t_vector * dt

    def get_dt_ns(self, f_muestreo_ghz, factor_escala_hw):
        """Retorna el diferencial de tiempo en nanosegundos para visualización UI."""
        return (1.0 / f_muestreo_ghz) / factor_escala_hw

    def get_pulso(self, idx, fuente="Muestra"):
        """Retorna el vector de voltaje de un pulso específico."""
        matriz = self.pulsos_fondo if fuente == "Background" else self.pulsos_muestra
        if matriz is None:
            raise ValueError(f"No hay datos cargados para la fuente: {fuente}")
            
        idx_seguro = int(np.clip(idx, 0, matriz.shape[0] - 1))
        return matriz[idx_seguro], idx_seguro

    def _procesar_matriz(self, matriz, kwargs_procesamiento):
        areas = []
        for i in range(matriz.shape[0]):
            area, _, _, _, _, _ = processing.calcular_area(
                t_ns=self.t_ns,
                v=matriz[i],
                **kwargs_procesamiento
            )
            if area > 0:
                areas.append(area)
        return areas

    def integrar_todos(self, kwargs_procesamiento):
        """
        Orquesta la integración de todos los pulsos usando los parámetros dados.
        kwargs_procesamiento debe ser un diccionario con los parámetros requeridos por processing.calcular_area.
        """
        if self.pulsos_muestra is None or self.t_ns is None:
            raise ValueError("Faltan datos de muestra o aplicar calibración antes de integrar.")
            
        self.areas_acumuladas = self._procesar_matriz(self.pulsos_muestra, kwargs_procesamiento)
        
        if self.pulsos_fondo is not None:
            self.areas_fondo = self._procesar_matriz(self.pulsos_fondo, kwargs_procesamiento)
        else:
            self.areas_fondo = []

    def exportar_csv(self, ruta):
        if not self.areas_acumuladas:
            raise ValueError("No hay datos integrados para exportar.")
            
        df = pd.DataFrame({"area_V_ns": self.areas_acumuladas})
        if self.areas_fondo:
            # Rellenar con NaN si las longitudes son distintas
            df_fondo = pd.DataFrame({"area_fondo_V_ns": self.areas_fondo})
            df = pd.concat([df, df_fondo], axis=1)
            
        df.to_csv(ruta, index=False)