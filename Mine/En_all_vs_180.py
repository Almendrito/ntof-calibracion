import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. Parámetros y Constantes Físicas
# ==========================================
m_n  = 1.008664
m_d  = 2.014101
m_He = 3.016029
Q    = 3.269 

# ==========================================
# 2. Función de la Cinemática
# ==========================================
def calcular_En(E_d_keV, angulo_deg):
    """Calcula la energía del neutrón en MeV dado E_d en keV y el ángulo en grados."""
    E_d_MeV = E_d_keV / 1000.0  # Conversión a MeV
    theta = np.radians(angulo_deg)
    
    factor_externo = E_d_MeV * (m_n * m_d) / ((m_n + m_He)**2)
    
    num = m_He * (m_n + m_He - m_d) * E_d_MeV + m_He * (m_He + m_n) * Q
    den = m_n * m_d * E_d_MeV
    
    interior_raiz = (num / den) - np.sin(theta)**2
    interior_raiz_seguro = np.where(interior_raiz >= 0, interior_raiz, np.nan)
    
    E_n = factor_externo * (np.cos(theta) + np.sqrt(interior_raiz_seguro))**2
    return E_n

# ==========================================
# 3. Cálculo de la Línea Base y Diferencias
# ==========================================
E_d_keV = np.linspace(0.1, 1000, 1000)

# Calculamos la energía a 180° que servirá como nuestra línea base
E_n_180 = calcular_En(E_d_keV, 180)

# Definimos los ángulos que queremos comparar contra 180°
angulos_a_comparar = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165]

# ==========================================
# 4. Formato y Visualización
# ==========================================
plt.figure(figsize=(8, 6))

# Iteramos sobre los ángulos, calculamos y restamos la base (180°)
for angulo in angulos_a_comparar:
    E_n_actual = calcular_En(E_d_keV, angulo)
    diferencia = E_n_actual - E_n_180
    
    plt.plot(E_d_keV, diferencia, label=f'{angulo}° - 180°')

# Configuraciones del gráfico
plt.title('Diferencia de energía respecto a la emisión en 180°', fontsize=16)
plt.xlabel('Deuteron energy $E_d$ (keV)', fontsize=16)
plt.ylabel('$\Delta E_n = E_n(\\theta) - E_n(180^\circ)$ (MeV)', fontsize=16)

plt.xlim(0, 1000)
plt.ylim(ymin=0) 

plt.tick_params(axis='both', direction='in', top=True, right=True)
plt.legend(title='Comparación', frameon=False, loc='upper left')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.show()