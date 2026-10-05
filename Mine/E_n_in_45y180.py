import matplotlib.pyplot as plt
import numpy as np


#Definir las constantes físicas e intervalos
m_n_SI = 1.675e-27  # mass of neutron in kg
m_d_SI = 3.344e-27  # mass of deuteron in kg

m_n_uma = 1.00866491588  # mass of neutron in atomic mass units
m_d_uma = 2.013553212745  # mass of deuteron
m_He3_uma = 3.016029320729  # mass of helium-3 in atomic mass units
Q = 3.27

E_d_Mev = np.linspace(0.1, 1000, 1000)  # Deuteron energy in MeV
E_d = E_d_Mev/1000  # Deuteron energy in MeV

angulos_grados = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]  # Angles in degrees


#Definir función para calcular la energía del neutrón en función de la energía del deuterón y el ángulo

def calcular_energia_neutron(E_d, angulo):
    theta_rad = np.radians(angulo)  # Convert angle to radians

    factor_constante = E_d * (m_n_uma / m_d_uma)/(m_He3_uma + m_n_uma)**2 

    fraccion_raiz_num = m_He3_uma*(m_n_uma + m_He3_uma - m_d_uma)*E_d + m_He3_uma*(m_He3_uma + m_n_uma)* Q 
    fraccion_raiz_den = m_n_uma*m_d_uma*E_d
    fraccion_raiz = fraccion_raiz_num / fraccion_raiz_den
    factor_variable =(np.cos(theta_rad) + np.sqrt(fraccion_raiz - np.sin(theta_rad)**2))**2  # Variable factor for energy calculation

    E_n = factor_constante * factor_variable  # Neutron energy in MeV
    return E_n

E_n_45 = calcular_energia_neutron(E_d, 45)  # Neutron energy for 45 degrees
E_n_180 = calcular_energia_neutron(E_d, 180)  # Neutron energy for 180 degrees

diferencia_energia = E_n_45 - E_n_180  # Difference in neutron energy between 45 and 180 degrees


plt.figure(figsize=(10, 6))

plt.plot(E_d, diferencia_energia, label='Difference in Neutron Energy', color='red')
plt.xlabel('Deuteron Energy (MeV)', fontsize = 16)
plt.ylabel('Difference in Neutron Energy (MeV)', fontsize = 16)
plt.title('Energy Distribution of Neutrons', fontsize = 18)
plt.xlim(0, 1)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()