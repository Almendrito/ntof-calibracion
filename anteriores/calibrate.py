import re

import matplotlib.pyplot as plt 
import numpy as np  


archivo = 'C:\\Users\\mateo\\OneDrive\\Escritorio\\A BG 5GaS 10k 02mV trigg att.xml'  # Replace with your actual file path
# Load the data from the text file
with open(archivo, 'r',encoding = 'utf-8') as f:
    contenido = f.read()

try:
    # 3. Pescamos el texto
    tiempo_texto = re.search(r'<Time[^>]*>(.*?)</Time>', contenido, re.DOTALL).group(1)
    voltaje_texto = re.search(r'<Data[^>]*>(.*?)</Data>', contenido, re.DOTALL).group(1)

    # 4. Limpieza a prueba de balas: 
    # Cambiamos las comas por espacios, separamos los números, y los forzamos a ser flotantes (decimales).
    # Esto ignora cualquier salto de línea, espacios extra o comas finales.
    tiempo = np.array([float(x) for x in tiempo_texto.replace(',', ' ').split()])
    voltaje = np.array([float(x) for x in voltaje_texto.replace(',', ' ').split()])

    # 5. ¡A graficar!
    plt.figure(figsize=(8, 5))
    
    # Truco extra: a veces el DRS4 guarda el voltaje invertido o en distintas unidades.
    # Si ves que el pulso va hacia abajo y quieres verlo positivo, puedes poner "-voltaje".
    plt.plot(tiempo, voltaje, color='tab:green', linewidth=1.5, label='Primer Pulso Extraído')
    
    plt.title('Pulso extraído del archivo XML del DRS4', fontsize=14)
    plt.xlabel('Tiempo (ns)', fontsize=12)
    plt.ylabel('Amplitud (mV)', fontsize=12)
    
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.grid(True, alpha=0.4)
    plt.legend()
    
    plt.show()

except AttributeError:
    print("Oops. No se encontraron las etiquetas <Time> o <Data> en el archivo.")
except ValueError as e:
    print(f"Error al convertir los números: {e}")