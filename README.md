# ntof-calibracion

Calibración en energía de los detectores de centelleo plástico (BC-408 con PMT XP2262/B) del diagnóstico nToF de Llampüdkeñ, a partir del borde Compton. Implementa el método NCC-SA de Kong et al. (2025): un modelo del espectro (Klein-Nishina convolucionado con una gaussiana), la selección automática del intervalo por correlación cruzada normalizada y el ajuste global por recocido simulado.

Resultado vigente del detector A con Co-60 a 1550 V: D = 1.251 ± 0.005 mV·ns/keVee, medido a la entrada del DRS4 con atenuador de 20 dB. El detalle y su historia están en `CLAUDE.md` y en `informe/`.

## Qué hay en cada carpeta

| Carpeta | Qué contiene |
| --- | --- |
| `analisis/` | Implementación vigente: lectura del DRS4 (XML y binario), pulso, área, espectro, NCC, recocido simulado, calibración e interfaz (`gui.py`) |
| `desde_cero/` | Reimplementación paso a paso con sus verificaciones |
| `replica_paper/` | Réplica literal del paper de Kong et al. y auditoría de sus erratas |
| `codigo/` | Implementación anterior, se conserva como referencia |
| `pulsos/` | Estudio del ancho y la forma de los pulsos |
| `verificaciones/` | Pruebas de partes del análisis |
| `simulacion/` | Copia de la simulación Geant4 del repositorio [ComptonSpectrum](https://github.com/Almendrito/ComptonSpectrum), la que generó las plantillas `.root`. Se queda aquí porque `analisis/informe.py` lee de ella la geometría y la constante de Birks. Los macros `co60.mac`, `cs137.mac` y `na22.mac` solo están en esta copia |
| `informe/` | Informe de calibración en LaTeX |
| `docs/` | Campaña de linealidad, formato binario del DRS4 y documentación del recocido simulado |
| `Mine/` | Energía del neutrón DD a 45° y 180° |
| `anteriores/` | Primeros análisis de abril de 2026: `Analizador/` (interfaz de espectros con muestra y fondo) y `calibrate.py` (lectura del XML del DRS4) |

## Datos

Los datos (15 GB de mediciones del DRS4 y las simulaciones `.root`) no están en este repositorio. Se quedan en OneDrive y el código los busca en la carpeta `datos/`. Para conectarla, una sola vez, en PowerShell:

```powershell
cd C:\Users\mateo\Repos\ntof-calibracion
New-Item -ItemType Junction -Path datos -Target "$env:OneDriveCommercial\Weas de la U\Física Nuclear\JcV\Calibración CE\datos"
```

Eso crea un acceso de carpeta (junction) que no ocupa espacio ni necesita permisos de administrador. Git lo ignora.

## Cómo correrlo

Necesita Python 3 con numpy, scipy y matplotlib. tkinter viene incluido en Python para Windows.

```powershell
python analisis/rutas.py     # muestra qué corridas encuentra en datos/
python analisis/gui.py       # interfaz de calibración
```

Los comandos de cada etapa están en `analisis/README.md`. Las figuras quedan en `analisis/figuras/` y no se suben a GitHub porque se regeneran.

## Referencia

Kong et al., "The Automatic Calibration Method of the Compton Edge Based on Normalized Cross-correlation and Simulated Annealing Algorithm", IEEE Transactions on Nuclear Science (2025). [arXiv:2507.16282](https://arxiv.org/abs/2507.16282)
