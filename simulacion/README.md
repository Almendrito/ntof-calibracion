# ComptonSpectrum

Simulación Geant4 11.4 del espectro Compton medido en un centellador plástico
**BC408** (base PVT), con análisis en Python.

## Física

El BC408 es de Z bajo (solo H y C), así que el efecto fotoeléctrico es
despreciable a cientos de keV: el gamma interactúa casi solo por **dispersión
Compton**. No hay fotopico de energía total, sino un **continuo Compton** que
termina en un **borde**. El electrón de retroceso deposita `T = Eγ − Eγ'`, con
máximo en la retrodispersión (180°):

```
T_max = 2·Eγ² / (m_e c² + 2·Eγ),    m_e c² = 0.511 MeV
```

| Fuente | Líneas (keV)   | Bordes Compton (keV) |
|--------|----------------|----------------------|
| Co-60  | 1173.2, 1332.5 | 963, 1118            |
| Cs-137 | 661.66         | 477                  |
| Na-22  | 511, 1274.5    | 341, 1062            |

El espectro **medido** difiere del ideal por: resolución del detector (redondea
el borde, efecto dominante en plástico), **scatter Compton múltiple** (rellena
por encima del borde), Doppler atómico (ensancha el borde, lo captura el modelo
LowEP de EM option4) y retrodispersión en el material circundante.

### De energía depositada a luz emitida (Birks)

El PMT integra **fotones**, no energía. Entre ambos está el quenching de Birks:
a lo largo de una traza muy ionizante, las moléculas excitadas se desexcitan sin
emitir, así que la luz por unidad de camino satura:

```
dL/dx = (dE/dx) / (1 + kB·dE/dx)
```

Esto hace que la respuesta del centellador sea **no lineal**: `L/E` depende de la
energía. Medido en esta simulación para el BC408 con `kB = 0.126 mm/MeV`:

| E electrón [keV] | 10 | 20 | 50 | 100 | 200 | 300 | 477 | 700 | 1000 | 1500 |
|---|---|---|---|---|---|---|---|---|---|---|
| L/E | 0.658 | 0.767 | 0.839 | 0.884 | 0.918 | 0.934 | 0.947 | 0.954 | 0.960 | 0.964 |

Consecuencias directas para la calibración:

- El borde Compton del Cs-137 cae en **463 keVee** de luz, no en 477 keV de
  energía. Calibrar `canal ↔ energía` contra `Edep` mete ese sesgo en la ganancia.
- Buena parte de la curvatura que las formas `quad`/`exp` de `calmodels.py`
  estaban absorbiendo **es Birks**, no el PMT. Ahora se predice en vez de ajustarse.
- El quenching se aplica **paso a paso**, nunca sobre el `Edep` total del evento:
  es no lineal en `dE/dx`, así que en un evento de scatter múltiple (varios
  electrones de baja energía) quenchar la suma da un resultado distinto y erróneo
  frente a sumar los quench individuales.

**`kB` no viene gratis de Geant4.** Su tabla interna
(`G4EmSaturation::InitialiseG4materials`) solo cubre 4 materiales — poliestireno,
BGO, lAr y PbWO4 — y deja `G4_PLASTIC_SC_VINYLTOLUENE` explícitamente comentado.
Sin fijar `kB` a mano el material queda con `kB = 0` y **no hay quenching en
absoluto**. Se fija en `DetectorConstruction::ApplyBirks()` y se ajusta con
`/det/birks`. El valor por defecto, 0.126 mm/MeV (= 0.0130 g/cm²/MeV), es el
típico de literatura para plásticos de la familia BC-400; **confírmalo contra tu
propio detector** si necesitas precisión absoluta, porque `kB` es el único
parámetro libre del modelo de luz.

## Diseño

- **Material:** `G4_PLASTIC_SC_VINYLTOLUENE` (base PVT del BC408, 1.032 g/cm³).
- **Geometría:** cilindro 2"×2" (r 2.54 cm, largo 5.08 cm) + carcasa de Al
  (pared lateral 10 mm, tapa frontal 6 mm hacia la fuente) + stub de PMT de
  vidrio. Grosores y volúmenes configurables por macro. La tapa frontal de Al
  atenúa el haz que entra y agrega dispersión/retrodispersión.
- **Física:** `G4EmStandardPhysics_option4` (Compton LowEP con perfil Doppler).
- **Scoring:** por evento, energía depositada (`Edep_keV`) **y luz emitida**
  (`Lvis_keVee`, con quenching de Birks paso a paso vía `G4EmSaturation`). La
  **resolución se aplica en el análisis** (convolución Gaussiana), así se puede
  re-suavizar sin re-simular.
- **Fuente:** Co-60 por defecto, gammas monoenergéticos limpios desde una fuente
  puntual. Muestrea una línea por evento; con `/gun/cascade true` emite ambas
  por evento (summing por coincidencia).

## Compilar

```bash
mkdir -p build && cd build
source /home/mateo/Software/geant4-install/bin/geant4.sh
cmake .. && make -j$(nproc)
```

## Simular (batch, multihilo)

```bash
cd build
./compton_spectrum co60.mac      # -> ComptonSpectrum_Co60_2M.root
./compton_spectrum cs137.mac
./compton_spectrum na22.mac
```

El nombre del archivo lleva los **primarios lanzados** en la corrida
(`ComptonSpectrum_<fuente>_<N>.root`): `3M`, `2.5M`, `20k`, `500`. Ese N es el
`beamOn`, **no** el número de filas del ntuple, que solo guarda los eventos con
depósito > 0. Con la geometría por defecto y haz `pencil`, alrededor del 35% de
los gammas de Cs-137 depositan algo, así que 3M primarios dan ~1,06M de filas.

El número tiene que salir igual en el master y en todos los hilos, o cada uno
abriría un archivo distinto y el merge del ntuple se rompería; por eso se toma de
`GetNumberOfEventToBeProcessed()` (el total de la corrida, que el worker hereda)
y no de los eventos que a cada hilo le tocaron.

Como consecuencia ya no hay un nombre fijo al que apuntar: las GUIs buscan
automáticamente el `.root` más reciente en `build/`, y los scripts de línea de
comandos aceptan comodines (`'../build/ComptonSpectrum_Cs137*.root'`, entre
comillas).

## Visualizar (WSL2)

```bash
cd build
./run_vis.sh
```

## Analizar y graficar

```bash
cd analysis
python3 plot_spectrum.py '../build/ComptonSpectrum_Co60*.root' --preset co60
```

Genera `compton_co60.png` con el espectro verdadero, el medido y la curva
Klein-Nishina, y marca los bordes. Opciones útiles:

- `--sigma KEV` : una sola Gaussiana de **sigma constante** (igual que la GUI).
  Sin `--sigma` usa el modelo dependiente de la energía `--res-a/--res-b/--res-c`
  con `FWHM/E = sqrt(a² + b²/E + c²/E²)`.
- `--lines 1173.2,1332.5 --weights 1,1` : líneas arbitrarias para la curva KN.
- `--bins --emax --output`.

## Visor de espectros (`spectrum_viewer.py`)

```bash
cd analysis
python3 spectrum_viewer.py                                   # y abres el archivo desde la GUI
python3 spectrum_viewer.py ../build/ComptonSpectrum_Cs137_3M.root   # o lo pasas directo
```

El visor simple: eliges uno o varios `.root` y los grafica como **cuentas vs
energía**. Es el de *mirar*, no el de calibrar; no ajusta resolución ni calibra
canales (para eso están las dos GUIs siguientes).

- **Varios archivos superpuestos**, un color por archivo: para comparar fuentes
  entre sí.
- **Observable**: `luz (Birks) [keVee]` (por defecto), `deposito [keV]`, o
  **`ambas`** superpuestas — línea sólida la luz, punteada el depósito. Esta
  última es la forma rápida de ver cuánto comprime el quenching: en Cs-137 el
  borde pasa de ~477 keV a ~450 keVee a simple vista.
- **Ejes**: número de bins, rango `E min`/`E max`, escala log en Y (con el límite
  inferior anclado al bin no vacío más bajo, para no comerse la cola del
  espectro), y trazo en escalones o continuo.
- **Normalización**: `cuentas` crudas, `por evento`, o `área = 1` — para comparar
  archivos con estadística distinta sin re-simular.
- **Marcas**: bordes Compton y energías de retrodispersión del preset, que se
  autodetecta por el nombre del archivo.
- **Guardar** la figura (PNG/PDF) o exportar el histograma a CSV.

## Interfaz interactiva (ajustar resolución a datos)

```bash
cd analysis
python3 gui_spectrum.py
```

GUI Tkinter para **ajustar en vivo el ancho** (resolución) del espectro simulado
y compararlo con tu espectro medido:

- Carga el `.root` simulado (autodetecta la fuente por el nombre del archivo) y
  elige la **observable** (`light` por defecto = `Lvis_keVee`, la luz con Birks;
  `edep` = depósito crudo). Ajustar la resolución contra la luz es lo correcto:
  el ensanchamiento que mides en tus datos es el de la carga del PMT.
- **Slider de sigma (keV)**: el detector convierte cada depósito en un pulso, así
  que la señal medida es la señal verdadera convolucionada con una sola Gaussiana
  de ancho constante. Al mover el slider se re-pliega el espectro al instante
  (matriz de respuesta que conserva las cuentas, sin ruido Monte Carlo). El panel
  muestra la FWHM (= 2.355·sigma) y la FWHM/E en cada borde Compton.
- **Carga tu espectro medido** (CSV o texto de dos columnas). Si está en
  canales, calíbralo con `gain` (keV/canal) y `offset`; escala vertical auto o
  manual.
- **Ajustar sigma a los datos**: auto-ajuste por mínimos cuadrados (scipy) del
  sigma para que el sim coincida con tu medido, desde una energía mínima
  configurable (evita la zona de umbral).
- Guarda la figura o exporta el espectro medido simulado a CSV.

`analysis/example_measured.csv` es un espectro de demostración sintético (Co-60,
sigma = 45 keV, calibrado a 2 keV/canal) para probar la carga, calibración y ajuste.

## Calibrar canal -> energía contra los datos del DRS4

```bash
cd analysis
python3 calibrate.py --sim '../build/ComptonSpectrum_Co60*.root' \
    --measured "/mnt/c/Users/mateo/OneDrive/Escritorio/B Co60 3k 1500V 0025 trigg att 20Db.xml" \
    --background "/mnt/c/Users/mateo/OneDrive/Escritorio/B BG 3k 1500V 0148 trigg att 20Db.xml" \
    --preset co60 --observable charge --channel 2
```

`calibrate.py` cierra el lazo entre la simulación y la medición real: el espectro
medido del DRS4 está en **carga (pC)** o **amplitud (mV)**, no en keV. Como en un
centellador plástico la luz (y por tanto la carga del pulso) es ~lineal con la
energía depositada, el canal se relaciona con la energía por una recta
`canal = gain·E + offset`. El script:

1. Construye el espectro medido desde el archivo DRS4 reusando `drs_spectrum.py`
   (histograma de carga/amplitud, resta de fondo opcional). También acepta un CSV
   de dos columnas (canal, cuentas) ya exportado, con `--measured archivo.csv`.
2. Pliega la simulación con la resolución del detector (Gaussiana de sigma).
3. Ajusta por mínimos cuadrados los parámetros de la calibración `canal<->E` y
   `sigma` (y el offset con `--free-offset`) para que la sim plegada coincida con
   el medido; el borde Compton es el ancla. La **forma** de la calibración se
   elige con `--cal-form linear|quad|exp` (ver más abajo).
4. Imprime la calibración `canal <-> energía`, el canal donde cae cada borde
   Compton, convierte canales arbitrarios (`--query 12,25,40`) y dibuja la
   comparación con doble eje (canal abajo, keV arriba). Exporta con `--export-csv`.

Opciones útiles: `--observable charge|amp`, `--channel`, `--polarity`,
`--norm real|live|none` (el XML de drsosc no guarda scaler por canal, así que
`live` suele caer a cuentas crudas; por eso el default es `real`), `--sigma0`,
`--gain0`, `--fit-cmin/--fit-cmax` (ventana de ajuste en canal), `--drs-dir`
(carpeta de `drs_spectrum.py`, por defecto el Escritorio). Empareja siempre una
medición y un fondo tomados con el **mismo trigger y atenuación**.

### GUI de calibración

```bash
cd analysis
python3 calibrate_gui.py
```

`calibrate_gui.py` es la misma calibración con interfaz Tkinter: seleccionas con
botones la simulación, el medido (DRS4 o CSV) y el fondo; ajustas observable,
canal, polaridad y ventana de integración; pulsas **Procesar medido** (el parseo
del XML grande corre en un hilo aparte para no congelar la ventana) y **Ajustar /
calibrado**. Eliges la **forma** de la calibración (lineal/quad/exp), el slider de
`sigma` re-dibuja en vivo, el cuadro de resultado muestra los parámetros y el RMS
del ajuste, y otro cuadro convierte cualquier canal a energía al instante. Exporta
el CSV calibrado o la figura.

Dibuja tres curvas (cada una con su casilla para mostrar/ocultar): **Medido**
(azul), **Depósito verdadero MC** (gris, el espectro simulado con el borde Compton
nítido, sin resolución) y **Sim plegada** (rojo, la simulación convolucionada con
la resolución, que es la que se compara con el medido). Las tres comparten la
misma calibración canal->keV y el mismo factor de escala. El CLI dibuja las mismas
tres curvas.

### Forma de la calibración canal <-> energía

La relación canal->energía se modela con `calmodels.py`, que ofrece tres formas
(todas invertibles, todas se reducen a la recta cuando la curvatura tiende a 0):

| forma    | `canal = f(E)`                  | cuándo |
|----------|---------------------------------|--------|
| `linear` | `a·E + b`                       | caso por defecto, estándar para centellador a HV fijo |
| `quad`   | `a·E² + b·E + c`                | curvatura suave (p.ej. saturación leve del PMT) |
| `exp`    | `off + G·(e^{B·E} − 1)/B`       | curvatura exponencial; `B→0` es la recta |

Se elige con `--cal-form` (CLI) o el combo "forma canal<->E" (GUI). El ajuste y
el reporte (RSS/RMS del residuo en la ventana) son idénticos para las tres.

**Normalización de intensidad.** No hace falta igualar tiempo de toma ni número de
eventos simulados: el ajuste aplica en cada paso un factor vertical
`A = Σ(modelo·y)/Σ(modelo·modelo)` (mínimos cuadrados) que pone la simulación a la
altura del medido. Además, al proyectar la sim de energía a canal se incluye el
**jacobiano `dE/dc`** (`project_to_channel`), de modo que la intensidad queda en la
misma base de densidad que el histograma medido. Para una recta `dE/dc = 1/gain`
es constante y lo absorbe `A` (no cambia nada); para `quad`/`exp` varía a lo largo
del espectro y corrige el sesgo de comparar histogramas con anchos de bin
desiguales.

**Ojo**:
con el espectro de **una sola fuente** no se puede distinguir la forma (el borde
único se confunde con la resolución: `gain` y `sigma` están correlacionados). Para
decidir de verdad si tu PMT a 1500 V es lineal o no, usa `linearity.py`.

### Resolución (sigma): ajustarla o fijarla

El `sigma` ensancha la simulación tanto como el detector ensancha tus datos. **No
sale del ancho temporal del pulso** (p.ej. 7 ns): ese ancho es la *velocidad* del
pulso (decaimiento del BC408, tránsito del PMT, electrónica) y no se convierte a
keV. La resolución en energía viene de la fluctuación de la *carga* (estadística de
fotoelectrones), y en un plástico (sin fotopico) se mide del **flanco del borde
Compton**: por eso ajustar sigma al borde es la forma de medirla, no un defecto.

Como la resolución es una propiedad del detector, puedes medirla una vez y luego
**fijarla**: ajusta sigma en una fuente limpia (Cs-137), lee el valor, y para las
demás fuentes usa `--fix-sigma --sigma0 <valor>` (CLI) o la casilla "fijar sigma"
(GUI), de modo que el ajuste mueva solo la calibración y deje sigma constante.
(Donde sí entran tus 7 ns es en la ventana de integración `--pre-ns/--post-ns`, que
debe cubrir todo el pulso para no perder carga y empeorar la resolución.)

### Curva de respuesta de luz (`light_response.py`)

```bash
# 1) barrido de electrones monoenergeticos (banco de pruebas de la respuesta pura)
cd build
for E in 10 20 50 100 200 300 477 700 1000 1500; do
  printf "/gun/particle e-\n/gun/energy $E keV\n/run/initialize\n/run/beamOn 20000\n" > e.mac
  ./compton_spectrum e.mac
done
cd ../analysis
python3 light_response.py curve '../build/ComptonSpectrum_E*keV.root'

# 2) cuanto comprime la luz al espectro Compton de una fuente real
python3 light_response.py edge ../build/ComptonSpectrum_Cs137_3M.root --preset cs137
```

`curve` dispara electrones desde el **centro** del centellador, donde se frenan
enteros (alcance de un e- de 1 MeV en PVT ≈ 4 mm, frente a 25.4 mm de radio), así
que `Lvis/E0` es la respuesta de luz pura, sin la mezcla de scatters múltiples del
espectro Compton.

`edge` recupera la transformación `E → L` por **emparejamiento de cuantiles**: como
el quenching es monótono, el cuantil `q` de la luz es la imagen del cuantil `q` de
la energía, así que la curva sale sin tener que localizar ningún borde. Esto
importa en Co-60, donde los dos bordes (963 y 1118 keV) se solapan y ningún
buscador de flancos los separa; ahí el script lo dice en vez de dar un número
falso. Las dos rutas coinciden (0.881 vs 0.884 a 100 keV; 0.945 vs 0.947 a ~500
keV), lo que valida el modelo de luz.

Dos advertencias que salen de estas medidas:

- El borde medido en `Edep` cae por **encima** del borde teórico (496 vs 477 keV en
  Cs-137). No es un error: el scatter múltiple rellena por encima del corte
  cinemático. Por eso conviene ajustar la **forma** simulada completa y no clavar
  el valor nominal del borde.
- El borde en luz cae en 463 keVee, un **-6.6%** respecto al mismo borde en energía.

### Chequeo de linealidad multi-fuente (`linearity.py`)

La ganancia del PMT (XP2262/B) es una ley de potencias muy empinada **del alto
voltaje**, pero a HV fijo eso es una constante: fija la *pendiente* de la recta
canal->energía, no su *curvatura*. Las desviaciones reales de la recta vienen de
la ley de Birks (baja energía) o de la saturación del PMT (alta energía), y se
comprueban con varios bordes Compton conocidos.

> **Con el modelo de luz activo, Birks ya no es un efecto que haya que ajustar
> aquí.** Si comparas los canales contra `Lvis_keVee` (el default), el quenching
> ya está dentro del modelo y cualquier curvatura que sobre apunta al PMT o a la
> electrónica, no al centellador. Contra `Edep_keV` los dos efectos se mezclan y
> el test pierde poder. Ajusta los bordes en **keVee**, no en keV.

```bash
# 1) localiza el canal del borde en cada fuente (independiente del modelo)
python3 linearity.py locate "Cs137.xml" --edge 477 --observable charge --channel 2
python3 linearity.py locate "Na22.csv"  --edge 1062
python3 linearity.py locate "Co60.csv"  --edge 1040
# 2) reúne los pares (E_keV  canal  etiqueta) en puntos.txt y ajusta las 3 formas
python3 linearity.py fit puntos.txt --output linealidad.png
```

`locate` halla el canal del borde como el punto de media altura del flanco de
bajada (no depende de ningún modelo). `fit` ajusta las tres formas a los puntos
`(E_borde, canal_borde)`, imprime el RMS del residuo de cada una y grafica
`canal vs E` con los residuos respecto a la recta: si los puntos se curvan, `quad`
o `exp` bajan mucho el residuo (hay no linealidad); si no, `linear` basta.
Necesitas al menos 4-5 bordes (Cs-137: 477; Na-22: 341 y 1062; Co-60: ~1040) para
que la comparación sea informativa.

**Rango de comparación**: arrastra horizontalmente sobre el gráfico (o escribe en
`canal min`/`canal max`) para elegir la zona donde se comparan medición y sim. El
rango queda sombreado, el escalado vertical y el ajuste de `gain`/`sigma` se
restringen a él (útil para anclar al borde Compton e ignorar el umbral de bajo
canal o la cola de scatter múltiple). El botón **Rango completo** lo limpia. En el
CLI es lo mismo con `--fit-cmin/--fit-cmax`.

## Comandos UI

```
/det/radius 2.54 cm      /det/halfZ 2.54 cm      # geometría (antes de /run/initialize)
/det/pmt true            /det/housing true
/det/housingWall 10 mm   /det/housingFront 6 mm  # grosor pared lateral / tapa frontal de Al
/det/birks 0.126                                  # constante de Birks en mm/MeV (0 = sin quenching)
/gun/preset co60|cs137|na22                       # fuente
/gun/energy 661.66 keV                            # una línea arbitraria
/gun/cascade true|false                           # summing por coincidencia
/gun/distance 10 cm                               # fuente -> cara frontal
/gun/beam pencil|cone|iso                         # divergencia del haz
/gun/particle gamma|e-                            # e-: banco de pruebas de la curva de luz L(E)
```

## Salida ROOT

Ntuple `ComptonSpectrum`, una fila por evento con depósito > 0:

| columna | significado |
|---|---|
| `EventID` | id del evento |
| `Edep_keV` | energía depositada en el centellador |
| `E0_keV` | energía primaria total del evento |
| `NCompton` | dispersiones Compton del fotón en el centellador |
| `NPhoto` | absorciones fotoeléctricas (~0 en plástico de Z bajo) |
| `FullAbs` | 1 si `Edep ≈ E0` (1%) |
| `Lvis_keVee` | **luz emitida**: energía visible tras Birks, sumada paso a paso |
| `QF` | factor de quenching del evento, `Lvis/Edep` (≤ 1) |

Histogramas `Edep` y `Lvis` (0–2600 keV) de conveniencia. Se lee con `uproot`
(no requiere instalar ROOT).

`Lvis_keVee` es la observable que hay que comparar con la carga del DRS4: el PMT
ve luz. Por eso `calibrate.py`, `calibrate_gui.py` y `plot_spectrum.py` usan la
luz **por defecto**, con `--sim-observable edep` / `--observable edep` para volver
al depósito crudo y aislar cuánto de la no linealidad viene de Birks. Los ROOT
generados antes de este cambio no tienen la columna: los scripts avisan y caen a
`Edep_keV`, pero conviene re-simular.
