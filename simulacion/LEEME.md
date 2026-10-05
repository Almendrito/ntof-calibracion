# La simulación de Geant4 que produjo los `.root` de `datos/`

Copia (21-sep-2026) de la fuente de `~/Proyectos/Geant4-Setup/ComptonSpectrum`
en WSL, que es la que generó `datos/ComptonSpectrum_Co60_3M.root` y
`datos/ComptonSpectrum_Cs137_10M.root`. Está acá para que el informe la lea
(`analisis/informe.geant4_geometria` saca de `src/DetectorConstruction.cc`
y `src/PhysicsList.cc` la geometría, el material, la constante de Birks y el
corte, sin tipearlos) y para que quien tenga que rehacerla o simular otra
fuente tenga el código al lado de los datos.

**Ojo con la otra carpeta.** En `JcV/ComptonSimulation/` (Windows) hay OTRA
simulación, más vieja: escribe un solo histograma en CSV, no tiene carcasa
ni PMT, emite iones de Cs-137 con decaimiento radiactivo y no guarda
`NCompton` ni `Lvis`. No es la que produjo los `.root` y el informe no la
lee. Hasta el 21-sep-2026 `informe.geant4_geometria` leía ESA por error;
los cuatro números que sacaba (densidad, Birks, radio, semialtura) coinciden
en las dos, así que no cambió ningún resultado, pero la descripción de la
geometría que se podía escribir a partir de ella estaba incompleta.

## Qué hay

| archivo | qué es |
|---|---|
| `src/DetectorConstruction.cc` | geometría, material, Birks, y la configuración de la fuente (mensajeros `/det/` y `/gun/`) |
| `src/PhysicsList.cc` | `G4EmStandardPhysics_option4` + `G4DecayPhysics`, corte de 0.1 mm |
| `src/PrimaryGeneratorAction.cc` | gammas monoenergéticos desde una fuente puntual en el eje; una línea por evento, sorteada por peso |
| `src/SteppingAction.cc` | acumula `Edep` y la luz con Birks paso a paso; cuenta `compt` y `phot` del fotón |
| `src/EventAction.cc`, `src/RunAction.cc` | el ntuple `ComptonSpectrum` (una fila por evento con depósito) y el nombre del archivo, que lleva el número de primarios |
| `co60.mac`, `cs137.mac`, `na22.mac` | los presets; el `beamOn` de los `.root` de `datos/` fue 3M y 10M (el nombre del archivo lo escribe el programa a partir del `beamOn`, no del macro) |
| `README.md` | el README original de WSL, con la tabla de `L/E` medida con electrones monoenergéticos y las herramientas de análisis de allá (que este proyecto NO usa: acá la plantilla la lee `analisis/simulacion.py`) |

Lo que NO se copió: `build/` (el ejecutable y los `.root`, que ya están en
`datos/`) y `analysis/` (scripts de Python de la etapa anterior al proyecto
`Calibración CE`).

## Compilar y correr (en WSL, Geant4 11.4.0)

```bash
source /home/mateo/Software/geant4-install/bin/geant4.sh
mkdir -p build && cd build && cmake .. && make -j$(nproc)
./compton_spectrum co60.mac        # -> ComptonSpectrum_Co60_<N>.root
./compton_spectrum cs137.mac
```

Los comandos de UI (`/det/radius`, `/det/housing`, `/gun/preset`,
`/gun/beam pencil|cone|iso`, `/gun/distance`, `/det/birks`, ...) están
listados en `README.md`, sección «Comandos UI». Los `.root` de `datos/` se
corrieron con los valores por defecto: haz `pencil`, fuente a 10 cm de la
tapa, carcasa y PMT construidos, kB = 0.126 mm/MeV.
