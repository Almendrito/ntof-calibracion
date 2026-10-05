# Plan del texto guía de la calibración

Este archivo es la lista de lo que `informe.tex` tiene que llegar a contener
para servir como **texto guía**: el documento que alguien que entra al
laboratorio lee para entender qué es calibrar un centellador plástico por
borde Compton, cómo se hizo acá, y con qué física y qué estadística se
sostiene cada paso. Se va completando sesión a sesión: cada bloque tiene
casillas, y una sesión de trabajo es uno o dos bloques, no más.

Fecha de este inventario: 20-sep-2026 (bloque D marcado el 21-sep-2026). `informe.tex` tiene 1842 líneas,
13 secciones y 4 apéndices, 13 figuras y 10 tablas generadas.

**Lector objetivo.** Estudiante de física al final del pregrado o al inicio de
una tesis, que sabe mecánica cuántica básica y estadística de primer curso,
que nunca vio un DRS4 ni una simulación de Geant4, y que va a tener que
repetir o extender esta calibración. El texto tiene que servirle sin tener que
abrir `CLAUDE.md` ni el código.

**Cómo se usa este archivo.**
- Antes de escribir, elegir un bloque de la sección 3 (el orden es una
  propuesta, no una obligación) y leer su lista de la sección 2.
- Al terminar, marcar las casillas y anotar la fecha al lado del bloque.
- Si al escribir aparece algo que no estaba en la lista, agregarlo acá antes
  de escribirlo en el `.tex`, para que el plan siga siendo el índice.
- Las reglas de la sección 4 valen para todo lo que se agregue.

---

## 0. Reestructuración del 24-sep-2026 (leer primero)

> La estadística que usa el código y falta en el informe tiene su propio plan,
> por sesiones: `PLAN_estadistica.md` (25-sep-2026).

El informe se reordenó para que cuente **qué se hizo, con qué lógica y cómo
se inspira en Kong et al.** Respaldo del anterior: `informe_v1_2026-09-24.tex`
y `.pdf`. Orden nuevo del cuerpo:

1. Introducción (problema, el NCC-SA en un párrafo, qué hace este trabajo,
   resultados en un párrafo, cómo leer + convenciones y etiquetas).
2. **El método de referencia y cómo lo adaptamos** (`sec:metodoref`): qué
   critica el paper, NCC-SA en 7 pasos con sus ecuaciones, `fig:flujo`
   (`flujo.tex`), `tab:adaptacion`, por qué la plantilla simulada, lo que el
   paper no trata.
3. **Montaje experimental y campaña** (`sec:montaje`): BC-408 + XP2262/B
   (`tab:hardware`), ley de ganancia, cadena de señal (`fig:cadena-senal`),
   atenuador, sesiones, lo no registrado.
4. La física del borde (incluye `sec:quenching`, nuevo).
5. Del pulso al espectro. 6. Geant4. 7. El ajuste (NCC, SA, barra,
   decisiones). 8. Resultados (calibraciones, ganancia con k = b/12, test A,
   detector D). 9. Resultados frente al paper. 10. **Conclusiones** (nuevo).
   11. Uso. 12. Limitaciones (sólo las abiertas).

Apéndices: A NCC/SA en detalle · B estadística y cierre · C controles del
observable · D clásicos medidos · E simulación · F auditoría del paper ·
G historia de errores · H glosario · I archivos · J reproducir · K ajustes.

Reglas nuevas: cada sección de método lleva `\enpaper`, `\adaptado` o
`\propio` al comienzo; lo que es historia ("una versión anterior…", "el
22-sep se midió…") va al apéndice G, no al cuerpo. El cuerpo quedó en ~45
páginas; si se quiere más corto, lo siguiente a mover es "Cómo se lee la
figura" del test A y la subsección de interacción de fotones.

Lo que esta reestructuración NO hizo (sigue abierto abajo): B6, B7, B8, C
con datos de laboratorio (actividad, distancia, umbral, carcasa), H
(propagación del error), I (formatos, software, corrida nueva).

---

## 1. Lo que hay hoy (inventario)

| Sección de `informe.tex` | Estado | Comentario |
|---|---|---|
| Cómo leer este informe | completa | mapa de preguntas, convenciones |
| Qué se midió | esqueleto | sesiones y fechas, qué es D. Sin hardware, sin fuentes, sin cadena de señal |
| Del pulso al espectro | completa | registro, onset, línea base, compuerta, control área/amplitud, espectro neto |
| Marco teórico | parcial | E_max, Klein-Nishina, modelo directo, NCC y SA muy bien. Falta toda la física de detección (interacción, centelleo, PMT, resolución) y la estadística |
| La cadena de análisis, paso a paso | completa | texto; falta diagrama y el cierre con verdad conocida |
| Las N calibraciones | completa | tabla, columna por columna, cómo mirar un ajuste |
| Test A (dos fuentes, mismo punto) | completa | |
| Comparación entre detectores (test B) | completa | |
| El detector D con la tierra mal | completa | |
| Comparación con el método publicado | completa | robustez, métricas, auditoría, erratas |
| Cómo usar estos números | completa | receta; falta propagación del error |
| Limitaciones y qué falta medir | completa | |
| Ap. Glosario | parcial | 19 símbolos; faltan los de estadística y los de la simulación |
| Ap. Qué archivo es cada cosa | completa | |
| Ap. Cómo reproducir | completa | |
| Ap. Ajustes de ejemplo | completa | |
| Bibliografía | **no existe** | el paper va en una nota al pie; no hay `.bib` |
| Montaje experimental | **no existe** | |
| Simulación Geant4 | completa | `sec:geant4` (21-sep-2026): qué se simula, qué se guarda, cómo se hace plantilla, qué falta, cómo correrla |

Figuras que existen: `integracion_pulso`, `pulsos_detectores`,
`area_a_espectro`, `teoria_ncc`, `teoria_sa`, `ajuste_A/B/C/D`,
`cociente_fuentes`, `comparacion_por_punto`, `ganancia_dv`,
`comparacion_detectores`. Todas las genera `analisis/informe.py`.

---

## 2. Lo que debe tener el texto guía (lista maestra)

Cada ítem dice qué tiene que decir, de dónde sale (sección de `CLAUDE.md`,
módulo, documento) y si hace falta una figura o tabla nueva (que se agrega en
`analisis/informe.py`, nunca a mano). Lo marcado **[LAB]** necesita un dato
que no está en el proyecto y hay que ir a buscarlo al laboratorio o a la
bitácora.

### A. Portada y orientación

- [x] Autoría y filiación reales (24-sep-2026).
- [x] (24-sep-2026: "Los resultados, en un párrafo" en la Introducción, y la sección Conclusiones) Resumen de una página con los números finales (D por detector y voltaje, b, estado del test A) para quien no va a leer más.
- [ ] Prefacio de medio párrafo: para quién es el texto, qué se asume que sabe, y en qué orden leerlo si se quiere solo la física, solo la receta o solo los resultados. La tabla de "Cómo leer" ya hace la mitad.
- [x] (24-sep-2026: etiquetas DEL PAPER / ADAPTADO / PROPIO) Convención tipográfica explícita: qué es un `aviso` (gris), qué es una decisión de análisis, qué es un resultado. Hoy los avisos están explicados; las decisiones no se distinguen.

### B. Marco teórico

Lo que hay (E_max, Klein-Nishina, modelo directo, NCC, SA) se conserva. Lo
que sigue va **antes** de "Dónde cae el borde", porque un lector nuevo llega
a la fórmula de E_max sin saber por qué el fotón y el electrón se reparten la
energía así.

**B1. Interacción de fotones con la materia** (hecho 21-sep-2026: `sec:interaccion`)
- [x] Los tres procesos con su escala en Z y E; el plástico como bolsa de electrones libres (Klein-Nishina × electrones/gramo reproduce el "incoherente" de XCOM al 0.9 % sobre 100 keV); camino libre medio y probabilidad de interactuar en el cilindro simulado (macros `LibreCamino*`, `ProbInteraccion*`).
- [x] Figura `teoria_interaccion`: BC-404 y NaI, de NIST XCOM (`datos/xcom/`, tablas con procedencia; citado como `\cite{xcom}`).
- [x] Fotopico: por qué lo tiene el NaI y no el plástico, con los cocientes medidos (`FotoSobreCompton*`, `CruceFotoCompton*`).

**B2. Cinemática Compton** (hecho 21-sep-2026: `sec:cinematica`)
- [~] E' y T con las leyes de conservación nombradas y el resultado; la deducción NO está línea por línea (se dice qué se conserva y se da la relación). Si el lector objetivo lo pide, agregar las tres líneas de álgebra.
- [x] Máximo en 180°, E_max, y la errata de la Ec. 1 del paper con el número que da (`EcUnoImpresa*`).
- [x] Pico de retrodispersión, con su energía por línea y el enlace al exceso de bajo canal.
- [x] Tabla `tab:cinematica` generada (E_γ, α, E_max, E_max/E_γ, retrodispersión, E' a 90°); figura `teoria_cinematica` (a).

**B3. Klein-Nishina** (hecho 21-sep-2026: `sec:kn`)
- [x] Supuestos (libre y en reposo), ligadura y Doppler frente a la resolución, con la comparación contra XCOM como verificación (`KNvsXcom`, `KNvsXcomBajo`).
- [x] Límite de Thomson (σ_T = 0.665 b) y la sección eficaz total analítica (`informe.kn_total`, verificada contra la integral numérica de la Ec. 2 a 1e-8).
- [x] El jacobiano escrito (dΩ = 2π d cos θ, dT/d cos θ) y por qué la distribución en energía CRECE hacia el borde.
- [x] Figura `teoria_cinematica` (b) y (c): dσ/dΩ contra Thomson, σ total contra energía. La dσ/dE por línea sigue en `teoria_ncc` (a).
- [ ] Citas de física de detección que faltan: Evans (*The Atomic Nucleus*) o Leo, para la deducción de Compton y la discusión del perfil Doppler; verificar antes de agregar.

**B4. Centelladores orgánicos**
- [ ] Mecanismo: niveles π del anillo bencénico, fluorescencia rápida, por qué el rendimiento es ~10^4 fotones/MeV y el tiempo de decaimiento es de ns (esto justifica la compuerta de 12 ns más adelante). **[LAB]** confirmar que el plástico es BC-404 (la simulación lo asume; `Analisis_BC404_Detector.docx` en `JcV/` lo dice) y tomar sus datos de hoja: rendimiento, τ, longitud de atenuación, densidad.
- [ ] Ley de Birks escrita: dL/dx = S (dE/dx) / (1 + kB dE/dx). Qué es kB, el valor usado en la simulación (0.126 mm/MeV, `ComptonSimulation/src/DetectorConstruction.cc`), y de dónde sale ese valor (citar).
- [ ] Qué es keVee, por qué para electrones el quenching es chico pero no nulo (Lvis/Edep = 0.863 a 50-100 keV, 0.954 sobre 700 keV, §7.17), y por qué eso hace que las dos escalas den D distintas en un 5 %. Hoy está repartido en tres avisos; acá va la explicación de una vez.
- [x] Figura Lvis/Edep contra Edep desde el `.root`: ya está, es el panel (c) de `simulacion_geant4` (hecho 21-sep-2026 con el bloque D). B4 la cita, no la repite.
- [ ] Por qué el test de linealidad se hace en keVee y no en keV (§7.46): en keV dos fuentes tienen que dar D distintas aunque la electrónica sea perfecta.

**B5. Fotomultiplicador y ganancia**
- [ ] Cómo funciona (fotocátodo, dinodos, multiplicación secundaria) y de dónde sale la ley D ∝ V^b: ganancia por dinodo δ ∝ V^k con k ~ 0.7-0.8, n dinodos, b = n·k. Con eso el b = 9.5 medido en el test B (§7.54) se vuelve interpretable (¿cuántos dinodos tiene el tubo?).
- [~] **[LAB]** modelo del PMT: ✔ **XP2262/B, 12 dínodos** (24-sep-2026); con eso k = b/12 ya está en el informe. Falta el divisor de tensión y el HV nominal. (Texto original:) modelo del PMT, número de dinodos, divisor de tensión, HV nominal. `CLAUDE.md` §2 y §7.45 dicen que ni siquiera está confirmado que sea un PMT (se dedujo del voltaje). Esto es bloqueante para B5 y B6.
- [ ] Por qué la ganancia depende tanto del voltaje (50 V son 30 % de ganancia) y qué implica para la reproducibilidad entre sesiones (test A, §7.55).
- [~] Qué es el atenuador y por qué su valor entra derecho en D (§7.45). El VALOR ya está: **20 dB en las 25 corridas** (§7.70), y con eso se puede escribir que las D del informe son "a la entrada del DRS4" y que la carga del PMT es 10x. Falta **[LAB]** el modelo del atenuador y dónde está exactamente en la cadena (antes o después del cable).

**B6. Resolución en energía**
- [ ] El modelo σ² = a·E + b·E² + c y qué proceso físico da cada término (estadística de fotoelectrones, no uniformidad y ganancia, ruido electrónico). Por qué el ajuste usa σ_rel como fracción del borde (un solo término, el estadístico) y qué se pierde con eso (§7.35: los términos extra se van a cero porque no hay palanca).
- [ ] Orden de magnitud esperado para plástico + PMT a 500-1000 keV, para que el lector sepa si 20-30 % FWHM es razonable o es un síntoma (§7.20 tiene la historia de cuando 43 % parecía real y era instrumental).
- [ ] Cómo el ancho de la compuerta entra en la resolución (§7.53, sección "Cuatro decisiones").

**B7. Qué deforma el espectro real, además de la resolución**
- [ ] Scattering múltiple: por qué llena la zona sobre el borde, cómo depende del tamaño del detector, cuánto es acá (25 % de los eventos, §7.17), y por qué no se puede absorber en σ (no es gaussiano). Hoy está en un ítem de la lista de "Por qué el borde no es un escalón"; merece figura.
- [~] Figura de la plantilla separada por número de scatterings: es el panel (b) de `simulacion_geant4` (21-sep-2026). Falta superponerle Klein-Nishina para cerrar el argumento; puede ir en la figura de B7 o agregarse a ese panel.
- [ ] Escape de electrones por las paredes y efecto del tamaño finito.
- [ ] Retrodispersión en el entorno (portafuente, mesa): el exceso de bajo canal (§7.32, §7.35). Qué se sabe, qué se descartó (no es umbral, no es fondo), qué falta medir.
- [ ] Apilamiento y cósmicos: por qué hay eventos a áreas enormes, cómo se ve en el fondo (§7.37), y por qué `RANGO_MAX` fijo (§7.7).
- [ ] La curva de eficiencia del umbral de trigger: por qué el espectro sube desde cero, y por qué eso es un piso duro para cualquier ventana (§7.48, `espectro.fin_del_umbral`).

**B8. Digitalización con el DRS4**
- [ ] Qué es un switched-capacitor array, por qué el muestreo no es uniforme (celda a celda) y por qué la integral tiene que usar el Δt real (la cadena ya lo dice; acá va el porqué).
- [ ] Rango dinámico (±0.5 V), saturación, número de celdas, tasa de muestreo, y qué es el tiempo muerto del DRSOSC (la placa dispara a 170 Hz y se guardan 71 ev/s, §7.34). Esto explica el tiempo vivo y el scaler.
- [ ] Trigger: umbral, retardo, por qué la posición del pulso en el registro cambia entre sesiones (48 ns contra 120 ns, §2) y por qué eso rompió el análisis viejo (§7.33). Enlazar con `docs/formato_binario_drs4.md`.
- [ ] Referencia al apéndice de formatos de archivo (XML del DRSOSC y binario).

**B9. Estadística del ajuste** (hecho 20-sep-2026: `\subsection{Estadística del ajuste}`, `sec:estadistica`)
- [x] Poisson en cada bin, la propagación con fondo restado y el piso en 1.
- [x] χ² de Neyman contra χ² de Pearson, con el sesgo MEDIDO con verdad conocida (tabla `tab:cierre`; `modelo.Problema(varianza="modelo")` es la opción Pearson, sólo para medirlo).
- [x] La amplitud se despeja (VARPRO, con la fórmula y la cita).
- [ ] La versión con nnls cuando hay dos columnas (A y C del exceso, §7.47): sólo mencionada; va con el modelo extendido (bloque E).
- [x] Curvatura contra perfil, con la figura (`teoria_estadistica` a) y los dos números de la corrida de ejemplo.
- [ ] El caso de c0 (21σ contra 2.1σ, §7.49) como ejemplo de cuando no coinciden: falta, va con el test del origen.
- [x] La correlación D-σ y la elipse; el escalado por sqrt(χ²/ndf) dicho como parche.
- [x] χ²/ndf y lo que no dice.
- [x] Sistemáticos: común contra punto a punto (cualitativo).
- [ ] Tabla del presupuesto de error de la corrida de referencia (`sistematicos.py`): falta generarla en `informe.py`.
- [x] Cierre con verdad conocida, por la cadena COMPLETA (`informe.cierre_conocido`, 30 realizaciones, Neyman y Pearson), con tabla y figura. Hallazgo: la barra reportada queda ~15 % corta porque la ventana se elige de nuevo en cada realización (ver `CLAUDE.md` §7.67).

**B10. NCC** (hecho 20-sep-2026)
- [x] Savitzky-Golay explicado (polinomio local de grado 3) y citado; el efecto del número de puntos MEDIDO sobre la corrida de ejemplo (`ncc.buscar_intervalo(n_pts=)`, macros `Npts*`).
- [x] Figura `teoria_ncc_detalle`: las dos curvas estandarizadas, una contra otra con ρ, y ρ/B contra el número de puntos. El ρ recalculado coincide con el del barrido (`\EjRhoRecalc` contra `\EjRho`).

**B11. Recocido simulado** (hecho 20-sep-2026)
- [x] Citas: Kirkpatrick 1983, Metropolis 1953, Geman y Geman 1984 (el teorema del enfriamiento logarítmico), Golub y Pereyra 1973.
- [x] Párrafo "Lo que se verificó del recocido, y dónde": `verificar_sa.py` y los dos PDF de `docs/`.

**B12. Los métodos clásicos y su sesgo** (hecho 20-sep-2026: `sec:clasicos`, antes del NCC)
- [x] Los tres métodos descritos; tabla `tab:clasicos` generada (sintético de Klein-Nishina del Cs-137, con y sin ruido, `informe.sesgo_clasicos`) y los tres sobre las 25 corridas reales contra la D del ajuste (`informe.clasicos_reales`, macros `ClasReal*`).
- [x] Figura `teoria_clasicos` (importa `codigo/etapas/clasicos.py`, no lo copia).
- [x] Nota histórica con Knox y Miller 1972 (el 89 % de la altura) y Dietze y Klein 1982 (Monte Carlo convolucionado), citas verificadas.
- [ ] Pendiente de física: explicar por qué sobre datos reales el sesgo (~10-12 %) es mayor que en el sintético (~4-6 %); el texto lo atribuye al scattering múltiple, que es plausible pero no está medido aparte. Se puede medir corriendo los clásicos sobre la plantilla de Geant4 convolucionada, sin ruido.

### C. Montaje experimental (sección nueva, va después de "Qué se midió")

- [ ] Diagrama de la cadena de señal: centellador, PMT, HV, atenuador, cable, DRS4, PC. Con los valores por sesión. **[LAB]**
- [ ] Los cuatro detectores: material, dimensiones, envoltura, acople óptico, PMT. Foto o esquema. **[LAB]**
- [ ] Las fuentes: nuclido, actividad nominal y fecha de referencia, encapsulado, distancia y posición respecto del detector, portafuente. **[LAB]** Importa para el exceso de bajo canal y para la simulación.
- [ ] Tabla de puntos de operación por corrida: HV, atenuador, umbral, retardo, tasa de muestreo, canal del DRS4, formato, fondo propio sí/no. `analisis/rutas.Corrida` ya tiene los campos y desde el 22-sep-2026 el atenuador está en las 25 (20 dB, §7.70); generarla desde ahí y marcar en la tabla lo que siga sin registrar —hoy: umbral en todas menos `co60`, y la distancia en todas— (§7.45).
- [ ] Qué se ve en la bitácora (`JcV/Bitacora`) de cada sesión, resumido: qué se cambió entre sesiones (la tierra del D, el retardo del trigger, la HV).
- [ ] Adquisición: cómo se configura el DRSOSC, dónde está el `drsosc.cfg`, qué es "Save as binary", y el `drs_daq.ps1` si se usó. Enlazar con el apéndice de formatos.

### D. La simulación de Geant4 (hecho 21-sep-2026: `sec:geant4`, entre el marco teórico y la cadena)

- [x] Geometría simulada, leída de la fuente: cilindro de 2 x 2 pulgadas de PVT (`G4_PLASTIC_SC_VINYLTOLUENE`, la base de BC-404/408), **con carcasa de Al (pared 10 mm, tapa 6 mm) y stub de PMT de vidrio**, en aire, fuente puntual a 10 cm con haz `pencil`. El "en vacío" y "sin carcasa" que decía este plan venían de la simulación equivocada (ver abajo). Dice explícitamente que NO está confirmado que coincida con el detector real, y cita el control de `w` sólo cualitativamente (el número, 0.796 ± 0.070 contra 0.751, es de `codigo/`; pasarlo a macro cuando el bloque E porte el modelo extendido).
- [x] Lista de física (`G4EmStandardPhysics_option4` + `G4DecayPhysics`), corte de 0.1 mm, Birks 0.126 mm/MeV con la ley escrita, citada (`birks1951`) y la trampa de que Geant4 no la trae para el PVT; aplicado paso a paso y por qué.
- [x] Qué se guarda por evento (las 8 columnas) y cuántos eventos: tabla `tab:simulacion` generada de los `.root` (primarios del nombre, filas, fracción que deposita, reparto por `NCompton`, fotoeléctricos, absorción total, `<L/Edep>` en tres zonas, bordes en luz).
- [x] Cómo se convierte en plantilla: 350 bines (`\SimBinsPlantilla`, leído de la firma de `simulacion.cargar`), interpolación lineal en ch/D (por eso no hay peine de aliasing en `analisis/`; el de §7.25 era del rebinning de `codigo/`), kernel normalizado, y **el binning medido** reajustando la corrida de ejemplo con 100/175/350/700 bines: D se mueve 0.06 % (`\BinsDesvMax`). Hallazgo lateral: la barra por curvatura sí depende del binning (0.0017 a 0.0031) porque la plantilla es lineal a trozos; queda dicho como otra razón para el perfil de Δχ².
- [x] La verificación de que la plantilla es de la fuente correcta (`lineas_primarias`) y por qué existe.
- [x] Figura `simulacion_geant4`: (a) Edep y Lvis del Co-60 con los bordes en las dos escalas, (b) Cs-137 por número de scatterings, (c) quenching evento a evento, (d) razón dato/modelo de la corrida de ejemplo en todo el rango (exceso y umbral). El panel de aliasing no se hizo porque en `analisis/` no ocurre (ver arriba).
- [x] Qué NO está simulado: entorno (→ exceso), umbral, PMT/electrónica, resolución, fondo/cósmicos/apilamiento, y si la geometría es la real.
- [x] Cómo correrla de nuevo: `simulacion/` (copia de la fuente de WSL, con `LEEME.md`), macros por fuente, nombre del archivo con el `beamOn`, qué tocar para otra fuente u otro detector.
- Lo que apareció al escribirlo (ya está en `CLAUDE.md` §7.69): **`informe.geant4_geometria` leía `JcV/ComptonSimulation/`, que es OTRA simulación** (CSV, sin carcasa, iones decayendo) y no la que produjo los `.root`; los cuatro números coincidían por casualidad de diseño. La real estaba sólo en WSL. Y el 2 % de eventos con `NCompton == 0` y sin fotoeléctrico son **electrones Compton nacidos en la tapa de Al** que entran al plástico (mediana 275 keV): la carcasa no es decorativa.
- [ ] Pendiente de física, para B7: cuánto de la forma sobre el borde cambia con el tamaño del detector (correr la simulación con `/det/radius` y `/det/halfZ` distintos). Sin eso "el scattering múltiple depende del tamaño" queda afirmado y no medido.

### E. La cadena de análisis (existe; agregar)

- [x] (24-sep-2026: `flujo.tex`, lado a lado con el del paper) Diagrama de flujo de una página: archivo, pulsos, áreas (+ caché), espectro neto (+ tiempo vivo), plantilla, NCC (ventana + semilla + cotas), SA multiarranque, pulido local, covarianza, D. Con el nombre del módulo de `analisis/` en cada caja.
- [ ] Una figura por corrida de referencia con las cuatro etapas (las pestañas 1 a 4 de `analisis/gui.py` ya las dibujan, §7.58): datos y plantilla, NCC, recocido, inicial contra final.
- [ ] Sección "Cierre con verdad conocida" (ver B9). Resultados de `registro.py --cierre` y de `codigo/ajuste/mcmc.py`.
- [x] (24-sep-2026: apéndice `ap:historia`) Sección "Errores que se cometieron y cómo se detectaron". Para un texto guía vale más que cualquier otra: (1) la línea base que se comía el pulso y el factor 3.4 en D (§7.33); (2) el tiempo de pared contra el tiempo vivo (§7.34); (3) el canal cableado que dio cero pulsos en silencio (§7.36); (4) la plantilla de la fuente equivocada que da un número bien formateado (§7.41); (5) las constantes de canal de la escala vieja escondidas en la GUI (§7.40); (6) los Co60.dat duplicados del detector B (ya está en "Limitaciones"). Una moraleja común: mirar pulsos crudos antes de refinar el ajuste (§9-bis.5).
- [ ] Modelo extendido: qué es (w, exceso exponencial), qué mide, por qué el informe no lo usa (§7.35, §7.47, "Cuatro decisiones"). Hoy se menciona como "una versión anterior"; contar que sigue disponible y qué da.
- [ ] Las dos implementaciones (`codigo/` y `analisis/`) y que coinciden al 0.2 % (§7.47). Es un control de software que el lector debería conocer.

### F. Resultados (existe)

- [ ] Nada estructural. Revisar que cada tabla tenga en el epígrafe cómo leerla (las de test A y test B lo tienen).
- [ ] Agregar la figura resumen del test A que hoy solo se cita por nombre (`resumen_linealidad.png`) como figura de página, o decir por qué no.

### G. Comparación con el paper (existe)

- [x] (24-sep-2026: apéndice `ap:auditoria`, 5 erratas + 2 ambigüedades + constantes) Las tres erratas están (Ec. 1, Ec. 12, SD/R² intercambiados) más el coseno. Agregar la ambigüedad Ec. 14 / Metropolis (§3 de `CLAUDE.md`, §7.22: con el costo NCC las amplitudes son indeterminables) porque afecta a quien quiera implementar el paper.
- [ ] Constantes del paper que no son transferibles (T_init = 2000, |ΔE| < 0.001, §7.22). Está en "La temperatura inicial no se puede copiar"; enlazar.
- [ ] Referencia a `replica_paper/` y a su GUI como la forma de ver el paper paso a paso.

### H. Uso y limitaciones (existe; agregar)

- [x] (24-sep-2026: `sec:error-energia` y apéndice `ap:propagacion`) Propagación del error a la energía convertida: σ_E/E = σ_D/D, más el sistemático, y qué hacer con la correlación D-σ si se quiere la resolución.
- [ ] Cómo interpolar a otro voltaje con la ley de potencia y con qué error.
- [~] (24-sep-2026: ejemplo de 1000 mV·ns con su error; falta el paso desde un pulso crudo) Ejemplo numérico completo, de un pulso hasta keV con su error.

### I. Apéndices

- [ ] **Glosario ampliado**: α, w, λ, C, T (temperatura del SA), k_m, k_f, tiempo vivo, scaler, onset, walk y jitter del trigger, VARPRO, χ² de Pearson y de Neyman, R̂ y n_eff (si entra el MCMC), keVee, quenching, kB.
- [ ] **Formatos de datos**: resumen de `docs/formato_binario_drs4.md` y de la estructura del XML (§2), con los tamaños y tiempos de parseo.
- [ ] **Software**: qué hace cada módulo de `analisis/`, las tres interfaces (`gui.py`, `visor.py`, `replica_paper/gui.py`) y para qué sirve cada una, los verificadores de `verificaciones/`. Una tabla, no prosa.
- [ ] **Cómo agregar una corrida nueva**: dónde copiarla (`datos/Detector X/<voltaje>/<fuente>`), cómo darla de alta en `rutas.py`, qué revisar antes de calibrar (`revisar.py`), cómo regenerar el informe. Es el procedimiento que un sucesor va a necesitar primero.
- [ ] **Historia**: por qué hay `codigo/` y `analisis/`, y qué cambió de la Etapa 8 a la 9 (el factor 3.4). Corto.
- [~] **Bibliografía**: `referencias.bib` existe (20-sep-2026) con 10 entradas verificadas (Kong, Metropolis, Kirkpatrick, Geman, Savitzky-Golay, Golub-Pereyra, Baker-Cousins, Dietze-Klein, Knox-Miller, Cowan) y `\bibliographystyle{unsrt}` antes de los apéndices; `informe.py` corre bibtex. Faltan las de física de detección, que entran con los bloques B1-B8. Candidatas, todas a verificar antes de citar:
  - Kong et al., IEEE TNS 2025, arXiv:2507.16282 (el método).
  - Klein y Nishina 1929, Z. Phys. (la sección eficaz).
  - Knoll, *Radiation Detection and Measurement*; Leo, *Techniques for Nuclear and Particle Physics Experiments* (interacción, centelladores, PMT, resolución).
  - Birks 1964, *The Theory and Practice of Scintillation Counting* (quenching).
  - Dietze y Klein 1982, NIM; Knox y Miller 1972, NIM (borde Compton en centelladores orgánicos, métodos clásicos).
  - Agostinelli et al. 2003, NIM A (Geant4); Allison et al. 2016 (Geant4 10).
  - Ritt, DRS4 (el chip); manual de la DRS4 Evaluation Board.
  - Metropolis et al. 1953; Kirkpatrick, Gelatt y Vecchi 1983 (SA).
  - Savitzky y Golay 1964, Anal. Chem.
  - Golub y Pereyra 1973 (VARPRO).
  - Cowan, *Statistical Data Analysis*, o Bevington (χ², Pearson/Neyman, covarianza).
  - NIST XCOM (secciones eficaces).
  - Hoja de datos del BC-404 y del PMT **[LAB]**.

---

## 3. Orden propuesto de sesiones

Cada fila es una sesión de trabajo. Las que dependen de datos del laboratorio
van después, para que lo que se puede escribir hoy no espere.

| # | Bloque | Entregable en `informe.tex` | Figuras/tablas nuevas en `informe.py` | Requiere | Hecho |
|---|---|---|---|---|---|
| 1 | B1 + B2 + bibliografía | subsecciones "Interacción de fotones" y "Cinemática Compton"; `referencias.bib` y `\bibliography` | secciones eficaces en BC-404; tabla E_γ / E_max / backscatter | nada | hecho 21-sep-2026 |
| 2 | B3 + B12 | jacobiano escrito; "Los métodos clásicos y su sesgo" | dσ/dΩ y dσ/dE; sesgo de los clásicos (de `clasicos.py`) | nada | hecho 20/21-sep-2026 |
| 3 | B9 + B10 + B11 | "Estadística del ajuste" completa; cierre con verdad conocida; NCC por dentro; citas del SA | tabla del cierre (Neyman/Pearson); `teoria_estadistica`, `teoria_ncc_detalle` | nada | hecho 20-sep-2026 |
| 4 | D | sección "La simulación de Geant4" | plantilla Edep/Lvis, descomposición por NCompton, quenching, dato/modelo; tabla de los `.root`; binning medido | leer `simulacion/src/` (la fuente REAL, de WSL), `.root` | hecho 21-sep-2026 |
| 5 | B4 + B7 | "Centelladores orgánicos y Birks"; "Qué deforma el espectro" | Klein-Nishina sobre el panel (b) de `simulacion_geant4`; escape y tamaño (barrido de `/det/radius`) | hoja de datos del BC-404 **[LAB]**; Lvis/Edep ya está en `simulacion_geant4`(c) | [ ] |
| 6 | E (diagrama, trampas, modelo ext.) | flujo de una página; "Errores que se cometieron" | diagrama; cuatro etapas de una corrida | nada | [ ] |
| 7 | C | sección "Montaje experimental" | cadena de señal; tabla de puntos de operación | PMT, fuentes, distancias **[LAB]**; atenuador ✔ 20 dB (§7.70) | [ ] |
| 8 | B5 + B6 + B8 | "Fotomultiplicador", "Resolución", "Digitalización" | ninguna o una (pulso con celdas no uniformes) | modelo del PMT **[LAB]** | [ ] |
| 9 | H + I (glosario, software, corrida nueva, historia) | propagación del error; apéndices | tabla de módulos | nada | [ ] |
| 10 | A + G + F + revisión | portada, resumen ejecutivo, consistencia de referencias cruzadas, índice de figuras | ninguna | todo lo anterior | [ ] |

Regla de cierre de cada sesión: `python analisis/informe.py` compila sin
errores ni referencias sin resolver, y el PDF se abre y se mira la sección
nueva a tamaño real (por §7.64 de `CLAUDE.md`, lo que se ve chico en pantalla
no se arregla subiendo la letra).

---

## 4. Reglas para no romper lo que ya funciona

Vienen de `informe/README.md` y de `CLAUDE.md` §7.64-7.66; se repiten acá
para que este archivo alcance.

- **Ningún número se tipea.** Si un párrafo nuevo necesita un número que
  sale de los datos, se agrega como macro en `analisis/informe.py`
  (`numeros`) y en el `.tex` va `\Macro`. Las constantes físicas (m_e c²,
  E_γ) ya son macro.
- **Ninguna figura se pega.** Se genera en `informe.py` con `figura(ancho,
  alto)`, letra del documento, PDF, sin `bbox_inches="tight"`, y el
  `width=` del `\includegraphics` es el mismo `ancho`.
- **Las figuras teóricas también salen de código** (`teoria_ncc`,
  `teoria_sa` son el modelo): así un cambio en `fisica.py` se ve en el
  informe.
- **Las tablas salen de `generado/`** y se traen con `\input`. El epígrafe
  vive en el `.tex`.
- **Avisos en gris** solo para lo que se puede malinterpretar. No para
  decorar.
- **Lo que es decisión de análisis** (compuerta de 12 ns, `inicio_rel`,
  sin exceso, sin fondo prestado) se dice como decisión, con su costo.
- **Las citas van al `.bib`**, no a notas al pie, y ninguna se agrega sin
  verificar la referencia (regla general del proyecto: no adivinar).
- Antes de escribir un tema de física, revisar si `CLAUDE.md` §7 ya tiene
  el número medido: el texto guía se ancla en lo medido acá, y la teoría
  general va para explicar ese número, no al revés.

---

## 5. Datos que hay que conseguir fuera del proyecto **[LAB]**

Bloquean los bloques C, B5, B6 y parte de B4. Conviene juntarlos en una sola
visita y anotarlos en la bitácora.

- [x] ✔ XP2262/B, 12 dínodos, confirmado por el usuario el 24-sep-2026. (Era: Modelo del fotodetector de cada detector (¿PMT? ¿cuál? ¿cuántos dinodos?) y HV nominal de operación. Sin esto la ley de ganancia no se puede interpretar (§7.45).
- [x] ✔ BC-408, cilindro de 2"×2", confirmado el 24-sep-2026; la carcasa sigue sin medir. (Era: Modelo y dimensiones del plástico de cada detector; si es BC-404 (lo que dice el informe y `Analisis_BC404_Detector.docx`) o BC-408 (lo que dice el README de la simulación; misma base PVT, misma densidad y composición, así que para el depósito da igual, pero hay que decir cuál es) y si mide 2 x 2 pulgadas como asume la simulación. También la carcasa: la simulación pone pared de Al de 10 mm y tapa de 6 mm ("medida", dice el comentario del `.cc`, pero no consta de qué detector), y un vidrio de 3 cm detrás. **Ojo:** desde el 21-sep el marco teórico cita ese tamaño (`\LadoDetector`, leído del `.cc` de Geant4) para la probabilidad de interacción; si el detector real mide otra cosa, hay que decirlo ahí.
- [x] Atenuador en cada sesión: **20 dB en todas**, para que el pulso del PMT no pase el rango de entrada del DRS4 (22-sep-2026, §7.70). Falta sólo dónde está exactamente en la cadena y de qué modelo es.
- [ ] Umbral y retardo del trigger por sesión (algunos están en el nombre del archivo o en `drsosc.cfg`).
- [ ] Fuentes: actividad y fecha de referencia, encapsulado, distancia al detector, portafuente y qué hay alrededor (para el exceso de bajo canal).
- [ ] Qué se cambió en el detector D entre el 1-sep y el 8-sep (la tierra), con palabras del que lo hizo.
- [ ] Foto o esquema del montaje.
