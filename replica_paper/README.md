# `replica_paper/` — réplica literal de Kong et al. (IEEE TNS 2025)

Réplica del método y del **análisis de resultados** de

> Kong, Zhang, Lin, Qiu, Chen, Wang. *The Automatic Calibration Method of the
> Compton Edge Based on Normalized Cross-correlation and Simulated Annealing
> Algorithm.* IEEE TNS 2025, arXiv:2507.16282.

**Regla de esta carpeta: acá sólo entra lo que el paper dice.** Plantilla
Klein-Nishina analítica (Ec. 2 y 10), no Monte Carlo. Nada de modelo extendido,
mezcla de scattering, término de exceso a bajo canal, χ² de Pearson iterado,
MCMC ni presupuesto de sistemáticos medido: todo eso vive en `codigo/` y
`analisis/` y es **extra** respecto del paper. Lo que en esta carpeta no está en
el paper aparece marcado `[NO está en el paper]` en el propio código.

Lo único prestado de afuera es el paso XML → `(canal, cuentas)`, que es
exactamente donde **arranca** el diagrama de flujo del paper (Fig. 2: *"Energy
Spectrum (Channel, Counts)"*). Ese préstamo está aislado en `datos.py`.

---

## Cómo correr

```bash
python replica_paper/gui.py              # ⭐ EL PAPER PASO A PASO, calculando en vivo
python replica_paper/ecuaciones.py       # redibuja la Fig. 1 del paper (validación de la Ec. 10)
python replica_paper/ncc.py              # Sección II-B sobre sintético con verdad conocida
python replica_paper/sa.py               # Sección II-A con las constantes del paper
python replica_paper/convencionales.py   # el sesgo de los dos métodos de comparación
python replica_paper/metodo.py           # la cadena de la Fig. 2, uno y dos bordes
python replica_paper/datos.py            # qué espectros hay
python replica_paper/auditoria.py        # AUDITORÍA de los números publicados
python replica_paper/correr.py --cierre  # cierre con verdad conocida del análisis de resultados
python replica_paper/correr.py --rapido  # nuestras corridas, sin los 10 bloques
python replica_paper/correr.py           # todo
python replica_paper/verificar_gui.py    # que la interfaz corra y no invente numeros
```

## Qué archivo cubre qué del paper

| archivo | del paper |
|---|---|
| `ecuaciones.py` | Ec. 1, 2, 3 + Tabla I, 8, 9, 10, 11, 12, 13, 14, 15, tal como se imprimen |
| `ncc.py` | Sección II-B: selección automática del intervalo por Pearson |
| `sa.py` | Sección II-A: recocido simulado, con `T_init = 2000`, `α = 0.95`, 100 iteraciones por temperatura, `Δx ~ U(−0.1δ, 0.1δ)`, corte por `T < 1e-6` o `\|ΔE\| < 0.001` en 10 ciclos |
| `convencionales.py` | los dos métodos de comparación (Gaussian FWHM, Difference) |
| `metodo.py` | el diagrama de flujo de la Fig. 2, de punta a punta |
| `resultados.py` | Sección II-C y III: recta con SD y R², similitud coseno, error relativo, reducción del error medio, robustez de intervalo |
| `auditoria.py` | reproduce Fig. 3, Fig. 5 y las erratas desde los propios números publicados |
| `datos.py` | **fuera del paper**: nuestras corridas → `(canal, cuentas)`, y los 10 bloques |
| `correr.py` | el informe completo |
| `verificar_gui.py` | **fuera del paper**: corre las 14 pestañas y compara sus números contra `informe.txt` |
| `gui.py` | **fuera del paper**: la interfaz. Una pestaña por paso, con el razonamiento al lado y los controles para mover lo que el paper deja sin especificar |

---

## `gui.py`: el paper paso a paso

`python replica_paper/gui.py` abre catorce pestañas, en el orden del propio
paper. Cada una tiene los controles a la izquierda, el gráfico al centro y, a
la derecha, cuatro secciones fijas: **qué dice el paper** (la cita literal),
**qué corre esta pestaña** (el módulo y la función), **lo que el paper no dice**
(las ambigüedades, que son exactamente los controles de la izquierda) y **qué
mirar**. Abajo a la izquierda, el resultado numérico recién calculado.

| pestaña | qué muestra |
|---|---|
| **· Mapa** | el diagrama de la Fig. 2 con la pestaña que calcula cada bloque, las ecuaciones con su función y su pestaña, las 5 erratas con su evidencia y las 9 ambigüedades con lo que se decidió |
| **0 · Datos** | donde arranca la Fig. 2: muestra, fondo y neto, con `k_m`, `k_f` y los bordes teóricos |
| **1 · Ec. 1** | la fórmula impresa contra la corregida, y los rótulos de la Fig. 3 encima |
| **2 · Ec. 2 y 10** | la Fig. 1 del paper recalculada, y el chequeo de que la Ec. 10 es la Ec. 2 cambiada de variable (las dos integran a la misma σ total) |
| **3 · Ec. 8-9** | la lineshape llevada a canal y convolucionada; el efecto de `pad`; el máximo y el 50% de altura, que **no** caen en `B` |
| **4 · Ec. 15** | el barrido de intervalos, S1 contra S2, el mapa de ρ y el barrido de σ de la plantilla |
| **5 · II-A** | el recocido con las constantes del paper, su historia de temperatura, y los 4 casos que muestran que `T_init = 2000` no es transferible |
| **6 · Fig. 2** | la cadena completa: NCC → SA → pulido, los tres métodos en la misma ventana, y el diagnóstico del solapamiento |
| **7 · Convencionales** | los dos métodos de comparación, y su sesgo medido contra un borde de posición conocida |
| **8 · Fig. 5** | los seis intervalos del paper escalados a nuestro borde, contra sus propios seis paneles |
| **9 · II-C** | las "diez mediciones", con la advertencia de qué prueban y qué no |
| **10 · Auditoría** | los números publicados contra sí mismos: SD/R² cruzados, rótulos cruzados, Ec. 1, Ec. 12, coseno, Fig. 5 |
| **11 · Sec. III** | la recta canal-energía, y por qué no se puede cerrar con lo medido |
| **12 · Cierre** | dos fuentes al mismo punto de operación con `D` conocida: la única parte donde se puede medir si el método acierta |

Dos reglas de la interfaz:

* **no calcula nada por su cuenta**: todo sale de llamar a los módulos de esta
  carpeta, así que no puede divergir de la réplica;
* **no cachea resultados**: cada `Recalcular` rehace el cálculo desde el
  espectro. Lo único que se reusa es el `.areas_gate.npz` de `analisis/area.py`,
  que es el *parseo* del archivo (clave = tamaño + fecha), no un resultado.

`python replica_paper/verificar_gui.py` chequea las dos: corre cada pestaña de
punta a punta sin `mainloop` (`leer_controles` → `trabajo` → `pintar`) y después
compara valor por valor contra `informe.txt`. Los 20 números coinciden.

---

## Lo que salió de auditar el paper

`python replica_paper/auditoria.py`

### 1. La Ec. 10 es correcta y la Fig. 1 se reproduce exacta

Con `A = 1`, `r₀ = 2.818e-13 cm` y `m₀c² = 511 keV`, la Ec. 10 da los valores del
eje de la Fig. 1 (en cm²/keV):

| E₀ [keV] | dσ/dE en E→0 | pico en el borde | el paper |
|---|---|---|---|
| 500 | 1.0199e-27 | 1.6795e-27 | 1.01 / 1.68 |
| 662 | 5.8179e-28 | 1.1253e-27 | 0.58 / 1.12 |
| 1000 | 2.5497e-28 | 6.5165e-28 | 0.255 / 0.645 |
| 2000 | 6.3742e-29 | 2.8483e-28 | 0.064 / 0.283 |
| 5000 | 1.0199e-29 | 1.0511e-28 | 0.010 / 0.11 |

### 2. La Ec. 1 impresa está mal, y el propio paper usa la corregida

`E_max = E_γ/(1 + 2E_γ/m_ec²)` da la energía del **fotón retrodispersado**. Para
Cs-137 daría 184.3 keV; la Fig. 3 del paper rotula **477**, que es lo que da
`E_max = 2αE_γ/(1+2α)`. Igual para 639, 963, 1061 y 1117 keV.

### 3. La Ec. 12 impresa está mal, y se puede medir en sus propios datos

Ec. 12: `CE_ch,2 = B₁ + D·ΔE` con *"ΔE = E₂ − E₁ energy difference between gamma
lines"*. En el panel de LaBr₃ los dos bordes del Co-60 caen en los canales 469 y
543, o sea **74 canales** de separación. Con `D` de la propia recta ajustada:

* ΔE de las líneas gamma (159.3 keV) → 76.0 canales
* ΔE de los **bordes** (154.7 keV) → 73.8 canales ✔

### 4. SD y R² están intercambiados en el recuadro de la Fig. 3, en los dos paneles

Los puntos rotulados de la Fig. 3 reproducen las rectas impresas a 4-5 cifras, lo
que además fija que **SD = RMS de los residuos** (ddof = 0):

| ajuste | ajuste de sus puntos | recuadro del paper |
|---|---|---|
| LaBr₃ fotopico (azul) | y = 2.0692x − 4.756, RMS 1.609198, R² 0.999974 | SD 2.183036, R² 0.999924 |
| LaBr₃ borde (rojo) | y = 2.0959x − 21.731, RMS 2.183036, R² 0.999924 | SD 1.609198, R² 0.999974 |
| NaI fotopico (azul) | y = 0.8525x − 58.882, RMS 1.981123, R² 0.999960 | SD 7.274524, R² 0.999055 |
| NaI borde (rojo) | y = 0.8397x − 40.183, RMS 7.274524, R² 0.999055 | SD 1.981123, R² 0.999960 |

Con la asignación correcta, la frase de la Sección III —*"the standard deviation
associated with the NCC-SA calibration consistently remains below that of the
full energy peak fitting method"*— **dice lo contrario de lo que muestran sus
propios datos**: 2.183 contra 1.609 en LaBr₃ y 7.275 contra 1.981 en NaI.

### 5. Siete rótulos cruzados en el panel de NaI

`Na22 (846, 662)` (662 keV es del Cs-137), `Cs137 (1050, 834)` (es del Mn-54),
`Mn54 (1442, 1173)` (es del Co-60), `Co60 (620, 477)` (es el borde del Cs-137),
etc. No cambian ningún ajuste —los pares canal-energía están bien— pero
atribuyen líneas a la fuente equivocada. El panel de LaBr₃ está bien.

### 6. La similitud coseno no es verificable

El paper reporta 99.9993% (LaBr₃) y 99.9999% (NaI) y **nunca dice entre qué
vectores**. Probadas cinco lecturas razonables (coeficientes (m, b); energías
reconstruidas en los canales medidos; en toda la escala; reconstruida contra
teórica para cada serie), ninguna reproduce los dos valores.

Y una advertencia que el cierre con verdad conocida vuelve concreta: entre dos
calibraciones cuyas **pendientes difieren un 3.0%**, la similitud coseno de las
energías reconstruidas da **0.99998812**. La métrica es insensible por
construcción; hay que reportar residuos en keV y la ordenada con su error.

### 7. Robustez de la Fig. 5, con número

Los seis paneles dan bordes 1521 / 1538 / 1540 / 1537 / 1540 / 1535 → dispersión
total **1.24%**; **0.33%** sacando el intervalo más angosto (`[1500,1600]`, el
único que no contiene el flanco entero), que es el que se aparta.

### 8. Cuatro inconsistencias internas más

1. La Ec. 14 define la loss como mínimos cuadrados; el párrafo de Metropolis dice
   que `E(x)` es *"the NCC-based matching error"*. Hay que elegir uno.
2. La Sección III dice que la robustez se probó *"calibrating EJ309 spectra"*; el
   epígrafe de la Fig. 5 dice BC408.
3. El borde del Co-60 de 1332.5 keV se rotula 1117 keV; la Ec. 1 corregida da
   1118.1.
4. `T_init = 2000` y `|ΔE_mean| < 0.001` son constantes **con unidades del
   costo**, y el costo no está fijado (punto 1). Medido en `sa.py`: con la Ec. 14
   sobre un espectro sintético, `T = 2000` acepta el **1.2%** de las propuestas;
   con el costo NCC (rango [0, 2]) acepta el **72%**. La constante no es
   transferible.

---

## Lo que el paper no especifica y hubo que decidir

Cada uno es un parámetro explícito en el código, con el default y el motivo.

| ambigüedad | dónde | qué se hizo |
|---|---|---|
| con qué σ se construye S2 (la plantilla del NCC) | `ncc.plantilla` | se barre una grilla y gana el de mayor ρ, que es el único criterio que el paper da |
| dónde empieza *"the Compton plateau"* y termina *"the edge minimum"* | `ncc.INI_REL`, `FIN_REL` | `[0.45·E_max, (1+4.5σ_rel)·E_max]`. **No es cosmético:** fija el mapa, porque `B = lo + (1−ini)/(fin−ini)·(hi−lo)` |
| largo mínimo de un intervalo candidato | `ncc.MIN_BINS` | 25 bins. Sin mínimo, Pearson encuentra correlación alta en un puñado de bins de cola vacía |
| cotas del SA (`δᵢ`, *"the parameter search range"*) | `metodo.COTAS_REL` | derivadas de la semilla del NCC |
| de dónde sale `x₀` | `metodo.semilla_por_barrido` | barrido grueso de (B₁, D, σ₁) con A₁, A₂ exactas por nnls. **Hace falta**: medido sobre sintético con verdad conocida, con 24 arranques dispersos al azar ni el SA ni el pulido local llegan al óptimo global |
| el lazo *"Optimal Parameter Solution? No → modify"* | `metodo._cerrar` | pulido local (`least_squares`) desde cada salida del SA, se queda con la de menor loss |
| con qué se suaviza para el método de la diferencia | `convencionales.VENTANA` | Savitzky-Golay de 21 bins (es lo que la Fig. 5 muestra como *"Smoothed Data"*) |
| qué región usan los métodos convencionales | `correr.tres_metodos` | la misma que elige el NCC. Es lo único que hace la comparación justa |
| qué es *"the maximum amplitude"* del método gaussiano | `convencionales.metodo_gauss` | la amplitud de la gaussiana ajustada, o sea borde = `μ + σ√(2ln2)` |

---

## Cierres con verdad conocida

Cada módulo trae el suyo (`python replica_paper/<archivo>.py`). Los dos que
importan:

**`metodo.py`** — un borde: NCC da el intervalo, el SA recupera `B` con **+0.6%**
y `σ` con −0.07%.

Dos bordes solapados, con la ventana y el σ de plantilla **fijos** (si no se
estaría midiendo también la elección de ventana) y sobre 6 realizaciones:

| σ/B | separación | Ec. 12 | error en B₁ | error en D |
|---|---|---|---|---|
| 0.08 | 2.01 σ | corregida | +2.29 ± 1.89% | **−13.60 ± 15.56%** |
| 0.08 | 2.01 σ | impresa | +2.14 ± 1.78% | −7.51 ± 13.69% |
| 0.04 | 4.01 σ | corregida | +0.12 ± 0.71% | **−0.57 ± 2.95%** |
| 0.04 | 4.01 σ | impresa | −0.20 ± 0.61% | −2.57 ± 2.66% |

Dos cosas de acá:

* **A 2σ de separación el modelo de dos bordes no sirve**: `D` sale sesgada −14%
  con ±15% de dispersión, con cualquiera de las dos versiones de la Ec. 12. Es
  exactamente el caso del Co-60 de este proyecto (2.2-2.5 σ).
* A 4σ sí funciona, y ahí aparece la errata: la diferencia entre las dos
  versiones es ~2%, del orden del 3% que predice `(159.3−154.7)/154.7`.

Donde la errata de la Ec. 12 se ve limpia es en los **datos publicados** (§3
arriba): 74 canales medidos contra 76.0 de la fórmula impresa y 73.8 de la
corregida.

**`correr.py --cierre`** — dos fuentes (Cs-137 y Co-60) al **mismo** punto de
operación, con `D` conocida, que es la situación que hoy no está medida. Recupera
la pendiente:

| método | m ajustada | error contra la verdad |
|---|---|---|
| NCC-SA | 0.7697 | **+0.06%** |
| Diferencia | 0.7475 | −2.82% |
| Gauss FWHM | 0.6390 | −16.93% |

O sea la afirmación central del paper —que los dos métodos convencionales sesgan
y el NCC-SA no— **se reproduce**, y sobre Co-60 (bordes solapados) el método
gaussiano se va +11.3% mientras el NCC-SA queda en +0.29%.

---

## Resultados sobre nuestras corridas

`python replica_paper/correr.py` → `replica_paper/informe.txt`

### Cs-137, Detector A @ 1600 V — una línea, el caso que sí funciona

NCC: intervalo `[579, 1389]` mV·ns, ρ = 0.99618 (σ_rel de S2 = 0.16).
SA: **B = 915.2 mV·ns**, σ = 115.5 (σ/B = 0.126 → 29.7% FWHM), χ²/ndf = 2.02.
Con el mapa por el origen de la Ec. 8, **D = 1.917 mV·ns/keV**.

| método | borde [mV·ns] | contra NCC-SA | SD entre 10 bloques |
|---|---|---|---|
| NCC-SA | 915.2 | — | 24.4 (2.67%) |
| Gauss FWHM | 951.7 | +4.0% | 24.3 (2.56%) |
| Diferencia | 829.0 | −9.4% | 139.4 (16.81%) |

Cada bloque tiene 1/10 de las cuentas, así que su dispersión estadística es
~√10 veces la de una medición completa; escalada, el NCC-SA queda en **0.84%** y
el método de la diferencia en 5.3%. Las diez mediciones del paper tienen cada
una la estadística completa **y** varían además en montaje, posición de la
fuente y ganancia, así que ninguna de las dos columnas es directamente su barra.

Robustez frente al intervalo (los seis de la Fig. 5, escalados al borde), con el
mismo criterio en los dos lados:

| | NCC-SA acá | el paper (sus 6 paneles) |
|---|---|---|
| todos los intervalos | 36.1% | 1.24% |
| sin el más angosto | **1.31%** | **0.33%** |

En los dos casos el intervalo más angosto es el que se aparta, y ahí está toda la
diferencia: el suyo tiene ~100 bins y el nuestro **30**. Sacándolo, la
reproducción es del mismo orden (1.3% contra 0.33%, cuatro veces peor con ~20
veces menos cuentas por bin).

El χ²/ndf = 2.02 no es del método: es que la plantilla Klein-Nishina no contiene
el exceso de cuentas a bajo canal que este detector tiene (documentado en
`CLAUDE.md` §7.32 y §7.35). El paper no modela nada de eso porque sus espectros
no lo muestran.

### Co-60, Detector A @ 1550 V — dos bordes solapados, el caso que no funciona

NCC: intervalo `[952, 1967]`, ρ = 0.99095 (σ_rel = 0.25).

| variante | B₁ | B₂ | D | A₂/A₁ | separación | χ²/ndf |
|---|---|---|---|---|---|---|
| Θ de la Ec. 11 (D libre) | 1148.4 | 1367.3 | 1.4149 | 1.55 | 2.16 σ | 1.37 |
| D = B₁/E_max,1 (lo que impone la Ec. 8) | 1124.2 | 1304.7 | **1.1669** | 6.26 | 2.46 σ | 1.18 |
| Ec. 12 impresa (ΔE de las gammas) | 1170.3 | 1363.8 | 1.2147 | — | — | 1.20 |

**D cambia un 21% según se deje libre o no**, con χ² casi igual, y A₂/A₁ toma
valores sin sentido físico (el Co-60 emite las dos líneas con intensidades
prácticamente iguales). El SD entre los 10 bloques es de **457 mV·ns sobre 1124,
un 41%** (12.9% escalado a la estadística completa), y la dispersión frente al
intervalo es del 15-46%. No es el
optimizador: es que Θ = {A₁, A₂, B₁, σ₁, σ₂, D} está sub-determinado a esta
resolución.

En cambio los dos métodos convencionales son *estables* acá (Gauss FWHM: SD 7.4,
dispersión de intervalo 1.3%) **y están sesgados**: el cierre con verdad conocida
los pone en −16.9% y −2.8% de pendiente. Poca dispersión no es exactitud, y esa
distinción es justamente lo que la Fig. 4 del paper mezcla en un solo gráfico.

### Comparación con `analisis/` [fuera del paper]

Para dar contexto, no como parte de la réplica:

| corrida | esta carpeta (paper puro, Klein-Nishina) | `analisis/` (con exceso modelado) |
|---|---|---|
| Cs-137 @ 1600 V | D = 1.917 | D = 1.9450 ± 0.0038 (−1.4%) |
| Co-60 @ 1550 V | D = 1.167 (D ligada) … 1.415 (D libre) | D = 1.2512 ± 0.0047 |

El Cs-137 cierra al 1.4%. El Co-60 no cierra, y la razón está arriba: por la ruta
del paper ese punto no es utilizable.

---

## Lo que NO se puede replicar con lo medido hoy

Está listado con el motivo al final de `correr.py`. En corto:

1. **Fig. 3 (recta canal-energía con SD y R²).** El paper la arma con 4 fuentes
   medidas con el mismo detector y la misma ganancia (5-6 bordes por recta).
   Nosotros tenemos **una fuente por punto de operación**: Co-60 a 1550 V y
   Cs-137 a 1600 V. Con un borde no hay recta; con dos, la recta pasa exacta y
   `SD = 0`, `R² = 1` no miden nada.
   Detalle estructural que conviene no perder: la Ec. 8 mapea canal a energía
   **por el origen**, así que todos los bordes de **una misma fuente** caen
   exactos sobre una recta por el origen. La recta sólo mide algo si los puntos
   vienen de fuentes **distintas**. Por eso el código agrupa por punto de
   operación y está listo para cuando existan Cs-137 y Co-60 al mismo voltaje.
2. **Fig. 4 (error relativo por método).** Depende de la recta del punto 1. La
   **barra** de error sí se calcula (10 bloques).
3. **La comparación contra la calibración por fotopico** y su similitud coseno:
   un centellador plástico no tiene fotopico. El paper tampoco la hace para
   BC408/EJ309, por lo mismo.
4. **Las 10 mediciones independientes.** Acá son 10 bloques de una misma corrida,
   o sea sólo estadística de conteo. La barra que sale es una **cota inferior**
   de la del paper, que incluye montaje, posición de la fuente y deriva.
5. **El segundo borde del Co-60 (1118.1 keV).** Ver abajo.

---

## El hallazgo propio de correr esto sobre nuestros datos

**Θ = {A₁, A₂, B₁, σ₁, σ₂, D} (Ec. 11) sobre-parametriza el modelo respecto de la
Ec. 8 del propio paper.** La Ec. 8 define el mapa canal-energía **por el
origen**, o sea `D = B₁/E_max,1` queda fijado por `B₁`. Dejar `D` libre además es
un grado de libertad que el modelo no tiene.

Se ve en los datos: sobre el Co-60, `D` cambia un **21%** (1.167 → 1.415) según
se ligue o no, con χ²/ndf 1.18 contra 1.37, y el cociente de amplitudes toma
valores sin sentido físico (A₂/A₁ = 1.5 o 6.3, cuando el Co-60 emite las dos
líneas con intensidades prácticamente iguales).

El cierre con verdad conocida de `metodo.py` dice cuándo el modelo de dos bordes
es utilizable y cuándo no, y la variable que manda es la **separación de los dos
bordes en unidades de σ**:

* a **2σ** (lo que se mide acá: 2.2-2.5 σ) `D` sale sesgada −14% con ±15% de
  dispersión → el parámetro no está determinado, y de hecho el SD de B₁ entre
  los 10 bloques del Co-60 es del **41%**;
* a **4σ** `D` sale −0.6% ± 3.0% → ahí sí sirve, y ahí es donde se distingue la
  errata de la Ec. 12.

Por eso al paper le funciona y acá no: en su BC408 (Fig. 5) el flanco da
σ/B ≈ 0.06, que para el Co-60 son ~2.7σ de separación **con estadística ~10
veces mayor**. La conclusión práctica es que el Co-60 no es una buena fuente de
calibración por esta ruta salvo con resolución claramente mejor; una fuente de
una sola línea (Cs-137, Mn-54) sí lo es, y ahí el método reproduce la robustez
que el paper reclama (1.31% de dispersión frente al intervalo, contra 0.33% del
paper con el mismo criterio y mucha más estadística).
