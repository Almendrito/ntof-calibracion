# Plan de `desde_cero/`: recrear la metodología sobre nuestros datos

Escrito el 27-sep-2026. Continúa desde la Etapa 0 (`lector.py`, ✅). Cada
etapa es **una sesión de Claude Code**: se pega el prompt de la etapa, se
revisa lo que entrega, se tacha la casilla y se anota la fecha.

## Qué se quiere lograr

Reconstruir, entendiendo cada paso, la cadena completa que va de un archivo del
DRS4 a la ganancia $D$ en mV·ns/keV, con los datos medidos de Co-60 y Cs-137:

1. primero **el método del paper tal cual** (Kong et al. 2025: plantilla
   Klein-Nishina, NCC, recocido simulado), aplicado a nuestros espectros;
2. después **las adaptaciones** que este proyecto mostró necesarias (plantilla
   de Geant4, barrido de un solo parámetro, χ² con pesos), midiendo cuánto
   cambia $D$ con cada una;
3. al final, **el lote completo** y los dos tests (ley de ganancia y acuerdo
   entre fuentes).

`analisis/` es el **oráculo**: no se importa en `desde_cero/` salvo en los
archivos `verificar_*.py`, que comparan contra él. `replica_paper/` es el
oráculo del paso 1 (el paper literal).

## Reglas (vienen de `README.md` y `CLAUDE.md` §4; valen para todas las etapas)

- Una etapa = un archivo = una idea. Cada archivo corre solo
  (`python desde_cero/<archivo>.py`) y deja una figura en `desde_cero/figuras/`.
- Nada entra sin comprobarse **dos veces**: contra el oráculo (mismo resultado
  sobre los mismos datos) y contra **datos fabricados con verdad conocida**.
- Lo mínimo que funcione: sin cachés, sin interfaz, sin barras de progreso,
  hasta que algo lo exija (si una etapa se vuelve lenta, se decide ahí).
- El porqué va en el código, al lado de la decisión. Comentarios en español.
- Al terminar cada etapa: actualizar la tabla de `README.md` y agregar al final
  de este archivo qué se encontró (sobre todo lo que no coincidió con el
  oráculo y por qué).
- Sin preguntas de comprensión al cierre (pedido explícito, `CLAUDE.md` §4).

## Corridas de trabajo

Siempre las mismas tres, para poder comparar etapa con etapa (claves del
catálogo `analisis/rutas.py`):

| clave | qué es | por qué |
|---|---|---|
| `cs137_1600` | Cs-137, detector A, 1600 V, XML, con fondo | la corrida de ejemplo del informe; un solo borde |
| `co60` | Co-60, detector A, 1550 V (junio), XML, con fondo | dos bordes solapados; registro largo |
| `a_co60_1600` | Co-60, detector A, 1600 V, binario (1-sep) | formato binario, sin fondo, registro corto tras el onset |

---

## Etapa 1 — `pulso.py`: mirar un pulso · [x] 27-sep-2026

**Qué resuelve.** De `(t, v)` de un evento, sacar línea base, ruido, polaridad,
pico, amplitud y onset. Todavía no se integra nada.

**Conceptos.** Mediana contra media para la línea base; desviación robusta
(1.4826·MAD); polaridad deducida, no supuesta; onset como cruce del 10 % de la
amplitud **retrocediendo desde el pico**; línea base relativa al onset (la
trampa de `CLAUDE.md` §7.33: una ventana fija se comía la subida del pulso).

**Oráculo.** `analisis/pulso.py`: `preparar`, `linea_base`, `detectar_pico`,
`detectar_onset`. **Verdad conocida:** pulsos fabricados (forma de la plantilla
de `analisis/.cache/plantilla_pulso_ago-sep.npz` × amplitud + ruido gaussiano +
jitter del onset) con onset y amplitud conocidos.

**Terminado cuando:** sobre 2000 eventos reales de cada corrida de trabajo,
onset y amplitud coinciden con el oráculo (onset a menos de un paso de muestreo,
amplitud a 1e-6 relativo); sobre los fabricados, el onset recuperado tiene
sesgo menor que 0.1 ns y la dispersión que corresponde al ruido.

```text
Etapa 1 de desde_cero. Lee desde_cero/README.md, desde_cero/PLAN.md (reglas y
Etapa 1) y CLAUDE.md §7.33. Crea desde_cero/pulso.py: dada una forma de onda
(t, v) de lector.eventos, calcula línea base (mediana de las muestras antes del
onset con una guarda de 4 ns), ruido (1.4826·MAD), polaridad (deducida de la
excursión máxima respecto de la mediana del registro), pico, amplitud y onset
(último cruce del 10 % de la amplitud retrocediendo desde el pico). Explica en
comentarios por qué mediana y no media, por qué 1.4826, y por qué se retrocede
desde el pico. Agrega desde_cero/verificar_pulso.py que compare contra
analisis/pulso.py sobre 2000 eventos de cs137_1600, co60 y a_co60_1600, y una
prueba con pulsos fabricados de onset y amplitud conocidos (usa la forma de
analisis/.cache/plantilla_pulso_ago-sep.npz). El demo del módulo dibuja 6
pulsos reales con línea base, onset y pico marcados. Reporta las diferencias
máximas con el oráculo y el sesgo en los fabricados.
```

## Etapa 2 — `area.py`: el observable · [x] 27-sep-2026

**Qué resuelve.** El número que se le asigna a cada pulso: el área en una
compuerta `[onset − pre, onset + post]` anclada al onset de cada evento.

**Conceptos.** Integral trapezoidal sobre el eje de tiempo real (el DRS4 no
muestrea uniforme); pulso promedio alineado y su curva de carga acumulada;
convergencia (dónde el pulso se hunde en el ruido); por qué una compuerta
**común** de 12 ns para comparar corridas (`CLAUDE.md` §7.53, §7.66); banderas
(saturado, truncado, sin pre-pulso, sin pulso) que se marcan y no se descartan
acá; el control área/amplitud contra amplitud.

**Oráculo.** `analisis/area.py:areas_de_archivo(ruta, metodo="global", pre=,
post=)` y `calibrar_ventana`. **Verdad conocida:** los pulsos fabricados de la
Etapa 1, cuya integral se conoce analíticamente, con registros truncados.

**Terminado cuando:** con el mismo `pre`/`post`, las áreas coinciden con el
oráculo evento a evento (diferencia relativa mediana < 1e-6); el control
área/amplitud sale plano sobre la amplitud mediana; sobre fabricados, el área
es la integral verdadera por la fracción de carga de la compuerta.

```text
Etapa 2 de desde_cero. Lee PLAN.md (Etapa 2) y CLAUDE.md §7.33, §7.52, §7.53 y
§7.66. Crea desde_cero/area.py usando desde_cero/pulso.py: área de cada evento
en la compuerta [onset − pre, onset + post], integral trapezoidal sobre el eje
de tiempo REAL, línea base restada, signo positivo, y banderas (saturado a 495
mV, compuerta fuera del registro, sin pre-pulso, sin pulso) que se marcan pero
no se descartan. Agrega una función que mida el pulso promedio alineado por
onset y su carga acumulada, y explica en el código por qué el informe usa una
compuerta común de 12 ns y no la que pide cada archivo. Crea
verificar_area.py: compara contra analisis/area.py:areas_de_archivo con el
mismo pre/post sobre 5000 eventos de cada corrida de trabajo, y valida con
pulsos fabricados de integral conocida, incluidos registros truncados. El demo
dibuja el pulso promedio con la compuerta, la carga acumulada, y el cociente
área/amplitud contra amplitud.
```

## Etapa 3 — `espectro.py`: del conteo al espectro neto · [x] 27-sep-2026

**Qué resuelve.** Histograma de áreas, tiempo vivo, resta de fondo y barra de
error por bin.

**Conceptos.** Binning común para muestra y fondo; tope del histograma por
percentil y no por máximo; tiempo vivo desde el scaler de la placa (recortado a
10 veces su mediana por las ráfagas) y no desde el reloj de pared
(`CLAUDE.md` §7.34); factores $k_m$, $k_f$; varianza de Poisson propagada
$k_m^2N_m + k_f^2N_f$ con piso 1; bines negativos que se conservan; el fin de
la subida del umbral como piso de cualquier ventana (§7.48).

**Antes de empezar:** `lector.eventos` entrega sólo `(t, v)`. Esta etapa
necesita también el reloj y el scaler de cada evento, así que hay que extender
el lector (y su verificación) sin romper lo que ya da.

**Oráculo.** `analisis/espectro.py:espectro_neto`, `tiempo_vivo`,
`fin_del_umbral`. **Verdad conocida:** muestra y fondo fabricados con tasas y
tiempos vivos conocidos.

**Terminado cuando:** el espectro neto y su error coinciden bin a bin con el
oráculo sobre `cs137_1600` y `co60`; el tiempo vivo coincide al 0.1 %; en el
fabricado, la resta devuelve el espectro de muestra verdadero dentro de su
barra.

```text
Etapa 3 de desde_cero. Lee PLAN.md (Etapa 3) y CLAUDE.md §7.34, §7.37 y §7.48.
Primero extiende desde_cero/lector.py para que pueda entregar también reloj y
scaler de cada evento (sin cambiar lo que ya devuelve por defecto) y amplía
verificar.py para comprobarlo contra analisis/drs_xml y drs_bin. Después crea
desde_cero/espectro.py: histograma con binning común para muestra y fondo
(tope = 1.6 × percentil 99.5 de la muestra, 500 bines), tiempo vivo desde el
scaler recortado a 10 veces su mediana, resta con k_m y k_f llevados al menor
tiempo vivo, varianza de Poisson propagada con piso 1, sin recortar negativos,
y el fin de la subida del umbral como primer máximo del espectro suavizado.
Explica en comentarios por qué el reloj de pared no sirve. Crea
verificar_espectro.py contra analisis/espectro.py sobre cs137_1600 y co60, y
una prueba con muestra y fondo fabricados de tasas conocidas. El demo dibuja
muestra, fondo escalado y neto con su barra.
```

## Etapa 4 — `fisica.py`: el borde y la forma teórica (paper, Ec. 1, 2 y 10) · [x] 27-sep-2026

**Qué resuelve.** La energía del borde Compton y la forma teórica del continuo,
sin nada del detector.

**Conceptos.** $E_{\max} = 2\alpha E_\gamma/(1+2\alpha)$ y la errata de la
Ec. 1 impresa; Klein-Nishina en ángulo (Ec. 2) y en energía del electrón
(Ec. 10) por el jacobiano; la comprobación de que las dos integran a la misma
sección eficaz total; el Co-60 como dos líneas.

**Oráculo.** `replica_paper/ecuaciones.py` (`e_max`, `dsigma_domega`,
`dsigma_de`). **Verdad conocida:** la sección eficaz de Thomson en el límite
de baja energía (0.665 b) y la integral de las dos formas.

**Terminado cuando:** bordes de 477.3, 963.4 y 1118.1 keV; las dos integrales
coinciden a 6 decimales para 500, 662, 1000 y 2000 keV; la figura reproduce la
Fig. 1 del paper.

```text
Etapa 4 de desde_cero. Lee PLAN.md (Etapa 4), la Sección II del paper
(docs/Paper ref.pdf, Ec. 1, 2 y 10) y CLAUDE.md §3 y §7.1. Crea
desde_cero/fisica.py con: E_max correcta y la impresa (para mostrar la errata),
dσ/dΩ de Klein-Nishina (Ec. 2), dσ/dE del electrón (Ec. 10) deducida por el
jacobiano (escribe la deducción en el docstring), y las líneas de Cs-137 y
Co-60. Verifica que integrar Ec. 2 en ángulo sólido y Ec. 10 en energía da la
misma sección eficaz total a 6 decimales, y el límite de Thomson. Crea
verificar_fisica.py contra replica_paper/ecuaciones.py. El demo redibuja la
Fig. 1 del paper para las tres líneas.
```

## Etapa 5 — `modelo.py`: el modelo directo (paper, Ec. 8 y 9) · [x] 27-sep-2026

**Qué resuelve.** Llevar la forma teórica al eje medido y borronearla con la
resolución: el espectro que el ajuste compara contra el dato.

**Conceptos.** Mapa por el origen $ch = D\,E$ (Ec. 8); kernel gaussiano
normalizado a suma 1 (Ec. 9) y por qué así $A$ y $\sigma$ quedan desacoplados;
**evaluar la convolución en una grilla más ancha que la ventana** (si no, el
ajuste achica $\sigma$, `CLAUDE.md` §7.14); el Co-60 con los dos bordes atados
por la misma $D$ y $\Delta E$ entre **bordes** (errata de la Ec. 12, §7.15).

**Oráculo.** `replica_paper/ecuaciones.py:modelo` y `modelo_dos_bordes`.
**Verdad conocida:** conservación de cuentas al convolucionar; el 50 % de
altura corriéndose del borde al crecer $\sigma$.

**Terminado cuando:** coincide con el oráculo a 1e-10; conserva el área; la
figura muestra el sesgo del 50 % de altura contra $\sigma$.

```text
Etapa 5 de desde_cero. Lee PLAN.md (Etapa 5), la Sección II-A del paper
(Ec. 8, 9, 11-13) y CLAUDE.md §7.2, §7.14 y §7.15. Crea desde_cero/modelo.py:
la forma de desde_cero/fisica.py llevada a canal con ch = D·E (mapa por el
origen), convolucionada con una gaussiana normalizada a suma 1 de ancho
σ_ch = σ_rel·D·E_ref, evaluada sobre una grilla extendida hacia abajo y hacia
arriba y recortada después. Para Co-60, las dos líneas con la misma D (así ΔE
entre bordes queda fijado solo; explica la errata de la Ec. 12). Explica por
qué el kernel suma 1 y por qué hace falta la grilla extendida. Crea
verificar_modelo.py contra replica_paper/ecuaciones.py y con pruebas de
conservación de cuentas. El demo dibuja Cs-137 y Co-60 para tres σ y marca el
punto del 50 % de altura contra el borde verdadero.
```

## Etapa 6 — `ncc.py`: la ventana automática (paper, Ec. 15, Sección II-B) · [x] 27-sep-2026

**Qué resuelve.** Elegir sin intervención el tramo del espectro que se ajusta,
y de paso la semilla de $D$ y las cotas de la búsqueda.

**Conceptos.** Coeficiente de Pearson y su invariancia a escala y pedestal;
suavizado Savitzky-Golay; el barrido del paper (inicio y fin libres, 100
puntos) contra el nuestro (sólo $B$, ventana fija por el origen, un punto por
bin) y por qué el primero elige ventanas degeneradas (`CLAUDE.md` §7.10, §7.11,
§7.27); semilla $D_0 = B/E_{\text{ref}}$ y cotas $[0.6, 1.6]\,D_0$.

**Oráculo.** `replica_paper/ncc.py` (el barrido del paper) y
`analisis/ncc.py:buscar_intervalo` (el nuestro). **Verdad conocida:** espectro
fabricado con el modelo de la Etapa 5 y $B$ conocido.

**Terminado cuando:** sobre el fabricado, $B$ recuperado a menos de un bin; sobre
las corridas reales, la ventana y $\rho$ coinciden con `analisis/ncc.py` usando
la misma plantilla; la figura muestra un caso donde el barrido libre del paper
elige una ventana degenerada.

```text
Etapa 6 de desde_cero. Lee PLAN.md (Etapa 6), la Sección II-B del paper y
CLAUDE.md §7.8-7.11 y §7.27. Crea desde_cero/ncc.py con: Pearson (Ec. 15),
suavizado Savitzky-Golay, el barrido del paper (inicio y fin libres,
reinterpolado a 100 puntos) y el nuestro (sólo B, ventana [0.45·B, f·B] con f
que deja entrar el último borde, un punto por bin), devolviendo ventana, ρ,
semilla D0 = B/E_ref y cotas. Usa la plantilla Klein-Nishina de modelo.py.
Explica por qué barrer sólo B y por qué no remuestrear a pocos puntos. Crea
verificar_ncc.py: con un espectro fabricado de B conocido, y contra
replica_paper/ncc.py y analisis/ncc.py sobre cs137_1600 y co60. El demo dibuja
ρ(B), la ventana ganadora y las dos curvas estandarizadas.
```

## Etapa 7 — `ajuste.py`: recocido simulado y la barra de error (paper, Sección II-A) · [x] 27-sep-2026

**Qué resuelve.** El mínimo global del $\chi^2$ en $(D, \sigma)$ y su barra.
Con esta etapa la réplica del paper sobre nuestros datos está completa.

**Conceptos.** $\chi^2$ con barras de Poisson (Neyman) y por qué pesar;
amplitud despejada en forma cerrada (VARPRO); Metropolis, enfriamiento
geométrico ($\alpha = 0.95$, 100 propuestas por temperatura, paso 0.1 del
ancho de la cota); $T_{\text{init}}$ **medida** y no copiada (§7.22);
multiarranque de 8 desde puntos al azar en las cotas; pulido con Nelder-Mead;
hessiano por diferencias finitas, $C = 2H^{-1}$, escalado por
$\sqrt{\chi^2/\text{ndf}}$; el cierre con verdad conocida.

**Oráculo.** `replica_paper/metodo.py:ajustar` (paper, Klein-Nishina) y
`replica_paper/informe.txt` (sus números). **Verdad conocida:** 30
realizaciones Poisson del modelo ajustado, pasadas por NCC + recocido.

**Terminado cuando:** sobre el fabricado, sesgo de $D$ compatible con cero y
cobertura cerca del 68 %; sobre `cs137_1600` y `co60`, $D$ con plantilla
Klein-Nishina coincide con `replica_paper` dentro de su barra; los 8 arranques
terminan en el mismo mínimo.

```text
Etapa 7 de desde_cero. Lee PLAN.md (Etapa 7), la Sección II-A del paper,
CLAUDE.md §7.21-7.23, §7.30, §7.31 y el apéndice B del informe
(informe/informe.tex, sec:est-chi2 y sec:est-incert). Crea desde_cero/ajuste.py:
χ² de Neyman en la ventana del NCC con la amplitud despejada en forma cerrada;
recocido con Metropolis, enfriamiento geométrico α = 0.95, 100 propuestas por
temperatura, paso 0.1 del ancho de la cota y T_init medida para una aceptación
inicial de ~80 %; multiarranque de 8 desde puntos uniformes en las cotas;
pulido con Nelder-Mead; barra por hessiano en diferencias finitas centradas,
C = 2H⁻¹, escalada por sqrt(χ²/ndf). Usa la plantilla Klein-Nishina (Etapas
4-6). Crea verificar_ajuste.py: cierre con verdad conocida (30 realizaciones
Poisson, sesgo y cobertura) y comparación contra replica_paper/metodo.py sobre
cs137_1600 y co60. El demo dibuja el ajuste, los residuos normalizados, el
costo por ciclo de los 8 arranques y el camino en el plano (D, σ).
```

## Etapa 8 — `plantilla.py`: de Klein-Nishina a Geant4 (lo que agrega este proyecto) · [x] 27-sep-2026

**Qué resuelve.** Cambiar la plantilla analítica por la simulación del detector
y medir cuánto cambia $D$.

**Conceptos.** Scattering múltiple (un cuarto de los eventos) y por qué corre el
borde; plantilla como histograma de `Edep_keV` interpolado en $ch/D$ (no
rebinneado: aliasing, §7.25); verificar que la simulación es de la fuente
correcta leyendo `E0_keV` (§7.41); escala de energía depositada contra luz
(quenching de Birks, §7.17).

**Oráculo.** `analisis/calibracion.py:calibrar` y
`analisis/figuras/calibraciones/resultados.json` (el resultado del informe).

**Terminado cuando:** con la plantilla de Geant4, $D$ de `cs137_1600` y
`co60` coincide con el oráculo dentro de su barra; el cambio respecto de
Klein-Nishina queda medido y explicado (del orden de +4 %, §7.18).

```text
Etapa 8 de desde_cero. Lee PLAN.md (Etapa 8), CLAUDE.md §7.17, §7.18, §7.25,
§7.41 y la sección sec:geant4 del informe. Crea desde_cero/plantilla.py: lee el
.root de la fuente (uproot, rama Edep_keV, 350 bines), verifica que las líneas
primarias (E0_keV) sean las de la fuente declarada y se niega a seguir si no,
y entrega la plantilla para modelo.py, que debe aceptarla en lugar de la
Klein-Nishina sin otros cambios. Repite el ajuste de la Etapa 7 con las dos
plantillas sobre cs137_1600 y co60 y tabula D, σ y χ²/ndf. Crea
verificar_plantilla.py contra analisis/figuras/calibraciones/resultados.json.
El demo superpone las dos plantillas y los dos ajustes, y separa la
simulación por número de dispersiones Compton.
```

## Etapa 9 — `lote.py`: todas las corridas y los dos tests · [x] 27-sep-2026

**Qué resuelve.** Calibrar las 28 corridas y hacer las dos pruebas del
informe.

**Conceptos.** Ley de ganancia $D \propto V^{b}$ ajustada en log-log con pivote
geométrico, contra la recta (test B); cociente Cs-137/Co-60 al mismo punto de
operación, corregido por quenching (0.980 en escala depositada), y la $D$
conjunta como media pesada (test A); propagación de errores.

**Oráculo.** `analisis/calibrar_todo.py` y su `resultados.json`;
`analisis/voltaje.py` y `analisis/registro.py`.

**Terminado cuando:** las 28 $D$ coinciden con `resultados.json` dentro de su
barra; $b$ por detector y fuente y los desacuerdos del test A coinciden con
las tablas del informe.

```text
Etapa 9 de desde_cero. Lee PLAN.md (Etapa 9), CLAUDE.md §7.46, §7.53-7.55 y
§7.72, y las secciones sec:detectores, sec:testA y el apéndice ap:propagacion
del informe. Crea desde_cero/lote.py: recorre el catálogo de corridas
(analisis/rutas.py sólo para las rutas), calibra cada una con la cadena de
desde_cero (Etapas 1-8, compuerta común de 12 ns, plantilla de Geant4), y
hace el test B (ley de potencia en log-log con pivote en la media geométrica,
contra la recta) y el test A (Cs/Co al mismo punto, corregido por quenching,
D conjunta como media pesada). Propaga los errores como en el apéndice
ap:propagacion. Crea verificar_lote.py contra
analisis/figuras/calibraciones/resultados.json y las tablas de
informe/generado/. El demo dibuja D(V) por detector y el cociente Cs/Co.
```

## Etapa 10 (opcional) — `clasicos.py`: la comparación del paper (Fig. 4 y 5) · [x] 27-sep-2026

**Qué resuelve.** Los métodos contra los que compara el paper (50 % de altura,
mínimo de la derivada, escalón gaussiano) y la robustez frente al intervalo.

**Oráculo.** `codigo/etapas/clasicos.py`, `replica_paper/convencionales.py`,
`replica_paper/resultados.py:robustez_intervalo`.

```text
Etapa 10 de desde_cero. Lee PLAN.md (Etapa 10), CLAUDE.md §7.44 y la sección
sec:clasicos del informe. Crea desde_cero/clasicos.py con los tres métodos
convencionales (ventana de suavizado derivada del ancho del flanco, región
anclada al borde del NCC) y el test de robustez de la Fig. 5 (seis intervalos
anidados). Mide su sesgo contra un borde fabricado de posición conocida, con y
sin ruido, y sobre las corridas de trabajo contra la D de la Etapa 8. Verifica
contra codigo/etapas/clasicos.py y replica_paper/. El demo dibuja el sesgo
contra σ/B.
```

---

## Orden y dependencias

| sesión | etapa | depende de | tiempo aproximado de cómputo |
|---|---|---|---|
| 1 | 1 · pulso | 0 | segundos |
| 2 | 2 · área | 1 | ~1 min por XML grande (parseo) |
| 3 | 3 · espectro | 0 (extensión), 2 | ~1 min |
| 4 | 4 · física | — | segundos |
| 5 | 5 · modelo | 4 | segundos |
| 6 | 6 · NCC | 3, 5 | segundos |
| 7 | 7 · ajuste | 6 | minutos (cierre: 30 realizaciones) |
| 8 | 8 · plantilla | 7 | segundos |
| 9 | 9 · lote | 8 | decenas de minutos (28 corridas) |
| 10 | 10 · clásicos | 6, 8 | minutos |

Las Etapas 4 y 5 no dependen de las 1-3: si una sesión de datos se atasca, se
puede adelantar la física.

## Bitácora

(Al cerrar cada etapa: fecha, qué se hizo, qué no coincidió con el oráculo y
por qué.)

- 27-sep-2026 — plan escrito. Etapa 0 ya estaba hecha (`lector.py`,
  `verificar.py`).
- 27-sep-2026 — **Etapa 1 cerrada** (`pulso.py`, `verificar_pulso.py`,
  `figuras/etapa1_pulsos.png`).
  - De paso: `lector._archivos_de_ejemplo` apuntaba a `datos/<xml>` y a `D:\`,
    que ya no existen desde §7.61, así que `verificar.py` no encontraba nada.
    Ahora `lector.CORRIDAS_TRABAJO` tiene las rutas de las tres corridas de
    trabajo (mismas claves que `analisis/rutas.py`, sin importarlo) y la
    Etapa 0 vuelve a verificar en cero.
  - **Contra el oráculo:** onset, amplitud, base, ruido y banderas idénticos
    (diferencia 0.0) en 2000 eventos de `cs137_1600`, `co60` y `a_co60_1600`.
    Es la misma cuenta en dos pasadas (mediana del registro → onset
    provisorio → base antes del onset con guarda de 4 ns → pico y onset
    definitivos); con la cola de respaldo en 80 ns, igual que el oráculo con
    `post=80`.
  - **Verdad conocida** (forma medida del pulso, eje temporal real del DRS4,
    ruido 1 mV, onset con jitter en 40-70 ns, polaridad al azar): polaridad
    100 %; onset interpolado con sesgo +0.002 ns a 100 mV y −0.005 a 300 mV;
    dispersión 0.061 / 0.021 ns contra 0.058 / 0.019 predichos por
    `sqrt(2/3)·σ/(A·f'(0))` (el `sqrt(2/3)` es el ruido promediado de la
    interpolación entre dos muestras).
  - **Lo que no cumple el criterio, y por qué:** a amplitud/ruido = 10 el
    onset sale **+0.36 ns tarde**. Separado a mano: +0.30 ns porque el flanco
    cruza el umbral (1σ) varias veces y retroceder desde el pico elige el
    último cruce, y +0.06 ns porque la amplitud sale sesgada hacia arriba
    (+1.4 mV, el máximo de muestras ruidosas). Buscar hacia adelante daría
    **−3.97 ns**: es la razón cuantitativa de retroceder desde el pico. El
    criterio de sesgo < 0.1 ns se exige desde amplitud/ruido ≥ 30 (a 30 da
    +0.064).
  - **Para la Etapa 2:** en los datos reales el 20-25 % de los eventos tiene
    amplitud/ruido < 30 (medianas 50, 112 y 81). Un onset tarde corre la
    compuerta y cambia el área de forma dependiente de la amplitud, o sea una
    no linealidad. Hay que medir cuánto con los fabricados, y en qué zona del
    espectro caen esos eventos (se espera: lejos del borde).
  - El onset "de muestra" (el que usa el oráculo para anclar la compuerta)
    tiene sesgo ≈ −dt/2 ≈ −0.1 ns por construcción. Se conserva igual que el
    oráculo; `t_onset_fino` queda disponible si la Etapa 2 lo necesita.
- 27-sep-2026 — **Etapa 2 cerrada** (`area.py`, `verificar_area.py`,
  `figuras/etapa2_area.png`). Compuerta `[onset − 2, onset + 12]` ns.
  - **Encontrado: la espiga del borde del DRS4 le gana a los pulsos chicos,
    también en `analisis/`.** Las primeras celdas del registro dan una espiga
    de 4 muestras (p99 de |v − mediana|: 60, 44, 29, 15 mV, después 2-4 mV).
    Con el recorte de 2 muestras de `analisis/pulso.py` quedan las de 29 y
    15 mV, y en un evento con pulso de menos de ~30 mV la espiga es el máximo:
    el onset cae en τ = 0, el pulso verdadero se ignora, área ≈ 0 y bandera
    SIN_PRE|TRUNCADO (o sea, descartado). Eran 30 de 5000 eventos en
    `cs137_1600` y 29 en `a_co60_1600`, todos de amplitud baja: un corte que
    depende de la amplitud, en la zona baja del espectro, lejos del borde.
    `pulso.RECORTE` pasó a **4**; con eso quedan 0. A los oráculos se les pasa
    `recorte=pulso.RECORTE` (ambos lo aceptan), así que la comparación sigue
    siendo exacta. `analisis/` NO se tocó.
  - **Contra el oráculo** (mismo pre/post/recorte, 5000 eventos, sin caché):
    áreas y banderas idénticas bit a bit fuera de los eventos SIN_PRE.
  - **Diferencia a propósito con el oráculo:** sin pre-pulso, la línea base
    sale de la cola. El oráculo la empieza en `onset + post + 4` = +16 ns con
    la compuerta de 12, donde el pulso vale todavía ~5 % y tiene una cola de
    ~1.5 %; acá en +84 ns (`pulso.COLA`), que no depende de la compuerta.
    Fabricados: −1.06 / −0.90 % el oráculo contra +0.03 / −0.16 % el nuestro
    (30 y 100 mV). Afecta a los 89 eventos SIN_PRE reales del Co-60 de junio
    (pulsos tempranos de verdad): +0.54 % de área mediana respecto del oráculo.
  - **Verdad conocida:** sin ruido, sesgo +0.02 % y dispersión 0.30 % (los
    bordes de la integral caen en muestras). Con ruido de 1 mV: sesgo del área
    +0.26 % a 30 mV, +0.03 % a 100, +0.02 % a 300; registros truncados: la
    bandera prende siempre y el área sale −33.9 % contra −33.9 % predicho.
  - **Cierra lo pendiente de la Etapa 1:** el onset tarde a bajo SNR sube el
    área **+1.2 % a SNR 10 y +2.2 % a SNR 5**. Es una no linealidad real, pero
    en pulsos de ≤ 10 mV, que con estas ganancias caen muy por debajo del
    borde Compton.
  - **El control área/amplitud miente a bajo SNR, y no por el área:** en
    fabricados el área sube 1.2 % a SNR 10 mientras área/amplitud baja
    **−10.8 %**, porque la amplitud (un máximo de muestras ruidosas) sale
    sesgada hacia arriba. La caída que se ve en el Co-60 real a 5 mV
    (6.1 contra 7.0 ns, ruido 0.44 mV → SNR ≈ 11) es eso. Por eso el control
    se lee sólo sobre la amplitud mediana.
  - **Datos reales:** la compuerta recoge 81.3 / 83.3 / 82.3 % de la carga a
    40 ns (§7.66 decía 81-83 %). Área/amplitud sobre la mediana varía 0.9 %
    en los dos Co-60 y **3.5 % en `cs137_1600`**, con un máximo en ~60 mV que
    no está explicado (§7.66 ya veía hasta 4.3 %). Queda anotado; no lo
    produce ninguno de los efectos medidos acá.
  - Descartes de pico a pico ≤ 3 mV: se marcan con SIN_PULSO en vez de
    tirarlos (el oráculo los tira); tramas repetidas sí se tiran (0 en las
    tres corridas).
- 27-sep-2026 — **Etapa 3 cerrada** (`espectro.py`, `verificar_espectro.py`,
  `figuras/etapa3_espectro.png`; `lector.eventos(..., con_reloj=True)` y
  `verificar.py` ampliado).
  - **Lector:** reloj y scaler idénticos al oráculo en los 5 archivos de las
    corridas de trabajo (XML con CHN2+Scaler1 y CHN1+Scaler0, y binario).
    `area.areas_de_archivo` guarda reloj y scaler de TODOS los eventos leídos
    (repetidos incluidos), porque el tiempo vivo cuenta disparos, no pulsos
    útiles.
  - **Contra el oráculo** (mismas áreas, archivos enteros): canales, neto,
    error, tiempo vivo y fin del umbral idénticos (diferencia 0.0). Se
    reproducen §7.34 y §7.37: Co-60 vivo 282.6 de 699.2 s (40.4 %),
    k_f = 0.2866; Cs-137 vivo 1.7 de 1372.8 s (0.12 %), k_f = 0.00064. Fin del
    umbral: 151.5 mV·ns en `cs137_1600` y 46.2 en `co60`. Con la compuerta de
    12 ns todas las áreas son más chicas que las de los números viejos de
    §7.48 (121-334), así que esos no se comparan.
  - **Verdad conocida** (fuente 150 Hz + fondo 10 Hz, muestra viva el 40 %,
    marcas en ráfagas de 200): tiempo vivo recuperado a ≤ 0.7 %; pulls del
    neto contra las cuentas de fuente esperadas, media −0.011 y ancho 1.022
    en 5 realizaciones × ~315 bines. O sea, la barra de error cubre lo que
    dice cubrir.
  - **Lo que se aprendió del scaler real, y conviene no olvidar:**
    - a tasa baja viene cuantizado en pasos de 10 Hz y la mitad de las
      lecturas son 0. La "mediana" que usa el recorte (sobre los valores > 0)
      da 20 Hz en el fondo de junio, que guarda 10 ev/s;
    - pesado por el tiempo **lee de menos**: 9.8 Hz contra 10.1 ev/s guardados
      en el fondo de junio, y 5.3 Hz contra 7.5 ev/s en el del 26-ago. La
      fracción viva sale 1.04 y 1.42, y el tope en 1 la deja en 100 %. **El
      "fondo 100 % vivo" lo sostiene el tope, no el scaler.** Físicamente es
      razonable (guardar 7-10 ev/s está muy por debajo de los ~70 ev/s que el
      DRSOSC escribe en el Co-60), pero no es una medición;
    - las espigas: el fondo de junio tiene 0.9 % de lecturas a más de 10 veces
      la mediana (hasta 90 kHz). Sin recorte daría 32 % vivo en vez de 100 %.
      Con recorte a 3× o a 10× da lo mismo, porque el tope absorbe el resto.
      El recorte no las elimina: con una fracción p de espigas los disparos
      se inflan en ~9p. Medido en fabricados: −1.6 % con p = 0.2 % contra
      −1.8 % predicho.
    - En las muestras (170-46000 Hz) nada de esto aparece: sin espigas y sin
      ceros.
  - **Un error del fabricado, no del módulo:** redondear los milisegundos
    aparte de los segundos hacía retroceder la marca 1 s cuando la fracción
    era ≥ 0.9995, y el tiempo vivo salía 0.5 % bajo sin espigas. Se redondea
    una sola vez, a ms enteros.
  - Costo: la Etapa 3 sobre las dos corridas con fondo parsea ~5 GB de XML
    en ~2 min. Todavía no pide caché.
- 27-sep-2026 — **Etapa 4 cerrada** (`fisica.py`, `verificar_fisica.py`,
  `figuras/etapa4_klein_nishina.png`).
  - Constantes CODATA 2018 (m_e c² = 510.99895 keV, r₀ = 2.8179403262e-13 cm).
    El paper no da m_e c²; `replica_paper` usa 511 y 2.818e-13 porque con eso
    se redibuja su Fig. 1. Para los bordes da igual (477.33 keV con los dos).
    Líneas de ENSDF: Cs-137 661.657 keV; Co-60 1173.228 y 1332.492 keV, con
    intensidades 0.9985 y 0.999826, o sea los dos bordes pesan lo mismo.
  - **Contra `replica_paper/ecuaciones.py`** (pasándole nuestras constantes):
    Ec. 1 (corregida e impresa), Ec. 2, Ec. 10 y Ec. 12 corregida iguales a
    redondeo de máquina (≤ 7e-16).
  - **Verdad conocida:** bordes 477.33 / 963.42 / 1118.10 keV; la Ec. 10
    escrita y la Ec. 2 pasada a energía por el jacobiano
    (2π m_e c²/E'² · dσ/dΩ, deducción en el docstring de `dsigma_de`) coinciden
    a 2e-15; la sección eficaz total por tres caminos (integral en ángulo,
    integral en energía, fórmula cerrada) coincide a 6 decimales en 500, 662,
    1000 y 2000 keV (0.289166, 0.256140, 0.211208, 0.146359 b); el límite de
    Thomson (0.665246 b) se recupera con la corrección de primer orden
    σ_T(1 − 2α) y un resto de segundo orden; la Fig. 1 del paper, leída del
    eje, coincide al 3 % (5 % en el pico de 5000 keV, que se lee "0.11").
  - Ec. 1 impresa: 184.3 / 209.8 / 214.4 keV, la energía del fotón
    retrodispersado. Ec. 12: 154.7 keV entre bordes contra 159.3 entre líneas.
  - Una afirmación mía que hubo que corregir con números: la fórmula cerrada
    de la σ total sí pierde precisión a energía baja por cancelación, pero
    recién abajo de ~10 eV (8e-8 a 10 eV, 4e-5 a 1 eV, 30 % a 0.01 eV). A
    cientos de keV es exacta a 1e-11.
  - Esta etapa no usa datos: la 5 lleva esta forma al eje medido.
- 27-sep-2026 — **Etapa 5 cerrada** (`modelo.py`, `verificar_modelo.py`,
  `figuras/etapa5_modelo.png`). Modelo `modelo(canales, A, D, sigma_rel, fuente)`.
  - **Decisiones que difieren del paper, a propósito:**
    - **A son cuentas.** La forma se normaliza a integral 1 en energía y cada
      canal recibe A·forma(ch/D)·paso/D. En el paper A absorbe el prefactor
      de la Ec. 10 y no significa nada; acá la suma del modelo es A, y eso es
      una prueba.
    - **Un solo σ en canales**, σ_ch = σ_rel·D·E_ref con E_ref el borde más
      bajo, para los dos bordes del Co-60. El paper tiene σ₁ y σ₂ (Ec. 11);
      físicamente σ₂/σ₁ ≈ √(1118/963) = 1.08, que acá se ignora. Motivo:
      §7.23 (con dos σ libres el ajuste los invertía). Si el 8 % importa, lo
      va a decir el χ² de la Etapa 7.
    - **Co-60 con una sola D**, así el segundo borde cae en D·E_max2 sin
      escribir la Ec. 12; las dos líneas pesan su intensidad ENSDF por su
      Klein-Nishina (aproximación de detector delgado).
  - **Contra `replica_paper`** (`modelo` y `modelo_dos_bordes` con la Ec. 12
    corregida, amplitudes convertidas, mismas constantes): iguales a ≤ 2.5e-15
    en Cs-137 y Co-60, a tres (D, σ_rel).
  - **Verdad conocida:** la convolución conserva cuentas (1e-16) para
    cualquier σ; la suma sin resolución da A a +0.13 % (sumar la forma en
    centros de bin); con σ_rel = 0.15, el 6.1 % de las cuentas cae en canales
    < 0 (energía ~0 borroneada), así que la suma sobre canales ≥ 0 es menor
    que A a propósito. Sin grilla extendida el modelo sale **51 % bajo** en
    el borde izquierdo de una ventana [0.45 B, 1.3 B] e idéntico en el medio:
    es §7.14. El 50 % de altura sobrestima el borde +0.8 / +2.0 / +3.5 /
    +4.9 / +7.5 % a σ_rel = 0.03 / 0.06 / 0.09 / 0.12 / 0.18, sin ruido (la
    columna "50 %" de §7.44, con ruido: +0.8 / +2.2 / +3.8 / +5.4 / +8.5).
  - Un defecto de mi propia prueba (d, el segundo borde): buscaba el último
    canal con cuentas > 1e-6 del máximo, que llega hasta la cola del kernel
    (4σ). La tolerancia pasó a 4σ + un paso.
- 27-sep-2026 — **Etapa 6 cerrada** (`ncc.py`, `verificar_ncc.py`,
  `figuras/etapa6_ncc.png`). Dos barridos con la plantilla Klein-Nishina de
  `modelo.py`: `barrido_B` (el nuestro: un parámetro, ventana
  [0.45 B, (E_max,último/E_ref + 3σ) B], Savitzky-Golay 21, un punto por bin)
  y `barrido_libre` (el del paper: inicio y fin libres, 100 puntos, sin
  suavizar). Devuelve la ventana, ρ, la semilla D0 = B/E_ref y las cotas
  [0.6, 1.6]·D0 y σ_rel ∈ [0.02, 0.30] para la Etapa 7.
  - **Contra los oráculos:** `barrido_libre` = `replica_paper/ncc.py` en la
    misma ventana y el mismo ρ (a 1e-12), en el Cs-137 real y en uno
    fabricado. Con el Co-60 no se compara, porque replica_paper arma la
    plantilla con una sola línea. `barrido_B` = `analisis/ncc.py` con la misma
    plantilla: el mismo B en los dos, y ρ a 5e-6 (el oráculo interpola la
    plantilla de una tabla y usa un kernel de 6σ).
  - **Verdad conocida** (espectros de `modelo.modelo` con Poisson, la grilla y
    las cuentas de la corrida real, σ verdadero 0.09): con la plantilla en el
    σ verdadero, B sale a +0.17 ± 0.51 bines (Cs-137) y −0.16 ± 0.66 (Co-60);
    el del paper, −1.3 y −0.55 bines. Con la plantilla a σ = 0.06 o 0.12, el
    nuestro se corre ≤ 2.1 bines y el del paper hasta 3.6.
  - **Sobre los datos reales:** Cs-137 B = 683.3 (D0 = 1.4315, ρ = 0.9942,
    ventana [307, 868]); Co-60 junio B = 960.7 (D0 = 0.9972, ρ = 0.9804);
    Co-60 binario B = 1263.6 (D0 = 1.3116, ρ = 0.9810). El barrido del paper
    da B a 0.3-1.4 % de esos. Las ventanas quedan todas por encima del fin
    del umbral (151, 46 y 125 mV·ns).
  - **No se reprodujo la ventana degenerada de §7.10,** que el plan pedía
    mostrar. Con el observable corregido y 500 bines, el barrido libre elige
    ventanas razonables incluso con un largo mínimo de 6 bines (en vez de
    25). Aquella salió del área mal integrada (espectro ancho, cola vacía
    larga, 800 bines de 10 mV·ns). No se fabricó un caso para mostrarla.
  - **Una hipótesis mía que los datos refutaron:** esperaba que el B del
    barrido libre dependiera más del σ de la plantilla (porque lo deduce del
    span supuesto, 1 + 4.5σ). En los datos reales, con σ entre 0.04 y 0.15,
    los dos se mueven varios %: nuestro 5.7 / 0.3 / 2.4 % y el del paper
    3.4 / 1.5 / 3.9 % (Cs-137 / Co-60 junio / Co-60 binario). Sobre Klein-Nishina
    pura el nuestro es más estable; sobre datos reales, no hay ganador.
  - **Por qué:** por debajo del borde el dato sube hacia canales bajos y la
    Klein-Nishina es plana (panel inferior de la figura). Es el scattering
    múltiple y el exceso, que la plantilla analítica no tiene. Con esa forma
    equivocada, qué σ "se parece más" no está bien definido. Consecuencia: el
    B del NCC es ventana y semilla, **no la calibración**. D sale del ajuste
    (Etapa 7) y el sesgo de forma lo corrige la plantilla de Geant4 (Etapa 8).
- 27-sep-2026 — **Etapa 7 cerrada** (`ajuste.py`, `verificar_ajuste.py`,
  `figuras/etapa7_ajuste.png`). χ² de Neyman en la ventana del NCC (nunca por
  debajo del fin del umbral), A despejada en forma cerrada, recocido
  (α = 0.95, 100 propuestas por ciclo, paso 0.1 de la cota, T_init medida para
  80 % de aceptación, corte relativo 1e-6 durante 10 ciclos), 8 arranques
  uniformes en las cotas, pulido Nelder-Mead, barra C = 2H⁻¹ escalada por
  √(χ²/ndf).
  - **Encontrado y corregido en `modelo.py` (afecta a las Etapas 5 y 6):**
    evaluar la forma en el CENTRO de cada bin (como el paper y
    `replica_paper`) hace el modelo discontinuo en D. Cuando el borde de
    Klein-Nishina cruza el centro de un bin, ese bin salta de 0 al pico. En el
    Cs-137 real, el χ² perfilado en D saltaba **22 unidades** a una fracción de
    barra del mínimo, y la curvatura daba una barra de 0.0004 que no
    significaba nada. `modelo.modelo(..., por_bin=True)` (default) integra la
    forma en cada bin con la acumulada F(E): continuo en D y suma exactamente
    A. Con eso la barra honesta es 0.0019, el perfil es una parábola, y D pasó
    de 1.4436 a 1.4378. `por_bin=False` queda para compararse con los
    oráculos (Etapa 5 con `replica_paper`, Etapa 6 con `analisis/ncc.py`), y
    las dos verificaciones vuelven a pasar. Con la plantilla integrada, el B
    del NCC se mueve −0.63 % (Cs-137) y +0.34 % (Co-60).
  - **Datos reales, Klein-Nishina:** Cs-137 D = 1.4378 ± 0.0023 (σ_rel 0.120,
    χ²/ndf 1.40); Co-60 junio 0.9931 ± 0.0036 (0.127, 1.50); Co-60 binario
    1.3135 ± 0.0036 (0.130, 2.07). corr(D, σ) = −0.5 a −0.7. Los residuos
    tienen estructura (Klein-Nishina no tiene scattering múltiple ni el
    exceso).
  - **Verificación:** los 8 arranques pulidos caen en el mismo D a 1e-7 (sin
    pulir, dispersión 1.3e-3: es el último paso del recocido, no mínimos
    distintos); la barra por curvatura coincide con el perfil Δχ² = 1
    (0.00193 / 0.00207 en Cs-137, 0.00291 / 0.00291 en Co-60). **Cierre con
    verdad conocida**, 30 realizaciones Poisson del ajuste del Cs-137 (~304
    cuentas por bin) por la cadena entera: sesgo −0.034 ± 0.029 %,
    dispersión/barra 1.04, cobertura 22/30 = 73 %. A esta estadística el sesgo
    de Neyman no se ve.
  - **Contra `replica_paper/metodo.py` sobre nuestros espectros** (no es la
    misma cuenta: Ec. 14 sin pesos, A en el recocido, σ de la plantilla por
    barrido, barrido libre; en el Co-60, σ₁, σ₂ y A₂/A₁ libres): Cs-137
    1.4289 (−0.62 %, ventana [534, 1003]); Co-60 **1.0758 (+8.3 %)**, ventana
    [755, 1670]. Lo del Co-60 es la sub-determinación de §7.23 con 5
    parámetros; el nuestro fija A₂/A₁ por la física (intensidades) y usa un σ.
  - Costo: ~25 s por corrida (~145 000 evaluaciones del modelo); el cierre
    de 30 realizaciones, 11 min.
- 27-sep-2026 — **Etapa 8 cerrada** (`plantilla.py`, `verificar_plantilla.py`,
  `figuras/etapa8_plantilla.png`). `Plantilla(fuente)` lee `Edep_keV` del .root
  (350 bines), verifica con `E0_keV` que las líneas simuladas sean las de la
  fuente declarada (si no, `ValueError`) y entra por
  `modelo.modelo(..., plantilla=)`; NCC y ajuste la pasan sin otros cambios.
  La plantilla se integra por bin del espectro con su acumulada, igual que la
  Klein-Nishina desde la Etapa 7.
  - **Contra el oráculo** (`resultados.json`, el resultado del informe):
    Cs-137 1.3563 ± 0.0025 / 1.3563 ± 0.0022 (0.0 barras); Co-60 junio
    0.9551 ± 0.0042 / 0.9558 ± 0.0022 (0.1); Co-60 binario 1.2631 ± 0.0037 /
    1.2617 ± 0.0037 (0.3). Las ventanas coinciden a ±3 mV·ns. Nuestra barra
    del Co-60 de junio es el doble de la del oráculo: está escalada por
    √(χ²/ndf) = 1.34 y sale de un χ² continuo en D (ver Etapa 7).
  - **El cambio:** Klein-Nishina / Geant4 − 1 = +6.0 % (Cs-137), +4.0 % y
    +4.0 % (Co-60). §7.18 decía +4 %. Rechazo de un .root de otra fuente: ok.
  - **Cierre con la plantilla de Geant4** (8 realizaciones del Cs-137): sesgo
    −0.10 ± 0.05 %, dispersión/barra 0.95. Pasa (< 3σ), pero está a 2σ; con
    30 realizaciones se sabría si hay un sesgo de ~0.1 %.
  - **Hallazgo de esta etapa: w, la fracción de scattering simple, no la fija
    el dato y mueve D ~6 %.** `plantilla.Mezcla(fuente, w)` =
    w·(NCompton = 1) + (1 − w)·(resto). La simulación tiene w = 0.730 (Cs-137)
    y 0.751 (Co-60).
    - Con la ventana del NCC (desde 0.45 B), el dato "prefiere" w ≈ 0.9-1.1:
      χ² 400 → 263 (Cs-137), 369 → 279, 503 → 254, entre w de la simulación y
      w = 1. Eso NO es scattering: los residuos de Klein-Nishina y de Geant4
      tienen la misma estructura (positivos en el borde inferior de la
      ventana), que es el exceso de bajo canal que el modelo no tiene
      (`gui.EXCESO = False` en la campaña).
    - Con la ventana desde 0.80 B, lejos del exceso, el χ² ya no distingue w:
      Co-60 junio 160 / 159 / 160, Co-60 binario 122 / 122 / 122, Cs-137
      125 / 125 / 141 (w de la simulación / 0.85 / 1). Pero D se mueve entre
      esos w **+4.4 %, +4.6 % y +5.4 %**. O sea el dato no mide w; lo fija la
      simulación, y un error en la geometría simulada entra directo en D.
      Sospechoso concreto, sin probar: según §7.69 la simulación tira un haz
      FINO sobre el eje del cilindro, que recorre los 5 cm enteros y maximiza
      el scattering múltiple; una fuente puntual real ilumina toda la cara.
    - Aun con el w de la simulación, D **deriva** al subir el inicio de la
      ventana de 0.45 B a 0.80 B: +1.1 % (Cs-137), +1.9 % (Co-60 junio),
      +1.8 % (Co-60 binario). Es el patrón de §7.32: el exceso sin modelar
      arrastra D. §7.35 y §7.47 lo cerraban modelando el exceso; la campaña
      del informe corre sin él.
    - Consecuencia para el informe (no se tocó): sus barras de D (~0.2-0.3 %)
      son solo estadísticas; estas dos dependencias (w y la ventana) son
      sistemáticos de 1-6 % que no están en ellas. Conviene medirlos en el
      lote (Etapa 9) o modelar el exceso antes.
- 27-sep-2026 — **Etapa 9 cerrada** (`lote.py`, `verificar_lote.py`,
  `resultados_lote.json`, `figuras/etapa9_lote.png`). Las 28 corridas por la
  cadena de desde_cero (compuerta de 12 ns, Geant4 Edep) en ~23 min; el
  resultado de cada corrida se guarda apenas termina en
  `resultados_lote.json` (es el entregable, no un caché escondido;
  `--rehacer` lo regenera y `--tests` corre solo los tests).
  - **Las fórmulas, solas:** con las D del oráculo, `test_B` y `test_A`
    reproducen `tabla_testB.tex` y `tabla_testA.tex` en todas las cifras
    impresas (b, su error, χ² de potencia y de recta, peor residuo, D
    conjunta, Cs/Co, desacuerdo, c₀ y su error, misma/distintas sesiones).
    Hubo que alinear una definición: el peor residuo es |D − pred|/D, no /pred.
    q(Cs)/q(Co) = 0.9799 sale de la simulación con la misma receta.
  - **Las D:** 24 de 28 dentro de 1 barra del oráculo (media −0.27 barras).
    Las otras 4 son TODAS las Cs-137 del detector D y todas salen bajas
    (−0.36 a −0.66 %). **Causa demostrada:** la espiga del borde del DRS4
    (Etapa 2). Recalibrando `d_cs137_1600` con el recorte de 2 muestras del
    oráculo da 0.4934, exactamente la suya. El oráculo pierde ahí **2.1 % de
    los pulsos** (1647 de 79672), todos chicos: el detector D es el de menor
    ganancia y tiene muchos pulsos por debajo de ~30 mV. **Es un error de
    `analisis/` que está en el informe:** D del detector D con Cs-137 sale
    +0.4-0.5 % alta (1.3-1.5 barras en las corridas buenas), y el desacuerdo
    del test A en D @ 1600 / 1650 V baja de +19.3 / +16.5 % a +18.5 / +15.9 %.
    Ninguna conclusión cambia. `analisis/` no se tocó.
  - **Test B con nuestras D:** A Co-60 9.45, A Cs-137 10.42, B 9.42 / 9.56,
    C 9.01 / 9.36, D 10.19 / 7.79. Todos a ≤ 0.4 barras del informe con la
    barra escalada por √(χ²/gdl). **Sin escalar, el b de C Co-60 difería 3.6
    barras** (9.01 contra 8.70): con χ² = 71-168 y 1 gdl la ley de potencia
    no describe esos tres puntos, b depende de los pesos relativos de las D,
    y la barra sin escalar de la tabla del informe (±0.05) es unas 12 veces
    corta. Lo mismo en D Co-60 (χ² 149-493, ±0.04 contra ±0.9 escalada).
    §7.73 ya decía que el err_b del informe no está escalado.
  - **Test A con nuestras D:** el mismo patrón que el informe, punto por
    punto dentro del error. A y B +9 a +15 % (A crece con el HV, B plano),
    C −8 a −11 %, D +18.5 / +15.9 / +2.0 % (misma sesión). Sigue sin cerrar:
    hace falta medir las dos fuentes en la misma sesión y el mismo día.
  - Lo que esta etapa NO midió, y queda de la Etapa 8: la dependencia de D
    con w (hasta ~6 %) y con el inicio de la ventana (1-2 %) no entra en
    ninguna barra de estas tablas.
- 27-sep-2026 — **Etapa 10 cerrada** (`clasicos.py`, `verificar_clasicos.py`,
  `figuras/etapa10_clasicos.png`). Cuatro métodos, porque el "Gaussian FWHM
  method" tiene dos lecturas en el proyecto: 50 % de altura, mínimo de la
  derivada (con la ventana de suavizado sacada del ancho del flanco, 0.75 ×
  caída 80 → 20 %), **escalón gaussiano** (erfc + fondo lineal, la lectura de
  `codigo/etapas`) y **gaussiana del paper** (gaussiana al pico, borde en
  μ + σ√(2 ln 2), la lectura literal de `replica_paper`). Región anclada al B
  del NCC ± 45 %.
  - **Contra los oráculos:** los cuatro, iguales a `codigo/etapas/clasicos.py`
    y `replica_paper/convencionales.py` a redondeo de máquina, en dos
    sintéticos y en las tres corridas de trabajo.
  - **Sesgo con verdad conocida** (Klein-Nishina de Cs-137, B conocido, 10
    realizaciones, σ/B = 0.03 / 0.06 / 0.09 / 0.12 / 0.18): 50 % +0.8 / +2.3 /
    +3.6 / +5.1 / +7.7 %; diferencia +0.1 / +1.6 / +1.3 / +2.6 / +4.7 %;
    escalón +2.2 / +4.0 / +5.4 / +6.5 / +7.6 %; gaussiana del paper +2.2 /
    +2.4 / +2.9 / +3.6 / +5.4 %. Todos sobrestiman y crecen con σ, del orden
    de §7.44 (la diferencia sale 0.6-2.7 puntos más baja que aquella tabla,
    que se midió sobre otro sintético: B = 600 canales y otro modelo
    directo). Sin ruido da casi lo mismo: el sesgo es del método.
  - **Sobre los datos, contra la D del ajuste con Geant4:** Cs-137 +10.1 a
    +12.7 % según el método (el sintético de Klein-Nishina a esa resolución,
    σ/B ≈ 0.11, da 2-6 %: el resto es la forma real, scattering múltiple incluido).
  - **Encontrado: en el Co-60 la comparación no está bien definida, y el
    informe la hace igual que la hacía yo al principio.** Los clásicos
    detectan la caída del conjunto de los dos bordes solapados, que la domina
    el SEGUNDO (1118 keV). Dividiendo por el PRIMERO (963 keV, `E_ref`, que es
    lo que hace `informe.clasicos_reales`) el sesgo sale **+19 a +23 %**, y
    un +16 % de eso es solo 1118/963. Contra el segundo borde da −3.6 a
    +6.2 %. Ninguna de las dos es "la" respuesta: un método que no sabe
    cuántos bordes hay no dice cuál midió. La mediana de +11.8 % que el
    informe da sobre 25 corridas mezcla Cs-137 (bien definida) con Co-60 (no),
    y conviene separarla por fuente. (El número exacto del informe para
    el Co-60 no se reprodujo: ancla la región en el B de su ajuste, no en el
    del NCC.)
  - **Robustez (Fig. 5), seis intervalos anidados:** nuestro ajuste con Geant4
    da el borde dentro de 1.1 % (Cs-137) y 1.0 % (Co-60) de un intervalo a
    otro. Los clásicos saltan: diferencia 14 % y 24 %, escalón 7.5 % y 21 %,
    en los intervalos anchos porque entran al exceso de bajo canal. 50 % y
    la gaussiana del paper solo funcionan en 2-3 de los intervalos (en los
    angostos no hay máximo adentro). Es el resultado cualitativo que el paper
    reporta para NCC-SA, ahora con esta cadena.

  **Con esto el plan de `desde_cero/` está completo (Etapas 0-10).** Lo que
  queda abierto, por orden de impacto en D, está en las bitácoras de las
  Etapas 8 y 9: (1) la fracción de scattering simple w, que el dato no fija y
  mueve D ~5 %; (2) la deriva de D con el inicio de la ventana (1-2 %) por el
  exceso de bajo canal sin modelar; (3) el error de la espiga del DRS4 en
  `analisis/` (+0.5 % en D del detector D con Cs-137); (4) la barra de b sin
  escalar en la tabla del test B; (5) el sesgo de los clásicos en el Co-60.
