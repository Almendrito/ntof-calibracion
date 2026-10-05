# El informe de calibración

`informe.pdf` sale de compilar `informe.tex`. Lo que hay acá se reparte en dos
grupos, y la diferencia importa:

| archivo | quién lo escribe |
|---|---|
| `informe.tex` | **vos**, a mano. Es el texto del informe. |
| `preambulo.tex` | **vos**, a mano. Paquetes, estilo, atajos, y las etiquetas de procedencia `\enpaper`, `\adaptado`, `\propio`. |
| `flujo.tex` | **vos**, a mano. El diagrama TikZ de las dos cadenas (paper y nuestra); es texto, no datos. |
| `informe_v1_2026-09-24.*` | respaldo del informe antes de la reestructuración del 24-sep-2026. |
| `generado/*.tex` | `python analisis/informe.py`. **Se reescribe entero** en cada corrida. |
| `figuras/*.pdf` | `python analisis/informe.py`. Idem. |
| `PLAN.md` | **vos**, a mano. Qué le falta al informe para ser el texto guía, bloque por bloque y sesión por sesión. Se lee antes de escribir. |

O sea que se puede editar el texto todo lo que haga falta y volver a correr el
script sin perder nada.

## Regenerar

```bash
python analisis/calibrar_todo.py              # calibra las corridas (si cambiaron los datos)
python analisis/calibrar_todo.py --resumen --control
python analisis/informe.py                    # figuras + tablas + numeros + PDF
python analisis/informe.py --sin-pdf          # sin compilar, para iterar rapido
```

A mano, si ya están las tablas: `pdflatex informe.tex`, `bibtex informe`,
y `pdflatex informe.tex` dos veces más (bibtex arma la bibliografía desde
`referencias.bib`; las dos pasadas siguientes resuelven las citas, el índice
y las referencias cruzadas). El script hace las cuatro.

Las secciones eficaces del marco teórico salen de `datos/xcom/xcom_*.txt`,
tablas descargadas de NIST XCOM con su procedencia (URL, fórmula, fecha) en
la cabecera. No se regeneran: son datos de referencia. Para otro material se
descarga otra tabla con el formulario "Compound" de XCOM y se guarda con la
misma cabecera; `informe.xcom(nombre)` la lee.

La sección «La simulación de Geant4» lee la geometría, el material, la
constante de Birks y el corte de la FUENTE de la simulación, que está copiada
en `simulacion/` (`informe.geant4_geometria`; el original vive en WSL, ver
`simulacion/LEEME.md`), y las estadísticas de los `.root` de `datos/`
(`informe.simulacion_resumen`, cacheado por tamaño y fecha del archivo). El
efecto del binning de la plantilla se MIDE reajustando la corrida de ejemplo
con 100 a 700 bines (`informe.binning_plantilla`, ~10 s, cacheado con la D
guardada como firma). Ojo: en `JcV/ComptonSimulation/` hay otra simulación
más vieja que NO es la de los `.root`; el informe no la lee.

`referencias.bib` se edita a mano y tiene su regla: **ninguna entrada se
agrega sin verificar** la revista, el volumen y las páginas. Se cita con
`\cite{clave}`.

La primera corrida después de cambiar el ajuste de la corrida de ejemplo
(detector A a 1600 V con Cs-137) tarda unos 4 minutos más: es el cierre con
verdad conocida (30 realizaciones Poisson por la cadena completa, dos veces
cada una, pesando con el dato y con el modelo). Queda en
`.cache/informe_cierre_<clave>.json` con la verdad como firma, así que sólo
se rehace si cambia el ajuste.

## Los números no se tipean

Cada número que el texto cita está en `generado/numeros.tex` como una macro, y
en el `.tex` se escribe la macro. Por ejemplo `\NCorridas` en vez de `25`, o
`\CocEsperado` en vez de `0.980`. La lista completa está al principio de ese
archivo. Si al escribir un párrafo hace falta un número que todavía no es
macro, conviene agregarlo en `analisis/informe.py` (función `numeros`) en vez
de escribirlo en el párrafo: así el informe no puede quedar desincronizado de
los datos.

La sección «Del pulso al espectro» lee además los pulsos de cada corrida
(`informe.pulsos_de`): mide la ventana sobre los primeros ~600 registros y
saca las áreas del `.npz` que la calibración ya dejó, con la misma
compuerta, así que no re-parsea nada; si el caché no estuviera, lo
regeneraría (80 s por XML). Sus números son macros como los demás.

Lo mismo con las tablas: salen de `analisis/informe.py` y se traen con
`\input{generado/tabla_...}`. El epígrafe y la posición sí viven en el `.tex`,
porque son decisiones del documento.

## Qué hay en cada figura

| archivo | qué muestra |
|---|---|
| `integracion_pulso.pdf` | del registro al área: pulsos crudos, un pulso con cada paso de la integración, el pulso promedio con la compuerta y la carga acumulada |
| `pulsos_detectores.pdf` | el pulso promedio de los cuatro detectores a 1600 V y su carga acumulada, normalizada a la carga a 40 ns |
| `area_a_espectro.pdf` | área contra amplitud (la forma no depende del tamaño) y los tres histogramas: muestra, fondo normalizado y neto |
| `teoria_ncc.pdf` | Klein-Nishina, la misma borroneada con la resolución medida, y el barrido del NCC |
| `teoria_sa.pdf` | los tres ingredientes del recocido simulado, con una función de juguete |
| `teoria_interaccion.pdf` | coeficientes másicos por proceso (fotoeléctrico, Compton, pares, Rayleigh) en BC-404 y en NaI, de NIST XCOM (`datos/xcom/`), con Klein-Nishina de electrón libre superpuesta |
| `teoria_cinematica.pdf` | energía del electrón contra ángulo para las tres líneas, Klein-Nishina en ángulo contra Thomson, y la sección eficaz total contra energía (analítica, numérica y XCOM) |
| `teoria_clasicos.pdf` | un borde sintético con la resolución de la corrida de ejemplo y dónde lo pone cada método clásico; el sesgo de los tres contra σ/B, con y sin ruido |
| `teoria_ncc_detalle.pdf` | el NCC por dentro: las dos curvas estandarizadas en la ventana ganadora, una contra otra, y ρ y B contra el número de puntos de remuestreo |
| `simulacion_geant4.pdf` | qué hay en la simulación: Edep contra Lvis del Co-60 con los bordes en las dos escalas, el Cs-137 separado por número de scatterings, el quenching de Birks evento a evento, y la razón dato/modelo del ajuste de ejemplo en todo el rango (el exceso de bajo canal y el umbral, que la simulación no tiene) |
| `teoria_estadistica.pdf` | el perfil de Δχ² contra la parábola de la curvatura, el plano (D, σ) con las dos elipses y el cierre, y la D recuperada realización por realización |
| `ajuste_A/B/C/D.pdf` | un ajuste por detector, a 1600 V: bordes en canal y espectro calibrado |
| `cociente_fuentes.pdf` | D(Cs-137)/D(Co-60) por punto (a 0.72 del ancho del texto) |
| `comparacion_por_punto.pdf` | los dos espectros de cada punto de operación, con una calibración común |
| `ganancia_dv.pdf` | D(V) con la ley de potencia y sus residuos |
| `comparacion_detectores.pdf` | ganancia, ganancia relativa, resolución y calidad de los ajustes |
| `resumen_ganancia.png`, `resumen_linealidad.png` | las figuras-resumen de `calibrar_todo`, tal cual; no entran en el PDF, el texto las cita por nombre |

## Tamaño y letra de las figuras

Las figuras se generan **del ancho que van a tener en la página**, con la
letra del documento y en PDF. Es lo que hace que no se vean chicas: una
figura de 15 pulgadas reducida a `\linewidth` (6.38 in en este documento)
pierde el 57% de su letra, y ninguna fuente grande lo compensa. La regla es
la inversa de lo intuitivo: figura chica, letra normal, y `\includegraphics`
sin escalar.

Los números viven en `analisis/informe.py` (sección 3) y salen de
`informe.tex` y `preambulo.tex`:

| qué | valor | de dónde |
|---|---|---|
| ancho del texto | (21.0 − 2·2.4) cm = 6.38 in | `a4paper`, `margin=2.4cm` |
| rótulos de ejes | 10 pt | el `\small` de los epígrafes, cuerpo 11 pt |
| ticks y leyendas | 9 pt | `\footnotesize` |
| tipografía | Latin Modern Roman | la misma `lmodern` del `.tex`, sacada del MiKTeX con `kpsewhich` |

Reglas que siguen de eso:

- **El `width=` del `\includegraphics` tiene que ser el mismo con que se
  generó la figura** (`figura(ancho=...)`): `\linewidth` para casi todas,
  `0.72\linewidth` para `cociente_fuentes`. Si se cambia uno hay que cambiar
  el otro, o la letra deja de medir lo que dice.
- Si cambian los márgenes o el cuerpo en `preambulo.tex`/`informe.tex`, hay
  que cambiar `MARGEN_CM`/`CUERPO_PT` en `informe.py` y regenerar.
- Nada de `bbox_inches="tight"`: recorta el margen, la figura queda más
  angosta que `\textwidth` y LaTeX la vuelve a estirar.
- Los títulos de panel no van dentro de la figura (el título es el
  epígrafe); los paneles llevan `(a)`, `(b)`, ... cuando el epígrafe los
  nombra así.
