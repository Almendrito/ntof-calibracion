# Campaña para cerrar la calibración y los dos tests de linealidad

Escrito el 28-ago-2026. Es la lista de lo que hay que MEDIR, en orden de
prioridad, y por qué cada corrida hace falta. Lo que se puede hacer sin medir ya
está hecho y vive en `analisis/`.

---

## De dónde se parte

| corrida | HV | fuente | D [mV·ns/keVee] | σ/B | χ²/ndf |
|---|---|---|---|---|---|
| junio 2026 | 1550 V | Co-60 | 1.2512 ± 0.0047 | 0.088 | 1.17 |
| 26-ago-2026 | 1600 V | Cs-137 | 1.9450 ± 0.0038 | 0.104 | 0.97 |

Dos calibraciones, cada una impecable por separado, y **ningún test posible**:
están a voltajes distintos y son de fuentes distintas, así que no hay ni dos
fuentes al mismo voltaje (test A) ni dos voltajes con la misma fuente (test B).

Y hay una señal de alarma. Si entre junio y agosto lo único que hubiera cambiado
fuera el HV, el cociente de las dos D daría el exponente de la ley de ganancia:

```
b = ln(1.9450 / 1.2512) / ln(1600 / 1550) = 13.9
```

Un PMT tiene b entre 6 y 10. **13.9 está fuera de rango**, así que entre las dos
corridas cambió algo más que 50 V: la distancia de la fuente, el tubo, o la
temperatura. **El atenuador queda descartado**: son 20 dB en las 25 corridas
(confirmado por el operador, 22-sep-2026; CLAUDE.md §7.70), y el test B lo había
anticipado —el punto de junio cae sobre la curva `D ∝ V^9.49` a −0.08%, cosa
imposible si la atenuación hubiera cambiado, porque el salto sería de un factor
10—. Es la razón principal por la que M1 encabeza la lista.

---

## Las corridas, en orden

### M1 — Cs-137 y Co-60 a 1550 V, en la misma sesión  ← el bloqueo

Con las dos fuentes al mismo voltaje y en la misma sesión aparece **el primer
test A real**: dos bordes independientes (477.3 y 963.4 keVee) sobre la misma
recta, y una ordenada al origen que se puede comparar con cero.

Repetir el Co-60 aunque ya exista el de junio: si el detector cambió, el de
junio no sirve, y no hay forma de saberlo sin volver a medirlo.

- fondo propio para cada una;
- misma distancia fuente-detector, anotada;
- ~20-30k eventos útiles alcanzan (ver "cuánta estadística").

### M2 — Barrido de HV con Cs-137: 1450, 1500, 1550, 1600, 1650, 1700 V

Es el test B. Los tres voltajes que ya existen (1600/1650/1700) **no alcanzan**,
y eso está calculado, no es una impresión: la separación entre una recta y la
ley de potencia, en el punto medio del rango, vale `(1/8)·b(b−1)·(ΔV/V̄)²`, o sea

| rango | span | separación (b = 8) | ¿decide? |
|---|---|---|---|
| 1600-1700 | 6.1% | 2.6% | no, es del tamaño del sistemático |
| 1500-1700 | 12.5% | 11% | sí |
| 1450-1700 | 15.9% | 18% | sí, con holgura |

`analisis/voltaje.py` lo verifica con datos sintéticos: con seis voltajes de
1450 a 1700 recupera b = 7.94 ± 0.07 sobre un valor verdadero de 8.00 y rechaza
la recta (χ² 10.6 contra 327); con sólo los tres de 1600 a 1700 el mismo código
avisa que **el span no alcanza para decidir**.

Detalles que cambian el resultado:

- **fondo propio en cada voltaje.** El de 1600 V no sirve a 1650: es otra
  ganancia.
- **umbral de trigger escalado con la ganancia** (`V_thr ∝ D(V)`). Si el umbral
  queda fijo en mV, el corte en energía se corre entre corridas y cada una queda
  contaminada distinto por abajo. El ajuste tiene un piso duro automático en el
  fin de la subida del trigger (`espectro.fin_del_umbral`), así que un umbral mal
  puesto no rompe el ajuste, pero se come parte del flanco.
- **repetir el primer voltaje al final de la sesión.** Un PMT deriva después de
  cambiar el HV: dejar 10-15 min de estabilización en cada punto y cerrar
  repitiendo el primero mide esa deriva en vez de suponerla.
- **verificar saturación en el punto más alto.** El DRS4 satura en ±500 mV. A
  1600 V hay 3 eventos saturados de 80322 (0.004%), así que sobra margen; a
  1700 V la ganancia es ~1.4 veces mayor y sigue sobrando, pero hay que mirar la
  bandera `MAL_SATURA` y no confiar en la extrapolación.

### M3 — Fuentes extra al voltaje de referencia

Con dos fuentes el test A tiene **un solo grado de libertad**: la recta con
ordenada libre pasa exacto por los dos puntos y lo único que se testea es si esa
ordenada es compatible con cero. Con cuatro fuentes hay residuos de verdad.

El paper de Kong et al. usa exactamente **Cs-137, Co-60, Mn-54 y Na-22**, que
dan seis bordes bien repartidos:

| fuente | línea γ [keV] | borde Compton [keV] |
|---|---|---|
| Na-22 | 511.0 | 340.7 |
| Cs-137 | 661.7 | 477.3 |
| Mn-54 | 834.8 | 639.5 |
| Na-22 | 1274.5 | 1061.7 |
| Co-60 | 1173.2 | 963.4 |
| Co-60 | 1332.5 | 1118.1 |

Mn-54 es la más limpia (una sola línea). Cada fuente nueva necesita **su
simulación Geant4 con esta geometría**: ya hay `.root` de Co-60 (3M) y de Cs-137
(10M), y `simulacion.coincide_fuente` verifica que la plantilla sea de la fuente
elegida antes de dejar calibrar. Sin simulación queda la ruta Klein-Nishina, con
su sesgo conocido de +4.3% por no tener scattering múltiple.

### M4 — Repeticiones en el punto de referencia (3 a 5)

El paper hace **diez mediciones independientes por par (centellador, fuente)** y
usa la dispersión entre ellas como medida de estabilidad. Acá la estadística del
ajuste da ±0.3%, y el sistemático es ~2-3%: la única forma de saber cuál es la
barra de error real de un punto es repetir la medición. Esa dispersión, y no
`err_D`, es lo que hay que usar en los residuos del test A.

### M5 — Qué es el exceso de bajo canal

Sigue siendo el sistemático más grande y ahora limita directamente los residuos
del test A. Dos corridas cortas al mismo HV lo resuelven:

- fuente al **doble de distancia** (si el exceso escala con la tasa, es
  dispersión en el material de alrededor);
- fuente **apantallada con plomo** (si no se mueve, es del fotodetector).

Hoy el exceso se modela con una exponencial de λ ≈ 380 mV·ns (Co-60) y ≈ 224
(Cs-137); cambiarla por una ley de potencia mueve D un 1.1% en el Co-60 y un
0.5% en el Cs-137. Eso es el sistemático de la forma, y no se cierra ajustando
mejor.

### M6 — Bitácora por corrida

El campo que faltaba y que hizo imposible comparar junio con agosto. Anotar, por
corrida:

```
fuente · HV [V] · atenuador [dB] · umbral [mV] · muestreo [GS/s]
distancia fuente-detector [mm] · temperatura [°C] · hora de inicio
número de eventos pedido · nombre del archivo
```

El catálogo ya tiene estos campos (`rutas.Corrida`) y `rutas.misma_cadena` avisa
cuando dos corridas que se van a combinar difieren en la cadena de señal, o
cuando el dato **no está registrado** (que no es lo mismo que ser igual).

✔ El atenuador ya está: **20 dB en las 25 corridas** (confirmado por el operador,
22-sep-2026). Sigue en la lista porque es el campo que más duele si alguna vez
cambia sin anotarse —el área se integra a la entrada del digitalizador, así que
20 dB dividen D por 10— y porque el valor hay que seguir anotándolo corrida a
corrida, no darlo por heredado. `rutas.misma_cadena` avisa si aparece uno
distinto.

---

## Cuánta estadística

Con 80k eventos el error estadístico de D es 0.25% y el sistemático 2-3%: **más
eventos no sirven de nada**. Unos 20-30k eventos útiles por punto alcanzan. El
tiempo de laboratorio rinde mucho más en más voltajes, más fuentes y más
repeticiones.

Ojo con el tiempo vivo: el DRSOSC escribe ~35 kB de XML por evento y no da
abasto. El Co-60 estuvo vivo el 40% y el Cs-137 de agosto el **0.13%** (46 kHz de
trigger real contra 58 eventos/s guardados). Subir el umbral de trigger es lo
que más rinde.

---

## Qué hace el análisis con cada cosa

| lo que se mide | qué habilita | dónde se procesa |
|---|---|---|
| M1 | test A, 2 fuentes, 1 gdl | `registro.linealidad`, `por_punto` |
| M2 | test B, b y el rechazo de la recta | `voltaje.ajustar`, `voltaje.colapso` |
| M3 | test A con residuos reales | `registro.linealidad` |
| M4 | la barra de error honesta de cada punto | `voltaje.puntos_del_registro` promedia |
| M5 | baja el sistemático dominante | `modelo.forma_exceso` (exp ↔ pot) |
| M6 | que las corridas sean comparables | `rutas.Corrida`, `rutas.misma_cadena` |

Con M1 y M2 hechas, `python analisis/registro.py --rehacer` y
`python analisis/voltaje.py` dan los dos tests sin tocar una línea de código.

---

## Criterios, fijados antes de mirar los resultados

- **Test A pasa** si la ordenada libre es compatible con cero dentro de 2σ *con
  el error total*, no con el estadístico, y los residuos de la recta por el
  origen quedan por debajo del sistemático punto a punto.
- **Test B** se reporta como el valor de **b** y el χ² de las dos hipótesis. Si
  el span es menor al 8%, la recta no se puede rechazar y hay que escribirlo.
- **La métrica del paper no sirve acá.** Kong et al. validan con *cosine
  similarity* 99.999% entre dos relaciones energía-canal. Entre vectores casi
  proporcionales esa métrica es insensible por construcción: dos calibraciones
  que difieren 5% en pendiente siguen dando ~99.99%. Reportar **residuos en keV**
  y la **ordenada con su error**.
