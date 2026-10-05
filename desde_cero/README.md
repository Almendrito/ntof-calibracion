# desde_cero — la cadena de análisis, construida por etapas

Reimplementación desde cero, chica y explicada, de lo que hace `analisis/`.
El objetivo no es tener otro programa que funcione: es **entender por qué cada
paso es como es**. `analisis/` se queda donde está y sirve de oráculo.

## Reglas

1. **Una etapa = un archivo = una idea.** Cada archivo corre solo
   (`python desde_cero/<archivo>.py`) y muestra lo que hace.
2. **Nada entra sin comprobarse.** Contra `analisis/` mientras exista, y contra
   datos fabricados con verdad conocida cuando haga falta.
3. **Lo mínimo que funcione.** Sin cachés, sin barras de progreso, sin interfaz,
   hasta que algo lo pida de verdad.
4. **El porqué va en el código**, al lado de la decisión que explica.

## Etapas

El plan detallado, con el prompt de Claude Code de cada etapa, está en
`PLAN.md`.

| # | archivo | qué resuelve | estado |
|---|---------|--------------|--------|
| 0 | `lector.py` | archivo del DRS4 → `(t, v)` por evento, en los dos formatos | ✅ |
| 1 | `pulso.py` | mirar un pulso: línea base, onset, amplitud | ✅ |
| 2 | `area.py` | el observable: integrar el pulso con compuerta → mV·ns | ✅ |
| 3 | `espectro.py` | histograma, tiempo vivo, resta de fondo, barra de error | ✅ |
| 4 | `fisica.py` | borde Compton y forma Klein-Nishina (Ec. 1, 2, 10) | ✅ |
| 5 | `modelo.py` | forma → canal, convolución (Ec. 8, 9) | ✅ |
| 6 | `ncc.py` | la ventana automática (Ec. 15) | ✅ |
| 7 | `ajuste.py` | recocido simulado, χ², barra de error → **D** | ✅ |
| 8 | `plantilla.py` | de Klein-Nishina a Geant4 | ✅ |
| 9 | `lote.py` | las 28 corridas, ley de ganancia y acuerdo entre fuentes | ✅ |
| 10 | `clasicos.py` | (opcional) los métodos convencionales y la robustez | ✅ |

Las etapas 1-3 son sobre el DATO y ahí estuvo el error más grande del proyecto
(el área mal integrada, `CLAUDE.md` §7.33). Las 4-7 son el método del paper;
la 8 y la 9, lo que este proyecto le agregó.

**Fuera del orden:** `ver_sa.py` mira cómo `analisis/` usa el recocido simulado,
sin reimplementarlo todavía. Es para entender el optimizador antes de que la
etapa 5 lo pida.

## Etapa 0 — leer los datos

```bash
python desde_cero/lector.py      # metadatos + figura con pulsos de los dos formatos
python desde_cero/verificar.py   # que dé lo MISMO que el lector de analisis/
```

La única función que hace falta desde afuera:

```python
from desde_cero.lector import eventos

for t, v in eventos(ruta, limite=200):   # t en ns, v en mV
    ...

# desde la Etapa 3, también el reloj y el scaler de la placa (tiempo vivo)
for t, v, reloj, scaler in eventos(ruta, con_reloj=True):
    ...
```

Cuatro decisiones, y por qué:

- **El formato se mira en los primeros 4 bytes, no en el nombre.** Los binarios
  del 1-sep no tienen extensión (`Co60`, `Cs137`, `B1700V`).
- **`eventos` es un generador.** El XML de Co-60 pesa 1.76 GB; leer 200 eventos
  cuesta 0.8 s y unos MB de RAM porque se corta apenas alcanza.
- **El canal se deduce del archivo.** Junio salió por `CHN2`, agosto por `CHN1`.
  Fijarlo a mano y errarle no da error: da cero pulsos, y el problema aparece
  cinco pasos después como un histograma vacío.
- **El eje temporal se reconstruye por evento** en el binario. El DRS4 es un
  anillo de 1024 condensadores: cada celda tiene su propio ancho (~5% de
  dispersión, vienen medidos en la cabecera) y el anillo se congela en una celda
  distinta cada vez (`T#`), así que el orden en que se suman los anchos cambia
  evento a evento. Se ve en la salida: `dt_min` ≠ `dt_max`.

Lo que ya encontró la comprobación:

- **Un error real, mío:** leía el "rango" del evento en el byte 14, que es la
  *hora*. También es un número chico y plausible, así que nada fallaba: los
  pulsos salían corridos ~16 mV. Lo cazó `verificar.py`.
- **Una rareza del archivo, no un error:** la celda 498 del Co-60 a 1600 V tiene
  ancho −0.0022 ns. Es 1 de 1024 y el 0.001% del registro; viene así de la
  calibración temporal del DRS4.
