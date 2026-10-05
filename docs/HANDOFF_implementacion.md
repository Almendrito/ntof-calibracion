# HANDOFF a Claude Code — Implementación etapas 4–7 (calibración Compton, método NCC-SA)

> Plan de implementación + **decisiones de método ya tomadas**. Pareja con `CLAUDE.md`
> (que tiene el contexto general, el método del paper y los hallazgos §7). Este archivo
> es el "qué construir y cómo decidir" para codificar directo.

---

## 0. Reglas de trabajo (mantener el estilo del proyecto)

- El usuario está **aprendiendo mientras construye**: explicar el PORQUÉ de cada patrón, no solo entregar código.
- **Por etapas**, cada una ejecutable y verificable sola.
- **Validar contra datos SINTÉTICOS con verdad conocida ANTES de correr sobre el Co-60 real.**
- Español. Comentarios en español. Anclar en el paper; marcar lo **[fuera del paper]**.

---

## 1. Estado actual

Hechas (funcionando): `spectrum.py` (Etapa 0), `physics.py` (Etapa 1), `model.py` (Etapa 2),
`fit_single.py` (Etapa 3). Datos: **solo Co-60 + fondo** (2 bordes solapados: 963.4 y 1118.1 keV).
`fit_single` usa `least_squares` local → **frágil** (cae en mínimos locales).

Faltan: **Etapa 4 (NCC)**, **Etapa 5 (SA)**, **Etapa 6 (dos bordes)**, **Etapa 7 (recta)**.

---

## 2. DECISIONES DE MÉTODO (ya tomadas — con justificación)

**D1 · Objetivo que minimiza el SA → mínimos cuadrados (Ec. 14), NO el error NCC.**
El paper es ambiguo (Ec. 14 dice LS; el texto de Metropolis dice "NCC-based error"). Decisión:
el SA minimiza **LS**. Motivo: ya verificamos (mapa de la loss, Etapa 3) que el mínimo GLOBAL
de la loss LS está en el borde real; el problema NUNCA fue el objetivo, sino el optimizador local.
El **NCC se usa SOLO para elegir el intervalo** (Etapa 4), no como objetivo del ajuste.

**D2 · Suavizado → Savitzky-Golay SOLO para el NCC, nunca para el ajuste final.**
El NCC es sensible a máximos espurios por ruido de Poisson (pocas cuentas). Suavizar los datos
antes de correlacionar. Pero el **ajuste (SA) corre sobre el neto CRUDO** (sin suavizar ni recortar),
para no destruir estadística ni sesgar. Semilla: ventana impar ~11–21, orden 3; ajustar mirando datos.

**D3 · Modelo del Co-60 → SIEMPRE dos bordes (Etapa 6).**
El ajuste de un borde (`fit_single`) sobre Co-60 da un punto SESGADO (cae entre los dos bordes,
σ inflado). Solo sirvió como bootstrap pedagógico. El resultado real usa el modelo compuesto.

**D4 · ΔE en la Ec. 12 → usar energías de BORDE (E_max2 − E_max1), NO de líneas gamma.**
Corrección al paper (análoga al erratum de la Ec. 1). El canal escala con la energía del
ELECTRÓN de retroceso (el borde), no con la del gamma. Entonces:
`CE_ch2 = B1 + D·(E_max2 − E_max1)`, con E_max1=963.4, E_max2=1118.1 keV → ΔE_borde=154.7 keV.

**D5 · Semilla de D → `D0 = B1 / E_max1`  [mV·ns / keV]  (pendiente, NO 1/E_max).**
D es "canales por keV" (la pendiente de calibración). Con el par aproximado (B1 ↔ 963.4 keV) y
origen asumido en 0, la semilla es B1/E_max1. Da al SA la escala correcta de arranque.

**D6 · Recta de calibración (Etapa 7) → con solo Co-60 son 2 puntos.**
Dos puntos determinan una recta EXACTA pero **no permiten validar linealidad** ni anclar energías
bajas (ambos puntos están juntos y altos, 963–1118 keV). Reportar pendiente+offset y ADVERTIR la
limitación. **Recomendación fuerte:** conseguir un Cs-137 (borde 477 keV) para un 3er punto y validar.
Decidir origen-vs-offset recién con los puntos reales en mano.

---

## 3. SPEC por etapa (qué codificar)

Reutiliza SIEMPRE `physics.py`, `model.py`, `spectrum.py`. No dupliques la lineshape ni el kernel.

### Etapa 4 — NCC · nuevo módulo `ncc.py`
Reemplaza el `intervalo=(lo,hi)` manual de `calibrar.py`.

Contrato sugerido:
```python
def seleccionar_intervalo(canales, cuentas, E_gamma,
                          ancho_ventana=..., suavizar=True, n_puntos=100) -> (lo, hi)
```
Algoritmo (Sección II-B):
1. Plantilla S2 = `dsigma_dE` muestreada en `n_puntos` (100) desde el plateau hasta el mínimo del borde.
2. (Opcional) pre-emborronar la plantilla con un σ nominal para que calce mejor con datos ya borrosos.
3. Suavizar `cuentas` (Savitzky-Golay, D2) → S1.
4. Para cada canal de inicio no-nulo: tomar una ventana de ancho `ancho_ventana`, **reinterpolar a 100
   puntos**, calcular Pearson (Ec. 15) contra S2, guardar ρ.
5. Devolver el intervalo (lo,hi) del ρ máximo.

Decisiones internas a probar: plantilla ideal vs pre-emborronada (empezar ideal); `ancho_ventana`
(tunable, ~pocos σ del borde). **Validar:** sintético con borde en B conocido → el intervalo debe
centrarse en B. La reinterpolación a 100 pts es OBLIGATORIA (invariancia al ancho del candidato).

### Etapa 5 — Simulated Annealing · nuevo módulo `optimize.py`
Mata la fragilidad del `least_squares`.

Contrato sugerido:
```python
def ajustar_sa(canales, cuentas, modelo_fn, p0, delta, bounds,
               T_init=2000, alpha=0.95, n_iter=100,
               T_min=1e-6, tol=1e-3, semilla=None) -> (params, loss)
```
Detalles (del paper, copiar tal cual):
- Objetivo: LS (Ec. 14) sobre el intervalo (D1). `E(x) = Σ (cuentas − modelo_fn(x))²`.
- Propuesta: `x_new = x + Δx`, con `Δx_i ~ U(-0.1·δ_i, 0.1·δ_i)`.
- Metropolis: aceptar con `P = exp(-ΔE / T_k)` (mejora → seguro; empeora → a veces).
- Enfriar geométrico: `T ← alpha·T` tras las `n_iter` iteraciones.
- Parar: `T < T_min` **o** `|ΔE_mean| < tol` durante 10 ciclos.
- `modelo_fn` **pluggable** (un borde o dos). Sembrar `p0` desde el NCC (edge ≈ centro del intervalo)
  y D con D5. Respetar `bounds` (A>0, σ>0, edge dentro del intervalo).

**Validar:** correr sobre el MISMO caso donde `least_squares` se atascaba (Etapa 3, B=5.0 desde un
arranque malo) → el SA debe recuperar el B real. Ese es el criterio de éxito.

### Etapa 6 — Modelo de dos bordes · extender `model.py`
Contrato sugerido:
```python
def forward_model_doble(canales, E_g1, E_g2, A1, A2, B1, sigma1, sigma2, D) -> S_comp
```
Lógica (Ecs. 11–13, con la corrección D4):
- Borde 1: `forward_model(canales, E_g1, A1, B1, sigma1)`.
- `E_max1 = compton_edge_energy(E_g1)`, `E_max2 = compton_edge_energy(E_g2)`.
- `CE_ch2 = B1 + D·(E_max2 − E_max1)`.
- Borde 2: `forward_model(canales, E_g2, A2, CE_ch2, sigma2)`.
- Devolver la suma.

Ajustar con el SA (Etapa 5), 6 parámetros Θ={A1,A2,B1,σ1,σ2,D}. Para Co-60: E_g1=1173.2, E_g2=1332.5.
**Validar:** reusar `_espectro_sintetico` de `fit_single.py` para generar 2 bordes con B1 y D conocidos
→ el ajuste debe recuperar ambos (sin el sesgo del modelo de un borde).

### Etapa 7 — Recta de calibración · nuevo módulo `calibrate.py`
Contrato sugerido:
```python
def ajustar_recta(puntos) -> {"pendiente", "offset", "R2", "energia_de_canal"}
```
- Puntos del ajuste Co-60: `(B1, 963.4)` y `(CE_ch2, 1118.1)`.
- Con 2 puntos: recta exacta (R² trivial). Reportar pendiente [keV/(mV·ns)] y offset, y **advertir**
  que no hay validación de linealidad (D6). Dejar `puntos` como lista extensible para más fuentes.
- Entregar la función inversa `energia_de_canal(ch)` = la calibración final.

---

## 4. Gotchas heredados (NO re-romper) — ver detalle en `CLAUDE.md §7`

1. **E_max** usa `2α` en el numerador (erratum de la Ec. 1). Ya correcto en `physics.py`.
2. Kernel gaussiano con **radio capado** a `(len(señal)-1)//2` (`np.convolve mode="same"` devuelve el
   largo del kernel si es más grande que la señal — pasa cuando σ se infla).
3. numpy 2.x: usar `np.trapezoid` vía `hasattr`.
4. Ajuste local frágil → esa es la razón del SA. El mínimo global SÍ está en el borde real.
5. El histograma debe tener **aire por encima del borde** (si no, el flanco se trunca y el fit falla).
6. `recortar_negativos=False` para AJUSTAR; `True` solo para GRAFICAR.

---

## 5. Flujo final objetivo (`calibrar.py` v2, un comando)

```
spectrum.espectro_neto(...)         # Etapa 0: XML → neto
   → ncc.seleccionar_intervalo(...) # Etapa 4: intervalo automático (reemplaza el manual)
   → optimize.ajustar_sa(forward_model_doble, ...)   # Etapas 5+6: ajuste robusto, 2 bordes
   → puntos (B1,963.4), (CE_ch2,1118.1)
   → calibrate.ajustar_recta(puntos)                 # Etapa 7: ganancia mV·ns ↔ keV
```

---

## 6. Fuera del código (bloqueante para calibrar BIEN)

Con solo Co-60 obtienes 2 puntos altos y juntos → ganancia aproximada, sin validar linealidad.
Conseguir **Cs-137** (borde 477 keV) —y/o Na-22— daría puntos en otras energías, permitiría anclar
la recta y comprobar la linealidad. El algoritmo corre igual con Co-60; la limitación es de datos,
no de código.
