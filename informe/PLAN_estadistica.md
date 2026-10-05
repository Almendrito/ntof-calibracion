# Plan: la estadística que usa el código y todavía no está en el informe

Inventario hecho el 25-sep-2026 recorriendo `analisis/` (la implementación
vigente) y comparándolo con `informe.tex`. La meta: **todo concepto estadístico
que el código usa tiene que poder encontrarse en el informe**, explicado al
nivel del lector objetivo de `PLAN.md` (fin de pregrado, estadística de primer
curso).

## Cómo se usa este archivo

- Una sesión = un bloque (o dos chicos). Al terminar: tachar las casillas
  (`- [x]`), poner la fecha al lado del título del bloque y anotar en qué
  sección quedó cada cosa.
- Estado de cada ítem: **✗** no está · **~** se menciona pero no se explica ·
  **✔** ya está (se lista sólo para no duplicarlo).
- Dónde va: el cuerpo lleva sólo el resumen y un puntero; el desarrollo va al
  **apéndice B** ("Estadística del ajuste y cierre"), que se propone
  reordenar por etapa de la cadena (ver bloque 0). Las fórmulas de
  propagación siguen en el apéndice C.
- Reglas de siempre (`PLAN.md` §4): ningún número tipeado (macro en
  `informe.py` o pgfmath en el `.tex`), ninguna figura pegada, ninguna cita
  sin verificar, y anclar cada concepto en un número de esta campaña.
- Al cerrar cada sesión: `python analisis/informe.py`, sin referencias sin
  resolver, y mirar la sección nueva en el PDF.

---

## Bloque 0 — Reordenar el apéndice B · ✔ 25-sep-2026

Hoy el apéndice B sigue el orden en que se escribió. Propuesta de
subsecciones, en el orden de la cadena, para que cada bloque de abajo tenga
dónde caer:

- [x] Hecho con las subsecciones que ya tienen contenido: B.1 Del conteo al
      espectro (`sec:est-conteo`) · B.2 El χ² y lo que dice (`sec:est-chi2`) ·
      B.3 La incertidumbre de los parámetros (`sec:est-incert`) ·
      B.4 Sistemáticos (`sec:est-sist`) · B.5 El cierre con verdad conocida
      (`sec:est-cierre`). Las que faltan (estimadores del pulso, optimización,
      combinar y comparar) se crean en la sesión de su bloque, para no dejar
      subsecciones vacías.
- [x] Lo que ya estaba se movió sin reescribirlo; "Lo que χ²/ndf no dice"
      pasó al final de B.2.

---

## Bloque 1 — El χ² y lo que dice (`modelo.py`, `calibracion.py`) · ✔ 25-sep-2026 (falta correr `informe.py`)

- [x] **✗ La distribución de χ²_min.** Si el modelo es correcto, χ²_min sigue
      una χ² con ndf grados de libertad: media ndf, desviación √(2·ndf). O sea
      que "χ²/ndf cerca de 1" quiere decir dentro de ±√(2/ndf): con ndf ≈ 300
      eso es ±0.08, y los χ²/ndf de 1.2 a 2.7 de la campaña están a muchas
      desviaciones. Agregar el valor p y cómo leerlo. Macro nueva: el número de
      desviaciones de la corrida de ejemplo.
- [x] **~ Grados de libertad con parámetros proyectados.** `modelo.Problema`
      cuenta ndf = bines de la ventana − (parámetros libres + amplitudes), o
      sea que A cuenta aunque se despeje. Decirlo y por qué.
- [x] **✗ Residuos normalizados ("pulls").** (dato − modelo)/σ por bin: si el
      modelo es correcto son N(0,1) sin estructura. Cómo se leen en las
      figuras de ajuste y qué forma deja el exceso de bajo canal. Evaluar si
      vale una figura (histograma de pulls de la corrida de ejemplo).
- [x] **~ χ² y verosimilitud de Poisson.** Hoy dice "para muchas cuentas es lo
      mismo que maximizar la verosimilitud". Escribir −2 ln L de Poisson y el
      χ² de cociente de verosimilitud de Baker y Cousins (ya citado) como la
      tercera opción sin sesgo, y por qué no se usa (con ~270 cuentas/bin el
      sesgo de Neyman es −0.04 %, tabla del cierre).
- [x] **✗ nnls.** Con el exceso modelado hay dos amplitudes (A y C) y se
      despejan con mínimos cuadrados con restricción de no negatividad
      (`modelo.amplitudes_optimas`): qué resuelve y por qué una amplitud
      negativa no es física. Cita a verificar: Lawson y Hanson (1974).
- [x] **~ Corte contra perfil del paisaje.** `Problema.mapa` dibuja un CORTE
      del χ² (el resto congelado en el óptimo), no el perfil. Decir la
      diferencia donde se muestre el paisaje.

Cómo quedó el bloque 1 (25-sep-2026), para la próxima sesión:

- Texto: B.2 tiene cinco párrafos nuevos ("Cuánto es «cerca de uno»", "Los
  grados de libertad", "La verosimilitud de Poisson", "Con dos amplitudes",
  "Los residuos normalizados") y B.3 uno ("Corte contra perfil"). La columna
  χ²/ndf de "Columna por columna" ahora dice cuánto es cerca de 1.
- Código: `informe.estadistica_chi2`, `figura_pulls` (→ `teoria_pulls.pdf`)
  y `macros_chi2` (17 macros: `EjNdf`, `EjChiSigmas`, `EjChiPvalor`,
  `EjChiNeyman/Pearson/BC`, `EjPull*`, `NdfMin/Max`, `ChiSigmasMin/Max`),
  enganchadas en `main`. `estadistica_chi2` avisa si el χ² reconstruido no
  coincide con el guardado.
- **Pendiente: correr `python analisis/informe.py`** en el PC. Cowork no
  puede leer los `.root` (no hay `uproot`), así que la figura de pulls y los
  números de Baker-Cousins y de los pulls salen como "?" en la vista previa
  (`\providecommand` en la cabecera de `informe.tex`). Los de ndf y valor p
  ya son los reales: 192 gdl, χ² = 393, 10.3σ, p ≈ 7·10⁻¹⁶; campaña 156–220
  gdl, 2 a 17σ.
- **Después de correrlo, revisar el párrafo de los pulls**: está escrito para
  cualquier resultado ("se leen juntos: …"). Con los números reales conviene
  decir cuál de los dos casos es (ruido de más sin estructura, o estructura
  con correlación entre vecinos).
- Cita nueva verificada: Lawson y Hanson (1974), `lawson1974`.

## Bloque 2 — Incertidumbre de los parámetros (`sa.errores_por_curvatura`) · ✔ 25-sep-2026 (falta correr `informe.py`)

- [x] **✗ El hessiano numérico.** Derivadas segundas por diferencias
      centradas; el paso se elige por parámetro (`Problema.pasos_curvatura`:
      un bin para c0, relativo para D). Por qué un paso más chico que el bin
      mide el quiebre de la plantilla interpolada y no la curvatura (el caso
      21σ contra 2σ ya está en el apéndice de historia; acá va el método).
- [x] **✗ Covarianza no definida positiva.** Qué significa (el punto no es un
      mínimo en alguna dirección, o el hessiano numérico es ruidoso), por qué
      pasa en el detector D, y el parche de `calibrar_todo.error_usable`
      (piso del 1 % de D, fila marcada). Es un parche y el informe tiene que
      decirlo donde aparezcan esas corridas.
- [x] **✗ Δχ² en dos dimensiones.** La elipse Δχ² = 1 del plano (D, σ) de la
      figura de estadística proyecta el ±1σ de cada parámetro, pero como
      región conjunta contiene el 39 %; la región del 68 % es Δχ² = 2.30.
      Aclararlo en el epígrafe y en el texto.
- [x] **~ Teorema de Wilks.** El test de la ordenada (`calibracion.test_ordenada`)
      convierte Δχ² entre modelos anidados en desviaciones como √Δχ². Enunciar
      el resultado (Δχ² ~ χ² con k grados de libertad, k = parámetros extra) y
      cuándo vale. Cita a verificar: Wilks (1938) o Cowan (ya citado).
- [ ] **✔** curvatura contra perfil, Δχ² = 1, escalado por √(χ²/ndf),
      correlación D–σ.

Cómo quedó el bloque 2 (25-sep-2026):

- Texto, todo en B.3: "Cómo se calcula la curvatura" (diferencias centradas,
  C = 2H⁻¹ y por qué el 2, el paso por parámetro), "Cuando la covarianza no
  sale" (definida positiva, el parche de `error_usable`), "La barra de uno y
  la región de dos" (39 % contra 68 %, 1 / 2.30 / 3.53) y "Comparar modelos
  anidados: el teorema de Wilks" (con las tres condiciones, incluida la del
  borde C ≥ 0). El epígrafe de la figura de estadística nombra la elipse
  nueva. Cita nueva verificada: Wilks (1938), `wilks1938`.
- Código: `modelo.PASO_CURV_REL` (antes un 1e-3 suelto en
  `pasos_curvatura`; mismo valor, ningún resultado cambia), la elipse de
  Δχ² = 2.30 en el panel (b) de `figura_estadistica`, y
  `informe.macros_incertidumbre` (`PasoCurvRel`, `DeltaChi{Uno,Dos,Tres}D`,
  `ProbElipseUno`, `NErrorParche`, `CierreDentroElipse`).
- Números ya medidos (no dependen de la plantilla): hoy ninguna corrida
  necesita el parche (las 29 covarianzas guardadas son definidas positivas),
  y 17 de las 30 realizaciones del cierre caen en la región del 68 %.
- Pendiente: correr `informe.py` para que la figura muestre la elipse de
  puntos (la vista previa todavía tiene la figura vieja).

## Bloque 3 — Cierre con verdad conocida y sistemáticos (`informe.cierre_conocido`, `sistematicos.py`)

- [ ] **✗ El error de la cobertura.** Con n = 30 realizaciones y p = 0.683, la
      cobertura esperada es binomial: 20.5 ± 2.5. El 18/30 que da Neyman está
      a 1σ de lo esperado; decirlo con el número para no leer de más el
      "algo corta". Macros: cobertura esperada y su σ.
- [ ] **✗ Estimadores muestrales del cierre.** Media, desviación estándar con
      n−1, error de la media (σ/√n) para el sesgo, y covarianza muestral para
      la elipse empírica. Y el nombre del procedimiento: muestreo Poisson del
      modelo ajustado = "toy Monte Carlo" / bootstrap paramétrico.
- [ ] **✗ Semidispersión.** `sistematicos._semidispersion` = (máx − mín)/2
      como error de una variación acotada. Alternativa (máx − mín)/√12 si se
      la supone uniforme; decir cuál se usa y por qué.
- [ ] **✗ Suma en cuadratura por alcance.** Los términos se suman en
      cuadratura separados en común y punto a punto
      (`sist_punto`, `sist_comun`, `hypot`), lo que supone independencia.
      Hoy el apéndice B sólo lo dice cualitativamente.
- [ ] **✗ Tabla de presupuesto de la corrida de referencia** (pendiente
      desde `PLAN.md` B9): generarla en `informe.py` desde
      `sistematicos.presupuesto`.
- [ ] **✗ Variantes que se miden y no se suman.** El criterio: una variante
      que los datos rechazan (χ²) o que la física descarta (el área por
      convergencia) no es un sistemático. Cita a verificar: Barlow,
      "Systematic errors: facts and fictions" (2002).

## Bloque 4 — Estimadores del pulso (`pulso.py`, `area.py`)

- [ ] **~ Mediana contra media** para la línea base: robustez (una muestra
      contaminada no la mueve) a cambio de algo de eficiencia con ruido
      gaussiano. Hoy se dice "es inmune"; falta el porqué.
- [ ] **~ MAD y el factor 1.4826.** Es 1/Φ⁻¹(3/4): hace que la MAD estime σ si
      el ruido es gaussiano. Hoy aparece el número sin la razón.
- [ ] **✗ Umbrales en unidades de ruido.** 5σ para "sin pulso", 3σ durante
      4 ns para la convergencia: probabilidad de un falso cruce por muestra y
      por qué exigir un tramo seguido baja esa probabilidad.
- [ ] **✗ Onset sub-muestra y jitter de muestreo.** Interpolación lineal del
      cruce del 10 % (`pulso.cruce_onset`); alinear por la muestra agrega un
      jitter uniforme de un paso, σ = Δt/√12, que ensancha el pulso promedio.
- [ ] **~ Promedio de pulsos.** El ruido del promedio baja como 1/√N; por qué
      se normaliza por amplitud y por qué se promedia con conteo por punto en
      vez de rellenar con ceros (`plantilla_pulso`).
- [ ] **✗ Error de la mediana por tramo** (`area.perfil`: 1.2533·σ/√n) en el
      control área/amplitud. Decir de dónde sale el 1.2533 (√(π/2)).
- [ ] **~ Percentiles como estimadores.** El rango del histograma (1.6 ×
      percentil 99.5) y el espacio tras el onset (percentil en eventos sobre
      la mediana): por qué un percentil y no el máximo.

## Bloque 5 — El ajuste de la plantilla del pulso (`pulso.ajustar_plantilla`, `plantilla_pulso`)

Existe en el código y tiene un lote corrido (`analisis/figuras/calibraciones_plantilla/`),
pero **no está en el informe ni en `CLAUDE.md`**. Antes de escribirlo hay que
decidir si entra como método o como control.

- [ ] **Decisión**: mirar los resultados del lote y anotar acá si el método
      entra al informe (y dónde) o queda como control.
- [ ] **✗ Mínimos cuadrados lineales en (A, b)** con solución cerrada para
      cada t0, y la búsqueda de t0 en grilla gruesa, fina y refinamiento
      parabólico.
- [ ] **✗ Por qué pesos uniformes por muestra** (el ruido es por celda del
      DRS4) y cuándo deja de ser óptimo (ruido correlacionado entre muestras:
      filtro óptimo con la covarianza del ruido).
- [ ] **✗ χ² reducido por evento** (ndf = n − 3) como bandera de apilamiento.
- [ ] **✗ Extrapolación de la cola** y su supuesto: forma independiente de la
      amplitud. El control es que área_plantilla/área_compuerta sea plano en
      amplitud.
- [ ] **✗ Alineación iterativa** de la plantilla (por t0 ajustado en vez de
      por cruce) y por qué reduce el ensanchamiento.

## Bloque 6 — Del conteo al espectro (`espectro.py`)

- [ ] **~ El tiempo vivo como estimador** y su incertidumbre (conteo de
      triggers); el recorte del scaler a 10× su mediana como estimador robusto
      frente a las ráfagas. Mostrar que su error no importa acá (k_f ≈ 6·10⁻⁴).
- [ ] **✗ Cuándo falla la aproximación gaussiana del χ².** Con pocas cuentas
      por bin (regla práctica de 5–10) el χ² deja de seguir su distribución;
      la campaña tiene una mediana de ~270 cuentas por bin en la ventana (macro `\CuentasBin`). Decir qué
      pasa en las colas de la ventana.
- [ ] **~ `fin_del_umbral`**: máximo local del espectro suavizado como
      estimador y su incertidumbre (un bin).
- [ ] **✔** Poisson por bin, varianza de la resta, piso 1, negativos
      conservados.

## Bloque 7 — Combinar mediciones y comparar modelos (`registro.py`, `voltaje.py`, `informe.py`)

- [ ] **✗ Recta contra potencia.** Tienen el mismo número de parámetros, así
      que Δχ² se compara directo; decir qué haría falta (test F, AIC) si no lo
      tuvieran.
- [ ] **✗ `voltaje.discriminacion`**: el desvío de la cuerda,
      (1/8)·b(b−1)·(ΔV/V̄)², deducido, y por qué con 1600–1700 V el test apenas
      separa las dos hipótesis.
- [ ] **✗ Sesgo de ajustar en logaritmo**: E[ln D] ≠ ln E[D]; con errores del
      0.2 % es despreciable, decir el orden de magnitud.
- [ ] **✗ `_pedestal_o_ganancia`**: media pesada de las ordenadas y su χ² de
      constancia, con el umbral heurístico "χ² < 2·gdl + 5"; compararlo con un
      valor p.
- [ ] **✗ El cierre de `voltaje._cierre` y `registro.cierre`** (verdad
      conocida para los tests A y B): mencionar que existen y qué prueban.
- [ ] **✗ Las métricas del paper, definidas.** R² = 1 − SS_res/SS_tot, tal como
      lo calcula `informe.metricas_paper` (con SS_tot centrado aunque la recta
      pase por el origen, que es una de las definiciones posibles y la que
      hace que valga 0.94 con dos puntos); la similitud coseno; la
      "dispersión" de la robustez, que es (máx − mín)/media.
- [ ] **✔** media pesada (D conjunta), mínimos cuadrados ponderados en forma
      cerrada, propagación (apéndice C), pivote geométrico.

## Bloque 8 — Lo que queda de NCC y recocido (`ncc.py`, `sa.py`)

- [ ] **~ ρ no es un test.** Con el suavizado, los puntos no son
      independientes (una ventana de Savitzky-Golay de 21 bines deja del orden
      de N/21 puntos efectivos); por eso ρ sirve para ordenar ventanas y no
      tiene un valor p.
- [ ] **✗ El kernel gaussiano discreto**: normalizado a suma 1 y truncado a
      n·σ; el salto del χ² cuando el radio cambia en un bin (comentario en
      `ncc._gaussiana`) y por qué eso rompe el hessiano.
- [ ] **~ Nelder–Mead**: qué es (símplex sin derivadas) y su criterio de
      corte. Cita a verificar: Nelder y Mead (1965).
- [ ] **✗ Cuántos arranques bastan**: si la cuenca del mínimo global ocupa una
      fracción p de la caja, la probabilidad de no caer nunca en ella es
      (1 − p)^n; con 8 arranques, qué p queda cubierto.
- [ ] **✗ Reproducibilidad**: semilla fija del generador (`semilla=7`) y el
      control de que recalibrar el lote reproduce las D exactas.
- [ ] **✔** Pearson, invariancia, Savitzky-Golay, Metropolis, Boltzmann,
      balance detallado, enfriamiento, T_init medida, tasa de aceptación.

## Bloque 9 — Transversal (al final)

- [ ] **Glosario**: MAD, percentil, valor p, ndf, pull, hessiano, cobertura,
      toy Monte Carlo, semidispersión, suma en cuadratura, nnls, Nelder–Mead,
      Wilks, R², similitud coseno.
- [ ] **Bibliografía**, verificando cada una antes de agregarla: Lawson y
      Hanson (nnls), Nelder y Mead, Wilks, Barlow (sistemáticos), la revisión
      de estadística del PDG, y James, *Statistical Methods in Experimental
      Physics*, como texto de referencia general junto a Cowan.

## Bloque 10 (opcional) — La implementación anterior (`codigo/`)

Sólo si el informe llega a citar números de `codigo/`. Hoy no cita el MCMC.

- [ ] MCMC de Metropolis–Hastings (`codigo/ajuste/mcmc.py`): propuesta con la
      covarianza del ajuste escalada por 2.4²/d, R̂ de Gelman–Rubin, tiempo de
      autocorrelación integrado y n_eff, mediana y percentiles del posterior.
- [ ] El χ² de Pearson iterado y el modelo extendido (fracción de un solo
      scattering w, exceso exponencial) como antecedente del presupuesto.
- [ ] Decidir: un párrafo en el apéndice de historia, o nada.

---

## Orden propuesto de sesiones

| # | Bloques | Por qué en este orden |
|---|---|---|
| 1 | 0 + 1 | Da la estructura y lo más usado: cómo juzgar un χ² |
| 2 | 2 | La barra de error es lo que se reporta en cada tabla |
| 3 | 3 | Cobertura y sistemáticos cierran la pregunta "¿el error es real?" |
| 4 | 4 + 6 | Estimadores del pulso y del conteo (menos riesgo, más texto) |
| 5 | 7 | Tests y métricas (resultados) |
| 6 | 5 | Depende de decidir si la plantilla del pulso entra |
| 7 | 8 + 9 | Cierre: detalles del NCC/SA, glosario y bibliografía |
| 8 | 10 | Opcional |
