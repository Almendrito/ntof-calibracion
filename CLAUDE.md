# CLAUDE.md — Calibración de centelladores plásticos por Compton edge (método NCC-SA)

> Archivo de contexto para retomar el proyecto sin partir de cero.
> Colócalo en la raíz del repositorio (Claude Code lo lee como contexto).

---

## 1. Qué es este proyecto

Replicar y adaptar el método de calibración de energía por **Compton edge** del paper
**Kong et al., "The Automatic Calibration Method of the Compton Edge Based on Normalized
Cross-correlation and Simulated Annealing Algorithm", IEEE TNS 2025 (arXiv:2507.16282)**,
para **calibrar la ganancia de un centellador plástico** leído con un SiPM/MPPC.

**Objetivo concreto:** convertir el eje medido —*área del pulso en mV·ns*— a *energía en keV*,
usando el borde Compton (los centelladores plásticos son de bajo Z, dominados por Compton,
y **no tienen fotopico**, así que la calibración clásica por fotopico no aplica).

> ⚠️ **25-ago-2026 — el número cambió por un factor 3.4.** Todo lo que dice
> `D ≈ 4.3 mV·ns/keV` en este archivo es de la Etapa 8 y está MAL: el área del pulso se
> calculaba con una línea base que se comía la subida del propio pulso (§7.33) y el fondo
> se restaba por tiempo de pared cuando la corrida estuvo viva el 40% (§7.34). El resultado
> vigente es **D_luz = 1.276 ± 0.023 mV·ns/keVee** (Etapa 9, `calibracion.py`). Las secciones
> viejas se dejan porque el MÉTODO sigue valiendo; los números de canal, no.

Contexto mayor: proyecto de **detección de neutrones / nToF** (tesis). Esta calibración es la
**rama de electrones/gammas**. La rama de protones de retroceso (neutrones) es no-lineal y va
aparte; NO mezclar ambas en esta recta.
⚠️ Corrección (24-ago-2026): decir que "Birks casi no muerde" en esta rama es falso. La
simulación da `Lvis/Edep` = 0.954 arriba de 700 keV y 0.863 a 50-100 keV (§7.17). Por eso la
calibración tiene DOS escalas que no hay que mezclar: keVee (luz) y keV (energía depositada),
separadas por un 4.5%. Ver `calibracion.py`.
⚠️ Corrección (25-ago-2026): que D_luz/D_edep coincida con 1/0.954 se venía vendiendo como
"control de consistencia independiente". **No lo es.** Lvis es Edep multiplicada por el
quenching, y arriba de 700 keV el quenching es plano, así que ajustar los mismos datos con
las dos plantillas TIENE que dar ese cociente. Lo que verifica es que el ajuste está
dominado por la zona donde el quenching satura — útil, pero mucho menos de lo que decía.

---

## 2. Hardware y datos

- **Detector:** ✔ **CONFIRMADO (24-sep-2026, por el usuario): plástico BC-408 cilíndrico de 2"×2" + fotomultiplicador XP2262/B (51 mm, 12 dínodos).** Con n = 12, el b = 7.7–10.4 del test B da k = b/12 ≈ 0.64–0.87 por etapa. Lo que sigue en este ítem es la historia de cómo se llegó a eso.
  (antes decía: centellador plástico + un fotodetector **sin confirmar**.)
  ⚠️ (28-ago-2026, §7.45) Todo este archivo decía "SiPM/MPPC MicroFC-30020", y **no puede
  ser**: el `1550` del nombre del archivo es el ALTO VOLTAJE, y un MicroFC se polariza a
  ~25 V. A 1550-1700 V es un PMT. Con eso se cae también la lectura de la Tabla I del paper
  que hacía este punto (`σ_gain = ENF·√N_pe` es un término de SiPM). Confirmar qué tubo es
  antes de modelar el ensanchamiento a partir de primeros principios.
- **Punto de operación** (§7.45): D no es del detector, es del detector A UN VOLTAJE.
  `analisis/rutas.Corrida` guarda `hv`, `atenuador`, `umbral` y `muestreo`, y es la clave
  por la que se agrupa. El catálogo tiene hoy **25 corridas**, por fecha de medición
  (`rutas.fecha`, leída del archivo): 1 del 23-jun (Co-60 a 1550 V), 3 del 26-ago (Cs-137 a
  1600/1650/1700 V), 12 del 1-sep, 3 del 2-sep y 6 del 8-sep, que son los cuatro detectores
  A-D a 1600/1650/1700 V en binario (§7.51, §7.61).
  ✔ **El atenuador SÍ está registrado: 20 dB en TODAS las corridas** (confirmado por el
  operador, 22-sep-2026; §7.70). Están para que el pulso del PMT no pase el rango de
  entrada del DRS4. El área se integra a la ENTRADA del digitalizador y **ningún módulo
  corrige por la atenuación**, así que todas las D de este archivo son "por mV·ns a la
  entrada del DRS4"; la carga a la salida del PMT es **10 veces** mayor (20 dB = factor 10
  en amplitud, y el área escala igual).
- **Eje de "canal":** área del pulso integrada, en **mV·ns** (no es un ADC channel clásico).
- **Datos disponibles:** las 25 corridas del catálogo, copiadas a `datos/Detector X/<V>/`
  (§7.61), así que el proyecto entero corre sin el disco externo. Archivos de junio
  (CHN2 es el canal usado):
  - Muestra: `A Co60 5GaS 50k 02mV trigg 1550 att.xml` — 1.76 GB, 50000 eventos, 699 s → 71.5 cps
  - Fondo:   `A BG 5GaS 10k 02mV trigg 1550 att.xml` — 353 MB, 10000 eventos, 986 s → 10.1 cps
- **Estructura del XML** (DRSOSC): `DRSOSC > Event > Board_2859 > CHN{n} > Data`, un tag por
  línea, `Data` = `"t,v"` con t en ns y v en mV; `Time` por evento da el reloj. Encoding
  ISO-8859-1. ⚠️ **El canal NO es fijo** (§7.36): junio = `CHN2` + `<Scaler1>`, la sesión del
  26-ago-2026 = `CHN1` + `<Scaler0>`. `spectrum.detectar_canal` lo saca del primer evento.
- **Y hay un SEGUNDO formato** (§7.51): la sesión del 1-sep-2026 salió con "Save as binary"
  del DRSOSC y **sin extensión en el nombre** (`Co60`, `Cs137`, `B1700V`). Formato byte a
  byte en `docs/formato_binario_drs4.md`. Las dos implementaciones lo leen y ninguna decide
  por el nombre: `analisis/drs.py` y `codigo/datos_io/spectrum.py` miran los primeros 4
  bytes (`DRS2`/`TIME`). ⚠️ En los diálogos de archivo hay que dejar el filtro en `*.*`,
  si no esas carpetas se ven vacías (§7.57).
- **Las dos corridas NO son intercambiables tal cual** (descubierto 25-ago-2026): se tomaron
  días distintos y con retardo de trigger distinto. El pulso está en **t ≈ 48 ns** en el
  Co-60 y en **t ≈ 120 ns** en el fondo. Con el área vieja eso significaba que muestra y
  fondo ni siquiera medían lo mismo (§7.33). El integrador con compuerta localiza el onset
  por evento, así que es inmune a eso.
- **Tiempo vivo** (§7.34): el `<Scaler1>` de cada evento da la tasa de trigger REAL de la
  placa. Co-60: 170 Hz reales contra 71.5 ev/s guardados → **40% vivo, 283 s de 699**.
  Fondo: 10 Hz contra 10.1 ev/s → **100% vivo**.
- **Zona útil del espectro (área vieja, ya no vale):** entre 0 y ~6000 mV·ns. El plateau Compton va
  de ~2000 a ~4200 y la caída llega hasta ~5800. Hay eventos sueltos hasta 46000 mV·ns
  (cósmicos/pile-up; el fondo tiene proporcionalmente más).
  Con el **área corregida** toda la física vive entre 0 y ~1900 mV·ns, y los dos bordes caen
  en **1229 y 1427 mV·ns**. `RANGO_MAX` pasa a 2500 con 500 canales (5 mV·ns por bin).
- **Co-60 → DOS bordes Compton solapados:** 963.4 keV (de γ=1173.2) y 1118.1 keV (de γ=1332.5).

---

## 3. El método (resumen físico)

Cadena conceptual del paper, en orden:

1. **E_max (borde teórico)** — energía conocida sin ajustar, dada la línea γ:
   `E_max = 2α·E_γ / (1 + 2α)`, con `α = E_γ / m_e c²` (m_e c² = 510.999 keV).
   ⚠️ La **Ec. 1 impresa del paper está mal** (da la energía del fotón retrodispersado,
   `E_γ/(1+2α)`); la correcta lleva `2α` en el numerador. Verificado contra los propios
   datos del paper (Fig. 3): 963 y 1118 keV para Co-60.
2. **Lineshape teórica dσ/dE_e** — Klein-Nishina (Ec. 2) pasada a energía del electrón vía
   Jacobiano → **Ec. 10**. Es la "plantilla" S2. Tiene un plateau, una acumulación afilada
   y un escalón (corte) en E_max.
3. **Forward model** — la lineshape llevada a espacio de canal (parámetro `B` = canal del
   borde) y convolucionada con una gaussiana de ancho `σ` (broadening del detector). Ec. 8-9.
   Parámetros: `A` (amplitud), `B`/`edge_ch` (canal del borde), `σ` (resolución).
4. **NCC (Pearson, Ec. 15)** — elige AUTOMÁTICAMENTE el intervalo del borde correlacionando
   la forma de los datos contra la plantilla. Inmune a offset y a escala → compara FORMA.
   Lo que el paper dice, literal (Sección II-B): genera **100 puntos de referencia (S2)**
   muestreando la lineshape "from the Compton plateau to the edge minimum", y **cada punto
   no nulo del dato es un posible inicio de intervalo**, con los intervalos candidatos
   **interpolados a 100 puntos (S1)**. O sea el paper barre inicio y fin LIBRES, que es el
   `modo="libre"` de `ncc.buscar_intervalo`. Nuestro default es `modo="origen"` (barrido 1D)
   y un punto por bin en vez de 100: son desviaciones NUESTRAS, justificadas en §7.10 y
   §7.27 con la propia Ec. 8 del paper y con medición. **[fuera del paper]**
5. **Simulated Annealing (SA)** — optimizador GLOBAL de {A, B, σ} que escapa mínimos locales
   (Metropolis + enfriamiento geométrico). Params del paper: T_init=2000, α∈[0.85,0.99]
   (usan 0.95), 100 iters/temperatura, corte T<1e-6 o |ΔE|<0.001 por 10 ciclos.
6. **Bordes solapados (Co-60)** — Θ = {A1, A2, B1, σ1, σ2, D}, con el 2º borde amarrado al 1º:
   `CE_ch,2 = B1 + D·ΔE`. **D es la pendiente de calibración** (canal↔energía) que buscamos.
   ⚠️ ΔE es la diferencia de los **bordes** (154.7 keV), no la de las líneas gamma (§7.15).
7. **Recta de calibración** — juntar puntos (canal, E_max) → ajuste lineal → ganancia mV·ns↔keV.

Ecuaciones clave del paper: 1 (E_max), 2 (Klein-Nishina), 8-9 (convolución+kernel),
10 (dσ/dE_e), 11-13 (solapados), 14 (loss mín. cuadrados), 15 (Pearson/NCC).

⚠️ **Inconsistencia del paper:** la Ec. 14 dice que SA minimiza mínimos cuadrados, pero el
texto de Metropolis dice que E(x) es "NCC-based matching error". Hay que ELEGIR uno al
implementar SA (para calibración de forma, NCC suele ser más robusto).

⚠️ **Dónde termina el paper.** Los 7 pasos de arriba son el paper y nada más. La réplica
literal vive en **`replica_paper/`**, cuya regla es "acá sólo entra lo que el paper dice"
(Klein-Nishina analítica, sin Geant4, sin MCMC, sin modelo extendido). Todo lo de `codigo/`
y `analisis/` es **extensión nuestra**, y el índice de qué es de quién está en §7.71. Si una
explicación mezcla las dos cosas sin decirlo, está mal escrita: eso es lo que §4 pide marcar
con **[fuera del paper]**.

---

## 4. Enfoque de trabajo (IMPORTANTE — mantenerlo)

El usuario (Mateo, física) está **APRENDIENDO mientras construye**. El estilo acordado:

- **Explicar el PORQUÉ** de cada patrón, no solo entregar código.
- **Construir por etapas**, cada una **ejecutable y verificable por sí sola**.
- **Validar contra datos SINTÉTICOS con verdad conocida** ANTES de correr sobre datos reales
  (generar espectro con B y σ conocidos, confirmar que el código los recupera). El sintético es
  solo un test unitario: el resultado que importa se saca de los datos REALES contra la
  plantilla teórica (Klein-Nishina), que es exactamente lo que hace el paper.
- **NO cerrar las etapas con preguntas de comprensión** (el usuario lo pidió explícitamente,
  24-ago-2026): avanzar en el desarrollo, explicar el porqué al pasar, sin interrogatorio.
- Idioma: **español**. Comentarios de código en español.
- Anclar explicaciones en el paper; marcar explícitamente lo que es **[fuera del paper]**.

---

## 5. Arquitectura del código

Módulos (cada uno corre standalone con `python3 <archivo>.py` y genera una figura demo):

| Archivo (dentro de `codigo/`) | Etapa | Contenido / funciones clave |
|------------------|-------|-----------------------------|
| `etapas/physics.py` | 1     | `compton_edge_energy(E_gamma)` (Ec.1 corregida), `dsigma_dE(E_e, E_gamma, A)` (Ec.10), **`FUENTES` + `NOTAS_FUENTE` + `Fuente`** (§7.39, §7.41: qué líneas gamma tiene lo medido, sus bordes, `E_ref`, `doble`, más `nota()` y `avisos()`). 9 presets: Am-241, Co-57, Ba-133, Cs-137, Mn-54, Na-22, Co-60, K-40, Tl-208. Pura teoría, trabaja en keV. |
| `etapas/model.py`   | 2     | `forward_model(canales, E_gamma, A, edge_ch, sigma, pad=True)` (lineshape→canal + convolución **sobre grilla extendida**, §7.14), `punto_mitad_altura(canales, y)` (método clásico 50%), `_kernel_gaussiano(...)` (gaussiana normalizada, radio capado). |
| `etapas/fit_single.py` | 3     | `ajustar_borde_unico(canales, cuentas, E_gamma, intervalo, p0)` (least_squares local — INGENUO/FRÁGIL), `_espectro_sintetico(...)` (generador de prueba). |
| `ajuste/sim_template.py` | 7     | Forward folding con la simulación Geant4: `cargar_simulacion(rama, filtro, bins, rango, arbol)` (uproot; `bins=175` para SA/MCMC, §7.26), `sigma_energia(E, alpha, beta)`, `modelo_simulado(...)` (matriz de respuesta, convolución de ancho VARIABLE), `ajustar_forward_folding(..., errores=, E_ref=)` (χ² ponderado + covarianza completa en `cov`; `E_ref` es sólo para reportar la resolución: `sigma_E_ref`/`fwhm_rel_ref`, antes `*_963`), `ajustar_pearson(..., var_fn)` (χ² de Pearson iterado, mata el sesgo de §7.31), `costo_sim_varpro(...)` (costo para el SA con A proyectada), `ajustar_sa(...)` (SA + pulido local), **`lineas_primarias(ruta, arbol)`** (qué gammas simula el `.root`, leyendo `E0_keV`) y **`comparar_con_fuente(...)`** (§7.41). |
| `rutas.py`          | —     | **Dónde vive cada cosa** (§7.43): `RAIZ`, `DATOS`, `FIGURAS`, `DOCS`, `CACHE`, los defaults `MUESTRA`/`FONDO`/`SIM`/`SIM_CS137`, y `dato()`/`figura()`/`cache()`/`relativa()`/`absoluta()`. Todo se resuelve desde `__file__`, así que el proyecto corre desde cualquier directorio. `python codigo/rutas.py` lista lo que ve y avisa si falta un dato. |
| `gui.py`            | —     | **Visor por etapas** (tkinter + matplotlib). 14 pestañas: la de **Archivos** (elige el `.root`, los dos `.xml` y la FUENTE; `TabArchivos`, `Datos.fijar_archivos/fijar_fuente/huella`, preferencia en `.cache/archivos.json`), 11 de análisis (controles + gráfico + explicación, todas parametrizadas por la fuente elegida), la de **★ Calibración** (que ahora CALCULA, §7.41) y la de **glosario de parámetros** (`GLOSARIO`, tabla filtrable con 26 entradas). Fuera de las pestañas: **`ajuste_extendido(...)`** (la receta de la Etapa 9 como función, la usan la 10 y la de Calibración), **`Datos.compat_simulacion/exigir_sim_compatible/componentes`**, y el **registro de calibraciones** (`cargar_registro`, `guardar_en_registro`, `recta_por_origen`, `.cache/calibraciones.json`). Recalcula de verdad, no muestra PNGs. Cálculos pesados en hilo aparte, con **barra de progreso** abajo (`App.progreso_iniciar/actualizar/terminar`, `Pestana.progreso(texto, frac)` y `Pestana.fase(a, b)`, `Datos.avisar`; §7.42). `activar_dpi_awareness()` + `App.configurar_nitidez()` para pantallas con escalado. |
| `calibracion.py`    | —     | **EL RESULTADO CONGELADO del Co-60 con este SiPM.** `D_LUZ`, `D_EDEP`, `ALPHA`, `area_a_energia(area, escala)`, `energia_a_area(...)`, `resolucion_sigma/fwhm_rel(E)`. Es lo que se importa desde otros proyectos. Para OTRA fuente u otro detector el número lo calcula la pestaña ★ Calibración y queda en `.cache/calibraciones.json` (§7.41); este archivo no se actualiza solo. |
| `etapas/sa.py`      | 5     | `recocido_simulado(costo, x0, cotas, ...)` (Metropolis + enfriamiento geométrico, `T_init="auto"`, `paso_adaptativo`), `sa_multiarranque(...)` (**usar esta**), `costo_minimos_cuadrados` (Ec.14), `costo_ncc` (Ec.15), `costo_dos_bordes_varpro` (4p con A1,A2 por nnls), `proyectar_amplitud(...)`. |
| `etapas/fit_double.py` | 6     | `modelo_dos_bordes(...)` (Ec.13), `componentes(...)`, `delta_E()` (Ec.12 corregida), `semilla_por_barrido(...)` (barrido grueso con A1,A2 exactas por nnls), `ajustar_dos_bordes(..., ligar_origen, ligar_sigmas)` → Θ={A1,A2,B1,σ1,σ2,D}. |
| `etapas/clasicos.py` | —  | **Los dos métodos convencionales contra los que compara el paper** (§7.44): `metodo_diferencia(...)` (mínimo de la derivada, con `ventana="auto"` deducida del ancho del flanco), `metodo_fwhm_gauss(...)` (ajuste de escalón borroneado `escalon_gaussiano` = erfc + fondo lineal → borde y FWHM), y `comparar(..., B_ref=)` que corre los tres junto al 50% de altura de `model.punto_mitad_altura`. No son el resultado: son la vara. |
| `etapas/ncc.py`     | 4     | `plantilla_referencia(E_gamma, sigma_rel, ...)` (S2 Klein-Nishina), `plantilla_simulada(centros, pesos, ...)` (S2 desde Geant4, §7.28), `pearson(s1,s2)` (Ec.15), `buscar_intervalo(..., modo="origen"\|"libre", plantilla=, E_ref=, n_puntos="auto")` → `lo,hi,rho,edge_ch,sigma,p0,D_seed,p0_sim`, `buscar_intervalo_simulado(..., E_ref=, E_ref2=)` (**la ruta de la Etapa 8**; con `E_ref2=None`, o sea fuente de una línea, la ventana se cierra como en Klein-Nishina). Selección AUTOMÁTICA del intervalo. |
| `ajuste/mcmc.py`    | 8     | `metropolis(log_post, x0, cov_prop, ...)` (propuesta = covarianza del ajuste × 2.4²/d), `correr_cadenas(...)` (arranques dispersos), `gelman_rubin` (R̂), `autocorr_integrada` (τ_int → n_eff), `resumen(...)`, `contornos_2d` / `elipse_covarianza` (dibujo sin `corner`), `log_posterior_sim(..., var_fn=)`, `posterior_calibracion(...)` (**el orquestador**). |
| `datos_io/spectrum.py` | 0     | **`formato(ruta)`/`es_binario(ruta)`** (XML del DRSOSC o binario del DRS4, por los primeros 4 bytes; el binario se delega en `analisis/drs_bin`, §7.57), **`pares_a_pulso(pares)`** (unifica los dos formatos en `(t, v)`), **`detectar_canal(ruta)`** (qué `CHN{n}`/`Scaler{n-1}` trae el archivo, §7.36), `espectro_neto(..., metodo="gate", canal=None)` (XML→espectro neto; devuelve `k_m`, `k_f`, `t_vivo_*`), **`_area_pulso_gate(...)`** (área con línea base pre-onset y compuerta, §7.33 — es el observable bueno), `_area_pulso(...)` (LEGACY, roto, solo para reproducir lo viejo), **`tiempo_vivo(...)`** (§7.34), `varianza_predicha(fondo, k_m, k_f)` (§7.31), `leer_areas_y_tiempo(...)` (escáner + caché `.areas_gate.npz`), `_iter_eventos(...)` (da también `Scaler1`), `suavizar(...)`. Los tres escalones del parseo aceptan `callback=` para la barra de progreso (§7.42). |
| `ajuste/modelo_ext.py` | 9     | **El modelo que mide los sistemáticos.** `cargar_componentes(...)` (separa Geant4 en 1-scattering y múltiple), `modelo_ext(canales, centros, h1, hm, p)` con Θ={A,D,α,w,C,λ,β,a_rel}, `sigma_energia(E,α,β,a_rel)` (3 términos), `ajustar(..., libres=[...], E_ref=)` (χ² de Pearson iterado, subconjunto de parámetros libre), `log_posterior(...)` para el MCMC. |
| `scripts/calibrar_v2.py` | run   | **EL SCRIPT BUENO (Etapa 9).** XML → área con compuerta → tiempo vivo → escalera de modelos → barrido de ventana → MCMC → `calibracion_v2.png`. `--rapido` salta el MCMC. |
| `documentacion.py`  | —    | Genera **`analisis_paso_a_paso.pdf`** (reportlab): el recorrido completo del análisis en orden de ejecución, con qué archivo y qué función hace cada cosa, las trampas ya resueltas y una lista priorizada de qué mejorar. Los números NO están tipeados: las constantes se leen de los módulos (`spectrum.GATE_PRE`, …), las figuras se generan de los datos reales y el resultado sale de correr el ajuste, así que el PDF no puede quedar desincronizado del código. `--sin-figuras` para iterar rápido sobre el texto. |
| `scripts/calibrar_ncc.py` | run   | Cadena de las Etapas 4-6: XML → neto → NCC → ajuste sembrado → `calibracion_ncc.png`. Superado por `calibrar_v2.py` (sus números salían del área rota). `--no-show` para no abrir ventanas. |
| `scripts/calibrar.py` | run   | Versión vieja e interactiva: pide el intervalo por teclado. Superada por `calibrar_ncc.py`. |
| `../verificaciones/lecturaxml.py` | —     | Código original del usuario (referencia histórica; superado por `spectrum.py`). |

Dependencias entre módulos (los imports son por paquete, `from etapas.physics import ...`):
`datos_io.spectrum` → (independiente) · `model` → `physics` ·
`fit_single` → `model`, `physics` · `ncc` → `physics`, `model` (+ `fit_single` solo en su demo) ·
scripts → `spectrum`, `fit_single`, `model`.

Requiere: `numpy`, `scipy`, `matplotlib`.

---

## 6. Convenciones y decisiones de diseño (respetar)

- **Energía siempre en keV.** El "canal" es área del pulso en **mV·ns**.
- **`A` absorbe el prefactor** `π r0²/(m_e c²)` de la Ec. 10: para NCC/ajuste solo importa la FORMA.
- **Kernel gaussiano NORMALIZADO** (suma 1): conserva área → `A` y `σ` quedan desacoplados
  (subir σ ensancha y baja sin cambiar cuentas totales).
- **Binning COMÚN** para muestra y fondo (mismo `range` fijo): sin esto, la resta de fondo es
  físicamente inválida (bins distintos = distintos mV·ns).
- **Restar línea base** de cada pulso antes de integrar el área (evita ancho artificial).
- **`recortar_negativos=False`** para AJUSTAR (neto real); `True` solo para GRAFICAR.
- **Eje de canales uniforme** es obligatorio para la convolución (`np.convolve`) y para el
  barrido vectorizado del NCC.
- **El mapa canal↔energía pasa por el ORIGEN** (`E_e = canal·E_max/B`, Ec. 8): área cero ↔
  energía cero, sin offset. Todo el código debe respetarlo — incluida la selección de intervalo
  del NCC (§7.10). Si algún día aparece un pedestal, se agrega explícitamente, no por la puerta
  de atrás dejando la ventana libre.
- **`RANGO_MAX` siempre FIJO** (8000 mV·ns), nunca `None`: ver §7.7.

---

## 7. Hallazgos y trampas YA resueltas (no re-romper)

1. **Erratum Ec. 1:** usar `2α` en el numerador para E_max (ver §3). Está bien en `physics.py`.
2. **`np.convolve(mode="same")` devuelve el largo del KERNEL si el kernel > señal** →
   `_kernel_gaussiano` capa el radio a `(len(señal)-1)//2`. Pasa cuando σ se infla.
3. **numpy 2.x no tiene `np.trapz`** → usar `np.trapezoid` vía `hasattr` (hecho en `spectrum.py`).
4. **El ajuste local (`fit_single`) es FRÁGIL:** `least_squares` cae en mínimos locales; el
   resultado depende del punto inicial (= centro del intervalo). Se mapeó el paisaje de la loss:
   el mínimo GLOBAL sí está en el B verdadero, pero desde un arranque malo el optimizador queda
   ~10% off con σ inflado. → **Esta es la motivación de la Etapa 5 (SA).** Mitigación interina:
   dar un intervalo centrado en el borde visible. Mejora barata pendiente: sembrar least_squares
   con un barrido grueso de B (global de pobre).
5. **Ajuste de UN borde sobre Co-60 (2 bordes) da un punto SESGADO** (queda entre los bordes,
   corrido al de mayor amplitud —963 keV—, con σ inflado). Punto sin sesgo → Etapa 6.
6. **El histograma debe extenderse por encima del borde**, o el flanco descendente se trunca y
   el ajuste falla. En datos reales suele haber eventos de mayor área (2º borde, pile-up,
   cósmicos) que dan aire; si no, subir `RANGO_MAX`.
7. **`RANGO_MAX=None` (auto) arruina el binning en datos reales.** Unos pocos eventos de área
   enorme (pile-up/cósmicos) estiran el `rango_max` a ~46000 mV·ns mientras toda la física vive
   entre 0 y ~6000: con 1024 bins eso deja el borde con ~130 bins y el flanco sub-muestreado.
   → Fijar `RANGO_MAX ≈ 8000` (deja aire sobre el borde, cumple §7.6). El hombro candidato a
   borde se ve entre ~2000 y ~5000 mV·ns en `espectro_neto.png`.
8. **El NCC (Etapa 4) también regala la SEMILLA del ajuste** [fuera del paper]: el intervalo
   ganador queda mapeado punto a punto contra la plantilla, así que de ahí salen B y σ. Sembrar
   `least_squares` con eso resuelve en la práctica la fragilidad de §7.4 (validado: error en B
   de -0.2 canales sobre 600, vs -26 con un intervalo malo elegido a mano). El SA (Etapa 5)
   sigue siendo la solución robusta definitiva, pero ya no es un bloqueador.
9. **El σ de la plantilla NCC casi no importa para B.** El paper no dice con qué σ construir S2
   (σ es justo lo que el SA ajusta). Barrido de `sigma_rel`: B apenas se mueve, pero ρ **sí**
   tiene un máximo claro, así que el propio NCC estima la resolución. Sintético: ρ = 0.960 →
   **0.9986** → 0.976 con el pico en σ_rel = 0.06 (verdadero 0.067).
10. **El barrido 2D del paper está SOBRE-PARAMETRIZADO respecto de su propio modelo.**
   [la desviación es NUESTRA: **fuera del paper**, aunque la justifique su Ec. 8] La Ec. 8
   usa `E_eff = E_e·B/E_max`: un mapa canal↔energía **por el origen** (área cero ↔ energía cero).
   Si el mapa pasa por el origen, fijado B la ventana YA está determinada: `[0.45·B, fin_rel·B]`.
   Barrer (inicio, fin) libres agrega un grado de libertad que el modelo no tiene. En datos
   reales eso se paga: el barrido libre eligió una ventana de 32 bins en la cola vacía con
   ρ=0.984, y hasta su candidato "razonable" (B≈3806) salía de una ventana de 360 mV·ns cuando
   la plantilla abarca 0.82·B ≈ 3100. → `ncc.buscar_intervalo(modo="origen")` barre **solo B**
   (1D). Es el default: más robusto, 7x más rápido, y sobre sintético da B exacto (600.0).
   El barrido del paper sigue disponible con `modo="libre"`.
11. **Pearson ignora la amplitud → la cola vacía gana si no se filtra.** Sobre datos casi nulos
   el z-score amplifica ruido puro hasta parecerse a cualquier forma. En `modo="libre"` hay que
   pasar `cuentas_crudas=` para activar el filtro `min_cuentas` (rechaza ventanas sin
   estadística). En `modo="origen"` el problema casi desaparece porque la ventana no puede
   encogerse: su ancho está amarrado a B.
12. **`ET.parse` es inviable con el XML de 1.76 GB** (~5·10⁷ tags `<Data>` → decenas de GB de
   RAM y minutos). `spectrum._iter_eventos` lo escanea como TEXTO (un tag por línea) y convierte
   cada pulso con `np.fromstring(",".join(pares), sep=",")`: **27 s** el Co-60, 5 s el fondo, con
   RAM plana. Además las áreas se cachean en `<archivo>.areas.npz` (clave = tamaño + mtime +
   umbral + baseline), así que las corridas siguientes son instantáneas.
13. ⚠️ **Superado por §7.33** (ese 0.18 era el área rota; con el área buena σ/B = 0.088).
   **σ_rel llega a 0.18 en el Co-60 real, y NO es (solo) el segundo borde.** ρ sube monótono
   0.69 (σ_rel=0.03) → **0.9918** (σ_rel=0.18) → baja. Predije que con dos bordes σ/B caería a
   ~0.10: **falso**. El modelo de dos bordes solo lo baja de 0.199 a 0.177 (42% FWHM). O sea
   el ensanchamiento es REAL, no un artefacto del solapamiento. Es esperable: la gaussiana del
   paper también se traga el **scattering Compton múltiple**, que en un centellador orgánico
   rellena la zona sobre el borde y no es gaussiano. Por eso la literatura de centelladores
   orgánicos suele usar función respuesta de Geant4. [fuera del paper]
14. **La convolución NO puede evaluarse solo dentro de la ventana de ajuste.** `np.convolve`
   rellena con ceros, así que el modelo pierde los eventos de más abajo del límite izquierdo
   que el broadening debería meter hacia adentro. El ajuste lo compensa **achicando σ**: en la
   prueba sintética de dos bordes daba σ un **18% bajo** y D un 3.9% corrido. → `forward_model`
   ahora tiene `pad=True` (default): extiende la grilla hasta el canal 0 y 4σ a la derecha,
   convoluciona, y recorta. Con eso σ pasa a -2.9% y D a +0.62%. Afecta también a la Etapa 3
   (el ajuste de un borde sembrado por NCC pasó a recuperar σ = 40.0 exacto).
15. **Erratum Ec. 12:** dice `ΔE = E2 - E1` "between gamma lines" (159.3 keV). Está mal: D
   convierte energía del ELECTRÓN a canal, y lo que se separa en canal son los BORDES, así que
   `ΔE = E_max,2 - E_max,1 = 154.7 keV`. Con el ΔE de los gammas, D sale ~3% corrido.
16. **En el Co-60 real, B1 es robusto pero D NO.** Barriendo variantes (6p / σ ligadas / origen
   ligado / 5 ventanas distintas): B1 ∈ [4153, 4198] (**±0.5%**), pero D ∈ [4.33, 5.05] (±8%).
   La causa: `B2 - B1 = 686` mV·ns mientras `σ1 = 741` mV·ns → los bordes están a **0.93 σ**.
   Por debajo de 1σ la descomposición interna no tiene información, aunque la posición del
   conjunto sí. → Para la calibración usar `ligar_origen=True` (D = B1/E_max,1), que da
   D ∈ [4.333, 4.344] en todas las variantes. Está justificado porque el test del origen
   sale limpio: offset en E=0 = -2.3% de B1.
17. **La simulacion** — `ComptonSpectrum_Co60_3M.root` (31 MB, TTree `ComptonSpectrum`,
   852079 eventos con Edep>0 de 3M generados). Ramas: `Edep_keV`, `Lvis_keVee` (luz con Birks
   ya aplicado), `E0_keV` (1173.2 con 436146 eventos y 1332.5 con 415933), `NCompton`, `NPhoto`,
   `FullAbs`, `QF`. Se lee con `uproot` (5.7.4 instalado). **25% de los eventos son scattering
   multiple** (145497 dobles, 35072 triples, ...), con Edep llegando hasta 1332.5 keV, muy por
   encima del borde de 1118. Quenching de electrones: `Lvis/Edep` = 0.863 a 50-100 keV y satura
   en **0.954** sobre 700 keV.
18. **El metodo del paper sobreestima D ~5% en este detector, y se puede demostrar.** Ajustando
   los MISMOS datos cambiando solo la plantilla: Geant4 completo (Edep) → D = 4.13; Geant4 con
   `NCompton==1` → 4.37; Klein-Nishina analitico → 4.34. Las dos ultimas coinciden porque
   Klein-Nishina **es** scattering simple. El scattering multiple mueve cuentas por encima del
   borde; un modelo que no lo tiene debe correr el borde hacia arriba en canal para explicarlas.
   → Con simulacion disponible, el forward folding le gana al metodo analitico. [fuera del paper]
   ⚠️ **Actualizado con el área corregida (§7.33):** el sesgo del método analítico es
   **+4.3%** (`fit_double` con dos bordes y origen ligado da 1.3157 contra 1.2618 del modelo
   extendido), y coincide casi exacto con lo que da fijar `w = 1` en `modelo_ext` (1.3048).
   O sea: el sesgo del método del paper ES el scattering múltiple, ni más ni menos.
   Y ahora los datos lo dicen solos: fijar `w = 1` empeora el ajuste en **Δχ² = +13**,
   mientras que fijar `w = 0.751` (la mezcla de la simulación) no lo empeora nada (Δχ² = 0).
19. **Ajustes SIN pesos = sin barras de error.** Todo hasta la Etapa 6 usaba mínimos cuadrados
   planos, que sobrepondera el plateau (muchas cuentas) y no permite cuantificar nada.
   `spectrum.espectro_neto` ahora devuelve `error` (Poisson propagado por la resta de fondo:
   `sqrt(N_m·k_m² + N_f·k_f²)`, con piso 1) y `ajustar_forward_folding` acepta `errores=` →
   χ²/ndf y covarianza. Resultado: χ²/ndf = 0.96 (597 gdl), `corr(D, alpha) = -0.70` (D y la
   resolución están correlacionadas, hay que reportarlas juntas). Presupuesto de error de D:
   estadística ±0.6%, ventana de ajuste ±1%, elección de plantilla ±2%.
   ⚠️ Dos correcciones de la Etapa 8: ese `error` es el χ² de **Neyman** y sesga D (§7.31), y
   el ±1% de ventana estaba subestimado porque la dependencia con el límite inferior es una
   deriva monótona, no dispersión (§7.32). El presupuesto vigente está en `calibracion.py`.
21. **El SA gana claro con 3 parámetros y se rompe con 6.** Un borde, arranque malo
   (B₀ = 0.70·B, σ₀ = 3σ): `least_squares` da B = 343 (error -257 sobre 600), el SA da
   **600.2** con cualquiera de los dos objetivos. Con los 6 parámetros crudos de la Ec. 11 no
   converge: el paso `0.1·δ_i` mezcla amplitudes (~5000) con D (~0.6) en el mismo vector.
   Proyectando A1, A2 fuera con nnls (VARPRO, 4 parámetros) + pulido local, encuentra el óptimo
   global en 2 de cada 3 semillas. **El SA es estocástico: usar `sa.sa_multiarranque`, una sola
   corrida no es un resultado.** El paso fijo del paper tampoco converge; `paso_adaptativo=True`
   (default, [fuera del paper]) escala el paso para mantener la aceptación cerca de 0.35.
22. **Dos constantes del paper no son transferibles: `T_init = 2000` y `|ΔE| < 0.001`.** Metropolis
   compara T con ΔE, así que T solo tiene sentido en las UNIDADES DEL COSTO. Medido: la
   calibración automática pide **T_init = 2.1e9** con la Ec. 14 (suma de cuadrados) y **3.19**
   con el error NCC. Nueve órdenes de magnitud. Lo mismo con el corte absoluto de 0.001. Por eso
   `recocido_simulado` usa `T_init="auto"` (calibra para ~80% de aceptación inicial) y tolerancia
   RELATIVA. Además, verificado: multiplicar A por 7.3 cambia el costo NCC en **exactamente 0**,
   así que con el objetivo NCC del texto de Metropolis las amplitudes son INDETERMINABLES; la
   Ec. 14 (mínimos cuadrados) es la lectura utilizable.
23. **El SA demuestra que Θ = {A1,A2,B1,σ1,σ2,D} está SUB-DETERMINADO para bordes solapados.**
   Es el hallazgo más importante de la Etapa 5, y llegó por hacer bien el trabajo del optimizador.
   Sobre los datos reales, cada libertad que el modelo tiene de más, el óptimo global la usa para
   algo sin sentido físico:
   - D libre (6p): χ²/ndf = **0.957**, mejor que el ajuste local (1.192)... con offset en E=0 de
     **-5928 mV·ns** y σ1/B1 = 0.44. Usa la componente de "963 keV" como fondo ancho.
   - origen ligado (5p): χ²/ndf = 0.959 con **D = 4.116** (el local sembrado daba 4.344, o sea
     el local NO estaba en el óptimo global)... pero σ1 = 1047 > σ2 = 601, invertido respecto de
     la física (σ crece con E).
   - origen + σ ligadas (4p): **A1 = 0**. Ajusta con una sola componente, cuando el Co-60 emite
     las dos líneas casi igual (la simulación: 436146 vs 415933 eventos, razón 0.95).
   → El ajuste local parecía sano solo porque la semilla lo mantenía en la zona física. La
   solución NO es cambiar de optimizador sino RESTRINGIR el modelo con física. Y la forma limpia
   de hacerlo es la plantilla Geant4 (§7.18), que fija de una sola vez las intensidades relativas
   de las dos líneas, el scattering múltiple y el quenching, dejando solo {A, D, alpha}.
   Consistencia: SA+origen da D = 4.12 y el forward folding 4.05-4.23; el local sembrado daba
   4.34. Las dos mejoras independientes empujan D en la misma dirección.
25. **Llevar la simulación a canal con σ→0 da un PEINE de picos, no el espectro ideal.**
   Es aliasing de binning: los bins de la simulación (2 keV ≈ 8.5 mV·ns con D≈4.2) no
   coinciden con los del histograma medido (10 mV·ns), así que unos canales acumulan dos
   bins y otros ninguno. Para dibujar "la simulación sin la resolución del detector" hay que
   ensanchar exactamente UN BIN: `modelo_simulado(..., alpha=0, beta=ancho_bin/D)`, que da
   σ_ch = ancho_bin. Con eso se ven los escalones limpios en ~3900 y ~4500 mV·ns.
24. **Tkinter en Windows sale borroso si no se declara DPI awareness.** Windows renderiza la
   ventana a 96 DPI y la ESTIRA como bitmap cuando la pantalla tiene escalado. Ninguna opción
   de matplotlib lo arregla: hay que llamar `ctypes.windll.shcore.SetProcessDpiAwareness(2)`
   **antes** de crear el root de Tk. Medido en este equipo: `winfo_fpixels('1i')` = **120**, o
   sea escalado al 125%. Después de activarlo hay que escalar a mano (`tk scaling` = dpi/72 y
   `figure.dpi` = 100·dpi/96), porque si no todo queda chico.
20. ⚠️ **DESMENTIDO por §7.33 — era instrumental.** Lo que sigue fue el razonamiento hasta
   el 25-ago-2026, y la conclusión ("revisar la ventana de integración y la colecta de luz")
   apuntaba al lugar correcto: el problema estaba en cómo se integraba el pulso, no en el
   detector. Con el área corregida la resolución es **20.6% FWHM**. Texto original:
   **La resolución de ~43% FWHM a 963 keV es REAL, no un artefacto del modelo.** Es alta hasta
   para plástico+SiPM. Pero las tres plantillas (Klein-Nishina, Geant4 completo, Geant4 de un
   solo scattering) piden 38-44%, así que no se puede echar la culpa al modelo. Antes de darla
   por buena conviene revisar la ventana de integración del pulso (205 ns a 5 GS/s: ¿entra toda
   la cola de recuperación del MicroFC-30020?) y la colecta de luz.
26. **El binning de la plantilla simulada es irrelevante, y eso es lo que hace posible la
   Etapa 8.** `modelo_simulado` cuesta 6.8 ms/eval con los 700 bins de 2 keV; el SA pide ~18000
   evaluaciones y el MCMC ~48000 por plantilla, o sea 11 minutos. Pero el kernel gaussiano tiene
   σ_E(963) = 178 keV, **90 veces más ancho que un bin**, así que rebinnear no cambia nada:

   | bins | ancho | D | α | χ²/ndf | eval |
   |------|-------|---|---|--------|------|
   | 700 | 2 keV | 4.2293 ± 0.0236 | 5.733 | 0.963 | 6.8 ms |
   | 350 | 4 keV | 4.2293 ± 0.0236 | 5.733 | 0.963 | 4.0 ms |
   | 175 | 8 keV | 4.2294 ± 0.0236 | 5.733 | 0.963 | 1.1 ms |
   | 70 | 20 keV | 4.2294 ± 0.0236 | 5.731 | 0.963 | 0.3 ms |

   → `bins=175` para todo lo iterativo (`Datos.simulacion(..., bins=175)` en la GUI). El default
   sigue en 700 por compatibilidad con las pestañas viejas.
27. **Interpolar a MÁS puntos en el NCC no agrega información; interpolar a POCOS infla ρ.**
   [**fuera del paper**: el paper fija 100 puntos, acá se usa un punto por bin]
   La ventana real tiene ~570 bins y el suavizado Savitzky-Golay de 21 bins deja sólo ~27 puntos
   INDEPENDIENTES, así que subir `n_puntos` sólo agrega puntos correlacionados:

   | n_puntos | ρ | B | D = B/963.4 |
   |----------|---|---|-------------|
   | 50 | 0.99402 | 4095 | 4.251 |
   | 100 (el default del paper) | 0.99268 | 4115 | 4.271 |
   | 200 | 0.99211 | 4130 | 4.287 |
   | 400 | 0.99193 | 4110 | 4.266 |
   | 571 (1 punto por bin) | 0.99191 | 4110 | 4.266 |

   B se mueve 0.4% en total y queda idéntico entre 400 y 571. Lo que sí pasa es que con pocos
   puntos ρ sale **inflado** (la interpolación promedia ruido), y ρ es justamente el criterio con
   que se comparan ventanas. → `n_puntos="auto"` = un punto por bin, sin inventar puntos.
   Restricción: n queda FIJO durante todo el barrido; si cambiara con B los ρ de ventanas
   distintas dejarían de ser comparables (misma trampa que §7.10).
28. **El NCC con plantilla SIMULADA funciona y da la semilla de D.** Cambiando sólo S2 (de
   Klein-Nishina a Geant4 llevado a energía con `modelo_simulado(D=1)`), el barrido de un
   parámetro da ρ = 0.992 y ventanas coherentes con lo que se elegía a mano:
   Lvis → [1850, 7856] con D₀ = 4.27; Edep → [1783, 7569] con D₀ = 4.11;
   Edep-1-scattering → [1893, 7280] con D₀ = 4.37. Los tres D₀ caen a ~1% del D final del ajuste.
   Con la plantilla simulada B ya no es "el borde" (hay dos): es el ancla del mapa por el
   origen, D = B/E_ref con E_ref = 963.4 keV. `fin_rel` también cambia: tiene que pasar el
   SEGUNDO borde (1118.1/963.4 + 3σ_rel), si no la ventana corta justo la parte que distingue
   una plantilla de otra.
29. **El α que prefiere el NCC (7.5) no es el α del χ² (5.9), y no es contradicción.** El NCC
   correlaciona forma sobre el espectro SUAVIZADO y SIN pesos, y además `fin_rel` depende de
   σ_rel, así que al barrer σ_rel está comparando ρ entre **ventanas distintas** — un confound
   del mismo tipo que §7.10. → el NCC entrega VENTANA y semilla de D; la resolución la
   determinan el χ² ponderado y el posterior. No usar el σ_rel del NCC como medida de resolución.
30. **El SA confirma que el óptimo de {A, D, α} es global.** Sobre las tres plantillas, partiendo
   sólo de la semilla del NCC (sin la grilla 40×24 que `ajustar_forward_folding` necesita para
   sembrarse), el SA llega a D = 4.2107 / 4.0224 / 4.3139 y el pulido local mueve el 4º decimal.
   O sea: acá el SA **no cambia el número**, y eso es el resultado — lo que sí cambiaba era el
   modelo sub-determinado de 6 parámetros de §7.23. Con VARPRO (A proyectada por mínimos
   cuadrados) el SA busca en 2 dimensiones y tarda ~50 s con 3 semillas.
31. **Estábamos usando el χ² de Neyman y eso sesga D hacia abajo.** [**fuera del paper**: la
   Ec. 14 es suma de cuadrados sin pesos] `espectro_neto` calcula la
   barra de error con las cuentas OBSERVADAS, así que un bin que fluctuó bajo recibe σ más
   chico → **más peso** → el ajuste se deja arrastrar hacia abajo. El sesgo va como 1/N por bin.
   Medido con verdad conocida (D = 4.20, plantilla Lvis):

   | cuentas/bin | Neyman (σ del dato) | Pearson (σ del modelo) |
   |-------------|---------------------|------------------------|
   | 2.6 | −1.92% | +0.67% |
   | 26 | −0.70% | +0.06% |
   | 258 | −0.05% | −0.02% |

   Y sobre **20 realizaciones** a 26 cuentas/bin: Neyman da sesgo −0.53% y cobertura del 68% de
   **8/20**; Pearson da −0.06% y **15/20** (lo que corresponde). Este experimento tiene ~15-25
   cuentas/bin en la ventana, o sea el sesgo era **comparable al error estadístico que
   reportábamos** (±0.6%). → `spectrum.varianza_predicha(fondo, k_m, k_f)` da
   `var = k_m·(modelo + fondo) + k_f·fondo`, y `sim_template.ajustar_pearson` itera 3 veces
   (la varianza depende de los parámetros). En el posterior hay que agregar el término
   `ln σ²(θ)`, que con varianza fija era constante y se podía ignorar.
   Sobre datos reales el efecto es +0.37% en D (4.2110 → 4.2264) y +8.7% en α.
32. **Hay un EXCESO de cuentas a bajo canal que la simulación no tiene.** ⚠️ Los números de
   canal de acá son de la escala vieja (§7.33): el exceso sigue existiendo pero ahora vive
   abajo de ~800 mV·ns, y desde la Etapa 9 se MODELA en vez de cortarse (§7.35), así que ya
   no manda el presupuesto de error. Lo que sigue vale como razonamiento, no como números.
   **Texto original:** Ajustando SÓLO arriba (ventana [3500, 7856], lejos de
   cualquier contaminación) y extrapolando el modelo hacia abajo, la razón dato/modelo es:

   | canal [mV·ns] | 1400 | 1800 | 2200 | 2600 | 3000 | 3400 | 3800 |
   |---------------|------|------|------|------|------|------|------|
   | dato/modelo | 1.41 | 1.17 | 1.10 | 0.97 | 0.95 | 1.00 | 0.99 |

   (y sigue creciendo hacia abajo: 4.8x en 1000-1200). Es un EXCESO, no un déficit, así que
   **no es el umbral de trigger**. Candidatos: dispersión en el material de alrededor y en el
   portafuente, que la simulación no modela; dark counts / afterpulsing del SiPM apilados.
   El fondo ya está restado, así que no es eso.

   Consecuencia directa: si la ventana de ajuste entra en esa zona, el modelo absorbe el exceso
   **ensanchando α y bajando D**. Barrido del límite INFERIOR (el superior no importa: D queda
   igual a 4 decimales entre hi = 7000 y 8000):

   | lo | 1600 | 1850 | 2000 | 2200 | 2500 | 2800 |
   |----|------|------|------|------|------|------|
   | D | 4.190 | 4.226 | 4.247 | 4.268 | 4.312 | 4.336 |
   | α | 6.67 | 6.38 | 6.22 | 6.05 | 5.72 | 5.55 |
   | FWHM(963) | 50.6% | 48.4% | 47.2% | 45.9% | 43.4% | 42.1% |

   Deriva MONÓTONA de 3.5%, no dispersión. Por eso el ±1% de ventana de §7.19 estaba
   subestimado: se había barrido 1800-2500 leyéndolo como ruido, cuando es una tendencia.
   → **La ventana tiene que arrancar arriba de ~2600 mV·ns** (inicio_rel ≈ 0.62 en vez del
   0.45 del paper).

   ⚠️ **Esto es un límite real del NCC**: compara FORMA, y no tiene manera de saber que a bajo
   canal hay física que la plantilla no contiene. Eligió inicio = 0.45·B = 1850 y se comió la
   zona sucia. La selección automática de ventana resuelve la subjetividad, NO la validez del
   modelo. `inicio_rel` es ahora un control de la pestaña 8 para poder verlo en vivo.

33. **EL AREA DEL PULSO ESTABA MAL CALCULADA, y era el error dominante de todo el
   proyecto.** `_area_pulso` promediaba la línea base sobre el **primer 20% del registro**
   (204 muestras = 40.4 ns) e integraba los 200 ns completos. Pero el pulso **arranca a los
   ~25 ns y tiene el pico a los 48**: la ventana de línea base se comía la subida del propio
   pulso. El sesgo `b → b + q_b/T_b` se multiplica por `T/T_b ≈ 5` al integrar todo el
   registro, así que:

   | integral/amplitud [ns] | mediana | std | eventos con signo invertido |
   |---|---|---|---|
   | línea base 0-12 ns (limpia) | −8.48 | 0.84 | 0 / 3639 |
   | línea base 0-40 ns (la vieja) | −7.71 | **11.33** | **954 / 3639 (26%)** |

   y el `abs()` reciclaba esos 954 como si fueran buenos. Peor todavía: el jitter del
   trigger mueve el pico entre 18 y 54 ns, así que **cada evento perdía una fracción
   distinta de su carga** → ensanchamiento multiplicativo puro. Y el fondo, tomado otro día
   con el pulso a 120 ns, tenía la línea base LIMPIA: muestra y fondo no medían lo mismo.
   → `spectrum._area_pulso_gate`: onset por evento (cruce del 10% retrocediendo desde el
   pico), línea base = mediana hasta 4 ns antes del onset, integral en
   `[onset−8, onset+55] ns` (`spectrum.GATE_PRE`/`GATE_POST`; acá decía −6, que
   nunca fue el valor del código). Verificado: el pulso medio devuelve el **100% de su área en
   80 ns** y `área/amplitud` es constante al 1% entre 20 y 220 mV, o sea la compuerta no
   corta nada y es lineal. Banderas de calidad: onset temprano (8% del Co-60), saturación
   del DRS4 a 500 mV (0.02% del Co-60, 2.3% del fondo), compuerta fuera del registro.
   **Consecuencias:** D pasó de 4.32 a 1.276 mV·ns/keVee, la resolución de 42.7% a 20.6%
   FWHM (§7.20 tenía razón: el ancho era instrumental), y el espectro pasó a tener un borde
   limpio en vez de una caída de 2000 mV·ns de ancho.

34. **La corrida de Co-60 estuvo VIVA solo el 40% del tiempo, y el fondo el 100%.**
   [**fuera del paper**: la Fig. 2 arranca en `(canal, cuentas)`, el tiempo vivo es aguas
   arriba] El
   `<Scaler1>` de cada evento es la tasa de trigger REAL de la placa: mediana **170 Hz** en
   el Co-60 contra 71.5 ev/s efectivamente guardados (el DRSOSC vuelca ~35 kB de XML por
   evento y no da abasto), y **10 Hz** contra 10.1 ev/s en el fondo. O sea T_vivo = 283 s de
   699 s de pared, y 986 s de 986. Normalizando por tiempo de PARED como se hacía, el fondo
   entraba con k = 0.709 cuando corresponde **k = 0.287**: sobre-resta de 2.5x.
   → `spectrum.tiempo_vivo(relojes, scalers, duracion, n_eventos)`, usado por defecto
   (`usar_tiempo_vivo=True`). Ojo: el scaler del fondo tiene ráfagas de hasta 90 kHz (ruido
   electrónico), por eso se recorta a 10x la mediana antes de integrar. Las marcas de
   `<Time>` NO sirven para esto: el DRSOSC las escribe en ráfagas (4124 eventos del Co-60
   comparten timestamp al ms).

35. **La mezcla de scattering y el exceso a bajo canal dejaron de ser sistemáticos "a mano"
   y pasaron a ser parámetros.** [**fuera del paper**: el modelo extendido no es de Kong
   et al., que usa la lineshape analítica sola] La Etapa 8 reportaba ±2.7% de sistemático comparando tres
   plantillas (Geant4 completo / 1 scattering / Klein-Nishina) y cortaba la ventana por
   arriba del exceso. `modelo_ext.py` mete las dos cosas adentro:
   `S = A·[w·h_1sc + (1−w)·h_multi] ⊗ G(σ_E) + C·exp(−ch/λ)`.
   Sobre los datos reales, ventana [300, 1900]: **w = 0.796 ± 0.070** contra el 0.751 de la
   simulación, o sea **los datos y la simulación coinciden a 0.6σ**: la geometría simulada
   reproduce bien cuánto scattering múltiple hay. χ²/ndf = 1.008 con 314 gdl,
   `corr(D, w) = +0.96`: la ignorancia sobre la plantilla ya está DENTRO del error de D, que
   es lo que se buscaba.
   Prueba de que w está realmente medido y no es decorativo: fijarlo en 0.751 no cambia el
   χ² (Δχ² = 0, D = 1.2537), y fijarlo en 1 —que es lo que asume la plantilla analítica del
   paper— lo empeora en Δχ² = +13 y sube D a 1.3048.
   Controles hechos:
   - **Cierre con verdad conocida** (20 realizaciones Poisson a la estadística real):
     sesgo en D **−0.10%**, rms 0.0116 contra error reportado 0.0127, cobertura 17/20.
   - **Barrido de la definición del área** (compuerta 45/55/70/90 ns, onset al 5/10/20%,
     guarda, con y sin descarte por calidad): D entre 1.2536 y 1.2799, o sea **±1.04%**. La
     de 45 ns es la única con χ² claramente peor (1.076): trunca el pulso.
   - **Ventana**: D varía 1.5% entre lo = 250 y lo = 800, **sin deriva monótona** (§7.32
     tenía 3.5% de deriva). El límite superior no mueve nada (0.1% entre 1600 y 2200).
   - **Ajustando SOLO la zona limpia [900, 1900] y sin componente de exceso**, con w libre:
     D = 1.2748 ± 0.0297 contra 1.2765 del modelo completo. O sea el exceso no arrastra a D.
   - Los términos `a·E` y `β` de la resolución se van a **cero**: con dos bordes a 963 y
     1118 keV no hay palanca para más que la √E.
   Lo que queda como sistemático dominante es la **FORMA** del exceso: cambiarla de
   exponencial a ley de potencia mueve D **+2.35%**, y las dos describen los datos igual de
   bien (χ²/ndf 1.008 contra 1.015). Eso no se cierra ajustando mejor: hay que medir el
   exceso (corrida con la fuente apantallada o corrida de lugar).

36. **El canal del DRSOSC no es una constante, y leerlo mal da CERO pulsos en silencio.**
   Las corridas de junio (Co-60, SiPM) salieron por **`CHN2` con `<Scaler1>`**; las de la
   sesión del **26-ago-2026** (Cs-137, "Detector A" a 1600 V, en el disco externo
   `D:\Sesión Calibración 26_08_2026\`) salieron por **`CHN1` con `<Scaler0>`**.
   `_iter_eventos` tenía `CHN2` cableado, así que recorría los 2.9 GB del Cs-137 sin
   encontrar un solo `<Data>` adentro del tag que buscaba: 80322 eventos leídos, **0 pulsos**,
   y el caché `.areas_gate.npz` guardaba ese vacío para siempre. Lo que se veía aguas abajo
   era un histograma de ceros o un `IndexError` cinco etapas más allá, nunca el problema.
   → Tres arreglos: `detectar_canal(ruta)` mira el primer evento y devuelve
   `(CHN{n}, Scaler{n-1})`; el canal entra en la **clave del caché** (si no, el `.npz` vacío
   se reusa aunque se arregle el código); y `leer_areas_y_tiempo` **levanta `ValueError`** si
   quedan 0 pulsos útiles, en vez de cachear el vacío. La pestaña Archivos muestra el canal
   detectado de cada XML al lado del archivo.
   Verificado: las áreas del Co-60 y del fondo viejos salen **bit a bit idénticas** después
   del refactor, y `calibrar_v2.py --rapido` reproduce D_luz = 1.2618 ± 0.0133.
   ⚠️ Como la clave del caché cambió de largo, la primera corrida re-parsea los XML una vez.

37. **La corrida de Cs-137 del 26-ago-2026 estuvo viva el 0.13%.** El `<Scaler0>` marca
   **46 kHz** de trigger real contra 58.5 ev/s guardados: T_vivo = **1.7 s de 1373 s** de
   pared. El fondo de esa sesión (2669 s) está 100% vivo, así que entra con
   **k_f = 0.00064** y es despreciable. El número es enorme y conviene confirmarlo con el
   umbral de trigger de esa corrida antes de usarlo para nada absoluto (tasas, eficiencias);
   para la FORMA del espectro —que es lo que calibra— no cambia nada, porque el fondo pesa
   casi cero de todos modos. Ojo también que el fondo de esa sesión tiene pulsos mucho más
   grandes que la muestra (p99 = 10064 contra 1022 mV·ns): sin fuente, el trigger lo disparan
   cósmicos.

38. **Primer número del Detector A (Cs-137, 1600 V), por la ruta analítica.** Con el canal
   arreglado, el espectro neto muestra un Compton limpio: plateau hasta ~780 mV·ns y caída
   que muere en ~1100. NCC con Klein-Nishina a E_γ = 661.657 keV (E_max = **477.3 keV**) da
   ventana [386, 1322], **ρ = 0.9991**, y el ajuste sembrado
   **B = 848 mV·ns → D ≈ 1.78 mV·ns/keV**, con σ/B = 11.9% (**28% FWHM** en el borde).
   Por la pestaña 6 de la GUI (binning de 5 mV·ns, ventana NCC [387, 1323]):
   **D = 1.7964**, B = 857.5, ρ = 0.9992, χ²/ndf = 1.907. El χ² alto es
   esperable: la ventana entra en el exceso de bajo canal, que la plantilla
   analítica no tiene (§7.35). El método clásico del 50% da 1.8802, **+4.7%**.
   Es PRELIMINAR por dos razones conocidas: la plantilla analítica sesga D **+4.3%** hacia
   arriba porque no tiene scattering múltiple (§7.18), y este ajuste no modela el exceso de
   bajo canal (§7.35). Para cerrarlo hace falta **una simulación Geant4 de Cs-137** con la
   geometría de este detector; el `.root` que hay es de Co-60 y no sirve como plantilla acá.
   Tampoco es el mismo detector que el de las Etapas 0-9 (aquel era el SiPM MicroFC-30020,
   este es "Detector A" a 1600 V), así que **D = 1.78 no contradice a D = 1.26**: son
   detectores distintos.

39. **La fuente dejó de estar cableada: `physics.Fuente` y el selector de la GUI.**
   Hasta acá `EG1, EG2 = 1173.2, 1332.5` vivían en `gui.py` y el Co-60 estaba
   metido en cada pestaña. Ahora la fuente es un objeto (`FUENTES` trae Co-60 y
   Cs-137; el combo "manual" acepta una o dos energías tipeadas) que se elige en
   la pestaña **Archivos**, se guarda en `.cache/archivos.json`, y de él salen
   `lineas`, `bordes`, `E_ref` (el ancla del mapa por el origen) y `doble`.
   Lo que cambia según la fuente:
   - **Pestaña 6** ajusta dos bordes solapados (Ec. 11-13) si hay dos líneas, y
     un solo borde (Ec. 8-10) si hay una; la solapa dice cuál está activo. Con
     una línea D = B/E_max sale directo y se muestra al lado el método clásico
     del 50% para comparar (en el Cs-137 medido: +4.7% de sesgo).
   - **NCC** (4 y 8) usa `E_gamma`, `E_ref` y `E_ref2` de la fuente. Con una
     sola línea no hay segundo borde que pasar, así que la ventana se cierra con
     el criterio de Klein-Nishina (borde + 4.5σ) en vez del de §7.28.
   - Las plantillas simuladas y la resolución se reportan en `E_ref`, no en un
     963.4 escrito a mano (`sigma_E_ref`, `fwhm_rel_ref`).
   **La simulación pasó a ser OPCIONAL**: sin `.root` la GUI hace la ruta
   analítica completa (pestañas 4 y 6) y las tres pestañas que la necesitan
   (7, 8, 10) muestran un aviso limpio en vez de reventar (`SinSimulacion`).
   Verificado que no hay regresión: con Co-60 + simulación, la pestaña 10 sigue
   dando **D = 1.2618 ± 0.0137** con w = 0.797, y `calibrar_v2.py --rapido`
   reproduce 1.2618 ± 0.0133.

40. **Dos números de la escala vieja seguían escondidos en la GUI, y uno rompía
   un ajuste.** Aparecieron al generalizar las ventanas:
   - `TabDosBordes` pasaba `B1_rango=(0.60*4357, 1.05*4357)`, o sea el barrido
     de semilla buscaba B1 entre 2614 y 4575 mV·ns cuando el borde real está en
     1229 (§7.33). Resultado: la pestaña 6 venía mostrando **D = 1.90** en vez
     de 1.32 desde el arreglo del área. Ahora el rango sale del NCC.
   - La pestaña 7 dibujaba `z = (can >= 1200) & (can <= 7500)`: con el área
     nueva, el panel izquierdo quedaba casi vacío. Ahora la ventana de dibujo se
     deriva de la de ajuste.
   Moraleja repetida (§9-bis.5): cuando cambia el OBSERVABLE hay que barrer
   todas las constantes de canal del proyecto, no sólo las del análisis.
   Por eso ahora las ventanas se piden a **`Datos.ventana_ajuste`** y a
   **`Datos.borde_estimado`** (NCC), que las ponen en escala solas: con el
   Co-60 dan (307, 1905) y (922, 1905), o sea los mismos 300-1900 y 900-1900 que
   estaban tipeados, pero calculados.

41. **La calibración entera se configura desde la fuente, y el último número
   escrito a mano era la pestaña del resultado** (27-ago-2026). La §7.39 había
   hecho configurable el MODELO (una línea o dos, E_ref, ventanas del NCC),
   pero la pestaña ★ Calibración seguía siendo un cartel con los números del
   Co-60 leídos de `calibracion.py`: con Cs-137 elegido mostraba un aviso y
   los mismos 1.261. Ahora **calcula sobre lo que está cargado**. Tres cosas
   que aparecieron al hacerlo:
   - **La plantilla Geant4 ES la fuente, y nadie lo estaba verificando.** Con
     Cs-137 elegido y el `.root` de Co-60 puesto, las pestañas 7, 8 y 10
     devolvían un D perfectamente formateado y sin sentido: la plantilla tiene
     sus bordes en 963 y 1118 keV y los datos uno solo en 477. No hay error en
     ningún lado, sólo un número equivocado. → `sim_template.lineas_primarias`
     lee la rama `E0_keV` (cuesta 0.2 s) y `Datos.exigir_sim_compatible` corta
     con el mismo aviso limpio de `SinSimulacion`; se puede destrabar a mano
     con «usar la simulación igual» en la pestaña Archivos. Es la misma clase
     de falla silenciosa que §7.36 (el `CHN2` cableado): lo que no se verifica,
     se cachea.
   - **El presupuesto de error ahora se MIDE en vivo en vez de citarse.** La
     pestaña corre las variantes sobre los datos cargados (barrido del límite
     inferior de la ventana, exceso exponencial contra ley de potencia, w
     fijo en 1 y en el de la simulación) y suma. Sobre el Co-60 reproduce la
     Etapa 9 casi exacto: estadística 1.08% (era 1.07), ventana 0.91% (0.75),
     forma del exceso 2.20% (2.35), total **2.8%** (2.9). El único término que
     se ARRASTRA es la definición del área (±1.04%, `ERR_AREA_REL`), porque
     cada variante obliga a re-parsear los XML.
   - **Un registro por fuente, `.cache/calibraciones.json`.** Una recta de
     calibración se hace con varias fuentes medidas en corridas distintas, así
     que el resultado de cada una se guarda con la clave `fuente|huella` y el
     panel derecho las superpone. Ojo con qué prueba eso: cada fuente ya
     entregó su D con el mapa por el origen, así que su punto cae exacto sobre
     su propia recta. **Lo que hay que mirar es la dispersión ENTRE fuentes**,
     y eso —no el ajuste conjunto— es el test de linealidad que sigue
     pendiente (§9-bis.3).
   De paso: `physics.FUENTES` pasó de 2 a 9 presets (Am-241 … Tl-208) con
   `NOTAS_FUENTE`, y `Fuente.avisos()` deriva las advertencias de las líneas
   en vez de tabularlas, así que valen también para la fuente "manual" (borde
   por debajo de 100 keV, más de dos líneas, bordes que se solapan). Y la
   receta de la Etapa 9 salió de la pestaña 10 a la función `ajuste_extendido`
   de `gui.py`, para que la pestaña de Calibración no tenga una segunda copia
   que pueda divergir. Verificado que no hay regresión: la 10 sigue dando
   **D = 1.2618 ± 0.0137** con w = 0.797 y χ²/ndf = 0.998, y
   `calibrar_v2.py --rapido` no se movió.

42. **Barra de progreso: la fracción se mide en BYTES, no en eventos** (27-ago-2026).
   El paso que parecía un cuelgue era siempre el mismo: `_iter_eventos`
   recorriendo 1.76 GB durante ~40 s sin decir nada. Ahora acepta
   `callback(fraccion, n_eventos)` y la cadena llega entera hasta la barra:
   `spectrum._iter_eventos` → `leer_areas_y_tiempo` → `espectro_neto` →
   `Datos.espectro` (que usa el gancho `Datos.avisar`, puesto por
   `Pestana.lanzar`) → cola → `App.progreso_actualizar`. Tres cosas que hubo
   que resolver:
   - **En bytes y no en eventos**, porque cuántos eventos tiene el XML no se
     sabe hasta terminarlo, mientras que `os.path.getsize` está desde el
     principio. `f.tell()` no sirve: en un archivo de texto iterado con `for`
     Python lo prohíbe. Se acumula `len(linea)`.
   - **CRLF.** El archivo se abre en modo texto, así que Python traduce `\r\n`
     a `\n` y `len(linea)` queda un byte corto POR LÍNEA. Con ~5·10⁷ líneas eso
     es enorme: la barra llegaría al 78% y ahí terminaría. Se lee `f.newlines`
     en el primer evento y se compensa. Verificado: la fracción cierra en
     100.0% y es lineal en el tiempo (19.8 / 39.7 / 59.8 / 79.9 / 100.0%).
   - **Los dos XML se reparten la barra por su TAMAÑO**, no mitad y mitad: la
     muestra pesa 5 veces el fondo (corte en 0.833), y con el reparto ingenuo
     el ETA salía el doble de lo que era.
   El resto de las pestañas reporta fracciones de grano grueso: la de
   Calibración cuenta los ajustes que va a correr ANTES de empezar (`[3/10]`,
   que es lo que hace posible un ETA honesto), y los MCMC de las pestañas 8 y
   10 usan el `callback(k, nk, paso, npas)` que `mcmc.metropolis` ya emitía 20
   veces por cadena. Cuando un paso no sabe cuánto falta manda `frac=None` y la
   barra queda indeterminada, que para lo que se pide —ver que el programa no
   se colgó— ya alcanza. `Pestana.fase(a, b)` mapea el 0-100% de un sub-paso
   al tramo que le toca, para que el parseo del XML no haga saltar la barra del
   resto de la tarea.

43. **Las carpetas, y por qué el código dejó de depender del directorio de
   trabajo** (27-ago-2026). La raíz tenía 46 archivos mezclados: 2 GB de XML,
   18 PNG, los PDF y 20 `.py`. Ahora hay `codigo/` (con `etapas/`, `datos_io/`,
   `ajuste/`, `scripts/`), `datos/`, `figuras/`, `docs/` y `verificaciones/`.
   Tres cosas que hubo que resolver y conviene no re-romper:
   - **`codigo/rutas.py` es el único que sabe dónde está cada cosa.** Resuelve
     todo desde `__file__`, así que ya no importa desde qué carpeta se ejecute
     (antes los `.xml` se nombraban relativos al cwd y sólo funcionaba estando
     parado en la carpeta del proyecto). Las figuras de las demos van a
     `figuras/`, los cachés a `.cache/`, y el PDF a la raíz.
   - **Cada módulo lleva 5 líneas de bootstrap** que agregan `codigo/` a
     `sys.path`. Hacen falta porque `python codigo/etapas/ncc.py` pone en
     `sys.path` la carpeta *del archivo*, no `codigo/`, y entonces
     `from etapas.physics import ...` no resolvería. Correr cada módulo solo es
     una característica del proyecto (§4), así que esto no es opcional.
   - **Los cachés de áreas siguen valiendo** porque `mv` conserva tamaño y
     mtime, que es de lo que está hecha la clave (`spectrum._clave_cache`), y
     los `.npz` se movieron junto a sus `.xml`. Lo mismo las cadenas de MCMC,
     cuya clave es `Datos.huella()` (nombre + tamaño + fecha, no la ruta).
     Lo único que se invalidó a propósito es `.cache/archivos.json`, que
     guardaba rutas viejas; ahora se guardan **relativas a la raíz**, así que
     mover la carpeta entera de lugar ya no las rompe.
   Verificado sin regresión: `calibrar_v2.py --rapido` da D_luz = 1.2618 ±
   0.0133 igual que antes, las demos de los 9 módulos corren, la GUI abre, y
   `documentacion.py` regenera el PDF — todo también desde un cwd de afuera
   del proyecto.

44. **Los tres métodos clásicos, y que los tres sobrestiman** (27-ago-2026).
   El paper compara NCC-SA contra dos métodos convencionales y esa comparación
   era la mitad que faltaba. Están en `etapas/clasicos.py`, corren en la
   pestaña 6 y en la de Calibración, y `python codigo/etapas/clasicos.py` los
   mide contra un borde sintético de posición CONOCIDA:

   | σ/B | 50% de altura | diferencia | FWHM gaussiano | diferencia SIN ruido |
   |-----|---------------|------------|----------------|----------------------|
   | 0.03 | +0.8% | +0.7% | +2.6% | +0.4% |
   | 0.06 | +2.2% | +1.5% | +4.7% | +1.4% |
   | 0.09 | +3.8% | +2.8% | +6.3% | +2.6% |
   | 0.12 | +5.4% | +4.4% | +7.5% | +3.8% |
   | 0.18 | +8.5% | +7.4% | +9.2% | +6.2% |

   **El sesgo es del MÉTODO, no de la estadística**: la última columna es el
   mismo método sobre el espectro sin ruido y va igual de +0.4% a +6.2%. La
   causa es común a los tres: la lineshape de Klein-Nishina no es un escalón,
   tiene una acumulación creciente antes del borde, y ni el 50% de altura, ni
   el mínimo de la derivada, ni un escalón con fondo lineal la contienen. Es la
   justificación cuantitativa de ajustar la lineshape entera.
   Sobre el Co-60 real: 50% → +11.3%, FWHM gaussiano → +12.5%, diferencia →
   +0.3% contra el modelo extendido. Ese +0.3% es en parte suerte (con dos
   bordes solapados el mínimo de la derivada cae entre los dos); el sintético
   dice que a esta resolución le corresponde ~+2.8%.
   Dos trampas que costaron encontrar:
   - **La ventana de suavizado no puede ser una constante.** La derivada del
     borde es un pozo cuyo ancho escala con la resolución; con 21 bins fijos, a
     σ/B = 0.09 aparecen tres mínimos locales degenerados (−29.1, −28.2, −28.1)
     entre los canales 580 y 619 y el ruido elige cuál gana: el sesgo salta de
     +1.5% a −3.3%. `_ventana_auto` la deriva de la caída 80%→20% del flanco.
     Ojo: medir ese ancho sobre la DERIVADA no sirve (con ruido el mínimo cae
     en una espiga y la ventana se achica justo cuando había que agrandarla);
     hay que medirlo sobre el espectro.
   - **La región NO es la ventana de ajuste.** El modelo extendido ajusta desde
     ~0.25·B porque CONTIENE el exceso de bajo canal; los clásicos no, y la
     caída del exceso es mucho más empinada que el borde. Pasándoles la ventana
     de ajuste (317-1963 mV·ns) el mínimo de la derivada se va al exceso y da
     **D = 0.363 en vez de 1.265, −71%**. Por eso `comparar` toma `B_ref` (el
     del NCC) y arma la región alrededor del borde.
   El factor 0.75 de `_ventana_auto` se eligió barriéndolo contra el caso sin
   ruido, o sea con el criterio de **reproducir el método**, no de minimizar el
   error: bajar el sesgo con un suavizado que casualmente cancela sería hacer
   trampa en la comparación.

45. **El "1550" del Co-60 es el ALTO VOLTAJE, no el trigger, y eso rompe la
   comparación entre corridas** (28-ago-2026, confirmado por el usuario). Tres
   consecuencias, todas grandes:
   - **El fotodetector NO puede ser el SiPM MicroFC-30020** que dice §2: un
     MicroFC se polariza a ~25 V. A 1550-1700 V es un PMT. Eso invalida además
     la lectura de la Tabla I del paper que hacía §2 (`σ_gain = ENF·√N_pe` es un
     término de SiPM). Queda como dato a confirmar en el laboratorio; el
     catálogo lo dice explícitamente en la nota de la corrida.
   - **Co-60 (1550 V) y Cs-137 (1600 V) están a ganancias distintas**, así que
     no son dos puntos de una misma recta. Si lo único que hubiera cambiado
     fuera el HV, el cociente de sus D daría el exponente de la ley de ganancia:
     `b = ln(1.9450/1.2512)/ln(1600/1550) = 13.9`, **fuera del rango 6-10 de un
     PMT**. O sea que entre junio y agosto cambió algo más que 50 V.
   - **El atenuador.** El área se integra a la ENTRADA del digitalizador, así
     que 20 dB dividen D por 10, y dos corridas con atenuación distinta no son
     comparables ni por asomo. Cuando se escribió esto el dato no constaba: el
     Co-60 de junio decía "att" en el nombre y de las de agosto no había
     registro. ✔ **RESUELTO (22-sep-2026): son 20 dB en todas** (§7.70), así que
     el riesgo era real pero no se materializó. Ni `codigo/` ni `analisis/`
     aplican el factor, y está bien: el observable es lo que entra al DRS4.
   → `rutas.Corrida` tiene ahora `hv`, `atenuador`, `umbral`, `muestreo` y
   `distancia`; `Corrida.punto()` es la etiqueta por la que agrupa el registro
   (detector + HV), y `rutas.misma_cadena` avisa cuando dos corridas que se van a
   combinar difieren en la cadena de señal **o cuando el dato no está
   registrado**, que no es lo mismo que ser igual.

46. **Las dos linealidades son cosas distintas y hay que decir cuál se está
   midiendo** (28-ago-2026). Vienen juntas en la pregunta "¿es lineal?":
   - **Test A, el eje** (`canal = D·E`): valida la calibración. Necesita varias
     FUENTES al mismo punto de operación. Con dos fuentes hay un solo grado de
     libertad (la ordenada); con cuatro hay residuos de verdad.
   - **Test B, la ganancia** (`D(V)`): caracteriza el fotodetector. Necesita la
     misma FUENTE a varios HV, y **la respuesta esperada no es una recta**: un
     PMT tiene `D ∝ V^b` con b entre 6 y 10.
   Dos cosas cuantitativas que salieron de plantearlo así:
   - **Con los tres voltajes que hay, el test B no se puede cerrar.** La
     separación entre la recta y la potencia, medida como desvío de la cuerda en
     el punto medio, vale `(1/8)·b(b−1)·(ΔV/V̄)²`: 2.6% entre 1600 y 1700 (o sea
     del tamaño del sistemático), 11% si se baja a 1500, 18% si se baja a 1450.
     `voltaje.py` lo verifica con verdad conocida: con seis voltajes de 1450 a
     1700 recupera b = 7.94 ± 0.07 sobre 8.00 y rechaza la recta (χ² 10.6 contra
     327); con sólo los tres de arriba avisa que el span no alcanza.
   - **El test A va en keVee, no en keV.** El quenching de Birks vale 0.9541 en
     el borde del Co-60 y **0.9349** en el del Cs-137: 2.0% de diferencia. En
     escala Edep dos fuentes tienen que dar D distintas aunque la electrónica sea
     perfectamente lineal; en Lvis no, porque la plantilla ya trae el quenching.
     Hacer el test en keV sería medir a Birks y llamarlo no linealidad.
   Y una que conviene no copiar del paper: Kong et al. validan con *cosine
   similarity* 99.999% entre dos relaciones energía-canal. Entre vectores casi
   proporcionales esa métrica es insensible por construcción (dos calibraciones
   que difieren 5% en pendiente siguen dando ~99.99%): hay que reportar residuos
   en keV y la ordenada con su error.

47. **En `analisis/` la ventana estaba comiendo el exceso, y eso explicaba el 3%
   de desacuerdo con `codigo/`** (28-ago-2026). `analisis/` daba D = 1.2230
   contra 1.2618, ocho veces el error estadístico. Barriendo `inicio_rel` sobre
   el Co-60 aparece una deriva **monótona** (1.2230 → 1.2430 entre 0.45 y 0.80)
   con el χ²/ndf bajando de 1.576 a 1.067: es §7.32 otra vez, en la
   implementación nueva. Se cerró portando el término de exceso:
   `S = A·[P(ch/D) ⊗ G(σ)] + C·exp(−(ch−lo)/λ)`, con λ como tercer parámetro no
   lineal, (A, C) proyectadas por **nnls** (VARPRO de dos columnas) y la ventana
   de ajuste extendida a 0.25·B. Resultado: **D = 1.2512**, sin deriva con la
   ventana (1.2512 con 0.25, 1.2524 con 0.35), y el desacuerdo con `codigo/` baja
   a **0.2%** contra su variante de `w` fija en el 0.751 de la simulación
   (1.2537), que es lo que este modelo hace. El resto del 3% era la libertad de
   `w`. La receta queda escrita en `calibracion.EXCESO` / `INICIO_AJUSTE_REL`,
   con la tabla del barrido al lado.
   **Presupuesto de error medido** (`analisis/sistematicos.py`, Co-60, keVee):
   estadística 0.38%, forma del exceso 0.56% (punto a punto), gate de
   integración 0.39% (común), ventana de ajuste 0.05% ← se volvió despreciable
   justo por modelar el exceso. **Total 0.78%.** Dos variantes se miden y NO se
   suman, por razones distintas: la plantilla de solo-1-scattering porque los
   datos la rechazan (χ²/ndf 1.180 contra 1.169, y D = 1.3025, o sea el +4.1%
   del método analítico del paper otra vez), y el método de área por
   convergencia porque se descarta **por física** —su ventana de integración va
   de 13 ns en los pulsos chicos a 42 en los grandes— y no por χ², que no los
   distingue (1.223 contra 1.169 con D a 0.9%).
   ⚠️ El barrido del gate de acá mueve posición y ancho con el onset fijo, y da
   0.39%; el de `codigo/` (§7.35) barría además la fracción del onset y los
   cortes de calidad y daba 1.04%. El de acá es un subconjunto, no una mejora.

48. **El umbral del trigger es un piso duro para cualquier ventana, y no
   respetarlo no sesga: rompe** (28-ago-2026). El Cs-137 a 1600 V no tiene ni
   una cuenta debajo de 121 mV·ns y sube hasta ~334: es la curva de eficiencia
   del discriminador, que ni la plantilla ni el exceso contienen. Arrancando la
   ventana en 229 el χ²/ndf salta a **4.26** y la matriz de errores sale
   singular; en 320 da 0.95. → `espectro.fin_del_umbral` lo detecta como el
   primer máximo local del espectro suavizado (la eficiencia CRECE y todo lo que
   hay debajo del borde DECRECE, así que donde deja de crecer una cosa le pasa el
   mando a la otra) y `modelo.armar` lo impone como piso. En el Co-60, tomado con
   el umbral en 0.2 mV, cae en 56 mV·ns y no restringe nada, que es lo correcto.

49. **Una barra de error puede estar mal por el paso de derivación, y el síntoma
   es que no coincide con el cociente de verosimilitud** (28-ago-2026).
   `sa.errores_por_curvatura` usaba `paso = 1e-3·|x|` para todos los parámetros.
   Para la ordenada c0 del mapa eso da 0.16 mV·ns, un 3% de un bin de 5.4, y como
   la plantilla se lleva a canal con `np.interp` —lineal a trozos sobre la grilla
   del histograma— lo que se mide ahí es el quiebre entre dos trozos y no la
   curvatura del χ². Daba **21σ donde Δχ² da 2.1σ**. Ahora el paso se pide por
   parámetro (`Problema.pasos_curvatura`; para c0, un bin). Control: para D el
   paso no cambia nada (err = 0.0040 con `paso_rel` 1e-2 y 0.0044 con 3e-4) y
   coincide con el perfil real de Δχ² = 1, que da 0.0040.
   Corolario general: **cuando la curvatura local y el cociente de verosimilitud
   no coinciden, gana el segundo**, y la diferencia es información y no ruido
   numérico. En el test del origen la brecha (2.1 contra 4.7 en el Co-60) dice
   que el χ² es mucho más plano lejos del mínimo de lo que su curvatura sugiere,
   o sea que c0 está degenerado con el exceso.

50. **El test del origen con UNA sola fuente da una cota, y de paso mide su
   propia debilidad** (28-ago-2026). Liberando la ordenada (`canal = D·E + c0`)
   sobre una corrida: Co-60 da c0 = −163 ± 35 mV·ns con Δχ² = 4.3 (2.1σ) y
   Cs-137 da +60 ± 34 con Δχ² = 0.5 (0.7σ). **Signos opuestos**, o sea no hay un
   pedestal consistente; y liberar c0 mueve D un 13% y un 7% respectivamente por
   una mejora de χ² casi nula, que es la degeneración c0-σ-exceso cuantificada.
   Sirve como cota ("si hubiera un pedestal grande, se vería") y no reemplaza al
   test con dos fuentes. Está en `calibracion.test_ordenada`.

51. **LOS DATOS NUEVOS NO SON XML: son el binario del DRS4** (01-sep-2026).
   La carpeta `D:\Sesión Calibración 26_08_2026\` tiene DOS sesiones adentro (las
   marcas de tiempo de cada evento lo dicen): el Cs-137 del detector A del 26-ago
   en XML, y del **01-sep los cuatro detectores A, B, C, D a 1600/1650/1700 V en
   formato BINARIO y sin extensión en el nombre** (`Co60`, `Cs137`, `B1700V`).
   Cs-137 en los cuatro, **Co-60 sólo en el A**, ninguna con fondo propio.
   → **`docs/formato_binario_drs4.md` tiene el formato byte a byte, un lector
   mínimo verificado y las trampas.** En el código no hay que hacer nada:
   `analisis/drs.py` elige el lector mirando los primeros 4 BYTES (`DRS2`/`TIME`),
   no la extensión, y `analisis/drs_bin.py` expone la misma interfaz que
   `drs_xml`. El catálogo pasó de 4 a 16 corridas.
   Dos cosas que se dedujeron del archivo en vez de cablearlas, por §7.36: si el
   scaler está o no (se prueban los dos tamaños de evento y se verifica que
   `EHDR` caiga donde predice), y cuántas placas y canales hay.

52. **En la sesión del 01-sep la compuerta de integración NO CABE en el
   registro** (01-sep-2026). El trigger quedó en t ≈ 123 ns de un registro de
   200 ns: después del onset sobran ~81 ns y el pulso necesita **83** para
   converger al 99.9%. Con la ventana que pide el pulso, el **84% de los eventos
   sale truncado**, y un evento truncado pierde una fracción de carga que depende
   de CUÁNDO llegó el pulso — justo el error que la compuerta anclada al onset
   existe para evitar (§7.33).
   → `area.calibrar_ventana` ahora mide el espacio disponible tras el onset **en
   los eventos de amplitud sobre la mediana** (el walk del trigger pone el onset
   tarde en los pulsos chicos y sesga el percentil crudo), descuenta 5 ns de
   margen y recorta `post`; devuelve `post_pedido`, `recortado`, `fraccion` y
   `frac_eventos_enteros`. Medido en el Co-60 a 1600 V: post 83 → 99.9% de carga
   y 16% de eventos útiles; post 80 → 99.5% y 71%; **post 75 → 98.0% y 94%**, con
   D corregida por carga idéntica (1.948 contra 1.949). Los eventos que se
   pierden son todos de amplitud < 20 mV. **Para la próxima toma: adelantar el
   retardo de trigger 40-50 ns.** No cambia nada en los XML viejos (ahí el pulso
   cae en t ≈ 48 ns y sobran 120-155 ns): `pipeline.py co60` sigue dando 1.2512.

53. **Para COMPARAR dos corridas hay que usar la MISMA compuerta** (01-sep-2026).
   El `post` que mide cada archivo depende de su propio pulso, y dos corridas con
   `post` distinto recogen fracciones distintas de la carga, lo que entra derecho
   en D. Las tres corridas de Cs-137 del 26-ago habían medido 55, 98 y 124 ns:

   | | b de D ∝ V^b | χ² | residuos |
   |---|---|---|---|
   | compuerta de cada archivo | 10.4 | 744 (1 gdl) | ±14% |
   | corregida por fracción de carga | 10.0 | — | ±2.3% |
   | **compuerta común de 75 ns** | **10.1** | **9.2 (1 gdl)** | **±1.5%** |

   → `pipeline.calibrar_completo(c, {"post_fijo": 75.0})`.

54. **TEST B CERRADO: la ganancia sigue una ley de potencia con b = 9.49**
   (01-sep-2026). Co-60 del detector A a 1550 (junio, XML) + 1600/1650/1700
   (1-sep, binario): **b = 9.49 ± 0.04, χ² = 7.6 con 2 gdl** contra 969 de la
   recta, residuos ±0.9%. **b cae dentro del 6-10 típico de un PMT**, que es la
   primera vez (§7.46 daba 13.9 combinando puntos que no eran comparables). De
   yapa: la corrida de junio cae sobre la misma curva a −0.08%, o sea entre junio
   y septiembre **no cambió el atenuador** ni nada de la cadena, que era la duda
   abierta de §7.45. Los otros detectores: B 9.06, C 9.07, D 7.68 (crudo).
   ✔ **Control cruzado (22-sep-2026):** el operador confirmó que fueron 20 dB en
   todas las sesiones (§7.70). O sea que la medición —que no sabía nada del
   atenuador— y la bitácora dicen lo mismo. Es la clase de coincidencia que vale
   la pena anotar: la ley de potencia se ajustó sin usar ese dato, y si la
   atenuación hubiera cambiado entre junio y septiembre el punto de junio habría
   caído un factor 10 afuera, no un 0.08%.

55. **TEST A: las dos fuentes NO comparten ganancia, y no es un pedestal**
   (01-sep-2026). Con compuerta común de 75 ns, detector A:

   | HV | Co-60 | Cs-137 | Cs/Co |
   |---|---|---|---|
   | 1600 V | 1.6983 ± 0.0041 | 1.8612 ± 0.0033 | 1.096 |
   | 1650 V | 2.2570 ± 0.0045 | 2.5237 ± 0.0057 | 1.118 |
   | 1700 V | 3.0152 ± 0.0064 | 3.4450 ± 0.0073 | 1.143 |

   Liberando la ordenada da +154 ± 5, +252 ± 7 y +407 ± 9 mV·ns. **Un pedestal
   del digitalizador es carga fija por evento: no depende del HV**, así que
   tendría que dar lo mismo a los tres voltajes. χ² = 613 con 2 gdl → no es
   pedestal. Traducido a voltaje con b = 9.49, el cociente equivale a **+0.97%,
   +1.18% y +1.41% de HV** (16-24 V): las dos sesiones no estaban al mismo punto
   de operación real. Descartado que sea del análisis: con plantilla analítica
   Klein-Nishina el desacuerdo sigue (+11.0%, +5.5%, +12.0%) y con el 50% de
   altura también (+6.9%, +4.8%, +8.5%).
   Tampoco es la atenuación: las dos sesiones tienen los mismos 20 dB (§7.70),
   y de todos modos un cambio de atenuador daría un factor fijo, no un cociente
   que crece con el HV (1.096 → 1.143).
   → **El test A no se cierra hasta medir Cs-137 y Co-60 en la MISMA sesión.**
   `registro.informe` ya no afirma "hay un pedestal" mirando un voltaje solo:
   junta las ordenadas de todos los puntos y separa las dos hipótesis
   (`_pedestal_o_ganancia`).

56. **`replica_paper/` tiene interfaz, y construirla cerró la única ecuación
   del paper que quedaba sin verificar** (02-sep-2026). `replica_paper/gui.py`:
   14 pestañas en el orden del paper (Mapa · Datos · Ec. 1 · Ec. 2 y 10 ·
   Ec. 8-9 · Ec. 15 · II-A · Fig. 2 · convencionales · Fig. 5 · II-C ·
   auditoría · Sec. III · cierre), con controles a la izquierda, gráfico al
   centro y, a la derecha, cuatro secciones fijas: **la cita literal del
   paper**, **qué función la implementa**, **lo que el paper NO dice** (que son
   exactamente los controles) y **qué mirar**. Dos reglas que la hacen
   confiable: no calcula nada por su cuenta —todo sale de llamar a los módulos
   de `replica_paper/`, así que no puede divergir de la réplica— y **no cachea
   resultados**: cada «Recalcular» rehace el cálculo desde el espectro. Lo único
   que se reusa es el `.areas_gate.npz` de `analisis/area.py`, que es el PARSEO
   del archivo (clave = tamaño + fecha), no un resultado.
   Lo que apareció al construirla:
   - **La Ec. 2 y la Ec. 10 cierran numéricamente, y nadie lo había medido.**
     `ecuaciones.dsigma_domega` estaba implementada explícitamente «para poder
     verificar que la Ec. 10 es su transformada» y esa verificación no existía
     en ningún lado. Integrando cada una por su propia variable (ángulo sólido
     y energía del electrón) dan la MISMA sección eficaz total a 6 decimales
     para E₀ = 500, 662, 1000 y 2000 keV (2.891788e-25, 2.562107e-25,
     2.112170e-25 y 1.463654e-25 cm²). O sea el jacobiano del paper está bien,
     y ahora la pestaña 2 lo recalcula en vivo en vez de afirmarlo.
   - **El costo real de «sin caché», medido:** las pestañas rápidas (Mapa,
     Ec. 1, Ec. 2/10, Ec. 8-9, auditoría) están abajo de 0.5 s y calculan
     solas al abrirlas; la Fig. 2 completa tarda 46 s sobre el Co-60 con el
     diagnóstico del solapamiento, la Fig. 5 tarda 127 s sobre el Cs-137
     (siete cadenas NCC-SA), la Sección III 46 s con las dos corridas y el
     cierre 52 s. Los 10 bloques son lo más caro y ya estaba medido por
     `correr.py`: 53 s en el Cs-137 y 453 s en el Co-60. Todas esas no
     calculan solas: hay que apretar «Recalcular», igual que en
     `codigo/gui.py`.
   - Dos cambios chicos en `correr.cierre()`, sin tocar lo que imprime: acepta
     `D_real`, `sigma_rel` y `semilla` (antes cableados) y devuelve
     `(analisis, rectas)`, para que la interfaz dibuje lo mismo que el informe
     dice en vez de regenerar los espectros por su cuenta.
   - **Un default que no se puede dejar librado al azar: `ligar_origen`.**
     `correr.tres_metodos` lo tiene en `True` (D = B₁/E_max,1, lo que impone la
     Ec. 8) y `metodo.ajustar_dos_bordes` en `False` (D libre, la Ec. 11 tal
     como se imprime). Una pestaña que llamara directo a `metodo.ajustar`
     mostraría B₁ = 1148.4 mientras el informe dice 1124.2, sin que nada
     falle. La interfaz usa `True` en todas menos en la pestaña 6, que es
     donde el punto ES poder destildarlo.
   - Trampa de Python que costó una corrida entera: `"50% de altura %8.1f"` es
     un format string de TRES conversiones, no de una, porque `% d` es un
     `%d` con flag de espacio. Hay un escáner de esto en la nota de abajo.
   - Y una de gráfico que era un error real de comparación: los seis
     intervalos de la Fig. 5 se numeran 1..6, pero el más angosto se saltea
     cuando queda por debajo de `min_bins`. Numerar 1..N los que sobreviven
     desalinea la curva contra los seis paneles del paper (nuestro «1» sería
     su «2»). Ahora se conserva el índice original y el salteado se marca.

   Cómo se verificó, dos pasos:
   - `smoke.py` instancia la `App` y, para cada pestaña, llama
     `leer_controles()`, `trabajo()` y `pintar()` en el hilo principal, sin
     `mainloop`. Es la forma barata de probar una GUI entera sin tocarla, y
     encontró el bug del format string.
   - `verificar.py` compara lo que devuelve la interfaz contra
     `replica_paper/informe.txt`, valor por valor. **Los 20 coinciden**
     (Cs-137: intervalo [579, 1389], ρ = 0.99618, B = 915.2, σ = 115.5,
     χ²/ndf = 2.02, Gauss 951.7, Diferencia 829.0; Co-60: [952, 1967],
     ρ = 0.99095, B₁ = 1124.2, B₂ = 1304.7, D = 1.1669, Gauss 1341.3,
     Diferencia 1110.9), y la robustez de la Fig. 5 sobre el Co-60 reproduce
     sus 46.2% / 1.3% / 88.2%. Es la prueba de que la interfaz no calcula
     nada por su cuenta.

57. **`codigo/` ya lee el binario, y el problema no era leerlo sino ENCONTRARLO**
   (02-sep-2026). `analisis/` sabía leer los binarios desde §7.51, pero
   `codigo/datos_io/spectrum.py` no, y las cuatro interfaces del proyecto
   abrían el diálogo de archivos filtrando `*.xml`. Con ese filtro las carpetas
   de la sesión del 1-sep **se ven vacías**, porque esos archivos no tienen
   extensión: la falla se presentaba como «no hay datos», no como un error.
   Tres arreglos, en orden de lo que costó encontrarlos:
   - **El filtro del diálogo.** Ahora «Datos del DRS4 (XML o binario) *.*» va
     PRIMERO en `codigo/gui.py`, `analisis/gui.py` y `pulsos/gui.py`. El
     formato lo decide el contenido (`spectrum.formato` / `drs.formato`, los
     primeros 4 bytes), nunca el nombre, y la etiqueta de cada archivo lo dice
     junto con el canal detectado — que es el diagnóstico que faltaba en
     `info_archivo`, gateado hasta ahora por `ruta.lower().endswith(".xml")`.
   - **`spectrum.py` delega, no duplica.** `detectar_canal` y `_iter_eventos`
     desvían a `analisis/drs_bin.py` cuando el archivo es binario; el parser
     byte a byte sigue existiendo en un solo lugar. Para que los consumidores
     no ramifiquen, `_iter_eventos` entrega en el lugar de `pares` un array
     (n, 2) en el binario y la lista de cadenas "t,v" en el XML, y
     **`spectrum.pares_a_pulso`** unifica los dos (la usan `gui.TabArea` y
     `documentacion.py`, que hacían el `np.fromstring` a mano).
   - **El fondo pasó a ser opcional.** Ninguna corrida del 1-sep tiene fondo
     propio. `espectro_neto(muestra, None)` devuelve neto = muestra con
     `k_f = 0` y `hay_fondo = False`, y la GUI lo dice en la etiqueta. Prestarle
     el fondo de otra sesión sería PEOR que no restar nada: el fondo depende del
     punto de operación (a 1700 V el mismo cósmico deja el triple de área), así
     que restar el de otro voltaje mete una forma ajena en el espectro.
   Además: la letra del disco externo **no es una propiedad del disco**, se la
   asigna Windows. `rutas._buscar_sesion` barre las unidades buscando la carpeta
   de la sesión (y respeta `CALIBRACION_SESION` si se quiere apuntar a una
   copia local), en vez de dar todo el catálogo por [NO MONTADO] cuando el
   disco entra como E:. Y `analisis/gui.py` expone por fin la **compuerta común**
   (`post_fijo`, §7.53), que es lo que hace comparables dos corridas y hasta
   ahora solo se podía pasar desde `pipeline.py`.
   Verificado con **`verificaciones/_verif_binario.py`**, que FABRICA un binario
   del DRS4 con verdad conocida (el disco externo no siempre está conectado) y
   comprueba que el eje temporal vuelve exacto, la tensión dentro de medio
   cuanto del ADC (0.0076 de 0.0153 mV), y que `codigo/spectrum` y
   `analisis/area` sacan **exactamente las mismas áreas** (diferencia relativa
   máxima 0.0e+00), que además cierran con la integral analítica del pulso a
   1.0000. Sin regresión sobre los XML: `calibrar_v2.py --rapido` sigue dando
   D_luz = 1.2618 ± 0.0133 y `pipeline.py co60` D = 1.2512 ± 0.0047.

   **Addendum (02-sep-2026): tres cosas que salieron al perseguir un cierre
   abrupto de `analisis/gui.py` que NO se pudo reproducir.**
   - **Un bug propio, ya arreglado:** `rutas._buscar_sesion` sondeaba de la A a
     la Z a ciegas. En Windows tocar una unidad extraíble presente pero vacía
     (lectora de tarjetas, DVD) tarda segundos y puede abrir el diálogo modal
     "No hay disco en la unidad", que sin consola parece un cuelgue. Ahora
     `rutas._unidades()` pide el mapa de bits a `GetLogicalDrives` y sólo
     sondea las que existen (import de `rutas`: 0.03 s).
   - **Una excepción en un callback de Tk era INVISIBLE.** Tk las imprime en
     stderr y sigue andando: sin consola no van a ninguna parte y la ventana
     queda a medio dibujar, que desde afuera se ve igual que un cierre. Ahora
     `App.report_callback_exception` las escribe en `.cache/gui_error.log` y
     las muestra en un cuadro.
   - **Un cierre abrupto no deja excepción**, así que hay que instrumentar:
     `faulthandler` vuelca la pila de todos los hilos en `.cache/gui_fallo.log`
     (el sospechoso típico con estos datos es un `EXCEPTION_IN_PAGE_ERROR` de
     Windows si el disco USB se duerme mientras `drs_bin` lee por `np.memmap`),
     y `App._traza` deja una línea por etapa con flush inmediato en
     `.cache/gui_traza.log`, cuya última línea dice hasta dónde llegó.
   Lo que SÍ se verificó sin fallar: las 16 corridas del catálogo por la cadena
   completa con todas las figuras, la corrida real con hilo + `mainloop` y las
   7 pestañas, seis variantes de "archivos sueltos" (binario, binario + sim,
   sim que no corresponde, XML de 2.8 GB con fondo y sim, un .txt y un .npz —
   las tres últimas dan error controlado), el parseo EN FRÍO del XML de 2.8 GB
   desde el USB con Tk vivo (86 s), y lanzar `python gui.py` parado dentro de
   `analisis/`.
   ⚠️ **Detector D da un ajuste degenerado** y conviene mirarlo aparte:
   χ²/ndf 3.3-3.7 contra ~1 de los demás, error de D nan o 0 (la covarianza sale
   no definida positiva, con `invalid value encountered in sqrt`), y D = 0.55 /
   0.69 / 0.86 contra 1.46-2.86 de B y C. No es un crash, pero tampoco es un
   número usable.

58. **`analisis/gui.py` muestra ahora los CUATRO estados de la cadena, y armar
   la tercera pestaña destapó que el paisaje del χ² estaba roto**
   (03-sep-2026). Las pestañas pasaron de 7 a 9 y las cuatro del medio son la
   cadena en el orden en que corre, con cuatro paneles cada una:
   **1 · Datos y simulación** (muestra/fondo, neto con la ventana, la plantilla
   de Geant4 en energía, y las dos superpuestas con la D ajustada),
   **2 · NCC e intervalo** (ρ(B), la ventana ganadora, las dos formas
   estandarizadas que Pearson compara, y el diagrama de dispersión de una
   contra otra), **3 · Semilla y recocido** (la semilla del NCC con las cotas,
   el costo ciclo a ciclo de los N arranques, T y aceptación, y el camino en el
   plano D-σ) y **4 · Inicial vs final** (los dos modelos sobre el dato, los
   residuos de los dos, cuánto se movió cada parámetro en % y en unidades de su
   propio error, y el corte del χ² sobre el segmento semilla→óptimo).
   Cuatro cosas que salieron de construirlas:
   - **`Problema.mapa` reventaba con el modelo por defecto, en silencio.**
     Barría (D, σ) dejando `lam` en `None`, así que `costo` hacía
     `not None > 0` → `TypeError`. Como el dibujo pasa por `Lienzo.dibujar`,
     que atrapa todo y escribe el error DENTRO del gráfico, la pestaña
     mostraba un cartel gris y nadie se enteraba. Ahora `mapa` pide `x_base` y
     congela el resto de los parámetros en el óptimo; lo que se dibuja es un
     CORTE del χ², no el perfilado (eso pediría reoptimizar en cada punto de
     la grilla).
   - **El detalle del recocido no salía de `calibrar` y había que sacarlo sin
     tocar el resultado.** `calibracion.calibrar` acepta un dict `detalles` y
     lo llena con la semilla completa, las cotas y las N corridas; `pipeline`
     lo deja en `salida["sa"]`. No entra en `r` porque `r` se serializa a JSON
     y esto son arrays. Y `sa.recocido` guarda ahora el CAMINO
     (`historia_x`, `historia_mejor_x`): una fila de floats por ciclo de
     temperatura, que es lo único que permite ver si el recocido cruzó de un
     valle a otro o bajó derecho.
   - **La pestaña 2 recalcula ρ y lo muestra al lado del ρ del barrido.** Los
     dos dan 0.97746 en el Co-60, o sea que la pestaña reconstruye bien lo que
     el barrido hizo. Es la misma idea que `replica_paper/verificar_gui.py`: la
     interfaz no tiene permiso de inventar números.
   - **Las escalas hay que elegirlas o no se ve nada.** Tres casos medidos: en
     el panel de "inicial vs final" el pico de bajo canal (880 cuentas contra
     150 del plateau) se lleva la escala y las dos curvas quedan pegadas; en el
     paisaje, el χ² de las paredes es ~100 veces el del valle y con el rango
     completo el valle queda de un solo color (se corta en 5 veces el mínimo);
     y en el costo contra ciclo hay que poner los límites a mano porque
     `costo` devuelve 1e30 fuera de lo físico y un solo valor así aplasta la
     escala logarítmica.
   Verificado con **`verificaciones/_verif_gui_pestanas.py`**, que corre la
   cadena y llama a cada método de figura DIRECTAMENTE (sin `Lienzo`, así que
   una excepción se ve) y deja los PNG en `figuras/gui_*.png`: las 8 dibujan
   sobre el Co-60 y sobre el Cs-137 —o sea también con fuente de UNA sola
   línea—. Sin regresión: `pipeline.py co60` y `calibracion.py` siguen dando
   D = 1.2512 ± 0.0047 (keVee) y 1.1973 ± 0.0038 (keV), y las demos de
   `sa.py` y `modelo.py` no se movieron.

59. **El alto de los gráficos de `analisis/gui.py` lo decidía la lista de
   etapas, y el espaciado se calculaba una sola vez** (03-sep-2026). Las cuatro
   pestañas de §7.58 son grillas de 2×2: lo que se ve depende del alto de cada
   panel. Medido con la ventana maximizada (1920×991 en este equipo):

   | | lienzo | fila de la grilla |
   |---|---|---|
   | lista de etapas arriba, `tight_layout` | 11.95 × 5.53 in | 2.76 in |
   | lista en la columna izquierda, `constrained` | 11.19 × 7.20 in | 3.60 in |

   Los cuatro cambios, en orden de cuánto aportan:
   - **La lista de etapas se movió a la columna izquierda.** Se lleva ~170 px de
     ALTO y sólo necesita ~430 de ancho; la columna izquierda tiene ese ancho y
     le sobra alto. Hay un interruptor «ver etapas» para plegarla, y NO se
     pliega sola al terminar: los tiempos por etapa son diagnóstico, y una
     ventana que se reacomoda sin que nadie la toque se lee como un error.
   - **`gui.ajustar(fig)` (motor `constrained`) en lugar de
     `fig.tight_layout()`.** `tight_layout` resuelve el espaciado UNA vez, con
     el tamaño que la figura tiene en ese momento, y el lienzo de Tk se estira
     con la ventana sin volver a llamarlo: los títulos se terminan pisando en
     cuanto la ventana no mide lo que medía al calcular. ⚠ Va como PRIMERA
     línea de cada figura, no como última: `fig.colorbar` mira si hay motor de
     layout y, si no lo hay, mete la barra en un gridspec anidado con
     `height_ratios = [0, 1, 0]`; encender el motor después hace que el layout
     divida por ese cero, y el `ZeroDivisionError` salta al DIBUJAR, lejos de
     la causa.
   - **La ventana arranca maximizada** y la tipografía se bajó a 9 pt (8 en los
     ejes), a nivel de módulo para que las figuras del verificador salgan
     iguales.
   - **Las tablas pasaron a 4 columnas** (unidades pegadas al nombre) y la
     columna del panel con tabla va 1.35 veces más ancha que la otra.

   Cuatro trampas de matplotlib que costaron encontrar, todas de la misma
   familia (fallan callado, deformando el dibujo en vez de dar error):
   - **`ax.table(bbox=...)` es (x0, y0, ANCHO, ALTO)**, no (x0, y0, x1, y1).
     Pasarle `(0, 0.09, 1, 1)` para «dejar 9% abajo» da una tabla de alto 1 que
     arranca en 0.09: se sale por arriba de los ejes y tapa el título.
   - **Un `ax.text` largo cuenta como decoración del eje**: constrained layout
     mide el tight bbox con el texto adentro y ACHICA el eje hasta que entre.
     Con una nota de 100 caracteres en una línea, la tabla quedó de un quinto
     del panel. Por eso `tabla` parte la nota en líneas de 46 caracteres y le
     reserva 0.048 de alto por línea (a 7.5 pt cada línea mide ~0.14 in y el
     eje de una grilla 2×2 en ventana chica mide ~2.9).
   - **Un título largo con `loc="left"` se sale del ancho de sus ejes** y se le
     monta al título del panel vecino; el motor no lo acorta ni lo parte en
     líneas. Los títulos de dos líneas están cortados a mano.
   - **matplotlib no achica el texto de una celda: lo RECORTA.** Y elegir el
     cuerpo según el ancho de la figura NO sirve: la figura se dibuja una vez,
     al terminar el cálculo, y después el lienzo la reescala sin volver a
     ejecutar ese código, así que lo que se mediría es el tamaño por defecto.
     Lo único que funciona en los dos tamaños es contenido corto.
   Verificado con `_verif_gui_pestanas.py` (Co-60 y Cs-137, las 8 pestañas) y
   con la ventana real a 1920×991 y a 1400×860, guardando lo que el lienzo
   dibuja de verdad. Sin regresión: `pipeline.py co60` sigue en D = 1.2512 ±
   0.0047 y `cs137_1600` en 1.9450 ± 0.0038.

60. **La interfaz se cerró de golpe «durante el SA», y el SA no tenía nada que
   ver** (03-sep-2026). El reporte fue «crashea al ejecutar el SA y sin mensaje
   del problema». Lo que dicen los dos logs que la GUI ya escribía:
   - `gui_traza.log` llega hasta `fin sa (2.1 s)`, `fin calib` y
     `--- terminado, D = 2.3984`. O sea que **el recocido terminó**, la
     calibración también, y el número salió. Se veía como «el SA» solamente
     porque la lista de etapas era lo último que se movía en pantalla.
   - `gui_fallo.log` tiene `Windows fatal exception: access violation`, y en el
     hilo principal la pila es **solo `mainloop`**: ningún frame de Python. Es
     decir que el proceso murió dentro del C de Tk procesando eventos, no
     adentro de un callback nuestro —si hubiera sido nuestro, `faulthandler`
     mostraría `_vaciar_cola → _mensaje → _terminado → ...`.
   El sospechoso que queda es lo que pasa justo ahí: `_terminado` pintaba las
   **siete** figuras de una, y las seis que no están al frente viven en un
   lienzo que Tk todavía no mapeó, o sea de **1×1 px** (medido). TkAgg vuelca
   cada render en una PhotoImage que después dibuja **Tk por su cuenta**, sin
   Python en el medio — exactamente el lugar donde apareció la violación de
   acceso.
   ⚠ **No está probado**: 60 vueltas pintando en lienzos de 1×1 no alcanzaron
   para reproducir la caída, ni tampoco tres calibraciones seguidas sobre las
   corridas B, C y D a 1700 V. Así que lo que se hizo son dos cosas distintas,
   y conviene no confundirlas:
   - **Sacar el trabajo frágil.** `Lienzo` dibuja ahora cuando la pestaña SE
     VE (`programar` / `asegurar` / `<Configure>`, con `MINIMO = 40 px`).
     Terminar una corrida cuesta UN dibujo en vez de siete, y no se blitea en
     ningún lienzo degenerado. `limpiar()` descarta el pendiente al empezar
     otra corrida, si no un dibujo viejo aparecería sobre un resultado nuevo.
   - **Hacer que un cierre así HABLE**, que es lo que se pidió. Eran tres
     agujeros: (a) `faulthandler` abría `gui_fallo.log` en modo `"w"`, así que
     **relanzar el programa borraba la única evidencia** —ahora se rota a
     `gui_fallo_anterior.log` y la ventana avisa al abrirse, con la pila en el
     panel de Resultado—; (b) la traza no decía nada del dibujado, así que un
     cierre posterior al cálculo parecía del cálculo —ahora hay una línea por
     cambio de pestaña y dos por figura, antes y después—; (c) el cartel gris
     de `Lienzo.dibujar` se comía la excepción sin dejar rastro —ahora también
     va a `gui_error.log`.
   Descartado por medición, no por intuición: el ajuste de esa corrida es sano
   (χ²/ndf = 1.05, sin NaN ni infinitos, límites de ejes finitos en las 8
   figuras) y los hilos de `fsspec`/`asyncio` que abre uproot **no se
   acumulan** entre corridas (quedan 4 hilos, corrida tras corrida).
   Moraleja para la próxima: cuando la ventana desaparece, la primera pregunta
   es qué dice la ÚLTIMA línea de `.cache/gui_traza.log`, y la segunda si
   `.cache/gui_fallo_anterior.log` tiene una pila. Las dos existen justamente
   para eso.

61. **Las 12 GB de la sesion se copiaron a `datos/`, ordenadas por detector y
   voltaje** (03-sep-2026). Hasta ahora las 16 corridas del catalogo vivian
   solo en un disco externo cuya letra ni siquiera es una propiedad del disco
   (§7.57), asi que la mitad del proyecto no se podia correr sin enchufarlo.
   La copia queda en `datos/Detector X/<voltaje>/<muestra>` y respeta lo que
   ya funcionaba, con tres arreglos:
   - **La carpeta del detector A a 1700 V se llamaba `1700`** y todas las
     demas `1700V`. En la copia se normalizo; `rutas._local` aplica la misma
     regla al buscar.
   - **`B1700V` no dice de que fuente es.** Es el Cs-137 del detector B a
     1700 V y en la copia se llama `Cs137`. Era el unico nombre opaco;
     `rutas._RENOMBRE_LOCAL` lo traduce.
   - **Se copia con `shutil.copy2`, que conserva la fecha de modificacion**,
     porque de eso depende el cache de areas... y ahi aparecio un bug: la
     clave de `analisis/area._clave` incluia la **ruta absoluta**, asi que
     copiar (o mover) un archivo invalidaba su cache aunque fuera byte a byte
     el mismo. Medido: `cs137_1600` desde la copia re-parseo 80 s y dio
     exactamente la misma D (1.9450 +- 0.0038), o sea trabajo puro al vacio.
     Y contradecia lo que §7.43 afirmaba del proyecto -que mover carpetas no
     costaba nada-, que era cierto para `codigo/spectrum` (clave = tamano +
     fecha) y falso para `analisis/area`. Se saco la ruta de la clave: ahora
     identifica al CONTENIDO. Costo una vez re-hacer los .npz de las 16
     corridas; de aca en adelante mover los datos es gratis.
   `rutas._ext` busca **primero en `datos/` y despues en el disco**, o sea que
   el catalogo entero funciona con el disco desconectado y sigue funcionando
   como antes si algun archivo no se copio. Cada carpeta lleva un `LEEME.txt`
   generado leyendo los archivos (no el nombre): que fuente es, que formato
   tiene y por que canal salio. Hace falta porque los binarios del DRS4 no
   tienen extension y en los dialogos de archivo hay que dejar el filtro en
   `*.*` para verlos (§7.57).
   ⚠ `datos/` esta dentro de OneDrive, asi que la copia se sincroniza a la
   nube: son ~14 GB en total. Si eso molesta, la escotilla es la variable de
   entorno `CALIBRACION_SESION`, que apunta la sesion a donde uno quiera.

62. **Los graficos de `analisis/gui.py` salian LITERALMENTE cortados, y no era
   el espaciado: faltaba un pedazo del dibujo** (08-sep-2026). TkAgg escala la
   figura por el *device pixel ratio* de la pantalla —en Windows,
   `tk scaling / (96/72)`, que con el `tk scaling = dpi/72` que declara la App
   (§7.24) vale **1.25**— pero lo aplica al recibir `<Map>`, o sea DESPUES de
   que `FigureCanvasTk.resize` ya fijo el tamano en pulgadas con el dpi viejo.
   Lo unico que hace para reconciliarlos es pedirle al widget que crezca, y
   adentro de un Notebook dentro de un Panedwindow ese pedido no se cumple: no
   llega ningun `<Configure>` nuevo, el buffer de Agg queda mas grande que la
   PhotoImage sobre la que se vuelca, y el blit **recorta en silencio**. Con la
   ventana en 1550x830: widget 1029x739 px contra un buffer de 1286x923, o sea
   que **el 20% de la derecha y el 20% de abajo no se dibujaban nunca**.
   ⚠ Un recorte asi NO se ve en un PNG guardado con `savefig`, que es lo que
   probaba `_verif_gui_pestanas.py` (§7.58): solo se ve en la ventana. De ahi
   **`verificaciones/_verif_gui_encuadre.py`**, que compara los DOS tamanos
   —el del render y el de la imagen— en cada pestana y a tres tamanos de
   ventana. Cuatro consecuencias, todas de la misma causa:
   - **El lienzo PIDE 360x240 px** (`Lienzo.PEDIDO`) y se conforma con el que
     le den; `Lienzo._encuadrar` le da a la figura el tamano exacto del widget
     con el dpi que tenga en ese momento. Sin lo primero el Panedwindow pedia
     2490x1194 px —mas que la pantalla— y en `pack` los hijos que se empaquetan
     al final se quedan sin lugar: con la ventana sin maximizar **desaparecian
     la barra de progreso de abajo y la barra de herramientas de matplotlib**.
     La barra de abajo se empaqueta ahora PRIMERO.
   - **El dpi de la figura va en 100**, el nominal. Estaba en `100·dpi/96` = 125
     para compensar el escalado de Windows, que es exactamente lo que TkAgg ya
     hace: el factor entraba dos veces, la figura terminaba a 156 dpi y todo se
     dibujaba un 25% mas grande de lo pedido, o sea con un **36% menos de
     superficie util**. (§7.24 sigue valiendo para `tk scaling`, que es de Tk;
     lo que sobraba era aplicarselo tambien a la figura.)
   - **La columna de controles se desplaza** cuando no entra: el Panedwindow
     reparte el ancho pero el alto lo impone la ventana, y los controles piden
     907 px contra los 830 de la ventana sin maximizar. La barra aparece solo
     si hace falta, y la decision se toma **en reposo** (`after_idle`) y no
     dentro del `<Configure>`, porque ahi `winfo_height` todavia puede ser el
     alto de antes.
   - **Doble clic sobre un panel lo agranda** a la figura entera (o el selector
     que esta al lado de la barra de herramientas); doble clic de nuevo vuelve.
     Aun sin recortes, cuatro paneles en 1008x697 px son ~500x350 cada uno. No
     se recalcula NADA: se guardan posiciones, visibilidad y anclas, se apaga el
     motor de layout y se reubica el panel elegido —lo unico que lo hace usable
     en el paisaje del χ², donde redibujar costaria ~6300 evaluaciones del
     modelo—. Los `twinx` se detectan por posicion (nacen con el mismo
     rectangulo que su panel) y las barras de color se anotan mientras se
     dibuja, envolviendo `fig.colorbar`: despues no hay como saber de quien es
     cada una.
   De paso cae media trampa de §7.59: «elegir el cuerpo de la tabla segun el
   ancho no sirve porque la figura se dibuja una vez y despues el lienzo la
   reescala» valia mientras la figura se construia con el tamano de fabrica.
   Ahora `_encuadrar` le pone el tamano real ANTES de llamar a la funcion de
   dibujo, asi que `tabla` saca el cuerpo del ancho del panel (5.8 a 8 pt) y
   dejo de recortar los encabezados con la ventana chica.
   `visor.py` hereda todo, porque usa el mismo `Lienzo`.

63. **El recorte de §7.62 se reporto de nuevo, NO se pudo reproducir, y la
   respuesta fue cerrar la mecanica en vez de perseguir el orden de eventos**
   (14-sep-2026). El reporte: con solo la pantalla del portatil los graficos
   siguen saliendo incompletos; al conectar un segundo monitor la ventana se
   "ajusta" y, devuelta al portatil, se ve bien. Lo que se midio en este
   mismo equipo (con el monitor externo conectado, que es lo unico que se
   pudo probar): arranque real maximizado por `mainloop`, las 7 pestanas,
   normal/maximizada, ida y vuelta al monitor externo (96 dpi contra 120), y
   con los dos Python instalados (3.13 y 3.14, ambos matplotlib 3.10 y Tk
   8.6.15). En TODOS los estados widget = render = imagen, sin recorte.
   Dato util de paso: Tk 8.6.15 **no cambia `tk scaling`** al cruzar de
   monitor (sigue en 1.668 con la ventana en el de 96 dpi), o sea que el
   *device pixel ratio* de TkAgg queda en 1.25 en los dos y el sistema es
   autoconsistente en pixeles.
   El recorte tiene UNA sola mecanica: Agg renderiza lo que mide la figura
   (pulgadas x dpi) y lo vuelca sobre una PhotoImage de tamano FIJO que mide
   lo que dijo el ultimo `<Configure>`; `Tk_PhotoPutBlock` no agranda una
   imagen con tamano fijado, la recorta por la derecha y por abajo. Lo que
   §7.62 hacia era mantener a los dos de acuerdo reaccionando a `<Configure>`
   y `<Map>`, y eso depende del orden en que Tk los entrega. Tres cosas:
   - **`gui._CanvasEncuadrado`**, el canvas de todo `Lienzo`: justo antes de
     renderizar compara widget, render e imagen; si no coinciden, corrige los
     dos que dependen de la figura y escribe en la traza
     `!!! reencuadre forzado en <pestana>: widget AxB, render CxD, imagen
     ExF, dpi, ratio`. Probado provocando a mano los dos estados rotos (la
     figura a 156 dpi sin `<Configure>`, y la imagen achicada): los dos
     quedan corregidos en el siguiente render. En el arranque real el
     contador `forzados` queda en 0, o sea la red no esta tapando nada hoy.
   - **La traza ahora dice con que tamano se dibujo cada cosa.** Al abrir:
     `entorno: pantalla, N monitor(es), tk scaling (px/in)`; cada vez que la
     ventana se asienta en un tamano, posicion o monitor nuevo: `ventana:
     WxH en (x, y), zoomed|normal, dpi del monitor`; y por lienzo, una linea
     `render <pestana>: widget, render, imagen, dpi, ratio` por cada tamano
     NUEVO (en reposo, 500 ms, para no escribir una por `<Configure>` al
     arrastrar). `App._traza` paso a ser `gui.trazar`, a nivel de modulo,
     para que los lienzos puedan escribirla.
   - **`_verif_gui_encuadre.py` prueba primero el arranque REAL**: maximizada
     de entrada y `mainloop` con pasos encadenados por `after`, y no
     `update()` en bucle con la ventana en `normal`, que era lo unico que
     probaba y no es lo que el usuario ve. Si hay mas de un monitor, pasea
     la ventana por el otro y la trae de vuelta. En cada estado exige
     widget = render = imagen y `forzados == 0`.
   Lo que sigue abierto: el caso con UN solo monitor no se pudo medir desde
   aca. Si vuelve a verse, las lineas `!!! reencuadre forzado` y `render` de
   `analisis/.cache/gui_traza.log` dicen que tamanos tenia cada cosa y en
   que ventana; si no hay ninguna `!!!` y los graficos igual se ven cortados,
   entonces NO es el blit y hay que mirar la ventana (`ventana:` contra
   `pantalla`), no la figura. De yapa: `App.destroy` cancela los `after`
   periodicos, que sobrevivian a la ventana y ensuciaban la salida de los
   verificadores con `invalid command name ..._vaciar_cola`.

64. **Las figuras del informe se veian chicas, y no era la letra: era el
   factor de escala** (14-sep-2026). `analisis/informe.py` las generaba de
   15 in de ancho con letra de 8-10 pt, y `informe.tex` las mete a
   `\linewidth`, que con `a4paper` y `margin=2.4cm` mide 16.2 cm = 6.38 in:
   escala 0.43, o sea 9 pt que quedan en 4. Subir la fuente no lo arregla,
   lo corre. La regla es la inversa de lo intuitivo: **la figura se hace del
   ancho que va a tener en la pagina, con letra normal, y se incluye sin
   escalar**; asi 1 pt en la figura es 1 pt en el documento. Lo que quedo:
   - **Seccion 3 de `informe.py`**: `ANCHO_TEXTO` derivado de los mismos
     numeros del `preambulo.tex` (A4, margen 2.4 cm) y `CUERPO_PT = 11`;
     `ESTILO` con rotulos de ejes a 10 pt (el `\small` de los epigrafes),
     ticks y leyendas a 9 (`\footnotesize`), nada bajo 8; `figura(ancho,
     alto)` en fraccion de `\textwidth`, `guardar` en PDF sin
     `bbox_inches="tight"` (recorta el margen y LaTeX vuelve a estirar), y
     tres ayudas: `etiqueta_panel` ("(a)" adentro, porque el titulo es el
     epigrafe), `pocos_ticks` y `leyenda_abajo` (`fig.legend(loc="outside
     lower center")`, lo unico que cabe con paneles de 2-3 in y 8 series).
   - **La misma tipografia que el documento.** El `.tex` usa `lmodern`;
     MiKTeX trae las Latin Modern en OTF pero no las instala en el sistema.
     `fuentes_del_documento` se las pide a `kpsewhich` y las registra en
     matplotlib (`font_manager.fontManager.addfont`), con `mathtext.fontset
     = "cm"` para la matematica. ⚠ `kpsewhich` escribe UTF-8: con
     `text=True` la ruta con acentos (`Física Nuclear`) sale rota y "no
     existe"; hay que decodificar los bytes a mano.
   - **`ajuste_*`, `ganancia_dv` y `cociente_fuentes` ya no se copian ni se
     recortan de PNG de 15 in**: se dibujan aca. El ajuste se reconstruye de
     lo guardado (`.cache/calibracion_<clave>_edep.json` + espectro cacheado
     → `modelo.Problema` → `prediccion(calibracion.parametros(r))`), sin
     recalibrar: el recocido es estocastico y la figura tiene que mostrar EL
     ajuste de las tablas. Verificado: chi2 y A reconstruidos coinciden con
     los guardados a todos los decimales en los cuatro detectores.
   - En `informe.tex` solo cambiaron las lineas `\includegraphics`: `.png` →
     `.pdf`, y el apendice de `0.9\linewidth` a `\linewidth` (los `ajuste_*`
     se generan a ancho completo). **El `width=` del `\includegraphics` y el
     `ancho=` de `figura()` tienen que coincidir**: `cociente_fuentes` va a
     0.72 en los dos.
   - Dos decisiones de dibujo que no estan en la GUI: el eje vertical de los
     `ajuste_*` se corta en 1.5 veces el maximo DENTRO de la ventana de
     ajuste (el pico del exceso, fuera de la ventana, llega a 5-6 veces eso en
     el Co-60 y aplastaba el flanco), y las leyendas van debajo de los
     paneles o se sacan cuando el epigrafe ya lo dice.
   Con eso `python analisis/informe.py` regenera las 10 figuras en PDF y
   compila (27 paginas); `informe/README.md` tiene la tabla de tamanos y las
   reglas para no re-romperlo.
   Y una de contenido, no de tamano: la figura 5 (`comparacion_por_punto`)
   decia "cuentas (norm.)" sin decir contra que, y no explicaba por que las
   lineas llenas no caen sobre las punteadas. Ahora el eje dice
   "cuentas / plateau" (cada espectro dividido por la mediana de sus cuentas
   en la mitad inferior de su ventana de ajuste), el texto de cada panel
   lleva "desac." como la tabla del test A, y el epigrafe y la subseccion
   "Como se lee la figura" (`sec:leer-porpunto`) dicen la mecanica: la llena
   esta en E*D_fuente/D_conjunta y la punteada en E, asi que coinciden solo
   si la fuente comparte la D conjunta; como esta es un compromiso entre las
   dos, cuando difieren una se corre a cada lado (Cs a la derecha y Co-60 a
   la izquierda en A y D, al reves en C). Trampa de LaTeX: con babel-spanish
   `$\sim 1\,\%$` da "Incompatible glue units"; el `\%` va fuera de `$`.

65. **La "semilla" del NCC NO es el punto de partida del recocido, y el
   informe lo decia como si lo fuera** (14-sep-2026). Al escribir en detalle
   que entrega el NCC (seccion 3.5 del informe, parrafos "Que entrega, paso a
   paso", "Como entra eso al recocido" y "En numeros") aparecio el matiz:
   `calibracion.calibrar` llama `sa.multiarranque(p.costo, cotas, n=8,
   semilla=7)` SIN `x0`, asi que los 8 arranques parten de puntos uniformes
   dentro de las cotas, a proposito (si todos arrancaran de D0, que terminen
   juntos no probaria nada). Lo que D0 = B/E_ref fija es (a) DONDE se busca,
   `COTA_D_REL = (0.6, 1.6)` veces D0; (b) el paso, 0.1 del ancho de la cota;
   (c) la referencia contra la que se mide cuanto corrigio el ajuste
   (`semilla_D`, `chi2_semilla`). sigma_rel y lambda no dependen de la
   semilla: `COTA_SIGMA_REL = (0.02, 0.30)`, `COTA_LAMBDA_MIN = 10`,
   `LAMBDA_SEMILLA_REL = 0.15`. Esos literales estaban sueltos dentro de
   `calibrar`; ahora son constantes de `analisis/calibracion.py` y el informe
   los cita por macro (`\CotaDMin`, `\SigmaRelNCC`, `\NBarrido`, `\NMarkov`,
   `\Alfa`, ...), leidos de las constantes y de las firmas con `inspect`, no
   copiados. La corrida de ejemplo del informe (A @ 1600 V, Cs-137) tiene sus
   propias macros `\Ej*` (B, rho, ventana, D0, cotas, D final, desvio), y
   `informe.corrida_ejemplo` es la unica definicion de cual es, para la
   figura 1(c) y para el texto. Ojo: en la campana `gui.EXCESO = False`, o sea
   modelo de DOS parametros y ventana de ajuste = ventana del NCC; el
   `pipeline.py co60` suelto usa `calibracion.EXCESO = True`, y por eso su
   ventana arranca en 0.25*B y la del lote en 0.45*B.

66. **El informe tiene una sección propia de los pulsos y su integración, y
   escribirla desarmó un número que se venía leyendo mal** (14-sep-2026).
   `informe.tex` §3 «Del pulso al espectro»: qué registra el DRS4, un pulso
   real con cada paso de `pulso.area_evento` marcado (`integracion_pulso.pdf`),
   cómo se mide la compuerta y por qué es corta, el control de área/amplitud
   contra amplitud, y de las áreas al espectro (binning, tiempo vivo, resta,
   error, fin del umbral), más `tabla_pulsos` con eventos, descartados,
   pared, vivo, registro tras el onset y carga por corrida. Todo sale de
   `informe.pulsos_de` (ventana medida + áreas del caché, misma compuerta que
   el lote) y va como macro; `espectro.RANGO_FACTOR`/`RANGO_PERCENTIL` pasaron
   a constantes de módulo para poder citarlas. Lo que apareció:
   - **La columna «Carga» (72-82 %) no compara detectores: compara
     registros.** Se normaliza a la carga del pulso REGISTRADO, y el registro
     termina a 47 ns del onset en el D del 8-sep, a ~81 en la sesión del 1-sep
     y a 130-158 en los XML. Normalizando a la carga acumulada hasta 40 ns
     (`informe.T_COMUN`, dentro del registro más corto), la compuerta de 12 ns
     recoge **81-83 % en las 25 corridas**: la forma del pulso es la misma en
     los cuatro detectores (pico a 4.1 ns del onset, 7.1 ns FWHM) y la
     compuerta común mide en todas la misma porción del mismo pulso. Ojo al
     dibujar la carga acumulada normalizada al máximo: `pulso_promedio`
     rellena con ceros fuera del registro, así que el D «converge» a 45 ns
     por el registro, no por el pulso.
   - **Área/amplitud contra amplitud es plano al 1-4 % por encima de la
     amplitud mediana** (mediana 2.7 %, máximo 4.3 % en la campaña) y cae
     4-5 % hacia los pulsos más chicos, donde la amplitud (una muestra, ~1 mV
     de ruido) y el onset al 10 % se miden peor. Los Co-60, con amplitudes
     hasta 400 mV, son los más planos (1 %).
   - `pre` no es constante: `calibrar_ventana` lo mide y va de 2.0 a 8.8 ns
     (el 8.8 es `a_co60_1700`). Es zona de línea base, aporta ruido y no
     señal, pero el informe lo dice en vez de escribir «2 ns».
   - Los descartes por bandera son 0.5-2.5 % salvo el Co-60 de junio (6.2 %,
     pulsos al principio del registro cuya compuerta se sale por ese extremo).
   Herramienta: para parchear archivos con `\\n` literales hay que usar un
   script en disco; en un heredoc de la consola de este entorno las barras
   dobles llegan simples y el `replace` no encuentra nada.

67. **El informe tiene marco teórico de estadística, los clásicos y el NCC por
   dentro, y el cierre con verdad conocida de la cadena COMPLETA de
   `analisis/`** (20-sep-2026). Son los bloques B9-B12 de `informe/PLAN.md`
   (el plan del texto guía, que dice qué le falta al informe sesión a
   sesión). Tres subsecciones nuevas en `informe.tex` (`sec:clasicos` antes
   del NCC, `sec:estadistica` después del SA, más párrafos en el NCC y citas
   en el SA), tres figuras (`teoria_clasicos`, `teoria_ncc_detalle`,
   `teoria_estadistica`), dos tablas (`tabla_clasicos`, `tabla_cierre`),
   `referencias.bib` con 10 entradas verificadas y bibtex en `informe.compilar`.
   Lo que apareció:
   - **`analisis/` no tenía cierre con verdad conocida de la cadena**: el
     `registro.py --cierre` prueba la recta del test A con puntos fabricados,
     no el ajuste. `informe.cierre_conocido` lo hace de verdad: toma el ajuste
     real de la corrida de ejemplo (A a 1600 V, Cs-137) como verdad, fabrica
     el espectro que predice (fuera de la ventana lleva el dato suavizado,
     para que el NCC y `fin_del_umbral` vean lo mismo que en una corrida
     real), lo muestrea 30 veces con Poisson, y pasa cada realización por
     NCC → recocido → pulido → curvatura, dos veces: con las barras del dato
     (Neyman) y con las del modelo (Pearson). ~2.5 s por ajuste, 4 min en
     total, cacheado en `.cache/informe_cierre_<clave>.json` con la verdad
     como firma. Resultado: sesgo −0.04 ± 0.03 % (Neyman) y +0.00 ± 0.03 %
     (Pearson) con 292 cuentas/bin; a esta estadística Neyman no importa.
   - **La barra reportada queda ~15 % corta**: dispersión 0.17 % contra 0.15 %
     reportado, cobertura 18/30. `corr(D, lo) = 0.85` entre realizaciones: la
     ventana se elige de nuevo en cada una (su límite inferior se mueve
     ±1 mV·ns, o sea un bin) y la curvatura, que se calcula con la ventana
     FIJA, no puede ver esa dispersión. La figura lo muestra con dos elipses
     (la de Δχ² = 1 del ajuste real y la empírica del cierre). El error
     reportado en `tab:calib` viene además escalado por sqrt(χ²/ndf), que
     sobre datos reales lo agranda más que eso.
   - **`modelo.Problema(varianza="modelo")`** es el χ² de Pearson iterado
     (3 vueltas de proyección de A con la varianza del modelo). Default
     `"dato"`, sin cambio de resultados: `pipeline.py co60` sigue en
     1.2512 ± 0.0047. Y `ncc.buscar_intervalo(n_pts=)` para medir el efecto
     del remuestreo: entre 12 y 400 puntos B se mueve 0.3 % y ρ va de 0.9828
     (12 puntos, inflado) a 0.9808 (un punto por bin).
   - **Los clásicos sobre las 25 corridas reales** (importando
     `codigo/etapas/clasicos.py` desde `informe.py`, con `codigo/` agregado a
     `sys.path`; el `import rutas` de ahí resuelve a `codigo/rutas.py` porque
     `analisis` se importa como paquete): mediana +11.8 % (50 % de altura),
     +10.2 % (derivada), +12.7 % (escalón gaussiano) contra la D del ajuste.
     Sobre el sintético de Klein-Nishina a la misma resolución dan 3-7 %. El
     informe atribuye la diferencia al scattering múltiple, que la plantilla
     tiene y Klein-Nishina no; es plausible pero NO está medido aparte
     (pendiente en `PLAN.md`, B12).
   - **El perfil de Δχ² coincide con la curvatura** en la corrida de ejemplo
     (0.00161 contra 0.00155): el valle es cuadrático. El escalado lleva la
     barra a 0.0022.
   Trampas de esta sesión: (a) las tablas anchas de 7 columnas se salen del
   texto y LaTeX sólo lo dice con un `Overfull \hbox` en el log; hay que
   mirar `grep Overfull informe.log` después de compilar (quedan 3, todas
   anteriores a esta sesión). (b) Una leyenda de matplotlib dentro de un
   panel de 2 in tapa lo que hay que ver; cuando el epígrafe ya dice qué es
   cada marca, se saca. (c) `\,` de LaTeX no existe en los rótulos de
   matplotlib fuera de modo matemático: sale el backslash literal.

68. **El marco teórico arranca ahora por la física** (21-sep-2026; bloques
   B1-B3 de `informe/PLAN.md`): `sec:interaccion` (los tres procesos, el
   plástico como bolsa de electrones libres, camino libre medio, por qué no
   hay fotopico), `sec:cinematica` (E' y T, E_max en 180°, la errata de la
   Ec. 1 con su número, el pico de retrodispersión, tabla `tab:cinematica`) y
   `sec:kn` (supuestos, Thomson, ligadura y Doppler, el jacobiano escrito).
   Lo que apareció:
   - **Las secciones eficaces son datos de NIST XCOM**, no números tipeados:
     `datos/xcom/xcom_C9H10.txt` y `xcom_NaI.txt`, descargadas con el
     formulario "Compound" (`POST` a `cgi-bin/Xcom/xcom3_2` con `Formula`,
     `Name` sin dígitos, `NumAdd=1`, `Output=on`; la tabla viene en el HTML,
     el botón "Download" de `data.pl` rechaza la consulta). Cabecera con
     URL, fórmula y fecha; `informe.xcom(nombre)` las lee. XCOM reparte
     C9H10 en H 0.0853 / C 0.9147 en peso, lo mismo que el material de
     `ComptonSimulation`. `xraylib` se probó y no sirve: llega a ~800 keV y
     no tiene producción de pares (quedó instalado, no se usa).
   - **Klein-Nishina de electrón libre × 3.26e23 electrones/g reproduce el
     "incoherente" de XCOM al 0.9 % entre 100 keV y 3 MeV**, y se separa
     hasta 35 % a 10 keV (ligadura). Es la verificación de que "libre y en
     reposo" vale acá, y está en la figura. La σ total analítica
     (`informe.kn_total`) coincide con la integral numérica de la Ec. 2 del
     paper a 1e-8; σ_T = 0.665 b, σ_KN(662) = 0.256 b.
   - **Números que salen solos y sirven de intuición**: en el plástico a
     662 keV el fotoeléctrico es 4e-5 veces el Compton (NaI: 0.13); el
     fotoeléctrico deja de dominar a 21 keV en el plástico y a 266 keV en
     NaI; camino libre medio 11.6 cm (662) y 16.2 cm (1332); probabilidad de
     interactuar en los 5.08 cm del cilindro simulado, 35 % y 27 %.
   - `informe.geant4_geometria()` lee densidad (1.032), Birks (0.126 mm/MeV),
     radio y semialtura (2.54 cm) del `DetectorConstruction.cc` con regex,
     para que el informe no los tipee. El tamaño del detector REAL sigue sin
     confirmarse (`PLAN.md`, sección 5).
   - Trampa de heredoc otra vez (§7.66): un patrón de `replace` con `\\`
     dentro de un `python - <<'EOF'` no coincide nunca; todo parche con
     barras va a un script en el scratchpad.

69. **El informe tiene sección propia de la simulación, y escribirla destapó
   que se estaba leyendo la simulación EQUIVOCADA** (21-sep-2026; bloque D de
   `informe/PLAN.md`). `informe.geant4_geometria` sacaba densidad, Birks,
   radio y semialtura de `JcV/ComptonSimulation/src/DetectorConstruction.cc`
   (Windows). Esa simulación escribe UN histograma en CSV, emite iones de
   Cs-137 con `G4RadioactiveDecayPhysics`, no tiene carcasa ni guarda
   `NCompton` ni `Lvis`: **no puede haber producido los `.root`**. La que sí
   los produjo es `~/Proyectos/Geant4-Setup/ComptonSpectrum` en WSL (Geant4
   11.4.0, según su `Geant4ConfigVersion.cmake`), y ahora está copiada en
   **`simulacion/`** (src, include, macros, CMake, README; `LEEME.md` dice
   qué hay y qué no). Los cuatro números coincidían en las dos, así que
   ningún resultado cambió; la descripción que se podía escribir, sí:
   - **La geometría real tiene carcasa de Al (pared 10 mm, tapa 6 mm hacia
     la fuente) y un cilindro de vidrio borosilicato de 3 cm detrás** como
     ventana de PMT, todo en aire; el plástico es `G4_PLASTIC_SC_VINYLTOLUENE`
     (el README de allá lo llama BC-408, este proyecto BC-404: misma base
     PVT, misma densidad y composición, así que para el depósito da igual,
     pero hay que confirmar cuál es el de verdad, `PLAN.md` §5). La fuente son
     **gammas monoenergéticos, una línea por evento sorteada por peso, haz
     `pencil` axial** a 10 cm de la tapa; no hay decaimiento ni betas. Los
     `.root` se corrieron con los defaults (los macros dicen `beamOn 2M`,
     pero el programa escribe el `beamOn` real en el nombre: 3M y 10M).
   - **El 2 % de eventos con `NCompton == 0` y sin fotoeléctrico son
     electrones Compton nacidos en la tapa de Al** que entran al plástico
     (23134 de 24212 en el Cs-137, mediana de Edep 275 keV). La carcasa no es
     decorativa y hay que decirlo en la tabla.
   - **El fotoeléctrico en el plástico sí hace fotopico, invisible**: el
     1.3 % de los eventos del Cs-137 tiene un `phot`, el 58 % de ellos tras
     UNA sola dispersión Compton (no "varias", que fue lo primero que se
     escribió), y son los 724 eventos de absorción total en 3.5 millones.
   - **El binning de la plantilla, medido en `analisis/`** reajustando A a
     1600 V con 100/175/350/700 bines: D se mueve **0.06 %** (1.3565 / 1.3565
     / 1.3563 / 1.3555). Es §7.26 de nuevo, ahora con la cadena vigente. Lo
     que sí se mueve es la barra por curvatura (0.0028 / 0.0031 / 0.0022 /
     0.0017): la plantilla se interpola linealmente y el paso de derivación
     ve los quiebres (§7.49). En `analisis/` NO aparece el peine de aliasing
     de §7.25 porque `plantilla_en_canal` interpola en vez de rebinnear.
   - Lo nuevo en `informe.py`: `simulacion_resumen` (estadísticas de los dos
     `.root`, cacheadas por tamaño + fecha), `binning_plantilla` (~10 s,
     cacheado con la D guardada como firma), `figura_simulacion` (Edep/Lvis
     del Co-60, Cs-137 por `NCompton`, quenching evento a evento, y la razón
     dato/modelo de la corrida de ejemplo en todo el rango: 1 en la ventana,
     1.7 por debajo, cero en el umbral), `tabla_simulacion` y ~50 macros
     `Sim*`/`Bins*`. Citas nuevas verificadas por CrossRef: `agostinelli2003`,
     `allison2016`, `birks1951`.
   Dos trampas de LaTeX: **los nombres de macro no admiten dígitos**
   (`\BinsD100` compila como `\BinsD` seguido de "100" y revienta con
   "Missing \begin{document}" dentro de `numeros.tex`, lejos de la causa), y
   una celda con `ComptonSpectrum_Cs137_10M.root` se sale 29 pt de la tabla
   sin más aviso que el `Overfull` del log (quedan 3, todos anteriores).

70. **El atenuador: 20 dB en TODAS las corridas, y por qué tenía que estar**
   (22-sep-2026, confirmado por el operador). Era la última incógnita de la
   cadena de señal (§7.45) y estuvo abierta un mes: el Co-60 de junio decía
   "att" en el nombre y de las 24 corridas restantes no constaba nada, así que
   `rutas.misma_cadena` avisaba "atenuador sin registrar" cada vez que dos
   corridas se iban a combinar.
   - **Para qué está.** El DRS4 recorta a **500 mV** (§7.33). Las amplitudes
     medidas del Co-60 a 1550 V tienen mediana 44.7 mV, p99 179.3 y máximo
     499.8 (o sea el 0.01% ya toca el techo). Multiplicadas por 10, que es lo
     que hay a la SALIDA del PMT, la mediana queda en 447 mV y el p99 en 1.8 V:
     **sin el atenuador el 47.3% de los eventos habría saturado**. No es una
     precaución nominal, es la diferencia entre medir el espectro y recortarlo
     por la mitad.
   - **Qué NO cambia: ningún número.** El área se integra a la ENTRADA del
     digitalizador y ningún módulo corrige por la atenuación, a propósito: el
     observable es lo que entra al DRS4. Verificado después del cambio:
     `pipeline.py co60` sigue dando **D = 1.2512 ± 0.0047** keVee con
     χ²/ndf = 1.169. Lo único que cambió es metadato (`Corrida.cadena()`, el
     JSON de `calibracion.py`) y que los avisos de `misma_cadena` se apagaron
     —por el motivo correcto: ahora el dato está, antes faltaba—.
   - **Cómo leer las D de este archivo.** Todas son "mV·ns a la entrada del
     DRS4 por keV". La carga del PMT es 10 veces mayor (20 dB = factor 10 en
     amplitud, y el área escala igual). Si alguna vez se compara con otro
     montaje, o se quiere la ganancia del tubo, ese factor hay que ponerlo a
     mano.
   - **De regalo, un control cruzado.** §7.54 había DEDUCIDO que la atenuación
     no cambió entre junio y septiembre, porque el punto de junio cae sobre la
     curva `D ∝ V^9.49` a −0.08%. Quien midió dice lo mismo, y las dos cosas
     son independientes: el ajuste no usó el dato del atenuador. Si hubiera
     cambiado, ese punto habría caído un factor 10 afuera.
   → `rutas.Corrida` lleva `atenuador=20` en las 25 corridas;
   `rutas.misma_cadena` se conserva tal cual para la corrida que algún día
   venga con otro valor.

71. **Qué es del paper y qué es nuestro** (22-sep-2026). Está acá porque la
   regla de §4 —marcar **[fuera del paper]**— se fue aplicando entrada por
   entrada y quedó dispersa; con 71 hallazgos ya no se puede reconstruir de
   memoria cuál es cuál. El detalle, con las citas literales, vive en
   `replica_paper/README.md`; esto es el índice.

   | | **Kong et al. 2025** | **nuestro** |
   |---|---|---|
   | observable | arranca en `(canal, cuentas)` — la Fig. 2 empieza ahí | área con compuerta anclada al onset (§7.33), tiempo vivo por scaler (§7.34), lector del binario del DRS4 (§7.51) |
   | plantilla | Klein-Nishina analítica (Ec. 2 y 10) | forward folding con Geant4 (§7.18), mezcla de scattering `w` ajustada (§7.35), escalas keVee/keV por Birks |
   | mapa canal↔energía | por el ORIGEN (Ec. 8) | igual, más el test de la ordenada como cota (§7.50) |
   | selección de intervalo | NCC/Pearson (Ec. 15), inicio y fin libres, todo reinterpolado a 100 puntos | barrido 1D por el origen (§7.10) y un punto por bin (§7.27) |
   | optimizador | SA + Metropolis, `T_init=2000`, `α=0.95`, 100 iteraciones por T, paso `U(−0.1δ, 0.1δ)` | `T_init="auto"` y corte relativo (§7.22), VARPRO por nnls y multiarranque. El **paso sigue siendo el fijo del paper** en `analisis/sa.py`: se midió contra el adaptativo y dan lo mismo (1.2219 ± 0.0042 contra 1.2197 ± 0.0046), así que no hay evidencia para cambiarlo. En `codigo/etapas/sa.py` sí es adaptativo por default (§7.21) |
   | función de costo | mínimos cuadrados SIN pesos (Ec. 14) — o el error NCC, según el párrafo | χ² ponderado con las barras del espectro neto (§7.19). La varianza por default es la del DATO (Neyman) en `analisis/`; la de Pearson (§7.31) está disponible y se usa para medir su propio sesgo. En `codigo/` el default es Pearson |
   | qué se ajusta | {A, B, σ} o Θ de la Ec. 11 | {A, D, σ_rel} más λ del exceso de bajo canal, con A y C proyectadas (§7.32, §7.35). Son 3 y no 6 porque la plantilla simulada fija A₂/A₁, B₂/B₁ y la mezcla de scattering |
   | incertidumbre | SD entre 10 mediciones repetidas | curvatura del χ² escalada por √(χ²/ndf) y presupuesto de sistemáticos medido variante a variante (§7.41, §7.47). El MCMC es de `codigo/` (Etapa 8); **`analisis/` no tiene MCMC** |
   | validación | recta canal-energía con SD, R² y similitud coseno | los dos tests separados: eje (A) y ganancia `D ∝ V^b` (B) (§7.46, §7.54) |
   | comparación | FWHM gaussiano y mínimo de la derivada | los mismos, medidos contra verdad conocida (§7.44) |

   ⚠️ La columna derecha **no es una sola implementación**: `codigo/` (Etapas
   0-10) y `analisis/` (Etapa 11) tomaron decisiones distintas en el
   optimizador, en la varianza del χ² y en el posterior. Cuando se cite un
   default hay que decir de cuál de las dos, porque no coinciden.

   Y dos cosas de `analisis/` que no entran en ninguna fila porque el paper no
   las trata en absoluto: el **piso de la ventana por el umbral del trigger**
   (`espectro.fin_del_umbral`, §7.48) y el **escalado del error por
   √(χ²/ndf)**, que es práctica estándar y no está en Kong et al.

   La regla práctica: **si algo de esta columna derecha hay que defender, se
   defiende con una medición de este proyecto, no citando al paper.** Y al
   revés, `replica_paper/` no puede incorporar nada de la derecha; ese es el
   punto de que exista.

72. **El Co-60 del detector B, y el primer detector que se porta bien**
   (22-sep-2026). Llegaron las tres corridas que faltaban (`datos/Detector B/
   {1600,1650,1700}V/Co60`, binario, CHN1, 80000 eventos cada una, medidas hoy
   13:34, 13:46 y 13:58). El catálogo pasa de 25 a **28 corridas**.
   Lo primero, por el incidente del 10-sep: se verificó que los tres archivos
   son **únicos byte a byte** —no son copias de los del detector C ni del
   Cs-137 del propio B— antes de meterlos al catálogo.
   - **Un caso especial que habría calibrado la corrida equivocada.** El bucle
     de `_SESION_BIN` tenía `_a = "B1700V" if (_det == "B" and _hv == 1700)`,
     porque el Cs-137 del B a 1700 V se llama así en el disco. Sin agregarle la
     fuente a la condición, el Co-60 del B a 1700 V se habría resuelto al
     archivo del **Cs-137** y habría salido un número perfectamente formateado
     de la corrida equivocada, sin que nada fallara. Es la misma familia que
     §7.36 y que las copias del 10-sep: lo que no se verifica, se cachea.
   - **Resultado** (escala Edep, compuerta común de 12 ns):

     | HV | Co-60 | Cs-137 (1-sep) | Cs/Co |
     |---|---|---|---|
     | 1600 V | 1.0033 ± 0.0017 | 1.0717 ± 0.0021 | 1.0682 ± 0.0028 |
     | 1650 V | 1.3472 ± 0.0021 | 1.4364 ± 0.0032 | 1.0663 ± 0.0029 |
     | 1700 V | 1.7780 ± 0.0033 | 1.9140 ± 0.0047 | 1.0765 ± 0.0033 |

     χ²/ndf 2.59 / 2.35 / 2.04, FWHM 32.4 / 31.1 / 32.0%, en línea con los
     otros Co-60 de la campaña.
   - **TEST B: `b = 9.44 ± 0.04`** con χ² = 4.6 (1 gdl) y residuos de ±0.25%,
     dentro del 6-10 de un PMT. Y acá está lo bueno: **el B es el único
     detector cuyas dos fuentes dan el mismo `b`**:

     | det | b (Co-60) | b (Cs-137) | diferencia |
     |---|---|---|---|
     | A | 9.37 ± 0.03 | 10.42 ± 0.04 | 21 σ |
     | **B** | **9.44 ± 0.04** | **9.56 ± 0.05** | **1.9 σ** |
     | C | 8.70 ± 0.05 | 9.34 ± 0.05 | 9.1 σ |
     | D | 10.10 ± 0.04 | 7.73 ± 0.05 | 37 σ |

   - **TEST A: sigue sin cerrar, pero el B falla DISTINTO que el A, y eso es
     información.** El cociente Cs/Co tendría que dar 0.9799 en escala Edep (es
     el quenching de Birks) y da 1.070. El desacuerdo es del mismo tamaño que
     el del detector A (§7.55), pero con otra forma: en el A **crece** con el
     voltaje (1.075 → 1.124, χ² = 108 contra constante) y en el B es **plano**
     (1.0698 ± 0.0017, χ² = 5.9 con 2 gdl). Un cociente plano es exactamente lo
     que deja una diferencia de ganancia FRACCIONAL fija entre las dos
     sesiones: 9.17% de exceso, que con b = 9.5 equivale a **+0.93% de HV, unos
     15 V**, del mismo orden que los 16-24 V que §7.55 dedujo para el A.
     Y el C se va para el **otro lado** (0.87-0.90), lo que descarta que sea un
     sesgo del análisis a favor de una fuente: un sesgo de plantilla tendría el
     mismo signo en los cuatro.
     Tampoco es la compuerta: rehaciendo el test A con la de 75 ns en vez de la
     de 12, el cociente del B se mueve **+1.07%, +0.80% y −0.08%** y sigue
     plano (1.0796 / 1.0748 / 1.0756). Contra un desacuerdo del 9%, eso no
     alcanza ni de lejos.
     → **El Co-60 del B es del 22-sep y su Cs-137 del 1-sep**, así que esto NO
     cierra el test A: sigue haciendo falta **las dos fuentes en la misma
     sesión**. Lo que sí cierra es la otra limitación, la de que el B tenía una
     sola fuente y no se le podía ni intentar.
   - **Dos cambios de código, ninguno numérico.** `rutas._SESION_BIN` tiene una
     fila más y una nota base por fila (la nota común cita la bitácora del
     1-sep, que no vale para una corrida de hoy); y `calibrar_todo.calibrar_una`
     guarda ahora la cadena de SEÑAL (`atenuador`, `cadena`, `umbral`,
     `muestreo`) en `resultados.json`, que decía "la cadena de medición, sin la
     cual D no significa nada" y guardaba sólo la compuerta: el atenuador
     quedaba en `null` aunque el catálogo ya lo tuviera (§7.70).
   - **Control:** recalibrar las 28 con `--rehacer` (2.2 min) reprodujo las 25
     anteriores **exactamente**, D por D, sin una sola fila distinta. El SA está
     sembrado (`semilla=7`), así que eso es lo que tiene que pasar; si alguna
     vez no pasa, hay un bug.
   - **El informe ya está regenerado** (`python analisis/informe.py`, 50
     páginas). Las macros se movieron solas: `\NCorridas` 25 → 28,
     `\NPuntosTestA` 9 → 12, `\NFigurasDiag` 150 → 168, y ninguna desapareció
     ni apareció, o sea que ninguna referencia del texto quedó colgada. Lo que
     sí hubo que tocar a mano son cinco frases que el generador no puede saber:
     la bitácora de sesiones (faltaba el 22-sep), el epígrafe que daba por
     hecho que la fila del B era la de "una sola fuente: sin test", las dos
     enumeraciones de detectores ("A y D" → "A, B y D"; "en A y en C" → "en A,
     en B y en C") y la limitación del detector B, que dejó de serlo. Y "los
     siete grupos" del test B pasó a `\NGruposTestB`, que es lo que la regla
     del informe pedía desde el principio: un número que cambia con los datos
     no se escribe en el párrafo.

---

## 8. Estado del proyecto (roadmap)

- [x] **Etapa 0** — `spectrum.py`: XML→espectro neto (binning común + baseline + norm. temporal).
- [x] **Etapa 1** — `physics.py`: E_max + lineshape Klein-Nishina (Ec. 10). Validado (963/1118 keV).
- [x] **Etapa 2** — `model.py`: forward model (convolución gaussiana). Muestra el sesgo del 50%.
- [x] **Etapa 3** — `fit_single.py`: ajuste de un borde con least_squares. Funciona pero FRÁGIL.
- [x] **Etapa 4** — `ncc.py`: selección automática del intervalo (Sección II-B, Ec. 15).
      Sintético (verdad conocida): B = 600.0 vs 600 real, ρ = 0.9986.
      **Datos reales de Co-60** (`calibrar_ncc.py`, 800 bins de 10 mV·ns, rango 0-8000):
      intervalo automático [1960, 7886], ρ = 0.9918, σ_rel = 0.18.
      → **B = 4375 mV·ns ↔ 963.4 keV**, ganancia **4.54 mV·ns/keV** (220 keV por 1000 mV·ns).
      NCC solo daba 4357 y el ajuste sembrado 4375 (0.4% de diferencia); el método clásico del
      50% da 4615 (sesgado). Número SESGADO por el segundo borde (§7.5): la Etapa 6 lo corrigió
      a 4185 mV·ns (-5.8%).
- [x] **Etapa 6** — `fit_double.py`: dos bordes solapados (Θ={A1,A2,B1,σ1,σ2,D}, Ec. 11-13).
      Sintético (verdad conocida): B1 -0.6%, B2 -0.4%, σ1 -2.9%, **D +0.62%**, offset ~0.
      **Datos reales de Co-60:** B1 = 4179 y B2 = 4865 mV·ns, σ1/B1 = 0.177, A2/A1 = 0.75.
      → **CALIBRACIÓN: D = 4.344 mV·ns/keV** (230.2 keV por 1000 mV·ns), con `ligar_origen=True`
      porque D libre está mal condicionado (§7.16). Puntos: (4185, 963.4) y (4857, 1118.1).
      El ajuste de UN borde daba 4435 mV·ns para 963 keV: **-5.8% de sesgo**, confirmado.
- [x] **Etapa 7** — Calibración vía forward folding con la simulación Geant4
      (`sim_template.py`): D_luz = 4.23 ± 0.10 mV·ns/keVee, D_edep = 4.05 ± 0.10,
      σ_E = 5.73·√E keV, χ²/ndf = 0.96 con 597 gdl. Ver §7.17-7.19.
      ⚠️ **Superada por la Etapa 8**: esos números salían con la ventana (2000, 8000) tipeada a
      mano y con χ² de Neyman. Los dos efectos se corrigieron y D subió un 2.2%.
- [x] **Etapa 8** — `mcmc.py` + pestaña 8: la cadena **NCC → SA → MCMC**, sin números a mano.
      El NCC elige la ventana correlacionando contra la simulación (§7.28), el SA confirma que
      el óptimo es global (§7.30) y el MCMC da el posterior (R̂ = 1.0001, n_eff ≈ 4000,
      aceptación 0.32). En el camino aparecieron los dos sistemáticos que faltaban: el sesgo de
      Neyman (§7.31) y el exceso de cuentas a bajo canal (§7.32).
      **RESULTADO (ventana [2600, 7900], χ² de Pearson, χ²/ndf = 0.918 con 528 gdl):**

      | plantilla | D (mediana posterior) | α [√keV] | FWHM(963) |
      |-----------|----------------------|----------|-----------|
      | Geant4 Lvis → **keVee** | **4.320** −0.028 +0.028 | 5.626 | 42.7% |
      | Geant4 Edep → **keV** | **4.133** −0.027 +0.025 | 5.711 | 43.4% |
      | Geant4 Edep, 1 scattering | 4.368 ± 0.027 | 5.273 | 40.1% |

      → **D_luz = 4.320 ± 0.028 (est) ± 0.122 (sist) mV·ns/keVee** (231.5 keVee por 1000 mV·ns)
      y **D_edep = 4.133 ± 0.026 (est) ± 0.118 (sist)** (242.0 keV por 1000 mV·ns), con
      corr(D, α) = **−0.78**. El cociente D_luz/D_edep = 1.045 = 1/0.954 reproduce el Birks de
      la simulación, que es un control de consistencia independiente.
      **El sistemático es 4x el estadístico**: lo que limita la calibración es cuánta física le
      falta a la simulación, no cuántos eventos se midieron.
- [x] **Etapa 5** — `sa.py`: Simulated Annealing (Sección II-A). Con 3 parámetros le gana
      claramente al optimizador local (B = 600.2 vs 343 desde el mismo arranque malo). Con los
      6 de la Ec. 11 hay que proyectar las amplitudes (VARPRO) y correr multiarranque.
      **Su aporte real fue diagnóstico:** demostró que Θ={A1,A2,B1,σ1,σ2,D} está
      sub-determinado para bordes solapados y que el ajuste local de la Etapa 6 no estaba en el
      óptimo global (§7.21-7.23). Dos constantes del paper (T_init=2000, |ΔE|<0.001) resultaron
      no ser transferibles entre objetivos (§7.22).
      **REPLICA DEL PAPER COMPLETA:** Etapas 0-7, las 15 ecuaciones. **Cinco erratas**,
      todas verificadas contra el PDF y reproducibles con `replica_paper/auditoria.py`
      (acá decía "3 erratas" hasta el 22-sep-2026, y una de las tres no era errata sino
      inconsistencia interna):
      1. **Ec. 1** sin el `2α` del numerador: da el fotón retrodispersado. La propia Fig. 3
         rotula 477 keV para Cs-137, que es lo que da la corregida (la impresa daría 184.3).
      2. **Ec. 12** con `ΔE` "between gamma lines" (159.3 keV) en vez de entre BORDES
         (154.7). En su panel de LaBr₃ los dos bordes del Co-60 están a 74 canales: la
         corregida predice 73.8 y la impresa 76.0.
      3. **SD y R² intercambiados** en los recuadros de la Fig. 3, en los DOS paneles. Con
         la asignación correcta, la frase de la Sección III ("the standard deviation
         associated with the NCC-SA calibration consistently remains below that of the full
         energy peak fitting method") **dice lo contrario de lo que muestran sus datos**:
         2.183 contra 1.609 en LaBr₃ y 7.275 contra 1.981 en NaI.
      4. **Rótulos cruzados en el panel de NaI de la Fig. 3**: `Na22 (846, 662)` es Cs-137,
         `Cs137 (1050, 834)` es Mn-54, `Mn54 (1442, 1173)` es Co-60, `Co60 (620, 477)` es el
         borde del Cs-137, y tres más. Los pares canal-energía están bien, así que no cambia
         ningún ajuste; lo que está mal es a qué fuente se atribuye cada línea. El panel de
         LaBr₃ está bien.
      5. **Fig. 5**: la Sección III dice que la robustez se probó con **EJ309** y el epígrafe
         de la figura dice **BC408**.
      Y dos cosas que NO son erratas y conviene no venderlas como tales: la ambigüedad
      Ec. 14 (mínimos cuadrados) contra el párrafo de Metropolis ("NCC-based matching
      error"), que es una inconsistencia interna a resolver al implementar; y la **similitud
      coseno no verificable** (reportan 99.9993% y 99.9999% sin decir nunca entre qué
      vectores; cinco lecturas razonables no reproducen esos números, y la métrica es
      insensible por construcción: dos calibraciones con pendientes 3% distintas dan
      0.99998812).

- [x] **Etapa 9** — `spectrum._area_pulso_gate` + `tiempo_vivo` + `modelo_ext.py` +
      `calibrar_v2.py`: se arregló el OBSERVABLE (§7.33-7.34) y los dos sistemáticos que
      quedaban pasaron a ser parámetros ajustados (§7.35).
      **RESULTADO VIGENTE** (ventana [300, 1900] mV·ns, χ² de Pearson, χ²/ndf = 0.977 con
      314 gdl, posterior de 4 cadenas × 9000 pasos, R̂ < 1.007, n_eff ≈ 1300):

      | cantidad | valor | antes (Etapa 8) |
      |----------|-------|-----------------|
      | **D_luz** [mV·ns/keVee] | **1.261 ± 0.014 (est) ± 0.034 (sist)** | 4.320 ± 0.028 ± 0.122 |
      | **D_edep** [mV·ns/keV] | **1.208 ± 0.014 ± 0.032** | 4.133 |
      | σ_E(E) | 2.79·√E keV → **21.1% FWHM** a 963 keV | 5.63·√E → 42.7% |
      | w (fracción 1-scattering) | 0.796 ± 0.070 (simulación: 0.751 ✓) | no existía |
      | χ²/ndf | 1.008 (314 gdl) | 0.918 (528 gdl) |
      | error total | 2.9% | 2.9% |

      Los bordes caen ahora en **1215 y 1410 mV·ns**, y están a 1.8σ uno del otro (antes,
      0.93σ).
      Presupuesto de error: estadística ±1.07% (ya incluye w y el exceso vía corr(D,w) = 0.96),
      ventana ±0.75%, definición del área ±1.04%, **forma del exceso ±2.35%** ← el que manda.
      ⚠️ El total sigue dando 2.9% como en la Etapa 8, pero eso NO quiere decir que no se haya
      ganado: aquel 2.9% estaba sobre un valor central equivocado y **no incluía ninguno de
      los dos sistemáticos que hoy dominan** (la definición del área no se había mirado nunca,
      y del exceso solo se sabía que había que esquivarlo). Lo que cambió es que el
      presupuesto ahora está completo y se sabe qué medir para bajarlo.

**Resuelto (24-ago-2026):** se eligió la Etapa 4 (NCC), y de yapa resolvió el otro camino: la
semilla que entrega el NCC blinda el ajuste local sin necesidad de un barrido aparte (§7.8).

- [x] **Etapa 10** — la interfaz se configura por la MUESTRA (27-ago-2026, §7.41):
      `physics.FUENTES` con 9 presets + `Fuente.nota/avisos`, verificación de que la
      plantilla Geant4 sea de la fuente elegida (`sim_template.lineas_primarias`),
      `gui.ajuste_extendido` como receta única de la Etapa 9, y la pestaña ★ Calibración
      calculando sobre los datos cargados con el presupuesto de error medido y un registro
      por fuente (`.cache/calibraciones.json`). Sin regresión: la 10 sigue dando
      D = 1.2618 ± 0.0137 y `calibrar_v2.py --rapido` no se movió.

- [x] **Etapa 11** — `analisis/`: reimplementación limpia de toda la cadena, sin
      importar nada de `codigo/`, y los dos tests de linealidad puestos donde se
      pueden correr (§7.45-7.50). El punto de operación (HV, atenuador, umbral)
      entra al catálogo y es la clave por la que agrupa el registro; el exceso de
      bajo canal se modela también acá, y con eso las dos implementaciones
      coinciden al 0.2%. **RESULTADO VIGENTE de `analisis/`:**

      | corrida | punto | escala | D | χ²/ndf | FWHM |
      |---|---|---|---|---|---|
      | Co-60 (junio) | A @ 1550 V | keVee | **1.2512 ± 0.0047** | 1.17 | 20.8% |
      | Co-60 (junio) | A @ 1550 V | keV | 1.1973 ± 0.0038 | 1.18 | — |
      | Cs-137 (26-ago) | A @ 1600 V | keVee | **1.9450 ± 0.0038** | 0.97 | 24.4% |
      | Cs-137 (26-ago) | A @ 1600 V | keV | 1.8304 ± 0.0035 | 0.92 | — |

      Los dos tests siguen **sin poder correrse**, y ahora se sabe exactamente
      por qué: hay una fuente por punto de operación y un voltaje por fuente
      (§7.45-7.46). Lo que falta es medición, priorizada en
      `docs/campana_linealidad.md`.

**Pendiente que SOLO se cierra midiendo (25-ago-2026):** qué es el exceso de cuentas a bajo
canal. Ya no obliga a cortar la ventana (se modela), pero su forma es el sistemático más
grande. Una corrida corta con la fuente apantallada, o movida de lugar, lo resuelve. Segundo
pendiente: una fuente más (Cs-137) para tener un segundo punto REAL de la recta (§9-bis.3).
Ese segundo punto ya tiene dónde vivir: el registro de la pestaña ★ Calibración guarda una
entrada por fuente y las superpone (§7.41). Lo que falta es la MEDICIÓN —y, para cerrarla
bien, una simulación Geant4 de esa fuente con esta geometría: sin ella la ruta que queda es
la analítica, con su +4.3% de sesgo conocido (§7.18, §7.38).

---

## 9. Cómo ejecutar (estado actual)

**Estructura de carpetas** (reorganizada el 27-ago-2026, §7.43):

```
Calibración CE/
├─ CLAUDE.md                   este archivo
├─ analisis_paso_a_paso.pdf    el documento del análisis (lo genera documentacion.py)
├─ analisis/                  ⭐ LA IMPLEMENTACIÓN VIGENTE (27/28-ago-2026)
│   ├─ rutas.py                catálogo de corridas CON punto de operación (HV, atenuador)
│   ├─ drs.py                  elige el lector por los primeros 4 BYTES (§7.51)
│   ├─ drs_xml · drs_bin       XML del DRSOSC · binario del DRS4, misma interfaz
│   ├─ pulso · area · espectro                archivo -> área -> espectro neto
│   ├─ revisar.py              revisión de una tanda antes de calibrar (§7.52)
│   ├─ visor.py                interfaz para VER lo calibrado (los dos tests)
│   ├─ fisica · simulacion · ncc · modelo · sa · calibracion    el ajuste
│   ├─ registro.py             calibraciones acumuladas + TEST A (canal vs energía)
│   ├─ voltaje.py              TEST B (ganancia vs alto voltaje)
│   ├─ pipeline.py  gui.py     la receta única, y la interfaz que la llama
│   ├─ calibrar_todo.py        calibra el catálogo entero -> figuras/calibraciones/
│   ├─ informe.py              ⭐ arma el informe: figuras A TAMAÑO DE INFORME y con
│   │                          su tipografía (§7.64), tablas y números -> informe/
│   └─ README.md               qué decide cada módulo y qué apareció al construirlo
├─ informe/                   ⭐ EL INFORME (LaTeX). informe.tex y preambulo.tex se
│   │                          editan a mano; generado/ (tablas, números) y
│   │                          figuras/ (PDF) los reescribe analisis/informe.py
│   ├─ PLAN.md                 ⭐ qué le falta al informe para ser TEXTO GUÍA,
│   │                          por bloques y sesiones (20-sep-2026); se lee antes de escribir
│   └─ README.md               qué es de quién, cómo regenerar, tamaños de figura
├─ replica_paper/             ⭐ RÉPLICA LITERAL del paper (1-sep-2026, §7.56)
│   │                          regla: acá solo entra lo que el paper dice.
│   │                          Klein-Nishina analítica, sin Geant4 ni MCMC.
│   ├─ ecuaciones.py           las 15 ecuaciones, y las dos versiones donde hay errata
│   ├─ ncc.py  sa.py           Sección II-B y II-A, con SUS constantes
│   ├─ convencionales.py       los dos métodos contra los que compara (Fig. 4)
│   ├─ metodo.py               el diagrama de la Fig. 2, de punta a punta
│   ├─ resultados.py           Sección II-C y III: recta, SD, R², coseno
│   ├─ auditoria.py            los números PUBLICADOS, contra sí mismos
│   ├─ datos.py                el único puente con analisis/ (archivo -> canal, cuentas)
│   ├─ correr.py               el informe completo -> informe.txt
│   ├─ gui.py                  ⭐ el paper PASO A PASO, calculando en vivo (§7.56)
│   ├─ verificar_gui.py       corre las 14 pestañas y las compara con informe.txt
│   └─ README.md               qué salió de auditar el paper, con números
├─ simulacion/                ⭐ LA FUENTE DE GEANT4 que produjo los .root de datos/
│                              (copia de ~/Proyectos/Geant4-Setup/ComptonSpectrum en
│                              WSL, §7.69). informe.py lee de acá la geometría. NO
│                              confundir con JcV/ComptonSimulation/, que es otra.
├─ codigo/                     la implementación anterior (Etapas 0-10), se conserva
│   ├─ rutas.py                DÓNDE está cada cosa (RAIZ, DATOS, FIGURAS, CACHE)
│   ├─ gui.py  calibracion.py  documentacion.py
│   ├─ etapas/    physics · model · fit_single · fit_double · ncc · sa
│   ├─ datos_io/  spectrum
│   ├─ ajuste/    sim_template · modelo_ext · mcmc
│   └─ scripts/   calibrar_v2 · calibrar_ncc · calibrar
├─ datos/                      los .xml del DRSOSC, sus cachés .npz y los .root
├─ figuras/                    todo lo que los módulos dibujan al correr solos
├─ docs/                       Paper ref.pdf, HANDOFF_implementacion.md,
│                              campana_linealidad.md (qué falta MEDIR),
│                              formato_binario_drs4.md (⭐ los datos SIN
│                              extensión: formato, lector y trampas, §7.51-53)
├─ verificaciones/             _verif_*.py, lecturaxml.py (histórico)
└─ .cache/                     preferencias de la GUI, cadenas MCMC, registro
```

Los comandos se corren **desde la raíz**, pero ya no hace falta: todo resuelve
sus rutas desde `__file__` (§7.43), así que funcionan desde cualquier carpeta.

```bash
pip install numpy scipy matplotlib uproot reportlab

# --- analisis/ : la implementación vigente ---
python analisis/pipeline.py co60               # LA CADENA COMPLETA por consola (~4 s con caché)
python analisis/gui.py                         # la interfaz que CORRE la calibración
python analisis/visor.py                       # la interfaz que la MUESTRA (4 pestañas:
                                               #   test A, test B, espectros, archivos)
python analisis/drs.py                         # qué formato tiene cada corrida (xml/binario)
python analisis/drs_bin.py                     # metadatos y pulsos de los binarios del DRS4
python analisis/revisar.py --todas             # revisión de la tanda antes de calibrar:
                                               #   compuerta, tiempo vivo, espectro, borde crudo
python analisis/calibracion.py                 # las dos corridas, las dos escalas, con figura
python analisis/registro.py                    # las calibraciones acumuladas y el TEST A
python analisis/registro.py --rehacer          # recalcula todo el catálogo y regraba
python analisis/registro.py --cierre           # el test A con verdad conocida
python analisis/sistematicos.py                # el presupuesto de error, medido variante a variante
python analisis/voltaje.py                     # el TEST B, D(V), con cierre y figura
python analisis/rutas.py                       # qué corridas ve, con su punto de operación
python analisis/calibrar_todo.py               # calibra el catálogo entero (resultados.json)
python analisis/informe.py                     # ⭐ figuras (PDF, a tamaño de informe) + tablas
                                               #   + números + pdflatex -> informe/informe.pdf
python analisis/informe.py --sin-pdf           #   lo mismo sin compilar

# --- replica_paper/ : la réplica literal del paper (§7.56) ---
python replica_paper/gui.py                    # ⭐ EL PAPER PASO A PASO, calculando en vivo
                                               #   14 pestañas: Mapa · Datos · Ec.1 · Ec.2/10 ·
                                               #   Ec.8-9 · Ec.15 (NCC) · II-A (SA) · Fig.2 ·
                                               #   convencionales · Fig.5 · II-C · auditoría ·
                                               #   Sec.III · cierre con verdad conocida
python replica_paper/auditoria.py              # los números PUBLICADOS, contra sí mismos
python replica_paper/correr.py                 # el informe completo (~8 min) -> informe.txt
python replica_paper/correr.py --rapido        # sin los 10 bloques (~1 min)
python replica_paper/correr.py --cierre        # dos fuentes al mismo punto, con D conocida
python replica_paper/ecuaciones.py             # redibuja la Fig. 1 del paper (valida la Ec. 10)
python replica_paper/metodo.py                 # la cadena de la Fig. 2 con verdad conocida
python replica_paper/verificar_gui.py          # que la GUI corra y NO invente numeros
python replica_paper/verificar_gui.py --rapido #   lo mismo sin las 4 pestanas lentas

# --- verificaciones/ ---
python verificaciones/_verif_gui_pestanas.py   # que las 8 pestañas de analisis/gui.py DIBUJEN
                                               #   (sin Tk: llama a cada figura directamente y
                                               #   deja los PNG en figuras/gui_*.png)
python verificaciones/_verif_gui_pestanas.py cs137_1600
python verificaciones/_verif_gui_encuadre.py   # que la VENTANA no recorte lo que dibuja
                                               #   (compara el tamaño del render con el de la
                                               #   imagen de Tk en cada pestaña; §7.62. Un
                                               #   recorte así no se ve en un PNG, solo en la
                                               #   ventana). Prueba primero el ARRANQUE REAL,
                                               #   maximizado y por mainloop, pasea por los
                                               #   monitores si hay más de uno, y después los
                                               #   3 tamaños en normal; exige además que
                                               #   `_CanvasEncuadrado` no haya corregido nada
                                               #   (§7.63)
python verificaciones/_verif_gui_encuadre.py cs137_1600
python verificaciones/_verif_binario.py        # que los DOS lectores del binario del DRS4
                                               #   den lo mismo, sobre un archivo fabricado
                                               #   con verdad conocida (no hace falta el
                                               #   disco externo)

# --- codigo/ : la implementación anterior ---
python codigo/scripts/calibrar_v2.py           # LA CADENA VIGENTE (Etapa 9), con MCMC (~1 min)
python codigo/scripts/calibrar_v2.py --rapido  # lo mismo sin MCMC (~10 s)
python codigo/ajuste/modelo_ext.py             # la escalera de modelos M0->M3, con figura
python codigo/gui.py                           # VISOR POR ETAPAS (lo mas comodo para explorar)
python codigo/documentacion.py                 # genera analisis_paso_a_paso.pdf (~10 s con cache)
python codigo/scripts/calibrar_ncc.py          # cadena de las Etapas 4-6 (numeros viejos)
python codigo/rutas.py                         # dice qué rutas ve y si los datos están
```

⚠️ **`codigo/` no corre hoy: no encuentra sus datos** (medido el 22-sep-2026). `codigo/rutas.py`
líneas 86-87 buscan `datos/A Co60 …xml` y `datos/A BG …xml`, pero al copiar la campaña
(§7.61) los XML de junio quedaron en **`datos/Detector A/`**, así que `calibrar_v2.py`
muere con `FileNotFoundError` antes de calcular nada. O sea que la frase "sin regresión:
`calibrar_v2.py --rapido` da 1.2618" que aparece en §7.41, §7.43 y §7.57 era cierta cuando
se escribió y **ya no se puede correr**. No es un problema de análisis —los archivos están
enteros y `analisis/` los usa— sino dos rutas de una línea cada una. Se deja anotado y sin
tocar: la implementación vigente es `analisis/`, y arreglarlo es cambiar `dato(...)` por
`dato("Detector A", ...)` en `codigo/rutas.py`.

⚠️ La primera corrida después del arreglo del área re-parsea los XML (~37 s) y escribe un
caché NUEVO, `*.areas_gate.npz`. El caché viejo `*.areas.npz` se conserva: es lo que usa
`metodo="legacy"` para reproducir el análisis anterior y ver el antes/después.

**`gui.py`** abre 14 pestañas (**Archivos** · Espectro · Klein-Nishina · Forward model ·
Un borde · NCC · Dos bordes · SA · Geant4 · NCC→SA→MCMC · **El área** · **Modelo ext.** ·
Calibración · ? Parámetros). Cada una tiene controles a la
izquierda, gráfico al centro y una explicación de qué mirar a la derecha. Las pestañas de SA,
Geant4 y NCC→SA→MCMC no calculan solas (tardan): hay que apretar "Recalcular". Casos
preparados para ver los problemas conocidos: RANGO_MAX=46000 en la pestaña 0 (§7.7), "pad"
destildado en la 2 (§7.14), modo "libre" en la 4 (§7.10), las ligaduras destildadas en la 6
(§7.23), y en la 8 destildar "Pearson" (§7.31) o bajar `inicio_rel` a 0.45 (§7.32).

La pestaña **Archivos** (la primera) es de dónde salen los datos y QUÉ se midió: elige el
`.xml` de la muestra, el `.xml` del fondo, el `.root` de la simulación (**opcional**) y la
**fuente** (Co-60 / Cs-137 / manual con una o dos energías), con un selector de TTree que se
llena leyendo el propio archivo (`inspeccionar_root`, solo cabeceras) y avisa si faltan las
ramas `Lvis_keVee` / `Edep_keV` / `NCompton`. "Aplicar y recargar" borra los cachés EN MEMORIA
de `Datos`, marca todas las demás pestañas como pendientes y deja una vista previa del
espectro neto con el diagnóstico: pulsos, tiempo de pared vs tiempo VIVO, `k_m`/`k_f`,
percentiles del área y el RANGO_MAX sugerido (1.6·p99). La selección se recuerda en
`.cache/archivos.json`. Detalles que importan:
- Los cachés de áreas (`*.areas_gate.npz`) viven al lado de cada XML y llevan la huella del
  archivo, así que cada corrida tiene el suyo y no hace falta invalidarlos a mano; si el XML
  elegido no tiene caché hay que tildar "parsear los XML" (~40 s la primera vez).
- La clave del caché de MCMC de la pestaña 8 lleva ahora `Datos.huella()` (nombre + tamaño +
  fecha de los tres archivos): sin eso, cambiar de datos y reusar un `.npz` viejo daría un
  posterior de otro experimento. Cambiar de archivos deja obsoletas las cadenas ya guardadas.
- Los números que citan las explicaciones (D = 1.276 mV·ns/keVee, bordes en 1229 y 1427, los
  defaults de RANGO_MAX = 2500 e `inicio_rel`) son de ESTE detector: con otros archivos hay
  que revisarlos.
- **Sin `.root`** (botón "Quitar la simulación") la GUI queda en la ruta analítica: las
  pestañas 4 y 6 calibran igual con Klein-Nishina, y 7, 8 y 10 avisan que falta la
  simulación en vez de fallar. Con `.root` cargado, el número bueno es el de la 10, que le
  mide el sesgo a la plantilla analítica (§7.18, §7.35).
- La **fuente** manda sobre qué modelo se ajusta (§7.39): con dos líneas, bordes solapados;
  con una, un solo borde. La solapa de la pestaña 6 cambia de nombre para decir cuál está
  activo. Las ventanas de ajuste ya no son constantes: salen del NCC
  (`Datos.borde_estimado` / `Datos.ventana_ajuste`) y se adaptan al detector. Hay 9 presets
  (Am-241, Co-57, Ba-133, Cs-137, Mn-54, Na-22, Co-60, K-40, Tl-208) más "manual" con una o
  dos energías; el informe muestra la nota de la fuente y los avisos derivados de sus líneas
  (§7.41).
- **La simulación tiene que ser de la misma fuente** (§7.41). La pestaña lee los gammas
  primarios del `.root` (`E0_keV`) y los compara con la fuente elegida: si no coinciden, lo
  dice con «⛔» y las pestañas 7, 8 y 10 quedan bloqueadas con el aviso limpio de
  `SinSimulacion` en vez de devolver un D equivocado. El checkbox «usar la simulación igual»
  destraba el bloqueo a propósito.

Las dos pestañas nuevas de la **Etapa 9**:
- **9 · El área** muestra el hallazgo de §7.33 en cuatro paneles: la forma media del pulso con
  las dos líneas base, la distribución de área/amplitud (donde se ve el 26% con el signo dado
  vuelta), el jitter del trigger contra la ventana de línea base vieja, y el espectro antes y
  después. Tarda ~1 s cada 2000 eventos porque relee los pulsos crudos.
- **10 · Modelo ext.** corre el modelo de §7.35 con w y el exceso libres. Los controles
  permiten fijar w en 1 o en 0.751 para ver, en vivo, lo que el análisis viejo llamaba
  "sistemático de plantilla", y tildar MCMC para ver el posterior D–w.

La pestaña 0 tiene ahora un selector **`gate` / `legacy`**: `legacy` reproduce el área rota, y
con RANGO_MAX = 8000 se ve el espectro viejo superpuesto al nuevo.

La pestaña **★ Calibración** ya no es un cartel: **corre la cadena entera sobre los archivos y
la fuente que están cargados** (§7.41). No calcula sola (tarda ~2 min con el barrido de
sistemáticos, 20 s sin él): hay que apretar "Recalcular". Elige la ruta según lo que haya —
modelo extendido con Geant4 si la simulación corresponde a la fuente, Klein-Nishina si no— y
arma el presupuesto de error MIDIÉNDOLO: barre el límite inferior de la ventana, cambia la
forma del exceso (exponencial ↔ ley de potencia) y fija w en 1 y en el de la simulación. El
único término arrastrado es la definición del área (±1.04%, §7.35). Cada resultado se guarda
en `.cache/calibraciones.json` con la clave `fuente|huella`, y el panel derecho superpone las
fuentes ya medidas: **la dispersión entre ellas es el test de linealidad**. Controles: escala
principal (Lvis→keVee o Edep→keV), ventana auto/manual, "barrer sistemáticos", "Exportar…" y
"Vaciar registro".

La pestaña **8 · NCC→SA→MCMC** encadena las tres etapas sobre la misma ventana: NCC con
plantilla Geant4 elige el intervalo y la semilla, el SA busca el óptimo global de {A, D, α},
y el MCMC (Metropolis-Hastings propio, `mcmc.py`) da el posterior con R̂, n_eff y la
correlación real. Corre las tres plantillas y las superpone en el plano D-α: **la distancia
entre las tres nubes es el error sistemático, y es ~6 veces el ancho de cada nube.** Tarda
~5 min la primera vez; las cadenas se cachean en `.cache/mcmc_*.npz`.

La pestaña **"? Parámetros"** es el glosario: qué es B, σ, D, ρ, Θ, α, ΔE, χ²/ndf, keVee vs
keV, etc., con unidades, el valor medido en este proyecto y dónde vive en el código. La
búsqueda normaliza griegas, así que escribir "sigma" encuentra las entradas escritas con σ
(ver `_ALIAS`/`_normalizar`; ojo que `_normalizar` pasa a minúscula ANTES de sustituir, por
eso las griegas van listadas en minúscula). Buscar "erratum" lista las erratas del paper.

**La barra de abajo (§7.42).** Al lado del texto de estado hay una barra de progreso.
Arranca indeterminada (va y viene: el programa está vivo) y se pasa sola a determinada,
con porcentaje y tiempo restante, en cuanto el paso que corre sabe cuánto falta. Los que
lo saben: el parseo de los XML (mide bytes leídos, así que es exacto), la pestaña
★ Calibración (cuenta los ajustes antes de empezar: `[3/10]`) y los MCMC de las pestañas
8 y 10 (cadena y paso). El ETA no aparece hasta tener algo de historia (3% y 2 s), porque
antes es una regla de tres sobre ruido.

⚠️ **Nitidez en Windows (§7.24):** sin `SetProcessDpiAwareness`, Windows renderiza la ventana
a 96 DPI y la estira como bitmap. Esta pantalla está a **120 px/pulgada (escalado 125%)**, así
que se veía borrosa. `activar_dpi_awareness()` se llama ANTES de crear la ventana Tk; después
`configurar_nitidez()` ajusta `tk scaling` y `figure.dpi`. Hay un selector de nitidez abajo a
la derecha que repinta en caliente (~0.8 s) sin recalcular, usando el `_ultimo` resultado
guardado en cada pestaña.

La primera corrida parsea los XML (~32 s en total) y cachea las áreas; las siguientes arrancan
al instante. Salida: diagnóstico del espectro, intervalo elegido por NCC, barrido de
sensibilidad en `sigma_rel`, punto de calibración, y `calibracion_ncc.png` (4 paneles: espectro
log, intervalo automático, S1-vs-S2 estandarizados, y el ajuste final).

Demos por etapa, cada una standalone (las figuras van a `figuras/`):
`python codigo/etapas/physics.py`, `.../model.py`, `.../fit_single.py`,
`.../ncc.py`, `.../fit_double.py`, `.../sa.py`, `python codigo/ajuste/sim_template.py`,
`.../mcmc.py`, `python codigo/calibracion.py`.

`python codigo/ajuste/mcmc.py` es el test con VERDAD CONOCIDA de la Etapa 8: genera datos con la
simulación a D = 4.20 y α = 5.70 con ruido de Poisson a la estadística real de este
experimento (26 cuentas/bin), y verifica que los tres parámetros caen dentro del 68% del
posterior. Ahí se ve también el sesgo de Neyman en vivo (§7.31).

---

## 9-bis. Sobre `HANDOFF_implementacion.md`

Documento de plan para las etapas 4-7, escrito contra el estado ANTERIOR (las da por
pendientes). **Ya está todo implementado**, así que es histórico: leerlo como registro de
decisiones, no como pendientes.

Se cumplieron D2 (Savitzky-Golay solo para NCC), D3 (siempre dos bordes), D4 (ΔE de bordes)
y D5 (semilla D0 = B1/E_max,1). Cuatro desviaciones, todas forzadas por los datos:

1. **Etapa 4** — el plan pide `seleccionar_intervalo(..., ancho_ventana=...)` con el ancho
   libre. Sobre datos reales eso elige ventanas degeneradas (§7.10). Se implementó
   `ncc.buscar_intervalo(modo="origen")`, que barre solo B.
2. **D1** confirmada y reforzada: con el objetivo NCC las amplitudes no son solo poco
   robustas, son INDETERMINABLES (§7.22).
3. **D6 subestima el problema**: no son "2 puntos que determinan una recta exacta". Los
   bordes están a 0.93σ (§7.16), así que la separación interna no da un segundo punto real.
   El Cs-137 hace falta para TENER un segundo punto, no solo para validar linealidad.
4. **El plan no contempla la simulación Geant4** (§7.17-7.18). Todo su §3 y §5 asume
   plantilla analítica; el forward folding la supera y el flujo de su §5 quedó obsoleto.
5. **Ni el plan ni el paper miran el OBSERVABLE.** Los dos dan por buena el área del pulso y
   se concentran en cómo ajustar el borde. Acá el error dominante estaba en el área (§7.33) y
   en el tiempo vivo (§7.34), o sea aguas arriba de todo lo que el plan discute. Moraleja para
   la próxima etapa: antes de refinar el ajuste, mirar unos cuantos pulsos crudos.

Nombres: el plan pide `optimize.py`, `calibrate.py` y `model.forward_model_doble`; quedaron
como `sa.py`, `calibracion.py` y `fit_double.modelo_dos_bordes`.

---

## 10. Referencia

Kong, Zhang, Lin, Qiu, Chen, Wang. *The Automatic Calibration Method of the Compton Edge Based
on Normalized Cross-correlation and Simulated Annealing Algorithm.* IEEE Transactions on Nuclear
Science, 2025. arXiv:2507.16282 [physics.ins-det].

73. **El informe se reestructuró para contar la lógica y la relación con el
   paper** (24-sep-2026). Pedido del usuario: que alguien entienda qué se hizo,
   por qué, y cómo se inspira en Kong et al. `informe/informe.tex` se armó de
   nuevo (respaldo: `informe_v1_2026-09-24.tex/.pdf`) con este orden:
   Introducción (problema, NCC-SA, qué hace este trabajo, resultados en un
   párrafo, cómo leer) → **El método de referencia y cómo lo adaptamos**
   (las 3 estrategias que el paper critica, NCC-SA en 7 pasos con sus Ec.,
   `flujo.tex` = diagrama TikZ lado a lado, `tab:adaptacion` = §7.71 como
   tabla, por qué la plantilla simulada no contradice al paper: su objeción
   es a la búsqueda en grilla, no a la simulación) → **Montaje** (BC-408 +
   XP2262/B, `tab:hardware`, cadena de señal en TikZ, por qué 20 dB, sesiones,
   lo no registrado) → Física (con subsección nueva de quenching y keVee) →
   Pulso → Geant4 → Ajuste (NCC/SA/barra condensados) → Resultados (calib,
   test B con k = b/12, test A, detector D) → Resultados frente al paper →
   **Conclusiones** (nuevo) → Uso → Limitaciones (sólo abiertas). A apéndices:
   NCC/SA en detalle, estadística y cierre, controles del observable, clásicos
   medidos, simulación, **auditoría con las 5 erratas** (antes el informe
   listaba 4 cosas), **historia de errores** (lo que era tono de bitácora en
   el cuerpo), glosario ampliado.
   - Etiquetas de procedencia `\enpaper`, `\adaptado`, `\propio`
     (preámbulo) al comienzo de cada sección de método: es la regla de §4
     ("marcar [fuera del paper]") hecha tipografía.
   - `\kMin`/`\kMax` se calculan en el `.tex` con pgfmath de `\BMin/12`,
     así que no se tipean.
   - Se corrigió: BC-404 → PVT/BC-408; "χ² de Pearson" en el párrafo del
     recocido cuando la campaña usa Neyman; "nueve puntos" y "25 corridas"
     tipeados y viejos; "4 sesiones" cuando eran 6 días; "geometría real"
     cuando sólo el tamaño es real.
   - ⚠️ La VM de Cowork no tiene `babel-spanish`: el `informe.pdf` que quedó
     es una vista previa compilada con nombres en castellano puestos a mano.
     El definitivo sale de `python analisis/informe.py` (MiKTeX).
   - Agregado después (24-sep-2026): `sec:dconj` (la D conjunta es la media de
     las D pesada por 1/σ², y por qué no se reporta como calibración),
     `sec:error-energia` (calibración vs resolución vs sistemáticos, con
     ejemplo: 1000 mV·ns en A@1600 Cs → σ_E ≈ 1.2 keV, calculado con pgfmath
     de `\EjErrEscalado` y `\DmilseisA`) y el apéndice `ap:propagacion`
     (fórmula por fórmula: σ_D escalado, 1/D, E, FWHM, D_conj y su χ², Cs/Co
     y desacuerdo, c0, b en log-log con pivote geométrico, interpolación en V,
     V_eq). Ojo: el err_b de `voltaje.ajustar` NO está escalado por
     √(χ²/ndf), y en varios grupos χ² ≫ ndf.
   - Pendiente en `informe.py`: la columna 1/D de `tabla_calibraciones` sale
     con un decimal (1.0, 0.8…), poca precisión para una columna que es la
     que se usa para convertir.

74. **La estadística del código, al informe, por bloques** (25-sep-2026).
   `informe/PLAN_estadistica.md` inventaría cada concepto estadístico que usa
   `analisis/` y dice si el informe lo explica. Bloques 0 y 1 hechos: el
   apéndice B quedó en subsecciones por etapa (`sec:est-conteo`, `-chi2`,
   `-incert`, `-sist`, `-cierre`) y B.2 explica la distribución de χ²_min
   (±√(2/ndf); la corrida de ejemplo está 10.3σ arriba con 192 gdl), los ndf
   con la amplitud contada, Baker-Cousins, nnls y los pulls. En `informe.py`:
   `estadistica_chi2`, `figura_pulls` (`teoria_pulls.pdf`), `macros_chi2`.
   Las macros nuevas tienen `\providecommand` de respaldo en la cabecera de
   `informe.tex` (salen "?" si `generado/` es viejo). Cowork no tiene
   `uproot`: la figura de pulls y sus números los genera el PC.
   ⚠️ Encontrado de paso: `pulso.py`/`area.py` tienen el método de área por
   **ajuste de plantilla** (`ajustar_plantilla`, `plantilla_pulso`,
   `metodo="plantilla"`) y un lote en `analisis/figuras/calibraciones_plantilla/`,
   sin entrada en este archivo. Si fue la sesión de Claude Code del 24-sep,
   conviene agregar acá qué midió.
   Bloque 2 (mismo dia): B.3 explica el hessiano por diferencias finitas y
   el paso por parametro (`modelo.PASO_CURV_REL`, nuevo nombre del 1e-3 que
   estaba suelto), la covarianza no definida positiva y el parche de
   `error_usable` (hoy no lo necesita ninguna corrida), Delta chi2 = 1 contra
   2.30 en dos dimensiones (elipse de puntos nueva en `teoria_estadistica`) y
   Wilks, con la condicion del borde C >= 0. `macros_incertidumbre`.
