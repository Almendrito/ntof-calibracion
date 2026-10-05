# `analisis/` — procesamiento de datos desde cero

Reimplementación limpia de la cadena **XML del DRS4 → espectro de energía**,
escrita el 27-ago-2026 sin importar nada de `codigo/`.

Va del XML crudo hasta la relacion energia-canal: espectro neto con barras de
error, seleccion automatica del intervalo (NCC), ajuste por recocido simulado, y
los dos tests de linealidad -el del eje y el de la ganancia contra el voltaje-.

```
python analisis/rutas.py             # qué corridas ve
python analisis/drs_xml.py           # E0 · leer el XML          -> figuras/00_*.png
python analisis/pulso.py             # E1 · la ventana           -> figuras/01_*.png
python analisis/area.py [N]          # E2 · los dos métodos      -> figuras/02_*.png
python analisis/espectro.py [N]      # E3 · el espectro          -> figuras/03_*.png

python analisis/hacer_espectro.py --corrida co60 --metodo ambos
python analisis/hacer_espectro.py --corrida cs137_1600 --limite 20000
python analisis/hacer_espectro.py --todas --limite 5000

python analisis/ver_espectros.py           # dibuja lo ya calculado, sin recalcular
python analisis/simulacion.py              # las plantillas Geant4 -> figuras/04_*.png
python analisis/ncc.py                     # el intervalo automatico -> figuras/05_*.png
python analisis/modelo.py                  # modelo y costo         -> figuras/06_*.png
python analisis/sa.py                      # el recocido            -> figuras/07_*.png
python analisis/calibracion.py             # LA CALIBRACION       -> figuras/08_*.png
python analisis/verificar_sa.py            # audita las decisiones de sa.py
python analisis/doc_sa.py                  # docs/recocido_simulado.pdf

python analisis/sistematicos.py            # el presupuesto de error, medido
python analisis/sistematicos.py --gates    # + el termino de definicion del area (relee los XML)
python analisis/voltaje.py                 # la ganancia contra el HV -> figuras/09_*.png
python analisis/gui.py                     # LA INTERFAZ
python analisis/pipeline.py co60           # la misma cadena, por consola
python analisis/registro.py                # las calibraciones acumuladas
python analisis/registro.py --rehacer      # recalcula TODO el catalogo y regraba
python analisis/registro.py --cierre       # el test de linealidad con verdad conocida
```

La primera pasada sobre un XML lo parsea (~40 s por GB) y cachea las áreas en
`analisis/.cache/`; las siguientes arrancan al instante. La clave del caché
incluye el archivo, el canal, el método y los parámetros de la ventana, así que
cambiar cualquiera de esos lo invalida solo. Los nombres son hashes:
`python analisis/area.py --cache` dice qué guardó cada uno.

---

## La cadena

| módulo | qué hace | la decisión que toma |
|---|---|---|
| `rutas.py` | catálogo de corridas y carpetas | nada; sólo dice dónde está cada cosa |
| `drs_xml.py` | XML → `(t, v)` en ns y mV | qué canal leer (lo **detecta**, no lo supone) |
| `pulso.py` | línea base, onset, integral acumulada | **cuál es la ventana de integración** |
| `area.py` | integrar todos los eventos | qué eventos se marcan y cuáles se usan |
| `espectro.py` | histograma, tiempo vivo, resta de fondo | cómo se normaliza la resta |
| `hacer_espectro.py` | los cuatro encadenados | — |
| `ver_espectros.py` | dibuja los `.npz` ya calculados | ninguna; no vuelve a tocar los XML |
| `fisica.py` | borde Compton y fuentes medidas | la Ec. 1 corregida (`2α` en el numerador) |
| `simulacion.py` | lee los `.root` de Geant4 | qué rama es la plantilla: `Lvis` (keVee) o `Edep` (keV) |
| `ncc.py` | selección automática del intervalo | **dónde se ajusta**, y la semilla de D |
| `modelo.py` | modelo directo `A·[P(ch/D) ⊗ G(σ)]` y χ² | qué se ajusta y con qué pesos |
| `sa.py` | recocido simulado, multiarranque, errores | cómo se busca el óptimo, y cómo se verifica |
| `calibracion.py` | **la relación energía-canal** | qué escala se reporta: keVee o keV |
| `verificar_sa.py` | contrasta las decisiones de `sa.py` | ninguna; audita |
| `doc_sa.py` | genera `docs/recocido_simulado.pdf` | ninguna; documenta |
| `pipeline.py` | **la receta única**, con progreso y cancelación | el orden de las etapas y cómo se pesa la barra |
| `gui.py` | la interfaz (9 pestañas; las 4 numeradas son la cadena) | ninguna; llama a `pipeline` |
| `registro.py` | calibraciones acumuladas por punto de operación | qué cuenta como el mismo punto, y el **test A** |
| `voltaje.py` | la ganancia contra el HV, `D(V)` | qué ley se ajusta, y si los datos **alcanzan** para decidir |
| `sistematicos.py` | el presupuesto de error de D | qué variante es incertidumbre y cuál está **descartada** |

⚠️ Los espectros parciales (los que salen con `--limite`, y los de las demos)
llevan `_lim<N>` o `demo_` en el nombre. Sin eso una corrida rápida de prueba
sobrescribe el resultado completo y no queda ninguna señal de que pasó — pasó
una vez con el Cs-137 y el `.npz` de 80 000 pulsos quedó reemplazado por uno de
8 000 sin ningún aviso.

### Cómo sale la ventana

Es el punto del proyecto y no hay ningún número tipeado:

1. Se leen ~600 pulsos, se les detecta el onset y se **alinean por él**
   (interpolando a una grilla común: el eje del DRS4 no es uniforme, así que no
   se puede alinear corriendo índices).
2. Se promedian normalizados por amplitud → la **forma** del pulso.
3. Sobre esa forma se calcula la **integral acumulada** `Q(t)` y se busca dónde
   empieza a crecer y dónde converge, con el ruido de la propia línea base como
   vara y exigiendo un tramo continuo de quietud (un solo punto plano no es
   convergencia: el pulso tiene un hombro).
4. Esos `(pre, post)` se aplican **anclados al onset de cada evento**.

Se hace sobre el promedio y no evento por evento por una razón medida: el ruido
del promedio de 400 pulsos es ~20 veces menor, así que ve la cola mucho más
abajo. Ver el hallazgo 4.

---

## Donde estan los datos

Las 16 corridas del catalogo estan copiadas dentro del proyecto, en
`datos/Detector X/<voltaje>/<muestra>` (12 GB, copiados del disco externo el
03-sep-2026). `rutas._ext` busca **primero ahi y despues en el disco**, asi que
todo corre con el disco desconectado; si algun archivo no esta copiado, sigue
usando el disco como antes.

Tres detalles de la copia que conviene no deshacer:

- la carpeta del detector A a 1700 V se llamaba `1700` en el disco y todas las
  demas `1700V`: en la copia esta normalizada, y `rutas._local` aplica la misma
  regla al buscar;
- `B1700V` era el unico archivo cuyo nombre no dice de que fuente es (es el
  Cs-137 del detector B a 1700 V): en la copia se llama `Cs137`, y
  `rutas._RENOMBRE_LOCAL` traduce;
- se copio con `shutil.copy2`, que **conserva la fecha de modificacion**,
  porque de eso depende el cache de areas. Ojo que hasta el 03-sep-2026 la
  clave de `area._clave` llevaba tambien la **ruta absoluta**, asi que la copia
  invalidaba el cache y obligaba a re-parsear 80 s por XML para dar el mismo
  numero (medido en `cs137_1600`: misma D = 1.9450 ± 0.0038). La ruta se saco
  de la clave -lo que identifica a estos datos es tamano + fecha, igual que en
  `codigo/datos_io/spectrum.py`-, asi que mover los datos de carpeta volvio a
  ser gratis. La cuenta de esa vez: rehacer los `.npz` de las 16 corridas.

Cada carpeta lleva un `LEEME.txt` escrito leyendo los archivos y no el nombre
(fuente, formato y canal detectado), porque los binarios del DRS4 no tienen
extension y en los dialogos hay que dejar el filtro en `*.*` para verlos.

Si conviene sacar los datos de OneDrive, la escotilla es la variable de entorno
`CALIBRACION_SESION`, que apunta la sesion a donde uno quiera.

## La interfaz

`python analisis/gui.py`. Panel de entradas a la izquierda (corrida del catálogo
**o archivos sueltos**, detector, fuente, escala keVee/keV, opciones del ajuste),
lista de las 7 etapas con su barra y su tiempo arriba, y pestañas de diagnóstico
abajo: Resultado, Ventana de integración, **1 · Datos y simulación**,
**2 · NCC e intervalo**, **3 · Semilla y recocido**, **4 · Inicial vs final**,
Θ y correlaciones, Paisaje del χ² y Registro.

**Las cuatro numeradas son la cadena en el orden en que corre**, un estado por
pestaña, y cada una tiene cuatro paneles:

1. **Datos y simulación.** Muestra y fondo con el binning común (log), el neto
   con la ventana de ajuste, la plantilla de Geant4 en su propio eje de energía
   y las dos superpuestas en energía con la D ajustada. Esa superposición es el
   control más directo de que la plantilla es de esta fuente: si los bordes no
   caen juntos, no hay ajuste que lo arregle (hallazgo 7.41 de CLAUDE.md, la
   falla silenciosa del `.root` de otra fuente). Se ve además el exceso de bajo
   canal como la diferencia entre el dato y la plantilla, y el escalón del
   umbral del trigger donde lo haya.
2. **NCC e intervalo.** La curva ρ(B) del barrido de un parámetro con una tabla
   de ρ a ±5% y ±10% de B (para leer si el máximo es agudo o una meseta), la
   ventana ganadora sobre el espectro —con la extensión hacia abajo que usa el
   ajuste marcada aparte—, las dos formas ya estandarizadas que Pearson compara
   de verdad, y el diagrama de dispersión de una contra la otra, que es la
   lectura geométrica de ρ. El ρ se **recalcula** en la pestaña y se muestra al
   lado del que devolvió el barrido: si no coinciden, la pestaña está
   reconstruyendo mal lo que el barrido hizo.
3. **Semilla y recocido.** De dónde arranca, cómo baja y dónde termina: la
   tabla de la semilla del NCC con las cotas y el óptimo, más una fila por
   arranque; el costo ciclo a ciclo de los N arranques con el ganador
   resaltado; temperatura y tasa de aceptación (si la aceptación arranca en ~1
   o en ~0, `T_init="auto"` calibró mal); y el camino en el plano (D, σ_rel),
   que es lo que distingue "bajó derecho" de "cruzó de un valle a otro". Ojo
   con lo que la tabla dice explícitamente: los arranques salen **al azar
   dentro de las cotas**, no de la semilla (`sa.multiarranque` con `x0=None`),
   y por eso que todos terminen en el mismo punto es evidencia de algo.
4. **Inicial vs final.** El modelo con la semilla del NCC y con el vector
   ajustado sobre los mismos datos, los residuos normalizados de los dos, la
   tabla de cuánto se movió cada parámetro (en % y en unidades de su propio
   error: si se movió menos que su barra, el ajuste no lo cambió de verdad) y
   un corte del χ² sobre el **segmento** que une los dos vectores. Que el
   mínimo de ese corte caiga en t = 1 es la comprobación más barata de que el
   ajuste no se quedó a medio camino ni se pasó de largo.

Lo que hace posible la pestaña 3 y la 4 es que `calibracion.calibrar` acepta un
dict `detalles` y lo llena con la semilla completa, las cotas y las N corridas
del recocido con su historia; `pipeline` lo deja en `salida["sa"]`. No entra en
el resultado porque el resultado se serializa a JSON y esto son arrays.
`sa.recocido` guarda ahora también el **camino** (`historia_x`), que cuesta una
fila de floats por ciclo de temperatura.

**El espacio de los gráficos, y las cuatro cosas que lo decidían.** Las cuatro
pestañas nuevas son grillas de 2×2, o sea que lo que se ve depende del alto que
tenga cada panel. Medido en esta pantalla (1920×991 con la ventana maximizada):

| | lienzo | fila de la grilla |
|---|---|---|
| lista de etapas cruzada arriba, `tight_layout` | 11.95 × 5.53 in | 2.76 in |
| lista en la columna izquierda, `constrained` | 11.19 × 7.20 in | 3.60 in |

- **La ventana arranca maximizada** (`state("zoomed")`). En 1500×900 cada panel
  queda de ~430×250 px, que es justo donde los rótulos empiezan a estorbar.
- **La lista de etapas vive en la columna izquierda**, no cruzada arriba del
  gráfico. Se lleva ~170 px de **alto** y sólo necesita ~430 de ancho, y la
  columna izquierda tiene ese ancho y le sobra alto: moverla le devuelve el 30%
  de altura al gráfico sin sacárselo a nada. El interruptor «ver etapas» de la
  barra de abajo la pliega, para ventanas chicas donde la columna no alcanza
  para las entradas y la lista a la vez. No se pliega sola al terminar: los
  tiempos por etapa son parte del diagnóstico, y una ventana que se reacomoda
  sin que nadie la toque se lee como un error.
- **`gui.ajustar(fig)` en vez de `fig.tight_layout()`**, y es la primera línea
  de cada figura, no la última. `tight_layout` resuelve el espaciado **una vez**,
  con el tamaño que la figura tiene en ese momento, y el lienzo de Tk se estira
  con la ventana sin volver a llamarlo: en una grilla de 2×2 los títulos y los
  rótulos se terminan pisando en cuanto la ventana no mide lo que medía al
  calcular. El motor `constrained` lo recalcula en cada *draw*. Va primero
  porque `fig.colorbar` mira si hay motor de layout: si no lo hay, mete la barra
  en un gridspec anidado con `height_ratios = [0, 1, 0]`, y encender el motor
  después hace que el layout divida por ese cero (`ZeroDivisionError` al
  dibujar, lejos de la causa).
- **Tipografía a nivel de módulo** (9 pt de cuerpo, 8 pt en los ejes) para que
  las figuras salgan iguales cuando las dibuja `_verif_gui_pestanas.py`, que no
  crea ventana.

Tres cosas que costaron encontrar y conviene no re-romper:

- **`ax.table(bbox=...)` es (x0, y0, ANCHO, ALTO)**, no (x0, y0, x1, y1). Pasarle
  `(0, 0.09, 1, 1)` para «dejar 9% abajo» da una tabla de alto 1 que arranca en
  0.09: se sale por arriba de los ejes y tapa el título.
- **Un `ax.text` largo cuenta como decoración del eje.** Constrained layout mide
  el *tight bbox* con el texto adentro y, si el texto es más ancho que el eje,
  **achica el eje** hasta que entre: con una nota de 100 caracteres en una línea
  la tabla quedaba de un quinto del panel. Por eso `tabla` parte la nota en
  líneas de 46 y le reserva 0.048 de alto por línea.
- **Un título largo con `loc="left"` se sale del ancho de sus ejes** y se le monta
  al título del panel de al lado; el motor no lo acorta ni lo parte. Los títulos
  de dos líneas están escritos así a mano, y el título de las tablas se dibuja
  adentro de los ejes.
- **matplotlib no achica el texto de una celda: lo RECORTA**, sin avisar (así se
  veían «punto INICIAI» y «A (VARPR»). Lo primero es contenido corto: unidades
  pegadas al nombre, cuatro columnas en vez de seis. Acá decía además que medir
  el ancho de la figura no servía «porque se dibuja una vez y después el lienzo
  la reescala»; eso valía mientras la figura se construía con el tamaño de
  fábrica, y dejó de valer con `Lienzo._encuadrar` (abajo), que le pone el
  tamaño real **antes** de llamar a la función de dibujo. `tabla` saca ahora el
  cuerpo del ancho del panel, entre 5.8 y 8 pt. `auto_set_font_size(True)` no
  sirve para esto: sólo sabe achicar, y una vez que encogió no vuelve a crecer
  aunque la ventana se agrande.

### Por qué los gráficos salían cortados (08-sep-2026)

El síntoma era «en las pestañas no caben todos los gráficos y se ven cortados».
No era el espaciado: **faltaba un pedazo del dibujo**. TkAgg escala la figura
por el *device pixel ratio* de la pantalla —en Windows, `tk scaling / (96/72)`,
que con el `tk scaling = dpi/72` que declara la App vale 1.25— pero lo aplica al
recibir `<Map>`, o sea **después** de que `FigureCanvasTk.resize` fijó el tamaño
en pulgadas con el dpi viejo. Lo único que hace para reconciliarlos es pedirle al
widget que crezca, y adentro de un Notebook dentro de un Panedwindow ese pedido
no se cumple: no llega ningún `<Configure>` nuevo, el buffer de Agg queda más
grande que la PhotoImage sobre la que se vuelca y el blit **recorta en
silencio**. Medido con la ventana en 1550×830:

| | widget | figura | buffer de Agg |
|---|---|---|---|
| antes | 1029×739 px | 8.23 × 5.91 in a 156 dpi | 1286×923 px ← se ve el 64% |
| ahora | 1008×697 px | 8.06 × 5.58 in a 125 dpi | 1008×697 px |

Un recorte así **no se ve en un PNG guardado con `savefig`**, que es lo que
prueba `_verif_gui_pestanas.py`: sólo se ve en la ventana. De ahí
`verificaciones/_verif_gui_encuadre.py`, que compara los dos tamaños —el del
render y el de la imagen— en cada pestaña y a tres tamaños de ventana.

Cuatro cosas salieron de eso, y las cuatro son la misma causa:

- **El lienzo PIDE un tamaño chico y fijo** (`Lienzo.PEDIDO`, 360×240) y se
  conforma con el que le den, y `Lienzo._encuadrar` le da a la figura el tamaño
  exacto del widget con el dpi que la figura tenga en ese momento. Sin lo
  primero, el Panedwindow terminaba pidiendo 2490×1194 px —más que la pantalla—
  y en `pack` los hijos que se empaquetan al final se quedan sin lugar: con la
  ventana sin maximizar **desaparecían la barra de progreso de abajo y la barra
  de herramientas de matplotlib**. La barra de abajo se empaqueta ahora primero.
- **El dpi de la figura va en 100**, el nominal. Estaba en `100·dpi/96` = 125
  para compensar el escalado de Windows, que es justo lo que TkAgg ya hace: el
  factor entraba dos veces, la figura terminaba a 156 dpi y todo se dibujaba un
  25% más grande de lo pedido, o sea con un **36% menos de superficie útil**.
- **La columna de controles se desplaza** cuando no entra. El Panedwindow
  reparte el ancho, pero el alto lo impone la ventana: los controles piden 907
  px y la ventana sin maximizar tiene 830, así que lo que sobraba quedaba
  cortado abajo y no había forma de llegar a él. La barra de desplazamiento
  aparece sólo cuando hace falta, y la decisión se toma **en reposo**
  (`after_idle`) y no dentro del `<Configure>`, porque ahí `winfo_height`
  todavía puede ser el alto de antes.
- **Doble clic sobre un panel lo agranda** a la figura entera; doble clic de
  nuevo, o el selector que está al lado de la barra de herramientas, vuelve a la
  grilla. Aun sin recortes, cuatro paneles en 1008×697 px son ~500×350 cada uno.
  No se recalcula nada: se guardan posiciones, visibilidad y anclas, se apaga el
  motor de layout y se reubica el panel elegido. Eso es lo único que lo hace
  usable en el paisaje del χ², donde volver a llamar a la función de dibujo
  costaría ~6300 evaluaciones del modelo. Los `twinx` se detectan por posición
  —nacen con el mismo rectángulo que su panel— y las barras de color se anotan
  mientras se dibuja, envolviendo `fig.colorbar`: después no hay cómo saber de
  quién es cada una (su eje sólo se distingue por la etiqueta `<colorbar>`).

**Cuando la ventana se cierra sola.** Un cierre por violación de acceso no deja
excepción: no hay `try/except` que lo vea ni mensaje que mostrar. Hay tres
archivos en `.cache/` para eso, y conviene mirarlos en este orden:

1. **`gui_traza.log`** — una línea por etapa, por cambio de pestaña y por figura
   dibujada (antes y después). **La última línea dice hasta dónde llegó.** Es lo
   que distingue «se cayó calculando» de «se cayó dibujando», y no es una
   distinción académica: el cierre del 03-sep-2026 se reportó como «crashea al
   ejecutar el SA» y la traza mostraba `fin sa`, `fin calib` y
   `--- terminado, D = 2.3984`. El recocido había terminado; lo que seguía era
   el dibujado.
2. **`gui_fallo_anterior.log`** — el volcado de `faulthandler` de la sesión que
   se cerró. Si el hilo principal aparece con **sólo `mainloop`** en la pila, el
   proceso murió dentro del C de Tk procesando eventos y no dentro de un
   callback nuestro. La ventana **avisa al abrirse** si la sesión anterior
   terminó así, y pega la pila en el panel de Resultado.
   ⚠ Antes esto no servía de nada: `faulthandler` abría `gui_fallo.log` en modo
   `"w"`, así que relanzar el programa —lo primero que uno hace— borraba la
   única evidencia. Ahora se rota antes de truncar.
3. **`gui_error.log`** — las excepciones de Python: las de los callbacks de Tk
   (que Tkinter se traga y manda a stderr, o sea a ninguna parte si no hay
   consola) y las de las figuras, que en pantalla son el cartel gris «no se
   pudo dibujar».

**Por qué las figuras se dibujan al abrir cada pestaña.** Una pestaña que Tk
todavía no mapeó mide **1×1 px**. Pintar ahí no le muestra nada a nadie y no es
gratis: TkAgg vuelca el render en una PhotoImage que después dibuja Tk por su
cuenta. `_terminado` pintaba las siete de una, seis de ellas en lienzos de 1×1.
Ahora `Lienzo.programar` deja el dibujo pendiente y `asegurar` lo hace cuando la
pestaña tiene tamaño de verdad (`MINIMO = 40 px`, más el `<Configure>` que llega
al mostrarla). Terminar una corrida cuesta un dibujo en vez de siete.

Que esto haya sido *la* causa de aquel cierre no está probado —60 vueltas
pintando en lienzos de 1×1 no lo reprodujeron, ni tres calibraciones seguidas
sobre B, C y D a 1700 V—, pero saca de la ventana trabajo frágil que no le
servía a nadie, y la próxima vez la traza dirá qué figura estaba en curso.

Las figuras se verifican sin abrir la ventana con
`python verificaciones/_verif_gui_pestanas.py`, que corre la cadena y llama a
cada método de figura directamente: dentro de `Lienzo.dibujar` una excepción se
convierte en un cartel gris en el propio gráfico, así que un error de
programación no se ve hasta que alguien abre esa pestaña. Así estuvo roto el
paisaje del χ²: `Problema.mapa` no sabía en qué valor congelar `lam` y tiraba
`TypeError`. Ahora `mapa` pide `x_base` y el corte se dibuja con el resto de los
parámetros en el óptimo.

**Θ y correlaciones** (antes «Parámetros del ajuste») es el vector Θ completo, no solo D: una tabla con
cada parámetro libre (D, σ_rel y, según el modelo, λ_exceso y/o c0) con su
valor, su error de curvatura ya escalado por √(χ²/ndf), y el error relativo;
las amplitudes A y C_exceso se listan aparte porque se proyectan exactas
(VARPRO/nnls) y no tienen curvatura. Al lado, la matriz de **correlación**
entre los parámetros libres: un σ_rel con error grande casi siempre viene con
una correlación fuerte con D (ver `corr_D_sigma` en `calibracion.calibrar`),
y la tabla sola no lo muestra. Abajo de la tabla van los tres indicadores de
**fiabilidad** que no son "un parámetro" pero sí dicen si el resultado es de
fiar: χ²/ndf (¿el modelo describe los datos?), la dispersión entre los N
arranques del SA (¿el óptimo es reproducible o cada arranque cae en un valle
distinto?) y ρ del NCC (¿la ventana automática correlaciona bien con la
plantilla?). Los datos salen de `r["libres"]/["x_libres"]/["err_libres"]/
["cov_libres"]`, que `calibrar()` guarda además de los campos con nombre
propio (`D`, `sigma_rel`, ...) para no romper lo que ya los lee.

**Por qué la barra está hecha así.** El parseo del XML es ~40 s y todo el resto
de la cadena junto ~3 s: el 95% del tiempo transcurre en un paso que, sin barra,
no dice nada. De ahí tres decisiones:

- la fracción del parseo se mide en **bytes leídos**, no en eventos (cuántos
  eventos tiene el archivo no se sabe hasta terminarlo);
- los pesos de la barra global se calculan **antes de empezar**, mirando el
  tamaño de cada XML y si sus áreas ya están cacheadas (`pipeline.repartir`).
  Los dos XML se reparten por bytes, no mitad y mitad: la muestra pesa 5 veces
  el fondo, y con el reparto ingenuo el ETA sale el doble;
- hay **una fila por etapa con su tiempo**. Una barra sola dice "va por el 30%";
  la lista dice además cuál etapa corre, cuánto tardó cada una y —si algo
  falla— dónde falló.

El cálculo va en un hilo aparte y se comunica por una cola; el hilo trabajador
nunca toca widgets ni matplotlib. Se puede cancelar a mitad del parseo.

**La receta vive en `pipeline.py`, no en la interfaz.** La GUI y la consola
llaman a la misma función. Si la GUI tuviera su copia, las dos divergirían: es
exactamente lo que pasó en el proyecto anterior.

**El registro es el pago de la campaña.** Cada calibración se guarda con su
detector. Cuando haya **dos fuentes distintas medidas con el mismo detector**,
la pestaña Registro ajusta la recta por el origen y con ordenada libre, y
reporta si la ordenada es compatible con cero. Ese es el test de linealidad, que
hoy no se puede hacer: los dos bordes del Co-60 están atados a la misma D por el
mapa por el origen (Ec. 12 del paper).

---

## Lo que apareció al construirlo

1. **Las dos corridas tienen ejes de tiempo distintos.** El Co-60 trae la
   calibración por celda del DRS4 (`dt` de 0.10 a 0.32 ns, 154% de variación);
   el Cs-137 del 26-ago trae un eje uniforme redondeado a dos valores
   (0.1950/0.1960), o sea se guardó **sin** esa corrección. Por eso toda
   integral va con `x=t` y nunca con un `dt` constante.

2. **El registro no empieza en cero.** `t[0]` vale −0.01 ns en el Co-60 y
   **27.6 ns** en el Cs-137. Todo se mide contra `tau = t - t[0]`. Además la
   primera y la última celda del anillo dan espigas de decenas de mV que no son
   señal: se recortan 2 muestras de cada extremo.

3. **La integral no converge a un plateau plano: hace techo y después baja.** En
   el Co-60 el techo está a 70 ns del onset y de ahí cae **1.66%**: hay un
   rebote de polaridad opuesta en la cola. Por eso la curva de convergencia se
   normaliza contra el techo y no contra el final del registro (si no pasa por
   101.7% y "converger al 100%" deja de significar algo).

4. **El método "convergencia evento por evento" corta con un umbral que depende
   de la amplitud, y eso deforma el espectro.** El criterio corta donde la señal
   se hunde en el ruido *de ese evento*, y en un pulso chico eso pasa antes.
   Ancho de la ventana por cuartil de amplitud, Co-60:

   | amplitud [mV] | 4–15 | 15–53 | 53–107 | 107–287 |
   |---|---|---|---|---|
   | ancho mediano [ns] | 13.1 | 26.7 | 35.1 | **42.1** |

   Consecuencia sobre la fracción de carga recogida (contra una ventana ancha de
   130 ns, bineado por amplitud):

   | corrida | global | convergencia |
   |---|---|---|
   | Co-60 | recoge 100.9%, deriva −16.6%* | 99.1%, deriva −9.9%* |
   | Cs-137 | recoge 97.1%, **deriva +4.1%** | 84.1%, **deriva +24.6%** |

   \* las derivas del Co-60 están dominadas por un pedestal de la ventana de
   referencia, no por el método; ver el hallazgo 5.

   **Veredicto: se usa el método global.** En el Cs-137 el método por
   convergencia mueve el borde ~15% hacia abajo y le mete un pico espurio a bajo
   canal; se ve directo en `figuras/espectro_cs137_1600.png`.

5. **Un cociente que se aparta de 1 no dice por qué; la diferencia sí.**
   Ajustando `area(ancha) − area(ventana) = c₀ + c₁·amplitud`:

   | corrida | c₀ (pedestal) | equivale a | término proporcional |
   |---|---|---|---|
   | Co-60 | −13.3 mV·ns | −0.19 mV sobre los 70 ns extra | +1.6% (el rebote) |
   | Cs-137 | +16.3 mV·ns | +0.51 mV sobre los 32 ns extra | +0.3% |

   Un término **constante** es un error de línea base multiplicado por los ns de
   más que integra la ventana ancha: el problema está en la *referencia*, no en
   la ventana que se está juzgando. Un término **proporcional** sí es forma. Sin
   esta separación, el Co-60 parecía tener una deriva de −16.6% que en realidad
   es un pedestal de la vara.

6. **El test que de verdad decide es la forma del pulso por cuartil de
   amplitud.** Fracción de la carga dentro de la ventana:

   | cuartil | Co-60 | Cs-137 |
   |---|---|---|
   | 1 (más chicos) | 96.86%* | 97.65% |
   | 2 | 99.95% | 96.55% |
   | 3 | 99.84% | 97.39% |
   | 4 (más grandes) | 99.53% | 98.02% |

   \* el primer cuartil del Co-60 está por debajo por *alineación*, no por
   forma: con amplitudes de 5–15 mV el onset es ruidoso y el promedio sale
   emborronado. Arriba de 15 mV la forma es constante al **0.4%**, o sea la
   ventana global recoge una fracción constante y el área es proporcional a la
   energía. Eso es lo único que la calibración necesita.

7. **El Cs-137 tiene una réplica del pulso a +35 ns del onset** (se ve en los
   pulsos crudos cerca de t = 105 ns y en la forma promedio), más un rizado de
   ~4 ns de período en la cola. Tiene toda la pinta de reflexión por
   desadaptación de impedancia. Por eso su integral no converge dentro del
   registro: la ventana sale de 100 ns y aun así deja **2.7%** de la carga
   afuera. No invalida la calibración —la fracción perdida es la misma para
   todos los pulsos (hallazgo 6)— pero conviene arreglarlo en el montaje.

8. **El tiempo vivo confirmado, independientemente.** Integrando el `<ScalerN>`:

   | corrida | pared | vivo | fracción | scaler |
   |---|---|---|---|---|
   | Co-60 | 699 s | **283 s** | 40.4% | 170 Hz |
   | fondo Co-60 | 986 s | 986 s | 100% | 20 Hz |
   | Cs-137 1600V | 1373 s | **2 s** | 0.1% | 46 350 Hz |
   | fondo Cs-137 | 2669 s | 2669 s | 100% | 10 Hz |

   El factor de resta del fondo del Co-60 es `k = 0.2866`; con tiempo de pared
   habría sido 0.7091, o sea se habría sobre-restado el fondo 2.5 veces.

9. **Elegir la ventana no era una decisión neutral: era elegir cuánto sesgo
   aceptar.** Con el modelo desnudo, D deriva de forma MONÓTONA con el arranque
   de la ventana, y el χ² mejora todo el tiempo:

   | `inicio_rel` | 0.45 | 0.55 | 0.65 | 0.75 | 0.80 |
   |---|---|---|---|---|---|
   | D (Co-60) | 1.2230 | 1.2302 | 1.2391 | 1.2422 | 1.2430 |
   | χ²/ndf | 1.576 | 1.372 | 1.165 | 1.090 | 1.067 |

   No es dispersión, es una tendencia: abajo hay cuentas que la simulación no
   tiene y el ajuste las absorbe ensanchando σ y bajando D. La salida no es
   esquivar la zona sino **modelarla**: `S = A·[P(ch/D) ⊗ G(σ)] + C·exp(−(ch−lo)/λ)`,
   con λ como tercer parámetro no lineal y (A, C) proyectadas por nnls. Con eso
   se ajusta desde 0.25·B y D deja de depender de la ventana (1.2512 con 0.25 y
   1.2524 con 0.35).

   Cierra además el desacuerdo con la implementación anterior: `analisis/` daba
   1.2230 contra 1.2618 de `codigo/`, un 3.1% que era ocho veces el error
   estadístico. La mitad era la ventana y la otra mitad la mezcla de scattering
   `w`, que `codigo/` ajusta y acá la fija la simulación: con `w` fijada en el
   0.751 de la simulación, `codigo/` da 1.2537 contra nuestro **1.2512, 0.2%**.

10. **El umbral del trigger es un piso duro, y ninguna ventana puede entrar
    ahí.** El Cs-137 a 1600 V no tiene ni una cuenta debajo de 121 mV·ns y sube
    hasta ~334: eso es la curva de eficiencia del discriminador, no física. Ni
    la plantilla ni la componente de exceso la contienen, así que ajustar ahí
    adentro no sesga, **rompe**: arrancando en 229 el χ²/ndf salta a 4.26 y la
    matriz de errores sale singular; arrancando en 320 da 0.95.
    `espectro.fin_del_umbral` lo detecta como el primer máximo local del
    espectro suavizado (la eficiencia sube, todo lo que hay debajo del borde
    baja, y donde deja de subir una cosa le pasa el mando a la otra) y
    `modelo.armar` lo impone como piso. En el Co-60, tomado con el umbral en
    0.2 mV, cae en 56 mV·ns y no restringe nada, que es lo correcto.

11. **La barra de error de un parámetro puede estar mal por elegir mal el paso
    de derivación.** El hessiano usaba `paso = 1e-3·|x|` para todos. Para la
    ordenada c0 del mapa eso da 0.16 mV·ns, un 3% de un bin de 5.4, y como la
    plantilla se lleva a canal con `np.interp` -o sea que es lineal a trozos
    sobre la grilla del histograma- lo que se mide ahí es el quiebre entre dos
    trozos, no la curvatura del χ². Daba **21σ donde el cociente de
    verosimilitud da 2.1**. Ahora el paso se pide por parámetro
    (`Problema.pasos_curvatura`) y para c0 vale un bin. Control: para D, cuya
    escala sí es su propio valor, el paso no cambia nada (0.0040 con 1e-2,
    0.0044 con 3e-4) y coincide con el perfil real de Δχ² = 1, que da 0.0040.

12. **El test del origen con una sola fuente da una cota, y de paso mide su
    propia debilidad.** Liberando la ordenada (`canal = D·E + c0`) sobre una
    sola corrida: Co-60 da c0 = −163 ± 35 mV·ns con Δχ² = 4.3 (2.1σ), Cs-137 da
    +60 ± 34 con Δχ² = 0.5 (0.7σ). **Signos opuestos**, o sea no hay un pedestal
    consistente. Y liberar c0 mueve D un 13% en el Co-60 y un 7% en el Cs-137
    por una mejora de χ² casi nula: eso es la degeneración entre c0, σ y el
    exceso de bajo canal, cuantificada. El test limpio sigue necesitando dos
    fuentes.

---

## Control cruzado contra la implementación anterior

Nada de `analisis/` importa nada de `codigo/`, así que la coincidencia es una
verificación real y no una tautología:

| cantidad | `analisis/` (nuevo) | `codigo/` (anterior) |
|---|---|---|
| tiempo vivo del Co-60 | 283 s de 699 (40.4%) | 283 s de 699 (40%) |
| `k_fondo` del Co-60 | 0.2866 | 0.287 |
| tiempo vivo del Cs-137 1600V | 2 s de 1373 (0.13%) | 1.7 s de 1373 (0.13%) |
| ventana de integración | medida: `[−2.2, +60.3]` ns, 62.5 ns de ancho | fijada a mano: `[−8, +55]`, 63 ns |
| borde del Co-60, 50% de altura | 1346 mV·ns | 1229 × 1.113 = 1368 mV·ns |
| D del Co-60 [mV·ns/keVee] | **1.2512** | 1.2537 con `w` fija en la simulación (1.2618 con `w` libre) |

La ventana la anterior la eligió a ojo y esta la mide; dan el mismo ancho por
caminos distintos. El borde coincide al 1.6%.

---

## Las dos linealidades, y qué falta para cerrarlas

Se confunden todo el tiempo y son cosas distintas:

| | qué valida | qué necesita | dónde |
|---|---|---|---|
| **Test A** — el eje, `canal = D·E` | la calibración | varias **fuentes** al mismo punto de operación | `registro.linealidad` |
| **Test B** — la ganancia, `D(V)` | el fotodetector | la misma **fuente** a varios HV | `voltaje.ajustar` |

Hoy **ninguno de los dos se puede hacer**: hay una fuente por punto de operación
(Co-60 a 1550 V, Cs-137 a 1600 V) y un voltaje por fuente. Lo que falta es
medición, y está en `docs/campana_linealidad.md` en orden de prioridad.

Dos cosas que conviene tener presentes antes de mirar cualquier resultado:

- **El test B no espera una recta.** Un PMT tiene `D ∝ V^b` con b entre 6 y 10.
  La pregunta es cuánto vale b y si los datos alcanzan para rechazar la recta;
  con los tres voltajes que hay (1600-1700, span 6%) la separación entre las dos
  hipótesis es del 2.6%, o sea del tamaño del sistemático: **no decide**, y el
  informe lo dice en vez de dibujar una recta.
- **El test A va en keVee, no en keV.** El quenching de Birks vale 0.9541 en el
  borde del Co-60 y 0.9349 en el del Cs-137 (`calibracion.quenching_simulacion`):
  un 2.0% de diferencia. En escala Edep dos fuentes tienen que dar D distintas
  aunque la electrónica sea perfectamente lineal; en Lvis no, porque la plantilla
  ya trae el quenching adentro.
