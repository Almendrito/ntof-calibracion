# El formato binario del DRS4 (los archivos sin extensión)

> Para sesiones futuras: **los datos nuevos no son XML.** Este documento dice qué
> son, cómo se leen, y qué trampas tienen. Verificado byte a byte contra los
> archivos de la sesión del 01-sep-2026.

---

## 1. Qué archivos son

Están en el disco externo, dentro de `D:\Sesión Calibración 26_08_2026\`. Ojo con
el nombre de esa carpeta: **adentro conviven dos sesiones**, y lo dicen las marcas
de tiempo que el propio DRS4 escribe en cada evento.

| sesión | qué | formato |
|---|---|---|
| **26-ago-2026** | Cs-137 del detector A a 1600/1650/1700 V, más el único fondo que existe (`BG.xml`, sólo a 1600 V) | XML |
| **01-sep-2026** | los cuatro detectores A, B, C, D a 1600/1650/1700 V. Cs-137 en los cuatro, **Co-60 sólo en el A**. Sin fondo | **binario** |

Los binarios **no tienen extensión**: se llaman `Co60`, `Cs137`, `B1700V`. Por eso
el despachador decide por los **primeros 4 bytes**, no por el nombre. Todos pesan
exactamente 167 044 112 bytes = cabecera de 4112 + 80 000 eventos de 2088.

La bitácora de la sesión está en `Bitácora.txt`, en esa misma carpeta.

---

## 2. Cómo leerlos en este proyecto

No hace falta tocar nada: `analisis/drs.py` elige el lector solo.

```python
from analisis import drs, rutas

c = rutas.corrida("a_co60_1600")          # está en el catálogo
print(drs.formato(c.muestra))             # 'binario'
print(drs.detectar_canal(c.muestra))      # ('CHN1', 'Scaler0')
drs.imprimir_metadatos(drs.metadatos(c.muestra))

for t, v, reloj, scaler in drs.iter_eventos(c.muestra, limite=10):
    ...            # t [ns], v [mV], reloj 'YYYY/MM/DD HH:MM:SS.mmm', scaler [Hz]
```

`drs_bin` expone **exactamente** la misma interfaz que `drs_xml`
(`detectar_canal`, `metadatos`, `iter_eventos`, `leer_pulsos`), así que `area`,
`pulso`, `espectro`, `pipeline` y la GUI no se enteran del formato.

* `analisis/drs.py` — despachador. Sniffea `DRS2`/`DRS8` en los primeros 4 bytes.
* `analisis/drs_bin.py` — el lector. `python analisis/drs_bin.py` imprime los
  metadatos de los binarios del catálogo y dibuja unos pulsos.
* `python analisis/drs.py` — lista qué formato tiene cada corrida del catálogo.
* `python analisis/revisar.py` — revisión completa de una tanda antes de calibrar.

---

## 3. El formato, byte a byte

Es el `Save as binary` del DRSOSC, el mismo que lee el `read_binary.cpp` que
viene con el software del DRS4.

### Cabecera (una vez, al principio)

| bytes | contenido |
|---|---|
| 4 | `"DRS2"` (versión del archivo; existe `"DRS8"`, no visto acá) |
| 4 | `"TIME"` |
| | **por cada placa:** |
| 2 | `"B#"` |
| 2 | número de serie, `uint16` little-endian (acá **2859**) |
| | **por cada canal de esa placa:** |
| 4 | `"C001"`, `"C002"`, … |
| 4096 | ancho de cada celda, 1024 × `float32` [ns] |

En nuestros archivos: 1 placa, 1 canal → cabecera de **4112 bytes**.

### Evento (se repite hasta el final del archivo)

| bytes | contenido |
|---|---|
| 4 | `"EHDR"` |
| 4 | número de serie del evento, `int32` |
| 16 | año, mes, día, hora, minuto, segundo, milisegundo, **rango**, 8 × `uint16` |
| | **por cada placa:** |
| 2 + 2 | `"B#"` + serie |
| 2 + 2 | `"T#"` + **celda de trigger**, `uint16` |
| | **por cada canal:** |
| 4 | `"C00n"` |
| 4 | **scaler** `int32` [Hz] — *puede no estar*, ver abajo |
| 2048 | muestras, 1024 × `uint16` |

En nuestros archivos: **2088 bytes por evento**, 80 000 eventos.

### Conversión a unidades físicas

```
v[mV] = (u16 / 65536 + rango/1000 - 0.5) * 1000
t[i]  = suma de ancho[(j + celda_de_trigger) % 1024]  para  j < i
```

* `rango` sale de la cabecera del evento (es el centro del rango en mV; **0** en
  nuestros archivos, o sea −500 a +500 mV). Se declara `uint16` pero representa
  un valor con signo: si es > 32767 hay que restarle 65536.
* El eje temporal **no es uniforme y cambia evento a evento** con la celda de
  trigger — igual que en el XML. No se puede cachear uno solo.
* Un par de anchos de celda salen levemente **negativos** (−0.002 ns): es la
  calibración del DRS4, no un error de lectura. La suma de los 1024 da 200.00 ns
  exactos → 5.12 GS/s.
* **La polaridad no se toca.** Los pulsos de estos archivos son negativos
  (amplitud mediana −84 mV), pero `pulso.py` deduce el signo de cada pulso.

---

## 4. Dos cosas que NO se dan por supuestas

`drs_bin.estructura()` las deduce del archivo en vez de cablearlas, y levanta
`ValueError` si no cierran, en vez de devolver basura:

1. **Si el scaler está o no.** Las versiones viejas del DRSOSC no lo escriben. Se
   prueban los dos tamaños de evento, se toma el que divide **exacto** el tamaño
   del archivo, y después se verifica que `"EHDR"` caiga donde ese tamaño predice
   en ocho eventos repartidos por todo el archivo.
2. **Cuántas placas y cuántos canales.** Se cuentan en la cabecera.

Es la misma lección que §7.36 de `CLAUDE.md` (el `CHN2` cableado que devolvía
**cero pulsos en silencio** sobre los XML del 26-ago): *lo que no se verifica, se
cachea*.

---

## 5. Lector mínimo, sin el proyecto

Si hace falta leer uno de estos archivos desde afuera:

```python
import numpy as np, struct

def leer_drs4_binario(ruta, n_max=None):
    """Devuelve (t[ns], v[mV], celda, scaler) por evento. Un canal, una placa."""
    with open(ruta, "rb") as f:
        cab = f.read(4112)
        assert cab[0:4] in (b"DRS2", b"DRS8") and cab[4:8] == b"TIME"
        serie = struct.unpack_from("<H", cab, 10)[0]
        ancho = np.frombuffer(cab, dtype="<f4", count=1024, offset=16)
        k = 0
        while True:
            ev = f.read(2088)
            if len(ev) < 2088 or ev[0:4] != b"EHDR":
                return
            rango = struct.unpack_from("<H", ev, 22)[0]
            rango = rango - 65536 if rango > 32767 else rango
            celda = struct.unpack_from("<H", ev, 30)[0]
            scaler = struct.unpack_from("<i", ev, 36)[0]
            u = np.frombuffer(ev, dtype="<u2", count=1024, offset=40)
            v = (u / 65536.0 + rango / 1000.0 - 0.5) * 1000.0
            t = np.concatenate(([0.0], np.cumsum(np.roll(ancho, -celda))[:-1]))
            yield t, v, celda, scaler
            k += 1
            if n_max and k >= n_max:
                return
```

Los 4112 y 2088 valen para **1 placa y 1 canal con scaler**. Con otra
configuración hay que recalcularlos con las tablas de §3 (o usar
`drs_bin.estructura`, que lo hace solo).

---

## 6. Trampas propias de estas corridas

### 6.1 La compuerta de integración no cabe en el registro

**Es el problema serio de la sesión del 01-sep.** El trigger quedó en
**t ≈ 123 ns de un registro de 200 ns**, así que después del onset sobran ~81 ns
y el pulso necesita **83 ns** para converger al 99.9% de su carga. Con la ventana
que pide el pulso, el **84% de los eventos sale truncado**.

Un evento truncado pierde una fracción de carga que depende de *cuándo llegó el
pulso*, que es exactamente el error que la compuerta anclada al onset existe para
evitar (§7.33).

`area.calibrar_ventana` lo resuelve solo: mide el espacio disponible tras el
onset **en los eventos de amplitud sobre la mediana** (el walk del trigger hace
que los pulsos chicos tengan el onset tarde y sesguen el percentil crudo),
descuenta 5 ns de margen y recorta `post`. Devuelve `post_pedido`, `recortado`,
`fraccion` (carga recogida) y `frac_eventos_enteros`. Medido en el Co-60 a 1600 V:

| post | carga recogida | eventos útiles | D corregida por carga |
|---|---|---|---|
| 83.0 ns (lo que pide el pulso) | 99.9% | 16% | — |
| 80 ns | 99.5% | 71% | 1.948 |
| **75 ns (con el margen)** | **98.0%** | **94%** | **1.949** |

Los eventos que se pierden son todos de amplitud < 20 mV: por encima de 50 mV la
truncación es 0%.

**Para la próxima toma: adelantar el retardo de trigger 40-50 ns.** Resuelve esto
de raíz.

### 6.2 Para COMPARAR corridas hay que usar la MISMA compuerta

El `post` que mide cada archivo depende de su propio pulso, y dos corridas con
`post` distinto recogen fracciones distintas de la carga — eso entra derecho en
D. Las tres corridas de Cs-137 del 26-ago habían medido 55, 98 y 124 ns (94.5% a
99.6% de carga):

| | ley de ganancia | χ² | residuos |
|---|---|---|---|
| compuerta de cada archivo | b = 10.4 | 744 (1 gdl) | ±14% |
| corregida por fracción de carga | b = 10.0 | — | ±2.3% |
| **compuerta común de 75 ns** | **b = 10.1** | **9.2 (1 gdl)** | **±1.5%** |

Se fuerza con `pipeline.calibrar_completo(c, {"post_fijo": 75.0})`.

### 6.3 Lo demás

* **Sin fondo.** Ninguna corrida binaria tiene fondo propio; el `BG.xml` del
  26-ago no sirve (otro día, otro formato y, para B/C/D, otro detector). El
  espectro sale **crudo** y `espectro_neto` lo dice.
* **Canal `C001` → `CHN1`.** `drs_bin.detectar_canal` traduce a la etiqueta que
  usa el resto del proyecto, para que la clave del caché y los informes hablen el
  mismo idioma. Los XML de junio son `CHN2`.
* **Tiempo vivo.** El Co-60 del 01-sep estuvo vivo el 43.8% (scaler 220 Hz,
  130 ev/s guardados). Las corridas de B, C y D estuvieron vivas el **0.3-0.6%**
  (scaler 25-55 kHz): el umbral del discriminador estaba muy abajo.
* **El detector D tiene 3.4× menos ganancia** que el A, como decía la bitácora.
  Su espectro queda tan comprimido que el umbral del trigger se come buena parte
  del continuo Compton.

---

## 7. Cómo verificar un archivo nuevo

```bash
python analisis/drs.py                  # qué formato tiene cada corrida
python analisis/drs_bin.py              # metadatos + pulsos de los binarios
python analisis/revisar.py <clave>      # revisión completa antes de calibrar
python analisis/visor.py                # pestaña 4 "Archivos", en interfaz
```

Lo que hay que mirar, en orden:

1. Que `estructura()` no levante `ValueError` — si el cuerpo no es múltiplo del
   tamaño de evento, el archivo está truncado o es otra configuración.
2. Que la **suma de los anchos de celda** dé el registro esperado (200 ns).
3. Que la fracción de **eventos que entran** enteros en la compuerta sea alta. Si
   no lo es, §6.1.
4. Que el **tiempo vivo** sea razonable. Con 0.3% el espectro sigue siendo válido
   en FORMA (es lo que calibra), pero las tasas absolutas no.

---

## 8. Dónde está cada cosa

| archivo | qué |
|---|---|
| `analisis/drs_bin.py` | el lector binario |
| `analisis/drs.py` | el despachador xml/binario |
| `analisis/drs_xml.py` | el lector de XML (sin cambios) |
| `analisis/rutas.py` | catálogo: 16 corridas, `SESION_EXTERNA`, `_ext()` |
| `analisis/revisar.py` | revisión de una tanda antes de calibrar |
| `analisis/visor.py` | interfaz para ver la calibración medida |
| `analisis/area.py` | `calibrar_ventana` (recorta la compuerta), `fraccion_de_carga` |
| `analisis/pipeline.py` | opción `post_fijo` para forzar la misma compuerta |
