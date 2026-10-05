"""
calibracion.py - EL RESULTADO: como convertir area de pulso [mV*ns] a energia.

Detector : centellador plastico + SiPM MicroFC-30020, leido con DRS4 (5 GS/s,
           1024 muestras = 200 ns), atenuador x10. Observable = area del pulso
           integrada en una compuerta alrededor del pulso, en mV*ns.
Fuente   : Co-60 (dos bordes Compton, 963.4 y 1118.1 keV), 283 s VIVOS de 699 s
           de pared, 45125 pulsos buenos, con fondo restado por tiempo vivo.
Fecha    : 25-ago-2026 (Etapa 9).

⚠️ EL NUMERO CAMBIO POR UN FACTOR 3.4 RESPECTO DE LA ETAPA 8
------------------------------------------------------------
Hasta la Etapa 8 este archivo decia D_LUZ = 4.32 mV*ns/keVee. Ese numero estaba
mal, y no por el metodo de ajuste sino por el OBSERVABLE: el area del pulso se
calculaba restando una linea base promediada sobre los primeros 40 ns del
registro, y el pulso arranca a los ~25 ns. O sea la linea base se comia la
subida del propio pulso. Con eso:

  - el 26% de los eventos salia con el area de signo INVERTIDO (y un abs() los
    reciclaba como buenos),
  - la dispersion de area/amplitud a energia fija era std = 11.3 ns contra 0.84
    con una linea base limpia,
  - y como el fondo se habia medido otro dia con el pulso en OTRA posicion del
    registro (120 ns en vez de 48), muestra y fondo ni siquiera compartian el
    observable, asi que la resta era invalida.

Ademas la resta de fondo se normalizaba por tiempo de PARED cuando la corrida de
Co-60 estuvo viva solo el 40% del tiempo (el DRSOSC no alcanza a volcar el XML a
170 Hz de trigger). Eso sobre-restaba el fondo por un factor 2.5.

Los dos arreglos estan en `spectrum.py` (`_area_pulso_gate` y `tiempo_vivo`).
Ver CLAUDE.md 7.33-7.35. Cualquier analisis anterior que use mV*ns hay que
rehacerlo: el factor entre la escala vieja y la nueva NO es una constante
limpia, porque la corrupcion dependia del jitter del trigger evento a evento.

DOS ESCALAS, no confundirlas
----------------------------
El SiPM ve LUZ, no energia. La simulacion Geant4 da Lvis/Edep = 0.954 sobre
700 keV para electrones (quenching de Birks), asi que:

  D_LUZ  = 1.261 mV*ns por keVee -> escala de LUZ (keV electron equivalente).
                                    Es la que sirve para la rama de neutrones:
                                    los protones de retroceso van en keVee.
  D_EDEP = 1.208 mV*ns por keV   -> escala de ENERGIA DEPOSITADA por gammas.

El cociente medido, 1.044, reproduce el 1/0.954 = 1.048 de la simulacion.

⚠️ Ojo con cuanto vale ese "control": las dos plantillas NO son independientes.
Lvis es Edep multiplicada por el quenching, y arriba de 700 keV el quenching es
plano en 0.954, asi que ajustar los mismos datos con una y con otra TIENE que
dar ese cociente. Verifica que el ajuste esta dominado por la zona donde el
quenching satura (que es lo que uno quiere, porque los bordes estan en 963 y
1118 keV), pero no es una comprobacion independiente de la calibracion.
El proyecto lo venia llamando "control de consistencia independiente"; no lo es.

Como se obtuvo (Etapa 9)
------------------------
1. AREA con compuerta: onset por evento (cruce del 10% de la amplitud
   retrocediendo desde el pico), linea base = mediana de las muestras PRE-onset
   (hasta 4 ns antes del inicio de la compuerta), integracion en
   [onset-8, onset+55] ns. El pulso medio devuelve el 100% de su area en 80 ns,
   asi que la compuerta no corta nada (verificado entre 20 y 220 mV de
   amplitud: area/amplitud constante al 1%).
   Los eventos sin pre-onset suficiente NO se descartan: se les toma la linea
   base del complemento de la compuerta y se marcan. Descartarlos seria un corte
   dependiente de la amplitud, porque el trigger tiene walk (los pulsos grandes
   cruzan el umbral ~12 ns antes que los chicos).
2. TIEMPO VIVO desde el `<Scaler1>` de la placa (tasa de trigger real).
3. MODELO EXTENDIDO (`modelo_ext.py`): la plantilla es la simulacion Geant4 con
   la mezcla scattering-simple/multiple AJUSTABLE (parametro w), mas una
   componente empirica para el exceso a bajo canal. Los dos sistematicos que la
   Etapa 8 trataba a mano son ahora parametros, y su ignorancia entra en el
   error de D por la correlacion (corr(D,w) = +0.95).
4. POSTERIOR por MCMC, 4 cadenas x 9000 pasos, R_hat < 1.002.

   chi2/ndf = 1.008 (314 gdl), ventana [300, 1900] mV*ns.

Lo que el arreglo del area resolvio de paso
-------------------------------------------
  - La resolucion cayo de 42.7% a 21.1% FWHM a 963 keV. El aviso de la Etapa 8
    ("43% es alto incluso para plastico+SiPM") tenia razon: era instrumental.
  - El exceso a bajo canal se achico y se corrio; ahora vive abajo de ~800
    mV*ns (0.65 del borde) y se puede MODELAR en vez de cortar, asi que el
    LIMITE de la ventana ya no es un sistematico dominante (D varia 1.5% entre
    lo=250 y 800, sin deriva monotona; antes derivaba 3.5%). Lo que si quedo
    dominando es la FORMA que se le asume a ese exceso.
  - El ruido electronico quedo medido: integrando una ventana VACIA de 75 ns el
    area da 0 +- 5 mV*ns, o sea 4 keV. No explica nada del ancho.

Incertidumbre
-------------
  estadistica (posterior MCMC, ya incluye w y el exceso) : +-0.0135  (1.07%)
  limite inferior de la ventana (250-800 mV*ns)          : +-0.0095  (0.75%)
  definicion del area (compuerta, onset, descarte)       : +-0.0131  (1.04%)
  forma del exceso (exponencial vs ley de potencia)      : +-0.0296  (2.35%)
  --> D_LUZ = 1.261 +- 0.036 mV*ns/keVee  (2.9% total)

Manda el sistematico, 2.7x el estadistico, y ahora se sabe CUAL: la forma que se
le asume al exceso a bajo canal. Exponencial y ley de potencia describen los
datos igual de bien (chi2/ndf 1.008 contra 1.015, indistinguibles) y difieren
2.35% en D. Eso NO se cierra ajustando mejor: hay que medir el exceso. Una
corrida con la fuente apantallada, o movida de lugar, dice si viene de
dispersion en el material de alrededor.

El barrido de definiciones del area (compuerta de 45/55/70/90 ns, onset al
5/10/20%, guarda, con y sin descarte por calidad) da D entre 1.2536 y 1.2799.
La de 45 ns es la unica con chi2 claramente peor (1.076): trunca el pulso.

⚠️ El total, 2.9%, es parecido al 2.9% que reportaba la Etapa 8, pero no
significa que no se haya ganado nada: aquel 2.9% estaba sobre un valor central
equivocado y no incluia NINGUNO de los dos sistematicos que hoy dominan (la
definicion del area no se habia mirado, y del exceso solo se sabia que habia que
esquivarlo). Lo que cambio es que ahora el presupuesto esta completo.

Corre solo:  python calibracion.py
"""

import numpy as np

# ---------------------------- LA CALIBRACION ----------------------------
# Mediana del posterior (Etapa 9, `calibrar_v2.py`). Entre corchetes, el valor
# de la Etapa 8, que salia del area con la linea base contaminada.
D_LUZ = 1.261          # mV*ns por keVee  (escala de luz; plantilla Lvis)  [antes 4.320]
D_LUZ_ERR_STAT = 0.0135  # posterior MCMC, semiancho del 68%
D_LUZ_ERR_SIS = 0.0337   # ventana + definicion del area + forma del exceso
D_LUZ_ERR = 0.036        # total

D_EDEP = 1.208         # mV*ns por keV    (energia depositada por gammas)  [antes 4.133]
D_EDEP_ERR_STAT = 0.0142
D_EDEP_ERR_SIS = 0.0323
D_EDEP_ERR = 0.035

# resolucion medida: sigma_E(E) = ALPHA * sqrt(E)  [keV], E en keV.
# El ajuste probo tambien los terminos a*E (colecta) y beta (electronica) y los
# dos se van a cero: con este alcance en energia la raiz sola alcanza.
ALPHA = 2.787          # sqrt(keV)   -> 21.1% FWHM a 963 keV               [antes 5.626]
ALPHA_ERR = 0.163

# Fraccion de scattering SIMPLE que piden los datos. La simulacion trae 0.751,
# o sea los datos y la simulacion COINCIDEN (0.6 sigma). Fijar w = 0.751 no
# empeora el ajuste (dchi2 = 0); fijar w = 1, que es lo que hace la plantilla
# analitica del paper, lo empeora en dchi2 = +13 y sube D un 3.4%.
W_SIMPLE = 0.796
W_SIMPLE_ERR = 0.070

# Correlaciones del posterior: si se propaga uno, van los tres.
CORR_D_ALPHA = -0.615
CORR_D_W = +0.962
# ------------------------------------------------------------------------


def area_a_energia(area_mVns, escala="luz"):
    """
    Convierte area de pulso [mV*ns] a energia.

    escala="luz"  -> keVee (keV electron equivalente). Default.
    escala="edep" -> keV de energia depositada por gammas.
    """
    D = D_LUZ if escala == "luz" else D_EDEP
    return np.asarray(area_mVns, float) / D


def energia_a_area(energia_keV, escala="luz"):
    """Inversa de area_a_energia."""
    D = D_LUZ if escala == "luz" else D_EDEP
    return np.asarray(energia_keV, float) * D


def resolucion_sigma(energia_keV):
    """sigma de la resolucion en energia [keV] a una energia dada."""
    return ALPHA * np.sqrt(np.maximum(np.asarray(energia_keV, float), 0.0))


def resolucion_fwhm_rel(energia_keV):
    """Resolucion relativa FWHM (fraccion, no %) a una energia dada."""
    E = np.maximum(np.asarray(energia_keV, float), 1e-9)
    return 2.355 * resolucion_sigma(E) / E


if __name__ == "__main__":
    print("=" * 68)
    print("CALIBRACION DEL CENTELLADOR PLASTICO + SiPM   (Etapa 9)")
    print("=" * 68)
    print(f"  D (luz)  = {D_LUZ:.3f} +- {D_LUZ_ERR_STAT:.3f} (est) "
          f"+- {D_LUZ_ERR_SIS:.3f} (sist) mV*ns/keVee"
          f"   -> {1000/D_LUZ:.0f} keVee por 1000 mV*ns")
    print(f"  D (edep) = {D_EDEP:.3f} +- {D_EDEP_ERR_STAT:.3f} (est) "
          f"+- {D_EDEP_ERR_SIS:.3f} (sist) mV*ns/keV "
          f"   -> {1000/D_EDEP:.0f} keV   por 1000 mV*ns")
    print(f"  resolucion: sigma_E = {ALPHA:.2f}*sqrt(E)  [keV]"
          f"   -> {100*resolucion_fwhm_rel(963.4):.1f}% FWHM a 963 keV")
    print(f"  mezcla de scattering: w = {W_SIMPLE:.3f} +- {W_SIMPLE_ERR:.3f} "
          f"simple  (la simulacion trae 0.751)")
    print(f"  corr(D,alpha) = {CORR_D_ALPHA:+.2f}   corr(D,w) = {CORR_D_W:+.2f}")

    print("\n  area [mV*ns] ->  E [keVee]   E [keV edep]   FWHM(E)")
    for a in (100, 250, 500, 750, 1000, 1215, 1410, 1800):
        e_luz = area_a_energia(a, "luz")
        e_dep = area_a_energia(a, "edep")
        print(f"  {a:9.0f}      {e_luz:9.1f}     {e_dep:9.1f}      "
              f"{100*resolucion_fwhm_rel(e_luz):5.1f}%")
    print(f"\n  los dos bordes Compton del Co-60 caen en "
          f"{energia_a_area(963.4):.0f} y {energia_a_area(1118.1):.0f} mV*ns")

    print("\n  --- que cambio respecto de la Etapa 8 ---")
    print("  D_LUZ paso de 4.320 a 1.261 mV*ns/keVee (factor 3.43). NO es un")
    print("  ajuste distinto: es que el AREA estaba mal calculada. La linea base")
    print("  se promediaba sobre los primeros 40 ns y el pulso empieza a los 25,")
    print("  asi que se restaba parte del propio pulso, multiplicada por 5 al")
    print("  integrar los 200 ns completos. El 26% de los eventos quedaba con el")
    print("  area de signo invertido. Ver CLAUDE.md 7.33.")
    print("  De paso cayeron dos avisos de la Etapa 8:")
    print("    - la resolucion (42.7% -> 21.1% FWHM) era instrumental, no fisica;")
    print("    - el sistematico de plantilla (2.7%) ya no se estima comparando")
    print("      tres ajustes: la mezcla w la miden los datos.")

    print("\n  AVISO: queda el exceso de cuentas a bajo canal, ahora abajo de")
    print("  ~800 mV*ns. Se modela con una exponencial y por eso ya no obliga a")
    print("  cortar la ventana, pero su FORMA es el sistematico mas grande que")
    print("  queda: exponencial y ley de potencia ajustan igual de bien y dan D")
    print("  con 2.35% de diferencia. Para cerrarlo hay que MEDIRLO: una corrida")
    print("  con la fuente apantallada o movida de lugar dice si es dispersion")
    print("  en el material de alrededor.")
