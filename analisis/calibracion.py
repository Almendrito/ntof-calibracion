"""
calibracion.py — la relacion energia-canal, que es el producto final.

El ajuste entrega D. Este modulo lo convierte en lo que se usa:

    energia [keV]  =  area [mV*ns] / D
    area   [mV*ns] =  D * energia

El mapa pasa por el ORIGEN, sin ordenada: area cero corresponde a energia cero.
No es una eleccion nuestra, es la Ec. 8 del paper (E_eff = E_e * B/E_max) y todo
el metodo depende de ella, incluida la seleccion de intervalo del NCC.

Sobre sigma y sus unidades
--------------------------
`sigma_rel` es ADIMENSIONAL. La convolucion se hace en canal con ancho

    sigma_ch = sigma_rel * D * E_ref     [mV*ns]

y como el mapa es lineal por el origen, eso es lo mismo que un ancho en energia

    sigma_E  = sigma_rel * E_ref         [keV]

Las dos son la misma cantidad, divididas o multiplicadas por D. Se reportan las
dos porque se usan para cosas distintas: sigma_ch para saber cuantos bines del
histograma abarca el borroneado, sigma_E para comparar con otras mediciones de
resolucion.

⚠️ Este modelo tiene una sola sigma CONSTANTE, asi que `sigma_E` es el ancho a
TODAS las energias, no solo en E_ref. Fisicamente la resolucion de un centellador
escala como sqrt(E) (fluctuacion del numero de fotoelectrones), asi que esto es
una aproximacion del modelo y no una medicion de la dependencia con la energia.
Mientras el ajuste viva en el entorno del borde, donde el rango de energias es
angosto, la diferencia es chica; extrapolar `sigma_E` a 100 keV seria inventar.

keVee contra keV
----------------
La plantilla `Lvis_keVee` es la LUZ, con el quenching de Birks ya aplicado, asi
que la D que sale de ella convierte area en keV **electron-equivalentes**. La
plantilla `Edep_keV` es la energia DEPOSITADA. Las dos calibraciones son
legitimas y distintas; hay que decir cual se esta usando. Su cociente reproduce
el quenching de la simulacion, pero eso NO es un control independiente: Lvis es
Edep multiplicada por el quenching, asi que el cociente tiene que dar eso.

Correr `python analisis/calibracion.py` deja `figuras/08_calibracion.png` y
`.cache/calibracion_<corrida>.json`.
"""

import json
import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, fisica, simulacion, modelo, sa


# --------------------------------------------------------------------------
# 0. La receta por defecto, en un solo lugar
# --------------------------------------------------------------------------
# El modelo incluye el exceso de bajo canal y ajusta desde 0.25*B. La decision
# esta medida, no es de gusto (Co-60, plantilla Lvis):
#
#   modelo                       ventana        D       sigma_rel  chi2/ndf
#   desnudo, inicio 0.45      [ 533, 1696]    1.2230     0.1256      1.576
#   desnudo, inicio 0.75      [ 898, 1713]    1.2422     0.1080      1.090
#   con exceso, ajuste 0.25   [ 296, 1696]    1.2512     0.0882      1.169
#
# El modelo desnudo tiene una deriva MONOTONA con el arranque de la ventana, o
# sea que ahi "elegir la ventana" es elegir cuanto sesgo se acepta. Con el
# exceso adentro la deriva se va y D deja de depender de esa eleccion (1.2512
# contra 1.2524 entre 0.25 y 0.35).
#
# Control cruzado contra la implementacion anterior (`codigo/`, modelo
# extendido con la mezcla de scattering `w` libre): daba D = 1.2618 con w
# ajustada y 1.2537 con w fijada en el 0.751 de la simulacion, que es
# exactamente lo que hace este modelo. Coincide al 0.2%.
EXCESO = True
INICIO_AJUSTE_REL = 0.25

# Lo que la semilla del NCC le entrega al recocido NO es un punto de partida
# (los arranques del multiarranque son uniformes dentro de las cotas, a
# proposito: ver `sa.multiarranque`) sino DONDE buscar: D entre 0.6 y 1.6
# veces D0 = B/E_ref. Generosas a proposito, porque D0 sale de un barrido sin
# pesos sobre el espectro suavizado y con un sigma supuesto, y puede venir
# sesgada unos por ciento; una cota estrecha dejaria el optimo afuera sin
# avisar. sigma_rel y lambda no dependen de la semilla: son fisica (una
# resolucion entre 2% y 30%) y geometria (un exceso que decae mas lento que
# media ventana es indistinguible de un fondo plano). Estan aca, con nombre,
# porque el informe los cita.
COTA_D_REL = (0.6, 1.6)          # cotas de D, relativas a la semilla D0
COTA_SIGMA_REL = (0.02, 0.30)    # cotas de sigma_rel
COTA_LAMBDA_MIN = 10.0           # [mV*ns] cota inferior de lambda (exceso)
LAMBDA_SEMILLA_REL = 0.15        # lambda0 = 0.15 * (ancho de la ventana)


# --------------------------------------------------------------------------
# 1. Conversion
# --------------------------------------------------------------------------
def a_energia(area, D):
    """area [mV*ns] -> energia [keV o keVee, segun la plantilla usada]."""
    return np.asarray(area, float) / D


def a_area(energia, D):
    """energia -> area [mV*ns]."""
    return np.asarray(energia, float) * D


def parametros(r):
    """
    El vector de parametros no lineales que corresponde a un resultado.

    Existe para que quien dibuje no tenga que saber si el ajuste tuvo o no
    componente de exceso: `p.prediccion(parametros(r))` da el modelo completo en
    los dos casos, mientras que `r["A"] * p.forma(...)` se olvidaria del exceso.
    """
    x = (r["D"], r["sigma_rel"])
    if r.get("exceso"):
        x = x + (r["lam_exceso"],)
    if r.get("libre_ordenada"):
        x = x + (r["c0"],)
    return x


# --------------------------------------------------------------------------
# 2. Calibrar una corrida
# --------------------------------------------------------------------------
def calibrar(corrida, rama="Lvis_keVee", metodo="global", n_arranques=8,
             verbose=True, esp=None, plantilla=None, inicio_rel=0.45,
             cb_ncc=None, cb_sa=None, exceso=None, tipo_exceso="exp",
             inicio_ajuste_rel="auto", ordenada=False, detalles=None,
             varianza="dato"):
    """
    Cadena completa: NCC -> recocido -> pulido -> errores -> relacion E/canal.

    `corrida` puede ser una clave del catalogo o una `rutas.Corrida` armada al
    vuelo; `esp` y `plantilla` permiten pasar lo ya calculado en vez de releerlo
    del disco (ver `modelo.armar`).

    `exceso` agrega la componente de bajo canal al modelo y libera un tercer
    parametro no lineal (lambda); con ella se baja tambien el arranque de la
    ventana (`inicio_ajuste_rel`), que es de lo que se trata: modelar el exceso
    para poder ajustar donde esta el flanco entero, en vez de esquivarlo
    recortando la ventana. Los dos vienen de EXCESO / INICIO_AJUSTE_REL; pasar
    `exceso=False` reproduce el modelo desnudo de tres parametros.

    `detalles`, si se le pasa un dict, se llena con lo que el ajuste usa por
    dentro y despues se tira: la semilla COMPLETA del NCC, las cotas y las N
    corridas del recocido con su historia. No entra en `r` a proposito, porque
    `r` se serializa a JSON y esto son arrays; y hace falta afuera porque es lo
    unico que permite mostrar el punto de partida, el camino y el punto final
    sin volver a correr el recocido.

    `varianza="modelo"` pesa el chi2 con la varianza PREDICHA (Pearson) en vez
    de la observada (Neyman, el default con que se calibra la campana). No
    cambia el resultado que se reporta: existe para medir el sesgo de Neyman
    con verdad conocida (ver `modelo.Problema`).

    `ordenada=True` libera la ordenada del mapa (`canal = D*E + c0`) en vez de
    imponer el origen. NO es como se calibra: es el TEST de la hipotesis del
    origen con una sola fuente, y hay que leerlo como cota (c0 se confunde con
    el exceso de bajo canal). Ver `test_ordenada`.

    Devuelve un dict con todo lo que hace falta para usar y para auditar la
    calibracion. El error de D se toma de la curvatura del chi2 y se ESCALA por
    sqrt(chi2/ndf): con chi2/ndf > 1 el criterio Delta_chi2 = 1 subestima las
    barras, porque el modelo no describe todo lo que hay en los datos. Es un
    parche estandar y hay que reportarlo como tal.
    """
    if exceso is None:
        exceso = EXCESO
    if inicio_ajuste_rel == "auto":
        inicio_ajuste_rel = INICIO_AJUSTE_REL if exceso else None

    c = rutas.resolver(corrida)
    p, semilla, ncc, e, _ = modelo.armar(c, metodo=metodo, rama=rama,
                                         verbose=False, esp=esp,
                                         plantilla=plantilla,
                                         inicio_rel=inicio_rel, cb_ncc=cb_ncc,
                                         exceso=exceso,
                                         tipo_exceso=tipo_exceso,
                                         inicio_ajuste_rel=inicio_ajuste_rel,
                                         ordenada=ordenada, varianza=varianza)
    cotas = [(COTA_D_REL[0] * semilla[0], COTA_D_REL[1] * semilla[0]),
             COTA_SIGMA_REL]
    if exceso:
        # lambda: escala de caida del exceso [mV*ns]. La cota superior la fija
        # el ancho de la ventana: un exceso que decae mas lento que eso es
        # indistinguible de un fondo plano y se lo come la plantilla.
        cotas.append((COTA_LAMBDA_MIN, max(0.5 * (p.hi - p.lo), 50.0)))
        semilla = tuple(semilla) + (LAMBDA_SEMILLA_REL * (p.hi - p.lo),)
    if ordenada:
        # c0 se acota a +-40% del borde: mas que eso no es un pedestal, es otro
        # ajuste. La semilla es 0, o sea la hipotesis que se esta testeando.
        B0 = semilla[0] * fisica.e_ref(rutas.resolver(corrida).fuente)
        cotas.append((-0.4 * B0, 0.4 * B0))
        semilla = tuple(semilla) + (0.0,)

    # multiarranque: una corrida del recocido es una muestra, no un resultado
    corridas = sa.multiarranque(p.costo, cotas, n=n_arranques, semilla=7,
                                callback=cb_sa)
    j = int(np.argmin([r["costo_pulido"] for r in corridas]))
    x = corridas[j]["x_pulido"]
    chi2 = corridas[j]["costo_pulido"]
    Ds_ms = np.array([r["x_pulido"][0] for r in corridas])

    if detalles is not None:
        detalles.update(semilla=tuple(float(v) for v in semilla),
                        cotas=[(float(a), float(b)) for a, b in cotas],
                        libres=list(p.libres),
                        chi2_semilla=float(p.costo(semilla)),
                        corridas=corridas, mejor=j, ncc=ncc)

    err, C = sa.errores_por_curvatura(p.costo, x, paso=p.pasos_curvatura(x))
    chi2_ndf = chi2 / p.ndf
    escala = np.sqrt(max(chi2_ndf, 1.0))

    D, s_rel = float(x[0]), float(x[1])
    E_ref = fisica.e_ref(c.fuente)
    bordes_E = fisica.bordes(c.fuente)

    r = dict(
        corrida=c.clave, fuente=c.fuente, detector=c.detector, rama=rama,
        metodo=metodo,
        # --- el punto de operacion: sin esto D no significa nada ---
        punto=c.punto(), hv=c.hv, atenuador=c.atenuador, cadena=c.cadena(),
        # --- el vector Theta completo, tal cual lo ve el optimizador ---
        # (D y sigma_rel se repiten abajo con nombre propio por compatibilidad
        # con lo que ya lee `informe`/el registro; esto es la version generica
        # que permite dibujar la tabla y la matriz de correlacion SIN saber de
        # antemano cuantos parametros tiene el modelo -2, 3 o 4, segun exceso/
        # ordenada).
        libres=list(p.libres),
        x_libres=[float(v) for v in x],
        err_libres=[float(v) for v in (err * escala)],
        cov_libres=(np.asarray(C, float) * escala ** 2).tolist(),
        # --- el resultado ---
        D=D, err_D=float(err[0] * escala),
        inv_D=1.0 / D, err_inv_D=float(err[0] * escala / D ** 2),
        # --- resolucion, en las dos unidades ---
        sigma_rel=s_rel, err_sigma_rel=float(err[1] * escala),
        sigma_E=s_rel * E_ref,                 # [keV o keVee]
        sigma_ch=s_rel * D * E_ref,            # [mV*ns]
        fwhm_E=2.3548 * s_rel * E_ref,
        fwhm_rel=2.3548 * s_rel,
        # --- donde caen los bordes ---
        E_ref=E_ref,
        bordes_E=list(bordes_E),
        bordes_ch=[D * b for b in bordes_E],
        err_bordes_ch=[float(err[0] * escala * b) for b in bordes_E],
        # --- calidad ---
        chi2=chi2, ndf=p.ndf, chi2_ndf=chi2_ndf, escala_error=float(escala),
        corr_D_sigma=float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])),
        dispersion_multiarranque=float(Ds_ms.std()),
        n_arranques=n_arranques,
        ventana=[float(p.lo), float(p.hi)],
        ventana_ncc=[float(ncc["lo"]), float(ncc["hi"])], rho=float(ncc["rho"]),
        semilla_D=float(semilla[0]),
        A=float(p.chi2_de(x)[1]),
        # --- el exceso de bajo canal (0 si no se modelo) ---
        exceso=bool(exceso), tipo_exceso=tipo_exceso if exceso else None,
        lam_exceso=float(x[2]) if exceso else None,
        err_lam_exceso=float(err[2] * escala) if exceso else None,
        C_exceso=float(p.C_exceso),
        # --- el test del origen (0 y sin error si no se libero) ---
        libre_ordenada=bool(ordenada),
        c0=float(x[p.libres.index("c0")]) if ordenada else 0.0,
        err_c0=float(err[p.libres.index("c0")] * escala) if ordenada else None,
    )
    if verbose:
        informe(r)
    # se devuelve tambien el dict del NCC: la interfaz dibuja el barrido de rho
    # como diagnostico, y recalcularlo seria repetir el trabajo caro del paso 5.
    return r, p, e, ncc


def test_ordenada(corrida, rama="Lvis_keVee", metodo="global", n_arranques=8,
                  esp=None, plantilla=None, verbose=True, **kw):
    """
    El test de la hipotesis del origen CON UNA SOLA FUENTE.

    Ajusta dos veces, con la ordenada fija en cero (como se calibra) y con la
    ordenada libre, y compara. Lo que hace posible el test es que el continuo
    Compton va de 0 a E_max: la forma entera lleva informacion sobre donde cae
    el cero del eje, aunque haya un solo borde.

    Es un test DEBIL y hay que decirlo: c0 se confunde con el exceso de bajo
    canal (los dos agregan cuentas abajo) y con sigma. Sirve como cota -"si
    hubiera un pedestal de mas de tanto, se veria"- y no reemplaza al test con
    dos fuentes, que es el limpio. Un Delta_chi2 chico NO prueba que c0 = 0:
    prueba que los datos no lo distinguen.
    """
    base = dict(rama=rama, metodo=metodo, n_arranques=n_arranques,
                esp=esp, plantilla=plantilla, verbose=False)
    base.update(kw)
    r0, p0, _, _ = calibrar(corrida, ordenada=False, **base)
    r1, p1, _, _ = calibrar(corrida, ordenada=True, **base)
    B = r0["D"] * r0["E_ref"]
    dchi2 = r0["chi2"] - r1["chi2"]
    # La significancia se toma del COCIENTE DE VEROSIMILITUD (sqrt del Delta
    # chi2 por un parametro), no de la curvatura local. Los dos no coinciden
    # -Co-60: 2.1 sigma contra 4.7- y la diferencia no es un detalle numerico:
    # dice que el chi2 es mucho mas PLANO lejos del minimo de lo que su
    # curvatura sugiere, o sea que c0 esta degenerado con el exceso de bajo
    # canal. Es la degeneracion que hace debil a este test, medida.
    n_sig = float(np.sqrt(max(dchi2, 0.0)))
    n_sig_curv = abs(r1["c0"]) / max(r1["err_c0"], 1e-9)
    res = dict(D_origen=r0["D"], err_D_origen=r0["err_D"],
               D_libre=r1["D"], err_D_libre=r1["err_D"],
               c0=r1["c0"], err_c0=r1["err_c0"], c0_rel=r1["c0"] / B,
               chi2_origen=r0["chi2"], chi2_libre=r1["chi2"],
               ndf_origen=r0["ndf"], ndf_libre=r1["ndf"],
               delta_chi2=dchi2, sigmas=n_sig, sigmas_curvatura=n_sig_curv,
               desplazamiento_D=abs(r1["D"] - r0["D"]) / r0["D"], borde=B,
               compatible=bool(n_sig < 2.0))
    if verbose:
        print("  TEST DEL ORIGEN con una sola fuente (%s, %s)"
              % (r0["corrida"], rama))
        print("     origen fijo  : D = %.4f +- %.4f   chi2 = %.1f (%d gdl)"
              % (r0["D"], r0["err_D"], r0["chi2"], r0["ndf"]))
        print("     ordenada libre: D = %.4f +- %.4f   chi2 = %.1f (%d gdl)"
              % (r1["D"], r1["err_D"], r1["chi2"], r1["ndf"]))
        print("     c0 = %+.1f +- %.1f mV*ns  =  %+.1f%% del borde (%.0f mV*ns)"
              % (r1["c0"], r1["err_c0"], 100 * res["c0_rel"], B))
        print("     Delta_chi2 = %.1f por 1 parametro mas  ->  %.1f sigma "
              "(la curvatura local diria %.1f: el chi2 es mas plano de lo que "
              "parece cerca del minimo)" % (dchi2, n_sig, n_sig_curv))
        print("     liberar c0 mueve D un %.1f%%: esa es la degeneracion"
              % (100 * res["desplazamiento_D"]))
        print("     -> %s" % (
            "compatible con el origen: liberar c0 no mejora el ajuste"
            if res["compatible"] else
            "NO compatible con cero: hay indicio de pedestal, revisar antes de "
            "usar la calibracion por el origen"))
        print("     (test debil: c0 se confunde con el exceso de bajo canal; "
              "el limpio necesita dos fuentes)")
    return res, r0, r1


def informe(r):
    u = "keVee" if r["rama"].startswith("Lvis") else "keV"
    print("-" * 76)
    print("  %s  (%s, %s)   plantilla %s  ->  escala en %s"
          % (r["corrida"], r["fuente"], r.get("punto", r["detector"]),
             r["rama"], u))
    if r.get("cadena"):
        print("  cadena de senal: %s" % r["cadena"])
    print("-" * 76)
    print("  RELACION ENERGIA-CANAL   (por el origen, sin ordenada)")
    print("     D          = %.4f +- %.4f  mV*ns/%s   (%.2f %%)"
          % (r["D"], r["err_D"], u, 100 * r["err_D"] / r["D"]))
    print("     1/D        = %.4f +- %.4f  %s/(mV*ns)"
          % (r["inv_D"], r["err_inv_D"], u))
    print("     o sea      : 1000 mV*ns  <->  %.1f %s" % (1000 * r["inv_D"], u))
    print("  BORDES COMPTON")
    for E, ch, ech in zip(r["bordes_E"], r["bordes_ch"], r["err_bordes_ch"]):
        print("     %7.1f %-5s  ->  %7.1f +- %.1f mV*ns" % (E, u, ch, ech))
    print("  RESOLUCION  (una sola sigma CONSTANTE en este modelo)")
    print("     sigma_rel  = %.4f +- %.4f   (adimensional)"
          % (r["sigma_rel"], r["err_sigma_rel"]))
    print("     sigma_E    = %.1f %s        sigma_ch = %.1f mV*ns"
          % (r["sigma_E"], u, r["sigma_ch"]))
    print("     FWHM       = %.1f %s   ->  %.1f %% en E_ref = %.1f %s"
          % (r["fwhm_E"], u, 100 * r["fwhm_rel"], r["E_ref"], u))
    print("  CALIDAD")
    print("     chi2/ndf   = %.3f  (%d gdl)   ->  errores escalados x%.2f"
          % (r["chi2_ndf"], r["ndf"], r["escala_error"]))
    print("     corr(D, sigma) = %+.2f   (reportarlos juntos)"
          % r["corr_D_sigma"])
    print("     multiarranque: %d corridas, dispersion en D de %.4f (%.3f %%)"
          % (r["n_arranques"], r["dispersion_multiarranque"],
             100 * r["dispersion_multiarranque"] / r["D"]))
    vn = r.get("ventana_ncc", r["ventana"])
    print("     ventana de ajuste [%.0f, %.0f] mV*ns   (NCC eligio [%.0f, %.0f],"
          " rho = %.5f, semilla D = %.4f)"
          % (r["ventana"][0], r["ventana"][1], vn[0], vn[1], r["rho"],
             r["semilla_D"]))
    if r.get("exceso"):
        print("  EXCESO DE BAJO CANAL  (modelado, forma '%s')" % r["tipo_exceso"])
        print("     lambda     = %.1f +- %.1f mV*ns   amplitud C = %.3g cuentas"
              % (r["lam_exceso"], r["err_lam_exceso"], r["C_exceso"]))


def guardar(r):
    # El metodo de area va en el nombre SOLO si no es el global. Sin el sufijo,
    # calibrar con `metodo="plantilla"` pisaria en silencio el JSON que lee
    # `informe.py`; con el sufijo siempre, los archivos del metodo por defecto
    # cambiarian de nombre y el informe dejaria de encontrarlos.
    sufijo = "" if r.get("metodo", "global") == "global" else "_" + r["metodo"]
    ruta = rutas.cache("calibracion_%s_%s%s.json"
                       % (r["corrida"], "lvis" if r["rama"].startswith("Lvis")
                          else "edep", sufijo))
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(r, f, indent=2, ensure_ascii=False)
    return ruta


def quenching_simulacion(ruta_root, e_ref, frac=0.7):
    """
    <Lvis/Edep> de la simulacion en la zona que domina el ajuste, definida como
    la energia por encima de `frac * e_ref`.

    El umbral tiene que ir atado al borde de CADA fuente, no fijo: con un umbral
    fijo de 700 keV el Cs-137 (borde en 477) no tiene un solo evento arriba y el
    promedio sale nan. Ademas el quenching de Birks no es constante -satura
    recien arriba de ~700 keV-, asi que el numero que corresponde comparar es el
    de la zona donde efectivamente esta el borde de esa fuente.
    """
    import uproot
    d = uproot.open(ruta_root)[simulacion.ARBOL].arrays(
        ["Edep_keV", "Lvis_keVee"], library="np")
    m = d["Edep_keV"] > frac * e_ref
    if m.sum() < 100:
        return float("nan")
    return float(np.mean(d["Lvis_keVee"][m] / d["Edep_keV"][m]))


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def _demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rutas.salida_utf8()
    claves = [k for k in rutas.PRINCIPALES
              if os.path.exists(rutas.cache("espectro_%s_global.npz" % k))]

    fig, axes = plt.subplots(len(claves), 3, figsize=(16.5, 4.4 * len(claves)),
                             squeeze=False)
    resumen = []

    for fila, clave in enumerate(claves):
        c = rutas.corrida(clave)
        print()
        print("=" * 76)
        print("%s  -  %s" % (clave, c.descripcion))
        print("=" * 76)

        r_luz, p, e, _ = calibrar(clave, rama="Lvis_keVee")
        guardar(r_luz)
        r_dep, _, _, _ = calibrar(clave, rama="Edep_keV")
        guardar(r_dep)

        q = quenching_simulacion(c.sim, r_luz["E_ref"])
        print()
        print("  Lvis vs Edep: D_luz/D_edep = %.4f ; el quenching medio de la "
              "simulacion arriba de %.0f keV es %.4f, o sea 1/q = %.4f"
              % (r_luz["D"] / r_dep["D"], 0.7 * r_luz["E_ref"], q, 1 / q))
        print("     ⚠ NO es un control independiente: Lvis ES Edep por el "
              "quenching, asi que el cociente TIENE que dar eso. Lo que verifica")
        print("       es que el ajuste esta dominado por la zona donde el "
              "quenching satura, que es util pero es menos de lo que parece.")
        resumen.append((clave, r_luz, r_dep))

        # ---------------------------------------------------------- figura --
        D, u = r_luz["D"], "keVee"

        ax = axes[fila][0]
        Emax = 1.35 * max(r_luz["bordes_E"])
        Es = np.linspace(0, Emax, 50)
        ax.plot(Es, a_area(Es, D), lw=1.8, color="C0",
                label="ch = %.4f · E" % D)
        ax.fill_between(Es, a_area(Es, D - r_luz["err_D"]),
                        a_area(Es, D + r_luz["err_D"]), color="C0", alpha=.25)
        for E, ch in zip(r_luz["bordes_E"], r_luz["bordes_ch"]):
            ax.plot([E], [ch], "o", color="C3", ms=7, zorder=4)
            ax.annotate("%.1f %s\n%.0f mV·ns" % (E, u, ch), xy=(E, ch),
                        xytext=(E * 0.42, ch * 1.02), fontsize=8, color="C3")
        ax.plot([0], [0], "s", color="k", ms=5)
        ax.annotate("el mapa pasa por el ORIGEN", xy=(0, 0),
                    xytext=(Emax * 0.05, a_area(Emax, D) * 0.12), fontsize=8)
        ax.set_title("%s: relacion energia-canal" % clave)
        ax.set_xlabel("energia [%s]" % u); ax.set_ylabel("area [mV·ns]")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

        ax = axes[fila][1]
        ax.step(e["canales"], e["neto"], where="mid", lw=1, color="0.55",
                label="espectro neto")
        ax.plot(p.canales, p.prediccion(parametros(r_luz)),
                lw=1.6, color="C1", label="modelo ajustado")
        tope = float(np.max(e["neto"][e["canales"] > 300]))
        for k, (E, ch, ech) in enumerate(zip(r_luz["bordes_E"],
                                             r_luz["bordes_ch"],
                                             r_luz["err_bordes_ch"])):
            ax.axvline(ch, color="C3", ls="--", lw=1.2)
            ax.axvspan(ch - ech, ch + ech, color="C3", alpha=.25)
            # escalonadas: con dos bordes a 190 mV*ns las etiquetas se pisan
            ax.annotate("%.0f mV·ns\n(%.0f %s)" % (ch, E, u), xy=(ch, 0),
                        xytext=(ch + 20, (0.85 - 0.22 * k) * tope),
                        fontsize=8, color="C3")
        ax.set_xlim(0, 1.55 * max(r_luz["bordes_ch"]))
        ax.set_title("donde caen los bordes, en mV·ns")
        ax.set_xlabel("area [mV·ns]"); ax.set_ylabel("cuentas")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

        ax = axes[fila][2]
        En = a_energia(e["canales"], D)
        ax.step(En, e["neto"] * D, where="mid", lw=1, color="0.55",
                label="medido, calibrado")
        ax.plot(a_energia(p.canales, D),
                p.prediccion(parametros(r_luz)) * D, lw=1.6,
                color="C1", label="modelo")
        for E in r_luz["bordes_E"]:
            ax.axvline(E, color="C3", ls="--", lw=1.2)
            ax.annotate("%.1f %s" % (E, u), xy=(E, 0),
                        xytext=(E + 8, 0.75 * np.max(e["neto"][e["canales"] > 300]) * D),
                        fontsize=8, color="C3", rotation=90)
        # ancho de la resolucion, dibujado a escala
        y0 = 0.35 * np.max(e["neto"][e["canales"] > 300]) * D
        ax.errorbar([r_luz["E_ref"]], [y0],
                    xerr=[[r_luz["sigma_E"]], [r_luz["sigma_E"]]],
                    fmt="o", color="C0", capsize=4,
                    label="±1σ = %.0f %s" % (r_luz["sigma_E"], u))
        ax.set_xlim(0, 1.55 * max(r_luz["bordes_E"]))
        ax.set_title("EL ESPECTRO CALIBRADO")
        ax.set_xlabel("energia [%s]" % u); ax.set_ylabel("cuentas por keV")
        ax.legend(fontsize=8); ax.grid(alpha=.3)

    # -------------------------------------------------------------- tabla --
    print()
    print("=" * 76)
    print("RESUMEN")
    print("=" * 76)
    print("  %-12s %-22s %-11s %-11s %-9s %s"
          % ("corrida", "detector", "D [keVee]", "D [keV]", "FWHM", "bordes [mV·ns]"))
    for clave, rl, rd in resumen:
        print("  %-12s %-22s %-11s %-11s %-9s %s"
              % (clave, (rl["detector"] or "")[:22],
                 "%.4f" % rl["D"], "%.4f" % rd["D"],
                 "%.1f%%" % (100 * rl["fwhm_rel"]),
                 ", ".join("%.0f" % b for b in rl["bordes_ch"])))
    print()
    print("  ⚠ Las dos corridas son de DETECTORES DISTINTOS: no son dos puntos")
    print("    de una misma recta y no se pueden combinar. Cada una es su propia")
    print("    calibracion.")
    print("  ⚠ Los dos bordes del Co-60 tampoco son dos puntos independientes:")
    print("    el mapa por el origen los ata a la misma D (es la Ec. 12 del")
    print("    paper). Para TESTEAR linealidad hace falta una segunda fuente")
    print("    medida con el MISMO detector.")

    fig.tight_layout()
    salida = rutas.figura("08_calibracion.png")
    fig.savefig(salida, dpi=130)
    print()
    print("figura -> %s" % rutas.relativa(salida))


if __name__ == "__main__":
    _demo()
