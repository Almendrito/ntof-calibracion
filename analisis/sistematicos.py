"""
sistematicos.py — el presupuesto de error de una calibracion, medido.

Por que hace falta
------------------
`calibracion.calibrar` devuelve `err_D`, que es estadistico (curvatura del chi2,
escalada por sqrt(chi2/ndf)). Sobre estas corridas vale ~0.3%, y las decisiones
de analisis mueven D varias veces mas que eso. Un numero con una barra que solo
cuenta los eventos es un numero con una barra equivocada.

Y para los tests de linealidad la distincion importa el doble, porque no todos
los sistematicos pesan igual:

  * un sistematico COMUN a todas las corridas (la escala de la plantilla, la
    definicion del area si es la misma en todas) se va en la PENDIENTE de la
    recta y no ensucia los residuos: mueve D pero no la linealidad;
  * uno PUNTO A PUNTO (la ventana, que la elige el NCC corrida por corrida; la
    forma del exceso, que puede ser distinta si cambia la geometria) es el que
    hay que comparar contra los residuos del test A.

Por eso cada variante se declara con su alcance y el informe suma las dos cosas
por separado.

Como se mide
------------
Ajustando otra vez con la decision cambiada y mirando cuanto se mueve D. Ninguna
de las variantes es "la mala": son elecciones defendibles, y la dispersion entre
ellas es la ignorancia que queda.

  ventana          arranque del ajuste en 0.25, 0.30, 0.35 y 0.45 de B
  forma del exceso exponencial contra ley de potencia
  gate del area    la ventana de integracion movida (solo con `gates=True`)

Dos variantes se miden pero NO se suman, y por razones distintas:

  * **metodo del area, global contra convergencia.** El chi2 no las distingue
    (sobre el Cs-137 dan 0.965 y 0.95 con D distintas en 5.4%), asi que por
    chi2 parecerian dos opciones igual de buenas. No lo son: el metodo por
    convergencia corta la integral con un umbral que depende de la amplitud, y
    esta MEDIDO que por eso su ventana va de 13 ns en los pulsos chicos a 42 en
    los grandes (hallazgo 4 del README). Un observable cuyo ancho de integracion
    depende de la amplitud no es proporcional a la energia, que es justo lo que
    la calibracion necesita. Se descarta por fisica, no por chi2, y meterlo como
    sistematico seria pagar dos veces por un error que ya se sabe cual es.
  * **mezcla de scattering.** La plantilla de scattering simple es la que asume
    el metodo analitico del paper; los datos la rechazan (sube el chi2). Se
    reporta para tener el numero.

Lo que SI es incertidumbre de la definicion del area es cuanto se mueve D al
mover la ventana de integracion dentro de lo razonable, y eso obliga a releer
los XML (~40 s por variante), asi que va detras de `gates=True`. Sin eso el
presupuesto NO tiene ese termino y el informe lo dice, en vez de omitirlo en
silencio.

Correr `python analisis/sistematicos.py` mide el presupuesto de las corridas que
tengan espectro cacheado. Tarda ~1 min por corrida.
"""

import os
import sys

import numpy as np

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analisis import rutas, espectro, simulacion, calibracion

# alcance de cada termino: "comun" se va en la pendiente de la recta de
# calibracion, "punto" es el que hay que comparar con los residuos del test A
ALCANCE = dict(ventana="punto", forma_exceso="punto", gate="comun",
               area_metodo="informativo", scattering="informativo")

# Variaciones del gate de integracion que se prueban con `gates=True`, como
# (pre, post) en ns. None = la ventana que mide `pulso.py`, que sobre el Co-60
# da (2.2, 60.3).
#
# OJO CON EL SIGNO: la ventana es `[onset - pre, onset + post]`, o sea que `pre`
# son los ns ANTES del onset y va POSITIVO. Poniendolo negativo la integral
# arranca DESPUES del onset y se come la subida del pulso: con pre = -5 el
# barrido daba D = 0.765 contra 1.251, un -39% que parecia un sistematico
# gigantesco y era el signo cambiado. Las variantes tienen que ser cambios
# defendibles de la ventana, no ventanas rotas.
GATES = (None, (2.2, 50.0), (2.2, 75.0), (5.0, 60.3), (1.0, 60.3))
#
# Medido sobre el Co-60 (Lvis), 28-ago-2026:
#
#   (pre, post) [ns]   medida    2.2/50    2.2/75    5.0/60.3   1.0/60.3
#   D                  1.2512    1.2457    1.2554    1.2506     1.2491
#
# Semidispersion 0.39%, y sin deriva: la ventana medida cae en el medio. O sea
# que la eleccion de la ventana no es donde vive la incertidumbre del area.
#
# ⚠ ALCANCE: esto varia la POSICION y el ANCHO de la ventana, con el onset y
# los cortes de calidad fijos. La implementacion anterior barria tambien la
# fraccion del onset (5/10/20%) y el descarte por banderas, y por eso su numero
# era mayor (1.04%). Los dos miden cosas distintas y el de aca es un subconjunto:
# no reemplazarlo sin decir que se amplio el barrido.


def _D(corrida, esp, plantilla, **kw):
    r, _, _, _ = calibracion.calibrar(corrida, esp=esp, plantilla=plantilla,
                                      verbose=False, n_arranques=5, **kw)
    return r


def _semidispersion(valores):
    """(max - min)/2, que es la convencion para una variacion acotada."""
    v = np.asarray(valores, float)
    return float((v.max() - v.min()) / 2.0)


def presupuesto(corrida, rama="Lvis_keVee", metodo="global", esp=None,
                plantilla=None, verbose=True, gates=False):
    """
    Mide el presupuesto de error de D para una corrida. Devuelve un dict.

    `esp` y `plantilla` se pasan ya cargados cuando se los tiene. `gates=True`
    agrega el termino de definicion del area releyendo los XML con la ventana de
    integracion movida: es el unico que necesita los datos crudos, y sin el el
    presupuesto queda INCOMPLETO (el informe lo avisa).
    """
    c = rutas.resolver(corrida)
    e = esp if esp is not None else espectro.cargar(
        "espectro_%s_%s" % (c.clave, metodo))
    pl = plantilla if plantilla is not None else simulacion.cargar(c.sim, rama)
    base = dict(rama=rama, metodo=metodo)

    r0 = _D(c, e, pl, **base)
    D0 = r0["D"]
    terminos, detalle = {}, {}

    # --- ventana ---------------------------------------------------------
    Ds = []
    for ini in (0.25, 0.30, 0.35, 0.45):
        Ds.append(_D(c, e, pl, inicio_ajuste_rel=ini, **base)["D"])
    terminos["ventana"] = _semidispersion(Ds)
    detalle["ventana"] = list(zip((0.25, 0.30, 0.35, 0.45), Ds))

    # --- forma del exceso ------------------------------------------------
    Ds = []
    for tipo in ("exp", "pot"):
        Ds.append(_D(c, e, pl, tipo_exceso=tipo, **base)["D"])
    terminos["forma_exceso"] = _semidispersion(Ds)
    detalle["forma_exceso"] = list(zip(("exp", "pot"), Ds))

    # --- definicion del area ---------------------------------------------
    otro = "convergencia" if metodo == "global" else "global"
    nombre = "espectro_%s_%s" % (c.clave, otro)
    if os.path.exists(rutas.cache("%s.npz" % nombre)):
        r2 = _D(c, espectro.cargar(nombre), pl, rama=rama, metodo=otro)
        # informativo: se descarta por fisica y no por chi2 (ver la cabecera)
        terminos["area_metodo"] = _semidispersion([D0, r2["D"]])
        detalle["area_metodo"] = [(metodo, D0, r0["chi2_ndf"]),
                                  (otro, r2["D"], r2["chi2_ndf"])]

    # --- gate de integracion (releyendo los XML) --------------------------
    if gates:
        Ds, hechos = [], []
        for g in GATES:
            kw = {} if g is None else dict(pre=g[0], post=g[1])
            try:
                eg = espectro.de_corrida(c.clave, metodo=metodo, verbose=False,
                                         **kw)
            except Exception as ex:
                detalle["gate"] = "no se pudo releer los XML: %s" % ex
                break
            Ds.append(_D(c, eg, pl, **base)["D"])
            hechos.append(g)
        else:
            terminos["gate"] = _semidispersion(Ds)
            detalle["gate"] = list(zip(hechos, Ds))
    else:
        detalle["gate"] = "no medido (pide releer los XML: gates=True)"

    # --- mezcla de scattering (informativo) ------------------------------
    try:
        pl1 = simulacion.cargar(c.sim, rama, solo_simple=True)
        r1 = _D(c, e, pl1, **base)
        terminos["scattering"] = _semidispersion([D0, r1["D"]])
        detalle["scattering"] = [("completa", D0, r0["chi2_ndf"]),
                                 ("solo 1 scattering", r1["D"], r1["chi2_ndf"])]
    except Exception as ex:                      # la rama NCompton puede faltar
        detalle["scattering"] = "no se pudo: %s" % ex

    suma = lambda tipo: float(np.sqrt(sum(
        v ** 2 for k, v in terminos.items() if ALCANCE.get(k) == tipo)))
    res = dict(corrida=c.clave, punto=c.punto(), rama=rama, D=D0,
               err_estadistico=r0["err_D"], terminos=terminos, detalle=detalle,
               sist_punto=suma("punto"), sist_comun=suma("comun"))
    res["sist_total"] = float(np.hypot(res["sist_punto"], res["sist_comun"]))
    res["err_total"] = float(np.hypot(res["err_estadistico"], res["sist_total"]))
    if verbose:
        informe(res)
    return res


def informe(r):
    u = "keVee" if r["rama"].startswith("Lvis") else "keV"
    print("-" * 76)
    print("  PRESUPUESTO DE ERROR   %s  (%s)   escala %s"
          % (r["corrida"], r["punto"], u))
    print("-" * 76)
    print("     D = %.4f mV*ns/%s" % (r["D"], u))
    print("     %-22s %10s %9s   %s" % ("termino", "en D", "relativo", "alcance"))
    print("     %-22s %10.4f %8.2f %%   %s"
          % ("estadistica", r["err_estadistico"],
             100 * r["err_estadistico"] / r["D"], "punto a punto"))
    for k, v in sorted(r["terminos"].items()):
        print("     %-22s %10.4f %8.2f %%   %s"
              % (k, v, 100 * v / r["D"], ALCANCE.get(k, "?")))
    print("     %-22s %10s %8s" % ("-" * 22, "-" * 10, "-" * 8))
    print("     %-22s %10.4f %8.2f %%   (lo que se compara con los residuos)"
          % ("sistematico punto", r["sist_punto"], 100 * r["sist_punto"] / r["D"]))
    print("     %-22s %10.4f %8.2f %%   (se va en la pendiente)"
          % ("sistematico comun", r["sist_comun"], 100 * r["sist_comun"] / r["D"]))
    print("     %-22s %10.4f %8.2f %%"
          % ("TOTAL", r["err_total"], 100 * r["err_total"] / r["D"]))
    if isinstance(r["detalle"].get("gate"), str):
        print("     ⚠ INCOMPLETO: falta el termino de definicion del area")
        print("       (%s)" % r["detalle"]["gate"])
    d = r["detalle"].get("area_metodo")
    if isinstance(d, list):
        print("     informativo, NO sumado: el otro metodo de area")
        print("        %-20s D = %.4f  chi2/ndf = %.3f" % d[0])
        print("        %-20s D = %.4f  chi2/ndf = %.3f" % d[1])
        print("        (el chi2 no los distingue; el de convergencia se")
        print("         descarta por fisica: su ventana depende de la amplitud)")
    d = r["detalle"].get("scattering")
    if isinstance(d, list):
        print("     informativo, NO sumado: con plantilla de solo 1 scattering")
        print("        %-20s D = %.4f  chi2/ndf = %.3f" % d[0])
        print("        %-20s D = %.4f  chi2/ndf = %.3f" % d[1])
        print("        (es lo que asume la plantilla analitica del paper; si el")
        print("         chi2 empeora, los datos la estan rechazando)")


def _demo():
    rutas.salida_utf8()
    gates = "--gates" in sys.argv
    if not gates:
        print("  (sin --gates: no se mide el termino de definicion del area,")
        print("   que obliga a releer los XML. El presupuesto queda incompleto.)")
        print()
    hechas = []
    for clave in sorted(rutas.CORRIDAS):
        c = rutas.corrida(clave)
        if not os.path.exists(rutas.cache("espectro_%s_global.npz" % clave)):
            continue
        if not (c.sim and os.path.exists(c.sim)):
            continue
        hechas.append(presupuesto(c, gates=gates))
        print()
    if not hechas:
        print("  no hay corridas con espectro cacheado")
        return
    print("=" * 76)
    print("RESUMEN")
    print("=" * 76)
    print("  %-14s %-14s %8s %9s %9s %9s"
          % ("corrida", "punto", "D", "estad.", "sist.", "total"))
    for r in hechas:
        print("  %-14s %-14s %8.4f %8.2f%% %8.2f%% %8.2f%%"
              % (r["corrida"], r["punto"], r["D"],
                 100 * r["err_estadistico"] / r["D"],
                 100 * r["sist_total"] / r["D"],
                 100 * r["err_total"] / r["D"]))
    print()
    razon = max(r["sist_total"] / max(r["err_estadistico"], 1e-12) for r in hechas)
    if razon > 1.5:
        print("  El sistematico llega a ser %.1f veces el estadistico: lo que"
              % razon)
        print("  limita esta calibracion no es cuantos eventos se midieron.")
    else:
        print("  Estadistico y sistematico son del mismo orden (razon maxima")
        print("  %.1f). Ojo: falta el termino de definicion del area en las" % razon)
        print("  corridas que se corrieron sin --gates.")


if __name__ == "__main__":
    _demo()
