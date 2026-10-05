"""
Que la ventana de `analisis/gui.py` NO recorte lo que dibuja.

Por que existe
--------------
Las figuras salian cortadas por el 20% de la derecha y el 20% de abajo, y no
habia ningun error: TkAgg escala la figura por el device pixel ratio de la
pantalla DESPUES de fijar su tamano en pulgadas, y adentro de un Notebook el
widget no puede crecer para compensar, asi que el buffer de Agg terminaba mas
grande que la PhotoImage sobre la que se vuelca y el blit recortaba en
silencio. Ver el docstring de `gui.Lienzo`.

Un recorte asi no se ve en una figura guardada con `savefig` -que es lo que
prueba `_verif_gui_pestanas.py`-: solo se ve en la ventana. De ahi que haga
falta esta verificacion aparte, que compara los DOS tamanos, el del render y
el de la imagen, en cada pestana y a varios tamanos de ventana.

De paso comprueba las tres cosas que el recorte tapaba o que se arreglaron con
el: que la barra de progreso de abajo este visible aunque la ventana sea
chica, que la columna de controles se pueda desplazar cuando no entra, y que
agrandar un panel (doble clic o selector) muestre un panel solo y vuelva.

Y desde el 14-sep-2026 prueba primero EL ARRANQUE REAL (7.63): la ventana
maximizada de entrada y el `mainloop` de verdad, no `update()` en bucle con la
ventana en `normal`, que era lo unico que se probaba y no es lo que el usuario
ve. Recorre las pestanas, y si hay mas de un monitor pasea la ventana por el
otro y la trae de vuelta (normal y maximizada), midiendo en cada estado que
widget, render e imagen coincidan y que `_CanvasEncuadrado` no haya tenido que
corregir nada (`forzados == 0`): si tuvo que hacerlo, el encuadre por eventos
fallo en ese camino, y la traza dice con que numeros.

    python verificaciones/_verif_gui_encuadre.py
    python verificaciones/_verif_gui_encuadre.py cs137_1600
"""

import os as _os
import sys as _sys

_RAIZ = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _RAIZ not in _sys.path:
    _sys.path.insert(0, _RAIZ)
if _os.path.dirname(_os.path.abspath(__file__)) not in _sys.path:
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))

import time
import types

from analisis import rutas
from analisis import gui as g
import _verif_gui_pestanas as V

TAMANOS = [(1920, 1000), (1550, 830), (1250, 720)]


def asentar(app, veces=3):
    """Deja que Tk termine de acomodar todo (y corra los `after_idle`)."""
    for _ in range(veces):
        app.update_idletasks()
        app.update()
        time.sleep(0.12)
    app.update()


def lienzos(app):
    return list(app.lienzos.items()) + [("resultado", app.lz_resultado)]


def _monitores():
    """[(handle, (l, t, r, b))] de los monitores, en px fisicos (win32)."""
    if _sys.platform != "win32":
        return []
    import ctypes
    from ctypes import wintypes
    u = ctypes.windll.user32
    rects = []
    PROC = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                              ctypes.POINTER(wintypes.RECT), ctypes.c_double)

    def cb(hmon, hdc, prect, data):
        r = prect.contents
        rects.append((hmon, (r.left, r.top, r.right, r.bottom)))
        return 1
    u.EnumDisplayMonitors(None, None, PROC(cb), 0)
    return rects


def _medir(app, nombre, lz, etiqueta, fallas):
    """Un lienzo visible: widget, render e imagen tienen que medir lo mismo."""
    c = lz.canvas
    (wid, ren, foto) = c._medidas()
    bien = (abs(ren[0] - wid[0]) <= 1 and abs(ren[1] - wid[1]) <= 1
            and foto == wid and c.forzados == 0)
    print("     %-30s %-12s widget %s render %s imagen %s forzados %d  %s"
          % (etiqueta, nombre, wid, ren, foto, c.forzados,
             "ok" if bien else "RECORTA"))
    if not bien:
        fallas.append("%s en %s: widget %s, render %s, imagen %s, forzados %d"
                      % (nombre, etiqueta, wid, ren, foto, c.forzados))


def arranque_real(salida):
    """
    La ventana como la ve el usuario: maximizada de entrada y con `mainloop`.

    Todo va encadenado con `after`, que es como llegan los eventos de verdad
    (un `update()` en bucle los entrega en otro orden y con la ventana en otro
    estado). Devuelve la lista de fallas.
    """
    print("\n  arranque real (maximizada, mainloop)")
    fallas = []
    app = g.App()
    pasos = []

    def luego(ms, f):
        app.after(ms, f)

    def siguiente():
        if pasos:
            pasos.pop(0)()

    def p_inicio():
        app.lift()
        app.attributes("-topmost", True)
        app._terminado(salida, salida["corrida"])
        luego(1500, siguiente)

    def p_resultado():
        print("     ventana %s, %s" % (app.geometry(), app.state()))
        _medir(app, "resultado", app.lz_resultado, "arranque", fallas)
        luego(200, siguiente)

    def hacer_pestana(nombre, lz):
        def p():
            app.nb.select(lz)
            luego(700, lambda: (_medir(app, nombre, lz, "pestana", fallas),
                                siguiente()))
        return p

    pasos.append(p_inicio)
    pasos.append(p_resultado)
    for nombre, lz in app.lienzos.items():
        pasos.append(hacer_pestana(nombre, lz))

    # el paseo por monitores, si hay mas de uno
    mons = _monitores()
    if len(mons) > 1:
        import ctypes
        u = ctypes.windll.user32
        lz = app.lienzos["comparacion"]

        def estado(etiqueta, y_despues):
            def p():
                hwnd = int(app.wm_frame(), 16)
                print("     %s: %s, dpi del monitor %d"
                      % (etiqueta, app.geometry(), u.GetDpiForWindow(hwnd)))
                _medir(app, "comparacion", lz, etiqueta, fallas)
                y_despues()
                luego(1200, siguiente)
            return p

        def otro_monitor():
            hwnd = int(app.wm_frame(), 16)
            actual = u.MonitorFromWindow(hwnd, 2)
            for hmon, (l, t, r, b) in mons:
                if hmon != actual:
                    return l + 60, t + 60
            return None

        def a_normal():
            app.state("normal")

        def al_otro():
            xy = otro_monitor()
            if xy:
                app.geometry("+%d+%d" % xy)

        def a_zoomed():
            app.state("zoomed")

        def a_casa():
            app.state("normal")
            l, t = mons[0][1][0], mons[0][1][1]
            app.geometry("+%d+%d" % (l + 60, t + 60))

        pasos.append(lambda: (app.nb.select(lz), luego(700, siguiente)))
        pasos.append(estado("maximizada, monitor de partida", a_normal))
        pasos.append(estado("normal, monitor de partida", al_otro))
        pasos.append(estado("normal, otro monitor", a_zoomed))
        pasos.append(estado("maximizada, otro monitor", a_casa))
        pasos.append(estado("normal, de vuelta", a_zoomed))
        pasos.append(estado("maximizada, de vuelta", lambda: None))
    else:
        pasos.append(lambda: (print("     (un solo monitor: sin paseo)"),
                              luego(100, siguiente)))

    pasos.append(app.destroy)
    luego(1500, siguiente)
    app.mainloop()
    return fallas


def main():
    rutas.salida_utf8()
    clave = _sys.argv[1] if len(_sys.argv) > 1 else "co60"
    print("=" * 74)
    print("encuadre de la ventana de analisis/gui.py sobre %s" % clave)
    print("=" * 74)

    salida, como = V.salida_de(clave)
    print("  %s   ->   D = %.4f" % (como, salida["resultado"]["D"]))
    if "ventana" not in salida:
        # sin el archivo crudo no hay pulsos que promediar: esa pestana
        # muestra su cartel de error a proposito, y deja una traza en
        # .cache/gui_error.log que no es una falla de esta verificacion
        print("  (sin archivo crudo: la pestana de la compuerta va a quedar "
              "en su cartel de error)")

    g.registro.guardar = lambda *a, **k: None      # no tocar el registro real
    fallas = arranque_real(salida)

    print("\n  ventana en normal, a tres tamanos (update en bucle)")
    app = g.App()
    app.state("normal")
    asentar(app)
    app._terminado(salida, salida["corrida"])
    asentar(app)
    for ancho, alto in TAMANOS:
        app.geometry("%dx%d+0+0" % (ancho, alto))
        asentar(app)
        print("\n  ventana %dx%d" % (ancho, alto))

        # (1) la barra de abajo tiene que estar, la ventana mida lo que mida
        barra = [w for w in app.winfo_children()
                 if w.winfo_class() == "TFrame"][-1]
        vista = bool(app.b_cancelar.winfo_ismapped())
        print("     barra de progreso visible: %s" % vista)
        if not vista:
            fallas.append("la barra de abajo no se ve a %dx%d" % (ancho, alto))

        # (2) la columna de controles se desplaza si no entra
        col = app._col_lienzo
        caja = col.bbox("all")
        sb = [w for w in col.master.winfo_children()
              if w.winfo_class() == "TScrollbar"][0]
        falta = bool(caja) and caja[3] > col.winfo_height()
        print("     columna: %d px de alto, contenido %d, scrollbar %s"
              % (col.winfo_height(), caja[3] if caja else 0,
                 "si" if sb.winfo_ismapped() else "no"))
        if falta != bool(sb.winfo_ismapped()):
            fallas.append("la scrollbar de la columna no corresponde a %dx%d"
                          % (ancho, alto))

        # (3) el render y la imagen tienen que medir lo mismo en TODAS
        for nombre, lz in lienzos(app):
            app.nb.select(lz if nombre != "resultado" else 0)
            asentar(app, 1)
            _medir(app, nombre, lz, "%dx%d" % (ancho, alto), fallas)

    # (4) agrandar un panel y volver, por las dos vias
    print("\n  agrandar un panel")
    app.geometry("1550x830+0+0")
    asentar(app)
    for nombre, lz in lienzos(app):
        app.nb.select(lz if nombre != "resultado" else 0)
        asentar(app, 1)
        n = len(lz._paneles)
        if n < 2:
            print("     %-12s %d panel, no aplica" % (nombre, n))
            continue
        visibles = []
        for k in range(n):
            # el doble clic de verdad: se le pasa el punto medio del panel
            x0, y0, x1, y1 = lz._paneles[k].bbox.extents
            lz._al_clic(types.SimpleNamespace(dblclick=True,
                                              x=(x0 + x1) / 2,
                                              y=(y0 + y1) / 2))
            asentar(app, 1)
            visibles.append(sum(1 for a in lz.fig.axes if a.get_visible()))
            if lz._solo != k:
                fallas.append("%s: doble clic sobre el panel %d no lo agrando"
                              % (nombre, k + 1))
            lz._al_clic(types.SimpleNamespace(dblclick=True, x=1, y=1))
            asentar(app, 1)
            if lz._solo is not None:
                fallas.append("%s: el segundo doble clic no volvio" % nombre)
        # y por el selector
        lz.cb_panel.current(1)
        lz._al_elegir_panel()
        asentar(app, 1)
        por_selector = lz._solo == 0
        lz.cb_panel.current(0)
        lz._al_elegir_panel()
        asentar(app, 1)
        print("     %-12s %d paneles, ejes visibles al agrandar %s, "
              "selector %s" % (nombre, n, visibles,
                               "ok" if por_selector else "FALLA"))
        if not por_selector:
            fallas.append("%s: el selector no agranda" % nombre)
        if any(v > 2 for v in visibles):
            fallas.append("%s: agrandar dejo mas de un panel a la vista"
                          % nombre)

    # una captura para mirar, con un panel agrandado
    lz = app.lienzos["comparacion"]
    app.nb.select(lz)
    asentar(app, 1)
    lz.elegir_panel(0)
    asentar(app, 1)
    destino = rutas.figura("gui_encuadre_%s.png" % clave)
    lz.fig.savefig(destino)
    print("\n  captura del panel agrandado -> %s" % rutas.relativa(destino))

    app.destroy()
    if fallas:
        print("\n  %d PROBLEMAS" % len(fallas))
        for f in fallas:
            print("     %s" % f)
        return 1
    print("\n  sin recortes en %d tamanos de ventana" % len(TAMANOS))
    return 0


if __name__ == "__main__":
    _sys.exit(main())
