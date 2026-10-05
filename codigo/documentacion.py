"""
documentacion.py - genera `analisis_paso_a_paso.pdf`: el recorrido completo del
analisis, paso por paso, contra el codigo que de verdad esta corriendo.

Por que existe
--------------
CLAUDE.md es un cuaderno de laboratorio: crece por el final, mezcla el hallazgo
con la correccion del hallazgo, y hay que leerlo entero para saber en que quedo
cada cosa. Este PDF es la otra vista: el ESTADO ACTUAL en orden de ejecucion,
del XML crudo hasta el numero final, diciendo en cada paso que archivo y que
funcion lo hacen.

Coherencia con el codigo
------------------------
Los numeros del documento NO estan tipeados: se calculan al generarlo. Las
constantes (compuerta, umbral de saturacion, cotas) se leen de los modulos, las
figuras salen de los datos reales, y el resultado sale de correr el ajuste. Si
alguien cambia `spectrum.GATE_POST`, el PDF de la proxima corrida lo dice.

Uso
---
    python documentacion.py                 # completo (~1 min)
    python documentacion.py --sin-figuras   # solo texto, para iterar rapido
"""

# --- para poder correr este archivo solo, desde cualquier carpeta ----------
# `python codigo/documentacion.py` pone en sys.path la carpeta DEL ARCHIVO, no
# `codigo/`, asi que los imports por paquete no resolverian. Estas lineas
# agregan `codigo/` una sola vez y son inocuas cuando el modulo se importa.
import os as _os
import sys as _sys
_CODIGO = _os.path.dirname(_os.path.abspath(__file__))
if _CODIGO not in _sys.path:
    _sys.path.insert(0, _CODIGO)
# --------------------------------------------------------------------------


import os
import sys
import time

import numpy as np


import rutas

# El PDF va a la RAIZ del proyecto (es el entregable, tiene que estar a
# mano); las figuras intermedias, al cache.
SALIDA = os.path.join(rutas.RAIZ, "analisis_paso_a_paso.pdf")
CARPETA_FIG = os.path.join(rutas.CACHE, "doc")


# ======================================================================
#  1) tipografia
# ======================================================================
def registrar_fuentes():
    """
    Usa las DejaVu que trae matplotlib, para tener griegas y flechas.

    Las Type1 estandar de reportlab (Helvetica y companiaa) estan en WinAnsi:
    no tienen sigma, alfa ni chi, que aparecen en casi todas las formulas de
    este proyecto. matplotlib instala DejaVu, que es Unicode completa, asi que
    se registran de ahi en vez de pedir una instalacion aparte.

    Devuelve (normal, negrita, italica, mono). Si algo falla, cae a Helvetica.
    """
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import matplotlib
    base = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
    pedido = [("DJV", "DejaVuSans.ttf"), ("DJV-B", "DejaVuSans-Bold.ttf"),
              ("DJV-I", "DejaVuSans-Oblique.ttf"), ("DJV-M", "DejaVuSansMono.ttf")]
    try:
        for nombre, archivo in pedido:
            pdfmetrics.registerFont(TTFont(nombre, os.path.join(base, archivo)))
        pdfmetrics.registerFontFamily("DJV", normal="DJV", bold="DJV-B",
                                      italic="DJV-I", boldItalic="DJV-B")
        return "DJV", "DJV-B", "DJV-I", "DJV-M"
    except Exception as e:
        print("  [aviso] sin DejaVu (%s): las griegas van a salir mal" % e)
        return "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Courier"


def estilos(fuentes):
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_JUSTIFY
    from reportlab.lib import colors
    normal, negrita, italica, mono = fuentes
    E = {}
    E["titulo"] = ParagraphStyle("titulo", fontName=negrita, fontSize=21,
                                 leading=26, spaceAfter=6)
    E["subtitulo"] = ParagraphStyle("subtitulo", fontName=normal, fontSize=11.5,
                                    leading=15, textColor=colors.HexColor("#444444"),
                                    spaceAfter=14)
    E["h1"] = ParagraphStyle("h1", fontName=negrita, fontSize=14.5, leading=18,
                             spaceBefore=16, spaceAfter=6,
                             textColor=colors.HexColor("#12395e"))
    E["h2"] = ParagraphStyle("h2", fontName=negrita, fontSize=11.5, leading=14,
                             spaceBefore=10, spaceAfter=3,
                             textColor=colors.HexColor("#1f6091"))
    E["cuerpo"] = ParagraphStyle("cuerpo", fontName=normal, fontSize=9.6,
                                 leading=13.4, alignment=TA_JUSTIFY,
                                 spaceAfter=5)
    E["lista"] = ParagraphStyle("lista", parent=E["cuerpo"], leftIndent=13,
                                bulletIndent=3, spaceAfter=3)
    E["codigo"] = ParagraphStyle("codigo", fontName=mono, fontSize=8.4,
                                 leading=11.4, leftIndent=10,
                                 textColor=colors.HexColor("#333333"),
                                 spaceBefore=3, spaceAfter=6)
    E["nota"] = ParagraphStyle("nota", parent=E["cuerpo"], fontSize=9.0,
                               leading=12.4, leftIndent=8, rightIndent=8,
                               spaceBefore=4, spaceAfter=4)
    E["pie"] = ParagraphStyle("pie", fontName=italica, fontSize=8.4,
                              leading=11, textColor=colors.HexColor("#555555"),
                              spaceBefore=2, spaceAfter=10)
    E["tabla"] = ParagraphStyle("tabla", fontName=normal, fontSize=8.6,
                                leading=11)
    E["tabla_b"] = ParagraphStyle("tabla_b", fontName=negrita, fontSize=8.6,
                                  leading=11)
    return E


# ======================================================================
#  2) helpers de contenido
# ======================================================================
def esc(t):
    """Escapa lo que reportlab interpreta como marcado."""
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


class Doc:
    """Acumulador de flowables con atajos legibles."""

    def __init__(self, E, fuentes):
        from reportlab.platypus import Spacer
        self.E = E
        self.f = fuentes
        self.items = []
        self._Spacer = Spacer

    # --- bloques ---
    def h1(self, texto, ancla=None):
        from reportlab.platypus import Paragraph
        t = esc(texto)
        if ancla:
            t = '<a name="%s"/>%s' % (ancla, t)
        self.items.append(Paragraph(t, self.E["h1"]))

    def h2(self, texto):
        from reportlab.platypus import Paragraph
        self.items.append(Paragraph(esc(texto), self.E["h2"]))

    def p(self, texto, estilo="cuerpo"):
        from reportlab.platypus import Paragraph
        self.items.append(Paragraph(marcar(texto), self.E[estilo]))

    def li(self, texto, vinieta="•"):
        from reportlab.platypus import Paragraph
        self.items.append(Paragraph(marcar(texto), self.E["lista"],
                                    bulletText=vinieta))

    def code(self, texto):
        from reportlab.platypus import Paragraph
        html = esc(texto).replace("\n", "<br/>").replace(" ", "&nbsp;")
        self.items.append(Paragraph(html, self.E["codigo"]))

    def nota(self, titulo, texto, color="#fff7e0", borde="#e0b040"):
        """Recuadro para las trampas y decisiones de diseño."""
        from reportlab.platypus import Paragraph, Table, TableStyle
        from reportlab.lib import colors
        from reportlab.lib.units import mm
        cuerpo = Paragraph("<b>%s</b> %s" % (esc(titulo), marcar(texto)),
                           self.E["nota"])
        t = Table([[cuerpo]], colWidths=[165 * mm])
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(color)),
            ("BOX", (0, 0), (-1, -1), 0.7, colors.HexColor(borde)),
            ("LEFTPADDING", (0, 0), (-1, -1), 7),
            ("RIGHTPADDING", (0, 0), (-1, -1), 7),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
        self.items.append(t)
        self.esp(4)

    def mejorar(self, texto):
        self.nota("Para mejorar:", texto, color="#eef6ec", borde="#79a86b")

    def tabla(self, filas, anchos=None, cabecera=True):
        from reportlab.platypus import Paragraph, Table, TableStyle
        from reportlab.lib import colors
        from reportlab.lib.units import mm
        datos = []
        for i, fila in enumerate(filas):
            est = self.E["tabla_b"] if (cabecera and i == 0) else self.E["tabla"]
            datos.append([Paragraph(marcar(str(c)), est) for c in fila])
        anchos = [a * mm for a in anchos] if anchos else None
        t = Table(datos, colWidths=anchos, repeatRows=1 if cabecera else 0)
        estilo = [("VALIGN", (0, 0), (-1, -1), "TOP"),
                  ("LINEBELOW", (0, 0), (-1, -1), 0.25,
                   colors.HexColor("#cccccc")),
                  ("LEFTPADDING", (0, 0), (-1, -1), 4),
                  ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                  ("TOPPADDING", (0, 0), (-1, -1), 3),
                  ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
        if cabecera:
            estilo += [("BACKGROUND", (0, 0), (-1, 0),
                        colors.HexColor("#eef2f6")),
                       ("LINEBELOW", (0, 0), (-1, 0), 0.7,
                        colors.HexColor("#8899aa"))]
        t.setStyle(TableStyle(estilo))
        self.items.append(t)
        self.esp(6)

    def figura(self, ruta, pie, ancho_mm=165):
        from reportlab.platypus import Image, Paragraph, KeepTogether
        from reportlab.lib.units import mm
        if not ruta or not os.path.exists(ruta):
            self.p("[figura no generada: se corrió con --sin-figuras]", "pie")
            return
        from PIL import Image as PILImage
        try:
            w, h = PILImage.open(ruta).size
        except Exception:
            w, h = 1400, 600
        ancho = ancho_mm * mm
        img = Image(ruta, width=ancho, height=ancho * h / w)
        cap = Paragraph(marcar(pie), self.E["pie"])
        self.items.append(KeepTogether([img, cap]))

    def esp(self, alto=6):
        self.items.append(self._Spacer(1, alto))

    def salto(self):
        from reportlab.platypus import PageBreak
        self.items.append(PageBreak())


def marcar(texto):
    """
    Marcado minimo: **negrita**, `codigo`, y el resto escapado.

    Se escapa PRIMERO y se sustituye despues, para que un `<Data>` del XML
    aparezca como texto y no como una etiqueta rota.
    """
    t = esc(texto)
    partes = t.split("`")
    t = "".join(p if i % 2 == 0
                else '<font face="DJV-M" size="8.6">%s</font>' % p
                for i, p in enumerate(partes))
    partes = t.split("**")
    t = "".join(p if i % 2 == 0 else "<b>%s</b>" % p
                for i, p in enumerate(partes))
    return t


# ======================================================================
#  3) los datos y los numeros del documento (se CALCULAN, no se tipean)
# ======================================================================
def reunir_numeros(rapido=False, avisar=print):
    """
    Corre la cadena sobre los datos de referencia y devuelve lo que el texto
    necesita citar. Todo lo que el PDF afirma numericamente sale de aca.
    """
    from datos_io import spectrum
    from etapas import physics
    import gui

    d = gui.Datos()
    n = {"fuente": d.fuente, "archivos": (d.archivo_muestra, d.archivo_fondo,
                                          d.archivo_sim)}
    avisar("  espectro neto ...")
    r = d.espectro()
    n["r"] = r
    n["borde"] = d.borde_estimado()
    n["iv"] = d.ventana_ajuste()
    n["constantes"] = {
        "GATE_PRE": spectrum.GATE_PRE, "GATE_POST": spectrum.GATE_POST,
        "GUARDA_BASE": spectrum.GUARDA_BASE, "FRAC_ONSET": spectrum.FRAC_ONSET,
        "V_SATURA": spectrum.V_SATURA,
    }
    n["hay_sim"] = d.hay_simulacion()
    n["compat"] = d.compat_simulacion()
    if n["hay_sim"] and n["compat"][0] and not rapido:
        avisar("  modelo extendido (Lvis) ...")
        n["fit"] = gui.ajuste_extendido(d, r, n["iv"], rama="Lvis_keVee")
        avisar("  modelo extendido (Edep) ...")
        n["fit_edep"] = gui.ajuste_extendido(d, r, n["iv"], rama="Edep_keV")
        n["lineas_sim"] = d.lineas_sim()
    n["datos"] = d
    return n


# ======================================================================
#  4) figuras, generadas de los datos reales
# ======================================================================
def _guardar(fig, nombre):
    os.makedirs(CARPETA_FIG, exist_ok=True)
    ruta = os.path.join(CARPETA_FIG, nombre + ".png")
    fig.savefig(ruta, dpi=155, bbox_inches="tight")
    return ruta


def figuras(n, rapido=False, avisar=print):
    """Todas las figuras del documento. Devuelve {nombre: ruta}."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5,
                         "axes.labelsize": 8.5, "legend.fontsize": 7.6,
                         "figure.autolayout": True})
    from datos_io import spectrum
    from etapas.physics import dsigma_dE, compton_edge_energy
    from etapas.model import forward_model, punto_mitad_altura

    F = {}
    fu, r = n["fuente"], n["r"]
    can, neto, err = r["canales"], r["neto"], r["error"]

    # ---------- 1) el pulso crudo y las dos lineas base ----------
    avisar("  figura: forma del pulso ...")
    N = 1500 if rapido else 4000
    P = np.zeros(1024)
    s_corta, s_larga, t_pico = [], [], []
    n_forma, dt = 0, 0.198
    for k, (_ts, pares, _sc) in enumerate(
            spectrum._iter_eventos(n["datos"].archivo_muestra)):
        if k >= N:
            break
        # `pares` llega como lista de cadenas "t,v" del XML y como array
        # (n, 2) del binario del DRS4; `pares_a_pulso` unifica los dos.
        t, v = spectrum.pares_a_pulso(pares)
        if t is None or t.size != 1024:
            continue
        tau = t - t[0]
        dt = float(np.median(np.diff(tau)))
        n_base = max(int(round(12.0 / dt)), 5)
        b_corta = np.median(v[:n_base])
        b_larga = v[:int(len(v) * 0.2)].mean()
        vc = v - b_corta
        i = int(np.argmin(vc))
        a = -vc.min()
        t_pico.append(tau[i])
        if a < 20 or not (8 < tau[i] < tau[-1] - 60):
            continue
        s_corta.append(np.trapezoid(vc, x=tau) / a)
        s_larga.append(np.trapezoid(v - b_larga, x=tau) / a)
        P += vc / a
        n_forma += 1
    P /= max(n_forma, 1)
    tt = np.arange(1024) * dt

    # Las ventanas se dibujan con las CONSTANTES DEL MODULO y con el onset
    # medido sobre el pulso medio, no con numeros a ojo: si alguien cambia
    # spectrum.GATE_POST, la figura del proximo PDF lo muestra.
    C = n["constantes"]
    i_pico = int(np.argmin(P))
    idx_on = np.flatnonzero(np.abs(P[:i_pico + 1]) < C["FRAC_ONSET"] * abs(P.min()))
    t_on = tt[int(idx_on[-1])] if idx_on.size else tt[i_pico]

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
    ax[0].plot(tt, P, color="C0", lw=1.4, label="pulso medio (normalizado)")
    ax[0].axvspan(0, max(t_on - C["GATE_PRE"] - C["GUARDA_BASE"], 0.1),
                  color="C2", alpha=0.22,
                  label="línea base nueva (pre-onset, con %.0f ns de guarda)"
                        % C["GUARDA_BASE"])
    ax[0].axvspan(0, 1024 * dt * 0.2, color="C3", alpha=0.12,
                  label="línea base vieja (20% del registro)")
    ax[0].axvspan(t_on - C["GATE_PRE"], t_on + C["GATE_POST"], color="C1",
                  alpha=0.13,
                  label="compuerta [onset−%.0f, onset+%.0f] ns"
                        % (C["GATE_PRE"], C["GATE_POST"]))
    ax[0].axvline(t_on, color="C1", lw=1.0, ls="--",
                  label="onset (%.0f%% de la amplitud)" % (100 * C["FRAC_ONSET"]))
    ax[0].axhline(0, color="0.6", lw=0.7)
    ax[0].set_xlabel("tiempo relativo al registro [ns]")
    ax[0].set_ylabel("amplitud normalizada")
    ax[0].set_title("Forma media del pulso (%d eventos)" % n_forma)
    ax[0].legend(loc="lower right")
    ax[0].grid(alpha=0.3)

    sc, sl = np.array(s_corta), np.array(s_larga)
    lim = np.percentile(np.abs(np.concatenate([sc, sl])), 99)
    bb = np.linspace(-lim, lim, 90)
    # ancho ROBUSTO (1.4826*MAD): la desviacion estandar de estas dos
    # distribuciones esta dominada por las colas, y lo que interesa comparar es
    # el nucleo, que es lo que se traduce en resolucion del espectro.
    mad = lambda x: 1.4826 * float(np.median(np.abs(x - np.median(x))))
    ax[1].hist(sl, bins=bb, color="C3", alpha=0.55,
               label="línea base vieja: ancho %.2f ns, %d%% con el signo dado "
                     "vuelta" % (mad(sl), round(100 * np.mean(sl > 0))))
    ax[1].hist(sc, bins=bb, color="C2", alpha=0.65,
               label="línea base nueva: ancho %.2f ns" % mad(sc))
    ax[1].axvline(0, color="k", lw=0.8)
    ax[1].set_xlabel("integral / amplitud  [ns]")
    ax[1].set_ylabel("eventos")
    ax[1].set_title("El mismo observable, medido de las dos formas")
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    F["pulso"] = _guardar(fig, "pulso")
    plt.close(fig)

    # ---------- 2) el espectro neto ----------
    avisar("  figura: espectro neto ...")
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.step(can, r["muestra"], where="mid", color="C0", lw=1.0, alpha=0.75,
            label="muestra (%s), escalada a tiempo vivo" % fu.nombre)
    ax.step(can, r["fondo"], where="mid", color="0.45", lw=1.0, alpha=0.8,
            label="fondo, escalado al mismo tiempo vivo")
    ax.step(can, neto, where="mid", color="C3", lw=1.3, label="neto = resta")
    ax.axvspan(n["iv"][0], n["iv"][1], color="C2", alpha=0.10,
               label="ventana de ajuste [%.0f, %.0f]" % n["iv"])
    for E in fu.bordes:
        pass
    ax.set_xlabel("área del pulso [mV·ns]")
    ax.set_ylabel("cuentas normalizadas")
    ax.set_title("Etapa 0 — espectro neto sobre binning común")
    ax.legend()
    ax.grid(alpha=0.3)
    F["espectro"] = _guardar(fig, "espectro")
    plt.close(fig)

    # ---------- 3) la lineshape teorica y el forward model ----------
    avisar("  figura: Klein-Nishina y forward model ...")
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
    E = np.linspace(1, 1.25 * max(fu.bordes), 3000)
    for Eg, c in zip(fu.lineas, ("C0", "C2", "C4")):
        y = dsigma_dE(E, Eg)
        Em = compton_edge_energy(Eg)
        ax[0].plot(E, y, color=c, lw=1.5,
                   label="E$_\\gamma$ = %.1f keV → borde %.1f" % (Eg, Em))
        ax[0].axvline(Em, color=c, ls="--", alpha=0.6)
    ax[0].set_xlabel("energía del electrón de retroceso [keV]")
    ax[0].set_ylabel("dσ/dE$_e$ (forma, u.a.)")
    ax[0].set_title("Etapa 1 — lineshape de Klein-Nishina (Ec. 10)")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    B = float(n["borde"]["edge_ch"])
    ejec = np.linspace(0, 1.6 * B, 1200)
    for sg, c in ((0.02 * B, "C0"), (0.09 * B, "C1"), (0.18 * B, "C3")):
        y = forward_model(ejec, fu.lineas[0], 1.0, B, sg)
        mitad = punto_mitad_altura(ejec, y)
        ax[1].plot(ejec, y / y.max(), color=c, lw=1.4,
                   label="σ/B = %.2f → 50%% en %.0f (%+.1f%%)"
                         % (sg / B, mitad, 100 * (mitad / B - 1)))
        ax[1].axvline(mitad, color=c, ls=":", lw=1.0)
    ax[1].axvline(B, color="k", lw=1.1, label="borde verdadero B = %.0f" % B)
    ax[1].set_xlabel("área del pulso [mV·ns]")
    ax[1].set_ylabel("modelo (normalizado)")
    ax[1].set_title("Etapa 2 — el método clásico del 50% está sesgado")
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    F["forward"] = _guardar(fig, "forward")
    plt.close(fig)

    # ---------- 4) el NCC ----------
    avisar("  figura: NCC ...")
    from etapas.ncc import plantilla_referencia
    suave = spectrum.suavizar(neto, ventana=21, orden=3)
    b = n["borde"]
    lo, hi = b["lo"], b["hi"]
    m = (can >= lo) & (can <= hi)
    # La plantilla se reconstruye con los parametros GANADORES del barrido, no
    # con unos cualesquiera: si no, la figura mostraria una comparacion que el
    # NCC nunca hizo.
    B_ncc = float(b["edge_ch"])
    s2 = plantilla_referencia(fu.lineas[0],
                              sigma_rel=float(b["sigma"]) / B_ncc,
                              n_puntos=int(m.sum()),
                              inicio_rel=lo / B_ncc, fin_rel=hi / B_ncc)[1]
    s1 = suave[m]
    z = lambda v: (v - v.mean()) / max(v.std(), 1e-12)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
    ax[0].step(can, neto, where="mid", color="0.75", lw=0.9, label="neto")
    ax[0].plot(can, suave, color="C0", lw=1.2, label="suavizado (Savitzky-Golay)")
    ax[0].axvspan(lo, hi, color="C2", alpha=0.14,
                  label="ventana elegida por NCC\n[%.0f, %.0f], ρ = %.4f"
                        % (lo, hi, b["rho"]))
    ax[0].axvline(b["edge_ch"], color="C3", ls="--",
                  label="B = %.0f mV·ns" % b["edge_ch"])
    ax[0].set_xlabel("área del pulso [mV·ns]")
    ax[0].set_ylabel("cuentas netas")
    ax[0].set_title("Etapa 4 — selección automática del intervalo")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    from etapas.ncc import pearson
    ax[1].plot(np.linspace(0, 1, len(s1)), z(s1), color="C0", lw=1.4,
               label="datos (estandarizados)")
    ax[1].plot(np.linspace(0, 1, len(s2)), z(s2), color="C3", lw=1.2, ls="--",
               label="plantilla S2 (σ/B = %.2f), ρ = %.4f"
                     % (float(b["sigma"]) / B_ncc, pearson(s1, s2)))
    ax[1].set_xlabel("posición dentro de la ventana")
    ax[1].set_ylabel("z-score")
    ax[1].set_title("Lo que compara Pearson (Ec. 15): la FORMA")
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    F["ncc"] = _guardar(fig, "ncc")
    plt.close(fig)

    # ---------- 5) las plantillas de Geant4 y el ajuste ----------
    if "fit" in n:
        avisar("  figura: plantillas y ajuste ...")
        from ajuste import modelo_ext as me
        f = n["fit"]
        centros, h1, hm, w_nat = f["comp"]
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
        ax[0].step(centros, h1 / h1.max(), where="mid", color="C3", lw=1.2,
                   label="1 scattering (Klein-Nishina)")
        ax[0].step(centros, hm / hm.max(), where="mid", color="C0", lw=1.2,
                   label="scattering múltiple")
        for E in fu.bordes:
            ax[0].axvline(E, color="0.4", ls="--", lw=0.9)
        for E in fu.lineas:
            ax[0].axvline(E, color="0.75", ls=":", lw=0.9)
        ax[0].set_xlabel("energía depositada / luz [keV]")
        ax[0].set_ylabel("normalizadas por separado")
        ax[0].set_title("Etapa 7 — las dos componentes de la simulación\n"
                        "(w = %.3f de scattering simple)" % w_nat)
        ax[0].legend()
        ax[0].grid(alpha=0.3)

        iv = f["intervalo"]
        zz = (can >= 0.25 * iv[0]) & (can <= min(1.15 * iv[1], can[-1]))
        ax[1].errorbar(can[zz], neto[zz], yerr=err[zz], fmt=".", ms=2, lw=0.4,
                       color="0.65", label="datos (neto)")
        ax[1].plot(can[zz], me.modelo_ext(can[zz], centros, h1, hm, f["p"]),
                   color="C0", lw=1.6, label="modelo: D = %.4f" % f["D"])
        p = f["p"].copy()
        ax[1].plot(can[zz], p[me.I_C] * np.exp(-can[zz] / p[me.I_LAM]), "k--",
                   lw=1.0, label="exceso a bajo canal")
        q = p.copy()
        q[me.I_C] = 0.0
        ax[1].plot(can[zz], me.modelo_ext(can[zz], centros, h1, hm, q),
                   color="C3", ls=":", lw=1.2, label="sólo Compton")
        ax[1].axvspan(iv[0], iv[1], color="C2", alpha=0.09)
        ax[1].set_xlabel("área del pulso [mV·ns]")
        ax[1].set_ylabel("cuentas netas")
        ax[1].set_title("Etapa 9 — modelo extendido, χ²/ndf = %.3f (%d gdl)"
                        % (f["chi2_red"], f["ndf"]))
        ax[1].legend()
        ax[1].grid(alpha=0.3)
        F["ajuste"] = _guardar(fig, "ajuste")
        plt.close(fig)

        # ---------- 6) residuos y recta ----------
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
        mm_ = (can >= iv[0]) & (can <= iv[1])
        mod = me.modelo_ext(can[mm_], centros, h1, hm, f["p"])
        ax[0].plot(can[mm_], (neto[mm_] - mod) / np.maximum(err[mm_], 1), ".",
                   ms=3, color="C0")
        ax[0].axhline(0, color="k", lw=0.9)
        ax[0].axhspan(-1, 1, color="C2", alpha=0.14)
        ax[0].set_ylim(-6, 6)
        ax[0].set_xlabel("área del pulso [mV·ns]")
        ax[0].set_ylabel("(dato − modelo)/σ")
        ax[0].set_title("Residuos normalizados")
        ax[0].grid(alpha=0.3)

        Ee = np.linspace(0, 1.25 * max(fu.bordes), 200)
        ax[1].plot(Ee, f["D"] * Ee, color="C0", lw=1.8,
                   label="luz: D = %.4f mV·ns/keVee" % f["D"])
        if "fit_edep" in n:
            ax[1].plot(Ee, n["fit_edep"]["D"] * Ee, color="C2", lw=1.5, ls="--",
                       label="edep: D = %.4f mV·ns/keV" % n["fit_edep"]["D"])
        for E in fu.bordes:
            ax[1].plot(E, f["D"] * E, "k*", ms=11)
            ax[1].annotate("borde %.0f keV → %.0f mV·ns" % (E, f["D"] * E),
                           (E, f["D"] * E), textcoords="offset points",
                           xytext=(-105, 7), fontsize=7)
        ax[1].set_xlabel("energía [keV]")
        ax[1].set_ylabel("área del pulso [mV·ns]")
        ax[1].set_title("La calibración: recta por el origen")
        ax[1].legend()
        ax[1].grid(alpha=0.3)
        F["recta"] = _guardar(fig, "recta")
        plt.close(fig)
    return F


# ======================================================================
#  5) el documento
# ======================================================================
def construir(n, F, E, fuentes):
    from datos_io import spectrum
    from etapas import physics
    d = Doc(E, fuentes)
    fu, r = n["fuente"], n["r"]
    C = n["constantes"]
    b = n["borde"]
    iv = n["iv"]
    fit = n.get("fit")
    fit_e = n.get("fit_edep")
    hoy = time.strftime("%d-%m-%Y")

    # ---------------------------------------------------------------- portada
    from reportlab.platypus import Paragraph
    d.items.append(Paragraph("Calibración de un centellador plástico por "
                             "Compton edge", E["titulo"]))
    d.items.append(Paragraph(
        "El análisis implementado, paso a paso, contra el código que corre.<br/>"
        "Método NCC + recocido simulado (Kong et al., IEEE TNS 2025, "
        "arXiv:2507.16282), adaptado a un centellador plástico leído con SiPM.<br/>"
        "Generado automáticamente el %s por <font face=\"DJV-M\">codigo/documentacion.py</font>."
        % hoy, E["subtitulo"]))

    d.p("Este documento describe **el estado actual** del análisis en orden de "
        "ejecución: qué le pasa a los datos desde el XML del osciloscopio hasta "
        "el número final, qué archivo y qué función hace cada cosa, qué "
        "decisiones de diseño hay detrás y qué quedaría por mejorar en cada "
        "paso. No es el cuaderno de laboratorio (eso es `CLAUDE.md`, que crece "
        "por el final y mezcla el hallazgo con su corrección posterior): acá "
        "está lo que hoy hace el código.")
    d.p("**Los números no están tipeados.** Se calculan al generar el PDF: las "
        "constantes se leen de los módulos, las figuras salen de los datos "
        "reales y el resultado sale de correr el ajuste. Si mañana cambia "
        "`spectrum.GATE_POST`, el documento de la próxima corrida lo dice.")

    d.h2("Sobre qué datos se generó este documento")
    d.tabla([
        ["", "archivo", "qué es"],
        ["Muestra", os.path.basename(n["archivos"][0]),
         "corrida con la fuente"],
        ["Fondo", os.path.basename(n["archivos"][1]),
         "corrida sin fuente, se resta"],
        ["Simulación", os.path.basename(n["archivos"][2]) or "(ninguna)",
         "plantilla de Geant4"],
        ["Fuente", fu.nombre,
         "líneas gamma: " + " y ".join("%.1f" % E for E in fu.lineas) + " keV"],
    ], anchos=[22, 78, 65])

    d.h2("El resultado, para tenerlo a mano")
    if fit:
        filas = [["cantidad", "valor", "de dónde sale"],
                 ["D (luz)", "%.4f ± %.4f mV·ns/keVee"
                  % (fit["D"], fit["err"].get("D", float("nan"))),
                  "`modelo_ext` con plantilla `Lvis_keVee`"],
                 ["D (edep)", "%.4f ± %.4f mV·ns/keV"
                  % (fit_e["D"], fit_e["err"].get("D", float("nan"))),
                  "lo mismo con `Edep_keV`"],
                 ["resolución", "σ_E = %.3f·√E keV → %.1f%% FWHM a %.0f keV"
                  % (fit["alpha"], 100 * fit["fwhm_rel_ref"], fit["E_ref"]),
                  "el parámetro α del mismo ajuste"],
                 ["mezcla w", "%.3f ± %.3f (la simulación trae %.3f)"
                  % (fit["w"], fit["err"].get("w", float("nan")), fit["w_nat"]),
                  "fracción de scattering simple, ajustada"],
                 ["bondad", "χ²/ndf = %.3f con %d grados de libertad"
                  % (fit["chi2_red"], fit["ndf"]), "χ² de Pearson iterado"]]
        d.tabla(filas, anchos=[24, 63, 78])
        d.p("El eje horizontal de todo el proyecto es el **área del pulso en "
            "mV·ns**, no un canal de ADC. `D` es la pendiente que lo convierte "
            "a energía, y el mapa pasa por el origen: área cero ↔ energía cero.")
    else:
        d.p("(Este documento se generó sin la simulación o con "
            "`--sin-figuras`: faltan los números del modelo extendido.)", "pie")

    d.salto()

    # ============================================================ panorama
    d.h1("0. El panorama en una página")
    d.p("El problema: un centellador plástico es de bajo Z, así que domina el "
        "efecto Compton y **no hay fotopico**. La calibración clásica —poner el "
        "pico de energía conocida en su canal— no aplica. Lo que sí hay es el "
        "**borde Compton**: la energía máxima que un gamma le puede dar a un "
        "electrón, que ocurre en retrodispersión y vale")
    d.code("E_max = 2·α·E_γ / (1 + 2·α),      α = E_γ / m_e c²   (m_e c² = %.3f keV)"
           % physics.MEC2)
    d.p("Esa energía es conocida sin ajustar nada. Si uno logra decir en qué "
        "área de pulso cae ese borde, tiene un punto de calibración. Todo el "
        "resto del análisis es responder esa pregunta con cuidado, porque el "
        "borde no es un pico: es un escalón borroneado por la resolución del "
        "detector y contaminado por física que la fórmula no contiene.")

    d.h2("La cadena completa")
    d.tabla([
        ["#", "paso", "archivo", "qué produce"],
        ["1", "XML → área del pulso", "`datos_io/spectrum.py`",
         "una área en mV·ns por evento"],
        ["2", "áreas → espectro neto", "`spectrum.espectro_neto`",
         "histograma con fondo restado y barras de error"],
        ["3", "física teórica", "`etapas/physics.py`",
         "E_max y la forma dσ/dE_e"],
        ["4", "forma → espacio de canal", "`etapas/model.py`",
         "el modelo que se compara con los datos"],
        ["5", "elegir la ventana", "`etapas/ncc.py`",
         "el intervalo de ajuste y la semilla, sin criterio humano"],
        ["6", "ajustar", "`etapas/fit_single.py`, `etapas/fit_double.py`",
         "D por la ruta analítica"],
        ["7", "optimizador global", "`etapas/sa.py`",
         "verificar que el óptimo es global"],
        ["8", "plantilla de Geant4", "`ajuste/sim_template.py`",
         "reemplaza la teoría analítica por la simulación"],
        ["9", "modelo extendido", "`ajuste/modelo_ext.py`",
         "**el resultado vigente**: sistemáticos como parámetros"],
        ["10", "posterior", "`ajuste/mcmc.py`",
         "barras de error y correlaciones reales"],
        ["11", "resultado", "`calibracion.py`, pestaña ★",
         "D, su presupuesto de error y el registro por fuente"],
    ], anchos=[8, 42, 45, 70])
    d.p("Los pasos 3 a 7 son la réplica del paper. Los pasos 8 a 11 son lo que "
        "hubo que agregarle para que el número fuera creíble en **este** "
        "detector, y son también donde está casi todo el error sistemático.")

    d.h2("Dónde vive cada cosa")
    d.code("Calibración CE/\n"
           "├─ analisis_paso_a_paso.pdf   este documento\n"
           "├─ CLAUDE.md                  el cuaderno de laboratorio\n"
           "├─ codigo/\n"
           "│   ├─ rutas.py               dónde está cada cosa\n"
           "│   ├─ gui.py  calibracion.py  documentacion.py\n"
           "│   ├─ etapas/    physics · model · fit_single · fit_double · ncc · sa\n"
           "│   ├─ datos_io/  spectrum\n"
           "│   ├─ ajuste/    sim_template · modelo_ext · mcmc\n"
           "│   └─ scripts/   calibrar_v2 · calibrar_ncc · calibrar\n"
           "├─ datos/       los .xml del DRSOSC, sus cachés .npz y los .root\n"
           "├─ figuras/     lo que los módulos dibujan al correr solos\n"
           "├─ docs/        el paper de referencia y notas históricas\n"
           "├─ verificaciones/  scripts sueltos de comprobación\n"
           "└─ .cache/      preferencias de la GUI, cadenas MCMC, registro")
    d.p("Cada módulo **corre solo** y deja su figura en `figuras/`: "
        "`python codigo/etapas/ncc.py` muestra la Etapa 5 sobre datos "
        "sintéticos con verdad conocida, y así con los demás. Es la forma "
        "prevista de tocar una etapa sin arrastrar el resto. Todas las rutas "
        "se resuelven desde `__file__` vía `codigo/rutas.py`, así que no "
        "importa desde qué carpeta se ejecute.")

    d.esp(4)

    # ============================================================ paso 1
    d.h1("1. Del XML crudo al área del pulso")
    d.p("**Archivo:** `datos_io/spectrum.py`. **Funciones:** `detectar_canal`, "
        "`_iter_eventos`, `_area_pulso_gate`, `tiempo_vivo`, "
        "`leer_areas_y_tiempo`.")
    d.p("El DRSOSC guarda un XML con un tag por línea y la estructura "
        "`DRSOSC > Event > Board_XXXX > CHN{n} > Data`, donde cada `Data` es "
        "un par «t,v» con t en ns y v en mV. El archivo de la muestra pesa "
        "1.76 GB y tiene del orden de 5·10⁷ tags `Data`.")

    d.h2("1.1 Leer el archivo sin morir en el intento")
    d.p("`ElementTree` construye un objeto por tag: decenas de GB de RAM y "
        "varios minutos. Como el formato del DRSOSC es rígido, "
        "`_iter_eventos` lo recorre **como texto**, línea por línea, y "
        "convierte cada pulso de una sola vez con "
        "`np.fromstring(\",\".join(pares), sep=\",\")`. La RAM queda plana y el "
        "archivo grande se lee en unos 30 s.")
    d.p("Dos cosas se detectan en vez de asumirse. `detectar_canal` mira el "
        "primer evento y devuelve qué `CHN{n}` y qué `Scaler{n-1}` trae: las "
        "corridas de junio salieron por `CHN2` y las de agosto por `CHN1`, y "
        "tener el canal cableado hacía que se recorrieran los 2.9 GB "
        "encontrando **cero** pulsos, en silencio. Además el resultado se "
        "cachea en `<archivo>.areas_gate.npz` con una huella que incluye "
        "tamaño, fecha, umbral y canal, así que la segunda corrida es "
        "instantánea.")
    d.nota("Trampa que ya costó cara:",
           "el caché guardaba también los resultados vacíos. Hoy "
           "`leer_areas_y_tiempo` levanta `ValueError` si quedan 0 pulsos "
           "útiles, en vez de dejar un `.npz` de ceros que reaparece cinco "
           "etapas más abajo como un `IndexError` incomprensible.")

    d.h2("1.2 El área: línea base fuera del pulso y compuerta")
    d.p("Este es **el paso más importante de todo el proyecto**, y el que más "
        "tiempo estuvo mal. `_area_pulso_gate` hace, por evento:")
    d.li("línea base provisoria = mediana de todo el registro (el pulso ocupa "
         "~30% de las muestras, así que la mediana no lo ve);")
    d.li("pico = máximo de |v|; **onset** = último cruce del %.0f%% de la "
         "amplitud antes del pico;" % (100 * C["FRAC_ONSET"]))
    d.li("compuerta = `[onset − %.0f ns, onset + %.0f ns]`;"
         % (C["GATE_PRE"], C["GATE_POST"]))
    d.li("línea base definitiva = **mediana de las muestras pre-onset**, "
         "dejando %.0f ns de guarda antes de que abra la compuerta;"
         % C["GUARDA_BASE"])
    d.li("área = integral trapezoidal dentro de la compuerta, con esa línea "
         "base restada.")
    d.p("Todo se mide contra `t − t[0]`, nunca contra t absoluto: el DRSOSC "
        "escribe el eje relativo a la celda de trigger y en ~2% de los eventos "
        "arranca corrido.")
    d.nota("Por qué esto importa tanto:",
           "la versión anterior promediaba la línea base sobre el **primer 20% "
           "del registro** (40 ns) e integraba los 200 ns completos. Pero el "
           "pulso arranca a los ~25 ns: la ventana de línea base se comía la "
           "subida del propio pulso. El sesgo se multiplica por T/T_b ≈ 5 al "
           "integrar todo, y como el jitter del trigger mueve el pico entre 18 "
           "y 54 ns, **cada evento perdía una fracción distinta de su carga** — "
           "un ensanchamiento multiplicativo puro. Una fracción grande de los "
           "eventos salía incluso con el área de **signo invertido**, y un "
           "`abs()` los reciclaba como buenos (la figura de abajo lo mide "
           "sobre estos datos). Arreglarlo cambió D por un factor 3.4 y bajó "
           "la resolución "
           "de 42.7% a ~21% FWHM. El panel derecho de la figura es esa "
           "comparación, sobre los datos reales.")
    d.figura(F.get("pulso"),
             "Izquierda: forma media del pulso con las dos definiciones de "
             "línea base y la compuerta. Derecha: el cociente integral/amplitud, "
             "que a resolución perfecta sería una constante. La distribución "
             "roja es lo que producía la línea base vieja.")

    d.h2("1.3 Banderas de calidad")
    d.p("`_area_pulso_gate` devuelve una bandera por evento, un bit por "
        "problema: **2** = saturación del digitalizador (|v| ≥ %.0f mV), "
        "**4** = el registro arranca dentro del pulso (queda truncado), "
        "**8** = no había pre-pulso suficiente y la línea base salió del "
        "complemento de la compuerta." % C["V_SATURA"])
    d.p("El bit 8 **no descarta**: el evento sirve, sólo tiene el área ~1.6% "
        "alta. La razón es sutil y vale la pena tenerla presente: descartar "
        "por posición del onset sería un corte **dependiente de la amplitud**, "
        "porque el trigger tiene walk (los pulsos grandes cruzan el umbral "
        "~12 ns antes que los chicos). Un corte así deforma el espectro, que "
        "es exactamente lo que no se puede hacer en una calibración.")
    d.mejorar("el pedestal del bit 8 (~6 mV·ns sobre la compuerta) se traduce "
              "casi 1:1 en D porque el mapa pasa por el origen. Hoy se marca y "
              "se acepta. Se podría corregir evento a evento en vez de "
              "tolerarlo, o cuantificar cuánto de D viene de ahí.")

    d.h2("1.4 Tiempo vivo")
    d.p("El `<Scaler>` de cada evento es la tasa de trigger **real** de la "
        "placa. `tiempo_vivo` la integra y la compara con los eventos "
        "efectivamente guardados:")
    d.code("fracción viva = eventos guardados / ∫ scaler dt\n"
           "T_vivo        = fracción viva × duración de pared")
    tvm, tvf = r["t_vivo_muestra"], r["t_vivo_fondo"]
    d.p("En estos datos la muestra estuvo viva **%.0f s de %.0f s de pared "
        "(%.0f%%)** y el fondo **%.0f s de %.0f s (%.0f%%)**: el DRSOSC vuelca "
        "~35 kB de XML por evento y no da abasto a tasa alta. Normalizar por "
        "tiempo de pared, que es lo que se hacía antes, sobre-restaba el fondo "
        "por un factor 2.5."
        % (tvm, r["t_muestra"], 100 * tvm / max(r["t_muestra"], 1e-9),
           tvf, r["t_fondo"], 100 * tvf / max(r["t_fondo"], 1e-9)))
    d.p("Las marcas de `<Time>` **no** sirven para esto: el DRSOSC las escribe "
        "en ráfagas y miles de eventos comparten el timestamp al milisegundo.")

    d.esp(4)

    # ============================================================ paso 2
    d.h1("2. El espectro neto")
    d.p("**Función:** `spectrum.espectro_neto`. Toma las áreas de las dos "
        "corridas y devuelve el histograma listo para ajustar.")
    d.li("**Binning común obligatorio.** Muestra y fondo se histograman con el "
         "mismo `range` fijo. Sin eso la resta es físicamente inválida: bins "
         "distintos serían mV·ns distintos.")
    d.li("**`RANGO_MAX` fijo, nunca automático.** Unos pocos eventos de área "
         "enorme (pile-up, cósmicos) estiran el máximo y dejan al borde con "
         "pocos bins. Hoy el default de la interfaz es 2500 mV·ns con 500 "
         "canales, y la pestaña de Archivos sugiere 1.6·p99 mirando los datos.")
    d.li("**Normalización por tiempo vivo**, con los pesos `k_m` = %.3f y "
         "`k_f` = %.3f con que entra cada corrida." % (r["k_m"], r["k_f"]))
    d.li("**Barras de error de Poisson propagadas por la resta:** "
         "`error = √(N_m·k_m² + N_f·k_f²)`, con piso 1 para que un bin vacío "
         "no tenga peso infinito.")
    d.li("**El descarte por calidad se aplica igual a muestra y fondo**, si no "
         "la resta deja de significar algo. Acá se descartaron %d eventos de "
         "la muestra y %d del fondo." % r["n_descartados"])
    d.li("**`recortar_negativos=False` para ajustar** (el neto real, con sus "
         "fluctuaciones negativas incluidas); `True` sólo para graficar.")
    d.figura(F.get("espectro"),
             "El espectro neto sobre el que se ajusta todo. La banda verde es "
             "la ventana que eligió el NCC (paso 5).")
    d.nota("Un detalle que sigue abierto:",
           "hay un **exceso de cuentas a bajo canal** que la simulación no "
           "tiene, y que no es el umbral de trigger (es un exceso, no un "
           "déficit). Candidatos: dispersión en el material de alrededor y en "
           "el portafuente, que la simulación no modela, o apilamiento de dark "
           "counts del SiPM. El fondo ya está restado, así que no es eso. "
           "Desde la Etapa 9 se **modela** en vez de esquivarlo (paso 9), pero "
           "la forma que se le asume es hoy el sistemático más grande.")

    d.esp(4)

    # ============================================================ paso 3
    d.h1("3. La física teórica: E_max y la lineshape")
    d.p("**Archivo:** `etapas/physics.py`. Es pura teoría: no toca datos y trabaja "
        "siempre en keV.")
    d.li("`compton_edge_energy(E_gamma)` — la Ec. 1 del paper, **corregida**.")
    d.li("`dsigma_dE(E_e, E_gamma, A)` — Klein-Nishina (Ec. 2) pasada a energía "
         "del electrón vía el jacobiano, o sea la Ec. 10. Vale 0 fuera de "
         "(0, E_max): arriba del borde no hay cinemática, y ese corte **es** el "
         "escalón que se busca.")
    d.li("`FUENTES`, `NOTAS_FUENTE` y la clase `Fuente` — qué líneas gamma "
         "tiene lo que se midió, sus bordes, el ancla `E_ref` y si el modelo "
         "de dos bordes aplica.")
    d.nota("Errata del paper (1 de 3):",
           "la Ec. 1 impresa, E_γ/(1+2α), es la energía del **fotón "
           "retrodispersado**, no la del electrón. La correcta lleva 2α en el "
           "numerador. Se verificó contra los propios datos del paper: con la "
           "fórmula corregida el Co-60 da 963 y 1118 keV, que es lo que "
           "muestra su Fig. 3.")
    d.p("En `A` se absorbe el prefactor π·r0²/(m_e c²) de la Ec. 10: para "
        "correlacionar y ajustar sólo importa la forma.")

    d.h2("La fuente es un parámetro, no una constante")
    d.p("Todo lo que depende de qué se midió sale de un objeto `Fuente`, que "
        "se elige en la interfaz. De él salen las líneas, los bordes, el ancla "
        "`E_ref` del mapa por el origen y si hay uno o dos bordes que ajustar. "
        "Con la fuente cargada acá (%s):" % fu.descripcion())
    filas = [["línea gamma [keV]", "borde Compton E_max [keV]",
              "área esperada [mV·ns]"]]
    for Eg, Em in zip(fu.lineas, fu.bordes):
        area = "%.0f" % (fit["D"] * Em) if fit else "—"
        filas.append(["%.1f" % Eg, "%.1f" % Em, area])
    d.tabla(filas, anchos=[45, 55, 55])
    d.figura(F.get("forward"),
             "Izquierda: la lineshape teórica, con el escalón en E_max. "
             "Derecha: el mismo borde visto con tres resoluciones, y por qué el "
             "método clásico del 50% de altura está sesgado.")

    d.esp(4)

    # ============================================================ paso 4
    d.h1("4. Del espacio de energía al espacio de canal")
    d.p("**Archivo:** `etapas/model.py`. **Función:** `forward_model(canales, E_gamma, "
        "A, edge_ch, sigma)`. Hace dos cosas:")
    d.code("1)  canal → energía:   E_e = canal · (E_max / edge_ch)\n"
           "2)  convolución:       S(canal) = dσ/dE_e(E_e) ⊗ G(σ)")
    d.p("El mapa es lineal **por el origen**: área cero ↔ energía cero, sin "
        "offset. Es la Ec. 8 del paper y una decisión de diseño que todo el "
        "código respeta. Si algún día aparece un pedestal real, se agrega "
        "explícitamente como parámetro, no por la puerta de atrás dejando "
        "libre la ventana de ajuste.")
    d.li("El kernel gaussiano está **normalizado a suma 1**, así que conserva "
         "área: subir σ ensancha y baja el pico sin cambiar las cuentas "
         "totales, y eso desacopla `A` de `σ`.")
    d.li("El radio del kernel se capa a `(len(señal)−1)//2`, porque "
         "`np.convolve(mode=\"same\")` devuelve el largo del **kernel** si el "
         "kernel es más largo que la señal — lo que pasa cuando σ se infla "
         "durante un ajuste que va mal.")
    d.nota("La trampa del padding:",
           "la convolución **no** puede evaluarse sólo dentro de la ventana de "
           "ajuste. `np.convolve` rellena con ceros, así que el modelo pierde "
           "los eventos de más abajo del límite izquierdo que el broadening "
           "debería meter hacia adentro, y el ajuste lo compensa **achicando "
           "σ**: en la prueba sintética de dos bordes daba σ un 18% bajo y D un "
           "3.9% corrido. Por eso `forward_model` tiene `pad=True` por "
           "defecto: extiende la grilla hasta el canal 0 y 4σ a la derecha, "
           "convoluciona, y recorta. En la interfaz se puede destildar para "
           "ver el efecto en vivo.")
    d.h2("Los tres métodos clásicos, y por qué no alcanzan")
    d.p("**Archivo:** `etapas/clasicos.py`. El paper no sólo propone NCC-SA: "
        "mide cuánto mejor es que dos métodos convencionales, y esa "
        "comparación está implementada:")
    d.li("**50% de altura** (`model.punto_mitad_altura`): se lee directo sobre "
         "los datos el canal donde el espectro cae a la mitad del máximo.")
    d.li("**Método de la diferencia** (`metodo_diferencia`): el borde es el "
         "canal donde la derivada es más negativa. La ventana de suavizado se "
         "deduce del ancho del propio flanco, porque con una constante el pozo "
         "de la derivada queda chato cuando el borde es ancho y el mínimo lo "
         "elige el ruido.")
    d.li("**FWHM gaussiano** (`metodo_fwhm_gauss`): se ajusta el flanco con un "
         "escalón borroneado por una gaussiana, `C/2·erfc((ch−B)/(√2σ))` más un "
         "fondo lineal. Es lo que se obtiene de convolucionar un escalón con la "
         "resolución, así que su punto de media altura cae en B por "
         "construcción y el mismo ajuste entrega el FWHM.")
    d.p("Si el borde fuera un escalón perfecto los tres darían lo mismo. Sobre "
        "un borde sintético de posición **conocida**, los tres **sobrestiman**, "
        "y el sesgo crece con la resolución:")
    d.tabla([["σ/B", "50% de altura", "diferencia", "FWHM gaussiano",
              "diferencia SIN ruido"],
             ["0.03", "+0.8%", "+0.7%", "+2.6%", "+0.4%"],
             ["0.06", "+2.2%", "+1.5%", "+4.7%", "+1.4%"],
             ["0.09", "+3.8%", "+2.8%", "+6.3%", "+2.6%"],
             ["0.12", "+5.4%", "+4.4%", "+7.5%", "+3.8%"],
             ["0.18", "+8.5%", "+7.4%", "+9.2%", "+6.2%"]],
            anchos=[18, 34, 30, 34, 40])
    d.p("La última columna es el mismo método sobre el espectro **sin ruido**: "
        "va igual de +0.4% a +6.2%, o sea **el sesgo es del método, no de la "
        "estadística**. La causa es común a los tres: la lineshape de "
        "Klein-Nishina no es un escalón, tiene una acumulación creciente antes "
        "del borde, y ninguno de los tres la contiene. Ésa es la justificación "
        "cuantitativa de ajustar la lineshape entera, que es lo que hace todo "
        "el resto de este documento.")
    d.nota("Una trampa al usarlos:",
           "hay que darles una región **alrededor del borde**, no la ventana de "
           "ajuste del modelo extendido. Esa ventana arranca en ~0.25·B a "
           "propósito porque ese modelo contiene el exceso de bajo canal; los "
           "clásicos no lo contienen, y la caída del exceso es mucho más "
           "empinada que el borde. Medido sobre el Co-60 real: con la ventana "
           "de ajuste el mínimo de la derivada se va al exceso y da D = 0.363 "
           "en vez de 1.265, un −71%. Por eso `comparar` toma el borde del NCC "
           "como ancla.")
    d.esp(4)

    # ============================================================ paso 5
    d.h1("5. Elegir la ventana de ajuste sin criterio humano (NCC)")
    d.p("**Archivo:** `etapas/ncc.py`. Es la Sección II-B del paper y, en la práctica, "
        "lo más útil que aporta: reemplaza «elegir a ojo dónde empieza y "
        "termina el borde» por un criterio reproducible.")
    d.p("La idea: se desliza una ventana sobre el espectro suavizado, se "
        "reinterpola a un número fijo de puntos y se calcula el **coeficiente "
        "de Pearson** (Ec. 15) contra la plantilla teórica. Pearson es inmune "
        "a offset y a escala, así que compara **forma**, no amplitud. Gana la "
        "ventana con ρ más alto.")
    d.p("Sobre estos datos: ventana **[%.0f, %.0f] mV·ns** con "
        "**ρ = %.4f**, y de ahí sale la estimación del borde "
        "**B = %.0f mV·ns**." % (b["lo"], b["hi"], b["rho"], b["edge_ch"]))
    d.figura(F.get("ncc"),
             "Izquierda: la ventana elegida sobre el espectro real. Derecha: "
             "los dos perfiles estandarizados que Pearson compara.")

    d.h2("Dos desviaciones del paper, ambas forzadas por los datos")
    d.p("**1. El barrido 2D del paper está sobre-parametrizado respecto de su "
        "propio modelo.** Si el mapa canal↔energía pasa por el origen (Ec. 8), "
        "fijado B la ventana ya está determinada: `[inicio_rel·B, fin_rel·B]`. "
        "Barrer inicio y fin libres agrega un grado de libertad que el modelo "
        "no tiene, y sobre datos reales se paga: el barrido libre eligió una "
        "ventana de 32 bins en la cola vacía con ρ = 0.984. Por eso el modo "
        "por defecto es `modo=\"origen\"`, que barre **sólo B** — más robusto, "
        "7× más rápido, y sobre sintético recupera B exacto. El barrido del "
        "paper sigue disponible con `modo=\"libre\"`.")
    d.p("**2. El número de puntos de la reinterpolación.** El paper usa 100. "
        "La ventana real tiene ~570 bins y el suavizado Savitzky-Golay de 21 "
        "bins deja sólo ~27 puntos independientes, así que subir el número no "
        "agrega información — pero **bajarlo infla ρ**, porque la "
        "interpolación promedia ruido. Y ρ es justamente el criterio con que "
        "se comparan ventanas. El default acá es `n_puntos=\"auto\"`: un punto "
        "por bin, sin inventar puntos. El número queda fijo durante todo el "
        "barrido, porque si cambiara con B los ρ dejarían de ser comparables.")
    d.nota("El límite real del NCC:",
           "compara forma, y no tiene manera de saber que a bajo canal hay "
           "física que la plantilla no contiene. Sobre estos datos eligió un "
           "inicio que se comía parte del exceso de bajo canal. **La selección "
           "automática resuelve la subjetividad, no la validez del modelo.** "
           "Es la razón por la que el paso 9 modela ese exceso en vez de "
           "confiar en que la ventana lo esquive.")
    d.p("De yapa, el NCC entrega la **semilla** del ajuste [fuera del paper]: "
        "la ventana ganadora queda mapeada punto a punto contra la plantilla, "
        "así que de ahí salen B, σ y D₀. Eso resuelve en la práctica la "
        "fragilidad del optimizador local.")

    d.esp(4)

    # ============================================================ paso 6
    d.h1("6. El ajuste analítico")
    d.p("**Archivos:** `etapas/fit_single.py` (un borde) y `etapas/fit_double.py` (dos "
        "bordes solapados). Cuál se usa lo decide la fuente elegida.")

    d.h2("6.1 Un solo borde")
    d.p("`ajustar_borde_unico` corre `least_squares` sobre {A, edge_ch, σ} "
        "dentro de la ventana. Con una sola línea gamma, D sale directo: "
        "`D = B / E_max`. Es el caso fácil, y por eso una fuente de una sola "
        "línea vale tanto para verificar.")
    d.nota("Es frágil por naturaleza:",
           "se mapeó el paisaje de la función de costo y el mínimo global sí "
           "está en el B verdadero, pero desde un arranque malo el optimizador "
           "queda ~10% off con σ inflado. Mitigación: sembrarlo con la semilla "
           "del NCC, que en la prueba con verdad conocida bajó el error en B de "
           "−26 canales a −0.2 sobre 600. La solución robusta es el paso 7.")

    d.h2("6.2 Dos bordes solapados")
    d.p("Con dos líneas gamma (Co-60, Na-22) el espectro tiene dos bordes que "
        "se superponen. `modelo_dos_bordes` implementa las Ec. 11-13, con el "
        "segundo borde amarrado al primero:")
    d.code("CE_ch,2 = B1 + D·ΔE          Θ = {A1, A2, B1, σ1, σ2, D}")
    d.nota("Errata del paper (2 de 3):",
           "la Ec. 12 define ΔE como la diferencia «entre las líneas gamma». "
           "Está mal: D convierte energía del **electrón** a canal, y lo que se "
           "separa en canal son los **bordes**. Para Co-60 son 154.7 keV y no "
           "159.3. Con el ΔE de los gammas, D sale ~3% corrido. En el código, "
           "`fit_double.delta_E`.")
    if fu.doble:
        sep = (fu.bordes[1] - fu.bordes[0]) / fu.bordes[0]
        d.p("Con la fuente cargada, los dos bordes están a %.0f%% uno de otro "
            "(%.1f y %.1f keV)." % (100 * sep, fu.bordes[0], fu.bordes[1]))
    d.p("**El problema de fondo:** en el Co-60 la separación entre los bordes "
        "es comparable a σ, así que la descomposición interna no tiene "
        "información aunque la posición del conjunto sí. Barriendo variantes, "
        "B1 se mueve ±0.5% pero D libre se mueve ±8%. Por eso se ajusta con "
        "`ligar_origen=True` (D = B1/E_max,1), que da D estable en todas las "
        "variantes, y está justificado porque el test del origen sale limpio: "
        "el offset extrapolado a E = 0 es −2.3% de B1.")
    d.mejorar("una fuente con los bordes bien separados (Na-22 los tiene en "
              "341 y 1062 keV) haría honesta la descomposición de dos "
              "componentes, en vez de tener que ligarla.")

    d.esp(4)

    # ============================================================ paso 7
    d.h1("7. El optimizador global (recocido simulado)")
    d.p("**Archivo:** `etapas/sa.py`. Es la Sección II-A del paper: Metropolis con "
        "enfriamiento geométrico, para escapar de mínimos locales.")
    d.p("Con 3 parámetros gana claro. Desde un arranque deliberadamente malo "
        "(B₀ = 0.70·B, σ₀ = 3σ), `least_squares` da B = 343 sobre un valor "
        "verdadero de 600, y el recocido da **600.2** con cualquiera de los dos "
        "objetivos.")
    d.nota("Dos constantes del paper NO son transferibles:",
           "Metropolis compara la temperatura con ΔE, así que T sólo tiene "
           "sentido **en las unidades del costo**. Medido acá, la calibración "
           "automática pide T_init ≈ 2·10⁹ con la Ec. 14 (suma de cuadrados) y "
           "≈ 3.2 con el error NCC: nueve órdenes de magnitud de diferencia "
           "contra el T_init = 2000 del paper. Lo mismo con su corte absoluto "
           "|ΔE| < 0.001. Por eso `recocido_simulado` usa `T_init=\"auto\"` "
           "(calibra para ~80% de aceptación inicial) y tolerancia relativa.")
    d.nota("Errata del paper (3 de 3):",
           "la Ec. 14 dice que el recocido minimiza mínimos cuadrados, pero el "
           "texto de Metropolis dice que el costo es «NCC-based matching "
           "error». Hay que elegir uno. Se verificó que con el objetivo NCC las "
           "amplitudes son **indeterminables** (multiplicar A por 7.3 cambia el "
           "costo en exactamente 0), así que la Ec. 14 es la lectura "
           "utilizable.")
    d.p("**Su aporte real fue diagnóstico.** Al buscar de verdad el óptimo "
        "global mostró que Θ = {A1, A2, B1, σ1, σ2, D} está **sub-determinado**: "
        "cada libertad de más el óptimo la usa para algo sin sentido físico "
        "(offset de −5928 mV·ns en E = 0; σ1 > σ2 cuando la resolución tiene "
        "que crecer con E; A1 = 0 cuando el Co-60 emite las dos líneas casi "
        "igual). La conclusión no fue cambiar de optimizador sino **restringir "
        "el modelo con física**, que es lo que hace el paso siguiente.")
    d.p("El recocido es estocástico: una sola corrida no es un resultado. Se "
        "usa `sa_multiarranque`.")

    d.esp(4)

    # ============================================================ paso 8
    d.h1("8. Forward folding con la simulación de Geant4")
    d.p("**Archivo:** `ajuste/sim_template.py`. Cambia la plantilla analítica por el "
        "espectro simulado, dejando el resto del esquema igual:")
    d.code("S_model(ch) = A · Σ_j  h_j · G(ch ; D·E_j , D·σ_E(E_j))\n"
           "σ_E(E) = √(α²·E + β²)          [α: estadística de fotoelectrones]")
    d.p("Qué contiene la simulación que la fórmula analítica no: scattering "
        "Compton **múltiple**, escape de los gammas dispersados (depende del "
        "tamaño del centellador), el quenching de Birks a baja energía y la "
        "eficiencia relativa de las dos líneas. Todo eso el ajuste analítico "
        "lo tiene que absorber inflando σ.")
    if n.get("lineas_sim"):
        d.p("La simulación cargada acá genera %s, y eso se **verifica** contra "
            "la fuente elegida antes de usarla."
            % ", ".join("%.1f keV (%.0f%%)" % (E, 100 * f)
                        for E, f in n["lineas_sim"]))
    d.nota("Un modo de falla silencioso que ahora está tapado:",
           "la plantilla **es** la fuente. Ajustar un espectro de Cs-137 contra "
           "una simulación de Co-60 no produce ningún error: produce un D "
           "perfectamente formateado y equivocado, porque la plantilla tiene "
           "sus bordes en otro lado. `sim_template.lineas_primarias` lee la "
           "rama `E0_keV` y la interfaz bloquea el forward folding si no "
           "coincide con la fuente elegida.")
    d.p("Dos escalas que **no hay que mezclar**: la rama `Lvis_keVee` es la luz "
        "visible con Birks ya aplicado —lo que ve el SiPM— y da keVee; la rama "
        "`Edep_keV` es la energía depositada y da keV. Su cociente reproduce el "
        "quenching de la simulación, pero eso **no** es un control "
        "independiente: Lvis es Edep multiplicada por el quenching, así que "
        "ajustar los mismos datos con las dos plantillas tiene que dar ese "
        "cociente.")
    d.p("Detalle de implementación que importa para la velocidad: el binning de "
        "la plantilla es **irrelevante** para el resultado, porque el kernel "
        "gaussiano es ~90 veces más ancho que un bin. Rebinnear de 700 a 175 "
        "bins deja D igual en la cuarta cifra y acelera 6×, que es lo que hace "
        "viables el recocido y el MCMC.")
    d.figura(F.get("ajuste"),
             "Izquierda: las dos componentes de la simulación, separadas por "
             "`NCompton`. El scattering múltiple deposita energía **por encima** "
             "del borde, y ése es exactamente el relleno que la plantilla "
             "analítica no tiene. Derecha: el ajuste del paso 9 sobre los datos.")

    d.esp(4)

    # ============================================================ paso 9
    d.h1("9. El modelo extendido: los sistemáticos como parámetros")
    d.p("**Archivo:** `ajuste/modelo_ext.py`. Es el resultado vigente y la idea "
        "central del análisis actual.")
    d.p("Hasta acá quedaban dos incertidumbres que se trataban **a mano**: "
        "cuánto scattering múltiple hay (se reportaba la diferencia entre "
        "plantillas como sistemático) y el exceso de bajo canal (se resolvía "
        "cortando la ventana, lo que tira estadística y deja una dependencia "
        "con dónde se corta). Las dos entran ahora al modelo:")
    d.code("S(ch) = A·[ w·h_1scatt + (1−w)·h_múltiple ] ⊗ G(σ_E)  +  C·exp(−ch/λ)\n"
           "σ_E(E) = √( a²·E² + α²·E + β² )")
    d.p("`w` es cuánto scattering simple **piden los datos**, y el término "
        "exponencial es el exceso. Θ = {A, D, α, w, C, λ, β, a}; se ajusta el "
        "subconjunto que uno declare libre.")
    if fit:
        d.p("Sobre estos datos, con la ventana [%.0f, %.0f]: **w = %.3f ± %.3f** "
            "contra **%.3f** de la simulación, o sea coinciden dentro de "
            "%.1fσ. El χ²/ndf es %.3f con %d grados de libertad."
            % (fit["intervalo"][0], fit["intervalo"][1], fit["w"],
               fit["err"].get("w", float("nan")), fit["w_nat"],
               abs(fit["w"] - fit["w_nat"]) / max(fit["err"].get("w", 1e-9), 1e-9),
               fit["chi2_red"], fit["ndf"]))
        d.p("Que la geometría simulada reproduzca cuánto scattering múltiple "
            "hay es un resultado en sí mismo. Y `w` está realmente medido, no "
            "es decorativo: fijarlo en el valor de la simulación no cambia el "
            "χ², mientras que fijarlo en 1 —que es lo que asume la plantilla "
            "analítica del paper— lo empeora en Δχ² ≈ +13 y sube D un 3.4%.")
    d.p("Lo que se gana: el error de D deja de ser «la dispersión entre tres "
        "ajustes que no se pueden comparar» y pasa a ser el ancho marginal del "
        "posterior, que ya incluye la ignorancia sobre `w` y sobre el exceso "
        "vía la correlación corr(D, w) ≈ +0.96. Eso sí es propagable.")
    d.p("Lo que **no** se gana: si la simulación no contiene un proceso, `w` no "
        "lo inventa. `w` mide cuánto pesa el scattering múltiple **simulado**, "
        "no valida la geometría.")
    d.nota("El χ² de Pearson, no el de Neyman:",
           "calcular la barra de error con las cuentas **observadas** sesga D "
           "hacia abajo, porque un bin que fluctuó bajo recibe σ más chico y "
           "por lo tanto más peso. Con la estadística de este experimento "
           "(~15-25 cuentas/bin) el sesgo era comparable al error estadístico "
           "que se reportaba. Medido sobre 20 realizaciones con verdad "
           "conocida: Neyman da −0.53% de sesgo y cobertura 8/20; Pearson da "
           "−0.06% y 15/20, que es lo que corresponde. Por eso "
           "`spectrum.varianza_predicha` usa la varianza del **modelo** y el "
           "ajuste itera tres veces, ya que esa varianza depende de los "
           "parámetros.")
    d.figura(F.get("recta"),
             "Izquierda: residuos normalizados del ajuste — planos, o sea el "
             "exceso no está tapando estructura del borde. Derecha: la recta "
             "de calibración resultante, con los bordes de esta fuente.")

    d.esp(4)

    # ============================================================ paso 10
    d.h1("10. El posterior (MCMC)")
    d.p("**Archivo:** `ajuste/mcmc.py`. Metropolis-Hastings propio, sin dependencias "
        "externas.")
    d.li("La propuesta usa la **covarianza del ajuste** escalada por 2.4²/d, "
         "que es el óptimo para un objetivo gaussiano en d dimensiones. Así la "
         "cadena viaja **a lo largo** de la degeneración D–α en vez de chocarse "
         "contra ella.")
    d.li("Las cadenas arrancan **dispersas**. No es un detalle: R-hat compara la "
         "varianza entre cadenas contra la interna, y si todas arrancan en el "
         "mismo punto R-hat da ~1 aunque ninguna haya explorado nada.")
    d.li("Diagnósticos: `gelman_rubin` (R-hat), `autocorr_integrada` (τ y por lo "
         "tanto n_eff) y la fracción de aceptación.")
    d.li("Cuando la varianza depende de los parámetros (Pearson), el posterior "
         "lleva el término `ln σ²(θ)`, que con varianza fija era constante y se "
         "podía ignorar.")
    d.p("`python codigo/ajuste/mcmc.py` corre el test con **verdad conocida**: genera datos "
        "con la simulación a D y α conocidos, con ruido de Poisson a la "
        "estadística real de este experimento, y verifica que los parámetros "
        "caen dentro del 68% del posterior. Ahí se ve también el sesgo de "
        "Neyman en vivo.")
    d.p("Además se hizo un **cierre con verdad conocida** de la cadena "
        "completa: 20 realizaciones Poisson a la estadística real dan sesgo en "
        "D de −0.10%, rms 0.0116 contra un error reportado de 0.0127, y "
        "cobertura 17/20.")

    d.esp(4)

    # ============================================================ paso 11
    d.h1("11. El resultado y su presupuesto de error")
    d.p("**Archivos:** `calibracion.py` (el número congelado del Co-60, que es "
        "lo que importan otros proyectos) y la pestaña **★ Calibración** de "
        "`gui.py`, que lo recalcula sobre los archivos y la fuente que estén "
        "cargados.")
    if fit:
        d.p("Sobre estos datos, la conversión es:")
        d.code("E [keVee] = área [mV·ns] / %.4f          "
               "área [mV·ns] = E [keVee] × %.4f\n"
               "σ_E(E) = %.3f · √E  keV        (E en keV)"
               % (fit["D"], fit["D"], fit["alpha"]))
    d.h2("De dónde sale cada término del error")
    d.p("La pestaña de Calibración **mide** el presupuesto en vez de citarlo: "
        "corre las variantes sobre los mismos datos y compara.")
    d.tabla([
        ["término", "cómo se mide", "cuánto da (Co-60, este detector)"],
        ["estadística", "ancho del posterior; ya incluye w y el exceso por "
         "correlación", "≈ 1.1%"],
        ["ventana", "barrido del límite inferior de la ventana de ajuste",
         "≈ 0.9%"],
        ["forma del exceso", "exponencial contra ley de potencia: describen los "
         "datos igual de bien", "**≈ 2.2-2.4% ← el que manda**"],
        ["definición del área", "barrido de compuerta (45/55/70/90 ns), onset "
         "(5/10/20%) y descarte", "≈ 1.0%"],
        ["total", "suma en cuadratura", "≈ 2.8-2.9%"],
    ], anchos=[32, 75, 58])
    d.p("El único término que se **arrastra** en vez de medirse en vivo es el "
        "de la definición del área, porque cada variante obliga a re-parsear "
        "los XML (~40 s cada una).")
    d.nota("Lo que esto dice:",
           "**el sistemático domina y ya se sabe cuál es.** La forma que se le "
           "asume al exceso de bajo canal no se cierra ajustando mejor: "
           "exponencial y ley de potencia dan χ²/ndf indistinguibles y D con "
           "~2.3% de diferencia. Hay que **medir** el exceso.")
    d.p("Cada resultado se guarda en `.cache/calibraciones.json` con una clave "
        "que incluye la fuente y la huella de los archivos, y la pestaña "
        "superpone las fuentes ya medidas. Ojo con qué prueba eso: cada fuente "
        "entregó su D con el mapa por el origen, así que su punto cae exacto "
        "sobre su propia recta. **Lo que hay que mirar es la dispersión entre "
        "fuentes**, y eso es el test de linealidad.")

    d.salto()

    # ============================================================ la interfaz
    d.h1("12. La interfaz")
    d.p("**Archivo:** `gui.py`. Catorce pestañas, una por etapa, cada una con "
        "controles a la izquierda, gráfico al centro y una explicación a la "
        "derecha. No muestra PNGs guardados: **recalcula**, que es lo que "
        "permite mover un parámetro y ver qué se rompe.")
    d.tabla([
        ["pestaña", "qué permite ver"],
        ["Archivos", "elegir muestra, fondo, simulación y **fuente**; "
         "diagnóstico de tiempo vivo y canal; verifica que la simulación "
         "corresponda a la fuente"],
        ["0 · Espectro", "el neto, con selector `gate`/`legacy` para ver el "
         "antes/después del arreglo del área"],
        ["1 · Klein-Nishina", "la lineshape teórica y los bordes"],
        ["2 · Forward model", "el efecto de σ y el sesgo del método del 50%; "
         "destildar «pad» muestra la trampa de la convolución"],
        ["3 · Un borde", "el ajuste local y su fragilidad"],
        ["4 · NCC", "la selección automática de ventana; modo «libre» muestra "
         "por qué el barrido 2D del paper falla"],
        ["6 · Dos bordes / Un borde", "el modelo que corresponda a la fuente; "
         "destildar las ligaduras muestra la sub-determinación"],
        ["5 · Recocido", "el óptimo global contra el local"],
        ["7 · Geant4", "las plantillas simuladas"],
        ["8 · NCC→SA→MCMC", "la cadena entera y el posterior; la distancia "
         "entre plantillas **es** el sistemático"],
        ["9 · El área", "el hallazgo del paso 1, en cuatro paneles"],
        ["10 · Modelo ext.", "el modelo del paso 9, con w y el exceso libres"],
        ["★ Calibración", "el resultado sobre los datos cargados, con el "
         "presupuesto de error medido y el registro por fuente"],
        ["? Parámetros", "glosario: qué es B, σ, D, ρ, Θ, α, keVee vs keV…"],
    ], anchos=[42, 123], cabecera=True)
    d.p("Los cálculos pesados corren en un hilo aparte y hay una barra de "
        "progreso abajo: indeterminada cuando el paso no sabe cuánto falta, y "
        "con porcentaje y tiempo restante cuando sí (el parseo de los XML mide "
        "bytes leídos; la pestaña de Calibración cuenta los ajustes antes de "
        "empezar; los MCMC reportan cadena y paso).")

    # ============================================================ resumen
    d.salto()
    d.h1("13. Qué revisar y qué mejorar")
    d.p("En orden de impacto sobre el número final.")
    d.h2("1. Medir el exceso de bajo canal — es el sistemático dominante")
    d.p("Hoy vale ~2.3% de D y no se cierra con más estadística ni con mejor "
        "ajuste: exponencial y ley de potencia describen los datos igual de "
        "bien. Lo que lo resuelve es **una corrida corta con la fuente "
        "apantallada, o movida de lugar**, que diga si viene de dispersión en "
        "el material de alrededor. Es la medición de mayor rendimiento "
        "pendiente.")
    d.h2("2. Una segunda fuente, para dejar de asumir la linealidad")
    d.p("Con una sola fuente la recta se traza con un punto y el origen. El "
        "registro por fuente ya existe y la interfaz superpone lo que haya; "
        "falta la medición. Para cerrarla bien hace falta también una "
        "**simulación de Geant4 de esa fuente** con esta geometría: sin ella "
        "queda la ruta analítica, con su sesgo conocido por no tener "
        "scattering múltiple.")
    d.h2("3. Recuperar el tiempo vivo")
    d.p("Con ~40% de tiempo vivo se tiró el 60% de la estadística. Bajar la "
        "tasa de trigger o acelerar el guardado da un factor √2.5 en el error "
        "estadístico sin tocar nada del análisis.")
    d.h2("4. La definición del área")
    d.p("Vale ~1.0% y hoy se arrastra del barrido de la Etapa 9 en vez de "
        "re-medirse. Sería razonable hacer que el barrido se pueda correr desde "
        "la interfaz (es lento porque re-parsea, pero se puede cachear una vez "
        "por definición). También queda abierto qué hacer con el pedestal de "
        "los eventos sin pre-pulso (bandera 8): hoy se marca y se acepta.")
    d.h2("5. Cosas menores pero reales")
    d.li("El término `a·E` y el `β` de la resolución se van a cero: con dos "
         "bordes tan cercanos no hay palanca para separar los tres términos. "
         "Con una fuente de energía muy distinta sí la habría.")
    d.li("El ajuste de dos bordes tiene que ligar el origen para ser estable. "
         "Con una fuente de bordes separados dejaría de hacer falta.")
    d.li("La simulación se cargó con la geometría actual; si cambia el montaje "
         "(portafuente, blindaje, distancia), `w` lo va a notar antes que "
         "cualquier otra cosa. Vale la pena mirarlo después de cada cambio "
         "mecánico.")

    d.esp(10)
    # `p` escapa el marcado, asi que aca va texto plano con el mini-formato de
    # `marcar` (backticks para codigo), no HTML crudo.
    d.p("Fin. Este documento se regenera con `python codigo/documentacion.py`; los "
        "números y las figuras se recalculan cada vez, así que refleja el "
        "código del día en que se corrió.", "pie")
    return d.items


# ======================================================================
#  6) armado
# ======================================================================
def pie_de_pagina(canvas, doc, fuentes):
    canvas.saveState()
    canvas.setFont(fuentes[0], 7.6)
    canvas.setFillGray(0.45)
    canvas.drawString(20 * 2.8346, 12 * 2.8346,
                      "Calibración por Compton edge — análisis paso a paso")
    canvas.drawRightString(190 * 2.8346, 12 * 2.8346, "%d" % doc.page)
    canvas.restoreState()


def main(sin_figuras=False, rapido=False):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate

    t0 = time.time()
    print("=" * 70)
    print("DOCUMENTACION DEL ANALISIS -> " + SALIDA)
    print("=" * 70)

    print("  reuniendo los numeros (se calculan, no se tipean) ...")
    n = reunir_numeros(rapido=rapido or sin_figuras)
    F = {} if sin_figuras else figuras(n, rapido=rapido)

    fuentes = registrar_fuentes()
    E = estilos(fuentes)
    print("  armando el PDF ...")
    doc = SimpleDocTemplate(SALIDA, pagesize=A4,
                            leftMargin=22 * mm, rightMargin=22 * mm,
                            topMargin=18 * mm, bottomMargin=20 * mm,
                            title="Calibracion por Compton edge - analisis paso a paso",
                            author="documentacion.py")
    contenido = construir(n, F, E, fuentes)
    doc.build(contenido,
              onFirstPage=lambda c, d_: pie_de_pagina(c, d_, fuentes),
              onLaterPages=lambda c, d_: pie_de_pagina(c, d_, fuentes))
    print("  listo: %s (%.0f kB, %.0f s)"
          % (rutas.relativa(SALIDA), os.path.getsize(SALIDA) / 1e3,
             time.time() - t0))


if __name__ == "__main__":
    main(sin_figuras="--sin-figuras" in sys.argv,
         rapido="--rapido" in sys.argv)
