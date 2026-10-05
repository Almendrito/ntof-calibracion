import numpy as np
from Gui.tabs.base_plot import BasePlotTab

class TabHistMuestra(BasePlotTab):
    """Pestaña 5: Espectro de Energía de la Muestra."""
    def __init__(self, parent):
        super().__init__(parent)
        self.bar_patches = []
        
        # Tooltip para el histograma
        self.annot = self.ax.annotate(
            "", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
            bbox=dict(boxstyle="round", fc="black", ec="white", alpha=0.8),
            arrowprops=dict(arrowstyle="->", color="white"),
            color="white", fontweight="bold"
        )
        self.annot.set_visible(False)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)

    def actualizar(self, areas_muestra, num_bins):
        self.limpiar_y_dibujar()
        if not areas_muestra:
            return

        _, _, self.bar_patches = self.ax.hist(
            areas_muestra, bins=num_bins,
            color="magenta", edgecolor="white", alpha=0.75, linewidth=0.4
        )
        self.ax.set_xlim(0, np.max(areas_muestra) * 1.1)
        
        self.aplicar_estilo_oscuro(
            titulo="Espectro de Energía (Solo Muestra)",
            xlabel="Área Integrada (V·ns)", ylabel="Cuentas"
        )
        self.canvas.draw()

    def on_hover(self, event):
        visible = self.annot.get_visible()
        if event.inaxes == self.ax:
            for bar in self.bar_patches:
                cont, _ = bar.contains(event)
                if cont:
                    x = bar.get_x() + bar.get_width() / 2
                    y = bar.get_height()
                    self.annot.xy = (x, y)
                    self.annot.set_text(f"Área: {x:.4f} V·ns\nCuentas: {int(y)}")
                    self.annot.set_visible(True)
                    self.canvas.draw_idle()
                    return
        if visible:
            self.annot.set_visible(False)
            self.canvas.draw_idle()


class TabHistComparacion(BasePlotTab):
    """Pestaña 6: Superposición de histogramas (Muestra vs Background)."""
    def actualizar(self, areas_muestra, areas_fondo, num_bins):
        self.limpiar_y_dibujar()
        if not areas_muestra or not areas_fondo:
            return

        mn = min(min(areas_muestra), min(areas_fondo))
        mx = max(max(areas_muestra), max(areas_fondo))
        bins_c = np.linspace(mn, mx, num_bins + 1)
        
        cm, _ = np.histogram(areas_muestra, bins=bins_c)
        cb, _ = np.histogram(areas_fondo, bins=bins_c)

        self.ax.hist(bins_c[:-1], bins_c, weights=cm.astype(float),
                     color="cyan", edgecolor="white", alpha=0.6, label="Muestra", linewidth=0.4)
        self.ax.hist(bins_c[:-1], bins_c, weights=cb.astype(float),
                     color="red", edgecolor="white", alpha=0.6, label="Background", linewidth=0.4)
                     
        self.ax.set_xlim(0, mx * 1.1)
        self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white")
        
        self.aplicar_estilo_oscuro(
            titulo="Comparación: Muestra vs Background",
            xlabel="Área Integrada (V·ns)", ylabel="Cuentas"
        )
        self.canvas.draw()


class TabHistNeto(BasePlotTab):
    """Pestaña 7: Resta directa de cuentas (Muestra - Background)."""
    def actualizar(self, areas_muestra, areas_fondo, num_bins):
        self.limpiar_y_dibujar()
        if not areas_muestra or not areas_fondo:
            return

        mn = min(min(areas_muestra), min(areas_fondo))
        mx = max(max(areas_muestra), max(areas_fondo))
        bins_c = np.linspace(mn, mx, num_bins + 1)
        
        cm, _ = np.histogram(areas_muestra, bins=bins_c)
        cb, _ = np.histogram(areas_fondo, bins=bins_c)

        neta = np.maximum(cm.astype(float) - cb.astype(float), 0)
        
        self.ax.bar(bins_c[:-1], neta, width=np.diff(bins_c), align="edge",
                    color="#00FF00", edgecolor="white", alpha=0.85, linewidth=0.4)
                    
        self.ax.set_xlim(0, mx * 1.1)
        self.aplicar_estilo_oscuro(
            titulo="Espectro Neto (Muestra − Background)",
            xlabel="Área Integrada (V·ns)", ylabel="Cuentas Netas"
        )
        self.canvas.draw()