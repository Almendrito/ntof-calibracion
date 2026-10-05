import numpy as np
from Gui.tabs.base_plot import BasePlotTab

class TabOriginal(BasePlotTab):
    def actualizar(self, t, v, offset, is_ns):
        self.limpiar_y_dibujar()
        xu = "ns" if is_ns else "muestra"
        self.ax.plot(t, v, color="yellow", linewidth=0.8)
        self.ax.axhline(offset, color="red", linestyle="--", label=f"Baseline: {offset:.4f} V")
        self.aplicar_estilo_oscuro("Pulso Original", f"Tiempo ({xu})", "Voltaje (V)")
        self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white", fontsize=8)
        self.canvas.draw()

class TabBaseline(BasePlotTab):
    def actualizar(self, t, v_corr, offset, is_ns):
        self.limpiar_y_dibujar()
        xu = "ns" if is_ns else "muestra"
        self.ax.plot(t, v_corr, color="cyan", linewidth=0.8)
        self.ax.axhline(0, color="white", linewidth=0.5, linestyle=":")
        self.aplicar_estilo_oscuro(f"Pulso Corregido – Offset: {offset:.4f} V", f"Tiempo ({xu})", "Voltaje Corregido (V)")
        self.canvas.draw()

class TabIntegracion(BasePlotTab):
    # CRÍTICO: Ahora acepta 8 argumentos (incluyendo v_teorico)
    def actualizar(self, t, v_corr, area, idx_i, idx_f, umbral_amp, is_ns, v_teorico=None):
        self.limpiar_y_dibujar()
        xu = "ns" if is_ns else "muestra"
        
        # 1. Dibujar pulso real
        self.ax.plot(t, v_corr, color="#00FF00", linewidth=0.8, label="Pulso real")

        # 2. Dibujar ajuste exponencial si existe
        if v_teorico is not None and idx_i is not None:
            # Graficamos la curva teórica superpuesta
            self.ax.plot(t[idx_i:], v_teorico[idx_i:], color="#FF00FF", 
                         linewidth=1.2, linestyle="-.", label="Ajuste Exponencial")

        # 3. Dibujar zona de integración
        if idx_i is not None and idx_f is not None:
            t_int = t[idx_i:idx_f + 1]
            v_int = v_corr[idx_i:idx_f + 1]
            self.ax.fill_between(t_int, v_int, color="orange", alpha=0.5, label="Área Integrada")
            self.ax.axvline(t[idx_i], color="white", linestyle="--", alpha=0.4)
            self.ax.axvline(t[idx_f], color="white", linestyle="--", alpha=0.4)

        self.ax.axhline(0, color="white", linewidth=0.5)
        self.aplicar_estilo_oscuro(f"Integración | Área = {area:.5f} V·ns", f"Tiempo ({xu})", "Voltaje (V)")
        self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white", fontsize=8)
        self.canvas.draw()