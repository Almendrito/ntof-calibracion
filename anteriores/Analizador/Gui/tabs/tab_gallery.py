import numpy as np
from Gui.tabs.base_plot import BasePlotTab

class TabGallery(BasePlotTab):
    """Pestaña 4: Maneja las visualizaciones de múltiples pulsos (Galería)."""
    
    def _corregir_y_orientar(self, matriz, baseline_muestras, polaridad):
        """Método auxiliar para corregir baseline y polaridad de un grupo de pulsos."""
        mats_corr = np.array([
            (v - (np.mean(v[:baseline_muestras]) if baseline_muestras > 0 else v[0]))
            for v in matriz
        ])
        return mats_corr if polaridad == "Positivo" else -mats_corr

    def dibujar_superposicion(self, matriz, t_ns, n_superponer, fuente, polaridad, baseline, metodo, idx_i=None, idx_f=None):
        self.limpiar_y_dibujar()
        n = min(n_superponer, matriz.shape[0])
        indices = np.random.choice(matriz.shape[0], n, replace=False)
        
        cmap = __import__("matplotlib").cm.get_cmap("plasma", n)
        
        for k, idx in enumerate(indices):
            v_plot = self._corregir_y_orientar(matriz[idx:idx+1], baseline, polaridad)[0]
            self.ax.plot(t_ns, v_plot, color=cmap(k), alpha=0.55, linewidth=0.7)

        if idx_i is not None and idx_f is not None:
            self.ax.axvspan(t_ns[idx_i], t_ns[idx_f], alpha=0.12, color="orange",
                            label=f"Ventana integración\n[{idx_i}:{idx_f}]")
            self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white", fontsize=8)

        self.aplicar_estilo_oscuro(
            titulo=f"{n} pulsos superpuestos – {fuente} [método: {metodo}]",
            xlabel="Tiempo (ns)", ylabel="Voltaje (V)"
        )
        self.canvas.draw()

    def dibujar_promedio(self, matriz, t_ns, fuente, polaridad, baseline, idx_i=None, idx_f=None):
        self.limpiar_y_dibujar()
        mats_corr = self._corregir_y_orientar(matriz, baseline, polaridad)
        promedio = np.mean(mats_corr, axis=0)
        sigma = np.std(mats_corr, axis=0)

        self.ax.fill_between(t_ns, promedio - sigma, promedio + sigma, color="cyan", alpha=0.25, label="± 1σ")
        self.ax.plot(t_ns, promedio, color="cyan", linewidth=1.5, label="Promedio")

        if idx_i is not None and idx_f is not None:
            self.ax.axvspan(t_ns[idx_i], t_ns[idx_f], alpha=0.15, color="orange", label="Ventana integración")

        self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white", fontsize=8)
        self.aplicar_estilo_oscuro(
            titulo=f"Pulso Promedio ± σ – {fuente} (N={matriz.shape[0]})",
            xlabel="Tiempo (ns)", ylabel="Voltaje (V)"
        )
        self.canvas.draw()

    def dibujar_comparacion(self, m_muestra, m_fondo, t_ns, n_superponer, polaridad, baseline, idx_i=None, idx_f=None):
        self.limpiar_y_dibujar()
        n = min(n_superponer, m_muestra.shape[0], m_fondo.shape[0])
        
        idx_m = np.random.choice(m_muestra.shape[0], n, replace=False)
        idx_f_rnd = np.random.choice(m_fondo.shape[0], n, replace=False)

        pulsos_m = self._corregir_y_orientar(m_muestra[idx_m], baseline, polaridad)
        pulsos_f = self._corregir_y_orientar(m_fondo[idx_f_rnd], baseline, polaridad)

        for v in pulsos_m:
            self.ax.plot(t_ns, v, color="cyan", alpha=0.35, linewidth=0.6)
        for v in pulsos_f:
            self.ax.plot(t_ns, v, color="red", alpha=0.35, linewidth=0.6)

        # Líneas fantasma para leyenda
        self.ax.plot([], [], color="cyan", label="Muestra", linewidth=1.5)
        self.ax.plot([], [], color="red", label="Background", linewidth=1.5)

        if idx_i is not None and idx_f is not None:
            self.ax.axvspan(t_ns[idx_i], t_ns[idx_f], alpha=0.10, color="orange", label="Ventana integración")

        self.ax.legend(facecolor="#111", edgecolor="white", labelcolor="white", fontsize=8)
        self.aplicar_estilo_oscuro(
            titulo=f"Muestra (cian) vs Background (rojo) – {n} pulsos c/u",
            xlabel="Tiempo (ns)", ylabel="Voltaje (V)"
        )
        self.canvas.draw()