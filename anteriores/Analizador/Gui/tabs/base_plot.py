from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class BasePlotTab(ttk.Frame):
    """
    Clase base para todas las pestañas que contienen un gráfico de Matplotlib.
    Configura el modo oscuro y el toolbar automáticamente.
    """
    def __init__(self, parent):
        super().__init__(parent)
        
        # 1. Crear figura de Matplotlib
        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.fig.patch.set_facecolor("#1a1a1a") # Fondo oscuro exterior
        
        # 2. Crear el eje principal (ax)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor("#101010")        # Fondo oscuro interior
        self.fig.tight_layout()
        
        # 3. Integrar en Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        
        # 4. Añadir la barra de herramientas (zoom, pan, guardar)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        
    def aplicar_estilo_oscuro(self, titulo, xlabel, ylabel):
        """Método de utilidad para aplicar textos en blanco rápidamente."""
        self.ax.set_title(titulo, color="white", fontsize=11)
        self.ax.set_xlabel(xlabel, color="white")
        self.ax.set_ylabel(ylabel, color="white")
        self.ax.tick_params(colors="white")
        
    def limpiar_y_dibujar(self):
        """Limpia el gráfico actual (útil para redibujar)"""
        self.ax.clear()
        self.ax.set_facecolor("#101010")