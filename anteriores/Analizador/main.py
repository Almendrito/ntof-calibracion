import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Importamos nuestro Modelo (Lógica)
from core.data_manager import DataManager
from core import processing

# Importamos nuestras Vistas (Interfaz)
from Gui.control_panel import ControlPanel
from Gui.tabs.tab_individual import TabOriginal, TabBaseline, TabIntegracion
from Gui.tabs.tab_gallery import TabGallery
from Gui.tabs.tab_spectra import TabHistMuestra, TabHistComparacion, TabHistNeto

class AnalizadorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Analizador Espectroscopía – BC404 / XP2262B [Arquitectura Modular]")
        self.root.geometry("1480x1000")

        # 1. Inicializar el Administrador de Datos
        self.data = DataManager()
        self.idx_pulso_actual = 0

        # 2. Construir la Interfaz Gráfica
        self._construir_ui()

    def _construir_ui(self):
        # Panel Izquierdo
        self.panel_control = ControlPanel(self.root, controller=self, width=280)
        self.panel_control.pack(side=tk.LEFT, fill=tk.Y)

        # Panel Derecho (Pestañas)
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=6, pady=8)

        # Instanciar las pestañas modulares
        self.tab_orig = TabOriginal(self.notebook)
        self.tab_base = TabBaseline(self.notebook)
        self.tab_int  = TabIntegracion(self.notebook)
        self.tab_gal  = TabGallery(self.notebook)
        self.tab_hm   = TabHistMuestra(self.notebook)
        self.tab_hc   = TabHistComparacion(self.notebook)
        self.tab_hn   = TabHistNeto(self.notebook)

        # Añadirlas al Notebook
        self.notebook.add(self.tab_orig, text=" 📈 1. Original ")
        self.notebook.add(self.tab_base, text=" ⚖️ 2. Baseline ")
        self.notebook.add(self.tab_int,  text=" ⚡ 3. Integración ")
        self.notebook.add(self.tab_gal,  text=" 🔭 4. Galería Pulsos ")
        self.notebook.add(self.tab_hm,   text=" 📊 5. Hist. Muestra ")
        self.notebook.add(self.tab_hc,   text=" 🔍 6. Muestra vs Bkg ")
        self.notebook.add(self.tab_hn,   text=" 🧲 7. Muestra Neta ")

    # ════════════════════════════════════════════════════════════════
    #  CONTROLADORES DE EVENTOS (Conectan UI con Datos)
    # ════════════════════════════════════════════════════════════════
    def cargar_muestra(self):
        ruta = filedialog.askopenfilename(title="Seleccionar Muestra", filetypes=[("NPZ", "*.npz")])
        if not ruta: return
        try:
            num_p, elim = self.data.cargar_muestra(ruta)
            
            # Actualizar lista visual
            self.panel_control.lista_archivos.delete(0, tk.END)
            for i in range(num_p):
                self.panel_control.lista_archivos.insert(tk.END, f"Pulso {i+1}")
                
            self.aplicar_conversion()
            msg = f"Muestra: {num_p} pulsos únicos." + (f"\nDuplicados eliminados: {elim}" if elim else "")
            self.panel_control.lbl_status.config(text=msg)
            messagebox.showinfo("Éxito", msg)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def cargar_fondo(self):
        ruta = filedialog.askopenfilename(title="Seleccionar Background", filetypes=[("NPZ", "*.npz")])
        if not ruta: return
        try:
            num_p, elim = self.data.cargar_fondo(ruta)
            self.aplicar_conversion()
            msg = f"Background: {num_p} pulsos únicos." + (f"\nDuplicados eliminados: {elim}" if elim else "")
            messagebox.showinfo("Éxito", msg)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def aplicar_conversion(self):
        try:
            # Mandar datos crudos al DataManager
            self.data.aplicar_conversion(
                v_scale=self.panel_control.v_scale.get(),
                v_offset=self.panel_control.v_offset.get(),
                y_ref=self.panel_control.y_ref.get(),
                f_muestreo_ghz=self.panel_control.f_muestreo_ghz.get(),
                factor_escala_hw=self.panel_control.factor_escala_hw.get()
            )
            # Actualizar etiqueta de tiempo
            dt = self.data.get_dt_ns(self.panel_control.f_muestreo_ghz.get(), self.panel_control.factor_escala_hw.get())
            self.panel_control.lbl_dt.config(text=f"→ dt real ≈ {dt*1000:.3f} ps/muestra")
            
            # Si hay un pulso seleccionado, recalcular visualización
            if self.panel_control.lista_archivos.curselection():
                self.procesar_datos_individual()
        except Exception as e:
            messagebox.showwarning("Advertencia", f"Error en conversión: {e}")

    def on_select_archivo(self, event):
        sel = self.panel_control.lista_archivos.curselection()
        if sel:
            self.idx_pulso_actual = sel[0]
            self.procesar_datos_individual()

    def procesar_datos_individual(self):
        if self.data.pulsos_muestra is None: return
        
        try:
            v_actual, _ = self.data.get_pulso(self.idx_pulso_actual, fuente="Muestra")
            t_plot = self.data.t_ns if self.data.t_ns is not None else self.data.t_vector
            is_ns = self.data.t_ns is not None
            
            parametros = self.panel_control.get_parametros_integracion()
            
            # Calcular matemáticas puras
            area, offset, idx_i, idx_f, v_corr, v_teorico = processing.calcular_area(
                t_ns=t_plot, v=v_actual, **parametros
            )
            
            self.panel_control.area_calculada.set(round(area, 6))

            # Mandar a dibujar a las pestañas
            self.tab_orig.actualizar(t_plot, v_actual, offset, is_ns)
            self.tab_base.actualizar(t_plot, v_corr, offset, is_ns)
            self.tab_int.actualizar(t_plot, v_corr, area, idx_i, idx_f, parametros["umbral_amplitud"], is_ns,v_teorico)
            
        except Exception as e:
            print(f"Error previsualizando pulso: {e}")

    def procesar_archivos_completos(self):
        try:
            parametros = self.panel_control.get_parametros_integracion()
            self.panel_control.lbl_status.config(text="Procesando datos...")
            self.root.update_idletasks()
            
            self.data.integrar_todos(parametros)
            
            bins = self.panel_control.num_bins.get()
            self.tab_hm.actualizar(self.data.areas_acumuladas, bins)
            
            if self.data.areas_fondo:
                self.tab_hc.actualizar(self.data.areas_acumuladas, self.data.areas_fondo, bins)
                self.tab_hn.actualizar(self.data.areas_acumuladas, self.data.areas_fondo, bins)
                self.notebook.select(self.tab_hn)
            else:
                self.notebook.select(self.tab_hm)
                
            self.panel_control.lbl_status.config(text=f"✓ Muestra: {len(self.data.areas_acumuladas)} pulsos válidos.")
        except Exception as e:
            messagebox.showerror("Error de Procesamiento", str(e))
            self.panel_control.lbl_status.config(text="Error en procesamiento.")

    def exportar_csv(self):
        ruta = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not ruta: return
        try:
            self.data.exportar_csv(ruta)
            messagebox.showinfo("Exportado", f"Archivo guardado en:\n{ruta}")
        except Exception as e:
            messagebox.showerror("Error", str(e))


if __name__ == "__main__":
    root = tk.Tk()
    app = AnalizadorApp(root)
    root.mainloop()