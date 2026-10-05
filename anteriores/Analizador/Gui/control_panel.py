import tkinter as tk
from tkinter import ttk
from Gui.custom_widgets import ScrollableFrame

class ControlPanel(ScrollableFrame):
    """
    Panel izquierdo de la interfaz que contiene todos los controles,
    entradas de parámetros y la lista de pulsos.
    """
    def __init__(self, parent, controller, width=270):
        super().__init__(parent, width=width)
        self.controller = controller  # Referencia a la clase principal para ejecutar acciones
        
        # ── Variables de control de Tkinter ──────────────────────────
        self.baseline_muestras  = tk.IntVar(value=200)
        self.polaridad          = tk.StringVar(value="Negativo")
        self.area_calculada     = tk.DoubleVar(value=0.0)
        self.num_bins           = tk.IntVar(value=100)

        # Calibración
        self.v_scale  = tk.DoubleVar(value=0.200)
        self.v_offset = tk.DoubleVar(value=0.0)
        self.y_ref    = tk.IntVar(value=127)

        # Escala temporal
        self.f_muestreo_ghz   = tk.DoubleVar(value=1.0)
        self.factor_escala_hw = tk.DoubleVar(value=10.0)

        # Integración
        self.metodo_integracion = tk.StringVar(value="ventana_fija")
        self.pre_pico_muestras  = tk.IntVar(value=20)
        self.post_pico_muestras = tk.IntVar(value=80)
        self.umbral_inicio      = tk.DoubleVar(value=0.02)
        self.umbral_fin         = tk.DoubleVar(value=0.005)
        self.umbral_amplitud    = tk.DoubleVar(value=0.01)

        # Visualización
        self.viz_n_superponer   = tk.IntVar(value=20)
        self.viz_fuente         = tk.StringVar(value="Muestra")
        
        self._construir_ui()

    def _sep(self):
        ttk.Separator(self.inner, orient="horizontal").pack(fill=tk.X, pady=6)

    def _construir_ui(self):
        pc = self.inner # Alias para acortar código
        
        # ── Título ────────────────────────────────────────────────────
        ttk.Label(pc, text="⚛ BC404 / XP2262B", font=("Courier", 12, "bold")).pack(pady=(8, 2))
        ttk.Label(pc, text="Espectroscopía Compton", font=("Courier", 9), foreground="gray").pack()
        self._sep()

        # ── Carga de archivos ─────────────────────────────────────────
        ttk.Label(pc, text="DATOS", font=("Arial", 9, "bold")).pack(anchor="w", padx=5)
        ttk.Button(pc, text="📂 Cargar Muestra (.npz)", 
                   command=self.controller.cargar_muestra).pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(pc, text="📂 Cargar Background (.npz)", 
                   command=self.controller.cargar_fondo).pack(fill=tk.X, padx=5, pady=2)
        self._sep()

        # ── Calibración ADC → V ───────────────────────────────────────
        f_cal = ttk.LabelFrame(pc, text="Calibración Hardware (ADC → V)", padding=5)
        f_cal.pack(fill=tk.X, padx=5, pady=3)
        for lbl, var in [("Escala (V/div):", self.v_scale),
                         ("Offset Vertical (V):", self.v_offset),
                         ("Y Reference (ADC):", self.y_ref)]:
            ttk.Label(f_cal, text=lbl).pack(anchor="w")
            ttk.Entry(f_cal, textvariable=var, width=12).pack(anchor="w")

        # ── Escala Temporal ───────────────────────────────────────────
        f_time = ttk.LabelFrame(pc, text="Escala Temporal", padding=5)
        f_time.pack(fill=tk.X, padx=5, pady=3)
        ttk.Label(f_time, text="f_muestreo nominal (GHz):").pack(anchor="w")
        ttk.Entry(f_time, textvariable=self.f_muestreo_ghz, width=12).pack(anchor="w")
        ttk.Label(f_time, text="Factor HW (×10 por defecto):").pack(anchor="w")
        ttk.Entry(f_time, textvariable=self.factor_escala_hw, width=12).pack(anchor="w")
        
        self.lbl_dt = ttk.Label(f_time, text="→ dt real ≈ -- ps/muestra", foreground="#aaaaaa")
        self.lbl_dt.pack(anchor="w")
        
        ttk.Button(f_time, text="⚡ Aplicar Calibración", 
                   command=self.controller.aplicar_conversion).pack(pady=(6, 2), fill=tk.X)
        self._sep()

        # ── Ajustes de análisis ───────────────────────────────────────
        f_aj = ttk.LabelFrame(pc, text="Ajustes de Análisis", padding=5)
        f_aj.pack(fill=tk.X, padx=5, pady=3)
        ttk.Label(f_aj, text="Baseline (muestras):").pack(anchor="w")
        ttk.Entry(f_aj, textvariable=self.baseline_muestras, width=12).pack(anchor="w")
        ttk.Label(f_aj, text="Umbral amplitud mín (V):").pack(anchor="w", pady=(5, 0))
        ttk.Entry(f_aj, textvariable=self.umbral_amplitud, width=12).pack(anchor="w")
        ttk.Label(f_aj, text="Polaridad:").pack(anchor="w", pady=(5, 0))
        ttk.OptionMenu(f_aj, self.polaridad, "Negativo", "Negativo", "Positivo").pack(anchor="w")

        # ── Método de integración ─────────────────────────────────────
        f_met = ttk.LabelFrame(pc, text="Método de Integración", padding=5)
        f_met.pack(fill=tk.X, padx=5, pady=3)
        
        ttk.Radiobutton(f_met, text="Ventana Fija (anti-afterpulse)",
                        variable=self.metodo_integracion, value="ventana_fija",
                        command=self._toggle_metodo).pack(anchor="w")
        ttk.Radiobutton(f_met, text="Umbral clásico",
                        variable=self.metodo_integracion, value="umbral",
                        command=self._toggle_metodo).pack(anchor="w")
        
        ttk.Radiobutton(f_met, text="Ajuste Exponencial (Automático)",
                        variable=self.metodo_integracion, value="ajuste_exponencial",
                        command=self._toggle_metodo).pack(anchor="w")

        self.f_vf = ttk.Frame(f_met)
        ttk.Label(self.f_vf, text="Muestras pre-pico:").pack(anchor="w")
        ttk.Entry(self.f_vf, textvariable=self.pre_pico_muestras, width=12).pack(anchor="w")
        ttk.Label(self.f_vf, text="Muestras post-pico:").pack(anchor="w", pady=(4, 0))
        ttk.Entry(self.f_vf, textvariable=self.post_pico_muestras, width=12).pack(anchor="w")

        self.f_umb = ttk.Frame(f_met)
        ttk.Label(self.f_umb, text="Umbral inicio (V):").pack(anchor="w")
        ttk.Entry(self.f_umb, textvariable=self.umbral_inicio, width=12).pack(anchor="w")
        ttk.Label(self.f_umb, text="Umbral cierre (V):").pack(anchor="w", pady=(4, 0))
        ttk.Entry(self.f_umb, textvariable=self.umbral_fin, width=12).pack(anchor="w")
        
        self._toggle_metodo() # Mostrar el frame correcto al inicio
        self._sep()

        # ── Galería y Preview ─────────────────────────────────────────
        # (Para mantener el código corto, aquí irían los botones de galería
        # que apuntan a self.controller.vista_superponer_n, etc.)
        f_prev = ttk.LabelFrame(pc, text="Preview: Efecto del Filtro", padding=5)
        f_prev.pack(fill=tk.X, padx=5, pady=3)
        ttk.Button(f_prev, text="👁 Previsualizar Integración", 
                   command=self.controller.procesar_datos_individual).pack(fill=tk.X, pady=(5, 2))
        
        f_area = ttk.Frame(f_prev)
        f_area.pack(fill=tk.X, pady=3)
        ttk.Label(f_area, text="Área actual (V·ns):").pack(side=tk.LEFT)
        ttk.Label(f_area, textvariable=self.area_calculada, 
                  foreground="#00ff88", font=("Courier", 9, "bold")).pack(side=tk.LEFT, padx=4)
        self._sep()

        # ── Acciones principales ──────────────────────────────────────
        ttk.Button(pc, text="📊 INTEGRAR TODOS", 
                   command=self.controller.procesar_archivos_completos).pack(fill=tk.X, padx=5, pady=4)
        ttk.Button(pc, text="💾 Exportar Áreas (.csv)", 
                   command=self.controller.exportar_csv).pack(fill=tk.X, padx=5, pady=2)
        self._sep()

        # ── Lista de pulsos ───────────────────────────────────────────
        ttk.Label(pc, text="LISTA DE PULSOS", font=("Arial", 9, "bold")).pack(anchor="w", padx=5)
        self.lista_archivos = tk.Listbox(pc, height=12, width=34, background="#1a1a1a", 
                                         foreground="white", selectbackground="#005577")
        self.lista_archivos.pack(fill=tk.X, padx=5, pady=3)
        self.lista_archivos.bind('<<ListboxSelect>>', self.controller.on_select_archivo)

        self.lbl_status = ttk.Label(pc, text="Sin datos", foreground="gray", wraplength=240)
        self.lbl_status.pack(padx=5, pady=5, anchor="w")

    def _toggle_metodo(self):
        if self.metodo_integracion.get() == "ventana_fija":
            self.f_vf.pack(fill=tk.X, pady=2)
            self.f_umb.pack_forget()
        else:
            self.f_umb.pack(fill=tk.X, pady=2)
            self.f_vf.pack_forget()
            
    def get_parametros_integracion(self):
        """Retorna un diccionario limpio con los parámetros para pasarlo al DataManager"""
        return {
            "baseline_muestras": self.baseline_muestras.get(),
            "polaridad": self.polaridad.get(),
            "umbral_amplitud": self.umbral_amplitud.get(),
            "metodo_integracion": self.metodo_integracion.get(),
            "pre_pico": self.pre_pico_muestras.get(),
            "post_pico": self.post_pico_muestras.get(),
            "umbral_inicio": self.umbral_inicio.get(),
            "umbral_fin": self.umbral_fin.get()
        }