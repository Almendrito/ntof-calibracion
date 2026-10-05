import tkinter as tk
from tkinter import ttk

class ScrollableFrame(ttk.Frame):
    """
    Frame genérico con scroll vertical interno.
    Ideal para paneles de control con muchas opciones.
    """
    def __init__(self, parent, width=260, **kwargs):
        super().__init__(parent, **kwargs)
        self._canvas = tk.Canvas(self, width=width, highlightthickness=0, background="#2b2b2b")
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=self._canvas.yview)
        self.inner = ttk.Frame(self._canvas)

        self.inner.bind(
            "<Configure>",
            lambda e: self._canvas.configure(scrollregion=self._canvas.bbox("all"))
        )
        self._canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self._canvas.configure(yscrollcommand=scrollbar.set)

        self._canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Scroll con rueda del ratón
        self._canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self._canvas.bind_all("<Button-4>", self._on_mousewheel)
        self._canvas.bind_all("<Button-5>", self._on_mousewheel)

    def _on_mousewheel(self, event):
        # Manejo compatible con Windows/Linux/macOS
        if event.num == 4:
            self._canvas.yview_scroll(-1, "units")
        elif event.num == 5:
            self._canvas.yview_scroll(1, "units")
        else:
            self._canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")