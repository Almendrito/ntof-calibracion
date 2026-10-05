"""
replica_paper — replica LITERAL de Kong et al., IEEE TNS 2025 (arXiv:2507.16282).

Regla de esta carpeta: aca solo entra lo que el paper dice. Nada de plantillas
Monte Carlo, modelo extendido, tiempo vivo como parametro, MCMC ni sistematicos
medidos: todo eso vive en `codigo/` y `analisis/` y es EXTRA respecto del paper.

Lo unico que se toma prestado de afuera es el paso XML -> (canal, cuentas), que
es justamente donde ARRANCA el diagrama de flujo del paper (Fig. 2: "Energy
Spectrum (Channel, Counts)"). Ese prestamo esta aislado en `datos.py` y marcado.
"""
