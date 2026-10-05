#!/usr/bin/env bash
# Lanza la visualizacion en WSL2. GALLIUM_DRIVER=d3d12 evita el driver ZINK de Mesa.
source /home/mateo/Software/geant4-install/bin/geant4.sh
GALLIUM_DRIVER=d3d12 ./compton_spectrum vis.mac
