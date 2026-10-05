#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
{
    RegisterPhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4DecayPhysics());

    // Corte de produccion fino: el rango CSDA de los electrones de retroceso de
    // ~1 MeV en plastico es de unos pocos mm; un corte de 0.1 mm sigue bien las
    // deltas cerca de la superficie sin penalizar demasiado el tiempo de CPU.
    SetDefaultCutValue(0.1 * mm);
}
