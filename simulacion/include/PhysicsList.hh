#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"

// -----------------------------------------------------------------------------
// PhysicsList
//
// EM "option4" (la mejor fisica EM de Geant4). Para gammas < 20 MeV usa
// G4LowEPComptonModel, que incluye el perfil Doppler de los electrones atomicos:
// esto ensancha intrinsecamente el borde Compton incluso antes de la resolucion
// del detector. Se agrega G4DecayPhysics por completitud (inocuo para gammas).
// -----------------------------------------------------------------------------
class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;
};

#endif
