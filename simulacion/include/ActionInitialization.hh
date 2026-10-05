#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

// -----------------------------------------------------------------------------
// ActionInitialization
//
// Registra RunAction en el master (para mezclar la salida) y, en cada hilo
// trabajador, el generador primario, EventAction, SteppingAction y RunAction.
// -----------------------------------------------------------------------------
class ActionInitialization : public G4VUserActionInitialization
{
public:
    explicit ActionInitialization(DetectorConstruction* det);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    DetectorConstruction* fDet = nullptr;
};

#endif
