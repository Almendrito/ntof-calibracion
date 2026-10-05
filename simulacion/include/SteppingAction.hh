#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"

class EventAction;
class G4EmSaturation;

// -----------------------------------------------------------------------------
// SteppingAction
//
// En cada step dentro del volumen "Scintillator" suma al EventAction dos cosas:
//
//   1. la energia depositada (Edep), y
//   2. la energia VISIBLE tras el quenching de Birks (la luz, en MeVee),
//      calculada por G4EmSaturation step a step.
//
// Si el track es un gamma, ademas cuenta la interaccion (compt / phot).
// -----------------------------------------------------------------------------
class SteppingAction : public G4UserSteppingAction
{
public:
    explicit SteppingAction(EventAction* ea);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

private:
    EventAction*          fEventAction = nullptr;
    const G4EmSaturation* fSaturation  = nullptr;   // propiedad de G4LossTableManager
};

#endif
