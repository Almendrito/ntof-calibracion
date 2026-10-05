#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4Gamma.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

SteppingAction::SteppingAction(EventAction* ea)
: fEventAction(ea),
  fSaturation(G4LossTableManager::Instance()->EmSaturation())
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const auto* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!preVol || preVol->GetName() != "Scintillator") return;

    // Energia depositada en el centellador (cualquier particula del step).
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.) {
        fEventAction->AddEdep(edep);

        // Luz emitida, via quenching de Birks:  dL/dx = (dE/dx) / (1 + kB dE/dx).
        // Se evalua PASO A PASO, nunca sobre el Edep total del evento: el
        // quenching es no lineal en dE/dx, asi que en un evento de scatter
        // multiple (varios electrones de baja energia) quenchar la suma da un
        // resultado distinto -y erroneo- frente a sumar los quench individuales.
        // La energia visible sale en MeV "equivalentes a electron" (MeVee): es la
        // luz salvo el factor constante fotones/MeVee del BC408.
        fEventAction->AddEvis(fSaturation->VisibleEnergyDepositionAtAStep(step));
    }

    // Diagnostico: contar las interacciones del foton en el centellador.
    // Con EM option4 el proceso sigue llamandose "compt"/"phot" (solo cambia el modelo).
    if (step->GetTrack()->GetParticleDefinition() == G4Gamma::Definition()) {
        const auto* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (proc) {
            const G4String& pname = proc->GetProcessName();
            if      (pname == "compt") fEventAction->AddCompton();
            else if (pname == "phot")  fEventAction->AddPhoto();
        }
    }
}
