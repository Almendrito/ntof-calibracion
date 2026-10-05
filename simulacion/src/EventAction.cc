#include "EventAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

void EventAction::BeginOfEventAction(const G4Event*)
{
    fEdep     = 0.;
    fEvis     = 0.;
    fNCompton = 0;
    fNPhoto   = 0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    if (fEdep <= 0.) return;   // solo eventos con deposito en el centellador

    // Energia primaria total del evento (suma sobre vertices/particulas).
    G4double e0 = 0.;
    for (G4int iv = 0; iv < event->GetNumberOfPrimaryVertex(); ++iv) {
        const auto* v = event->GetPrimaryVertex(iv);
        for (G4int ip = 0; ip < v->GetNumberOfParticle(); ++ip)
            e0 += v->GetPrimary(ip)->GetKineticEnergy();
    }
    const G4int fullAbs = (e0 > 0. && std::fabs(fEdep - e0) < 0.01 * e0) ? 1 : 0;

    auto* am = G4AnalysisManager::Instance();
    am->FillNtupleIColumn(0, event->GetEventID());
    am->FillNtupleDColumn(1, fEdep / keV);
    am->FillNtupleDColumn(2, e0 / keV);
    am->FillNtupleIColumn(3, fNCompton);
    am->FillNtupleIColumn(4, fNPhoto);
    am->FillNtupleIColumn(5, fullAbs);
    am->FillNtupleDColumn(6, fEvis / keV);
    am->FillNtupleDColumn(7, (fEdep > 0.) ? fEvis / fEdep : 0.);
    am->AddNtupleRow(0);

    am->FillH1(0, fEdep / keV);
    am->FillH1(1, fEvis / keV);
}
