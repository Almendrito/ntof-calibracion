#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

// -----------------------------------------------------------------------------
// EventAction
//
// Acumula por evento, sobre todos los steps dentro del volumen "Scintillator":
//
//   fEdep : energia depositada total.
//   fEvis : energia visible tras el quenching de Birks, es decir la LUZ emitida
//           en MeVee (la suma de los quench de cada step, no el quench de la
//           suma). Multiplicada por los fotones/MeVee del BC408 da el numero de
//           fotones de centelleo; como ese factor es constante, MeVee ya es la
//           luz en las unidades utiles para calibrar contra la carga del DRS4.
//
// mas los contadores de diagnostico. Al final del evento, si hubo deposito,
// llena una fila del ntuple y los histogramas.
// -----------------------------------------------------------------------------
class EventAction : public G4UserEventAction
{
public:
    EventAction() = default;
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddEdep(G4double e) { fEdep += e; }
    void AddEvis(G4double e) { fEvis += e; }
    void AddCompton()        { ++fNCompton; }
    void AddPhoto()          { ++fNPhoto; }

private:
    G4double fEdep     = 0.;
    G4double fEvis     = 0.;
    G4int    fNCompton = 0;
    G4int    fNPhoto   = 0;
};

#endif
