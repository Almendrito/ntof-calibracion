#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class G4Run;

// -----------------------------------------------------------------------------
// RunAction
//
// Define el ntuple y el histograma del espectro, abre/cierra el archivo ROOT y
// mezcla los ntuples de los hilos.
//
// El archivo se llama ComptonSpectrum_<fuente>_<N>.root, donde N son los
// primarios LANZADOS en la corrida (no las filas del ntuple, que solo cuenta los
// eventos con deposito > 0). Tenerlo en el nombre evita confundir corridas de
// distinta estadistica, que es facil porque el espectro se ve igual y solo
// cambia el ruido.
// -----------------------------------------------------------------------------
class RunAction : public G4UserRunAction
{
public:
    explicit RunAction(DetectorConstruction* det);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

private:
    static G4String FormatCount(G4int n);   // 3000000 -> "3M", 20000 -> "20k"
    G4String FileName(const G4Run* run) const;

    DetectorConstruction* fDet = nullptr;
};

#endif
