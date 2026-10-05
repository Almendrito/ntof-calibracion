#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

// -----------------------------------------------------------------------------
// PrimaryGeneratorAction
//
// Dispara gammas desde una fuente puntual en el eje z, leyendo la configuracion
// (energias, pesos, cascada, distancia, tipo de haz) desde DetectorConstruction.
//   - cascade=false: muestrea una linea por evento segun los pesos.
//   - cascade=true : emite todas las lineas en el mismo evento (cascada real,
//                    p.ej. los dos gammas de Co-60), produciendo summing.
// -----------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    explicit PrimaryGeneratorAction(DetectorConstruction* det);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    G4ThreeVector SampleDirection(const G4String& beamType, G4double srcZ) const;
    void FireGamma(G4Event* event, G4double energy,
                   const G4ThreeVector& pos, const G4ThreeVector& dir);

    G4ParticleGun* fGun = nullptr;
    DetectorConstruction* fDet = nullptr;
};

#endif
