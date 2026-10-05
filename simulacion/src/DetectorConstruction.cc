#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4IonisParamMat.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cstdio>

DetectorConstruction::DetectorConstruction()
: fScintRadius(2.54 * cm),   // 2" de diametro
  fScintHalfZ(2.54 * cm),    // 2" de largo (semieje = 2.54 cm)
  fBuildPMT(true),
  fBuildHousing(true),
  fHousingWall(10. * mm),    // pared lateral de Al medida
  fHousingFront(6. * mm),    // tapa frontal de Al (hacia la fuente)
  fBirks(0.126 * mm / MeV),   // BC408/PVT, valor de literatura (ver README)
  fCascade(false),
  fSourceDistance(10. * cm),
  fBeamType("pencil"),
  fSourceLabel("Co60"),
  fParticle("gamma")
{
    SetPreset("co60");   // llena fEnergies / fWeights por defecto
    DefineCommands();
}

DetectorConstruction::~DetectorConstruction()
{
    delete fGeoMessenger;
    delete fSrcMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();
    auto* air   = nist->FindOrBuildMaterial("G4_AIR");
    // Base PVT del BC408 (= BC400 / EJ-200): H 8.5% / C 91.5% en masa, 1.032 g/cm3.
    // Los fluors que distinguen al BC408 cambian la luz, no la fisica del deposito:
    // por eso el DEPOSITO se simula con el PVT desnudo y la LUZ sale del quenching
    // de Birks aplicado paso a paso (ver SteppingAction).
    auto* pvt   = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    ApplyBirks();
    auto* glass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");   // ventana de PMT (borosilicato)
    auto* alu   = nist->FindOrBuildMaterial("G4_Al");

    // --- Mundo ---
    G4double worldHalf = std::max(fScintRadius, fScintHalfZ) + 30. * cm;
    auto* worldSolid = new G4Box("World", worldHalf, worldHalf, worldHalf + fSourceDistance);
    auto* worldLV = new G4LogicalVolume(worldSolid, air, "World");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV,
                                      "World", nullptr, false, 0, true);

    // --- Centellador BC408 (PVT): cilindro centrado en el origen, eje z ---
    // El haz entra por la cara frontal en z = -fScintHalfZ.
    auto* scintSolid = new G4Tubs("Scintillator", 0., fScintRadius, fScintHalfZ, 0., twopi);
    auto* scintLV = new G4LogicalVolume(scintSolid, pvt, "Scintillator");
    auto* scintVis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.45));
    scintVis->SetForceSolid(true);
    scintLV->SetVisAttributes(scintVis);
    new G4PVPlacement(nullptr, G4ThreeVector(), scintLV,
                      "Scintillator", worldLV, false, 0, true);

    // --- Carcasa de Al opcional: pared lateral (cilindro) + tapa frontal ---
    // Grosores medidos: pared lateral fHousingWall, tapa frontal fHousingFront.
    // La tapa frontal queda entre la fuente y el centellador (atenua el haz que
    // entra y agrega dispersion/retrodispersion). La pared lateral cubre el largo
    // del centellador; la tapa cubre toda la cara frontal del bote (radio incluido
    // el grosor de la pared).
    if (fBuildHousing) {
        auto* alVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.3));

        // Pared lateral: tubo de fScintRadius a fScintRadius + fHousingWall.
        auto* wallSolid = new G4Tubs("Wall", fScintRadius, fScintRadius + fHousingWall,
                                     fScintHalfZ, 0., twopi);
        auto* wallLV = new G4LogicalVolume(wallSolid, alu, "Wall");
        new G4PVPlacement(nullptr, G4ThreeVector(), wallLV, "Wall", worldLV, false, 0, true);
        wallLV->SetVisAttributes(alVis);

        // Tapa frontal (hacia la fuente, en -z): disco de Al que cubre todo el frente.
        auto* capSolid = new G4Tubs("FrontCap", 0., fScintRadius + fHousingWall,
                                    0.5 * fHousingFront, 0., twopi);
        auto* capLV = new G4LogicalVolume(capSolid, alu, "FrontCap");
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -(fScintHalfZ + 0.5 * fHousingFront)),
                          capLV, "FrontCap", worldLV, false, 0, true);
        capLV->SetVisAttributes(alVis);
    }

    // --- Stub de PMT (vidrio borosilicato) detras del centellador (+z) ---
    // Principal contribuyente al pico de retrodispersion por detras del detector.
    if (fBuildPMT) {
        const G4double pmtHalfZ = 1.5 * cm;
        auto* pmtSolid = new G4Tubs("PMT", 0., fScintRadius, pmtHalfZ, 0., twopi);
        auto* pmtLV = new G4LogicalVolume(pmtSolid, glass, "PMT");
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fScintHalfZ + pmtHalfZ + 0.6 * mm),
                          pmtLV, "PMT", worldLV, false, 0, true);
        auto* pmtVis = new G4VisAttributes(G4Colour(0.85, 0.8, 0.2, 0.3));
        pmtVis->SetForceSolid(true);
        pmtLV->SetVisAttributes(pmtVis);
    }

    return worldPV;
}

// -----------------------------------------------------------------------------
// Presets de fuente. Energias y pesos relativos de los fotones por decaimiento.
// -----------------------------------------------------------------------------
void DetectorConstruction::SetPreset(const G4String& name)
{
    fEnergies.clear();
    fWeights.clear();
    if (name == "co60" || name == "Co60" || name == "Co-60") {
        fEnergies = {1173.2 * keV, 1332.5 * keV};   // bordes Compton: 963 y 1118 keV
        fWeights  = {1.0, 1.0};                       // ~1 gamma de cada uno por decaimiento
        fSourceLabel = "Co60";
    } else if (name == "cs137" || name == "Cs137" || name == "Cs-137") {
        fEnergies = {661.66 * keV};                   // borde Compton: 477 keV
        fWeights  = {1.0};
        fSourceLabel = "Cs137";
    } else if (name == "na22" || name == "Na22" || name == "Na-22") {
        fEnergies = {511.0 * keV, 1274.5 * keV};      // bordes: 341 y 1062 keV
        fWeights  = {2.0, 1.0};                        // dos 511 (aniquilacion) por cada 1274.5
        fSourceLabel = "Na22";
    } else {
        G4cerr << "[DetectorConstruction] preset desconocido: '" << name
               << "'. Use co60 | cs137 | na22." << G4endl;
    }
}

void DetectorConstruction::SetBirks(G4double kB_mm_per_MeV)
{
    fBirks = kB_mm_per_MeV * mm / MeV;
}

void DetectorConstruction::SetSingleEnergy(G4double e)
{
    fEnergies = {e};
    fWeights  = {1.0};
    char buf[32];
    std::snprintf(buf, sizeof(buf), "E%.0fkeV", e / keV);
    fSourceLabel = buf;
}

// Fija la constante de Birks en el material del centellador.
//
// G4EmSaturation la lee del material via GetIonisation()->GetBirksConstant() y,
// si vale 0, no aplica quenching alguno. Su tabla interna NO incluye
// G4_PLASTIC_SC_VINYLTOLUENE (esta comentado en G4EmSaturation.cc), asi que sin
// esta llamada la luz simulada seria exactamente proporcional a la energia
// depositada y el efecto que queremos estudiar desapareceria.
void DetectorConstruction::ApplyBirks()
{
    auto* pvt = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    if (!pvt) {
        G4cerr << "[DetectorConstruction] no se encontro el material del centellador; "
               << "no se pudo fijar la constante de Birks." << G4endl;
        return;
    }
    pvt->GetIonisation()->SetBirksConstant(fBirks);
    G4cout << "[DetectorConstruction] BC408: kB = " << fBirks * MeV / mm << " mm/MeV"
           << "  (= " << fBirks * MeV / mm * 0.1 * pvt->GetDensity() / (g / cm3)
           << " g/cm2/MeV)";
    if (fBirks <= 0.) G4cout << "   <-- kB = 0: SIN quenching, luz proporcional a Edep";
    G4cout << G4endl;
}

void DetectorConstruction::UpdateGeometry()
{
    // Solo valido despues de /run/initialize (ver Guia de Geant4).
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::DefineCommands()
{
    // --- Geometria: /det/ ---
    fGeoMessenger = new G4GenericMessenger(this, "/det/", "Geometria y luz del detector BC408");
    fGeoMessenger->DeclarePropertyWithUnit("radius", "cm", fScintRadius,
        "Radio del centellador (fijar antes de /run/initialize)");
    fGeoMessenger->DeclarePropertyWithUnit("halfZ", "cm", fScintHalfZ,
        "Semilargo del centellador, mitad del largo (antes de /run/initialize)");
    fGeoMessenger->DeclareProperty("pmt", fBuildPMT,
        "Construir el stub de PMT (true/false)");
    fGeoMessenger->DeclareProperty("housing", fBuildHousing,
        "Construir la carcasa de Al (true/false)");
    fGeoMessenger->DeclarePropertyWithUnit("housingWall", "mm", fHousingWall,
        "Grosor de la pared lateral de Al (antes de /run/initialize)");
    fGeoMessenger->DeclarePropertyWithUnit("housingFront", "mm", fHousingFront,
        "Grosor de la tapa frontal de Al hacia la fuente (antes de /run/initialize)");
    // kB se pasa como numero desnudo en mm/MeV: Geant4 no tiene esa unidad en su
    // tabla (solo la inversa, MeV/cm), asi que DeclarePropertyWithUnit no sirve.
    fGeoMessenger->DeclareMethod("birks", &DetectorConstruction::SetBirks,
        "Constante de Birks kB en mm/MeV (antes de /run/initialize). 0 = sin quenching");
    fGeoMessenger->DeclareMethod("update", &DetectorConstruction::UpdateGeometry,
        "Reconstruir la geometria (solo despues de /run/initialize)");

    // --- Fuente: /gun/ ---
    fSrcMessenger = new G4GenericMessenger(this, "/gun/", "Fuente gamma");
    fSrcMessenger->DeclareMethod("preset", &DetectorConstruction::SetPreset,
        "Preset de fuente: co60 | cs137 | na22");
    fSrcMessenger->DeclareMethodWithUnit("energy", "keV", &DetectorConstruction::SetSingleEnergy,
        "Fijar una sola linea gamma de la energia dada (reemplaza la lista)");
    fSrcMessenger->DeclareProperty("cascade", fCascade,
        "true: emite todas las lineas por evento (cascada/summing); false: muestrea una");
    fSrcMessenger->DeclarePropertyWithUnit("distance", "cm", fSourceDistance,
        "Distancia de la fuente puntual a la cara frontal del centellador");
    fSrcMessenger->DeclareProperty("beam", fBeamType,
        "Emision: pencil (haz axial) | cone (cono que cubre el detector) | iso (4pi)");
    fSrcMessenger->DeclareProperty("particle", fParticle,
        "gamma: fuente externa (espectro Compton) | e-: electrones desde el centro "
        "del centellador, para medir la curva de respuesta de luz L(E)");
    fSrcMessenger->DeclareProperty("label", fSourceLabel,
        "Etiqueta para el nombre del archivo de salida");
}
