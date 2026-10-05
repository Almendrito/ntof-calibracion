#include "RunAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

#include <cstdio>

RunAction::RunAction(DetectorConstruction* det)
: fDet(det)
{
    auto* am = G4AnalysisManager::Instance();
    am->SetDefaultFileType("root");
    am->SetNtupleMerging(true);   // mezcla los ntuples de todos los hilos en un archivo
    am->SetVerboseLevel(1);

    // -------------------------------------------------------------------------
    // Ntuple: una fila por evento con deposito > 0 en el centellador.
    //
    // Edep_keV  : energia depositada en el centellador (la observable clave).
    // E0_keV    : energia primaria total del evento (linea muestreada, o suma en cascada).
    // NCompton  : numero de dispersiones Compton del foton en el centellador.
    // NPhoto    : numero de absorciones fotoelectricas (~0 en plastico de Z bajo).
    // FullAbs   : 1 si Edep ~ E0 (1%), absorcion casi total del evento.
    // Lvis_keVee: LUZ emitida = energia visible tras el quenching de Birks,
    //             sumada step a step. Es la observable que hay que comparar con
    //             la carga del DRS4: el PMT ve luz, no energia.
    // QF        : factor de quenching del evento, Lvis/Edep (<= 1). Cuanto se
    //             aparta de 1 es la no linealidad intrinseca del centellador.
    // -------------------------------------------------------------------------
    am->CreateNtuple("ComptonSpectrum", "Espectro Compton en BC408");
    am->CreateNtupleIColumn("EventID");    // 0
    am->CreateNtupleDColumn("Edep_keV");   // 1
    am->CreateNtupleDColumn("E0_keV");     // 2
    am->CreateNtupleIColumn("NCompton");   // 3
    am->CreateNtupleIColumn("NPhoto");     // 4
    am->CreateNtupleIColumn("FullAbs");    // 5
    am->CreateNtupleDColumn("Lvis_keVee"); // 6
    am->CreateNtupleDColumn("QF");         // 7
    am->FinishNtuple(0);

    // Histograma de conveniencia (espectro "verdadero", sin resolucion).
    // 0..2600 keV cubre incluso la suma en cascada de Co-60 (2505 keV).
    am->CreateH1("Edep", "Energia depositada;E_{dep} [keV];Cuentas", 1300, 0., 2600.);
    am->CreateH1("Lvis", "Luz emitida (Birks);L [keVee];Cuentas", 1300, 0., 2600.);
}

// Compacta el numero de eventos para que quepa en el nombre del archivo:
// 3000000 -> "3M", 2500000 -> "2.5M", 20000 -> "20k", 500 -> "500".
G4String RunAction::FormatCount(G4int n)
{
    char buf[32];
    if (n >= 1000000 && n % 1000000 == 0)      std::snprintf(buf, sizeof(buf), "%dM", n / 1000000);
    else if (n >= 1000000 && n % 100000 == 0)  std::snprintf(buf, sizeof(buf), "%.1fM", n / 1.0e6);
    else if (n >= 1000 && n % 1000 == 0)       std::snprintf(buf, sizeof(buf), "%dk", n / 1000);
    else                                        std::snprintf(buf, sizeof(buf), "%d", n);
    return {buf};
}

// Nombre del archivo de salida. DEBE dar exactamente lo mismo en el master y en
// todos los workers: si difiriera, cada hilo abriria un archivo distinto y el
// merge del ntuple se romperia. Por eso se usa GetNumberOfEventToBeProcessed()
// (el total de la corrida, que el worker hereda del master) y no el numero de
// eventos que a ese hilo le toco procesar.
G4String RunAction::FileName(const G4Run* run) const
{
    G4String name = "ComptonSpectrum_" + fDet->GetSourceLabel();
    const G4int n = run ? run->GetNumberOfEventToBeProcessed() : 0;
    if (n > 0) name += "_" + FormatCount(n);
    return name + ".root";
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    // Fuerza la inicializacion de G4EmSaturation en el master, ANTES de que los
    // hilos worker den un solo step. Su constructor solo se inicializa del todo
    // si la tabla de materiales ya existe; hacerlo aqui (geometria ya construida,
    // workers aun sin arrancar) evita depender de ese orden y evita que varios
    // hilos redimensionen sus vectores estaticos a la vez.
    if (IsMaster()) {
        auto* sat = G4LossTableManager::Instance()->EmSaturation();
        sat->InitialiseG4Saturation();
        sat->DumpBirksCoefficients();
    }

    auto* am = G4AnalysisManager::Instance();
    am->OpenFile(FileName(run));
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    auto* am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();

    if (IsMaster()) {
        G4cout << "\n====================================================="
               << "\n  Espectro Compton BC408 - simulacion terminada"
               << "\n  Fuente        : " << fDet->GetSourceLabel()
               << "\n  Primarios     : " << run->GetNumberOfEvent()
               << "\n  Archivo ROOT  : " << FileName(run)
               << "\n=====================================================\n"
               << G4endl;
    }
}
