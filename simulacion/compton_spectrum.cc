#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv)
{
    // Modo visual: sin argumentos o cuando se pasa explicitamente vis.mac.
    // Cualquier otro .mac se ejecuta en modo batch (sin UI).
    bool visualMode = (argc == 1) ||
                      (argc == 2 && G4String(argv[1]) == "vis.mac");
    G4UIExecutive* ui = visualMode ? new G4UIExecutive(argc, argv) : nullptr;

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    auto* det = new DetectorConstruction();
    runManager->SetUserInitialization(det);
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(det));

    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    auto* UImanager = G4UImanager::GetUIpointer();
    if (!ui) {
        UImanager->ApplyCommand("/control/execute " + G4String(argv[1]));
    } else {
        G4String macFile = (argc == 2) ? G4String(argv[1]) : "vis.mac";
        UImanager->ApplyCommand("/control/execute " + macFile);
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
    return 0;
}
