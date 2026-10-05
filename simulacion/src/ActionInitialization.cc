#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* det)
: fDet(det)
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction(fDet));
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction(fDet));

    auto* eventAction = new EventAction();
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(eventAction));
    SetUserAction(new RunAction(fDet));
}
