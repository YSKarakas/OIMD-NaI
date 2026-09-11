// scint_optical -- scintillator light-collection simulation.
//
// Everything about a run comes from the macro passed on the command line, so a
// run is fully described by that file plus the recorded environment. There is no
// interactive mode and no hidden default: a missing material specification is an
// error rather than a silent fallback.
//
//   scint_optical <macro>
//
// The run manager is deliberately serial. Parallelism in this project lives at
// the level of whole runs (scint.runner executes many configurations at once),
// which keeps each individual result bit-for-bit reproducible from its seed.

#include "Actions.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include <iostream>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: scint_optical <macro>\n";
    return 2;
  }

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  runManager->SetUserInitialization(new scint::DetectorConstruction());
  runManager->SetUserInitialization(new scint::PhysicsList());
  runManager->SetUserInitialization(new scint::ActionInitialization());

  G4UImanager* ui = G4UImanager::GetUIpointer();
  const G4int status = ui->ApplyCommand(G4String("/control/execute ") + argv[1]);

  delete runManager;
  return status == 0 ? 0 : 1;
}
