#include <TPrimary.h>

#include "Physics/Hadronization/PythiaSingleton.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

ClassImp(PythiaSingleton)

PythiaSingleton * PythiaSingleton::fgInstance = 0;

//____________________________________________________________________________
PythiaSingleton::PythiaSingleton()
{
    // Constructor
    if (fgInstance) {
      // TODO raise assertion?
      LOG("PythiaSingleton", pERROR) <<
          "Instance of PythiaSingleton already exists";
      return;
    }

    fPythia = new Pythia8::Pythia();
}
//____________________________________________________________________________
PythiaSingleton::~PythiaSingleton()
{
    // Destructor
    if (fgInstance) {
        fgInstance->Delete();
        delete fgInstance;
        fgInstance = 0;
    }
    delete fPythia;
}
//____________________________________________________________________________
PythiaSingleton* PythiaSingleton::Instance() 
{
    return fgInstance ? fgInstance : (fgInstance = new PythiaSingleton());
}

