#include <TPrimary.h>

#include "Physics/Hadronization/PythiaSingleton.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

using namespace Pythia8;

ClassImp(PythiaSingleton)

PythiaSingleton * PythiaSingleton::fgInstance = 0;

//____________________________________________________________________________
Pythia8::LHAup_Genie::LHAup_Genie()
{
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillInit(int idA, int idB, double eA, double eB)
{
    fillBeam("A", idA, eA);
    fillBeam("B", idB, eB);
    fillProcess();
    // Done.
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::clearInit()
{
    fillBeam("A", 0, -1.);
    fillBeam("B", 0, -1.);
    fillProcess();
    // Done.
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::setInit()
{
    // Call the routine that does the job.
    if (!fillHepRup()) return false;
    // Store beam and strategy info.
    setBeamA(idBeamAStore, eBeamAStore, pdfgAStore, pdfsAStore);
    setBeamB(idBeamBStore, eBeamBStore, pdfgBStore, pdfsBStore);
    setStrategy(-4);
    // Store process info. Protect against vanishing cross section.
    double xsec = max( 1e-10, xsecStore);
    addProcess( idProcStore, xsec, xsecerrStore, xwgtmaxStore);
    // Done.
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::setEvent(int idProcIn)
{
    // Call the routine that does the job.
    if (!fillHepEup()) return false;
    // Store process info.
    setProcess(idProcIn, xwgtStore, scaleStore, aqedStore, aqcdStore);
    // Store particle info.
    for (int ip = 0; ip < nParticlesStore; ++ip) {
        addParticle(idPartStore[ip], statusPartStore[ip], mother1PartStore[ip],
          mother2PartStore[ip], colPartStore[ip], acolPartStore[ip],
          pPartStore[ip][0], pPartStore[ip][1], pPartStore[ip][2], pPartStore[ip][3],
          pPartStore[ip][4], vtimPartStore[ip], spinPartStore[ip], scaleStore);
    }
    // Store x values (here E = pup_out[ip][3]), but note incomplete info.
    setPdf( idPartStore[0], idPartStore[1], pPartStore[0][3]/eBeamAStore,
        pPartStore[1][3]/eBeamBStore, 0., 0., 0., false);
    // Done.
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillBeam(string beam, int id, double e, int pdfg, int pdfs)
{
    if (beam == "A") {
        idBeamAStore = id; 
        eBeamAStore  = e;
        pdfgAStore   = pdfg;
        pdfsAStore   = pdfs;
    } else {
        idBeamBStore = id; 
        eBeamBStore  = e;
        pdfgBStore   = pdfg;
        pdfsBStore   = pdfs;
    }
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillProcess(int idProcIn, double xsec, double xsecerr,
    double xwgtmax)
{
    idProcStore  = idProcIn;    
    xsecStore    = xsec;
    xsecerrStore = xsecerr;
    xwgtmaxStore = xwgtmax;
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillEventInfo( int nParticles, double xwgt, double scale,
    double aqed, double aqcd)
{
    nParticlesStore = nParticles;
    xwgtStore = xwgt;
    scaleStore = scale;
    aqedStore = aqed;
    aqcdStore = aqcd;
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillNewParticle(int idPart, int statusPart, vector<double> pPart,
    int mother1Part, int mother2Part, int colPart, int acolPart, double vtimPart,
    double spinPart)
{
    idPartStore.push_back(idPart);
    statusPartStore.push_back(statusPart);
    if (mother1Part == 0 && statusPart>0) mother1Part = 1;
    mother1PartStore.push_back(mother1Part);
    if (mother2Part == 0 && statusPart>0) mother2Part = 2;
    mother2PartStore.push_back(mother2Part);
    if (colPart  == 0 && abs(idPart) < 10 && idPart > 0) colPart = 1;
    colPartStore.push_back(colPart);
    if (acolPart == 0 && abs(idPart) < 10 && idPart < 0) acolPart = 1;
    acolPartStore.push_back(acolPart);
    vtimPartStore.push_back(vtimPart);
    spinPartStore.push_back(spinPart);
    pPartStore.push_back(pPart);
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::clearEvent()
{
    nParticlesStore = 0;
    xwgtStore = scaleStore = aqedStore = aqcdStore = 0.; 
    idPartStore.clear();
    statusPartStore.clear();
    mother1PartStore.clear();
    mother2PartStore.clear();
    colPartStore.clear();
    acolPartStore.clear();
    vtimPartStore.clear();
    spinPartStore.clear();
    pPartStore.clear();
    return true;
}
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
    fEventReader = new Pythia8::LHAup_Genie();
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
    delete fEventReader;
}
//____________________________________________________________________________
PythiaSingleton* PythiaSingleton::Instance() 
{
    return fgInstance ? fgInstance : (fgInstance = new PythiaSingleton());
}

