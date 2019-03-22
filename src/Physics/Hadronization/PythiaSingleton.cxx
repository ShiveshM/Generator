#include <TPrimary.h>

#include "Physics/Hadronization/PythiaSingleton.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

using namespace std;

ClassImp(PythiaSingleton)

PythiaSingleton * PythiaSingleton::fgInstance = 0;

//____________________________________________________________________________
// Set the two beam momentum deviations and the beam vertex.
// Note that momenta are in units of GeV and vertices in mm,
// always with c = 1, so that e.g. time is in mm/c.
void Pythia8::GBeamShape::pick() {
  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Set beam A transverse momentum deviation by a two-dimensional Gaussian.
  if (allowMomentumSpread) {
    //    deltaPzA = -9.999;
    deltaPzA = sigmaPzA - eMax;
  } 
}
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
bool Pythia8::LHAup_Genie::fillEventInfo(int nParticles, double xwgt,
    double scale, double aqed, double aqcd)
{
    nParticlesStore = nParticles;
    xwgtStore = xwgt;
    scaleStore = scale;
    aqedStore = aqed;
    aqcdStore = aqcd;
    return true;
}
//____________________________________________________________________________
bool Pythia8::LHAup_Genie::fillNewParticle(int idPart, int statusPart,
    Pythia8::Vec4 pPart, int mother1Part, int mother2Part, int colPart,
    int acolPart, double vtimPart, double spinPart)
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
    vector<double> temp;
    temp.push_back(pPart.px());
    temp.push_back(pPart.py());
    temp.push_back(pPart.pz());
    temp.push_back(pPart.e());
    temp.push_back(pPart.mCalc());
    pPartStore.push_back(temp);
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
      LOG("PythiaSingleton", pERROR) <<
          "Instance of PythiaSingleton already exists";
      return;
    }

    fPythia = new Pythia8::Pythia();
}
//____________________________________________________________________________
PythiaSingleton::~PythiaSingleton()
{
    // // Destructor
    if (beamMap.empty()) {
      if (fPythia) {
        delete fPythia;
        fPythia = 0;
      }
      if (fEventReader) {
        delete fEventReader;
        fEventReader = 0;
      }
      if (fBeamShape) {
        delete fBeamShape;
        fBeamShape = 0;
      }
    }
    else {
      map< pair<int, int>, pair<Pythia8::Pythia*,
        pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> > >::iterator itr;
      for (itr = beamMap.begin(); itr != beamMap.end(); ++itr) {
        if (itr->second.first) {
            delete itr->second.first;
            itr->second.first = 0;
        }
        if (itr->second.second.first) {
            delete itr->second.second.first;
            itr->second.second.first = 0;
        }
        if (itr->second.second.second) {
            delete itr->second.second.second;
            itr->second.second.second = 0;
        }
      }
    }

    if (fgInstance) {
        delete fgInstance;
        fgInstance = 0;
    }
}
//____________________________________________________________________________
bool PythiaSingleton::BeamConfigExists(int beamA, int beamB) 
{
    map< pair<int, int>, pair<Pythia8::Pythia*,
      pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> > >::iterator itr;
    itr = beamMap.find(pair<int, int>(beamA, beamB));
    if (itr != beamMap.end()) {
        return true;
    }
    else {
        return false;
    }
}
//____________________________________________________________________________
void PythiaSingleton::InitializeBeam(int beamA, int beamB, double eMax) 
{
    // Assign current pythia object to the map if it's empty.
    if (beamMap.empty()) {
        fEventReader = new Pythia8::LHAup_Genie();
        fBeamShape = new Pythia8::GBeamShape(eMax);
        pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> fUtil(fEventReader, fBeamShape);
        beamMap[pair<int, int>(beamA, beamB)] =
          pair<Pythia8::Pythia*, pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> >
          (fPythia, fUtil);
        return;
    }

    map< pair<int, int>, pair<Pythia8::Pythia*,
      pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> > >::iterator itr;
    itr = beamMap.find(pair<int, int>(beamA, beamB));
    if (itr != beamMap.end()) {
        // If the beam config exists, assign existing pythia object to fPythia.
        fPythia      = itr->second.first;
        fEventReader = itr->second.second.first;
        fBeamShape   = itr->second.second.second;
    }
    else {
        // Otherwise, create a new Pythia object.
        fPythia = new Pythia8::Pythia(
            fPythia->settings, fPythia->particleData, false
        );
        fEventReader = new Pythia8::LHAup_Genie();
        fBeamShape = new Pythia8::GBeamShape(eMax);
        pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> fUtil(fEventReader, fBeamShape);
        beamMap[pair<int, int>(beamA, beamB)] = pair<Pythia8::Pythia*,
          pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> >
          (fPythia, fUtil);
    }

    return;
}
//____________________________________________________________________________
