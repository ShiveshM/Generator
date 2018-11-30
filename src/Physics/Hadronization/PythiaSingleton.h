#ifndef _PYTHIA_SINGLETON_H_
#define _PYTHIA_SINGLETON_H_

// Avoid the inclusion of dlfcn.h by Pythia.h that CINT is not able to process
#ifdef __CINT__
#define _DLFCN_H_
#define _DLFCN_H
#endif

#include <TObject.h>
#include <TPrimary.h>

#include "Pythia8/Pythia.h"
using namespace Pythia8;


namespace Pythia8 {

//==========================================================================

// A derived class with initialization information from the HEPRUP
// Fortran commonblock and event information from the HEPEUP one.

class LHAup_Genie : public LHAup {

public:

  // Constructor.
  LHAup_Genie() {}

  // Fill init information from external program.
  bool fillInit(int idA, int idB, double eA, double eB) {
    fillBeam("A", idA, eA);
    fillBeam("B", idB, eB);
    fillProcess();
    // Done.
    return true;
  }
  // Clear init information from external program.
  bool clearInit() {
    fillBeam("A", 0, -1.);
    fillBeam("B", 0, -1.);
    fillProcess();
    // Done.
    return true;
  }

  // Translate stored init information to Pythia bookkeeping.
  bool setInit() {
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

  // Translate stored event information to Pythia bookkeeping.
  bool setEvent(int idProcIn = 0) {
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

  // Internal bookkeeping of beam information.
  int idBeamAStore, idBeamBStore, pdfgAStore, pdfgBStore, pdfsAStore, pdfsBStore;
  double eBeamAStore, eBeamBStore;
  bool fillBeam(string beam, int id, double e, int pdfg = -1, int pdfs = -1) {
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

  // Internal bookkeeping of process/cross section information.
  int idProcStore;
  double xsecStore, xsecerrStore, xwgtmaxStore;
  bool fillProcess(int idProcIn = 0, double xsec = 1., double xsecerr = 0.,
    double xwgtmax = 0.) {
    idProcStore  = idProcIn;    
    xsecStore    = xsec;
    xsecerrStore = xsecerr;
    xwgtmaxStore = xwgtmax;
    return true;
  }

  // Internal bookkeeping of event information.
  int nParticlesStore;
  double xwgtStore, scaleStore, aqedStore, aqcdStore; 
  vector<int> idPartStore, statusPartStore, mother1PartStore, mother2PartStore,
    colPartStore, acolPartStore;
  vector<double> vtimPartStore, spinPartStore;
  vector< vector<double> > pPartStore;

  // Fill information on current event.
  bool fillEventInfo( int nParticles, double xwgt, double scale,
    double aqed = 0.00729735, double aqcd = 0.13) {
    nParticlesStore = nParticles;
    xwgtStore = xwgt;
    scaleStore = scale;
    aqedStore = aqed;
    aqcdStore = aqcd;
    return true;
  }

  // Attach new particle to the event.
  bool fillNewParticle(int idPart, int statusPart, vector<double> pPart,
    int mother1Part = 0, int mother2Part = 0, int colPart = 0,
    int acolPart = 0, double vtimPart = -1., double spinPart = -9.){
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

  // Clear event information.
  bool clearEvent(){
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

  Settings* settingsPtr;
  void setPointers(Settings* settings) {settingsPtr = settings;}

protected:

  // User-written routine that does the intialization and fills heprup_out.
  virtual bool fillHepRup() {return true;}

  // User-written routine that does the event generation and fills hepeup_out.
  virtual bool fillHepEup() {return true;}

private:

  // Store beam energies to calculate x values.
//  double eBeamA, eBeamB;
  LHAscales scalesNow;

};

}

class Pythia;

class PythiaSingleton : public TObject{

public:
    PythiaSingleton();
    virtual ~PythiaSingleton();

    static PythiaSingleton * Instance ();
    Pythia8::Pythia *      Pythia8     () {return fPythia;}
    Pythia8::LHAup_Genie * EventReader () {return fEventReader;}

private:
    static PythiaSingleton * fgInstance;   ///< singleton instance
    Pythia8::Pythia      *   fPythia;      ///< PYTHIA8 instance
    Pythia8::LHAup_Genie *   fEventReader; ///< LHAup instance

ClassDef(PythiaSingleton,1)
};

#endif    // _PYTHIA_SINGLETON__H_
