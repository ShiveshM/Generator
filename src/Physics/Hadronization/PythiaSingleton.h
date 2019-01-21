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
    LHAup_Genie();

    // Fill init information from external program.
    bool fillInit(int idA, int idB, double eA, double eB);

    // Clear init information from external program.
    bool clearInit();

    // Translate stored init information to Pythia bookkeeping.
    bool setInit();

    // Translate stored event information to Pythia bookkeeping.
    bool setEvent(int idProcIn = 0);

    // Internal bookkeeping of beam information.
    int idBeamAStore, idBeamBStore, pdfgAStore, pdfgBStore, pdfsAStore, pdfsBStore;
    double eBeamAStore, eBeamBStore;
    bool fillBeam(string beam, int id, double e, int pdfg = -1, int pdfs = -1);

    // Internal bookkeeping of process/cross section information.
    int idProcStore;
    double xsecStore, xsecerrStore, xwgtmaxStore;
    bool fillProcess(int idProcIn = 0, double xsec = 1., double xsecerr = 0.,
        double xwgtmax = 0.);

    // Internal bookkeeping of event information.
    int nParticlesStore;
    double xwgtStore, scaleStore, aqedStore, aqcdStore; 
    vector<int> idPartStore, statusPartStore, mother1PartStore, mother2PartStore,
        colPartStore, acolPartStore;
    vector<double> vtimPartStore, spinPartStore;
    vector< vector<double> > pPartStore;

    // Fill information on current event.
    bool fillEventInfo(int nParticles, double xwgt, double scale,
        double aqed = 0.00729735, double aqcd = 0.13);

    // Attach new particle to the event.
    bool fillNewParticle(int idPart, int statusPart, vector<double> pPart,
        int mother1Part = 0, int mother2Part = 0, int colPart = 0,
        int acolPart = 0, double vtimPart = -1., double spinPart = -9.);

    // Clear event information.
    bool clearEvent();

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
    ~PythiaSingleton();

    static PythiaSingleton * Instance ();
    Pythia8::Pythia *      Pythia8     () {return fPythia;}
    Pythia8::LHAup_Genie * EventReader () {return fEventReader;}
    bool BeamConfigExists(int beamA, int beamB);
    void InitializeBeam  (int beamA, int beamB);

private:
    static PythiaSingleton * fgInstance;   ///< singleton instance

    std::map< std::pair<int, int>,
              std::pair<Pythia8::Pythia*, Pythia8::LHAup_Genie*> > beamMap;
    Pythia8::Pythia * fPythia;             ///< PYTHIA8 instance
    Pythia8::LHAup_Genie * fEventReader;   ///< LHAup instance

ClassDef(PythiaSingleton,1)
};

#endif    // _PYTHIA_SINGLETON__H_
