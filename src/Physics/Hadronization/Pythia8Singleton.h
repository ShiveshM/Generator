//____________________________________________________________________________
/*!

\class    Pythia8Singleton

\brief    Provides access to the PYTHIA8 instance.

\author   Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
          Queen Mary University of London

\created  September 21, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _PYTHIA8_SINGLETON_H_
#define _PYTHIA8_SINGLETON_H_

// Avoid the inclusion of dlfcn.h by Pythia.h that CINT is not able to process
#ifdef __CINT__
#define _DLFCN_H_
#define _DLFCN_H
#endif

#include <TObject.h>
#include <TPrimary.h>

#include "Framework/Conventions/GBuild.h"

using namespace std;

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"

namespace Pythia8 {
//____________________________________________________________________________
class GBeamShape : public BeamShape {

public:

  // Constructor.
  GBeamShape(double eMax) : BeamShape() {this->eMax = eMax;}

  // Initialize beam parameters.
  // In this particular example we will reuse the existing settings names
  // but with modified meaning, so init() in the base class can be kept.
  // virtual void init( Settings& settings, Rndm* rndmPtrIn);

  // Set the two beam momentum deviations and the beam vertex.
  // Note that momenta are in units of GeV and vertices in mm,
  // always with c = 1, so that e.g. time is in mm/c.
  virtual void pick();

  double eMax;
  double pZspread;
  virtual void setSpread(double pIn) {pZspread = pIn;}
  
};
//____________________________________________________________________________
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
    bool fillNewParticle(int idPart, int statusPart, Pythia8::Vec4 pPart,
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
//____________________________________________________________________________
class Pythia;

#endif

class Pythia8Singleton : public TObject{

public:
    Pythia8Singleton();
    virtual ~Pythia8Singleton();

    static Pythia8Singleton * Instance () {
        return fgInstance ? fgInstance : (fgInstance = new Pythia8Singleton());
    }
#ifdef __GENIE_PYTHIA8_ENABLED__
    Pythia8::Pythia *      Pythia8     () {return fPythia;}
    Pythia8::LHAup_Genie * EventReader () {return fEventReader;}
    Pythia8::GBeamShape *  GBeamShape  () {return fBeamShape;}
    bool BeamConfigExists(int beamA, int beamB);
    void InitializeBeam  (int beamA, int beamB, double eMax);
#endif

private:
    static Pythia8Singleton * fgInstance;  ///< singleton instance
#ifdef __GENIE_PYTHIA8_ENABLED__
    map< pair<int, int>, pair<Pythia8::Pythia*,
      pair<Pythia8::LHAup_Genie*, Pythia8::GBeamShape*> > > beamMap;
    Pythia8::Pythia      * fPythia;      ///< PYTHIA8 instance
    Pythia8::LHAup_Genie * fEventReader; ///< LHAup instance
    Pythia8::GBeamShape  * fBeamShape;   ///< GBeamShape instance
#endif

ClassDef(Pythia8Singleton,1)
};

#endif    // _PYTHIA8_SINGLETON__H_
