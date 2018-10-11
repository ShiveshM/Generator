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

class Pythia;

class PythiaSingleton : public TObject{

public:
    PythiaSingleton();
    virtual ~PythiaSingleton();

    static PythiaSingleton * Instance ();
    Pythia8::Pythia * Pythia8  () {return fPythia;}

private:
    static PythiaSingleton * fgInstance;  ///< singleton instance
    Pythia8::Pythia * fPythia;    ///< PYTHIA8 instance

ClassDef(PythiaSingleton,1)
};

#endif    // _PYTHIA_SINGLETON__H_
