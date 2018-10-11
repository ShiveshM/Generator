//____________________________________________________________________________
/*!

\class    genie::PythiaDecayer

\brief    Interface to PYTHIA particle decayer. \n
          The PythiaDecayer is a concrete implementation of the DecayModelI
          interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 20, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PYTHIA_DECAYER_I_H_
#define _PYTHIA_DECAYER_I_H_

#include "Physics/Decay/DecayModelI.h"
#include "Physics/Hadronization/PythiaSingleton.h"

namespace genie {

class PythiaDecayer : public DecayModelI {

public:

  PythiaDecayer();
  PythiaDecayer(string config);
  virtual ~PythiaDecayer();

  // implement the DecayModelI interface  
  bool           IsHandled      (int pdgc)                      const;
  void           Initialize     (void)                          const;
  TClonesArray * Decay          (const DecayerInputs_t & inp)   const;
  double         Weight         (void)                          const;
  void           InhibitDecay   (int pdgc, TDecayChannel * dc=0) const;
  void           UnInhibitDecay (int pdgc, TDecayChannel * dc=0) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig             (void);
  double SumBR                  (int pdgc) const;
  int    FindPythiaDecayChannel (int pdgc, TDecayChannel* dc) const;
  bool   MatchDecayChannels     (int pdgc, int ichannel, TDecayChannel * dc) const;

  mutable PythiaSingleton * fPythia8;  ///< PYTHIA8 wrapper class
  mutable double fWeight;
//bool fForceDecay;
};

}         // genie namespace

#endif    // _PYTHIA_DECAYER_H_
