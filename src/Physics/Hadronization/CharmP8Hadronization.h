//____________________________________________________________________________
/*!

\class    genie::CharmP8Hadronization

\brief    Provides access to the PYTHIA6 hadronization models. \n
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

          Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
          Queen Mary University of London

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CHARM_P8_HADRONIZATION_H_
#define _CHARM_P8_HADRONIZATION_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Interaction.h"

#include "Physics/Hadronization/Pythia8Singleton.h"
class TF1;

namespace genie {

class Spline;
class FragmentationFunctionI;

class CharmP8Hadronization : public EventRecordVisitorI {

public:
  CharmP8Hadronization();
  CharmP8Hadronization(string config);
  virtual ~CharmP8Hadronization();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  //
  void Configure(const Registry & config);
  void Configure(string config);


private:

  void           LoadConfig          (void);
  void           Initialize          (void)                                    const ;
  TClonesArray * Hadronize           (const Interaction* )                     const ;
  int            GenerateCharmHadron (int nupdg, double EvLab)                 const ;

  double         Weight              (void)                                    const ;

  mutable TGenPhaseSpace fPhaseSpaceGenerator; ///< a phase space generator

  // Configuration parameters
  //
  bool                           fCharmOnly;   ///< don't hadronize non-charm blob
  TF1 *                          fCharmPT2pdf; ///< charm hadron pT^2 pdf
  const FragmentationFunctionI * fFragmFunc;   ///< charm hadron fragmentation func
  Spline *                       fD0FracSpl;   ///< nu charm fraction vs Ev: D0
  Spline *                       fDpFracSpl;   ///< nu charm fraction vs Ev: D+
  Spline *                       fDsFracSpl;   ///< nu charm fraction vs Ev: Ds+
  double                         fD0BarFrac;   ///< nubar \bar{D0} charm fraction
  double                         fDmFrac;      ///< nubar D- charm fraction
  mutable Pythia8Singleton *     fPythia;      ///< remnant (non-charm) hadronizer
};

}         // genie namespace

#endif    // _CHARM_P8_HADRONIZATION_H_

