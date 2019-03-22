//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TMath.h>
#include <TH1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Physics/Decay/DecayModelI.h"
#include "Physics/Hadronization/PythiaHadronization.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Hadronization/FragmRecUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
PythiaHadronization::PythiaHadronization() :
HadronizationModelBase("genie::PythiaHadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaHadronization::PythiaHadronization(string config) :
HadronizationModelBase("genie::PythiaHadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaHadronization::~PythiaHadronization()
{
}
//____________________________________________________________________________
void PythiaHadronization::Initialize(void) const
{
  fPythia8 = PythiaSingleton::Instance();
  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  pythia8->settings.resetAll();

  // sync GENIE/PYTHIA8 seed number
  RandomGen::Instance();
}
//____________________________________________________________________________
TClonesArray * 
  PythiaHadronization::Hadronize(
         const Interaction * interaction) const
{
  LOG("PythiaHad", pNOTICE) << "Running PYTHIA hadronizer";

  if(!this->AssertValidity(interaction)) {
     LOG("PythiaHad", pERROR) << "Returning a null particle list!";
     return 0;
  }

  // get kinematics / init-state / process-info

  const Kinematics &   kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  assert(target.HitQrkIsSet()); 

  double W = kinematics.W();

  int  probe       = init_state.ProbePdg();
  int  hit_nucleon = target.HitNucPdg();
  int  hit_quark   = target.HitQrkPdg();
  int  out_lepton  = interaction->FSPrimLeptonPdg();
  bool from_sea    = target.HitSeaQrk();

  LOG("PythiaHad", pNOTICE)
          << "Hit nucleon pdgc = " << hit_nucleon << ", W = " << W;
  LOG("PythiaHad", pNOTICE)
            << "Selected hit quark pdgc = " << hit_quark
                           << ((from_sea) ? "[sea]" : "[valence]");

  // check hit-nucleon assignment, input neutrino & interaction type
  bool isp  = pdg::IsProton           (hit_nucleon);
  bool isn  = pdg::IsNeutron          (hit_nucleon);
  bool isv  = pdg::IsNeutrino         (probe);
  bool isvb = pdg::IsAntiNeutrino     (probe);
//bool isl  = pdg::IsNegChargedLepton (probe);
//bool islb = pdg::IsPosChargedLepton (probe);
  bool iscc = proc_info.IsWeakCC      ();
  bool isnc = proc_info.IsWeakNC      ();
  bool isdm = proc_info.IsDarkMatter  ();
  bool isem = proc_info.IsEM          ();
  bool isu  = pdg::IsUQuark           (hit_quark);
  bool isd  = pdg::IsDQuark           (hit_quark);
  bool iss  = pdg::IsSQuark           (hit_quark);
  bool isub = pdg::IsAntiUQuark       (hit_quark);
  bool isdb = pdg::IsAntiDQuark       (hit_quark);
  bool issb = pdg::IsAntiSQuark       (hit_quark);

  //
  // Generate the quark system (q + qq) initiating the hadronization
  //

  int final_quark = 0; // leading quark (hit quark after the interaction)

  // Figure out the what happens to the hit quark after the interaction
  if (isnc || isem || isdm) {
    // NC, EM
    final_quark = hit_quark;
  } else {
    // CC
    if      (isv  && isd ) final_quark = kPdgUQuark;
    else if (isv  && iss ) final_quark = kPdgUQuark;
    else if (isv  && isub) final_quark = kPdgAntiDQuark;
    else if (isvb && isu ) final_quark = kPdgDQuark;
    else if (isvb && isdb) final_quark = kPdgAntiUQuark;
    else if (isvb && issb) final_quark = kPdgAntiUQuark;
    else {
      LOG("PythiaHad", pERROR)
        << "Not allowed mode. Refused to make a final quark assignment!";
      return 0;
    }
  }//CC

  //
  // PYTHIA -> HADRONIZATION
  //

  LOG("PythiaHad", pNOTICE)
        << "Fragmentation System: "
        << "Hit q = " << hit_quark << ", Final q = " << final_quark;

  // Kinematics.
  RefFrame_t rf = kRfLab; // LAB frame
  TLorentzVector probe_p4       = *init_state.GetProbeP4(rf);
  TLorentzVector hit_nucleon_p4 = target.HitNucP4();
  TLorentzVector out_lepton_p4  = kinematics.FSLeptonP4();

  Pythia8::Vec4 probeV4 = Pythia8::Vec4(
      probe_p4.Px(), probe_p4.Py(), probe_p4.Pz(), probe_p4.E()
      );
  Pythia8::Vec4 hitNucV4 = Pythia8::Vec4(
      hit_nucleon_p4.Px(), hit_nucleon_p4.Py(), hit_nucleon_p4.Pz(), hit_nucleon_p4.E()
      );
  Pythia8::Vec4 outLepV4 = Pythia8::Vec4(
      out_lepton_p4.Px(), out_lepton_p4.Py(), out_lepton_p4.Pz(), out_lepton_p4.E()
      );

  // Neutrino energy in LAB frame.
  double eNu = probeV4.e();

  // Maximum neutrino energy allowable.
  double eMax = 1e5; // GeV

  // Setup Pythia object.
  bool beamConfigExists = fPythia8->BeamConfigExists(probe, hit_nucleon);
  fPythia8->InitializeBeam(probe, hit_nucleon, eMax);
  Pythia8::Pythia      * pythia8     = fPythia8->Pythia8();
  Pythia8::LHAup_Genie * eventReader = fPythia8->EventReader();
  Pythia8::GBeamShape  * beamShape   = fPythia8->GBeamShape();

  // Boost to CMS of probe and hit nucleon.
  Pythia8::RotBstMatrix toCMS;
  toCMS.toCMframe(probeV4, hitNucV4);
  probeV4.rotbst(toCMS);
  hitNucV4.rotbst(toCMS);
  outLepV4.rotbst(toCMS);

  // Calculate hit quark momentum.
  double Q2 = -(probe_p4 - out_lepton_p4).M2();
  double x  = Q2 / (2 * hit_nucleon_p4.M() * (probe_p4.E() - out_lepton_p4.E()));
  TLorentzVector hit_quark_p4   = TLorentzVector(
      0., 0., -hit_nucleon_p4.M() * x, hit_nucleon_p4.M() * x
  );
  Pythia8::Vec4 hitQrkV4 = Pythia8::Vec4(
      hit_quark_p4.Px(), hit_quark_p4.Py(), hit_quark_p4.Pz(), hit_quark_p4.E()
  );

  // Calculate final quark 4 momentum.
  Pythia8::Vec4 finQrkV4 = probeV4 + hitQrkV4 - outLepV4;
  double mom3_2 = pow(finQrkV4.px(),2)+pow(finQrkV4.py(),2)+pow(finQrkV4.pz(),2);
  double finQrkMass_2 = pow(pythia8->particleData.m0(final_quark),2);
  finQrkV4.e(sqrt(mom3_2+finQrkMass_2));

  // Determine how jetset treats un-stable particles appearing in hadronization

  bool pi0_decflag = pythia8->particleData.canDecay(kPdgPi0);
  bool K0_decflag  = pythia8->particleData.canDecay(kPdgK0);
  bool K0b_decflag = pythia8->particleData.canDecay(kPdgAntiK0);
  bool L0_decflag  = pythia8->particleData.canDecay(kPdgLambda);
  bool L0b_decflag = pythia8->particleData.canDecay(kPdgAntiLambda);
  bool Dm_decflag  = pythia8->particleData.canDecay(kPdgP33m1232_DeltaM);
  bool D0_decflag  = pythia8->particleData.canDecay(kPdgP33m1232_Delta0);
  bool Dp_decflag  = pythia8->particleData.canDecay(kPdgP33m1232_DeltaP);
  bool Dpp_decflag = pythia8->particleData.canDecay(kPdgP33m1232_DeltaPP);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("PythiaHad", pDEBUG) << "Original decay flag for pi0           =  " << pi0_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for K0            =  " << K0_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for \bar{K0}      =  " << K0b_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for Lambda        =  " << L0_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for \bar{Lambda0} =  " << L0b_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for D-            =  " << Dm_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for D0            =  " << D0_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for D+            =  " << Dp_decflag;
  LOG("PythiaHad", pDEBUG) << "Original decay flag for D++           =  " << Dpp_decflag;
#endif

  pythia8->particleData.mayDecay(kPdgPi0,              false ); // don't decay pi0
  pythia8->particleData.mayDecay(kPdgK0,               false ); // don't decay K0
  pythia8->particleData.mayDecay(kPdgAntiK0,           false ); // don't decay \bar{K0}
  pythia8->particleData.mayDecay(kPdgLambda,           false ); // don't decay Lambda0
  pythia8->particleData.mayDecay(kPdgAntiLambda,       false ); // don't decay \bar{Lambda0}
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaM,  true  ); // decay Delta-
  pythia8->particleData.mayDecay(kPdgP33m1232_Delta0,  true  ); // decay Delta0
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaP,  true  ); // decay Delta+
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaPP, true  ); // decay Delta++

  // -- hadronize --
  if (!beamConfigExists) {
    // Set up Pythia to read hard scattering from external provider (via
    // the Les Houches Event Accord functionalities)
    pythia8->readString("ProcessLevel:all = on");
    pythia8->readString("Beams:frameType = 5");
    pythia8->readString("Check:beams = off");
    pythia8->readString("TimeShower:QEDshowerByL  = off");
    pythia8->readString("LesHouches:setLeptonMass = 2");
    pythia8->readString("LesHouches:setQuarkMass  = 2");
    pythia8->readString("LesHouches:matchInOut    = off");
    pythia8->readString("LesHouches:mRecalculate   = 1e-6");  
    pythia8->readString("Beams:allowMomentumSpread = on");
    pythia8->readString("BeamRemnants:primordialKT = on");
    pythia8->readString("BeamRemnants:primordialKTsoft = 0.0");
    pythia8->readString("BeamRemnants:primordialKThard = 0.0");    
    pythia8->readString("Print:quiet = on");

    // Hand BeamShape pointer to Pythia.
    pythia8->setBeamShapePtr(beamShape);

    // Set up a global instance of LHAup.
    eventReader->setPointers(&pythia8->settings);
    eventReader->fillInit(
        probe, hit_nucleon, eMax, pythia8->particleData.m0(hit_nucleon)
    );
    eventReader->setInit();
    pythia8->setLHAupPtr(eventReader);

    // Initialize Pythia.
    pythia8->init();
  }
  else {
    pythia8->event.reset();
  }

  // The variable neutrino energy can be handled using the BeamShape interface.
  // I have to reinitialize the BeamShape in each event, but this is only
  // because I want to set a parameter that represents the reduction in energy
  // for the incoming neutrino.
  pythia8->settings.parm("Beams:sigmaPzA", eNu);
  beamShape->init(pythia8->settings, &pythia8->rndm);

  // This should set the LHA event using fortran common blocks
  eventReader->clearEvent();
  eventReader->fillEventInfo(4, 1.0, 10.0);

  // Incoming particles.
  eventReader->fillNewParticle(probe,     -1, probeV4);
  eventReader->fillNewParticle(hit_quark, -1, hitQrkV4);

  // Outgoing particles.
  eventReader->fillNewParticle(out_lepton,  1, outLepV4);
  eventReader->fillNewParticle(final_quark, 1, finQrkV4);
  eventReader->setEvent(); 

  // Now call Pythia to process information.
  pythia8->next();

  // List the event information
  pythia8->event.list();
  pythia8->stat();

  // restore pythia decay settings so as not to interfere with decayer 
  pythia8->particleData.mayDecay(kPdgPi0,             pi0_decflag);
  pythia8->particleData.mayDecay(kPdgK0,              K0_decflag);
  pythia8->particleData.mayDecay(kPdgAntiK0,          K0b_decflag);
  pythia8->particleData.mayDecay(kPdgLambda,          L0_decflag);
  pythia8->particleData.mayDecay(kPdgAntiLambda,      L0b_decflag);
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaM, Dm_decflag);
  pythia8->particleData.mayDecay(kPdgP33m1232_Delta0, D0_decflag);
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaP, Dp_decflag);
  pythia8->particleData.mayDecay(kPdgP33m1232_DeltaPP,Dpp_decflag);

  // get record
  Pythia8::Event &fEvent = pythia8->event;
  int numpart = fEvent.size();
  // assert(numpart>0);
  if (numpart == 0) {
     LOG("PythiaHad", pERROR) << "Returning a null particle list!";
     return 0;
  }

  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", numpart);
  particle_list->SetOwner(true);

  // Offset the initial particles.
  int ioff = -1; int ioffUpper = -1;
  for (int i = 0; i < numpart; ++i) {
    if (fEvent[i].status() > 80 && !(fEvent[i].id() == out_lepton)) {
      ioff = fEvent[i].mother1();
      ioffUpper = fEvent[i].mother2();
      break;
    }
  }
  assert(ioff > -1);

  // Boost into hadronic CM frame.
  Pythia8::RotBstMatrix toHadCMS;
  toHadCMS.toCMframe(fEvent[ioff].p(), fEvent[ioffUpper].p());

  for (int i = ioff; i < numpart; ++i) {
    /*
     * Convert Pythia8 status code to Pythia6
     * Initial quark has a pythia6 status code of 12
     * The initial quark and the fragmented particles have a pythia6 code of 11 (kIStNucleonTarget)
     * Final state particles have a positive pythia8 code and a pythia6 code of 1 (kIStStableFinalState)
     */

    // Modify mother indexes.
    int mother1 = fEvent[i].mother1();
    int mother2 = fEvent[i].mother2();
    if (mother1 < ioff) mother1 = -1;
    else mother1 -= ioff;
    if (mother2 < ioff) mother2 = -1;
    else mother2 -= ioff;

    // Modify daughter indexes.
    int daughter1 = fEvent[i].daughter1();
    int daughter2 = fEvent[i].daughter2();
    if (daughter1 == 0) daughter1 = -1;
    else daughter1 -= ioff;
    if (daughter2 == 0) daughter2 = -1;
    else daughter2 -= ioff;

    GHepStatus_t gStatus;
    gStatus = (fEvent[i].status()>0) ? kIStStableFinalState : kIStNucleonTarget;

    LOG("PythiaHad", pDEBUG)
        << "Adding final state particle pdgc = " << fEvent[i].id()
        << " with status = " << gStatus;

    if (fEvent[i].status() > 0){
      if( pdg::IsQuark  (fEvent[i].id()) || 
              pdg::IsDiQuark(fEvent[i].id()) ) {
        LOG("PythiaHad", pERROR)
            << "Hadronization failed! Bare quark/di-quarks appear in final state!";
        particle_list->Delete();
        delete particle_list;
        return 0;            
      }
    }

    // Boost back to LAB frame.
    Pythia8::Vec4 p = fEvent[i].p();
    p.rotbst(toHadCMS);

    new((*particle_list)[i]) GHepParticle(
            fEvent[i].id(),
            gStatus,
            mother1,
            mother2,
            daughter1,
            daughter2,
            p.px(),               // [GeV/c]
            p.py(),               // [GeV/c]
            p.pz(),               // [GeV/c]
            p.e(),                // [GeV]
            fEvent[i].xProd(),    // [mm]
            fEvent[i].yProd(),    // [mm]
            fEvent[i].zProd(),    // [mm]
            fEvent[i].tProd());   // [mm/c]
  }

  utils::fragmrec::Print(particle_list);
  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * 
   PythiaHadronization::SelectParticles(
            const Interaction * interaction) const
{
// Works the opposite way (compared with the KNO hadronization model)
// Rather than having this method as one of the hadronization model components,
// we extract the list of particles from the fragmentation record after the
// hadronization has been completed.

  TClonesArray * particle_list = this->Hadronize(interaction);

  if(!particle_list) return 0;

  bool allowdup=true;
  PDGCodeList * pdgcv = new PDGCodeList(allowdup);
  pdgcv->reserve(particle_list->GetEntries());

  GHepParticle * particle = 0;
  TIter particle_iter(particle_list);

  while ((particle = (GHepParticle *) particle_iter.Next())) 
  {
    if (particle->Status()==kIStStableFinalState) pdgcv->push_back(particle->Pdg());
  }
  particle_list->Delete();
  delete particle_list;

  return pdgcv;
}
//____________________________________________________________________________
TH1D * PythiaHadronization::MultiplicityProb(
     const Interaction * interaction, Option_t * opt) const
{
// Similar comments apply as in SelectParticles()

  if(!this->AssertValidity(interaction)) {
     LOG("PythiaHad", pWARN) 
                << "Returning a null multipicity probability distribution!";
     return 0;
  }
  double maxmult   = this->MaxMult(interaction);
  TH1D * mult_prob = this->CreateMultProbHist(maxmult);

  const int nev=500;
  GHepParticle * particle = 0;

  for(int iev=0; iev<nev; iev++) {

     TClonesArray * particle_list = this->Hadronize(interaction);
     double         weight        = this->Weight();

     if(!particle_list) { iev--; continue; }

     int n = 0;
     TIter particle_iter(particle_list);
     while ((particle = (GHepParticle *) particle_iter.Next())) 
     {
       if (particle->Status()==kIStStableFinalState) n++;
     }   
     particle_list->Delete();
     delete particle_list;
     mult_prob->Fill( (double)n, weight);
  }

  double integral = mult_prob->Integral("width");
  if(integral>0) {
    // Normalize the probability distribution
    mult_prob->Scale(1.0/integral);
  } else {
    SLOG("PythiaHad", pWARN) << "probability distribution integral = 0";
    return mult_prob;
  }

  string option(opt);

  bool apply_neugen_Rijk = option.find("+LowMultSuppr") != string::npos;
  bool renormalize       = option.find("+Renormalize")  != string::npos;

  // Apply the NeuGEN probability scaling factors -if requested-
  if(apply_neugen_Rijk) {
    SLOG("KNOHad", pINFO) << "Applying NeuGEN scaling factors";
     // Only do so for W<Wcut
     const Kinematics & kinematics = interaction->Kine();
     double W = kinematics.W();
     if(W<fWcut) {
       this->ApplyRijk(interaction, renormalize, mult_prob);
     } else {
        SLOG("PythiaHad", pDEBUG)
              << "W = " << W << " < Wcut = " << fWcut
                                << " - Will not apply scaling factors";
     }//<wcut?
  }//apply?

  return mult_prob;
}
//____________________________________________________________________________
double PythiaHadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
void PythiaHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaHadronization::LoadConfig(void)
{
  // the configurable PYTHIA parameters used here are the ones used by NUX 
  // (see A.Rubbia's talk @ NuINT-01)
  // The defaults are the values used by PYTHIA
  // Use the NUX config set to set the tuned values as used in NUX.

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();

  GetParam( "PYTHIA-SSBarSuppression", fSSBarSuppression ) ;
  GetParam( "PYTHIA-GaussianPt2",      fGaussianPt2      ) ;
  // TODO: find PYTHIA8 equivalent of this parameter
  // GetParam( "PYTHIA-NonGaussianPt2Tail", fNonGaussianPt2Tail  ) ;
  GetParam( "PYTHIA-RemainingEnergyCutoff", fRemainingECutoff ) ;

  // fPythia->SetPARJ(23, fNonGaussianPt2Tail);
  pythia8->settings.parm("StringFlav:probStoUD", fSSBarSuppression);
  pythia8->settings.parm("Diffraction:primKTwidth", fGaussianPt2);
  pythia8->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  GetParam( "Wcut", fWcut ) ;

  // decayer
  fDecayer = 0;
  if( GetConfig().Exists("Decayer") ) {
     fDecayer = dynamic_cast<const DecayModelI *> (this->SubAlg("Decayer"));
     assert(fDecayer);
  }

  // Load NEUGEN multiplicity probability scaling parameters Rijk
   //neutrinos
   GetParam( "DIS-HMultWgt-vp-CC-m2",  fRvpCCm2  ) ;
   GetParam( "DIS-HMultWgt-vp-CC-m3",  fRvpCCm3  ) ;
   GetParam( "DIS-HMultWgt-vp-NC-m2",  fRvpNCm2  ) ;
   GetParam( "DIS-HMultWgt-vp-NC-m3",  fRvpNCm3  ) ;
   GetParam( "DIS-HMultWgt-vn-CC-m2",  fRvnCCm2  ) ;
   GetParam( "DIS-HMultWgt-vn-CC-m3",  fRvnCCm3  ) ;
   GetParam( "DIS-HMultWgt-vn-NC-m2",  fRvnNCm2  ) ;
   GetParam( "DIS-HMultWgt-vn-NC-m3",  fRvnNCm3  ) ;
   //Anti-neutrinos
   GetParam( "DIS-HMultWgt-vbp-CC-m2", fRvbpCCm2 ) ;
   GetParam( "DIS-HMultWgt-vbp-CC-m3", fRvbpCCm3 ) ;
   GetParam( "DIS-HMultWgt-vbp-NC-m2", fRvbpNCm2 ) ;
   GetParam( "DIS-HMultWgt-vbp-NC-m3", fRvbpNCm3 ) ;
   GetParam( "DIS-HMultWgt-vbn-CC-m2", fRvbnCCm2 ) ;
   GetParam( "DIS-HMultWgt-vbn-CC-m3", fRvbnCCm3 ) ;
   GetParam( "DIS-HMultWgt-vbn-NC-m2", fRvbnNCm2 ) ;
   GetParam( "DIS-HMultWgt-vbn-NC-m3", fRvbnNCm3 ) ;

  LOG("PythiaHad", pDEBUG) << GetConfig() ;
}
//____________________________________________________________________________
bool PythiaHadronization::AssertValidity(const Interaction * interaction) const
{
  // check that there is no charm production 
  // (GENIE uses a special model for these cases)
  if(interaction->ExclTag().IsCharmEvent()) {
     LOG("PythiaHad", pWARN) << "Can't hadronize charm events";
     return false;
  }
  // check the available mass
  double W = utils::kinematics::W(interaction);
  if(W < this->Wmin()) {
     LOG("PythiaHad", pWARN) << "Low invariant mass, W = " << W << " GeV!!";
     return false;
  }

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  if( ! target.HitQrkIsSet() ) {
     LOG("PythiaHad", pWARN) << "Hit quark was not set!";
     return false;
  }

  int  probe       = init_state.ProbePdg();
  int  hit_nucleon = target.HitNucPdg();
  int  hit_quark   = target.HitQrkPdg();
//bool from_sea    = target.HitSeaQrk();

  // check hit-nucleon assignment, input neutrino & weak current
  bool isp  = pdg::IsProton           (hit_nucleon);
  bool isn  = pdg::IsNeutron          (hit_nucleon);
  bool isv  = pdg::IsNeutrino         (probe);
  bool isvb = pdg::IsAntiNeutrino     (probe);
  bool isdm = pdg::IsDarkMatter         (probe);
  bool isl  = pdg::IsNegChargedLepton (probe);
  bool islb = pdg::IsPosChargedLepton (probe);
  bool iscc = proc_info.IsWeakCC      ();
  bool isnc = proc_info.IsWeakNC      ();
  bool isdmi = proc_info.IsDarkMatter  ();
  bool isem = proc_info.IsEM          ();
  if( !(iscc||isnc||isem||isdmi) ) {
    LOG("PythiaHad", pWARN) 
       << "Can only handle electro-weak interactions";
    return false;
  }
  if( !(isp||isn) || !(isv||isvb||isl||islb||isdm) ) {
    LOG("PythiaHad", pWARN) 
      << "Invalid initial state: probe = " 
      << probe << ", hit_nucleon = " << hit_nucleon;
    return false;
  }

  // assert that the interaction mode is allowed
  bool isu  = pdg::IsUQuark     (hit_quark);
  bool isd  = pdg::IsDQuark     (hit_quark);
  bool iss  = pdg::IsSQuark     (hit_quark);
  bool isub = pdg::IsAntiUQuark (hit_quark);
  bool isdb = pdg::IsAntiDQuark (hit_quark);
  bool issb = pdg::IsAntiSQuark (hit_quark);

  bool allowed = (iscc && isv  && (isd||isub||iss))  ||
                 (iscc && isvb && (isu||isdb||issb)) ||
                 (isnc && (isv||isvb) && (isu||isd||isub||isdb||iss||issb)) ||
                 (isdmi && isdm && (isu||isd||isub||isdb||iss||issb)) ||
                 (isem && (isl||islb) && (isu||isd||isub||isdb||iss||issb));
  if(!allowed) {
    LOG("PythiaHad", pWARN) 
      << "Impossible interaction type / probe / hit quark combination!";
    return false;
  }

  return true;
}
//____________________________________________________________________________
/*
void PythiaHadronization::SwitchDecays(int pdgc, bool on_off) const
{
  LOG("PythiaHad", pNOTICE)
     << "Switching " << ((on_off) ? "ON" : "OFF")
                     << " all PYTHIA decay channels for particle = " << pdgc;

  int flag     = (on_off) ? 1 : 0;
  int kc       = fPythia->Pycomp(pdgc);
  int first_ch = fPythia->GetMDCY(kc,2);
  int last_ch  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  for(int ich = first_ch; ich < last_ch; ich++) fPythia->SetMDME(ich,1,flag);
}
*/
//____________________________________________________________________________
/*
void PythiaHadronization::HandleDecays(TClonesArray * plist) const
{
// Handle decays of unstable particles if requested through the XML config.
// The default is not to decay the particles at this stage (during event
// generation, the UnstableParticleDecayer event record visitor decays what
// is needed to be decayed later on). But, when comparing various models
// (eg PYTHIA vs KNO) independently and not within the full MC simulation
// framework it might be necessary to force the decays at this point.

  if(!fDecayer) {
    LOG("PythiaHad", pWARN) << "No decayer was specified!";
    return;
  }

  this->SwitchDecays(kPdgLambda,     true); // decay Lambda
  this->SwitchDecays(kPdgAntiLambda, true); // decay \bar{Lambda}
  this->SwitchDecays(kPdgSigmaP,     true); // decay Sigma+
  this->SwitchDecays(kPdgSigma0,     true); // decay Sigma0
  this->SwitchDecays(kPdgSigmaM,     true); // decay Sigma-
  this->SwitchDecays(kPdgAntiSigmaP, true); // decay Sigma+
  this->SwitchDecays(kPdgAntiSigma0, true); // decay Sigma0
  this->SwitchDecays(kPdgAntiSigmaM, true); // decay Sigma-
  this->SwitchDecays(kPdgXi0,        true); // decay Xi0
  this->SwitchDecays(kPdgXiM,        true); // decay Xi-
  this->SwitchDecays(kPdgAntiXi0,    true); // decay \bar{Xi0}
  this->SwitchDecays(kPdgAntiXiP,    true); // decay \bar{Xi+}
  this->SwitchDecays(kPdgOmegaM,     true); // decay Omega-
  this->SwitchDecays(kPdgAntiOmegaP, true); // decay \bar{Omega+}

  int mstj21 = fPythia->GetMSTJ(21);
  fPythia->SetMSTJ(21,1); 
  fPythia->SetMSTJ(22,2);                  
  fPythia->SetPARJ(71,100);                  

  //-- loop through the fragmentation event record & decay unstables
  int idecaying   = -1; // position of decaying particle
  GHepParticle * p =  0; // current particle

  TIter piter(plist);
  while ( (p = (GHepParticle *) piter.Next()) ) {
     idecaying++;
     GHepStatus_t status = p->Status();
     int pdg    = p->Pdg();

     bool decay_it = (status<10) && 
                     ( pdg == kPdgLambda ||
                       pdg == kPdgAntiLambda ||
                       pdg == kPdgSigmaP ||
                       pdg == kPdgSigma0 ||
                       pdg == kPdgSigmaM ||
                       pdg == kPdgAntiSigmaP ||
                       pdg == kPdgAntiSigma0 ||
                       pdg == kPdgAntiSigmaM ||
                       pdg == kPdgXi0 ||
                       pdg == kPdgXiM ||
                       pdg == kPdgAntiXi0 ||
                       pdg == kPdgAntiXiP ||
                       pdg == kPdgOmegaM  ||
                       pdg == kPdgAntiOmegaP );

     // bother for final state particle only
     if(decay_it) {

          LOG("PythiaHad", pINFO)
                     << "Decaying particle with pdgc = " << p->Pdg();

          DecayerInputs_t dinp;

          TLorentzVector p4;
          p4.SetPxPyPzE(p->Px(), p->Py(), p->Pz(), p->Energy());

          dinp.PdgCode = p->Pdg();
          dinp.P4      = &p4;

          TClonesArray * decay_products = fDecayer->Decay(dinp);
          if(decay_products) {
                  //--  mark the parent particle as decayed & set daughters
                  p->SetStatus(kIStNucleonTarget);

                  int nfp = plist->GetEntries();          // n. fragm. products
                  int ndp = decay_products->GetEntries(); // n. decay products

                  p->SetFirstDaughter ( nfp );          // decay products added at
                  p->SetLastDaughter  ( nfp + ndp -1 ); // the end of the fragm.rec.

                  //--  add decay products to the fragmentation record
                  GHepParticle * dp = 0;
                  TIter dpiter(decay_products);

                  while ( (dp = (GHepParticle *) dpiter.Next()) ) {
                    if(dp->Status()>10) continue;
                    dp->SetFirstMother(idecaying);
                    new ( (*plist)[plist->GetEntries()] ) GHepParticle(*dp);
                  }

                  //-- clean up decay products
                  decay_products->Delete();
                  delete decay_products;
           }

     } // KS < 10 : final state particle (as in PYTHIA LUJETS record)
  } // particles in fragmentation record

  fPythia->SetMSTJ(21,mstj21); // restore mstj(21)
}
*/
//____________________________________________________________________________

