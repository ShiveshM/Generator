//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

         Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
         Queen Mary University of London
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TH1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h" 
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Hadronization/Pythia8Hadronization.h"
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
Pythia8Hadronization::Pythia8Hadronization() :
EventRecordVisitorI("genie::Pythia8Hadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadronization::Pythia8Hadronization(string config) :
EventRecordVisitorI("genie::Pythia8Hadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadronization::~Pythia8Hadronization()
{

}
//____________________________________________________________________________
void Pythia8Hadronization::Initialize(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = Pythia8Singleton::Instance();
  Pythia8::Pythia * pythia8 = fPythia->Pythia8();
  pythia8->settings.resetAll();

  // sync GENIE/PYTHIA8 seed number
  RandomGen::Instance();
#endif
}
//____________________________________________________________________________
void Pythia8Hadronization::ProcessEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  TClonesArray * particle_list = this->Hadronize(interaction);

  if(! particle_list ) {
    LOG("Pythia8Hadronization", pWARN) << "Got an empty particle list. Hadronizer failed!";
    LOG("Pythia8Hadronization", pWARN) << "Quitting the current event generation thread";
    
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    
    return;
   }

  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = is_nucleus ? kIStHadronInTheNucleus : kIStStableFinalState ;
  

  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);


  GHepParticle * neutrino  = event->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  GHepParticle * p = 0;
  TIter particle_iter(particle_list);
  while ((p = (GHepParticle *) particle_iter.Next())) 
  {

    int pdgc = p -> Pdg() ;
        
    // set the proper status according to a number of things:
    // interaction on a nucleaus or nucleon, particle type
    GHepStatus_t ist = ( p -> Status() ==1 ) ? istfin : kIStDISPreFragmHadronicState;
    // not that this status here ^ is the pythia state

    // handle gammas, and leptons that might come from internal pythia decays
    // mark them as final state particles
    bool not_hadron = ( pdgc == kPdgGamma ||
			pdg::IsNeutralLepton(pdgc) ||
			pdg::IsChargedLepton(pdgc) ) ;

    if( not_hadron )  { ist = kIStStableFinalState; }
    p -> SetStatus( ist ) ;
 
    p->SetFirstMother(mom + p->FirstMother() );
    p->SetLastMother( -1 );
    
    // In Pythia6 having no daughters means daughter == 0 hnce the following check
    int ifd = (p->FirstDaughter() <= 0 ) ? -1 : mom  + p->FirstDaughter();
    int ild = (p->LastDaughter()  <= 0 ) ? -1 : mom  + p->LastDaughter();
    p->SetFirstDaughter(ifd);
    p->SetLastDaughter (ild);
    
    // the Pythia particle position is overridden
    p -> SetPosition( vtx ) ; 

    event->AddParticle(*p);
  }

  particle_list -> Delete() ; 
  delete particle_list ;
}
//____________________________________________________________________________
TClonesArray * 
  Pythia8Hadronization::Hadronize(
         const Interaction * interaction) const
{
  LOG("Pythia8Had", pNOTICE) << "Running PYTHIA hadronizer";

  if(!this->AssertValidity(interaction)) {
     LOG("Pythia8Had", pERROR) << "Returning a null particle list!";
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

  LOG("Pythia8Had", pNOTICE)
          << "Hit nucleon pdgc = " << hit_nucleon << ", W = " << W;
  LOG("Pythia8Had", pNOTICE)
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

  int  final_quark = 0; // leading quark (hit quark after the interaction)

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
      LOG("Pythia8Had", pERROR)
        << "Not allowed mode. Refused to make a final quark assignment!";
      return 0;
    }
  }//CC

  //
  // PYTHIA -> HADRONIZATION
  //

#ifdef __GENIE_PYTHIA8_ENABLED__
  LOG("Pythia8Had", pNOTICE)
        << "Fragmentation / Init System: "
        << "Hit q = " << final_quark << ", Final q = " << final_quark;

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
  bool beamConfigExists = fPythia->BeamConfigExists(probe, hit_nucleon);
  fPythia->InitializeBeam(probe, hit_nucleon, eMax);
  Pythia8::Pythia      * pythia8     = fPythia->Pythia8();
  Pythia8::LHAup_Genie * eventReader = fPythia->EventReader();
  Pythia8::GBeamShape  * beamShape   = fPythia->GBeamShape();

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

  bool pi0_decflag = fPythia->Pythia8()->particleData.canDecay(kPdgPi0);
  bool K0_decflag  = fPythia->Pythia8()->particleData.canDecay(kPdgK0);
  bool K0b_decflag = fPythia->Pythia8()->particleData.canDecay(kPdgAntiK0);
  bool L0_decflag  = fPythia->Pythia8()->particleData.canDecay(kPdgLambda);
  bool L0b_decflag = fPythia->Pythia8()->particleData.canDecay(kPdgAntiLambda);
  bool Dm_decflag  = fPythia->Pythia8()->particleData.canDecay(kPdgP33m1232_DeltaM);
  bool D0_decflag  = fPythia->Pythia8()->particleData.canDecay(kPdgP33m1232_Delta0);
  bool Dp_decflag  = fPythia->Pythia8()->particleData.canDecay(kPdgP33m1232_DeltaP);
  bool Dpp_decflag = fPythia->Pythia8()->particleData.canDecay(kPdgP33m1232_DeltaPP);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for pi0           =  " << pi0_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for K0            =  " << K0_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for \bar{K0}      =  " << K0b_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for Lambda        =  " << L0_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for \bar{Lambda0} =  " << L0b_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for D-            =  " << Dm_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for D0            =  " << D0_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for D+            =  " << Dp_decflag;
  LOG("Pythia8Had", pDEBUG) << "Original decay flag for D++           =  " << Dpp_decflag;
#endif

  fPythia->Pythia8()->particleData.mayDecay(kPdgPi0,              false ); // don't decay pi0
  fPythia->Pythia8()->particleData.mayDecay(kPdgK0,               false ); // don't decay K0
  fPythia->Pythia8()->particleData.mayDecay(kPdgAntiK0,           false ); // don't decay \bar{K0}
  fPythia->Pythia8()->particleData.mayDecay(kPdgLambda,           false ); // don't decay Lambda0
  fPythia->Pythia8()->particleData.mayDecay(kPdgAntiLambda,       false ); // don't decay \bar{Lambda0}
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaM,  true  ); // decay Delta-
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_Delta0,  true  ); // decay Delta0
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaP,  true  ); // decay Delta+
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaPP, true  ); // decay Delta++

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
  fPythia->Pythia8()->particleData.mayDecay(kPdgPi0,             pi0_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgK0,              K0_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgAntiK0,          K0b_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgLambda,          L0_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgAntiLambda,      L0b_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaM, Dm_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_Delta0, D0_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaP, Dp_decflag);
  fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaPP,Dpp_decflag);

// get record
  Pythia8::Event &fEvent = pythia8->event;
  int numpart = fEvent.size();
  // assert(numpart>0);
  if (numpart == 0) {
     LOG("PythiaHad", pERROR) << "Returning a null particle list!";
     return 0;
  }
#endif

  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", numpart);
  particle_list->SetOwner(true);

#ifdef __GENIE_PYTHIA8_ENABLED__
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

  // Hadronic 4vec
  TLorentzVector p4Had = kinematics.HadSystP4();

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  for (int i = 1; i < numpart; ++i) {
    /*
     * Convert Pythia8 status code to Pythia6
     * Initial quark has a pythia6 status code of 12
     * The initial quark and the fragmented particles have a pythia6 code
     * of 11 (kIStNucleonTarget)
     * Final state particles have a positive pythia8 code and a pythia6 code of
     * 1 (kIStStableFinalState)
     * The fragmentation products are generated in the hadronic CM frame
     * where the z>0 axis is the \vec{phad} direction. For each particle 
     * returned by the hadronizer:
     * - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
     * - rotate its 3-momentum from LAB' to LAB
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

    TLorentzVector p4o(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
    p4o.Boost(beta); 
    TVector3 p3 = p4o.Vect();
    p3.RotateUz(unitvq); 
    TLorentzVector p4(p3,p4o.Energy());

    // insert the particle in the list
    new((*particle_list)[i]) GHepParticle(
            fEvent[i].id(),     // pdg
            gStatus,            // status
            mother1,            // first parent
            mother2,            // second parent
            daughter1,          // first daughter
            daughter2,          // second daughter
            p4.Px(),            // px [GeV/c]
            p4.Py(),            // py [GeV/c]
            p4.Pz(),            // pz [GeV/c]
            p4.Energy(),        // e  [GeV]
            fEvent[i].xProd(),  // x  [mm]
            fEvent[i].yProd(),  // y  [mm]
            fEvent[i].zProd(),  // z  [mm]
            fEvent[i].tProd()   // t  [mm/c]
    );
  }
#endif

  utils::fragmrec::Print(particle_list);

  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * 
   Pythia8Hadronization::SelectParticles(
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
TH1D * Pythia8Hadronization::MultiplicityProb( const Interaction * interaction ) const
{
// Similar comments apply as in SelectParticles()

  if(!this->AssertValidity(interaction)) {
     LOG("Pythia8Had", pWARN) 
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
    SLOG("Pythia8Had", pWARN) << "probability distribution integral = 0";
    return mult_prob;
  }

  return mult_prob;
}
//____________________________________________________________________________
double Pythia8Hadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
bool Pythia8Hadronization::AssertValidity(const Interaction * interaction) const {

  // check that there is no charm production
  // (GENIE uses a special model for these cases)
  if(interaction->ExclTag().IsCharmEvent()) {
    LOG("Hadronization", pWARN) << "Can't hadronize charm events";
    return false;
  }
  // check the available mass
  double W = utils::kinematics::W(interaction);
  if(W < this->Wmin()) {
    LOG("Hadronization", pWARN) << "Low invariant mass, W = "
				<< W << " GeV!!";
    return false;
  }

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  if( ! target.HitQrkIsSet() ) {
    LOG("Hadronization", pWARN) << "Hit quark was not set!";
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
    LOG("Hadronization", pWARN)
      << "Can only handle electro-weak interactions";
    return false;
  }
  if( !(isp||isn) || !(isv||isvb||isl||islb||isdm) ) {
    LOG("Hadronization", pWARN)
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
    LOG("Hadronization", pWARN)
      << "Impossible interaction type / probe / hit quark combination!";
    return false;
  }

  return true;
}
//____________________________________________________________________________
double Pythia8Hadronization::MaxMult(const Interaction * interaction) const
{
  double W = interaction->Kine().W();

  double maxmult = TMath::Floor(1 + (W-kNeutronMass)/kPionMass);
  return maxmult;
}
//____________________________________________________________________________
TH1D * Pythia8Hadronization::CreateMultProbHist(double maxmult) const
{
  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  TH1D * mult_prob = new TH1D("mult_prob",
			      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
  mult_prob->SetDirectory(0);

  return mult_prob;
}
//____________________________________________________________________________
double Pythia8Hadronization::Wmin(void) const {

  return (kNucleonMass+kPionMass);
}
//____________________________________________________________________________
void Pythia8Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadronization::LoadConfig(void)
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  // the configurable PYTHIA parameters used here are the ones used by NUX 
  // (see A.Rubbia's talk @ NuINT-01)
  // The defaults are the values used by PYTHIA
  // Use the NUX config set to set the tuned values as used in NUX.

  GetParam( "PYTHIA-SSBarSuppression", fSSBarSuppression ) ;
  GetParam( "PYTHIA-GaussianPt2",      fGaussianPt2      ) ;
  GetParam( "PYTHIA-NonGaussianPt2Tail", fNonGaussianPt2Tail  ) ;
  GetParam( "PYTHIA-RemainingEnergyCutoff", fRemainingECutoff ) ;

  GetParam( "PYTHIA-DiQuarkSuppression", fDiQuarkSuppression ) ;
  GetParam( "PYTHIA-LightVMesonSuppression", fLightVMesonSuppression ) ;
  GetParam( "PYTHIA-SVMesonSuppression", fSVMesonSuppression ) ;
  GetParam( "PYTHIA-Lunda", fLunda ) ;
  GetParam( "PYTHIA-Lundb", fLundb ) ;
  GetParam( "PYTHIA-LundaDiq", fLundaDiq ) ;

  fPythia->Pythia8()->settings.parm("StringFlav:probQQtoQ", fDiQuarkSuppression);
  fPythia->Pythia8()->settings.parm("StringFlav:mesonUDvector", fLightVMesonSuppression);
  fPythia->Pythia8()->settings.parm("StringFlav:mesonSvector", fSVMesonSuppression);
  fPythia->Pythia8()->settings.parm("StringZ:aLund", fLunda);
  fPythia->Pythia8()->settings.parm("StringZ:bLund", fLundb);
  fPythia->Pythia8()->settings.parm("StringZ:aExtraDiquark", fLundaDiq);

  fPythia->Pythia8()->settings.parm("StringFlav:probStoUD", fSSBarSuppression);
  fPythia->Pythia8()->settings.parm("Diffraction:primKTwidth", fGaussianPt2);
  fPythia->Pythia8()->settings.parm("StringPT:enhancedFraction", fNonGaussianPt2Tail);
  fPythia->Pythia8()->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  double Wcut, Wmin ;
  GetParam( "Wcut",            Wcut ) ;
  GetParam( "KNO2PYTHIA-Wmin", Wmin ) ;

  if ( Wcut > Wmin ) {
    LOG("Pythia8Hadronization", pFATAL) << "Wcut value too high and in conflict with the KNO2PYTHIA-Wmin" ;
    LOG("Pythia8Hadronization", pFATAL) << "Wcut = " << Wcut ;
    LOG("Pythia8Hadronization", pFATAL) << "KNO2PYTHIA-Wmin = " << Wmin ;
  }

  LOG("Pythia8Hadronization", pDEBUG) << GetConfig() ;
#endif
}
