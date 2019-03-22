//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 21, 2009 - CA
   Remove SyncSeeds() function. Now the GENIE/PYTHIA8 random number generator
   seed is synchonized at the genie::RandomGen() initialization.
 @ Oct 02, 2009 - CA
   Re-organize code and implement the `UnInhibitDecay(int,TDecayChannel*)
   const' and `InhibitDecay(int,TDecayChannel*) const' methods.
   Test/fix the code to match a ROOT TDecayChannel to a PYTHIA8 decay channel.
   Decay() returns null if decay is inhibited or if the sum{branching ratios}
   for all enabled decay channels is non-positive. In case of inhibited decay
   channels, a weight is calculated as w = 1./sum{BR for enabled channels}.
 @ Feb 04, 2010 - CA
   Comment out (unused) code using the fForceDecay flag

*/
//____________________________________________________________________________

#include <vector>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TDecayChannel.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/Decay/PythiaDecayer.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"

using std::vector;

using namespace genie;

//____________________________________________________________________________
PythiaDecayer::PythiaDecayer() :
DecayModelI("genie::PythiaDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::PythiaDecayer(string config) :
DecayModelI("genie::PythiaDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::~PythiaDecayer() 
{ 
}
//____________________________________________________________________________
bool PythiaDecayer::IsHandled(int code) const
{
// does not handle requests to decay baryon resonances
  
  if( utils::res::IsBaryonResonance(code) ) {
     LOG("PythiaDec", pINFO)
       << "This algorithm can not decay particles with PDG code = " << code;
     return false;
  } else return true;
}
//____________________________________________________________________________
void PythiaDecayer::Initialize(void) const
{
  fWeight = 1.;

  fPythia8 = PythiaSingleton::Instance();
  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  pythia8->settings.resetAll();

  pythia8->readString("ProcessLevel:all = off");
  pythia8->readString("Print:quiet      = on");

  // sync GENIE/PYTHIA8 seeds
  RandomGen::Instance();
  pythia8->init();
}
//____________________________________________________________________________
TClonesArray * PythiaDecayer::Decay(const DecayerInputs_t & inp) const
{
  fWeight = 1.; // reset weight

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();

  int pdgc = inp.PdgCode;

  if ( ! this->IsHandled(pdgc) ) return 0;
  
  bool md = pythia8->particleData.canDecay(pdgc);
  if(not md) {
    LOG("PythiaDec", pNOTICE)
       << (PDGLibrary::Instance())->Find(pdgc)->GetName() 
       << " decays are inhibited!";
    return 0;
  }

  double sumbr = this->SumBR(pdgc);
  if(sumbr <= 0) {
    LOG("PythiaDec", pNOTICE)
       << "The sum of enabled "
       << (PDGLibrary::Instance())->Find(pdgc)->GetName() 
       << " decay channel branching rations is non-positive!";
    return 0;
  }

  fWeight = 1./sumbr; // update weight to account for inhibited channels

  double E  = inp.P4->Energy();
  double M  = pythia8->particleData.m0(pdgc);
  double pz = Pythia8::sqrtpos(E*E - M*M);

  pythia8->event.reset();

  pythia8->event.append(pdgc, 11, 0, 0, 0., 0., pz, E, M);
  pythia8->next();

  // List event information
  pythia8->event.list();
  pythia8->stat();

  //-- get decay products
  Pythia8::Event &fEvent = pythia8->event;
  int numpart = fEvent.size();
   int ioff = 0;
  if (fEvent[0].id() == 90) ioff = -1;
   TClonesArray * pl = new TClonesArray("genie::GHepParticle", numpart);
  if(!pl) return 0;
   for (int i = 1; i < numpart; ++i) {
    /*
     * Convert Pythia8 status code to Pythia6
     * Decayed/fragmented particle has a pytahi6 code of 11 (kIStNucleonTarget)
     * Final state particles have a negative pythia8 code and a pythia6 code of 1 (kIStStableFinalState)
     */
    GHepStatus_t gStatus = (fEvent[i].status()>0) ? kIStStableFinalState : kIStNucleonTarget;
    new((*pl)[i]) GHepParticle(
            fEvent[i].id(),
            gStatus,
            fEvent[i].mother1()   + ioff,
            fEvent[i].mother2()   + ioff,
            fEvent[i].daughter1() + ioff,
            fEvent[i].daughter2() + ioff,
            fEvent[i].px(),       // [GeV/c]
            fEvent[i].py(),       // [GeV/c]
            fEvent[i].pz(),       // [GeV/c]
            fEvent[i].e(),        // [GeV]
            fEvent[i].xProd(),    // [mm]
            fEvent[i].yProd(),    // [mm]
            fEvent[i].zProd(),    // [mm]
            fEvent[i].tProd());   // [mm/c]
  }

  //-- transfer ownership and return
  pl->SetOwner(true);    
  return pl;
}
//____________________________________________________________________________
double PythiaDecayer::Weight(void) const 
{
  return fWeight; 
}
//____________________________________________________________________________
void PythiaDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return; 

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();

  if(!dc) {
    LOG("PythiaDec", pINFO)
       << "Switching OFF ALL decay channels for particle = " << pdgc;
    pythia8->particleData.mayDecay(pdgc, false);
    return;
  }

  LOG("PythiaDec", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdgc;

  int ichannel = this->FindPythiaDecayChannel(pdgc, dc);
  if(ichannel != -1) {
    Pythia8::ParticleDataEntry * fPDE = pythia8->particleData.particleDataEntryPtr(pdgc);
    fPDE->channel(ichannel).onMode(0); // switch-off
  }
}
//____________________________________________________________________________
void PythiaDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return; 

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  Pythia8::ParticleDataEntry * fPDE = pythia8->particleData.particleDataEntryPtr(pdgc);

  if(!dc) {
    LOG("PythiaDec", pINFO)
      << "Switching ON all PYTHIA decay channels for particle = " << pdgc;

    pythia8->particleData.mayDecay(pdgc, true);

    // loop over pythia decay channels
    int size_channels = fPDE->sizeChannels();
    for (int ichannel = 0; ichannel < size_channels; ++ichannel) {
      fPDE->channel(ichannel).onMode(1); // switch-on
    }
    return;
  }//!dc

  LOG("PythiaDec", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdgc;

  int ichannel = this->FindPythiaDecayChannel(pdgc, dc);
  if(ichannel != -1) {
    fPDE->channel(ichannel).onMode(1); // switch-on
  }
}
//____________________________________________________________________________
double PythiaDecayer::SumBR(int pdgc) const
{
// Sum of branching ratios for enabled channels
//
  double sumbr=0.;

  bool has_inhibited_channels=false;

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  Pythia8::ParticleDataEntry * fPDE = pythia8->particleData.particleDataEntryPtr(pdgc);
  int size_channels = fPDE->sizeChannels();

  // loop over pythia decay channels
  for (int ichannel = 0; ichannel < size_channels; ++ichannel) {
     bool enabled = (fPDE->channel(ichannel).onMode() == 1);
     if (!enabled) { 
       has_inhibited_channels = true; 
     } else {
       sumbr += fPDE->channel(ichannel).bRatio();
     }

/*
     LOG("PythiaDec", pDEBUG) 
        << "ich: " << ichannel << ", " << ((enabled)? "ON" : "OFF")
        << ", br = " << fPDE->channel(ichannel).bRatio()
        << ", curr_sum{br}[enabled] = " << sumbr;
*/
  }

  if(!has_inhibited_channels) return 1.;

  LOG("PythiaDec", pINFO) << "Sum{BR} = " << sumbr;

  return sumbr;
}
//____________________________________________________________________________
int PythiaDecayer::FindPythiaDecayChannel(int pdgc, TDecayChannel* dc) const
{	
  if(!dc) return -1;

  bool found_match = false;

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  Pythia8::ParticleDataEntry * fPDE = pythia8->particleData.particleDataEntryPtr(pdgc);
  int size_channels = fPDE->sizeChannels();

  // loop over pythia decay channels
  for (int ichannel = 0; ichannel < size_channels; ++ichannel) {

     // does the  current pythia channel match the input TDecayChannel?
     LOG("PythiaDec", pINFO)
         << "\nComparing PYTHIA's channel = " << ichannel
         << " with TDecayChannel = " << dc->Number();
 
     found_match = this->MatchDecayChannels(pdgc, ichannel, dc);
     if(found_match) {
         LOG("PythiaDec", pNOTICE)
            << " ** TDecayChannel id = " << dc->Number()
            << " corresponds to PYTHIA8 channel id = " << ichannel;
         return ichannel;
     }//match?
  }//loop pythia decay ch.

  LOG("PythiaDec", pWARN)
     << " ** No PYTHIA8 decay channel match found for "
     << "TDecayChannel id = " << dc->Number();

  return -1;
}
//____________________________________________________________________________
bool PythiaDecayer::MatchDecayChannels(int pdgc, int ichannel, TDecayChannel* dc) const
{
  // num. of daughters in the input TDecayChannel & the input PYTHIA ichannel
  int nd = dc->NDaughters();

  Pythia8::Pythia * pythia8 = fPythia8->Pythia8();
  Pythia8::ParticleDataEntry * fPDE = pythia8->particleData.particleDataEntryPtr(pdgc);

  int py_nd = 0;
  for (int i = 0; i < 8; i++) {
     if(fPDE->channel(ichannel).product(i) != 0) py_nd++;
  }

  LOG("PythiaDec", pDEBUG)
    << "NDaughters: PYTHIA = " << py_nd << ", ROOT's TDecayChannel = " << nd;

  if(nd != py_nd) return false;
  
  //
  // if the two channels have the same num. of daughters, then compare them
  //
 
  // store decay daughters for the input TDecayChannel
  vector<int> dc_daughter(nd); 
  for (int i = 0; i < nd; i++) {
     dc_daughter[i] = dc->DaughterPdgCode(i);
  }
  // store decay daughters for the input PYTHIA's ichannel
  vector<int> py_daughter(nd); 
  for (int i = 0; i < py_nd; i++) {
      py_daughter[i] = fPDE->channel(ichannel).product(i);
  }

  // sort both daughter lists 
  sort( dc_daughter.begin(), dc_daughter.end() ); 
  sort( py_daughter.begin(), py_daughter.end() );

  // compare  
  for(int i = 0; i < nd; i++) {
    if(dc_daughter[i] != py_daughter[i]) return false;
  }
  return true;
}
//____________________________________________________________________________
void PythiaDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PythiaDecayer::LoadConfig(void)
{
// Read configuration options or set defaults
/*
  // check whether we are asked to force the decay / default = false
  fForceDecay = fConfig->GetBoolDef("ForceDecay", false);
*/
}
//____________________________________________________________________________
