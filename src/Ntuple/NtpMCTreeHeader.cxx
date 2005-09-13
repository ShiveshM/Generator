//____________________________________________________________________________
/*!

\class   genie::NtpMCTreeHeader

\brief   MINOS-style Ntuple Class to hold an output MC Tree Header

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TObjString.h>

#include "Ntuple/NtpMCTreeHeader.h"

using std::string;
using namespace genie;

ClassImp(NtpMCTreeHeader)

using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCTreeHeader & hdr)
  {
     hdr.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCTreeHeader::NtpMCTreeHeader() :
TNamed("header","GENIE output tree header")
{
  this->Init();
}
//____________________________________________________________________________
NtpMCTreeHeader::NtpMCTreeHeader(const NtpMCTreeHeader & hdr) :
TNamed("header","GENIE output tree header")
{
  this->Copy(hdr);
}
//____________________________________________________________________________
NtpMCTreeHeader::~NtpMCTreeHeader()
{

}
//____________________________________________________________________________
void NtpMCTreeHeader::Fill(NtpMCFormat_t fmt)
{
  format =  fmt;

  //-- read/store all GENIE environmental variables that have been set
  string envvar = "";
  envvar = (gSystem->Getenv("GEVGL")    ? gSystem->Getenv("GEVGL")    : "");
  env.gevgl.SetString(envvar.c_str());
  envvar = (gSystem->Getenv("GSPLOAD")  ? gSystem->Getenv("GSPLOAD")  : "");
  env.gspload.SetString(envvar.c_str());
  envvar = (gSystem->Getenv("GSPSAVE")  ? gSystem->Getenv("GSPSAVE")  : "");
  env.gspsave.SetString(envvar.c_str());
  envvar = (gSystem->Getenv("GMSGCONF") ? gSystem->Getenv("GMSGCONF") : "");
  env.gmsgconf.SetString(envvar.c_str());

  //-- need to add code to save the MC job configuration from AlgConfigPool
}
//____________________________________________________________________________
void NtpMCTreeHeader::PrintToStream(ostream & stream) const
{
  stream << "NtpRecord format: "
                 << NtpMCFormat::AsString(this->format) << endl;
/*
  stream << "Stored MC job configuration: " << endl;
  TIter listiter(&this->config);
  TObjString * entry = 0;
  while( (entry = (TObjString *)listiter.Next()) ) {
     stream << entry->GetString().Data() << endl;
  }
*/
}
//____________________________________________________________________________
void NtpMCTreeHeader::Copy(const NtpMCTreeHeader & hdr)
{
  this->format = hdr.format;
/*
  this->config.Delete();
  TIter listiter(&hdr.config);
  TObjString * entry = 0;
  while( (entry = (TObjString *)listiter.Next()) ) {
    this->config.Add(entry);
  }
*/
}
//____________________________________________________________________________
void NtpMCTreeHeader::Init(void)
{
  this->format = kNFUndefined;
/*
  this->config.Delete();
*/
}
//____________________________________________________________________________

