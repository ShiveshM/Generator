//____________________________________________________________________________
/*!

\class   genie::ROOTGeomAnalyzer

\brief   A ROOT/GEANT Geometry Analyzer

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created May 24, 2005

*/
//____________________________________________________________________________

#ifndef _ROOT_GEOMETRY_ANALYZER_H_
#define _ROOT_GEOMETRY_ANALYZER_H_

#include <string>

#include <TGeoManager.h>

#include "EVGDrivers/GeomAnalyzerI.h"
#include "PDG/PDGUtils.h"

class TGeoVolume;
class TGeoMaterial;
class TGeoElement;

using std::string;

namespace genie {

class ROOTGeomAnalyzer : public GeomAnalyzerI {

public :

  ROOTGeomAnalyzer(string filename);
 ~ROOTGeomAnalyzer();

  // analyzer configuration options
  void SetScannerNPoints(int np) { fNPoints = np; };
  void SetScannerNRays  (int nr) { fNRays   = nr; };

  // implement the GeomAnalyzerI interface

  const PDGCodeList &    ListOfTargetNuclei    (void);
  const PathLengthList & ComputeMaxPathLengths (void);

  const PathLengthList &
           ComputePathLengths
             (const TLorentzVector & x, const TLorentzVector & p);

  const TVector3 &
           GenerateVertex
             (const TLorentzVector & x, const TLorentzVector & p, int tgtpdg);

private:

  void   Initialize              (string filename);
  void   BuildListOfTargetNuclei (void);
  int    GetTargetPdgCode        (const TGeoMaterial * const m) const;
  int    GetTargetPdgCode        (const TGeoElement  * const e) const;
  double ComputeMaxPathLengthPDG (double* XYZ, double* direction, int pdgc);


  int              fMaterial;               ///< [input] selected material for vertex
  TGeoManager *    fGeometry;               ///< [input] detector geometry
  int              fNPoints;                ///< max path length scanner: points/surface [def:200]
  int              fNRays;                  ///< max path length scanner: rays/point [def:200]
  TVector3 *       fCurrVertex;             ///< current generated vertex
  PathLengthList * fCurrPathLengthList;     ///< current list of path-lengths
  PathLengthList * fCurrMaxPathLengthList;  ///< current list of Max path-lengths
  PDGCodeList *    fCurrPDGCodeList;        ///< current list of target nuclei
};

}      // genie namespace

#endif // _ROOT_GEOMETRY_ANALYZER_H_
