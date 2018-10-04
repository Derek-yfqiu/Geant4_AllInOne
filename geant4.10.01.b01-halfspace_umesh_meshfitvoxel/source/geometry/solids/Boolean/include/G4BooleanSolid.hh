//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4BooleanSolid.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4BooleanSolid
//
// Class description:
//
// Abstract base class for solids created by boolean operations
// between other solids.

// History:
//
// 10.09.98 V.Grichine, created
//
// --------------------------------------------------------------------
#ifndef G4BOOLEANSOLID_HH
#define G4BOOLEANSOLID_HH

#include "G4DisplacedSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

class HepPolyhedronProcessor;

class G4BooleanSolid : public G4VSolid
{
  public:  // with description
 
    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB   );

    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB,
                          G4RotationMatrix* rotMatrix,
                    const G4ThreeVector& transVector    );

    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB , 
                    const G4Transform3D& transform   );

    virtual ~G4BooleanSolid();

    virtual const G4VSolid* GetConstituentSolid(G4int no) const;
    virtual       G4VSolid* GetConstituentSolid(G4int no);
      // If Solid is made up from a Boolean operation of two solids,
      // return the corresponding solid (for no=0 and 1).
      // If the solid is not a "Boolean", return 0.

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    virtual G4GeometryType  GetEntityType() const;
    virtual G4Polyhedron* GetPolyhedron () const;

    std::ostream& StreamInfo(std::ostream& os) const;

    inline G4int GetCubVolStatistics() const;
    inline G4double GetCubVolEpsilon() const;
    inline void SetCubVolStatistics(G4int st);
    inline void SetCubVolEpsilon(G4double ep);

    inline G4int GetAreaStatistics() const;
    inline G4double GetAreaAccuracy() const;
    inline void SetAreaStatistics(G4int st);
    inline void SetAreaAccuracy(G4double ep);

    G4ThreeVector GetPointOnSurface() const;
    //qiu
    inline void    SetPolyhedron(G4Polyhedron* aPolyhedron) {fpPolyhedron = aPolyhedron;};
    //set the polyhedron representing this solid;
    //ATTENTION:because it is difficult to generate the polyhedron by this code,
    //for displaying this solid the polyhedron should be set
    //this solid takes the ownership of this polyhedron

  public:  // without description

    G4BooleanSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4BooleanSolid(const G4BooleanSolid& rhs);
    G4BooleanSolid& operator=(const G4BooleanSolid& rhs);
      // Copy constructor and assignment operator.

  protected:
  
    G4Polyhedron* StackPolyhedron(HepPolyhedronProcessor&,
                                  const G4VSolid*) const;
      // Stack polyhedra for processing. Return top polyhedron.

    inline G4double GetAreaRatio() const;
      // Ratio of surface areas of SolidA to total A+B

  protected:
  
    G4VSolid* fPtrSolidA;
    G4VSolid* fPtrSolidB;

    mutable G4double fAreaRatio; // Calculation deferred to GetPointOnSurface()

  private:

    G4int    fStatistics;
    G4double fCubVolEpsilon;
    G4double fAreaAccuracy;
    G4double fCubicVolume;
    G4double fSurfaceArea;

    mutable G4Polyhedron* fpPolyhedron;

    G4bool  createdDisplacedSolid;
      // If & only if this object created it, it must delete it

} ;

#include "G4BooleanSolid.icc"

#endif
