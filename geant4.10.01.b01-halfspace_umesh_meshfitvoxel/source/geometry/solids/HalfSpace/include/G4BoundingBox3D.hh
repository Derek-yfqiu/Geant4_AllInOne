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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4BoundingBox3D
//
// Class description:
// 
// Definition of a generic solid's bounding box in the 3D space.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BoundingBox3D_h
#define __G4BoundingBox3D_h 1

//#include "G4Ray.hh"
#include "G4ThreeVector.hh"
//#include "G4ThreeVector.hh"

class G4BoundingBox3D
{

public:  // with description

  G4BoundingBox3D();
  G4BoundingBox3D(const G4ThreeVector&);
  G4BoundingBox3D(const G4ThreeVector&, const G4ThreeVector&);
  ~G4BoundingBox3D();
    // Constructors & destructor.

  G4BoundingBox3D(const G4BoundingBox3D& right);
  G4BoundingBox3D& operator=(const G4BoundingBox3D& right);
    // Copy constructor and assignment operator.

  void Init(const G4ThreeVector&);
  void Init(const G4ThreeVector&, const G4ThreeVector&);
  void Extend(const G4ThreeVector&);
  void Margin (const G4double & aMargin);
    // To create/extend the bounding box

  inline G4ThreeVector GetBoxMin() const;
  inline G4ThreeVector GetBoxMax() const;
  inline G4double GetDistance() const;
  inline void SetDistance(G4double distance0);
    // Accessors.

  G4int Inside(const G4ThreeVector&) const;
    // Returns 1 if the point is inside and on the bbox.
    // Returns 0 if the point is outside the bbox.

  inline G4ThreeVector GetMiddlePoint() const;
  inline G4double      GetSize() const;
  inline G4double      GetDX() const;
  inline G4double      GetDY() const;
  inline G4double      GetDZ() const;


public:

  inline G4int GetTestResult() const;
  G4int Test(const G4ThreeVector& aPoint, const G4ThreeVector& aVector);

  static const G4BoundingBox3D space;

private:

  G4int BoxIntersect(const G4ThreeVector&,
             const G4ThreeVector&,
             const G4ThreeVector&) const;

  G4double DistanceToIn(const G4ThreeVector&,
            const G4ThreeVector&) const;
    
private:

  G4ThreeVector box_min;
  G4ThreeVector box_max;
  G4double distance;

  G4int test_result;

  G4ThreeVector MiddlePoint;
  G4double      size; //the radius for the bounding sphere
  G4ThreeVector GeantBox;
  G4double kCarTolerance;
  G4double      Dx, Dy, Dz; //Delta dimension
};

#include "G4BoundingBox3D.icc"

#endif
