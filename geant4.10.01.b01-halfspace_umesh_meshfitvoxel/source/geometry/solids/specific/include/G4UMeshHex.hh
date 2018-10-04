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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P.                       *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
//
// $Id: G4UMeshHex.hh 76263 2013-11-08 11:41:52Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4UMeshHex
//
// Class description:
//
//   A G4UMeshHex is a hexahedron
//

// History:
// -------
// 28.03.2014 - Y.Qiu(Karlsruhe Institute of Technology)
// --------------------------------------------------------------------
#ifndef G4UMESHHEX_HH
#define G4UMESHHEX_HH

#include "G4UMeshElement1st.hh"

using namespace std;

class G4UMeshHex : public G4UMeshElement1st
{

  public:  // with description

    G4UMeshHex(const G4String& pName,
          G4ThreeVector anchor,
          G4ThreeVector p2,
          G4ThreeVector p3,
          G4ThreeVector p4,
          G4ThreeVector p5,
          G4ThreeVector p6,
          G4ThreeVector p7,
          G4ThreeVector p8
          );

    virtual ~G4UMeshHex();

    G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;


  public:   // without description

    const char* CVSHeaderVers()
      { return "$Id: G4UMeshHex.hh 76263 2013-11-08 11:41:52Z gcosmo $"; }
    const char* CVSFileVers()
      { return CVSVers; }


  private:
    static const char CVSVers[];


    bool            checkPosVolume() {return true; /*Not implemented yet*/};
    //check if the volume is positive. Not imiplemented
    void            formFaces();
    //form the faces using id of nodes in fNodeList
    void            calVolume();
    //calculate the volume size and assign to fCubicVolume
    void            revOrientation();
    //reverse the orientation of the nodes, try to make volume to positive


};

#endif

