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
// * and NASA under contract number NNG04CT05P                        *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// $Id: G4UMeshTet.cc 76263 2013-11-08 11:41:52Z gcosmo $
//
// class G4UMeshTet
//
// Implementation for G4UMeshTet class
//
// History:
//
//  20140328 - Yuefeng Qiu Created
//
// --------------------------------------------------------------------

#include "G4UMeshTet.hh"

const char G4UMeshTet::CVSVers[]="$Id: G4UMeshTet.cc *Not generated* $";

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"

#include "G4ThreeVector.hh"

#include <cmath>

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4UMeshTets are created.
// A Tet has all of its geometrical information precomputed

G4UMeshTet::G4UMeshTet(const G4String& pName,
                 G4ThreeVector anchor,
                 G4ThreeVector p2,
                 G4ThreeVector p3,
                 G4ThreeVector p4)
  : G4UMeshElement1st(pName)
{

    //push the nodes into the node list
    fNodeList.push_back(anchor);
    fNodeList.push_back(p2);
    fNodeList.push_back(p3);
    fNodeList.push_back(p4);

    //pre-processing the all necessary data
    preProc();
}



//form the faces using id of nodes in fNodeList
void G4UMeshTet::formFaces()
{
//                                  z
//                                .
//                              ,/
//                             /
//                           3
//                         ,/|`\
//                       ,/  |  `\
//                     ,/    '.   `\
//                   ,/       |     `\
//                 ,/         |       `\
//                0-----------'.--------2 --> y
//                 `\.         |      ,/
//                    `\.      |    ,/
//                       `\.   '. ,/
//                          `\. |/
//                             `1
//                                `\.
//                                  ` x


        //form the face for the hexahedron
        //node for each face is order in anticlockwise when see from inside
        vector <G4int> Face1, Face2, Face3, Face4;
        Face1.push_back(0);
        Face1.push_back(1);
        Face1.push_back(2);

        Face2.push_back(0);
        Face2.push_back(2);
        Face2.push_back(3);

        Face3.push_back(0);
        Face3.push_back(3);
        Face3.push_back(1);

        Face4.push_back(1);
        Face4.push_back(3);
        Face4.push_back(2);

        fFaceNodesList.push_back(Face1);
        fFaceNodesList.push_back(Face2);
        fFaceNodesList.push_back(Face3);
        fFaceNodesList.push_back(Face4);
}

//calculate the volume size and assign to fCubicVolume
void  G4UMeshTet::calVolume()
{
    FoundNagtvVol = false;
    G4double theVolume =  calTetraVolume( fNodeList[0],
            fNodeList[1],
            fNodeList[2],
            fNodeList[3]);
    if (theVolume < 0) FoundNagtvVol = true;

    fCubicVolume = theVolume ;
}

//reverse the orientation of the nodes, try to make volume to positive
void G4UMeshTet::revOrientation()
{
    //swap node 1 <-> 2
    G4ThreeVector aBuffer = fNodeList[1];
    fNodeList[1] = fNodeList[2];
    fNodeList[2] = aBuffer;

}



//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UMeshTet::~G4UMeshTet()
{
}



//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4UMeshTet::GetEntityType() const
{
  return G4String("G4UMeshTet");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UMeshTet::Clone() const
{
  return new G4UMeshTet(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4UMeshTet::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
  << "    *** Dump for solid - " << GetName() << " ***\n"
  << "    ===================================================\n"
  << " Solid type: G4UMeshTet\n"
  << " Parameters: \n" ;
  for (unsigned int i=0; i<fNodeList.size(); i++) {
      os << "p"<<i<<"\t: "<<fNodeList[i]/mm<<" mm \n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}
