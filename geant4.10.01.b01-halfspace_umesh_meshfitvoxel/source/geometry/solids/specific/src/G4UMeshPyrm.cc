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
// $Id: G4UMeshPyrm.cc 76263 2013-11-08 11:41:52Z gcosmo $
//
// class G4UMeshPyrm
//
// Implementation for G4UMeshPyrm class
//
// History:
//
//  20140328 - Yuefeng Qiu Created
//
// --------------------------------------------------------------------

#include "G4UMeshPyrm.hh"

const char G4UMeshPyrm::CVSVers[]="$Id: G4UMeshPyrm.cc *Not generated* $";

//#include "G4VoxelLimits.hh"
//#include "G4AffineTransform.hh"

//#include "G4VPVParameterisation.hh"

//#include "Randomize.hh"

//#include "G4VGraphicsScene.hh"
//#include "G4Polyhedron.hh"
//#include "G4VisExtent.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <cmath>

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4UMeshPyrms are created.
// A Tet has all of its geometrical information precomputed

G4UMeshPyrm::G4UMeshPyrm(const G4String& pName,
                 G4ThreeVector anchor,
                 G4ThreeVector p2,
                 G4ThreeVector p3,
                 G4ThreeVector p4,
                 G4ThreeVector p5
                 )
  : G4UMeshElement1st(pName)
{

    //push the nodes into the node list
    fNodeList.push_back(anchor);
    fNodeList.push_back(p2);
    fNodeList.push_back(p3);
    fNodeList.push_back(p4);
    fNodeList.push_back(p5);

    //pre-processing the all necessary data
    preProc();
}



//form the faces using id of nodes in fNodeList
void G4UMeshPyrm::formFaces()
{
    /*
//                                  z
//                                .
//                              ,/
//                             /
//                           4
//                         ,/|`\
//                       ,/ /|  `\
//                     ,/  / '.   `\
//                   ,/  /    |     `\
//                 ,/  /      |       `\
//                0---/-------'.--------3 --> y
//             ,/    /         |      ,/
//            ,/    /          |    ,/
//           ,/   /            '. ,/
//          ,/  /               |/
//          1-/--------'.--------2
//                                `\.
//                                  ` x
    */
        //form the face for the hexahedron
        //node for each face is order in anticlockwise when see from inside
        vector <G4int> Face1, Face2, Face3, Face4, Face5;
        Face1.push_back(0);
        Face1.push_back(1);
        Face1.push_back(2);
        Face1.push_back(3);

        Face2.push_back(0);
        Face2.push_back(4);
        Face2.push_back(1);

        Face3.push_back(1);
        Face3.push_back(4);
        Face3.push_back(2);

        Face4.push_back(3);
        Face4.push_back(2);
        Face4.push_back(4);

        Face5.push_back(0);
        Face5.push_back(3);
        Face5.push_back(4);


        fFaceNodesList.push_back(Face1);
        fFaceNodesList.push_back(Face2);
        fFaceNodesList.push_back(Face3);
        fFaceNodesList.push_back(Face4);
        fFaceNodesList.push_back(Face5);
}

//calculate the volume size and assign to fCubicVolume
void  G4UMeshPyrm::calVolume()
{
    // decomposition hexa into tetra
    const  int nbTet = 2;
    const static int vtab[nbTet][4] = {
      // hexahedron
        { 0, 1, 3, 4 },
        { 1, 2, 3, 4 },};

    double theVolume = 0.0;
    FoundNagtvVol = false;
    for (int i = 0; i <nbTet  ; i++) {
        G4double aVolume =  calTetraVolume( fNodeList[ vtab[i][0] ],
                fNodeList[ vtab[i][1] ],
                fNodeList[ vtab[i][2] ],
                fNodeList[ vtab[i][3] ]);
        if (aVolume < 0) FoundNagtvVol = true;
      theVolume +=aVolume;
    }
    fCubicVolume = theVolume ;
}

//reverse the orientation of the nodes, try to make volume to positive
void G4UMeshPyrm::revOrientation()
{
    //swap node 1 <-> 3
    G4ThreeVector aBuffer = fNodeList[1];
    fNodeList[1] = fNodeList[3];
    fNodeList[3] = aBuffer;
}



//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UMeshPyrm::~G4UMeshPyrm()
{
}



//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4UMeshPyrm::GetEntityType() const
{
  return G4String("G4UMeshPyrm");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UMeshPyrm::Clone() const
{
  return new G4UMeshPyrm(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4UMeshPyrm::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
  << "    *** Dump for solid - " << GetName() << " ***\n"
  << "    ===================================================\n"
  << " Solid type: G4UMeshPyrm\n"
  << " Parameters: \n" ;
  for (unsigned int i=0; i<fNodeList.size(); i++) {
      os << "p"<<i<<"\t: "<<fNodeList[i]/mm<<" mm \n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}