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
// $Id: G4Hex.cc 76263 2013-11-08 11:41:52Z gcosmo $
//
// class G4Hex
//
// Implementation for G4Hex class
//
// History:
//
//  20040903 - Marcus Mendenhall, created G4Hex
//  20041101 - Marcus Mendenhall, optimized constant dot products with
//             fCdotNijk values
//  20041101 - MHM removed tracking error by clipping DistanceToOut to 0
//             for surface cases
//  20041101 - MHM many speed optimizations in if statements
//  20041101 - MHM changed vdotn comparisons to 1e-12 instead of 0.0 to
//             avoid nearly-parallel problems
//  20041102 - MHM Added extra distance into solid to DistanceToIn(p,v)
//             hit testing
//  20041102 - MHM added ability to check for degeneracy without throwing
//             G4Exception
//  20041103 - MHM removed many unused variables from class
//  20040803 - Dionysios Anninos, added GetPointOnSurface() method
//  20061112 - MHM added code for G4VSolid GetSurfaceArea()
//  20100920 - Gabriele Cosmo added copy-ctor and operator=()
//
// --------------------------------------------------------------------

#include "G4Hex.hh"

const char G4Hex::CVSVers[]="$Id: G4Hex.cc *Not generated* $";

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
// in which a very large number of G4Hexs are created.
// A Tet has all of its geometrical information precomputed

G4Hex::G4Hex(const G4String& pName,
                 G4ThreeVector anchor,
                 G4ThreeVector p2,
                 G4ThreeVector p3,
                 G4ThreeVector p4,
                 G4ThreeVector p5,
                 G4ThreeVector p6,
                 G4ThreeVector p7,
                 G4ThreeVector p8
                 )
  : G4VSolid(pName), fpPolyhedron(0)
{

    //push the nodes into the node list
    fNodeList.push_back(anchor);
    fNodeList.push_back(p2);
    fNodeList.push_back(p3);
    fNodeList.push_back(p4);
    fNodeList.push_back(p5);
    fNodeList.push_back(p6);
    fNodeList.push_back(p7);
    fNodeList.push_back(p8);

    //pre-processing the all necessary data
    preProc();
}

//////////////////////////////////////////////////////////////////////////
//
// preprocessing the solid
void  G4Hex::preProc()
{
    //check if the volume positive
//    if (!checkPosVolume()) {
//        G4Exception("G4Hex::preProc()", "GeomSolids0002", FatalException,
//                    "Volume is negative!.");
//        return;
//    }

    //calculate the boundary box
    fBoundaryBox.calBdBox(fNodeList);

    //form the faces
    formFaces();

    //checking volume if found negative
    calVolume();
    if (FoundNagtvVol) {
        revOrientation(); //try to reverse orientation
        calVolume(); //calculate the volume again
        if (FoundNagtvVol) {
            G4Exception("G4Hex::preProc()", "GeomSolids0002", FatalException,
                        "Volume is negative, and cannot fix by reversing orientation!.");
            return;
        }
    }


    //calculate the face normal
    calFaceNormal();
    //calculat the face normal dot face center
    calFaceNdotCenter();

    //check the if 4 point are on the same surface within tolerance
    fTol = kCarTolerance;
    if (!checkCreaseFace()) {
        G4cerr << "WARNING: found an elment with crease face! "<<G4endl;
    }

    //calculate the volume and area
//    calVolume();
    calArea();

    //calculate the bary center
    calCenter();

    //calculat the bounding sphere radius (fMaxSize)
    calBoundingSphereRadius();
}

//form the faces using id of nodes in fNodeList
void G4Hex::formFaces()
{
    /*
    //    ^ z
    //    |
    //    4----------7
    //    |\         |\
    //    | \    f2  | \
    //    |  \  f5   |  \
    //    |   5------+---6
    //    | f3|      | f4|
    //    0---+------3 - | ---> y
    //     \  |    f6 \  |
    //      \ |   f1   \ |
    //       \|         \|
    //        1----------2
    //         \
    //          x
    */
        //form the face for the hexahedron
        //node for each face is order in anticlockwise when see from inside
        vector <G4int> Face1, Face2, Face3, Face4, Face5,Face6;
        Face1.push_back(0);
        Face1.push_back(1);
        Face1.push_back(2);
        Face1.push_back(3);

        Face2.push_back(4);
        Face2.push_back(7);
        Face2.push_back(6);
        Face2.push_back(5);

        Face3.push_back(0);
        Face3.push_back(4);
        Face3.push_back(5);
        Face3.push_back(1);

        Face4.push_back(3);
        Face4.push_back(2);
        Face4.push_back(6);
        Face4.push_back(7);

        Face5.push_back(0);
        Face5.push_back(3);
        Face5.push_back(7);
        Face5.push_back(4);

        Face6.push_back(2);
        Face6.push_back(1);
        Face6.push_back(5);
        Face6.push_back(6);

        fFaceNodesList.push_back(Face1);
        fFaceNodesList.push_back(Face2);
        fFaceNodesList.push_back(Face3);
        fFaceNodesList.push_back(Face4);
        fFaceNodesList.push_back(Face5);
        fFaceNodesList.push_back(Face6);
}


//calculate face normal
//should after fFaceNode is set
void G4Hex::calFaceNormal()
{
    //we just pickup two edge and make the cross produce
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {
        vector <G4int> aFaceNodes = fFaceNodesList[i];
        if (aFaceNodes.size() < 3) {
            G4Exception("G4Hex::calFaceNormal()", "GeomSolids0002", FatalException,
                        "Face is formed at least by three points.");
            return;
        }
        //two edge which pointing to the same node
        G4ThreeVector aFirstEdge = fNodeList[aFaceNodes[1]] - fNodeList[aFaceNodes[0]]; //vector 0->1
        G4ThreeVector aSecondEdge = fNodeList[aFaceNodes[1]] - fNodeList[aFaceNodes[2]];//vector  2->1
        G4ThreeVector aNormal = aFirstEdge.cross(aSecondEdge);  //cross product, a normal is facing outside
        fFaceNormalList.push_back(aNormal.unit());
    }
}

//calculate the face normal dot face center point
//should after face normal caculation
void   G4Hex::calFaceNdotCenter()
{
    if (fFaceNodesList.size() != fFaceNormalList.size()) {
        G4Exception("G4Hex::calFaceNdotCenter()", "GeomSolids0002", FatalException,
                    "Face normal is not calculated");
        return;
    }

    //calculate the face centers
    G4ThreeVectorList aCenterList;
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {
        vector <G4int> aFaceNodes = fFaceNodesList[i];
        G4ThreeVector aCenter = G4ThreeVector(); //initial with 0,0,0
        //sum all nodes and average
        for (unsigned int j=0; j<aFaceNodes.size(); j++)
            aCenter += fNodeList[aFaceNodes[j]];
        aCenter /= G4double(aFaceNodes.size());
        aCenterList.push_back(aCenter);
    }
    //calculate the center dot normal
    for (unsigned int i=0; i<fFaceNormalList.size() ; i++)
        fCdotNList.push_back(aCenterList[i].dot(fFaceNormalList[i]));
}

//calculate the distance of a point to the face indexed by the a id
G4double  G4Hex::calDistanceToFace(const G4ThreeVector & aPoint, G4int aFaceIdx) const
{
    //we don't check the normal list and other list any more
    //calculate the distance to the face
    return aPoint.dot(fFaceNormalList[aFaceIdx])-fCdotNList[aFaceIdx];
}
//calculate the volume size and assign to fCubicVolume
void  G4Hex::calVolume()
{
    // decomposition hexa into tetra
    const static int vtab[5][4] = {
      // hexahedron
      { 1, 4, 3, 0 },
      { 4, 1, 6, 5 },
      { 1, 3, 6, 2 },
      { 4, 6, 3, 7 },
      { 1, 4, 6, 3 }};

    double theVolume = 0.0;
    FoundNagtvVol = false;
    for (int i = 0; i <5  ; i++) {
        G4double aVolume =  calTetraVolume( fNodeList[ vtab[i][0] ],
                fNodeList[ vtab[i][1] ],
                fNodeList[ vtab[i][2] ],
                fNodeList[ vtab[i][3] ]);
        if (aVolume < 0) FoundNagtvVol = true;
      theVolume +=aVolume;
    }
    fCubicVolume = theVolume ;
}

//calculate the volume of a tetrahedron
double G4Hex::calTetraVolume(const G4ThreeVector & n1,
                             const G4ThreeVector & n2,
                             const G4ThreeVector & n3,
                             const G4ThreeVector & n4)
{
  G4double x1 = n1.x();
  G4double y1 = n1.y();
  G4double z1 = n1.z();

  G4double x2 = n2.x();
  G4double y2 = n2.y();
  G4double z2 = n2.z();

  G4double x3 = n3.x();
  G4double y3 = n3.y();
  G4double z3 = n3.z();

  G4double x4 = n4.x();
  G4double y4 = n4.y();
  G4double z4 = n4.z();

  G4double Q1 = -(x1-x2)*(y3*z4-y4*z3);
  G4double Q2 =  (x1-x3)*(y2*z4-y4*z2);
  G4double R1 = -(x1-x4)*(y2*z3-y3*z2);
  G4double R2 = -(x2-x3)*(y1*z4-y4*z1);
  G4double S1 =  (x2-x4)*(y1*z3-y3*z1);
  G4double S2 = -(x3-x4)*(y1*z2-y2*z1);

  return (Q1+Q2+R1+R2+S1+S2)/6.0;
}

//calculate the area and assign into fSurfaceArea
void  G4Hex::calArea()
{
    //fetch all face, and trianglize them
    G4double theArea = 0.0;
    for (unsigned int i=0; i<fFaceNodesList.size(); i++)
        theArea += calAreaOfFace(i);
    //return
    fSurfaceArea = theArea;
}

//calculaet the area of the face with id= aFaceIdx
G4double G4Hex::calAreaOfFace(G4int aFaceIdx)
{
    G4double theArea = 0.0;
    vector <G4int> aFaceNodes = fFaceNodesList[aFaceIdx];
    G4ThreeVector aPoint1 = fNodeList[aFaceNodes[0]]; //the first point
    for (unsigned int i=0; i<aFaceNodes.size()-2; i++) {
        G4ThreeVector aPoint2 = fNodeList[aFaceNodes[i+1]]; // the second point
        G4ThreeVector aPoint3 = fNodeList[aFaceNodes[i+2]]; // the third point
        theArea += calTriangleArea(aPoint1, aPoint2, aPoint3);
    }
    return theArea;
}

//calculate the area of a triangle
//see G4Tet or vtkCellIntegrator::IntegrateTriangle
G4double  G4Hex::calTriangleArea(const G4ThreeVector & n1,
                                const G4ThreeVector & n2,
                                const G4ThreeVector & n3) const
{
    G4ThreeVector aVec1 = n2 - n1;
    G4ThreeVector aVec2 = n3 - n1;
    return aVec1.cross(aVec2).mag() / 2;
}

//calculte the barycenter (cell centroid)
void  G4Hex::calCenter()
{
    G4ThreeVector aCenter = G4ThreeVector(); //default 0,0,0
    for (unsigned int i=0; i<fNodeList.size(); i++) {
        aCenter += fNodeList[i];
    }
    aCenter /= G4double(fNodeList.size());
    fCentroid = aCenter;
}
//calculte the bounding sphere radius, a bounding
//sphere is a sphere inclose all vertices
void   G4Hex:: calBoundingSphereRadius()
{
    G4double aRadius = DBL_MAX;
    for (unsigned int i=0; i<fNodeList.size(); i++) {
        G4double tmpDouble = (fNodeList[i] - fCentroid).mag();
        if (tmpDouble < aRadius)
            aRadius = tmpDouble;
    }
    fMaxSize = aRadius;
}




//////////////////////////////////////////////////////////////////////////
//
// checking the crease face
// return false if Crease face found
bool G4Hex::checkCreaseFace()
{
    for (unsigned int i=0; i< fFaceNodesList.size(); i++)
    {
//        if (fFaceNodesList[i].size() < 3) {
//            G4Exception("G4Hex::checkCreaseFace()", "GeomSolids0002", FatalException,
//                        "Face is formed at least by three points.");
//            return false;
//        }
        //check other nodes begin at the forth node
        vector <G4int> aFaceNodes = fFaceNodesList[i];
        for (unsigned int j=3; j<aFaceNodes.size(); j++) {
            G4ThreeVector aNode = fNodeList[aFaceNodes[j]];  //get the other nodes
            G4double aDistance = calDistanceToFace(aNode, i);
            if (aDistance > fTol) return false; //if distance larger than tolerance, failed
        }
    }
    return true;
}

//reverse the orientation of the nodes, try to make volume to positive
void G4Hex::revOrientation()
{
    //swap node 1 <-> 3
    G4ThreeVector aBuffer = fNodeList[1];
    fNodeList[1] = fNodeList[3];
    fNodeList[3] = aBuffer;
    //swap node 5 <-> 7
    aBuffer = fNodeList[5];
    fNodeList[5] = fNodeList[7];
    fNodeList[7] = aBuffer;
}


//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Hex::G4Hex( __void__& a )
  : G4VSolid(a), fCentroid(0,0,0), fTol(0.), fMaxSize(0.),
    fCubicVolume(0.), fSurfaceArea(0.), fpPolyhedron(0), FoundNagtvVol(false)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Hex::~G4Hex()
{
  delete fpPolyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Hex::G4Hex(const G4Hex& rhs)
  : G4VSolid(rhs),
    fNodeList(rhs.fNodeList), fFaceNodesList(rhs.fFaceNodesList),
    fFaceNormalList(rhs.fFaceNormalList), fCdotNList(rhs.fCdotNList),
    fBoundaryBox(rhs.fBoundaryBox), fCentroid(rhs.fCentroid),
    fTol(rhs.fTol),fMaxSize(rhs.fMaxSize),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fpPolyhedron(0), FoundNagtvVol(false)
{
}


///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Hex& G4Hex::operator = (const G4Hex& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fNodeList =          rhs.fNodeList;
   fFaceNodesList =     rhs.fFaceNodesList;
   fFaceNormalList=     rhs.fFaceNormalList;
   fCdotNList=          rhs.fCdotNList;
   fBoundaryBox =       rhs.fBoundaryBox;
   fCentroid =          rhs.fCentroid;
   fTol =               rhs.fTol;
   fMaxSize =           rhs.fMaxSize;
   fCubicVolume =       rhs.fCubicVolume;
   fSurfaceArea =       rhs.fSurfaceArea;
   fpPolyhedron =       0;
   FoundNagtvVol=       false;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// CheckDegeneracy

//G4bool G4Hex::CheckDegeneracy( G4ThreeVector anchor,
//                               G4ThreeVector p2,
//                               G4ThreeVector p3,
//                               G4ThreeVector p4 )
//{
//  G4bool result;
//  G4Hex *object=new G4Hex("temp",anchor,p2,p3,p4,&result);
//  delete object;
//  return result;
//}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Hex::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Hex::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                    G4double& pMin, G4double& pMax) const
{
  G4double xMin,xMax;
  G4double yMin,yMax;
  G4double zMin,zMax;

  if (pTransform.IsRotated())
  {
      //get the transformed node list
    G4ThreeVectorList aNodeList ;
    for (unsigned int i=0; i<fNodeList.size(); i++)
        aNodeList.push_back(pTransform.TransformPoint(fNodeList[i]));

    //calculate the new boundary box
    BdBox aBdBox;
    aBdBox.calBdBox(aNodeList);
    xMin    = aBdBox.XMin;
    xMax    = aBdBox.XMax;
    yMin    = aBdBox.YMin;
    yMax    = aBdBox.YMax;
    zMin    = aBdBox.ZMin;
    zMax    = aBdBox.ZMax;
  }
  else
  {
    G4double xoffset = pTransform.NetTranslation().x() ;
    xMin    = xoffset + fBoundaryBox.XMin;
    xMax    = xoffset + fBoundaryBox.XMax;
    G4double yoffset = pTransform.NetTranslation().y() ;
    yMin    = yoffset + fBoundaryBox.YMin;
    yMax    = yoffset + fBoundaryBox.YMax;
    G4double zoffset = pTransform.NetTranslation().z() ;
    zMin    = zoffset + fBoundaryBox.ZMin;
    zMax    = zoffset + fBoundaryBox.ZMax;
  }

  if (pVoxelLimit.IsXLimited())
  {
    if ( (xMin > pVoxelLimit.GetMaxXExtent()+fTol) || 
         (xMax < pVoxelLimit.GetMinXExtent()-fTol)  )  { return false; }
    else
    {
      xMin = std::max(xMin, pVoxelLimit.GetMinXExtent());
      xMax = std::min(xMax, pVoxelLimit.GetMaxXExtent());
    }
  }

  if (pVoxelLimit.IsYLimited())
  {
    if ( (yMin > pVoxelLimit.GetMaxYExtent()+fTol) ||
         (yMax < pVoxelLimit.GetMinYExtent()-fTol)  )  { return false; }
    else
    {
      yMin = std::max(yMin, pVoxelLimit.GetMinYExtent());
      yMax = std::min(yMax, pVoxelLimit.GetMaxYExtent());
    }
    }

    if (pVoxelLimit.IsZLimited())
    {
      if ( (zMin > pVoxelLimit.GetMaxZExtent()+fTol) ||
           (zMax < pVoxelLimit.GetMinZExtent()-fTol)  )  { return false; }
    else
    {
      zMin = std::max(zMin, pVoxelLimit.GetMinZExtent());
      zMax = std::min(zMax, pVoxelLimit.GetMaxZExtent());
    }
  }

  switch (pAxis)
  {
    case kXAxis:
      pMin=xMin;
      pMax=xMax;
      break;
    case kYAxis:
      pMin=yMin;
      pMax=yMax;
      break;
    case kZAxis:
      pMin=zMin;
      pMax=zMax;
      break;
    default:
      break;
  }

  return true;
} 

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Hex::Inside(const G4ThreeVector& p) const
{
    //calculate the distance of this point to each face
//    G4bool isOutside = false;  //default: not outside
    G4bool isInside = true;    //default: is inside
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {  //each face
        G4double aDistance = calDistanceToFace(p, i);
//        isOutside = isOutside || (aDistance > fTol);
        if (aDistance > fTol) return kOutside;//check if outside one face
        isInside = isInside && (aDistance < -fTol);   //check if inside all face
    }

    if (isInside) return kInside;
    else return kSurface;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned.
// This assumes that we are looking from the inside!

G4ThreeVector G4Hex::SurfaceNormal( const G4ThreeVector& p) const
{
    //calculate the distance to each face
    vector <G4double> aDistanceList ;
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {
        aDistanceList.push_back(std::fabs(calDistanceToFace(p, i)));
    }
    const G4double delta = 0.5*kCarTolerance;
    G4ThreeVector sumnorm(0., 0., 0.);
    G4int noSurfaces=0;
    //see if the point is on the face
    for (unsigned int i=0; i<aDistanceList.size(); i++) {
        if (aDistanceList[i] <= delta) {
            noSurfaces ++; //on how many face: points might at the corner
            sumnorm += fFaceNormalList[i]; //sum the normal
        }
    }
    //check on face, and reutrn the normal
    if( noSurfaces > 0 )  {
        if( noSurfaces == 1 )
            return sumnorm;
        else
            return sumnorm.unit();
    }
    else { //not on the face, calculate the nearst face
        G4double tmpMin = DBL_MAX;
        G4int aFaceIdx = 0; //face index
        for (unsigned int i=0; i<aDistanceList.size(); i++) {
            if (aDistanceList[i] <tmpMin) {
                tmpMin = aDistanceList[i];
                aFaceIdx = i;
            }
        }
        return fFaceNormalList[aFaceIdx];
    }
}
///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
// All this is very unrolled, for speed.

G4double G4Hex::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{
    G4ThreeVector vu(v.unit()), hp;
    G4double vdotn, t, tmin=kInfinity;

    G4double extraDistance=10.0*fTol; // a little ways into the solid

    //loop each face
    for (unsigned int i=0; i<fFaceNodesList.size(); i++)
    {
        vdotn = -vu.dot(fFaceNormalList[i]); //-(ray vector dot face normal) = -cos(theta); here cos(theta) <0
        if (vdotn > 1e-12) // this is a candidate face, since it is pointing at us
        {
            t = (calDistanceToFace(p, i) / vdotn);// #  distance to intersection= calDistanceToFace * (-cos(theta))
            if ((t >= -fTol) && (t < tmin)) {
                hp=p+vu*(t+extraDistance); // a little beyond point of intersection, for detecting if inside
                G4bool isInside = true;
                for (unsigned int j=0; j<fFaceNodesList.size(); j++) {
                    if (j == i) continue; //skip current face
                    isInside = isInside && (calDistanceToFace(hp, j) < 0.0);   //check if exactly inside all face
                }
                if (isInside) tmin = t;  //update the nearst intersect distance
            }
        }
    }

    return std::max(0.0,tmin); //no less than 0.0
}

//////////////////////////////////////////////////////////////////////////
// 
// Approximate distance to tet.
// returns distance to sphere centered
// - If inside return 0

G4double G4Hex::DistanceToIn(const G4ThreeVector& p) const
{
  G4double dd=(p-fCentroid).mag() - fMaxSize - fTol;
  return std::max(0.0, dd);
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.

G4double G4Hex::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm, G4ThreeVector *n) const
{
    G4ThreeVector vu(v.unit());
    G4double  vdotn, TMin = kInfinity;
    vector <G4double > aTList;
    aTList.assign(fFaceNodesList.size(), kInfinity);

    //calculate the distance to each face
    for (unsigned int i=0; i<fFaceNodesList.size(); i++)     {
        vdotn=vu.dot(fFaceNormalList[i]);
        if (vdotn > 1e-12)  // #we're heading towards this face, so it is a candidate
            //be in mind the distance of a point is negative when inside the solid
            //therefore we should correct it with "-"
            aTList[i] = -calDistanceToFace(p, i) / vdotn;
    }

    G4int aFaceIdx =0;
    for (unsigned int i=0; i<aTList.size(); i++) {
        if (aTList[i]< TMin) {
            TMin = aTList[i];
            aFaceIdx = i;
        }
    }


    if ((TMin == kInfinity || TMin < -fTol))
    {
      DumpInfo();
      std::ostringstream message;
      message << "No good intersection found or already outside!?" << G4endl
              << "p = " << p / mm << "mm" << G4endl
              << "v = " << v  << G4endl;
      G4Exception("G4Hex::DistanceToOut(p,v,...)", "GeomSolids1002",
                  JustWarning, message);
      if(validNorm)
        *validNorm=false; // flag normal as meaningless

    }
    else if(calcNorm && n)
    {
        *n=fFaceNormalList[aFaceIdx];
        if(validNorm) { *validNorm=true; }
    }

    return std::max(TMin,0.0); // avoid TMin<0.0 by a tiny bit
                             // if we are right on a face
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0

G4double G4Hex::DistanceToOut(const G4ThreeVector& p) const
{
    G4double aDistMin = kInfinity;
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {
        //be in mind the distance of a point is negative when inside the solid
        //therefore we should correct it with "-"
        G4double aDistance = -calDistanceToFace(p, i);
        if (aDistance < aDistMin) aDistMin = aDistance;
    }

  // if any one of these is negative, we are outside,
  // so return zero in that case
    return (aDistMin < fTol)? 0:aDistMin;
}

////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Note: Caller has deletion responsibility

G4ThreeVectorList*
G4Hex::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
  G4ThreeVectorList* vertices = new G4ThreeVectorList();

  if (vertices)  {
      for (unsigned int i=0; i<fNodeList.size(); i++)
          vertices->push_back(pTransform.TransformPoint(fNodeList[i]));
  }
  else  {
    DumpInfo();
    G4Exception("G4Hex::CreateRotatedVertices()",
                "GeomSolids0003", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }
  return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Hex::GetEntityType() const
{
  return G4String("G4Hex");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Hex::Clone() const
{
  return new G4Hex(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Hex::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
  << "    *** Dump for solid - " << GetName() << " ***\n"
  << "    ===================================================\n"
  << " Solid type: G4Hex\n"
  << " Parameters: \n" ;
  for (unsigned int i=0; i<fNodeList.size(); i++) {
      os << "p"<<i<<"\t: "<<fNodeList[i]/mm<<" mm \n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


////////////////////////////////////////////////////////////////////////
//
// GetPointOnFace
//
// Auxiliary method for get point on surface
// aFaceArea: return area of this face to avoid repeat calculation

G4ThreeVector G4Hex::GetPointOnFace(G4int aFaceIdx, G4double & aFaceArea) const
{
    //trianglize the face, and get a random point in each triangle
    vector <G4ThreeVector > aRandPointList;
    vector <G4double> aAreaList;

    vector <G4int> aFaceNodes = fFaceNodesList[aFaceIdx];
    G4ThreeVector aPoint1 = fNodeList[aFaceNodes[0]]; //the first point
    for (unsigned int i=0; i<aFaceNodes.size()-2; i++) {
        G4ThreeVector aPoint2 = fNodeList[aFaceNodes[i+1]]; // the second point
        G4ThreeVector aPoint3 = fNodeList[aFaceNodes[i+2]]; // the third point
        aRandPointList.push_back(GetPointOnTriangle(aPoint1,aPoint2, aPoint3));
        aAreaList.push_back( calTriangleArea( aPoint1, aPoint2, aPoint3));
    }

    //use area as probibility function, get a point among these random points
    vector <G4double> aAreaProbTable;
    aAreaProbTable.push_back(0.0); //first should be 0.
    for (unsigned int i=0; i<aAreaList.size(); i++)
        //aAreaProbTable contains: 0.0. A1, A1+A2, A1+A2+A3, ...
        aAreaProbTable.push_back(aAreaProbTable[i] + aAreaList[i]);
    aAreaList.clear(); //no use any more

    aFaceArea = aAreaProbTable[aAreaProbTable.size() -1]; //return the total area of this face, avoid recalculation
    G4double aNeedle = RandFlat::shoot(0.,aFaceArea); //0~ total area of this face
    for (unsigned int i=1; i<aAreaProbTable.size(); i++) {  //start from 1
        if (aNeedle >= aAreaProbTable[i-1] && aNeedle <aAreaProbTable[i])
            return aRandPointList[i-1];  //return i-1, not i
    }
    return aRandPointList[aAreaProbTable.size() -2];  //in case aNeedle = aFaceArea
}


G4ThreeVector G4Hex::GetPointOnTriangle(G4ThreeVector p1, G4ThreeVector p2,
                                    G4ThreeVector p3) const
{
  G4double lambda1,lambda2;
  G4ThreeVector v, w;
//  v = p3 - p1;
//  w = p1 - p2;
  v = p3 - p1; //we consider p1 as base
  w = p2 - p1;

  lambda1 = RandFlat::shoot(0.,1.);
//  lambda2 = RandFlat::shoot(0.,lambda1);//this produce non-uniform random points
  lambda2 = RandFlat::shoot(0.,1.);
  if ((lambda1 + lambda2) > 1) {
      lambda1 = 1 - lambda1;
      lambda2 = 1 - lambda2;
  }
//  return (p2 + lambda1*w + lambda2*v);
  return (p1 + lambda1*w + lambda2*v); //p1 as base
}


////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Hex::GetPointOnSurface() const
{
    //get a random point in each face
    vector <G4ThreeVector > aRandPointList;
    vector <G4double> aAreaList;
    for (unsigned int i=0; i<fFaceNodesList.size(); i++) {
        G4double aArea = 0.0; //for getting area
        aRandPointList.push_back(GetPointOnFace(i, aArea));
        aAreaList.push_back(aArea);
    }

    //use area as probibility function, get a point among these random points
    vector <G4double> aAreaProbTable;
    aAreaProbTable.push_back(0.0); //first should be 0.
    for (unsigned int i=0; i<aAreaList.size(); i++)
        //aAreaProbTable contains: 0.0. A1, A1+A2, A1+A2+A3, ...
        aAreaProbTable.push_back(aAreaProbTable[i] + aAreaList[i]);
    aAreaList.clear(); //no use any more

    G4double aTotalArea = aAreaProbTable[aAreaProbTable.size() -1]; // the total area of this solid, avoid recalculation
    G4double aNeedle = RandFlat::shoot(0.,aTotalArea); //0~ total area of this face
    for (unsigned int i=1; i<aAreaProbTable.size(); i++) {  //start from 1
        if (aNeedle >= aAreaProbTable[i-1] && aNeedle <aAreaProbTable[i])
            return aRandPointList[i-1];  //return i-1, not i
    }
    return aRandPointList[aAreaProbTable.size() -2];  //in case aNeedle = aTotalArea

}

////////////////////////////////////////////////////////////////////////
//
// GetVertices

std::vector<G4ThreeVector> G4Hex::GetVertices() const
{
    //we know they are the same type
    return static_cast <std::vector<G4ThreeVector> > (fNodeList);
}

////////////////////////////////////////////////////////////////////////
//
// GetCubicVolume

G4double G4Hex::GetCubicVolume()
{
  return fCubicVolume;
}

////////////////////////////////////////////////////////////////////////
//
// GetSurfaceArea

G4double G4Hex::GetSurfaceArea()
{
  return fSurfaceArea;
}

// Methods for visualisation

////////////////////////////////////////////////////////////////////////
//
// DescribeYourselfTo

void G4Hex::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////////////////////////
//
// GetExtent

G4VisExtent G4Hex::GetExtent() const
{
  return G4VisExtent (fBoundaryBox.XMin, fBoundaryBox.XMax,
                      fBoundaryBox.YMin, fBoundaryBox.YMax,
                      fBoundaryBox.ZMin, fBoundaryBox.ZMax);
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron

G4Polyhedron* G4Hex::CreatePolyhedron () const
{
  G4Polyhedron *ph=new G4Polyhedron;
  typedef G4int G4int4[4];
  typedef G4double G4double3[3];

  G4double3 * aPointArray = new G4double3 [fNodeList.size()];
  //put the nodes into the point list
  for (unsigned int i=0; i<fNodeList.size(); i++) {
      aPointArray[i][0] = fNodeList[i].x();
      aPointArray[i][1] = fNodeList[i].y();
      aPointArray[i][2] = fNodeList[i].z();
  }

  //we fist trianglize all faces
  vector < vector <G4int> > aFacetList;

  for (unsigned int i=0; i<fFaceNodesList.size(); i++)
  {
      vector <G4int> aFaceNodes = fFaceNodesList[i];
      vector <G4int> aFacet;
      aFacet.resize(4);
      aFacet[0] = aFaceNodes[0] + 1; //the first point, +1 because required by G4Polyhedron
      aFacet[3] = 0; // 0 indicates this is a triangle facet
      for (unsigned int j=0; j<aFaceNodes.size()-2; j++) {
          aFacet[1] = aFaceNodes[j+1] +1; // the second point, +1 because required by G4Polyhedron
          aFacet[2] = aFaceNodes[j+2] +1; // the third point, +1 because required by G4Polyhedron
          aFacetList.push_back(aFacet);
      }
  }

  G4int4 * aFacetArray = new G4int4 [aFacetList.size()];
  for (unsigned int i=0; i<aFacetList.size(); i++) {
      aFacetArray[i][0] = aFacetList[i][0];
      aFacetArray[i][1] = aFacetList[i][1];
      aFacetArray[i][2] = aFacetList[i][2];
      aFacetArray[i][3] = aFacetList[i][3];
  }

  ph->createPolyhedron(fNodeList.size(), aFacetList.size(),
                       aPointArray, aFacetArray);
  delete [] aPointArray;
  delete [] aFacetArray;

  return ph;
}

////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron

G4Polyhedron* G4Hex::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}


void BdBox::calBdBox(G4ThreeVectorList & aNodeList)
{
     XMin = DBL_MAX,  YMin = DBL_MAX, ZMin = DBL_MAX;
     XMax = DBL_MIN,  YMax = DBL_MIN, ZMax = DBL_MIN;
    for (unsigned int i=0; i< aNodeList.size(); i++) {
        G4ThreeVector aPoint =  aNodeList[i];
        if (aPoint.x() < XMin) XMin = aPoint.x();
        if (aPoint.y() < YMin) YMin = aPoint.y();
        if (aPoint.z() < ZMin) ZMin = aPoint.z();
        if (aPoint.x() > XMax) XMax = aPoint.x();
        if (aPoint.y() > YMax) YMax = aPoint.y();
        if (aPoint.z() > ZMax) ZMax = aPoint.z();
    }
}

