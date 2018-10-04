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
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4UMeshElement1st
//
// Class description:
//
//   A G4UMeshElement1st is a base class for first-order unstructured mesh
//   elements.
//

// History:
// -------
// 28.03.2014 - Y.Qiu(Karlsruhe Institute of Technology)
// --------------------------------------------------------------------
#ifndef G4UMESHELEMENT1ST_HH
#define G4UMESHELEMENT1ST_HH

#include "G4VSolid.hh"

using namespace std;

struct BdBox {
    G4double    XMin;
    G4double    XMax;
    G4double    YMin;
    G4double    YMax;
    G4double    ZMin;
    G4double    ZMax;
    G4double    getDx() {return (XMax - XMin)/2;} //half length of the box
    G4double    getDy() {return (YMax - YMin)/2;}
    G4double    getDz() {return (ZMax - ZMin)/2;}
    void        calBdBox(G4ThreeVectorList & aNodeList);
    //calculate the boundary box with a list of nodes
};


class G4UMeshElement1st : public G4VSolid
{

  public:  // with description

    G4UMeshElement1st(const G4String& pName );

     ~G4UMeshElement1st();

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;
    // Methods for solid

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;

    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=false,
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;

    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    virtual G4GeometryType GetEntityType() const;

    virtual G4VSolid* Clone() const;

    virtual std::ostream& StreamInfo(std::ostream& os) const;

    G4ThreeVector GetPointOnSurface() const;

    // Functions for visualization

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4Polyhedron* GetPolyhedron      () const;

  public:   // without description

    G4UMeshElement1st(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UMeshElement1st(const G4UMeshElement1st& rhs);
    G4UMeshElement1st& operator=(const G4UMeshElement1st& rhs);
      // Copy constructor and assignment operator.

    const char* CVSHeaderVers()
      { return "$Id: G4UMeshElement1st.hh 76263 2013-11-08 11:41:52Z gcosmo $"; }
    const char* CVSFileVers()
      { return CVSVers; }
//    void PrintWarnings(G4bool flag)
//      { warningFlag=flag; }
//    static G4bool CheckDegeneracy(G4ThreeVector anchor,
//                                  G4ThreeVector p2,
//                                  G4ThreeVector p3,
//                                  G4ThreeVector p4);
    std::vector<G4ThreeVector> GetVertices() const;
      // Return the four vertices of the shape.

  protected:  // with description

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

  private:

    static const char CVSVers[];

protected:

    G4ThreeVector GetPointOnFace(G4int aFaceIdx, G4double &aFaceArea) const;
    G4ThreeVector GetPointOnTriangle(G4ThreeVector p1, G4ThreeVector p2,
                                 G4ThreeVector p3) const;
    void            preProc();
    //pre-processing the solid, should be called in derived class constructor

    // ATTENTION! the following four method is PURE virtual,
    // SHOULD BE implenment in derive class
    virtual bool            checkPosVolume() {return true; /*Not implemented yet*/};
    //check if the volume is positive. Not imiplemented
    virtual void            formFaces() {};
    //form the faces using id of nodes in fNodeList
    virtual void            calVolume()  {};
    //calculate the volume size and assign to fCubicVolume
    virtual void            revOrientation()  {};
    //reverse the orientation of the nodes, try to make volume to positive


    bool            checkCreaseFace();
    //check if any face is Crease
    void            calFaceNormal();
    //calculate the face normal
    void            calFaceNdotCenter();
    //calculate the face normal dot face center point

    G4double        calDistanceToFace(const G4ThreeVector &aPoint, G4int aFaceIdx) const;
    //calculate the distance of a point to the face indexed by the a id


    double          calTetraVolume(const G4ThreeVector & n1,
                                 const G4ThreeVector & n2,
                                 const G4ThreeVector & n3,
                                 const G4ThreeVector & n4);
    //calculate the volume size of a tetrahedron
    void            calArea();
    //calculate the area and assign into fSurfaceArea
    G4double calAreaOfFace(G4int aFaceIdx);
    //calculaet the area of the face with id= aFaceIdx
    G4double calTriangleArea(const G4ThreeVector & n1,
                                    const G4ThreeVector & n2,
                                    const G4ThreeVector & n3) const;
    //calculate the area of a triangle

    void            calCenter();
    //calculte the barycenter (centroid)
    void            calBoundingSphereRadius();
    //calculte the bounding sphere radius, a bounding
    //sphere is a sphere inclose all vertices






  protected:

    G4ThreeVectorList fNodeList;
    //vector container to contain nodes
    vector <vector <G4int > > fFaceNodesList;
    //2D container to contain node id for each face
    //fFaceNodesList.size() == number of faces
    G4ThreeVectorList fFaceNormalList;
    //to contain normal of each face
    vector <G4double> fCdotNList;
    //face center dot face normal, for convenience of calculation
    BdBox fBoundaryBox;
    //boundary box
    G4ThreeVector fCentroid;
    //centroid
    G4double  fTol, fMaxSize;
    //tolerance and Radius of the bounding sphere
    G4double fCubicVolume, fSurfaceArea;

    mutable G4Polyhedron* fpPolyhedron;

    G4bool  FoundNagtvVol;
    //for checking negative volume
    //if one of the splitted Tetrahedron has negative volume
    //then FoundNagtvVol is true

};

#endif

