#ifndef G4HALFSPACESOLID_HH
#define G4HALFSPACESOLID_HH

#include "G4VSolid.hh"
#include "G4HalfSpaceSurface.hh"
#include "G4BoundingBox3D.hh"


class G4HalfSpaceSolid : public G4VSolid
{
public:
    G4HalfSpaceSolid(const G4String& name);
    G4HalfSpaceSolid(const G4String& name, std::vector<G4HalfSpaceSurface*> aSurfaceList,
                     const G4ThreeVector& BBoxLowerPoint, const G4ThreeVector& BBoxHigherPoint,
                     const G4double & Volume = 0.0, const G4double & SurfaceArea = 0.0  );
    //constructor for the solid, including a list of surface, boundary box,
    //and if applicable, volume and surface area

    virtual ~G4HalfSpaceSolid();


//    virtual void Initialize();
//    // ???

    G4bool      checkValidity();
    //check if the solid is valid: logical bounded solid
    const G4BoundingBox3D * getBoundaryBox() const {return m_BBox;};
    void        setBoundaryBox(const G4ThreeVector& aLowerPoint, const G4ThreeVector& aHigherPoint  );
    void        setBoundaryBox( G4BoundingBox3D * aBBox){m_BBox = aBBox;};

    //set the boundary for the solid
    //ATTENTION: this is a mendatory function for building a solid!
    inline G4bool      isBoundaryBox () const {return m_BBox != NULL; };
    //return if boundary if set or no

    void        setVolume(const G4double & aVolume);
    void        setSurfaceArea(const G4double & aArea);
    //method to set volume and surface area from outside of this class

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis              pAxis      ,
                           const G4VoxelLimits&     pVoxelLimit,
                           const G4AffineTransform& pTransform ,
                           G4double&                pMin       ,
                           G4double&                pMax        ) const;
    // Calculate the minimum and maximum extent of the solid, when under the
    // specified transform, and within the specified limits. If the solid
    // is not intersected by the region, return false, else return true.

    virtual EInside Inside(const G4ThreeVector &aPoint) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector&aPoint) const;
    // Calculates the normal of the surface at a point on the surface
    // The sense of the normal depends on the sense of the surface.
    //Note: we assume the point is on the surface, therefore check for the surface
    //within tolerences. if outside of tolerence, choose the closest surface and calculate
    //the normal

    virtual G4double DistanceToIn(const G4ThreeVector&aPoint) const;
    // Calculates the shortest distance ("safety") from a point
    // outside the solid to any boundary of this solid.
    // Return 0 if the point is already inside.

    virtual G4double DistanceToIn(const G4ThreeVector &aPoint,
                                  const G4ThreeVector &aVector) const;
    // Calculates the distance from a point Pt outside the solid
    // to the solid's boundary along a specified direction vector V.
    // Note: Intersections with boundaries less than the tolerance must
    //       be ignored if the direction is away from the boundary.

    virtual G4double DistanceToOut(const G4ThreeVector&aPoint) const;
    // Calculates the shortest distance ("safety") from a point inside the
    // solid to any boundary of this solid.
    // Return 0 if the point is already outside.

    virtual G4double DistanceToOut(register const G4ThreeVector& aPoint,
                                   register const G4ThreeVector& aVector,
                                   const G4bool  calcNorm=false ,
                                   G4bool        *validNorm=0   ,
                                   G4ThreeVector *n=0             ) const;
    // Calculates the distance from a point inside the solid to the solid`s
    // boundary along a specified direction vector.
    // Return 0 if the point is already outside.
    // Note: If the shortest distance to a boundary is less than the
    //       tolerance, it is ignored. This allows for a point within a
    //       tolerant boundary to leave immediately.


    virtual G4String GetEntityType() const;
    // Returns identifier for solid type entity.

    virtual G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.
    // The caller has responsibility for ownership.

    virtual std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

    void DescribeYourselfTo (G4VGraphicsScene& scene) const;
    // Dispatch function which identifies the solid to the graphics scene.

    G4Polyhedron* CreatePolyhedron () const;
    // Create a G4Polyhedron/...  (It is the caller's responsibility
    // to delete it).  A null pointer means "not created".
    virtual G4Polyhedron* GetPolyhedron () const;
    // Smart access function - creates on request and stores for future
    // access.  A null pointer means "not available".

    inline void    SetPolyhedron(G4Polyhedron* aPolyhedron) {fpPolyhedron = aPolyhedron;};
    //set the polyhedron representing this solid;
    //ATTENTION:because it is difficult to generate the polyhedron by this code,
    //for displaying this solid the polyhedron should be set
    //this solid takes the ownership of this polyhedron

    virtual G4double GetCubicVolume();
      // return the volume
    // we check if the volume is set or not,
    //if not, calculaet and cache it

    virtual G4double GetSurfaceArea();
    // return the surface area
  // we check if the area is set or not,
  //if not, calculaet and cache it

    virtual G4ThreeVector GetPointOnSurface() const;
      // Returns a random point located on the surface of the solid.
      // Points returned are not necessarily uniformly distributed.

public:  // without description

  G4HalfSpaceSolid(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4HalfSpaceSolid(const G4HalfSpaceSolid& rhs);
  G4HalfSpaceSolid& operator=(const G4HalfSpaceSolid& rhs);
    // Copy constructor and assignment operator.

protected:

  G4ThreeVectorList* CreateRotatedVertices(const G4AffineTransform&) const;
  void               Reset() const ;
  //reset this solid, activate all the surfaces
  void               setCarTolerance(const G4double & aCarTolerance) const;
  G4bool             Outside(const G4ThreeVector &aPoint, const G4int & skipIdx = -1) const;
  //return true if the point outside one of the surface
  //skipIdx >= 0 to specify a surface to skip


private:
    std::vector <G4HalfSpaceSurface*> m_SurfaceList;
    //a list of surface which form this convex solid
    G4BoundingBox3D  *           m_BBox;
    //bounadry box of the solid
//    G4bool                       m_isActive;
//    //the solid is active or not


    G4int    fStatistics;
    G4double fCubVolEpsilon;
    G4double fAreaAccuracy;
    G4double fCubicVolume; //solid volume, if = 0 then it has to be calcualted
    G4double fSurfaceArea; //surface area, if = 0 then it has to be calcualted
    // Statistics, error accuracy and cached value for volume and area.

    mutable G4Polyhedron* fpPolyhedron;




};

#endif // G4HALFSPACESOLID_HH
