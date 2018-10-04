#ifndef G4HALFSPACESURFACE_HH
#define G4HALFSPACESURFACE_HH

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"


#include <vector>
//using namespace std;

enum  G4HSSurfType {
    G4HSPlane,
    G4HSSphere,
    G4HSQuadric,
    G4HSCylinder,
    G4HSCone,
    G4HSTorus
};
//for denoting axises the surface parallel or on
enum G4HSAxis {
    AxisX,
    AxisY,
    AxisZ
};

class G4HalfSpaceSurface
{
public:
    G4HalfSpaceSurface();

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint = G4ThreeVector()) const = 0;
    //return the normal of the surface at the location of a point
    //except plane, all other surface should use a point as input
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                      std::vector <G4ThreeVector> & IntersectPoints,
                                      std::vector <G4double>  &     IntersectDistances) = 0;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances

    virtual EInside         Inside(const G4ThreeVector & aPoint) const = 0;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //ATTENTION: the position is decided by the solid, therefore the position should be
    //adjust by m_Sense


    virtual G4double        HowNear( const G4ThreeVector& aPoint )const   ;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value

    virtual G4HSSurfType    getType() const = 0;
    //return surface type

    virtual G4String        printMe() const = 0;



    inline G4int            getSense() {return m_Sense;};
    //return the sense of the surface
    inline void             setSense(const G4int & aSense) {m_Sense = aSense;};
    //set the sense
    void                    Reset();
    //reset this surface
    inline void             setActive(const G4bool & isActive) {m_isActive = isActive;};
    inline G4bool           getIsActive() {return m_isActive;};
    //get and set the active status

    inline std::vector <G4ThreeVector>  getIntersectPoints() {return m_IntersectPoints ; };
    inline std::vector <G4double>  getIntersectDistances() {return m_IntersectDistances ; };

    inline void             setTolerance(const G4double & aTolerance)
                            {kCarTolerance = aTolerance;};

    G4HalfSpaceSurface(const G4HalfSpaceSurface & right);
    //copy constructor
    G4HalfSpaceSurface& operator=(const G4HalfSpaceSurface & right);
    //copy operator



protected:
    G4int                   m_Sense;
    // sense to indicate which half-space it represent:  1--positive side, -1-- negative side
//    G4String                m_Type;
//    //surface type: plane, cylinder, cone, ...
//    G4double m_Distance;
//    //use to get the intersect distance

    std::vector <G4ThreeVector>  m_IntersectPoints;
    //restore the intersect points after the Intersect calculation
    std::vector <G4double>       m_IntersectDistances;
    //restore the intersect distance after the Intersect calculation
    G4double kCarTolerance, kAngTolerance;
    //tolerance of coordinates and angle

    G4bool                  m_isActive;



};

#endif // G4HALFSPACESURFACE_HH
