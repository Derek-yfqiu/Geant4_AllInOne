#ifndef G4HALFSPACEELLIPTICALTORUS_HH
#define G4HALFSPACEELLIPTICALTORUS_HH

#include "G4HalfSpaceTorus.hh"

class G4HalfSpaceEllipticalTorus : G4HalfSpaceTorus
{
public:
    G4HalfSpaceTorus(const G4ThreeVector & Center, const G4ThreeVector & Axis,
                     const G4double & MaxRadius, const G4double & MinRadius1,
                     const G4double & MinRadius2);

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint) const ;
    //return the normal of the surface at the location of a point
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual EInside         Inside(const G4ThreeVector & aPoint) const ;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //ATTENTION: the position is decided by the solid, therefore the position should be
    //adjust by m_Sense

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value

    virtual  G4HSSurfType   getType() const {return G4HSEllipticalTorus;}
    //return surface type

    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec) ;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances
    //ATTENTION: Remenber to Reset before using this method!
    //algorithm : from G4BREPToroidalSurface


protected:
    G4double                m_r2;
    //the second small radius

};

#endif // G4HALFSPACEELLIPTICALTORUS_HH
