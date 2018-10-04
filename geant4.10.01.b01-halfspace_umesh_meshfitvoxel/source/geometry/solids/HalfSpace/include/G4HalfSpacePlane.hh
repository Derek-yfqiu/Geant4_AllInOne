#ifndef G4HALFSPACEPLANE_HH
#define G4HALFSPACEPLANE_HH

#include "G4HalfSpaceSurface.hh"
#include "G4Plane.hh"

/*!
 * \brief The G4HalfSpacePlane class
 * half space plane,
 * \sa An Introduction to Ray Tracing
 */
class G4HalfSpacePlane : public G4HalfSpaceSurface
{

public:
//    G4HalfSpacePlane();
    G4HalfSpacePlane(const G4double aA, const G4double aB,
                     const G4double aC,const G4double aD,
                     const G4int    Sense = 1);
    // equation of the plane:
    // A*X+B*Y+C*Z - D=0
    G4HalfSpacePlane(G4HSAxis theAxis, const G4double aD,
                     const G4int    Sense = 1);
    //Planes which normal to one axis,
    //theAxis : normal to axis: AxisX,AxisY,AxisZ
    //aD Distance from the origin
    //Sense: on which half-space

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint = G4ThreeVector()) const ;
    //return the normal of the surface at the location of a point
    //except plane, all other surface should use a point as input
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                      std::vector <G4ThreeVector> & IntersectPoints,
                                      std::vector <G4double>  &     IntersectDistances) ;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances

    virtual EInside         Inside(const G4ThreeVector & aPoint) const ;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //ATTENTION: the position is decided by the solid, therefore the position should be
    //adjust by m_Sense

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const ;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value

    G4HSSurfType            getType() const {return G4HSPlane;}
    //return surface type

    G4String                printMe() const {
                                       std::stringstream strstr;
                                       strstr<< "Plane. D="<<m_Plane.d <<"\t A="<<m_Plane.a<<
                                                       "\t B="<<m_Plane.b<<"\t C="<<m_Plane.c<<std::endl;
                                       return strstr.str();}



    G4HalfSpacePlane(const G4HalfSpacePlane & right);
    //copy constructor
    G4HalfSpacePlane& operator=(const G4HalfSpacePlane & right);
    //copy operator
    G4HalfSpacePlane operator-();
    //reverse Sense operator

protected:
    G4Plane                 m_Plane;
    //store parameters for this plane
    G4ThreeVector           m_Normal;
    //plane have a unique normal everywhere


};

#endif // G4HALFSPACEPLANE_HH
