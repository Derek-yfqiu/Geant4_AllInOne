#ifndef G4HALFSPACESPHERE_HH
#define G4HALFSPACESPHERE_HH
#include "G4HalfSpaceSurface.hh"

/*!
 * \brief The G4HalfSpaceSphere class
 *  half-space sphere surface
 *  \sa An Introduction to Ray Tracing
 */
class G4HalfSpaceSphere : public G4HalfSpaceSurface
{
public:
    G4HalfSpaceSphere(const G4double & CenterX,
                      const G4double & CenterY,
                      const G4double & CenterZ,
                      const G4double & Radius,
                      const G4int    Sense = 1);
    //equation: (X-CenterX)^2 + (Y-CenterY)^2 + (Z-CenterZ)^2 = R^2

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint) const ;
    //return the normal of the surface at the location of a point
    //except plane, all other surface should use a point as input
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                      std::vector <G4ThreeVector> & IntersectPoints,
                                      std::vector <G4double>  &     IntersectDistances) ;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances
    //ATTENTION: Remenber to Reset before using this method!
    //algorithm : see An Introduction to Ray Tracing

    virtual EInside         Inside(const G4ThreeVector & aPoint) const ;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //ATTENTION: the position is decided by the solid, therefore the position should be
    //adjust by m_Sense

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const ;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value

    virtual G4HSSurfType    getType() const {return G4HSSphere;}
    //return surface type

    G4String                printMe() const {std::stringstream strstr;
                                             strstr << "Sphere. Radius="<<m_Radius<<"\t Center="<<m_Center.X
                                                    <<",/t"<<m_Center.Y<<",/t"<<m_Center.Z<<std::endl;
                                            return strstr.str();}


    G4HalfSpaceSphere(const G4HalfSpaceSphere & right);
    //copy constructor
    G4HalfSpaceSphere& operator=(const G4HalfSpaceSphere & right);
    //copy operator
    G4HalfSpaceSphere operator-();
    //reverse Sense operator

protected:

    G4ThreeVector           m_Center;
    //the sphere center
    G4double                m_Radius;
    //the radius
    G4double                m_Radius2;
    //square of the radius, for calculation convenience



};

#endif // G4HALFSPACESPHERE_HH
