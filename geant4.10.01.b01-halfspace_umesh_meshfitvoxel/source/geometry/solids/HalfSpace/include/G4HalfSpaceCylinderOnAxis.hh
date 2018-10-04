#ifndef G4HALFSPACECYLINDERONAXIS_HH
#define G4HALFSPACECYLINDERONAXIS_HH

#include "G4HalfSpaceCylinder.hh"

class G4HalfSpaceCylinderOnAxis : public G4HalfSpaceCylinder
{
public:
    G4HalfSpaceCylinderOnAxis(G4HSAxis thePllAxis, const G4double & Radius,
                        const G4int Sense = 1);
    //cylinder which axis is on coordinate axis
    //thePllAixs: which axis the cylinder is on
    //Radius:
    //Sense: which half-space the surface represent

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint) const ;
    //return the normal of the surface at the location of a point
    //except plane, all other surface should use a point as input
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const ;
    //calculate the distance to the surface
    //if the point is "inside" the surface, return a positive value
    //Should be implement by subclass

    virtual EInside         Inside(const G4ThreeVector & aPoint) const ;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //we substitute the point into the quadric equation,
    //if the result >0 : Inside the half-space
    //if result<0: Outside the half-space
    //if result=0: on the surface

    virtual   G4HSSurfType  getType() const {return G4HSCylinder;}
    //return surface type

    G4HalfSpaceCylinderOnAxis(const G4HalfSpaceCylinderOnAxis & right);
    //copy constructor
    G4HalfSpaceCylinderOnAxis& operator=(const G4HalfSpaceCylinderOnAxis & right);
    //copy operator
    G4HalfSpaceCylinderOnAxis operator-();
    //reverse Sense operator

protected:
    //Aq is the same with G4HalfSpaceCylinder

    virtual G4double        Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector);
    //method to calcaulate the Bq,
    //Bq = 2* (A*X0*Xd + B*(X0*Yd+ Xd*Y0) + C*(X0*Zd+ Xd*Z0) + D*Xd +
    //         E*Y0*Yd + F*(Y0*Zd + Yd*Z0) + G*Yd + H*Z0*Zd + I*Zd)
    virtual G4double        Cq(const G4ThreeVector & aPoint);
    //Cq= A*X0^2 + 2*B*X0*Y0 + 2*C*X0*Z0 + 2*D*X0
    //  + E*Y0^2 +  2*F*Y0*Z0 + 2*G*Y0
    //    + H*Z0^2 + 2*I*Z0 + J

};

#endif // G4HALFSPACECYLINDERONAXIS_HH
