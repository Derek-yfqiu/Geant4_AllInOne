#ifndef G4HALFSPACECONE_HH
#define G4HALFSPACECONE_HH
#include "G4HalfSpaceQuadric.hh"


class G4HalfSpaceCone  : public G4HalfSpaceQuadric
{
public:
//    G4HalfSpaceCone(G4HSAxis thePllAxis, const G4double & T2, const G4double &XYZpi,
//                    const G4int Sense = 1);
    //cone which axis is on the coordinate
    //thePllAixs: which axis the cylinder is on
    //T2: see equation, represent the slope of the line on cone surface
    //XYZpi: see equation,
    //Sense: which half-space the surface represent
    G4HalfSpaceCone();
    G4HalfSpaceCone( G4HSAxis thePllAxis, const G4double & T2,
                     const G4double & Xpi, const G4double & Ypi,
                    const G4double & Zpi, const G4int Sense = 1);
    //cone which axis is parallel with the coordinate
    //thePllAixs: which axis the cylinder is parallel
    //T2: see equation, represent the slope of the line on cone surface
    //Xpi, Ypi, Zpi; see equation,
    //Sense: which half-space the surface represent

    virtual G4ThreeVector   Normal(G4ThreeVector aPoint) const ;
    //return the normal of the surface at the location of a point
    //except plane, all other surface should use a point as input
    //ATTENTION: the normal should point outward of the solid, therefore should
    //be adjust by m_Sense

    virtual EInside         Inside(const G4ThreeVector & aPoint) const ;
    //judge the position of the point if inside/outside/on the surface
    //return the position kInside/kOutside/kSurface
    //we substitute the point into the quadric equation,
    //if the result >0 : Inside the half-space
    //if result<0: Outside the half-space
    //if result=0: on the surface

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const ;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value
    //Should be implement by subclass


    virtual   G4HSSurfType  getType() const {return G4HSCone;}
    //return surface type

    G4String                printMe() const {std::stringstream strstr;
                                             strstr << "Cone. Axis="<< char(88/*decimal value of X*/+int(m_PllAxis))
                                                    <<"\t T2="<<m_T2<<"\t Xpi="<<m_Xpi<<"\t Ypi="<<m_Ypi<<"\t Zpi="<<m_Zpi<<std::endl;
                                              return strstr.str();}

    G4HalfSpaceCone(const G4HalfSpaceCone & right);
    //copy constructor
    G4HalfSpaceCone& operator=(const G4HalfSpaceCone & right);
    //copy operator
    G4HalfSpaceCone operator-();
    //reverse Sense operator

protected:
    virtual G4double        Aq(const G4ThreeVector & aVector);
    //method to calcuate the Aq. the aVector is the vector of the ray
    //Aq = A*X^2 +  2*B*X*Y + 2*C*X*Z + E*Y^2 + 2*F*Y*Z  + H*Z^2

    virtual G4double        Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector);
    //method to calcaulate the Bq,
    //Bq = 2* (A*X0*Xd + B*(X0*Yd+ Xd*Y0) + C*(X0*Zd+ Xd*Z0) + D*Xd +
    //         E*Y0*Yd + F*(Y0*Zd + Yd*Z0) + G*Yd + H*Z0*Zd + I*Zd)
    virtual G4double        Cq(const G4ThreeVector & aPoint);
    //Cq= A*X0^2 + 2*B*X0*Y0 + 2*C*X0*Z0 + 2*D*X0
    //  + E*Y0^2 +  2*F*Y0*Z0 + 2*G*Y0
    //    + H*Z0^2 + 2*I*Z0 + J
protected:
    G4HSAxis m_PllAxis; //the axis which the cylinder paralleled
    G4double m_T2;
    G4double m_Xpi, m_Ypi, m_Zpi;

};

#endif // G4HALFSPACECONE_HH
