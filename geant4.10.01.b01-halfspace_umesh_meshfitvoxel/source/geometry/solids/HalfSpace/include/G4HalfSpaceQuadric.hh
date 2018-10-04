#ifndef G4HALFSPACEQUADRIC_HH
#define G4HALFSPACEQUADRIC_HH
#include "G4HalfSpaceSurface.hh"


class G4HalfSpaceQuadric : public G4HalfSpaceSurface
{
public:
    G4HalfSpaceQuadric();
    G4HalfSpaceQuadric (const G4double & aA,
                        const G4double & aB,
                        const G4double & aC,
                        const G4double & aD,
                        const G4double & aE,
                        const G4double & aF,
                        const G4double & aG,
                        const G4double & aH,
                        const G4double & aI,
                        const G4double & aJ,
                        const G4int    Sense = 1);
    //general Quadric equation:
    //    A*X^2 + 2*B*X*Y + 2*C*X*Z  + 2*D*X
    //  + E*Y^2 + 2*F*Y*Z            + 2*G*Y
    //  + H*Z^2                      + 2*I*Z   + J =0


    void SetParams(const G4double & aA,
                   const G4double & aB,
                   const G4double & aC,
                   const G4double & aD,
                   const G4double & aE,
                   const G4double & aF,
                   const G4double & aG,
                   const G4double & aH,
                   const G4double & aI,
                   const G4double & aJ);

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

    virtual G4double        HowNear( const G4ThreeVector& aPoint)const;
    //calculate the distance to the surface
    //if the point is inside the surface, return a positive value
    //Implement by subclass if possible

    virtual  G4HSSurfType   getType() const {return G4HSQuadric;}
    //return surface type

    G4String                printMe() const {std::stringstream strstr;
                                             strstr << "Quadric. A="<<A<<",\t B="<<B
                                                       <<",\t C="<<C
                                                         <<",\t D="<<D
                                                           <<",\t E="<<E
                                                             <<",\t F="<<F
                                                               <<",\t G="<<G
                                                                 <<",\t H="<<H
                                                                   <<",\t I="<<I
                                                                     <<",\t J="<<J<<std::endl;
                                            return strstr.str();}


    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                      std::vector <G4ThreeVector> & IntersectPoints,
                                      std::vector <G4double>  &     IntersectDistances) ;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances
    //ATTENTION: Remenber to Reset before using this method!
    //algorithm : see An Introduction to Ray Tracing

    G4HalfSpaceQuadric(const G4HalfSpaceQuadric & right);
    //copy constructor
    G4HalfSpaceQuadric& operator=(const G4HalfSpaceQuadric & right);
    //copy operator
    G4HalfSpaceQuadric operator-();
    //reverse Sense operator

protected:
    virtual G4double        Aq(const G4ThreeVector & aVector);
    //method to calcuate the Aq. the aVector is the vector of the ray
    //Aq = A*X^2 +  2*B*X*Y + 2*C*X*Z + E*Y^2 + 2*F*Y*Z  + H*Z^2
    //The children class can use more simple way to calculate it

    virtual G4double        Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector);
    //method to calcaulate the Bq,
    //Bq = 2* (A*X0*Xd + B*(X0*Yd+ Xd*Y0) + C*(X0*Zd+ Xd*Z0) + D*Xd +
    //         E*Y0*Yd + F*(Y0*Zd + Yd*Z0) + G*Yd + H*Z0*Zd + I*Zd)
    virtual G4double        Cq(const G4ThreeVector & aPoint);
    //Cq= A*X0^2 + 2*B*X0*Y0 + 2*C*X0*Z0 + 2*D*X0
    //  + E*Y0^2 +  2*F*Y0*Z0 + 2*G*Y0
    //    + H*Z0^2 + 2*I*Z0 + J

protected:

    //see general euqation, these parametere is not set in constructor
    G4double A;
    G4double E;
    G4double H;
    G4double B;
    G4double C;
    G4double F;
    G4double D;
    G4double G;
    G4double I;
    G4double J;


};

#endif // G4HALFSPACEQUADRIC_HH
