#ifndef G4HALFSPACETORUS_HH
#define G4HALFSPACETORUS_HH
#include "G4HalfSpaceSurface.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
/*!
 * \brief The G4HalfSpaceTorus class
 *  torus with circle cross section
 */
class G4HalfSpaceTorus : public G4HalfSpaceSurface
{
public:
    G4HalfSpaceTorus(const G4ThreeVector & Center, const G4ThreeVector & Axis,
                     const G4double & MaxRadius, const G4double & MinRadius,
                     const G4int Sense = 1);
    //contructor for circular torus
    //Center: the torus center
    //Axis: the axis
    //MaxRadius: the large radius
    //MinRadius: the small radius
    //ATTENTION: by default, the torus has center (0,0,0) and axis (0,0,1)

    G4HalfSpaceTorus(const G4ThreeVector & Center, const G4ThreeVector & Axis,
                     const G4double & MaxRadius, const G4double & MinRadius1,
                     const G4double & MinRadius2,
                     const G4int Sense = 1);
    //contructor for elliptical torus
    //contructor for circular torus
    //Center: the torus center
    //Axis: the axis
    //MaxRadius: the large radius
    //MinRadius1: the small ellipse radius in radial direction
    //MinRadius2: the small ellipse radius in axis direction

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

    virtual  G4HSSurfType   getType() const {return G4HSTorus;}
    //return surface type

    G4String                printMe() const {std::stringstream strstr;
                                             strstr << "Torus. CenterX="<< m_Center.X<<
                                                       "\t CenterY=" <<m_Center.Y<<
                                                       "\t CenterZ=" <<m_Center.Z<<
                                                       "\t AxisX=" <<m_Axis.X<<
                                                       "\t AxisY=" <<m_Axis.Y<<
                                                       "\t AxisZ=" <<m_Axis.Z<<
                                                       "\t Outer Radius="<<m_R<<
                                                       "\t Inner Radius1="<<m_r1<<
                                                       "\t Inner Radius2="<<m_r2<<std::endl;
                                            return strstr.str();}


    virtual G4int           Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                      std::vector <G4ThreeVector> & IntersectPoints,
                                      std::vector <G4double>  &     IntersectDistances) ;
    //calculat the intersection of the ray with the surface
    //return the number of intersection points and update m_IntersectPoints and m_IntersectDistances
    //ATTENTION: Remenber to Reset before using this method!
    //algorithm : from G4BREPToroidalSurface

    G4HalfSpaceTorus(const G4HalfSpaceTorus & right);
    //copy constructor
    G4HalfSpaceTorus& operator=(const G4HalfSpaceTorus & right);
    //copy operator
    G4HalfSpaceTorus operator-();
    //reverse Sense operator

protected:
    G4int SolveQuartic(G4double c[], G4double s[]);
    G4int SolveCubic  (G4double c[], G4double s[]);
    G4int SolveQuadric(G4double c[], G4double s[]);
      // Algorithms for solving quadratic, cubic and quartic equations.
    G4int IsZero(G4double x) const;

protected:

    G4ThreeVector           m_Center;
    //center Point of the torus
    G4ThreeVector           m_Axis;
    //rotation axis of the torus
    G4double                m_R;
    //the large radius
    G4double                m_r1;
    //the first small radius
    G4double                m_r2;
    //the second small radius
    G4AffineTransform           m_Transfrom;
    //transforming default torus with center (0,0,0) axis (0,0,1)
    //to m_Center and m_Axis
    G4AffineTransform           m_InvertTransfrom;
    //Inverted transformation from m_Transfrom
    //calculate in advance to avoid recalculate in functions
    const G4double     EQN_EPS;
    //for judge if zero

};

#endif // G4HALFSPACETORUS_HH
