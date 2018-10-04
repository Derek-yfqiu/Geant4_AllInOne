#include "G4HalfSpaceEllipticalTorus.hh"

G4HalfSpaceEllipticalTorus::G4HalfSpaceEllipticalTorus(
        const G4ThreeVector & Center, const G4ThreeVector & Axis,
        const G4double & MaxRadius, const G4double & MinRadius1,
        const G4double & MinRadius2)
{
}

G4int  G4HalfSpaceEllipticalTorus::Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec)
{
    // Variables. Should be optimized later...
    G4ThreeVector Base = aPoint;   // Base of the intersection ray
    G4ThreeVector DCos = aVec.unit();     // Direction cosines of the ray
    G4int	     nhits =0;		      // Number of intersections
    G4double   rhits[4];		      // Intersection distances
    G4double   p, a0, b0;	      // Related constants
    G4double   f, l, t, g1, q, m1, u;   // Ray dependent terms
    G4double   C[5];		      // Quartic coefficients

    //calculate the relate constances
    p =( m_r1 * m_r1) / (m_r2 * m_r2);
    a0 = 4 * m_R * m_R;
    b0 = (m_R * m_R) - ( m_r1 * m_r1);

    //	Compute ray dependent terms.
    f = 1. - DCos.y()*DCos.y();
    g1 = f + p * DCos.y()*DCos.y();
    l = 2. * (Base.x()*DCos.x() + Base.z()*DCos.z());
    t = Base.x()*Base.x() + Base.z()*Base.z();
    q = a0 / (g1*g1);



}
