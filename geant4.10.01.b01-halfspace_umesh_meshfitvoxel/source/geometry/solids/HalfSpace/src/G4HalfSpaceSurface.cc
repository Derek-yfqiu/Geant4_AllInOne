#include "G4HalfSpaceSurface.hh"
#include "G4GeometryTolerance.hh"

G4HalfSpaceSurface::G4HalfSpaceSurface()
{
    m_Sense = 1; //default positive sense
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
    m_isActive = true;
}

void G4HalfSpaceSurface::Reset()
{
    m_isActive = true;
    m_IntersectPoints.clear();
    m_IntersectDistances.clear();
}

G4double  G4HalfSpaceSurface:: HowNear( const G4ThreeVector&  )const
{
    //default no calculation of the distance
    return 0.;
}
//copy constructor
G4HalfSpaceSurface::G4HalfSpaceSurface(const G4HalfSpaceSurface & right)
{
    //G4HalfSpaceSurface general
    m_Sense =right.m_Sense;
    m_isActive = true;
    kCarTolerance = right.kCarTolerance;
    kAngTolerance = right.kAngTolerance;
}

//copy operator
G4HalfSpaceSurface& G4HalfSpaceSurface::operator=(const G4HalfSpaceSurface & right)
{
    if(&right == this) return *this;
    //G4HalfSpaceSurface general
    m_Sense =right.m_Sense;
    m_isActive = true;
    kCarTolerance = right.kCarTolerance;
    kAngTolerance = right.kAngTolerance;

    return *this;
}
