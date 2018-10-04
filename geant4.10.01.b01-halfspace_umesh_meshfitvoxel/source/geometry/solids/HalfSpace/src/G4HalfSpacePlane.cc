#include "G4HalfSpacePlane.hh"

//G4HalfSpacePlane::G4HalfSpacePlane()
//{
//}

/*!
 * \brief G4HalfSpacePlane::G4HalfSpacePlane
 *  constructor with parameter
 * \param aA parameter for x
 * \param aB parameter for y
 * \param aC parameter for z
 * \param aD distance to origin
 */
G4HalfSpacePlane::G4HalfSpacePlane(const G4double aA, const G4double aB,
                 const G4double aC, const G4double aD, const G4int Sense)
{
    m_Plane.a = aA;
    m_Plane.b = aB;
    m_Plane.c = aC;
    m_Plane.d = aD;
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpacePlane::G4HalfSpacePlane", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
    //the normal should be changed according to the sense of the plane!
    //if the sense is 1, the normal should pointing to outside (negative side)
    //therefore it should multiply with -1.
    m_Normal = G4ThreeVector (m_Plane.a, m_Plane.b, m_Plane.c).unit() * (-m_Sense)  ;

}

G4HalfSpacePlane::G4HalfSpacePlane(G4HSAxis theAxis, const G4double aD,
                 const G4int    Sense )

{
    switch (theAxis)
    {
    case AxisX:
        m_Plane.a = 1.;
        m_Plane.b = 0.;
        m_Plane.c = 0.;
        m_Plane.d = aD;
        break;
    case AxisY:
        m_Plane.a = 0.;
        m_Plane.b = 1.;
        m_Plane.c = 0.;
        m_Plane.d = aD;
        break;
    case AxisZ:
        m_Plane.a = 0.;
        m_Plane.b = 0.;
        m_Plane.c = 1.;
        m_Plane.d = aD;
        break;
    }
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpacePlane::G4HalfSpacePlane", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
    //the normal should be changed according to the sense of the plane!
    //if the sense is 1, the normal should pointing to outside (negative side)
    //therefore it should multiply with -1.
    m_Normal = G4ThreeVector (m_Plane.a, m_Plane.b, m_Plane.c).unit() * (-m_Sense)  ;

}

/*!
 * \brief G4HalfSpacePlane::Normal
 *  return the normal of the plane
 * \param aPoint not applicable here
 * \return a unit vector of the normal
 */
G4ThreeVector   G4HalfSpacePlane::Normal(G4ThreeVector) const
{   //this function have no checking
   return m_Normal;
}

/*!
 * \brief G4HalfSpacePlane::Inside
 *  judge the position of the point, if the point on the solid side
 * \param aPoint the point to be judge
 * \return return the position
 */
EInside   G4HalfSpacePlane::Inside(const G4ThreeVector & aPoint) const
{
    //f(x,y,z) < 0 : inside; >0: outside; = 0: on surface
//    G4double aDist = (m_Plane.a * aPoint.x() +
//            m_Plane.b * aPoint.y() +
//            m_Plane.c * aPoint.z() + m_Plane.d )
//            * m_Sense; //adjust the sign according to the solid
    G4double aDist = HowNear(aPoint);
    if (aDist > kCarTolerance) //positive means inside
        return kInside;
    else if (aDist < -kCarTolerance)  //negative means outside
        return kOutside;
    else
        return kSurface;
}

/*!
 * \brief G4HalfSpacePlane::HowNear
 *  calculate the distance to the surface
 *  if the point is inside the surface, return a positive value
 * \param aPoint the point to be judge
 * \return the distance
 */
G4double  G4HalfSpacePlane::HowNear(const G4ThreeVector &aPoint )const
{
    //the sign of the distance should take the sense into account
    return  ( aPoint.x()*m_Plane.a + aPoint.y()*m_Plane.b +
            aPoint.z()*m_Plane.c - m_Plane.d ) *m_Sense;
}

/*!
 * \brief Intersect
 *  calculat the intersection of the ray with the surface
 *  return the number of intersection points and update m_IntersectPoints and m_IntersectDistances
 *  \sa  An Introduction to Ray Tracing, chapter 2
 * \param aRay a ray
 * \return  number of intersections
 */
G4int G4HalfSpacePlane::Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                  std::vector<G4ThreeVector> &IntersectPoints, std::vector<G4double> &IntersectDistances)
{
    G4ThreeVector aUnitVec = aVec.unit() ; //normalize the vector

    //calculat the Normal*aVec to see if parallel
    G4double aVd = aUnitVec.x() * m_Plane.a +
            aUnitVec.y() * m_Plane.b +
            aUnitVec.z() * m_Plane.c ;
    if (std::fabs(aVd) >= kCarTolerance)  //if not parallel
    {
        //calculate the intersect distance
        G4double aV0 = -(aPoint.x() * m_Plane.a +
                         aPoint.y() * m_Plane.b +
                         aPoint.z() * m_Plane.c - m_Plane.d);
        G4double aT = aV0/ aVd;  //because the ray vector is normalized, this T is distance

        //see if the ray travel against the plane
        //calculate the solution //from BREP/G4FPlane
        G4double solx,soly,solz;
        solx = aPoint.x() + aT* aUnitVec.x();
        soly = aPoint.y() + aT* aUnitVec.y();
        solz = aPoint.z() + aT* aUnitVec.z();

        // solve tolerance problem
        if( (aT*aUnitVec.x() >= -kCarTolerance/2) && (aT*aUnitVec.x() <= kCarTolerance/2) )
          solx = aPoint.x();

        if( (aT*aUnitVec.y() >= -kCarTolerance/2) && (aT*aUnitVec.y() <= kCarTolerance/2) )
          soly = aPoint.y();

        if( (aT*aUnitVec.z() >= -kCarTolerance/2) && (aT*aUnitVec.z() <= kCarTolerance/2) )
          solz = aPoint.z();

        G4bool xhit = (aUnitVec.x() < 0 && solx <= aPoint.x()) || (aUnitVec.x() >= 0 && solx >= aPoint.x());
        G4bool yhit = (aUnitVec.y() < 0 && soly <= aPoint.y()) || (aUnitVec.y() >= 0 && soly >= aPoint.y());
        G4bool zhit = (aUnitVec.z() < 0 && solz <= aPoint.z()) || (aUnitVec.z() >= 0 && solz >= aPoint.z());

        if( xhit && yhit && zhit ) {
//            m_IntersectPoints.push_back(G4ThreeVector(solx, soly, solz));
//            m_IntersectDistances.push_back(aT);
            IntersectPoints.push_back(G4ThreeVector(solx, soly, solz));
            IntersectDistances.push_back(aT);
            return 1;
        }
        else return 0; //not intersection

    }
    return 0;
}

//copy constructor
G4HalfSpacePlane::G4HalfSpacePlane(const G4HalfSpacePlane & right)
    :G4HalfSpaceSurface(right)
{
    //for this class
    m_Plane  = right.m_Plane;
    m_Normal = right.m_Normal;

}

//copy operator
G4HalfSpacePlane& G4HalfSpacePlane::operator=(const G4HalfSpacePlane & right)
{
    if(&right == this) return *this;
    //G4HalfSpaceSurface general
    G4HalfSpaceSurface::operator =(right);

    //for this class
    m_Plane  = right.m_Plane;
    m_Normal = right.m_Normal;
    return *this;
}

//reverse Sense operator
G4HalfSpacePlane G4HalfSpacePlane::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpacePlane aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}

