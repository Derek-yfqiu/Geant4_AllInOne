#include "G4HalfSpaceSphere.hh"

G4HalfSpaceSphere::G4HalfSpaceSphere(const G4double & CenterX,
                                     const G4double & CenterY,
                                     const G4double & CenterZ,
                                     const G4double & Radius,
                                     const G4int    Sense)
{
    m_Center.set(CenterX, CenterY, CenterZ) ;
    if (Radius > 0.0) m_Radius = Radius;
    else m_Radius = 0.0;
    m_Radius2 = m_Radius * m_Radius;
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceSphere::G4HalfSpaceSphere", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}

G4ThreeVector   G4HalfSpaceSphere::Normal(G4ThreeVector aPoint) const
{
    //if the point is the center, then no way to calcualte the normal
    if (aPoint == m_Center) return G4ThreeVector(0., 0., 0.);
    G4ThreeVector aNormal = m_Center - aPoint;
    if (m_Sense == 1) return aNormal.unit() ;
    else return -aNormal.unit();
}

EInside G4HalfSpaceSphere::Inside(const G4ThreeVector & aPoint) const
{
    G4double aDist = HowNear(aPoint);
    if (aDist > kCarTolerance) //positive means inside
        return kInside;
    else if (aDist < -kCarTolerance)  //negative means outside
        return kOutside;
    else
        return kSurface;
}

G4double  G4HalfSpaceSphere::HowNear( const G4ThreeVector& aPoint)const
{
 //calculate the distance to the sphere surface
    G4double aDist = (m_Center - aPoint).mag() - m_Radius;
    //the points "Inside" has positive value;
    //when Sense = 1, the half-space it represent is outer space, therefore
    //point has positive distance to the sphere surface return positive value;
    //vice versa
    return aDist * m_Sense;
}

/*!
 * \brief G4HalfSpaceSphere::Intersect
 *  intersection between a ray and the sphere surface
 * \param aPoint the Origin of the ray
 * \param aVec the direction of the ray
 * \return the number of intersections
 * \sa  Book: An Introduction to Ray Tracing
 * \todo the tolerance treatment is not accurate enough
 */
G4int G4HalfSpaceSphere::Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                   std::vector<G4ThreeVector> &IntersectPoints, std::vector<G4double> &IntersectDistances)
{
    //check first
    if (m_Center.x() >= kInfinity || m_Center.y() >= kInfinity || m_Center.z() >= kInfinity ||
            m_Radius <= kCarTolerance) {
        G4Exception(" G4HalfSpaceSphere::Intersect", "GeomSolids1003",
                    FatalException, "The Center is Infinit far away or the radius is 0!");
    }
    //preparation
    G4ThreeVector aUnitVec = aVec.unit();
    G4double Tol2 = kCarTolerance * kCarTolerance;
    //Calculate the distance from ray origin to sphere center
    G4ThreeVector aOC = m_Center - aPoint;
    G4double aLoc2 = aOC.mag2(); //square distance

    //details see Book: An Introduction to Ray Tracing, and also the detail design document
    if (aLoc2 >= m_Radius2) //if the origin of the ray is outside the sphere
    {
        //Calculate the distance from sphere center to the ray (point-line distance)
        G4double aTca = aOC.dot(aUnitVec);
        //if aTca < 0 the ray is going way -- no intersection
        //Tolerance treament: if aTca <= kCarTolerance and >= -kCarTolerance,
        //the ray origin is on the tolerance layer, and the ray just slip -- no intersection
        if (aTca <= kCarTolerance)
            return 0;
        //calculate the half-length of the ray segment inside the sphere surface
        //we calculate the square distance first
        G4double aThc2 = m_Radius2 - (aLoc2 - aTca* aTca);
        //if aThc2 <0, no real solution for the intersect points
        //tolerance treament: if the ray just slip, --no intersection
        if (aThc2 <= Tol2)
            return 0;

        //now calculate the intersect points and distance
        G4double aThc = sqrt(aThc2);
        G4double aT1 = aTca - aThc; //solution 1
        G4double aT2 = aTca + aThc; // solution 2
        G4int aNbSln = 0; //number of valid solution

        //check if intersect distance is valid
        if (aT1 >= kCarTolerance && aT1 < kInfinity) {
//            m_IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT1,
//                                                      aPoint.y() + aUnitVec.y()*aT1,
//                                                      aPoint.z() + aUnitVec.z()*aT1));
//            m_IntersectDistances.push_back(aT1);
            IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT1,
                                                      aPoint.y() + aUnitVec.y()*aT1,
                                                      aPoint.z() + aUnitVec.z()*aT1));
            IntersectDistances.push_back(aT1);
            aNbSln++;
        }
        if (aT2 >= kCarTolerance && aT2 < kInfinity) {
//            m_IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT2,
//                                                      aPoint.y() + aUnitVec.y()*aT2,
//                                                      aPoint.z() + aUnitVec.z()*aT2));
//            m_IntersectDistances.push_back(aT2);
            IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT2,
                                                      aPoint.y() + aUnitVec.y()*aT2,
                                                      aPoint.z() + aUnitVec.z()*aT2));
            IntersectDistances.push_back(aT2);
            aNbSln++;
        }
        return aNbSln;
    }
    else  //inside
    {
        G4double aTca = aOC.dot(aUnitVec);
        G4double aThc2 = m_Radius2 - (aLoc2 - aTca* aTca);
        if (aThc2 <= Tol2)  //on the tolerance layer
            return 0;
        //only the solution 2 is valid when point is inside
        G4double aT2 = aTca + sqrt(aThc2); // solution 2
        if (aT2 >= kCarTolerance && aT2 < kInfinity) {
//            m_IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT2,
//                                                      aPoint.y() + aUnitVec.y()*aT2,
//                                                      aPoint.z() + aUnitVec.z()*aT2));
//            m_IntersectDistances.push_back(aT2);
            IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT2,
                                                      aPoint.y() + aUnitVec.y()*aT2,
                                                      aPoint.z() + aUnitVec.z()*aT2));
            IntersectDistances.push_back(aT2);
            return 1;
        }
        else return 0;
    }
}
//copy constructor
G4HalfSpaceSphere::G4HalfSpaceSphere(const G4HalfSpaceSphere & right)
    :G4HalfSpaceSurface(right)
{
    //for this class
    m_Center = right.m_Center;
    m_Radius = right.m_Radius;
    m_Radius2 =right.m_Radius2;

}

//copy operator
G4HalfSpaceSphere& G4HalfSpaceSphere::operator=(const G4HalfSpaceSphere & right)
{
    if(&right == this) return *this;
    //G4HalfSpaceSurface general
    G4HalfSpaceSurface::operator =(right);

    //for this class
    m_Center = right.m_Center;
    m_Radius = right.m_Radius;
    m_Radius2 =right.m_Radius2;
    return *this;
}

//reverse Sense operator
G4HalfSpaceSphere G4HalfSpaceSphere::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceSphere aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}
