#include "G4HalfSpaceQuadric.hh"

G4HalfSpaceQuadric::G4HalfSpaceQuadric()
{
    A = 0.0;
    E = 0.0;
    H = 0.0;
    B = 0.0;
    C = 0.0;
    F = 0.0;
    D  = 0.0;
    G  = 0.0;
    I  = 0.0;
    J  = 0.0;
}

G4HalfSpaceQuadric::G4HalfSpaceQuadric (const G4double & aA,
                                        const G4double & aB,
                                        const G4double & aC,
                                        const G4double & aD,
                                        const G4double & aE,
                                        const G4double & aF,
                                        const G4double & aG,
                                        const G4double & aH,
                                        const G4double & aI,
                                        const G4double & aJ,
                                        const G4int    Sense)
{
    SetParams(aA, aB, aC, aD, aE, aF, aG, aH, aI, aJ);
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceQuadric::G4HalfSpaceQuadric", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");

}

void G4HalfSpaceQuadric::SetParams(const G4double & aA,
               const G4double & aB,
               const G4double & aC,
               const G4double & aD,
               const G4double & aE,
               const G4double & aF,
               const G4double & aG,
               const G4double & aH,
               const G4double & aI,
               const G4double & aJ)
{
    A =  aA;
    E =  aE;
    H =  aH;
    B =  aB;
    C =  aC;
    F =  aF;
    D  =  aD;
    G  =  aG;
    I  =  aI;
    J  =  aJ;
}

G4ThreeVector   G4HalfSpaceQuadric::Normal(G4ThreeVector aPoint) const
{
    G4ThreeVector aNormal;
    aNormal.setX(A*aPoint.x() + B*aPoint.y() + C*aPoint.z() + D); //*2 is ignored
    aNormal.setY(B*aPoint.x() + E*aPoint.y() + F*aPoint.z() + G);//*2 is ignored
    aNormal.setZ(C*aPoint.x() + F*aPoint.y() + H*aPoint.z() + I);//*2 is ignored

    //take ellipse as example, the calculated normal is pointing outward
    //if m_sense =1, the half-space it represent is outside world,
    //To make normal pointing "Outside", the normal should be reversed
    if (m_Sense == 1) return -aNormal.unit();
    else return aNormal.unit();
}

EInside G4HalfSpaceQuadric::Inside(const G4ThreeVector & aPoint) const
{
    //we substitute the point into the quadric equation,
    //if the result >0 : Inside the half-space
    //if result<0: Outside the half-space
    //if result=0: on the surface
    G4double aResult =A*aPoint.x()*aPoint.x() +
                    2*B*aPoint.x()*aPoint.y() +
                    2*C*aPoint.x()*aPoint.z() +
                    2*D*aPoint.x() +
                      E*aPoint.y()*aPoint.y() +
                    2*F*aPoint.y()*aPoint.z() +
                    2*G*aPoint.y() +
                      H*aPoint.z()*aPoint.z() +
                    2*I*aPoint.z()   +
                      J;
//no a good way    G4double aResult = HowNear(aPoint);
    aResult = aResult * m_Sense; //adjust by sense
    //2016-03-24 multiplying 100 is because the  "aResult" is not the actual distance to the surface
    //instead it is just an estimation by substitute the point coordinate to the surface equation.
    //it encounter failure in the situation that the point is actually on the surface
    //therefore here we enlarge the tolerance to make it pass
    //2016-03-29 the tolerance is *1000 to handle some failure on this surface type
    if (aResult > (kCarTolerance*1000)) //positive means inside half-space
        return kInside;
    else if (aResult < (-kCarTolerance*1000))  //negative means outside half-space
        return kOutside;
    else
        return kSurface;
}

G4double        G4HalfSpaceQuadric::HowNear( const G4ThreeVector& )const
{
    //because the computation of the isotropic distance to General qualdirc is expensive
    //we use the most conservative value -- 0;
    return 0.;
}

G4double  G4HalfSpaceQuadric::Aq(const G4ThreeVector & aVector)
{
    return    A*aVector.x()*aVector.x() +
            2*B*aVector.x()*aVector.y() +
            2*C*aVector.x()*aVector.z() +
              E*aVector.y()*aVector.y() +
            2*F*aVector.y()*aVector.z() +
              H*aVector.z()*aVector.z();
}

G4double   G4HalfSpaceQuadric:: Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector)
{
     return 2* (A*aPoint.x()*aVector.x() +
                B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
                C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
                D*aVector.x() +
                E*aPoint.y()*aVector.y() +
                F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
                G*aVector.y() +
                H*aPoint.z()*aVector.z() +
                I*aVector.z());
}

G4double        G4HalfSpaceQuadric:: Cq(const G4ThreeVector & aPoint)
{
   return    A*aPoint.x()*aPoint.x() +
           2*B*aPoint.x()*aPoint.y() +
           2*C*aPoint.x()*aPoint.z() +
           2*D*aPoint.x()+
             E*aPoint.y()*aPoint.y() +
           2*F*aPoint.y()*aPoint.z() +
           2*G*aPoint.y() +
             H*aPoint.z()*aPoint.z() +
           2*I*aPoint.z() +
             J ;
}


G4int G4HalfSpaceQuadric:: Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                     std::vector<G4ThreeVector> &IntersectPoints, std::vector<G4double> &IntersectDistances)
{
    G4ThreeVector aUnitVec = aVec.unit(); //normalize
    //calculate parameters
    G4double aAq = Aq(aUnitVec);
    G4double aBq = Bq(aPoint, aUnitVec);
    G4double aCq = Cq(aPoint);
    //check the Aq
    if (aAq == 0) {
        G4double aT = - aCq / aBq ;
        //deal with tolerance, < Tolerence then ignore this intersect
        if (aT < kCarTolerance || aT >= kInfinity)
            return 0;
//        m_IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT,
//                                                  aPoint.y() + aUnitVec.y()*aT,
//                                                  aPoint.z() + aUnitVec.z()*aT));
//        m_IntersectDistances.push_back(aT);// the aT is just the distance, because aUnitVec is normalized
        IntersectPoints.push_back(G4ThreeVector(aPoint.x() + aUnitVec.x()*aT,
                                                  aPoint.y() + aUnitVec.y()*aT,
                                                  aPoint.z() + aUnitVec.z()*aT));
        IntersectDistances.push_back(aT);// the aT is just the distance, because aUnitVec is normalized

        return 1;
    }
    else  {
        G4double aSquare = aBq*aBq - 4*aAq*aCq;
        //if aSquare < 0no intersect
        if (aSquare <0) return 0;
        //if aSquare = 0, slip the surface, also no intersection
        //Tolerance treatment: t2-t1 = 2*sqrt(aSquare) should < 2*kCarTolerance
        //therefore aSquare should > kCarTolerance^2
        if (aSquare < kCarTolerance * kCarTolerance) return 0;

        //now calculate the intersect points and distance
        G4double aSqrt = sqrt(aSquare);
        G4double aT1 = (-aBq - aSqrt)/(2*aAq); //solution 1 //might be negative
        G4double aT2 = (-aBq + aSqrt)/(2*aAq); // solution 2
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
}


//copy constructor
G4HalfSpaceQuadric::G4HalfSpaceQuadric(const G4HalfSpaceQuadric & right)
    :G4HalfSpaceSurface(right)
{
    //for this class
    A =  right.A;
    E =  right.E;
    H =  right.H;
    B =  right.B;
    C =  right.C;
    F =  right.F;
    D  = right.D;
    G  = right.G;
    I  = right.I;
    J  = right.J;

}

//copy operator
G4HalfSpaceQuadric& G4HalfSpaceQuadric::operator=(const G4HalfSpaceQuadric & right)
{
    if(&right == this) return *this;
    //G4HalfSpaceSurface general
    G4HalfSpaceSurface::operator =(right);

    //for this class
    A =  right.A;
    E =  right.E;
    H =  right.H;
    B =  right.B;
    C =  right.C;
    F =  right.F;
    D  = right.D;
    G  = right.G;
    I  = right.I;
    J  = right.J;
    return *this;
}

//reverse Sense operator
G4HalfSpaceQuadric G4HalfSpaceQuadric::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceQuadric aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}
