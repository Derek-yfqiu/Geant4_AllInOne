#include "G4HalfSpaceCone.hh"


//G4HalfSpaceCone::G4HalfSpaceCone( G4HSAxis thePllAxis, const G4double & T2,
//                const G4double & XYZpi, const G4int Sense )
//{
//    m_PllAxis = thePllAxis;
//    if (T2 > 0) m_T2 = T2;
//    else m_T2 = 0;

//    switch (m_PllAxis)
//     {
//     case AxisX:
//        m_Xpi = XYZpi;  m_Ypi = 0.0; m_Zpi = 0.0;
//        A= -T2; D= T2*XYZpi; E= 1;  /* G=-Ypi;*/    H=1;    /*I=-Zpi;*/    J= /*Ypi*Ypi + Zpi*Zpi*/- T2*XYZpi*XYZpi;
//        //B=C=F=G=I=0
//         break;
//     case AxisY:
//        m_Xpi = 0.0;  m_Ypi = XYZpi; m_Zpi = 0.0;
//        A= 1;   /*D= -Xpi;*/  E= -T2; G= T2*XYZpi; H= 1;   /*I= -Zpi;*/   J= /*Xpi*Xpi +  Zpi*Zpi*/ - T2*XYZpi*XYZpi;
//        //B=C=F=D=I=0
//         break;
//     case AxisZ:
//        m_Xpi = 0.0;  m_Ypi = 0.0; m_Zpi = XYZpi;
//        A= 1;   /*D= -Xpi; */ E= 1;   /*G=-Ypi;*/    H= -T2; I= T2*XYZpi; J= /*Xpi*Xpi + Ypi*Ypi*/ -T2*XYZpi*XYZpi;
//        //B=C=F=D=G=0
//         break;
//     }
//    if (Sense == 1 || Sense == -1) m_Sense = Sense;
//    else   G4Exception("G4HalfSpaceCone::G4HalfSpaceCone", "GeomSolids1003",
//                             FatalException, "The sense should be 1 or -1!");
//}

G4HalfSpaceCone::G4HalfSpaceCone()
{
}

G4HalfSpaceCone::G4HalfSpaceCone( G4HSAxis thePllAxis, const G4double & T2,
                                  const G4double & Xpi, const G4double & Ypi,
                                  const G4double & Zpi, const G4int Sense)
{
    m_PllAxis = thePllAxis;
    if (T2 > 0) m_T2 = T2;
    else m_T2 = 0;
    m_Xpi = Xpi;
    m_Ypi = Ypi;
    m_Zpi = Zpi;

    switch (m_PllAxis)
     {
     case AxisX:
        A= -T2; D= T2*Xpi; E= 1;   G=-Ypi;    H=1;    I=-Zpi;    J= Ypi*Ypi + Zpi*Zpi - T2*Xpi*Xpi;
        //B=C=F=0
         break;
     case AxisY:
        A= 1;   D= -Xpi;   E= -T2; G= T2*Ypi; H= 1;   I= -Zpi;   J= Xpi*Xpi +  Zpi*Zpi - T2*Ypi*Ypi;
        //B=C=F=0
         break;
     case AxisZ:
        A= 1;   D= -Xpi;   E= 1;   G=-Ypi;    H= -T2; I= T2*Zpi; J= Xpi*Xpi + Ypi*Ypi -T2*Zpi*Zpi;
        //B=C=F=0
         break;
     }
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceCone::G4HalfSpaceCone", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}


G4ThreeVector   G4HalfSpaceCone::Normal(G4ThreeVector aPoint) const
{
    G4ThreeVector aNormal;
    switch (m_PllAxis)
     {
     case AxisX:
        aNormal.setX(A*aPoint.x() + /*B*aPoint.y() + C*aPoint.z() + */D ); //*2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z()*/ + G);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() + I);//H=1; *2 is ignored
         break;
     case AxisY:
        aNormal.setX(  aPoint.x() + /*B*aPoint.y() + C*aPoint.z() +*/ D); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ E*aPoint.y() + /*F*aPoint.z() +*/ G );//*2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() + I);//H=1; *2 is ignored
         break;
     case AxisZ:
        aNormal.setX(  aPoint.x() /*+ B*aPoint.y() + C*aPoint.z()*/ + D); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z()*/ + G);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y()  +*/ H*aPoint.z() + I);//*2 is ignored
         break;
     }
    if (m_Sense == 1) return -aNormal.unit();
    else return aNormal.unit();
    //the normal is possible to be (0,0,0)
}

EInside G4HalfSpaceCone::Inside(const G4ThreeVector & aPoint) const
{
    G4double aResult =0.;

    switch (m_PllAxis)
     {
     case AxisX:
        //B=C=F=0
        aResult =A*aPoint.x()*aPoint.x() +
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
                            2*D*aPoint.x() +
                                aPoint.y()*aPoint.y() +  //E=1
        //                    2*F*aPoint.y()*aPoint.z() +
                            2*G*aPoint.y() +
                                aPoint.z()*aPoint.z() +  //H=1
                            2*I*aPoint.z()   +
                              J;
        break;
     case AxisY:
        //B=C=F=0
        aResult =  aPoint.x()*aPoint.x() + //A=1
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
                            2*D*aPoint.x() +
                              E*aPoint.y()*aPoint.y() +
        //                    2*F*aPoint.y()*aPoint.z() +
                            2*G*aPoint.y() +
                                aPoint.z()*aPoint.z() + //H=1
                            2*I*aPoint.z()   +
                              J;
        break;
     case AxisZ:
        //B=C=F=0
        aResult =  aPoint.x()*aPoint.x() + //A=1
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
                            2*D*aPoint.x() +
                                aPoint.y()*aPoint.y() +  //E=1
        //                    2*F*aPoint.y()*aPoint.z() +
                            2*G*aPoint.y() +
                              H*aPoint.z()*aPoint.z() +
                            2*I*aPoint.z()   +
                              J;
        break;
     }

//no a good way    G4double aResult = HowNear(aPoint);
    aResult = aResult * m_Sense; //adjust by sense

    /*
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
        */

    // **we use flexiale way to judge the position of the point**
    if ((abs(aResult) > 100000*kCarTolerance ) || (abs(aResult) < 100*kCarTolerance )) {
        //for those "envident" results, we use rough calculation results
        if (aResult > (kCarTolerance*100)) //positive means inside half-space
            return kInside;
        else if (aResult < (-kCarTolerance*100))  //negative means outside half-space
            return kOutside;
        else
            return kSurface;

    }
    else {
        //else we calculate the actually distance to the cone surface
        aResult = HowNear(aPoint);
        if (aResult > (kCarTolerance*100)) //positive means inside half-space
            return kInside;
        else if (aResult < (-kCarTolerance*100))  //negative means outside half-space
            return kOutside;
        else
            return kSurface;
    }
}

G4double        G4HalfSpaceCone::HowNear( const G4ThreeVector& aPoint)const
{
    //see General Design document for the detail
    G4double aDist = 0.0;
    //2016-03-30 should be absolute value, otherwise result has no "sense"
    G4double DX = abs(aPoint.x() - m_Xpi);//aPoint.x() - m_Xpi;
    G4double DY = abs(aPoint.y() - m_Ypi);//aPoint.y() - m_Ypi;
    G4double DZ = abs(aPoint.z() - m_Zpi);//aPoint.z() - m_Zpi;
    G4double R2, R ,O,  N, P ;
//take the cone parallel to Z-axis as example
//    R^2 = t2 * (Z-Z’)2
//    Distance from point to cone axis O= Sqrt(Dx*Dx + Dy*Dy),
//    N= O – R
//    M/N = L/P
//    P= sqrt(R^2 + L^2)
//    We get M= N* L/ sqrt(R^2 + L^2)
//    This is also applicable when point is inside (negative value)

    switch (m_PllAxis)
    {
    case AxisX:
        R2 = m_T2* DX*DX;
        R = sqrt(R2);
        O = sqrt(DY*DY + DZ*DZ);
        N = O - R;
        //here L = DX
        P = sqrt(R2 + DX*DX);
        aDist = N * DX / P;
        break;
    case AxisY:
        R2 = m_T2* DY*DY;
        R = sqrt(R2);
        O = sqrt(DX*DX + DZ*DZ);
        N = O - R;
        //here L = DY
        P = sqrt(R2 + DY*DY);
        aDist = N * DY / P;
        break;
    case AxisZ:
        R2 = m_T2* DZ*DZ;
        R = sqrt(R2);
        O = sqrt(DX*DX + DY*DY);
        N = O - R;
        //here L = DZ
        P = sqrt(R2 + DZ*DZ);
        aDist = N * DZ / P;
        break;
    }

    //the point inside the cone is < 0;
    //if m_Sense = 1, the cone represents the outside half-space,
    //therefore it should multiply with Sense
    return aDist * m_Sense;//??
}

G4double  G4HalfSpaceCone::Aq(const G4ThreeVector & aVector)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        aResult = A*aVector.x()*aVector.x() +
//                2*B*aVector.x()*aVector.y() +
//                2*C*aVector.x()*aVector.z() +
                    aVector.y()*aVector.y() + //E=1
//                2*F*aVector.y()*aVector.z()  +
                    aVector.z()*aVector.z(); //H=1
         break;
     case AxisY:
        aResult =   aVector.x()*aVector.x() + //A=1
//                2*B*aVector.x()*aVector.y() +
//                2*C*aVector.x()*aVector.z() +
                  E*aVector.y()*aVector.y() +
//                2*F*aVector.y()*aVector.z()  +
                    aVector.z()*aVector.z(); //H=1
         break;
     case AxisZ:
        aResult =   aVector.x()*aVector.x() + //A=1
//                2*B*aVector.x()*aVector.y() +
//                2*C*aVector.x()*aVector.z() +
                    aVector.y()*aVector.y() + //E=1
//                2*F*aVector.y()*aVector.z()  +
                  H*aVector.z()*aVector.z();
         break;
     }

    return aResult;
}

G4double   G4HalfSpaceCone:: Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        aResult = 2* (  A*aPoint.x()*aVector.x() +
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
                      D*aVector.x() +
                        aPoint.y()*aVector.y() + //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
                      G*aVector.y() +
                        aPoint.z()*aVector.z() +  //H=1
                      I*aVector.z());
         break;
     case AxisY:
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
                      D*aVector.x() +
                      E*aPoint.y()*aVector.y() +
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
                      G*aVector.y() +
                        aPoint.z()*aVector.z() +  //H=1
                      I*aVector.z());
         break;
     case AxisZ:
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
                      D*aVector.x() +
                        aPoint.y()*aVector.y() + //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
                      G*aVector.y() +
                      H*aPoint.z()*aVector.z() +
                      I*aVector.z());
         break;
     }
    return aResult;
}

G4double        G4HalfSpaceCone:: Cq(const G4ThreeVector & aPoint)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        aResult = A*aPoint.x()*aPoint.x() +
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() +//E=1
//                2*F*aPoint.y()*aPoint.z() +
                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
                2*I*aPoint.z() +
                  J ;
         break;
     case AxisY:
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
                2*D*aPoint.x()+
                  E*aPoint.y()*aPoint.y() +
//                2*F*aPoint.y()*aPoint.z() +
                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
                2*I*aPoint.z() +
                  J ;
         break;
     case AxisZ:
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() + //E=1
//                2*F*aPoint.y()*aPoint.z() +
                2*G*aPoint.y() +
                  H*aPoint.z()*aPoint.z() +
                2*I*aPoint.z() +
                  J ;
         break;
     }
    return aResult;

}



//copy constructor
G4HalfSpaceCone::G4HalfSpaceCone(const G4HalfSpaceCone & right)
    :G4HalfSpaceQuadric(right)
{
    m_PllAxis = right.m_PllAxis;
    m_T2 = right.m_T2;
    m_Xpi = right.m_Xpi;
    m_Ypi = right.m_Ypi;
    m_Zpi = right.m_Zpi;
}

//copy operator
G4HalfSpaceCone& G4HalfSpaceCone::operator=(const G4HalfSpaceCone & right)
{
    if(&right == this) return *this;

    //base class
    G4HalfSpaceQuadric::operator =(right);

    m_PllAxis = right.m_PllAxis;
    m_T2 = right.m_T2;
    m_Xpi = right.m_Xpi;
    m_Ypi = right.m_Ypi;
    m_Zpi = right.m_Zpi;
    return *this;
}

//reverse Sense operator
G4HalfSpaceCone G4HalfSpaceCone::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceCone aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}
