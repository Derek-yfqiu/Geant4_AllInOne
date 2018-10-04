#include "G4HalfSpaceConeOnAxis.hh"

G4HalfSpaceConeOnAxis::G4HalfSpaceConeOnAxis(G4HSAxis thePllAxis,
                                             const G4double & T2,
                                             const G4double &XYZpi,
                                             const G4int Sense )
{
    m_PllAxis = thePllAxis;
    if (T2 > 0) m_T2 = T2;
    else m_T2 = 0;

    switch (m_PllAxis)
     {
     case AxisX:
        m_Xpi = XYZpi;  m_Ypi = 0.0; m_Zpi = 0.0;
        A= -T2; D= T2*XYZpi; E= 1;  /* G=-Ypi;*/    H=1;    /*I=-Zpi;*/    J= /*Ypi*Ypi + Zpi*Zpi*/- T2*XYZpi*XYZpi;
        //B=C=F=G=I=0
         break;
     case AxisY:
        m_Xpi = 0.0;  m_Ypi = XYZpi; m_Zpi = 0.0;
        A= 1;   /*D= -Xpi;*/  E= -T2; G= T2*XYZpi; H= 1;   /*I= -Zpi;*/   J= /*Xpi*Xpi +  Zpi*Zpi*/ - T2*XYZpi*XYZpi;
        //B=C=F=D=I=0
         break;
     case AxisZ:
        m_Xpi = 0.0;  m_Ypi = 0.0; m_Zpi = XYZpi;
        A= 1;   /*D= -Xpi; */ E= 1;   /*G=-Ypi;*/    H= -T2; I= T2*XYZpi; J= /*Xpi*Xpi + Ypi*Ypi*/ -T2*XYZpi*XYZpi;
        //B=C=F=D=G=0
         break;
     }
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceCone::G4HalfSpaceCone", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}

G4ThreeVector   G4HalfSpaceConeOnAxis::Normal(G4ThreeVector aPoint) const
{
    G4ThreeVector aNormal;
    switch (m_PllAxis)
     {
     case AxisX:
        //B=C=F=G=I=0
        aNormal.setX(A*aPoint.x() + /*B*aPoint.y() + C*aPoint.z() + */D ); //*2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z() + G*/);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() /*+ I*/);//H=1; *2 is ignored
         break;
     case AxisY:
        //B=C=F=D=I=0
        aNormal.setX(  aPoint.x()/* + B*aPoint.y() + C*aPoint.z() + D*/); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ E*aPoint.y() + /*F*aPoint.z() +*/ G );//*2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() /*+ I*/);//H=1; *2 is ignored
         break;
     case AxisZ:
        //B=C=F=D=G=0
        aNormal.setX(  aPoint.x() /*+ B*aPoint.y() + C*aPoint.z() + D*/); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z() + G*/);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y()  +*/ H*aPoint.z() + I);//*2 is ignored
         break;
     }
    if (m_Sense == 1) return -aNormal.unit();
    else return aNormal.unit();
    //the normal is possible to be (0,0,0)
}

EInside G4HalfSpaceConeOnAxis::Inside(const G4ThreeVector & aPoint) const
{
    G4double aResult =0.;  //bad initiation
    switch (m_PllAxis)
     {
     case AxisX:
        //B=C=F=G=I=0
        aResult =A*aPoint.x()*aPoint.x() +
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
                            2*D*aPoint.x() +
                                aPoint.y()*aPoint.y() +  //E=1
        //                    2*F*aPoint.y()*aPoint.z() +
//                            2*G*aPoint.y() +
                                aPoint.z()*aPoint.z() +  //H=1
//                            2*I*aPoint.z()   +
                              J;
        break;
     case AxisY:
        //B=C=F=D=I=0
        aResult =  aPoint.x()*aPoint.x() + //A=1
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
//                            2*D*aPoint.x() +
                              E*aPoint.y()*aPoint.y() +
        //                    2*F*aPoint.y()*aPoint.z() +
                            2*G*aPoint.y() +
                                aPoint.z()*aPoint.z() + //H=1
//                            2*I*aPoint.z()   +
                              J;
        break;
     case AxisZ:
        //B=C=F=D=G=0
        aResult =  aPoint.x()*aPoint.x() + //A=1
        //                    2*B*aPoint.x()*aPoint.y() +
        //                    2*C*aPoint.x()*aPoint.z() +
//                            2*D*aPoint.x() +
                                aPoint.y()*aPoint.y() +  //E=1
        //                    2*F*aPoint.y()*aPoint.z() +
//                            2*G*aPoint.y() +
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

G4double   G4HalfSpaceConeOnAxis:: Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        //B=C=F=G=I=0
        aResult = 2* (  A*aPoint.x()*aVector.x() +
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
                      D*aVector.x() +
                        aPoint.y()*aVector.y() + //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
//                      G*aVector.y() +
                        aPoint.z()*aVector.z()   //H=1
                      /*I*aVector.z()*/);
         break;
     case AxisY:
        //B=C=F=D=I=0
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
//                      D*aVector.x() +
                      E*aPoint.y()*aVector.y() +
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
                      G*aVector.y() +
                        aPoint.z()*aVector.z()   //H=1
                      /*I*aVector.z()*/);
         break;
     case AxisZ:
        //B=C=F=D=G=0
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
//                      D*aVector.x() +
                        aPoint.y()*aVector.y() + //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
//                      G*aVector.y() +
                      H*aPoint.z()*aVector.z() +
                      I*aVector.z());
         break;
     }
    return aResult;
}

G4double        G4HalfSpaceConeOnAxis:: Cq(const G4ThreeVector & aPoint)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        //B=C=F=G=I=0
        aResult = A*aPoint.x()*aPoint.x() +
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() +//E=1
//                2*F*aPoint.y()*aPoint.z() +
//                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
//                2*I*aPoint.z() +
                  J ;
         break;
     case AxisY:
        //B=C=F=D=I=0
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
//                2*D*aPoint.x()+
                  E*aPoint.y()*aPoint.y() +
//                2*F*aPoint.y()*aPoint.z() +
                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
//                2*I*aPoint.z() +
                  J ;
         break;
     case AxisZ:
        //B=C=F=D=G=0
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
//                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() + //E=1
//                2*F*aPoint.y()*aPoint.z() +
//                2*G*aPoint.y() +
                  H*aPoint.z()*aPoint.z() +
                2*I*aPoint.z() +
                  J ;
         break;
     }
    return aResult;

}



//copy constructor
G4HalfSpaceConeOnAxis::G4HalfSpaceConeOnAxis(const G4HalfSpaceConeOnAxis & right)
    :G4HalfSpaceCone(right)
{
}

//copy operator
G4HalfSpaceConeOnAxis& G4HalfSpaceConeOnAxis::operator=(const G4HalfSpaceConeOnAxis & right)
{
    if(&right == this) return *this;

    //base class
    G4HalfSpaceCone::operator =(right);

    return *this;
}

//reverse Sense operator
G4HalfSpaceConeOnAxis G4HalfSpaceConeOnAxis::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceConeOnAxis aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}
