#include "G4HalfSpaceCylinderOnAxis.hh"


G4HalfSpaceCylinderOnAxis::G4HalfSpaceCylinderOnAxis(G4HSAxis thePllAxis, const G4double & Radius,
                                                     const G4int Sense)
{
    m_PllAxis = thePllAxis;
    m_XYZ1 = 0.0;
    m_XYZ2 = 0.0;
    if (Radius > 0.0) m_Radius = Radius;
    else m_Radius = 0.0;

    switch (m_PllAxis)
     {
     case AxisX:
        E=1; H=1; J= -Radius*Radius;
         break;
     case AxisY:
        A=1; H=1; J= -Radius*Radius;
         break;
     case AxisZ:
        A=1; E=1; J= -Radius*Radius;
         break;
     }
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceCylinder::G4HalfSpaceCylinder", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}


G4ThreeVector   G4HalfSpaceCylinderOnAxis::Normal(G4ThreeVector aPoint) const
{
    G4ThreeVector aNormal;
    switch (m_PllAxis)
     {
     case AxisX:
        aNormal.setX(/*A*aPoint.x() + B*aPoint.y() + C*aPoint.z() + D*/ 0.0); //*2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z() + G*/);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() /*+ I*/);//H=1; *2 is ignored
         break;
     case AxisY:
        aNormal.setX(  aPoint.x() /*+ B*aPoint.y() + C*aPoint.z() + D*/); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() + E*aPoint.y() + F*aPoint.z() + G*/ 0.0);//*2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() +*/   aPoint.z() /*+ I*/);//H=1; *2 is ignored
         break;
     case AxisZ:
        aNormal.setX(  aPoint.x() /*+ B*aPoint.y() + C*aPoint.z() + D*/); //A=1; *2 is ignored
        aNormal.setY(/*B*aPoint.x() +*/ aPoint.y() /*+ F*aPoint.z() + G*/);//E=1; *2 is ignored
        aNormal.setZ(/*C*aPoint.x() + F*aPoint.y() + H*aPoint.z() + I*/0.0);//*2 is ignored
         break;
     }
    if (m_Sense == 1) return -aNormal.unit();
    else return aNormal.unit();
    //the normal is possible to be (0,0,0)
}

G4double        G4HalfSpaceCylinderOnAxis::HowNear( const G4ThreeVector& aPoint)const
{

    G4double aDist = 0.0;
    switch (m_PllAxis)
    {
    case AxisX:
        aDist  = sqrt(aPoint.y()*aPoint.y() + aPoint.z()*aPoint.z()) - m_Radius;
        break;
    case AxisY:
        aDist  = sqrt(aPoint.x()*aPoint.x() + aPoint.z()*aPoint.z()) - m_Radius;
        break;
    case AxisZ:
        aDist  = sqrt(aPoint.x()*aPoint.x() + aPoint.y()*aPoint.y()) - m_Radius;
        break;
    }
    //When point is inside the cylinder, the calculate distance is negative
    //But when m_Sense =1; the cylinder represents the outside half-space,
    //therefore it should multiply with Sense
    return aDist* m_Sense;
}

EInside  G4HalfSpaceCylinderOnAxis::Inside(const G4ThreeVector & aPoint) const
{
    G4double aResult =0.;
    switch (m_PllAxis)
    {
    case AxisX:
        aResult =/*A*aPoint.x()*aPoint.x() +
                            2*B*aPoint.x()*aPoint.y() +
                            2*C*aPoint.x()*aPoint.z() +
                            2*D*aPoint.x() +*/
                              aPoint.y()*aPoint.y() + //E=1
//                            2*F*aPoint.y()*aPoint.z() +
//                            2*G*aPoint.y() +
                              aPoint.z()*aPoint.z() + //H=1
//                            2*I*aPoint.z()   +
                              J;        //A=B=C=D=F=0
        break;
    case AxisY:
        aResult =aPoint.x()*aPoint.x() + //A=1
//                            2*B*aPoint.x()*aPoint.y() +
//                            2*C*aPoint.x()*aPoint.z() +
//                            2*D*aPoint.x() +
//                              E*aPoint.y()*aPoint.y() +
//                            2*F*aPoint.y()*aPoint.z() +
//                            2*G*aPoint.y() +
                              aPoint.z()*aPoint.z() + //H=1
//                            2*I*aPoint.z()   +
                              J;        //B=C=E=F=G=0
        break;
    case AxisZ:
        aResult =aPoint.x()*aPoint.x() + //A=1
//                            2*B*aPoint.x()*aPoint.y() +
//                            2*C*aPoint.x()*aPoint.z() +
//                            2*D*aPoint.x() +
                              aPoint.y()*aPoint.y() + //E=1
//                            2*F*aPoint.y()*aPoint.z() +
//                            2*G*aPoint.y() +
//                              H*aPoint.z()*aPoint.z() +
//                            2*I*aPoint.z()   +
                              J;        //B=C=F=H=I=0
        break;
    }

//no a good way    G4double aResult = HowNear(aPoint);
    aResult = aResult * m_Sense; //adjust by sense
    //2016-03-24 multiplying 100 is because the  "aResult" is not the actual distance to the surface
    //instead it is just an estimation by substitute the point coordinate to the surface equation.
    //it encounter failure in the situation that the point is actually on the surface
    //therefore here we enlarge the tolerance to make it pass
    if (aResult > (kCarTolerance*100)) //positive means inside half-space
        return kInside;
    else if (aResult < (-kCarTolerance*100))  //negative means outside half-space
        return kOutside;
    else
        return kSurface;
}

G4double   G4HalfSpaceCylinderOnAxis:: Bq(const G4ThreeVector & aPoint, const G4ThreeVector & aVector)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        aResult = 2* (/*A*aPoint.x()*aVector.x() +*/
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
//                      D*aVector.x() +
                        aPoint.y()*aVector.y() + //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
//                      G*aVector.y() +
                        aPoint.z()*aVector.z()   //H=1
                      /*I*aVector.z()*/);
         break;
     case AxisY:
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
//                      D*aVector.x() +
//                      E*aPoint.y()*aVector.y() +
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
//                      G*aVector.y() +
                        aPoint.z()*aVector.z()   //H=1
                      /*I*aVector.z()*/);
         break;
     case AxisZ:
        aResult = 2* (  aPoint.x()*aVector.x() + //A=1
//                      B*(aPoint.x()*aVector.y()+ aVector.x()*aPoint.y())+
//                      C*(aPoint.x()*aVector.z()+ aVector.x()*aPoint.z()) +
//                      D*aVector.x() +
                        aPoint.y()*aVector.y()  //E=1
//                      F*(aPoint.y()*aVector.z() + aVector.y()*aPoint.z())+
//                      G*aVector.y()
//                      H*aPoint.z()*aVector.z() +
                      /*I*aVector.z()*/);
         break;
     }
    return aResult;
}

G4double        G4HalfSpaceCylinderOnAxis:: Cq(const G4ThreeVector & aPoint)
{
    G4double aResult = 0.0;
    switch (m_PllAxis)
     {
     case AxisX:
        aResult =/* A*aPoint.x()*aPoint.x() +*/
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
//                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() +//E=1
//                2*F*aPoint.y()*aPoint.z() +
//                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
//                2*I*aPoint.z() +
                  J ;
         break;
     case AxisY:
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
//                2*D*aPoint.x()+
//                  E*aPoint.y()*aPoint.y() +
//                2*F*aPoint.y()*aPoint.z() +
//                2*G*aPoint.y() +
                    aPoint.z()*aPoint.z() +//H=1
//                2*I*aPoint.z() +
                  J ;
         break;
     case AxisZ:
        aResult =   aPoint.x()*aPoint.x() + //A=1
//                2*B*aPoint.x()*aPoint.y() +
//                2*C*aPoint.x()*aPoint.z() +
//                2*D*aPoint.x()+
                    aPoint.y()*aPoint.y() + //E=1
//                2*F*aPoint.y()*aPoint.z() +
//                2*G*aPoint.y() +
//                  H*aPoint.z()*aPoint.z() +
//                2*I*aPoint.z() +
                  J ;
         break;
     }
    return aResult;

}



//copy constructor
G4HalfSpaceCylinderOnAxis::G4HalfSpaceCylinderOnAxis(const G4HalfSpaceCylinderOnAxis & right)
    :G4HalfSpaceCylinder(right)
{
}

//copy operator
G4HalfSpaceCylinderOnAxis& G4HalfSpaceCylinderOnAxis::operator=(const G4HalfSpaceCylinderOnAxis & right)
{
    if(&right == this) return *this;

    //base class
    G4HalfSpaceCylinder::operator =(right);

    return *this;
}

//reverse Sense operator
G4HalfSpaceCylinderOnAxis G4HalfSpaceCylinderOnAxis::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceCylinderOnAxis aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}
