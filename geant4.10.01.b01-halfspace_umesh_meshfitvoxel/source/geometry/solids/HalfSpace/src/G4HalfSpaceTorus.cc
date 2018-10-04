#include "G4HalfSpaceTorus.hh"
#include "G4PhysicalConstants.hh"
#include "G4JTPolynomialSolver.hh"

G4HalfSpaceTorus::G4HalfSpaceTorus(const G4ThreeVector &Center,
                                   const G4ThreeVector &Axis,
                                   const G4double &MaxRadius,
                                   const G4double &MinRadius,
                                   const G4int Sense ) :EQN_EPS(1e-9)
{
    m_Center = Center;
    m_Axis = Axis;
    m_R  = MaxRadius;
    m_r1 = MinRadius;
    m_r2 = MinRadius;  // the same as m_r1
    //initiate the rotation with the axis and 0 degree angle
    G4RotationMatrix aRotate = G4RotationMatrix (Axis, 0.);
    //Transformation for transfrom torus in origin (0,0,0) and axis (0,0,1)
    //to the position which locate by axis and center
    m_Transfrom = G4AffineTransform(aRotate, Center);
    m_InvertTransfrom = m_Transfrom.Invert();
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceCylinder::G4HalfSpaceCylinder", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}

G4HalfSpaceTorus::G4HalfSpaceTorus(const G4ThreeVector & Center, const G4ThreeVector & Axis,
                 const G4double & MaxRadius, const G4double & MinRadius1,
                 const G4double & MinRadius2,const G4int Sense ):EQN_EPS(1e-9)
{
    m_Center = Center;
    m_Axis = Axis;
    m_R  = MaxRadius;
    m_r1 = MinRadius1;
    m_r2 = MinRadius2;
    //initiate the rotation with the axis and 0 degree angle
    G4RotationMatrix aRotate = G4RotationMatrix (Axis, 0.);
    //Transformation for transfrom torus in origin (0,0,0) and axis (0,0,1)
    //to the position which locate by axis and center
    m_Transfrom = G4AffineTransform(aRotate, Center);
    m_InvertTransfrom = m_Transfrom.Invert();
    if (Sense == 1 || Sense == -1) m_Sense = Sense;
    else   G4Exception("G4HalfSpaceCylinder::G4HalfSpaceCylinder", "GeomSolids1003",
                             FatalException, "The sense should be 1 or -1!");
}

G4ThreeVector  G4HalfSpaceTorus:: Normal(G4ThreeVector aPoint) const
{
    //Algorithm from Graphic Gems II
    //first transform the point from global coordinate to default torus coordinate
    G4ThreeVector bPoint = m_InvertTransfrom.TransformPoint(aPoint);
    G4double aD = sqrt(bPoint.x()*bPoint.x() + bPoint.y()*bPoint.y());
    G4double aF = 2*(aD - m_R) /(aD * m_r1 * m_r1);
    G4ThreeVector aNormal;
    aNormal.setX(bPoint.x() * aF);
    aNormal.setY(bPoint.y() * aF);
    aNormal.setZ(2 * bPoint.z() /(m_r2 * m_r2));
    aNormal = m_Transfrom.TransformAxis(aNormal);

    //adjust the normal according to the sense
    if (m_Sense == 1) return -aNormal.unit();
    else return aNormal.unit();
}

EInside   G4HalfSpaceTorus::Inside(const G4ThreeVector & aPoint) const
{
    //first transform the point from global coordinate to default torus coordinate
    G4ThreeVector bPoint = m_InvertTransfrom.TransformPoint(aPoint);
    G4double   rho, a0, b0;	      // Related constants
    rho = ( m_r1 * m_r1) / (m_r2 * m_r2);
    a0  = 4. * m_R*m_R;
    b0  = m_R*m_R - m_r1*m_r1;
    //the torus equation is separated
    G4double aTmpEq = bPoint.x()* bPoint.x() + bPoint.y() * bPoint.y() +
            rho* bPoint.z()* bPoint.z() + b0;
    G4double bTmpEq = a0 * (bPoint.x()* bPoint.x() + bPoint.y() * bPoint.y());

    G4double aResult = aTmpEq*aTmpEq - bTmpEq;
    aResult = aResult * m_Sense; //adjust by sense
    //2016-03-24 multiplying 100 is because the  "aResult" is not the actual distance to the surface
    //instead it is just an estimation by substitute the point coordinate to the surface equation.
    //it encounter failure in the situation that the point is actually on the surface
    //therefore here we enlarge the tolerance to make it pass
    //2016-03-29 the tolerance is *1000 to handle some failure on this surface type
    if (aResult > (kCarTolerance*1000)) //positive means "inside half-space", actually point is outside the torus
        return kInside;
    else if (aResult < (-kCarTolerance*1000))  //negative means outside half-space
//    if (aResult > kCarTolerance) //positive means "inside half-space", actually point is outside the torus
//        return kInside;
//    else if (aResult < -kCarTolerance)  //negative means "outside half-space"
        return kOutside;
    else
        return kSurface;
}


G4double        G4HalfSpaceTorus::HowNear( const G4ThreeVector& )const
{//too expensive to calcualte, we use the conservative way
    return 0.0;
}

//take from G4BREPSolid/G4ToroidalSurface::Intersect
G4int  G4HalfSpaceTorus::Intersect(const G4ThreeVector & aPoint, const G4ThreeVector & aVec,
                                   std::vector<G4ThreeVector> &IntersectPoints, std::vector<G4double> &IntersectDistances)
{
    // Variables. Should be optimized later...
    G4ThreeVector Base = aPoint;   // Base of the intersection ray
    G4ThreeVector DCos = aVec.unit();     // Direction cosines of the ray
    G4int	     nhits =0;		      // Number of intersections
    G4double   rhits[4];		      // Intersection distances
//    G4double   hits[4] = {0.,0.,0.,0.}; // Ordered intersection distances
    G4double   rho, a0, b0;	      // Related constants
    G4double   f, l, t, g1, q, m1, u;   // Ray dependent terms
    G4double   C[5];		      // Quartic coefficients

    //we transform the point in global coordinate to
    //the coordinate of default torus : center : 0,0,0; axis: 0,0,1
    Base = m_InvertTransfrom.TransformPoint(Base); //translate + rotate
    DCos = m_InvertTransfrom.TransformAxis(DCos);  //only rotate

    //ATTENTION: The default axis of this torus is Y-Axis
    //We switch to Z-axis to make it better for transformation!
/*
    //	Compute constants related to the torus.
    //qiu!! I think G4ToroidalSurface has a bug here!!
    //rho = a^2/b^2, for circular torus, a = b, therefore rho = 1 !!
//    G4double rnorm = m_R - m_r1;
//    rho = m_r1*m_r1 / (rnorm*rnorm);
    rho = ( m_r1 * m_r1) / (m_r2 * m_r2);
    a0  = 4. * m_R*m_R;
    b0  = m_R*m_R - m_r1*m_r1;

    //	Compute ray dependent terms.
    f = 1. - DCos.y()*DCos.y();
    l = 2. * (Base.x()*DCos.x() + Base.z()*DCos.z());
    t = Base.x()*Base.x() + Base.z()*Base.z();
    g1 = f + rho * DCos.y()*DCos.y();
    q = a0 / (g1*g1);
    m1 = (l + 2.*rho*DCos.y()*Base.y()) / g1;
    u = (t +    rho*Base.y()*Base.y() + b0) / g1;

    //	Compute the coefficients of the quartic.

    C[4] = 1.0;
    C[3] = 2. * m1;
    C[2] = m1*m1 + 2.*u - q*f;
    C[1] = 2.*m1*u - q*l;
    C[0] = u*u - q*t;
*/

    //	Compute constants related to the torus.
    //qiu!! I think G4ToroidalSurface has a bug here!!
    //rho = a^2/b^2, for circular torus, a = b, therefore rho = 1 !!
//    G4double rnorm = m_R - m_r1;
//    rho = m_r1*m_r1 / (rnorm*rnorm);
    rho = ( m_r1 * m_r1) / (m_r2 * m_r2);
    a0  = 4. * m_R*m_R;
    b0  = m_R*m_R - m_r1*m_r1;

    //	Compute ray dependent terms.
    f = 1. - DCos.z()*DCos.z();
    g1 = f + rho * DCos.z()*DCos.z();
    l = 2. * (Base.x()*DCos.x() + Base.y()*DCos.y());
    t = Base.x()*Base.x() + Base.y()*Base.y();
    q = a0 / (g1*g1);
    m1 = (l + 2.*rho*DCos.z()*Base.z()) / g1;
    u = (t +    rho*Base.z()*Base.z() + b0) / g1;

    //	Compute the coefficients of the quartic.

    //!!!!!!!!!!!!ATTENTION !!!!!!!!!!!!
    //the following code is switch from Graphic Gems solver to G4JTPolynomialSolver
    //BECAUSE the  Graphic Gems solver does not give correct answer!!

//    C[4] = 1.0;
//    C[3] = 2. * m1;
//    C[2] = m1*m1 + 2.*u - q*f;
//    C[1] = 2.*m1*u - q*l;
//    C[0] = u*u - q*t;
    C[0] = 1.0;
    C[1] = 2. * m1;
    C[2] = m1*m1 + 2.*u - q*f;
    C[3] = 2.*m1*u - q*l;
    C[4] = u*u - q*t;
//    //	Use quartic root solver found in "Graphics Gems" by Jochen Schwarze.
//        nhits = SolveQuartic (C,rhits);
    //!! use G4JTPolynomialSolver
    //Original from  G4Torus::TorusRootsJT
    G4double  srd[4] /* real solution*/, si[4] /*image solution*/ ;
    G4int num;
    G4JTPolynomialSolver  torusEq;
    num = torusEq.FindRoots( C, 4, srd, si );

    for (int i = 0; i < num; i++ )     {
      if( si[i] == 0. ) {  // store real roots (image part = 0)
          rhits[nhits] = srd[i];
          nhits++;
      }
    }

//no need    //	SolveQuartic returns root pairs in reversed order.
//    m1 = rhits[0]; u = rhits[1]; rhits[0] = u; rhits[1] = m1;
//    m1 = rhits[2]; u = rhits[3]; rhits[2] = u; rhits[3] = m1;

    //push the intersect points and distance into the array
    if (nhits == 0) return 0;
    else
    {
        G4int aNbSln = 0; //number of valid solution
        for (int i=0; i<nhits; i++) {
            G4double aT = rhits[i]; //get a real solution
            if (aT >= kCarTolerance && aT < kInfinity) {
                G4ThreeVector bPoint = G4ThreeVector(Base.x() + DCos.x()*aT,
                                                     Base.y() + DCos.y()*aT,
                                                     Base.z() + DCos.z()*aT);
                //transform it again to global coordinate
                bPoint = m_Transfrom.TransformPoint(bPoint);
//                m_IntersectPoints.push_back(bPoint);
//                m_IntersectDistances.push_back(aT);
                IntersectPoints.push_back(bPoint);
                IntersectDistances.push_back(aT);
                aNbSln++;
            }
        }
        return aNbSln;
    }
}


//copy constructor
G4HalfSpaceTorus::G4HalfSpaceTorus(const G4HalfSpaceTorus & right)
    :G4HalfSpaceSurface(right), EQN_EPS(1e-9)
{

    //for this class
    m_Center = right.m_Center;
    m_Axis = right.m_Axis;
    m_R = right.m_R;
    m_r1 = right.m_r1;
    m_r2 = right.m_r2;
    m_Transfrom = right.m_Transfrom;
    m_InvertTransfrom = right.m_InvertTransfrom;
//    EQN_EPS = right.EQN_EPS;

}

//copy operator
G4HalfSpaceTorus& G4HalfSpaceTorus::operator=(const G4HalfSpaceTorus & right)
{
    if(&right == this) return *this;
    //G4HalfSpaceSurface general
    G4HalfSpaceSurface::operator =(right);

    //for this class
    m_Center = right.m_Center;
    m_Axis = right.m_Axis;
    m_R = right.m_R;
    m_r1 = right.m_r1;
    m_r2 = right.m_r2;
    m_Transfrom = right.m_Transfrom;
    m_InvertTransfrom = right.m_InvertTransfrom;
//    EQN_EPS = right.EQN_EPS;
    return *this;
}

//reverse Sense operator
G4HalfSpaceTorus G4HalfSpaceTorus::operator-()
{
//    m_Sense = - m_Sense;
//    return *this;
    //we should not change the sense of current surface
    G4HalfSpaceTorus aCopy (*this);
    aCopy.setSense(- aCopy.getSense());
    return aCopy;  //?? realy works?
}

G4int G4HalfSpaceTorus::SolveQuartic(G4double cc[], G4double ss[]  )
{
  // From Graphics Gems I by Jochen Schwartz

  G4double  coeffs[ 4 ];
  G4double  z, u, v, sub;
  G4double  A, B, C, D;
  G4double  sq_A, p, q, r;
  G4int     i, num;

    // Normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0

  A = cc[ 3 ] / cc[ 4 ];
  B = cc[ 2 ] / cc[ 4 ];
  C = cc[ 1 ] / cc[ 4 ];
  D = cc[ 0 ] / cc[ 4 ];

  //  substitute x = y - A/4 to eliminate cubic term:
  // x^4 + px^2 + qx + r = 0

  sq_A = A * A;
  p    = - 3.0/8 * sq_A + B;
  q    = 1.0/8 * sq_A * A - 1.0/2 * A * B + C;
  r    = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D;

  if (IsZero(r))
  {
    // no absolute term: y(y^3 + py + q) = 0

    coeffs[ 0 ] = q;
    coeffs[ 1 ] = p;
    coeffs[ 2 ] = 0;
    coeffs[ 3 ] = 1;

    num = SolveCubic(coeffs, ss);

    ss[ num++ ] = 0;
  }
  else
  {
    // solve the resolvent cubic ...
    coeffs[ 0 ] = 1.0/2 * r * p - 1.0/8 * q * q;
    coeffs[ 1 ] = - r;
    coeffs[ 2 ] = - 1.0/2 * p;
    coeffs[ 3 ] = 1;

    (void) SolveCubic(coeffs, ss);

    // ... and take the one real solution ...
    z = ss[ 0 ];

    // ... to Build two quadric equations
    u = z * z - r;
    v = 2 * z - p;

    if (IsZero(u))
      u = 0;
    else if (u > 0)
      u = std::sqrt(u);
    else
      return 0;

    if (IsZero(v))
      v = 0;
    else if (v > 0)
      v = std::sqrt(v);
    else
      return 0;

    coeffs[ 0 ] = z - u;
    coeffs[ 1 ] = q < 0 ? -v : v;
    coeffs[ 2 ] = 1;

    num = SolveQuadric(coeffs, ss);

    coeffs[ 0 ]= z + u;
    coeffs[ 1 ] = q < 0 ? v : -v;
    coeffs[ 2 ] = 1;

    num += SolveQuadric(coeffs, ss + num);
  }

  // resubstitute

  sub = 1.0/4 * A;

  for (i = 0; i < num; ++i)
    ss[ i ] -= sub;

  return num;
}


G4int G4HalfSpaceTorus::SolveCubic(G4double cc[], G4double ss[]  )
{
  // From Graphics Gems I bu Jochen Schwartz
  G4int     i, num;
  G4double  sub;
  G4double  A, B, C;
  G4double  sq_A, p, q;
  G4double  cb_p, D;

  // Normal form: x^3 + Ax^2 + Bx + C = 0
  A = cc[ 2 ] / cc[ 3 ];
  B = cc[ 1 ] / cc[ 3 ];
  C = cc[ 0 ] / cc[ 3 ];

  //  substitute x = y - A/3 to eliminate quadric term:
  //	x^3 +px + q = 0
  sq_A = A * A;
  p = 1.0/3 * (- 1.0/3 * sq_A + B);
  q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);

  // use Cardano's formula
  cb_p = p * p * p;
  D = q * q + cb_p;

  if (IsZero(D))
  {
    if (IsZero(q)) // one triple solution
    {
      ss[ 0 ] = 0;
      num = 1;
    }
    else // one single and one G4double solution
    {
      G4double u = std::pow(-q,1./3.);
      ss[ 0 ] = 2 * u;
      ss[ 1 ] = - u;
      num = 2;
    }
  }
  else if (D < 0) // Casus irreducibilis: three real solutions
  {
    G4double phi = 1.0/3 * std::acos(-q / std::sqrt(-cb_p));
    G4double t = 2 * std::sqrt(-p);

    ss[ 0 ] =   t * std::cos(phi);
    ss[ 1 ] = - t * std::cos(phi + CLHEP::pi / 3);
    ss[ 2 ] = - t * std::cos(phi - CLHEP::pi / 3);
    num = 3;
  }
  else // one real solution
  {
    G4double sqrt_D = std::sqrt(D);
    G4double u = std::pow(sqrt_D - q,1./3.);
    G4double v = - std::pow(sqrt_D + q,1./3.);

    ss[ 0 ] = u + v;
    num = 1;
  }

  // resubstitute
  sub = 1.0/3 * A;

  for (i = 0; i < num; ++i)
    ss[ i ] -= sub;

  return num;
}


G4int G4HalfSpaceTorus::SolveQuadric(G4double cc[], G4double ss[] )
{
  // From Graphics Gems I by Jochen Schwartz
  G4double p, q, D;

  // Normal form: x^2 + px + q = 0
  p = cc[ 1 ] / (2 * cc[ 2 ]);
  q = cc[ 0 ] / cc[ 2 ];

  D = p * p - q;

  if (IsZero(D))
  {
    ss[ 0 ] = - p;
    return 1;
  }
  else if (D < 0)
  {
    return 0;
  }
  else if (D > 0)
  {
    G4double sqrt_D = std::sqrt(D);

    ss[ 0 ] =   sqrt_D - p;
    ss[ 1 ] = - sqrt_D - p;
    return 2;
  }

  return 0;
}

G4int G4HalfSpaceTorus::IsZero(G4double x) const
{
  if((x) > -EQN_EPS && (x) < EQN_EPS)
    return 1;
  else return 0;
}
