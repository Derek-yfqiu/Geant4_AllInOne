//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BoundingBox3D.cc
//
// ----------------------------------------------------------------------

#include "G4BoundingBox3D.hh"
#include "geomdefs.hh"
#include "G4GeometryTolerance.hh"

const G4BoundingBox3D G4BoundingBox3D::
          space( G4ThreeVector(-kInfinity, -kInfinity, -kInfinity),
         G4ThreeVector(+kInfinity, +kInfinity, +kInfinity)  );

/////////////////////////////////////////////////////////////////////////////

G4BoundingBox3D::G4BoundingBox3D()
{
  distance = 0;
  test_result = 0;
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4BoundingBox3D::G4BoundingBox3D(const G4ThreeVector& p1, const G4ThreeVector& p2)
{
  Init(p1, p2);
}

G4BoundingBox3D::G4BoundingBox3D(const G4ThreeVector& p)
{
  Init(p);
}

G4BoundingBox3D::~G4BoundingBox3D()
{
}

G4BoundingBox3D::G4BoundingBox3D(const G4BoundingBox3D& right)
  : box_min(right.box_min), box_max(right.box_max),
    distance(right.distance), test_result(right.test_result),
    MiddlePoint(right.MiddlePoint), size (right.size), GeantBox(right.GeantBox),
    kCarTolerance(right.kCarTolerance), Dx(right.Dx), Dy(right.Dy), Dz(right.Dz)
{
}

G4BoundingBox3D& G4BoundingBox3D::operator=(const G4BoundingBox3D& right)
{
  if (&right == this) return *this;
  box_min  = right.box_min;
  box_max  = right.box_max;
  distance = right.distance;
  test_result = right.test_result;
  MiddlePoint = right.MiddlePoint;
  size  = right.size;
  GeantBox = right.GeantBox;
  kCarTolerance = right.kCarTolerance;
  
  return *this;
}

void G4BoundingBox3D::Init(const G4ThreeVector& p1, const G4ThreeVector& p2)
{
  // L. Broglia
  // Maybe temporary
  // Create a BBox bigger than the reality

  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  box_min.setX( std::min(p1.x(), p2.x()) - kCarTolerance );
  box_min.setY( std::min(p1.y(), p2.y()) - kCarTolerance );
  box_min.setZ( std::min(p1.z(), p2.z()) - kCarTolerance );
  box_max.setX( std::max(p1.x(), p2.x()) + kCarTolerance );
  box_max.setY( std::max(p1.y(), p2.y()) + kCarTolerance );
  box_max.setZ( std::max(p1.z(), p2.z()) + kCarTolerance );
  
  // Calc half spaces
  GeantBox = (box_max - box_min)*0.5;
  MiddlePoint = (box_min + box_max)*0.5;
  size = (box_min - box_max).mag()/2; //size of the bounding sphere
  Dx = box_max.x() - box_min.x();
  Dy = box_max.y() - box_min.y();
  Dz = box_max.z() - box_min.z();
  test_result = 0;
  distance = 0;
}


void G4BoundingBox3D::Init(const G4ThreeVector& p)
{
  box_min= box_max= MiddlePoint= p;
  GeantBox= G4ThreeVector(0, 0, 0);
  test_result = 0;
  distance= 0;
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}


/////////////////////////////////////////////////////////////////////////////

void G4BoundingBox3D::Extend(const G4ThreeVector& p)
{
  
  // L. Broglia
  // Maybe temporary
  // Create a BBox bigger than the reality

  if (p.x() < box_min.x()) 
    box_min.setX( p.x() - kCarTolerance );
  else if (p.x() > box_max.x())
    box_max.setX( p.x() + kCarTolerance );
 
  if (p.y() < box_min.y()) 
    box_min.setY( p.y() - kCarTolerance );
  else if (p.y() > box_max.y()) 
    box_max.setY( p.y() + kCarTolerance );

  if (p.z() < box_min.z())
    box_min.setZ( p.z() - kCarTolerance );
  else if (p.z() > box_max.z())
    box_max.setZ( p.z() + kCarTolerance );

  // L. Broglia
  // Now re-calculate GeantBox and MiddlePoint
  GeantBox    = (box_max - box_min)*0.5;
  MiddlePoint = (box_min + box_max)*0.5;
  size = (box_min - box_max).mag()/2; //size of the bounding sphere
  Dx = box_max.x() - box_min.x();
  Dy = box_max.y() - box_min.y();
  Dz = box_max.z() - box_min.z();

  
}

void  G4BoundingBox3D::Margin (const G4double & aMargin)
{
    //Y.Qiu
    box_min.set(box_min.x()-aMargin, box_min.y()-aMargin, box_min.z()-aMargin);
    box_max.set(box_max.x()+aMargin, box_max.y()+aMargin, box_max.z()+aMargin);
}


////////////////////////////////////////////////////////////////////////////


G4int G4BoundingBox3D::Test(const G4ThreeVector &aPoint, const G4ThreeVector &aVector)
{
//  const G4ThreeVector&  tmp_ray_start = aPoint;
//  const G4ThreeVector& tmp_ray_dir   = aVector;

//  G4ThreeVector  ray_start = tmp_ray_start ;
//  G4ThreeVector ray_dir   = tmp_ray_dir   ;
    G4ThreeVector  ray_start = aPoint ;
    G4ThreeVector ray_dir   = aVector   ;

  G4double rayx,rayy,rayz;
  rayx = ray_start.x();
  rayy = ray_start.y();
  rayz = ray_start.z();

  // Test if ray starting point is in the bbox or not
  if((rayx < box_min.x()) || (rayx > box_max.x()) ||
     (rayy < box_min.y()) || (rayy > box_max.y()) ||		
     (rayz < box_min.z()) || (rayz > box_max.z())   )
  {
    // Outside, check for intersection with bbox
    
    // Adapt ray_starting point to box

    const G4ThreeVector ray_start2 = G4ThreeVector( ray_start - MiddlePoint );
    distance = DistanceToIn(ray_start2, ray_dir);

    if(!distance)
      test_result = 0; // Miss
    else
      test_result = 1; // Starting point outside box & hits box
  }
  else
  {
    // Inside
    // G4cout << "\nRay starting point Inside bbox.";
    test_result = 1;
    distance = 0;
  }

  return test_result;
}

///////////////////////////////////////////////////////////////////////////////


// Does an intersection exist?
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possiblity of an intersection)


G4int G4BoundingBox3D::BoxIntersect(const G4ThreeVector&      ,
                    const G4ThreeVector&  p   ,
                    const G4ThreeVector& v    ) const
{
  G4double safx, safy, safz;
  G4double  fdx,  fdy,  fdz;
  
  fdx = GeantBox.x();    
  fdy = GeantBox.y();    
  fdz = GeantBox.z();

  safx=std::fabs(p.x())-fdx;   // minimum distance to x surface of shape
  safy=std::fabs(p.y())-fdy;
  safz=std::fabs(p.z())-fdz;
  
  // Will we Intersect?
  // If safx/y/z is >=0 the point is outside/on the box's x/y/z extent.
  // If both p.X()/y/z and v.X()/y/z repectively are both positive/negative,
  // travel is in a G4ThreeVec away from the shape.

  if ( ( (p.x()*v.x()>=0.0 ) && safx>0.0 ) || 
       ( (p.y()*v.y()>=0.0 ) && safy>0.0 ) ||
       ( (p.z()*v.z()>=0.0 ) && safz>0.0 )    )
    return 0; // No intersection  	
  else
    return 1; // Possible intersection
}

///////////////////////////////////////////////////////////////////////////////


// Distance to in
// Calculate distance to box from outside - return kBig if no intersection
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possiblity of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is Inside the exact
// shape.

//G4double G4BoundingBox::distance_to_in(const G4ThreeVector& gbox, const G4ThreeVector& p, const G4ThreeVec& v) const
G4double G4BoundingBox3D::DistanceToIn(const G4ThreeVector& p,
                       const G4ThreeVector& v) const
{
    G4double safx,  safy,  safz,  snxt = 0;  // snxt = default return value
    G4double smin, sminx, sminy, sminz;
    G4double smax, smaxx, smaxy, smaxz;
    G4double stmp;
    G4double kBig = 10e20;
    G4double fdx,fdy,fdz;
    
    fdx = GeantBox.x();        
    fdy = GeantBox.y();        
    fdz = GeantBox.z();    

    safx = std::fabs(p.x())-fdx;   // minimum distance to x surface of shape
    safy = std::fabs(p.y())-fdy;
    safz = std::fabs(p.z())-fdz;

    // Will we Intersect?
    // If safx/y/z is >=0 the point is outside/on the box's x/y/z extent.
    // If both p.X()/y/z and v.X()/y/z repectively are both positive/negative,
    // travel is in a G4ThreeVec away from the shape.

    if ( ( ( p.x()*v.x()>=0.0 ) && safx>0.0) || 
	 ( ( p.y()*v.y()>=0.0 ) && safy>0.0) || 
	 ( ( p.z()*v.z()>=0.0 ) && safz>0.0)    )
      return snxt;   	
    
    // Compute min / max distance for x/y/z travel:
    if (safx<0.0)
    {
      // Inside x extent => Calc distance until trajectory leaves extent
      sminx=0.0;
      if (v.x()) 
	smaxx = fdx/std::fabs(v.x()) - p.x()/v.x();
      else
	smaxx = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.x()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = std::fabs(v.x());
	sminx = safx/stmp;
	smaxx = (fdx+std::fabs(p.x()))/stmp;
      }
    }
    
    if (safy<0.0)
    {
      // Inside y extent => Calc distance until trajectory leaves extent
      sminy=0.0;
      if (v.y()) 
	smaxy = fdy/std::fabs(v.y()) - p.y()/v.y();
      else
	smaxy = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.y()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = std::fabs(v.y());
	sminy = safy/stmp;
	smaxy = (fdy+std::fabs(p.y()))/stmp;
      }
    }
    
    if (safz<0.0)
    {
      // Inside z extent => Calc distance until trajectory leaves extent
      sminz=0.0;
      if (v.z()) 
	smaxz = fdz/std::fabs(v.z()) - p.z()/v.z();
      else 
	smaxz = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.z()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = std::fabs(v.z());
	sminz = safz/stmp;
	smaxz = (fdz+std::fabs(p.z()))/stmp;
      }
    }

    // Find minimum allowed Dist given min/max pairs
    if (sminx>sminy) 
      smin = sminx; // MAX(sminx,sminy,sminz)
    else 
      smin = sminy;
    
    if (sminz>smin) 
      smin=sminz;

    if (smaxx<smaxy) 
      smax = smaxx; // MIN(smaxx,smaxy,smaxz)
    else 
      smax = smaxy;
    
    if (smaxz<smax) 
      smax = smaxz;

    // If smin <= kCarTolerance then only clipping `tolerant' Area
    // -> no intersection
    
    if ((smin>0.) && (smin<=smax))  { snxt=smin; }
    
    return snxt;
}


///////////////////////////////////////////////////////////////////////////////

G4int G4BoundingBox3D::Inside(const G4ThreeVector& Pt) const
{
  if( ( Pt.x() >= box_min.x() && Pt.x() <= box_max.x() ) &&
      ( Pt.y() >= box_min.y() && Pt.y() <= box_max.y() ) &&
      ( Pt.z() >= box_min.z() && Pt.z() <= box_max.z() )    )
    return 1;
  else
    return 0;
}
