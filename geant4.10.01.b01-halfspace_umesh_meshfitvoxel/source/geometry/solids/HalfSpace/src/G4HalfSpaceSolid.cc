#include "G4HalfSpaceSolid.hh"

#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "Randomize.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

G4HalfSpaceSolid::G4HalfSpaceSolid(const G4String &name):
    G4VSolid(name)
{
    m_BBox = NULL;
//    m_isActive = true;
}

/*!
 * \brief G4HalfSpaceSolid::G4HalfSpaceSolid
 *  create a G4HalfSpaceSolid using a list of surface,
 *  the sense of these surface should be correctly set
 * \param name name of the solid
 * \param aSurfaceList a vector container with pointers
 *  to the surfaces, this solid takes the ownership of these pointer!
 */
G4HalfSpaceSolid::G4HalfSpaceSolid(const G4String& name, std::vector<G4HalfSpaceSurface*> aSurfaceList,
                                   const G4ThreeVector& BBoxLowerPoint, const G4ThreeVector& BBoxHigherPoint ,
                                   const G4double & Volume , const G4double & SurfaceArea  ):
    G4VSolid(name), m_SurfaceList(aSurfaceList)
{
    m_BBox = NULL;
    setBoundaryBox(BBoxLowerPoint, BBoxHigherPoint);
    setVolume(Volume);
    setSurfaceArea(SurfaceArea);
    checkValidity(); //check the solid error

    //we use a explicit tolerance larger than the default 1e-9,
    //because we usually handle CAD geometries.
    kCarTolerance = 1e-7 ;
//    m_isActive = true;
}


G4HalfSpaceSolid::G4HalfSpaceSolid( __void__& a )
    : G4VSolid(a), m_BBox (NULL),
   fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.),
   fCubicVolume(0.), fSurfaceArea(0.), fpPolyhedron(0)
{
}


G4HalfSpaceSolid::G4HalfSpaceSolid(const G4HalfSpaceSolid& rhs)
    : G4VSolid(rhs), m_SurfaceList (rhs.m_SurfaceList),
    fStatistics(rhs.fStatistics), fCubVolEpsilon(rhs.fCubVolEpsilon),
    fAreaAccuracy(rhs.fAreaAccuracy), fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea)
{
    *m_BBox = *rhs.m_BBox;
    *fpPolyhedron = *rhs.fpPolyhedron;
}

G4HalfSpaceSolid& G4HalfSpaceSolid::operator = (const G4HalfSpaceSolid& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   if (m_BBox != NULL) delete m_BBox;
   *m_BBox = *rhs.m_BBox;
   m_SurfaceList = rhs.m_SurfaceList;
   fStatistics= rhs.fStatistics; fCubVolEpsilon= rhs.fCubVolEpsilon;
   fAreaAccuracy= rhs.fAreaAccuracy; fCubicVolume= rhs.fCubicVolume;
   fSurfaceArea= rhs.fSurfaceArea;
   *fpPolyhedron = *rhs.fpPolyhedron;

   return *this;
}


G4HalfSpaceSolid::~G4HalfSpaceSolid()
{
    //because we take the ownership, so we need to delete
    for (unsigned  int i=0; i<m_SurfaceList.size(); i++)
        delete m_SurfaceList[i];
    m_SurfaceList.clear();

    if (m_BBox != NULL) delete m_BBox;
    if (fpPolyhedron != 0) delete fpPolyhedron;

}

G4bool G4HalfSpaceSolid::checkValidity()
{
    if (!isBoundaryBox()) {
        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");
    return false;
    }
    //not implemented yet
//    m_isActive = true;
    return true;
}


/*!
 * \brief G4HalfSpaceSolid::setBoundaryBox
 *  set the Boundary box,
 *  the boundary box is not calculated in this class, therefore
 *  it is mendatory to set a boundary box for this solid!
 * \param aLowerPoint  the lower point of the bounadry box
 * \param aHigherPoint the higher point of the boundary box
 */
void    G4HalfSpaceSolid::setBoundaryBox(const G4ThreeVector& aLowerPoint, const G4ThreeVector& aHigherPoint  )
{
// no need    if (aLowerPoint.x() >= aHigherPoint.x() ||
//            aLowerPoint.y() >= aHigherPoint.y() ||
//            aLowerPoint.z() >= aHigherPoint.z() ) {
//        G4Exception("G4HalfSpaceSolid::setBoundaryBox", "GeomSolids1003",
//                    JustWarning, "The second point should be a higher point! boudnary box not set!");
//        return;
//    }
    if (m_BBox != NULL ) delete m_BBox;
    m_BBox = new G4BoundingBox3D (aLowerPoint, aHigherPoint);
}

void        G4HalfSpaceSolid::setVolume(const G4double & aVolume)
{
    if (aVolume > 0) fCubicVolume = aVolume;
    else fCubicVolume = 0.0;
}

void        G4HalfSpaceSolid::setSurfaceArea(const G4double & aArea)
{
    if (aArea > 0) fSurfaceArea = aArea;
    else fSurfaceArea = 0.0;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4HalfSpaceSolid::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
}


//from  G4BREPSolid::CalculateExtent
G4bool G4HalfSpaceSolid::CalculateExtent(const EAxis              pAxis      ,
                       const G4VoxelLimits&     pVoxelLimit,
                       const G4AffineTransform& pTransform ,
                       G4double&                pMin       ,
                       G4double&                pMax        ) const
{
    if (!isBoundaryBox())
        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");
    G4ThreeVector Min = m_BBox->GetBoxMin();
    G4ThreeVector Max = m_BBox->GetBoxMax();

    if (!pTransform.IsRotated())
      {
        // Special case handling for unrotated boxes
        // Compute x/y/z mins and maxs respecting limits, with early returns
        // if outside limits. Then switch() on pAxis
        //
        G4double xoffset,xMin,xMax;
        G4double yoffset,yMin,yMax;
        G4double zoffset,zMin,zMax;

        xoffset=pTransform.NetTranslation().x();
        xMin=xoffset+Min.x();
        xMax=xoffset+Max.x();
        if (pVoxelLimit.IsXLimited())
        {
          if (xMin>pVoxelLimit.GetMaxXExtent()
              ||xMax<pVoxelLimit.GetMinXExtent())
          {
            return false;
          }
          else
          {
            if (xMin<pVoxelLimit.GetMinXExtent())
            {
              xMin=pVoxelLimit.GetMinXExtent();
            }
            if (xMax>pVoxelLimit.GetMaxXExtent())
            {
              xMax=pVoxelLimit.GetMaxXExtent();
            }
          }
        }

        yoffset=pTransform.NetTranslation().y();
        yMin=yoffset+Min.y();
        yMax=yoffset+Max.y();
        if (pVoxelLimit.IsYLimited())
        {
          if (yMin>pVoxelLimit.GetMaxYExtent()
              ||yMax<pVoxelLimit.GetMinYExtent())
          {
            return false;
          }
          else
          {
            if (yMin<pVoxelLimit.GetMinYExtent())
            {
              yMin=pVoxelLimit.GetMinYExtent();
            }
            if (yMax>pVoxelLimit.GetMaxYExtent())
            {
              yMax=pVoxelLimit.GetMaxYExtent();
            }
          }
        }

        zoffset=pTransform.NetTranslation().z();
        zMin=zoffset+Min.z();
        zMax=zoffset+Max.z();
        if (pVoxelLimit.IsZLimited())
        {
          if (zMin>pVoxelLimit.GetMaxZExtent()
              ||zMax<pVoxelLimit.GetMinZExtent())
          {
            return false;
          }
          else
          {
            if (zMin<pVoxelLimit.GetMinZExtent())
            {
              zMin=pVoxelLimit.GetMinZExtent();
            }
            if (zMax>pVoxelLimit.GetMaxZExtent())
            {
              zMax=pVoxelLimit.GetMaxZExtent();
            }
          }
        }

        switch (pAxis)
        {
          case kXAxis:
            pMin=xMin;
            pMax=xMax;
            break;
          case kYAxis:
            pMin=yMin;
            pMax=yMax;
            break;
          case kZAxis:
            pMin=zMin;
            pMax=zMax;
            break;
          default:
            break;
        }
        pMin-=kCarTolerance;
        pMax+=kCarTolerance;

        return true;
      }
    else
      {
        // General rotated case - create and clip mesh to boundaries

        G4bool existsAfterClip=false;
        G4ThreeVectorList *vertices;

        pMin=+kInfinity;
        pMax=-kInfinity;

        // Calculate rotated vertex coordinates
        //
        vertices=CreateRotatedVertices(pTransform);
        ClipCrossSection(vertices,0,pVoxelLimit,pAxis,pMin,pMax);
        ClipCrossSection(vertices,4,pVoxelLimit,pAxis,pMin,pMax);
        ClipBetweenSections(vertices,0,pVoxelLimit,pAxis,pMin,pMax);

        if ( (pMin!=kInfinity) || (pMax!=-kInfinity) )
        {
          existsAfterClip=true;

          // Add 2*tolerance to avoid precision troubles
          //
          pMin-=kCarTolerance;
          pMax+=kCarTolerance;
        }
        else
        {
          // Check for case where completely enveloping clipping volume.
          // If point inside then we are confident that the solid completely
          // envelopes the clipping volume. Hence set min/max extents according
          // to clipping volume extents along the specified axis.
          //
          G4ThreeVector clipCentre(
                  (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
                  (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
                  (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);

          if (Inside(pTransform.Inverse().TransformPoint(clipCentre))!=kOutside)
          {
            existsAfterClip=true;
            pMin=pVoxelLimit.GetMinExtent(pAxis);
            pMax=pVoxelLimit.GetMaxExtent(pAxis);
          }
        }
        delete vertices;
        return existsAfterClip;
      }
}


G4ThreeVectorList*
G4HalfSpaceSolid::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
    if (!isBoundaryBox())
        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");
    G4ThreeVector Min = m_BBox->GetBoxMin();
    G4ThreeVector Max = m_BBox->GetBoxMax();


  G4ThreeVectorList *vertices;
  vertices=new G4ThreeVectorList();

  if (vertices)
  {
    vertices->reserve(8);
    G4ThreeVector vertex0(Min.x(),Min.y(),Min.z());
    G4ThreeVector vertex1(Max.x(),Min.y(),Min.z());
    G4ThreeVector vertex2(Max.x(),Max.y(),Min.z());
    G4ThreeVector vertex3(Min.x(),Max.y(),Min.z());
    G4ThreeVector vertex4(Min.x(),Min.y(),Max.z());
    G4ThreeVector vertex5(Max.x(),Min.y(),Max.z());
    G4ThreeVector vertex6(Max.x(),Max.y(),Max.z());
    G4ThreeVector vertex7(Min.x(),Max.y(),Max.z());

    vertices->push_back(pTransform.TransformPoint(vertex0));
    vertices->push_back(pTransform.TransformPoint(vertex1));
    vertices->push_back(pTransform.TransformPoint(vertex2));
    vertices->push_back(pTransform.TransformPoint(vertex3));
    vertices->push_back(pTransform.TransformPoint(vertex4));
    vertices->push_back(pTransform.TransformPoint(vertex5));
    vertices->push_back(pTransform.TransformPoint(vertex6));
    vertices->push_back(pTransform.TransformPoint(vertex7));
  }
  else
  {
    G4Exception("G4BREPSolid::CreateRotatedVertices()", "GeomSolids0003",
                FatalException, "Out of memory - Cannot allocate vertices!");
  }
  return vertices;
}

EInside G4HalfSpaceSolid::Inside(const G4ThreeVector& aPoint) const
{
    if (!isBoundaryBox())
        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");

    //check if outside boundary box
    if( !m_BBox->Inside(aPoint) )
      return kOutside;

//no need to     Reset(); //reset the solid

    //check for each surface, if one is outside , then return kOutside
    std::vector <EInside > aResultList ;
    for (unsigned int i=0; i<m_SurfaceList.size();i++) {
        EInside aResult = m_SurfaceList[i]->Inside(aPoint);
        if (aResult == kOutside)
            return kOutside;
        //else
        aResultList.push_back(aResult);
    }

    //check if on one of the surface, if yes then return kSurface
    for (unsigned int i=0; i<aResultList.size(); i++) {
        if (aResultList[i] == kSurface)
            return kSurface;
    }
    //else  kInside
    return kInside;
}

G4bool G4HalfSpaceSolid::Outside(const G4ThreeVector& aPoint, const G4int & skipIdx) const
{
//    if (!isBoundaryBox())
//        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
//                    FatalException, "The boundary box is not set!");

    //check if outside boundary box
    if( !m_BBox->Inside(aPoint) )
      return true;

    //check for each solid, if one is outside , then return kOutside
    for (unsigned int i=0; i<m_SurfaceList.size() ;i++) {
        if (i == skipIdx) continue;
        if (m_SurfaceList[i]->Inside(aPoint)== kOutside)
            return true;
    }
    return false;
}


/*!
 * \brief SurfaceNormal
 *  This function calculates the normal of the surface at a point on the
    surface. If the point is not on the surface the result is the normal of closest surface,
    and the result is not garantee to be correct!
    Note : the sense of the normal depends on the sense of the surface.
 * \return
 */
G4ThreeVector G4HalfSpaceSolid::SurfaceNormal(const G4ThreeVector& aPoint) const
{
    //Because the HowNear() method dose not return reliable value (quadric and torus)
    //Here we use Inside method instead
    //detail explaination check the General Design document
    /*
    //first try to find the surface within tolerance
    vector <G4double > aDistList;
    G4double aDist;
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        aDist = m_SurfaceList[i]->HowNear(aPoint);
        if (aDist >= -kCarTolerance && aDist <= kCarTolerance ) {
            return m_SurfaceList[i]->Normal(aPoint);
        }
        aDistList.push_back(aDist);
    }
    //second try to enlarge the tolerance
    G4double newTolerance = 100* kCarTolerance;
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        aDist = m_SurfaceList[i]->HowNear(aPoint);
        if (aDist >= -newTolerance && aDist <= newTolerance ) {
            return m_SurfaceList[i]->Normal(aPoint);
        }
    }
    //third try to find the surface with smallest distance
    aDist = kInfinity;
    G4int aIdx = -1;
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        //the aDistList has the same element as m_SurfaceList
        if (aDistList[i] < aDist) {
            aDist = aDistList[i];
            aIdx = i;
        }
    }
    if (aIdx != -1)
        return m_SurfaceList[aIdx]->Normal(aPoint);
    return G4ThreeVector(0.,0.,0.);
    */
    //first try with default tolerance
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        EInside aPosition =  m_SurfaceList[i]->Inside(aPoint);
        if (aPosition == kSurface) {
            return m_SurfaceList[i]->Normal(aPoint);
        }
    }
    //if failed, try to enlarge the tolerance and check again
    //2014-05-26 the tolerance is enlarged because many failed in
    //G4HalfSpaceTorus::Inside()
//    setCarTolerance(kCarTolerance * 100);
    setCarTolerance(kCarTolerance * 10000); //apprx. 1e-5
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        EInside aPosition =  m_SurfaceList[i]->Inside(aPoint);
        if (aPosition == kSurface) {
            setCarTolerance(kCarTolerance ); //recovered the tolerance
            return m_SurfaceList[i]->Normal(aPoint);
        }
    }
    //still failed? No idea now
    setCarTolerance(kCarTolerance ); //recovered the tolerance
    return G4ThreeVector(0.,0.,0.);
}

G4double  G4HalfSpaceSolid::DistanceToIn(const G4ThreeVector& aPoint) const
{
    if (!isBoundaryBox())
        G4Exception("G4HalfSpaceSolid::CalculateExtent", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");
    //just calculate the approximate distance to the boundary box
//    G4ThreeVector aBBoxCenter = m_BBox->GetMiddlePoint();
//    G4double      aBBoxSize = m_BBox->GetSize();//the radius of bounding sphere of the bbox
//    G4double aDist = (aPoint - aBBoxCenter).mag() - aBBoxSize ;
//    return aDist > kCarTolerance ? aDist: 0.0; //if less than tolerance, return 0.0

    // Returns largest perpendicular distance to the closest x/y/z sides of
    // the box, which is the most fast estimation of the shortest distance to box
    G4double SafeX, SafeY, SafeZ, safe = 0.0 ;
    G4ThreeVector boxMin = m_BBox->GetBoxMin();
    G4ThreeVector boxMax = m_BBox->GetBoxMax();
    G4ThreeVector midPoint = m_BBox->GetMiddlePoint();
    if (aPoint.x() <= midPoint.x()) SafeX = boxMin.x() - aPoint.x();
    else SafeX = aPoint.x() - boxMax.x();
    if (aPoint.y() <= midPoint.y()) SafeY = boxMin.y() - aPoint.y();
    else SafeY = aPoint.y() - boxMax.y();
    if (aPoint.z() <= midPoint.z()) SafeZ = boxMin.z() - aPoint.z();
    else SafeZ = aPoint.z() - boxMax.z();

    if (SafeX > safe) { safe = SafeX ; }
    if (SafeY > safe) { safe = SafeY ; }
    if (SafeZ > safe) { safe = SafeZ ; }
    return  safe;
}

G4double G4HalfSpaceSolid::DistanceToIn( const G4ThreeVector& aPoint,
                               const G4ThreeVector& aVector) const
{
//no need    Reset(); //reset the solid
    //for each surface, calculate the intersect points,
    //if not intersect, deactivate the surface
    //!!! COMMENT OUT beacuse error in parallel calculation!!
//    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
//        G4int aNbIntersect = m_SurfaceList[i]->Intersect(aPoint, aVector);
//        if (aNbIntersect == 0)
//            m_SurfaceList[i]->setActive(false);
//    }
    //!!!

    //for each intersected surface, get the intersect points
    //and check if real on the boundary of the solid
    //if smallest then the shortest distance, udpate it
    //ATTENTION, when on the tolerance layer, the direction should be check!
    G4double theMinDist = kInfinity;
    std::vector <G4ThreeVector > aIntersectPoints;
    std::vector <G4double > aIntersectDistances;
    G4double aDist;
    G4ThreeVector aNorm;
    G4double aTestStep = 10.0* kCarTolerance; //2014-08-15 add a test step
    G4ThreeVector aTestPoint;
    G4ThreeVector aUnitVec = aVector.unit() ; //normalize the vector
    G4HalfSpaceSurface * theEnterSurface = NULL; //add qiu 2016-03-24 fetch the enter surface to check the normal
    //first test if current point is on the surface of this solid, in case two solid stick together
    if (!Outside( aPoint, -1 ))    { //The point is on-surface or inside -1 means not skip any surface
        aDist = 0;  //using 0, so that the current aPoint is kept
        aTestPoint = aPoint + aUnitVec * (aDist + aTestStep);  //add a test step to avoid the ray just scratch the solid
        if (!Outside( aTestPoint, -1 )) {  //test the testpoint if inside/on-surface the solid
            //the aPoint is within the boundary, we check if the direction is
            //opposite with surface normal, is no then ignore
            for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
                if (m_SurfaceList[i]->Inside(aPoint)== kSurface) {
                    //we found the surface which the aPoint is located
                    aNorm = m_SurfaceList[i]->Normal(aPoint);
                    if (aNorm * aVector < 0)
                        //directly return this distance because no other nearer intersect points can be found.
                        return aDist;
                }
            }
        }
    }
    //it could happend that the above calculation is failed in the test step, so we have to calculate also the other points.
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        //!!! modified here because of parallel calculatoin
        aIntersectPoints.clear();
        aIntersectDistances.clear();
        G4int aNbIntersect = m_SurfaceList[i]->Intersect(aPoint, aVector, aIntersectPoints, aIntersectDistances);
        if (aNbIntersect >0){
            //        if (m_SurfaceList[i]->getIsActive()){
            //get the intersect points and distance
            //            aIntersectPoints = m_SurfaceList[i]->getIntersectPoints();
            //            aIntersectDistances = m_SurfaceList[i]->getIntersectDistances();
            for (unsigned int j=0; j< aIntersectPoints.size(); j++ ) {
                if (!Outside( aIntersectPoints[j], i )) {  //skip the current surface i
                    //these two array for intersects should be consistent!
                    aDist= aIntersectDistances[j];
                    if (aDist < theMinDist) {
                        //2014-08-15 add a test step to avoid the ray just scratch the solid
                        aTestPoint = aPoint + aUnitVec * (aDist + aTestStep);
                        if (!Outside( aTestPoint, i )) {  //test the testpoint if inside/on-surface the solid
                            if (aDist > kCarTolerance)
                                theMinDist = aDist;  //update the smallest distance
                            else {
                                //the aPoint is within the boundary, we check if the direction is
                                //opposite with surface normal, is no then ignore
                                aNorm = m_SurfaceList[i]->Normal(aPoint);
                                if (aNorm * aVector < 0)
                                    theMinDist = aDist;
                            }
                        }

                    }
                }
            }
        }
    }



    return theMinDist;
}


G4double G4HalfSpaceSolid::DistanceToOut(const G4ThreeVector& aPoint) const
{
    //calculate the distance to all the surface, return the smallest one
    //because the solid is convex, if outside one surface, then it is outside
    G4double theMinDist = kInfinity;
    for (unsigned int i=0; i<m_SurfaceList.size(); i++)     {
        //for some of the surface e.g. General Quadric and Torus
        //the HowNear method return a 0. value. for this case it is conservative value
        //it also means that the return value is not useful for ComputeSafety() in Navigator
        G4double aDist = m_SurfaceList[i]->HowNear(aPoint);
        if (aDist <= kCarTolerance) return 0.0;  // distance <0 means outside
        else if (aDist < theMinDist) //update the shortest distance
            theMinDist = aDist;
    }

    return theMinDist;
}


G4double G4HalfSpaceSolid::DistanceToOut(register const G4ThreeVector& aPoint,
                               register const G4ThreeVector& aVector,
                               const G4bool  calcNorm ,
                               G4bool        *validNorm   ,
                               G4ThreeVector *n             ) const
{
//no need    Reset(); //reset the solid
    //for each surface, calculate the intersect points,
    //if not intersect, deactivate the surface
    //!!! COMMENT OUT beacuse error in parallel calculation!!
//    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
//        G4int aNbIntersect = m_SurfaceList[i]->Intersect(aPoint, aVector);
//        if (aNbIntersect == 0)
//            m_SurfaceList[i]->setActive(false);
//    }
    //!!!

    //for each intersected surface, get the intersect points
    //and check if real on the boundary of the solid
    //if smallest then the shortest distance, udpate it
    //ATTENTION, when on the tolerance layer, the direction should be check!
    G4double theMinDist = kInfinity;
    std::vector <G4ThreeVector > aIntersectPoints;
    std::vector <G4double > aIntersectDistances;
    G4double aDist;
    G4ThreeVector aNorm;
    G4double aTestStep = 10.0* kCarTolerance; //2014-08-15 add a test step
    G4ThreeVector aTestPoint;
    G4ThreeVector aUnitVec = aVector.unit() ; //normalize the vector
    G4ThreeVector theExitPoint (kInfinity,kInfinity,kInfinity); //to obtain the exit point
    G4HalfSpaceSurface * theExitSurface = NULL; //add qiu 2014-08-15 fetch the exit surface to check the normal
//    G4String IntersectSurfInfo ;  //information of the surface which the particle intersected.
    for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
        aIntersectPoints.clear();
        aIntersectDistances.clear();
        G4int aNbIntersect = m_SurfaceList[i]->Intersect(aPoint, aVector, aIntersectPoints, aIntersectDistances);
        if (aNbIntersect >0){
//        if (m_SurfaceList[i]->getIsActive()){
            //get the intersect points and distance
//            aIntersectPoints = m_SurfaceList[i]->getIntersectPoints();
//            aIntersectDistances = m_SurfaceList[i]->getIntersectDistances();
            for (unsigned int j=0; j< aIntersectPoints.size(); j++ ) {
                if (!Outside( aIntersectPoints[j], i )) {  //skip the current surface i
                    //these two array for intersects should be consistent!
                    aDist= aIntersectDistances[j];
                    if (aDist < theMinDist) {
                        //2014-08-15 add a test step to avoid the ray just scratch the solid
                        //!!!! THIS IS DANGEROUS!!
                        //!!!! ERROR happen here, when a point is on surface and the ray goes inside, this
                        //!!!! will never suceed!

//ERROR                        aTestPoint = aPoint + aUnitVec * (aDist + aTestStep);
//ERROR                         if (Outside( aTestPoint, i )) {  //test the testpoint if inside/on-surface the solid
                            if (aDist > kCarTolerance) {
                                theExitPoint = aIntersectPoints[j];
                                theMinDist = aDist;  //update the smallest distance
                                theExitSurface = m_SurfaceList[i];//add qiu 2014-08-15
                            }
                            else {
                                //the aPoint is within the boundary, we check if the direction is
                                //opposite with surface normal, is no then ignore
                                aNorm = m_SurfaceList[i]->Normal(aPoint);
                                if (aNorm * aVector > 0) { //opposite to DistanceToIn
                                    theExitPoint = aIntersectPoints[j];
                                    theMinDist = aDist;
                                    theExitSurface = m_SurfaceList[i];//add qiu 2014-08-15
                                }
                            }
//ERROR                       }
                    }
                }
            }
        }
    }

    //calculate the exit normal if neccessary
    if (calcNorm && theExitSurface != NULL) {
        //add qiu 2014-08-15
        G4HSSurfType theExitSurfaceType = theExitSurface->getType();
        switch (theExitSurfaceType)
        {
        case G4HSPlane:
            //for plane, the solid should lays on its behind
            *validNorm = true; //this side is convex
            if (theExitPoint.x() != kInfinity)
                *n = SurfaceNormal(theExitPoint);
            break;
        case G4HSTorus: case G4HSQuadric:
            //for torus, we could not decide, because both side can be concave
            //for General quadric, althougt we can decide for Ellipsoid, but other cannot
            // therefore, we make it conservative.
            *validNorm = false;
            if (theExitPoint.x() != kInfinity)
                *n = SurfaceNormal(theExitPoint);
            break;
        case G4HSSphere: case G4HSCone: case G4HSCylinder:
            // for sphere, cylinder and cone, if the solid lays on the positive side,
            // then it is concave, otherwise should be convex
            *validNorm = theExitSurface->getSense() == 1 ? false : true;
            if (theExitPoint.x() != kInfinity)
                *n = SurfaceNormal(theExitPoint);
            break;
        default: //never goes here
            *validNorm = false;
            if (theExitPoint.x() != kInfinity)
                *n = SurfaceNormal(theExitPoint);
            break;
        }
        //else remain n, what so ever

//        *validNorm = true; //because this is a convex solid
//        if (theExitPoint.x() != kInfinity)
//            *n = SurfaceNormal(theExitPoint);
    }

// for checking if the particle is lost

    std::stringstream outmsg;
    outmsg<< "##INFO##:Going out from solid: "<<this->GetName()<<" !!"<<G4endl;
    if (theMinDist >= kInfinity){
        if (theExitSurface != NULL)
            G4cout << theExitSurface->printMe(); //print the exit surface.
        G4cout << outmsg.str() ;
        G4cout << "##WARNING##:A particle is lost in solid: "<<this->GetName()<<" !!"<<G4endl;
        G4cout << std::setprecision (17) << "Location: ("<<aPoint.x()<<","<<aPoint.y()<<","<<aPoint.z()<<")\t"
               << "Direction:("<<aVector.x()<<","<<aVector.y()<<","<<aVector.z()<<")"<<std::setprecision<<G4endl;
    }
    else{
//        G4cout << "##INFO##:Going out from solid: "<<this->GetName()<<" !!"<<G4endl;

    }

    return theMinDist;
}


G4String G4HalfSpaceSolid::GetEntityType() const
{
    return "G4HalfSpaceSolid";
}

G4VSolid* G4HalfSpaceSolid::Clone() const
{
    return new G4HalfSpaceSolid(*this);
}

std::ostream& G4HalfSpaceSolid::StreamInfo(std::ostream& os) const
{
    os << "-----------------------------------------------------------\n"
       << "    *** Dump for solid - " << GetName() << " ***\n"
       << "    ===================================================\n"
       << " Solid type: " << GetEntityType() << "\n"
       << " Parameters: \n"
       << "   Number of surfaces: " << m_SurfaceList.size() << "\n"
       << "-----------------------------------------------------------\n";

    return os;
}

void G4HalfSpaceSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4HalfSpaceSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}

G4Polyhedron* G4HalfSpaceSolid::CreatePolyhedron () const
{
  // Approximate implementation, just a box ...

  G4ThreeVector Min = m_BBox->GetBoxMin();
  G4ThreeVector Max = m_BBox->GetBoxMax();

  return new G4PolyhedronBox (Max.x(), Max.y(), Max.z());
}


G4double G4HalfSpaceSolid::GetCubicVolume()
{
    if (fCubicVolume == 0) {
        if (!isBoundaryBox())
            G4Exception("G4HalfSpaceSolid::GetCubicVolume", "GeomSolids1003",
                        FatalException, "The boundary box is not set!");
        fCubicVolume = G4VSolid::GetCubicVolume();
    }
    return fCubicVolume;
}

G4double G4HalfSpaceSolid::GetSurfaceArea()
{
    if (fSurfaceArea == 0) {
        if (!isBoundaryBox())
            G4Exception("G4HalfSpaceSolid::GetSurfaceArea", "GeomSolids1003",
                        FatalException, "The boundary box is not set!");
        fSurfaceArea = G4VSolid::GetSurfaceArea();
    }
    return fSurfaceArea;
}



G4ThreeVector G4HalfSpaceSolid::GetPointOnSurface() const
{
    if (!isBoundaryBox())
        G4Exception("G4HalfSpaceSolid::GetPointOnSurface", "GeomSolids1003",
                    FatalException, "The boundary box is not set!");

//no need    Reset(); // to clear cache intersect points in last calcualtion!
    //generate a random ray inside the boundary box
    G4ThreeVector BoxMin = m_BBox->GetBoxMin();
    G4ThreeVector BoxMax = m_BBox->GetBoxMax();

    G4ThreeVector aResult ;
    std::vector <G4ThreeVector > aIntersectPoints;
    std::vector <G4double > aIntersectDistances;
    while (1)
    {
        aIntersectPoints.clear();
        aIntersectDistances.clear();
        //a random point inside the bbox, and a random direction
        G4ThreeVector aPoint = G4ThreeVector(G4RandFlat::shoot(BoxMin.x(), BoxMax.x()),
                                             G4RandFlat::shoot(BoxMin.y(), BoxMax.y()),
                                             G4RandFlat::shoot(BoxMin.z(), BoxMax.z()));
        G4ThreeVector aVector = G4ThreeVector(G4RandFlat::shoot(-1.0, 1.0),
                                              G4RandFlat::shoot(-1.0, 1.0),
                                              G4RandFlat::shoot(-1.0, 1.0));
        //calculate the intersect with surface
        std::vector <G4ThreeVector > aValidPoints;

        for (unsigned int i=0; i<m_SurfaceList.size(); i++) {
            G4int aNbIntersect = m_SurfaceList[i]->Intersect(aPoint, aVector,aIntersectPoints,aIntersectDistances);
            if (aNbIntersect != 0) {
//                vector <G4ThreeVector > aIntersectPoints = m_SurfaceList[i]->getIntersectPoints();
                for (unsigned int j=0; j< aIntersectPoints.size(); j++ ) {
                    if (!Outside( aIntersectPoints[j], i )) {  //skip the current surface i
                        aValidPoints.push_back(aIntersectPoints[j]);
                    }
                }
            }
        }
        if (aValidPoints.size() == 0)  //no intersection
            continue;
        else if (aValidPoints.size() == 1) { //only one intersection
            aResult = aValidPoints[0];
            break;
        }
        else {
        //randomly choosing one point from it
        G4int aIdx = (G4int) G4RandFlat::shoot(0.0, aValidPoints.size());
        aResult = aValidPoints[aIdx];
        break;
        }
    }

    return aResult;
}


/*!
 * \brief G4HalfSpaceSolid::Reset
 *  reset all the surfaces
 */
void        G4HalfSpaceSolid::Reset() const
{
    for (unsigned int i=0; i<m_SurfaceList.size(); i++)
        m_SurfaceList[i]->Reset();
}

void        G4HalfSpaceSolid::setCarTolerance(const G4double &aCarTolerance) const
{
    for (unsigned int i=0; i<m_SurfaceList.size(); i++)
        m_SurfaceList[i]->setTolerance(aCarTolerance);
}
