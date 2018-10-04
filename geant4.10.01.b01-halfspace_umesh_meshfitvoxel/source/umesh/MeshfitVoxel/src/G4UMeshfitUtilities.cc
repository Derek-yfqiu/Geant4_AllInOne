#include "G4UMeshfitUtilities.hh"

G4UMeshfitBoundaryBox::G4UMeshfitBoundaryBox
(G4double Xmin, G4double Ymin,  G4double Zmin, G4double Xmax, G4double Ymax, G4double Zmax)
{
    m_Xmin = Xmin;
    m_Ymin = Ymin;
    m_Zmin = Zmin;
    m_Xmax = Xmax;
    m_Ymax = Ymax;
    m_Zmax = Zmax;
}

G4UMeshfitBoundaryBox::G4UMeshfitBoundaryBox
(G4ThreeVector lowerPoint, G4ThreeVector upperPoint )
{
    m_Xmin = lowerPoint.getX();
    m_Ymin = lowerPoint.getY();
    m_Zmin = lowerPoint.getZ();
    m_Xmax = upperPoint.getX();
    m_Ymax = upperPoint.getY();
    m_Zmax = upperPoint.getZ();
}

//generate the boundary box by a list of nodes
G4UMeshfitBoundaryBox::G4UMeshfitBoundaryBox
(std::vector <G4ThreeVector> &aNodeList )
{
    G4double tmpXmin = 1e99,  tmpYmin = 1e99, tmpZmin = 1e99;
    G4double tmpXmax = -1e99,  tmpYmax = -1e99, tmpZmax = -1e99;

//    for (std::vector <G4ThreeVector>::iterator it= aNodeList.begin();
//         it != aNodeList.end(); ++it) {
    //        G4ThreeVector aPoint = *it;
    for (unsigned int i=0; i< aNodeList.size(); i++) {
        G4ThreeVector aPoint =  aNodeList[i];
        if (aPoint.x() < tmpXmin) tmpXmin = aPoint.x();
        if (aPoint.y() < tmpYmin) tmpYmin = aPoint.y();
        if (aPoint.z() < tmpZmin) tmpZmin = aPoint.z();
        if (aPoint.x() > tmpXmax) tmpXmax = aPoint.x();
        if (aPoint.y() > tmpYmax) tmpYmax = aPoint.y();
        if (aPoint.z() > tmpZmax) tmpZmax = aPoint.z();
    }
    m_Xmin = tmpXmin;
    m_Ymin = tmpYmin;
    m_Zmin = tmpZmin;
    m_Xmax = tmpXmax;
    m_Ymax = tmpYmax;
    m_Zmax = tmpZmax;
}

/*!
 * \brief isInside
 *  return if the point is inside the boundary box
    if a point is lying on the slice surface Xmax, Ymax, Zmax
    it is consider NOT inside because normaly this point will be
    search in next voxel
 * \param aPoint
 * \return
 */
G4bool  G4UMeshfitBoundaryBox::isInside (const G4ThreeVector & aPoint ) const
{
    //using >= for the Min and < for the Max
    if (aPoint.x() >= m_Xmin && aPoint.x() < m_Xmax &&
            aPoint.y() >= m_Ymin && aPoint.y() < m_Ymax &&
            aPoint.z() >= m_Zmin && aPoint.z() < m_Zmax)
        return true;
    else return false;
}

/*!
 * \brief G4UMeshfitBoundaryBox::expandSize
 *  expand the boundary box, the min=min-aMargin and max= max+aMargin
    it need to expand the boundary box then adding a tolerance
 * \param aMargin the dimension to expand, can be negative
 */
void  G4UMeshfitBoundaryBox::expandSize(const G4double  aMargin)
{
    m_Xmin -= aMargin;
    m_Ymin -= aMargin;
    m_Zmin -= aMargin;
    m_Xmax += aMargin;
    m_Ymax += aMargin;
    m_Zmax += aMargin;
}

void G4UMeshfitBoundaryBox::printMyself ()
{
    G4cout<< "Boundray Box: X="<< m_Xmin << "\t to \t "<< m_Xmax <<G4endl
          << "              Y="<< m_Ymin << "\t to \t "<< m_Ymax <<G4endl
          << "              Z="<< m_Zmin << "\t to \t "<< m_Zmax <<G4endl;
}



G4UMeshfitBoundingSphere::G4UMeshfitBoundingSphere(G4double Radius, G4ThreeVector & Center)
{
    m_Radius = Radius;
    m_Center = Center;
}
G4UMeshfitBoundingSphere::G4UMeshfitBoundingSphere(std::vector <G4ThreeVector> &aPointList)
{
    m_Radius = 0.0;
    m_Center = G4ThreeVector(); //(0,0,0)
    //calculate the center point by averaging all the points
    if (aPointList.empty()) return;
    G4double tmpX = 0.0, tmpY = 0.0, tmpZ = 0.0;
    for (unsigned int i=0; i<aPointList.size(); i++) {
        tmpX += aPointList[i].x();
        tmpY += aPointList[i].y();
        tmpZ += aPointList[i].z();
    }
    tmpX /= aPointList.size();
    tmpY /= aPointList.size();
    tmpZ /= aPointList.size();
    m_Center.set(tmpX, tmpY, tmpZ);

    //calculte the largest distance to all points, therefore the radius
    G4double tmpRadius = 0.0;
    G4double aTmpMod = 0.0;
    for (unsigned int i=0; i<aPointList.size(); i++) {
        G4ThreeVector aTmpVector = G4ThreeVector(aPointList[i].x() - m_Center.x(),
                                                 aPointList[i].y() - m_Center.y(),
                                                 aPointList[i].z() - m_Center.z());
        aTmpMod = aTmpVector.mag();  //get the modulus of this vector
        tmpRadius = tmpRadius >  aTmpMod ? tmpRadius : aTmpMod; //get the large one
    }
    m_Radius = tmpRadius;
}


G4double  G4UMeshfitBoundingSphere::DistanceToCeneter (const G4ThreeVector & aPoint) const
{
    return G4ThreeVector(aPoint.x()-m_Center.x(),
                         aPoint.y()-m_Center.y(),
                         aPoint.z()-m_Center.z()).mag();
}

void G4UMeshfitBoundingSphere:: expandSize(const G4double aMargin)
{
    m_Radius += aMargin;
}


