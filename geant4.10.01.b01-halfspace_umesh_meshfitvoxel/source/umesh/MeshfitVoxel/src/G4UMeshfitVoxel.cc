#include "G4UMeshfitVoxel.hh"
#include <math.h>
#include <algorithm>
#include <G4String.hh>

//static member initiation
G4int  G4UMeshfitVoxel::SubdivideCriteria ;
G4int  G4UMeshfitVoxel::VerboseLevel;
vector <G4bool>     G4UMeshfitVoxel::CheckList;

G4UMeshfitVoxel::G4UMeshfitVoxel()
{
    Init();
}

G4UMeshfitVoxel::G4UMeshfitVoxel(G4LogicalVolume * theMotherLogicVolume)
{
    Init();
    m_MotherLogicVol = theMotherLogicVolume;
    m_Tolerance = theMotherLogicVolume->GetSolid()->GetTolerance();
}

G4UMeshfitVoxel::G4UMeshfitVoxel(G4LogicalVolume * theMotherLogicVolume,
                G4int theVoxelLevel,
                G4UMeshfitBoundaryBox * theVoxelBoundaryBox)
{
    Init();
    m_MotherLogicVol = theMotherLogicVolume;
    m_Level = theVoxelLevel;
    m_VoxelBoundary = theVoxelBoundaryBox;
    m_Tolerance = theMotherLogicVolume->GetSolid()->GetTolerance();

}

G4UMeshfitVoxel::~G4UMeshfitVoxel()
{
    //delete boundaries
    if (m_RealBoundary != 0)
        delete m_RealBoundary;
    if (m_VoxelBoundary != 0)
        delete m_VoxelBoundary;
    if (m_RigidBoundingSphere != 0)
        delete m_RigidBoundingSphere;
    if (m_LooseBoundingSphere != 0)
        delete m_LooseBoundingSphere;
    //delete children
    if (!m_ChildrenVoxels.empty()) {
        for (vector< G4UMeshfitVoxel * >::iterator it= m_ChildrenVoxels.begin();
             it != m_ChildrenVoxels.end(); ++it) {
            if (*it != 0) delete *it;
        }
    }
    m_ChildrenVoxels.clear();

}


/*!
 * \brief setElementList
 *set the element list, with the Physcial Volume index inside
 *the array of the mother logical volume
 * \param aElementList a container which has the index of these element
 *
 */
void G4UMeshfitVoxel:: setElementList (vector<G4int> & aElementList)
{
    m_ElementList.clear();
    for (vector<G4int>::iterator it= aElementList.begin(); it != aElementList.end(); ++it )
        m_ElementList.push_back(*it);
}
void  G4UMeshfitVoxel::setVerboseLevel(G4int aVerboseLevel)
{
    VerboseLevel  = aVerboseLevel;
}

void   G4UMeshfitVoxel:: setSubdivideCriteria(G4int exceedNoOfElm)
{
    SubdivideCriteria  = exceedNoOfElm;
}

/*!
 * \brief G4UMeshfitVoxel::setCheckList
 *  set list of boolean, marks for all daugther PV in mother LV if they are
    calculated (for safetyDistance or Step), so don't do repeat calculations
    the size() of this container SHOULD identical with GetNoDaughter
 * \param aCheckList
 */
void  G4UMeshfitVoxel::setCheckList (vector <G4bool> & aCheckList)
{
    CheckList = aCheckList;
}

vector <G4bool>  G4UMeshfitVoxel::getCheckList ()
{
    return CheckList;
}


/*!
 * \brief G4UMeshfitVoxel::sortToVoxels
 * sort the elements to children voxles
 * the sorting is perform by Recursion
 * if the number of elements in this voxel is larger than \a SubdivideCriteria, then
 * it will carried out Subdivision into children voxels
 * \param isConservative     if isConservative, random point will be generated on the solid surface, and test
   them also for which voxel contains them
 * \return
 */
bool  G4UMeshfitVoxel::sortToVoxels (const G4bool isConservative)
{
    //calculate the real boundary box
    G4UMeshfitBoundaryBox * aBdrBox = calRealBoundaryBox();
    if (aBdrBox != 0) setRealBoundaryBox(aBdrBox);
    else return false;
    //expand the size using aTolerance/2 so that the node on the boundary surface can be found
    //IMPORTANT: Only expand the big mother voxel, because otherwise will inter an infinite loop!
    if (m_Level == 0) m_RealBoundary->expandSize(m_Tolerance/2);
    if (VerboseLevel >=3 ) m_RealBoundary->printMyself();

    //calculte the bounding spheres
    calBoundingSphere();
    m_RigidBoundingSphere->expandSize(m_Tolerance/2);
    m_LooseBoundingSphere->expandSize(m_Tolerance/2);


    //check if needs to create children voxel
    if (m_ElementList.size() > SubdivideCriteria)
        createChildrenVoxels();

    if (VerboseLevel >=3 )
        G4cout << "Voxel level "<<m_Level<<" has "<< m_ElementList.size() << " elements."<<G4endl;
    if (hasChildrenVoxels())
    {
        if (VerboseLevel >=3 )
            G4cout << "                  "<< m_ChildrenVoxels.size() << " voxles."<<G4endl;
        //Loop the element list and sort them into children voxels
        vector <G4int> aPartlyElementList;
        for (vector <G4int>::iterator itElm= m_ElementList.begin(); itElm!=m_ElementList.end(); ++itElm)
        {
            G4int aElmNo = *itElm;
            G4Tet * aTetElm = getElement(aElmNo);
            vector <G4ThreeVector> aPointList = aTetElm->GetVertices(); //get the node list
            vector <G4UMeshfitVoxel* >aVoxelList; //for containing voxels these points within
            //if conservative (more careful), generate 1000 random points on the solid surface and check them also
            if (isConservative){
                for (int i=0; i<1000; i++ )  {
                    aPointList.push_back(aTetElm->GetPointOnSurface());
                }
            }
            //find voxels contains these points
            for (vector <G4ThreeVector>::iterator itPoint =aPointList.begin(); itPoint != aPointList.end(); ++itPoint )
            {
                G4ThreeVector aPoint = *itPoint;
                //first check if inside the boundary box
                if (!m_RealBoundary->isInside(aPoint)) continue;
                //find the voxel
                G4UMeshfitVoxel* aVoxel = findChildVoxelwithPoint(aPoint);
                if (aVoxel == 0){
                    G4ExceptionDescription msg;
                    msg <<  "Point not found: \t" << aPoint.x() << "\t " << aPoint.y() << "\t " << aPoint.z() << G4endl;
                    G4Exception("G4UMeshfitVoxel::ComputeSafety",
                                "umesh", FatalException, msg);
                    return false;
                }
                //remove the repeat voxels: some element will fully contain in one voxel
                G4bool isRepeat = false;
                for (vector <G4UMeshfitVoxel* >::iterator itVoxel=aVoxelList.begin();
                     itVoxel !=aVoxelList.end(); ++itVoxel ){
                    if (aVoxel == *itVoxel) isRepeat = true;
                }
                if (!isRepeat)
                    aVoxelList.push_back(aVoxel);
            }
            //remove the repeat voxels: some element will fully contain in one voxel
            //        for (vector <G4UMeshfitVoxel* >::iterator itVoxel=aVoxelList.begin();
            //             itVoxel !=aVoxelList.end(); ++itVoxel ){
            //            for (vector <G4UMeshfitVoxel* >::reverse_iterator itVoxelRev= aVoxelList.rbegin();
            //                 itVoxelRev != itVoxel; ++itVoxelRev){
            //                if (*itVoxel == *itVoxelRev) //if repeat
            //                    aVoxelList.erase(itVoxelRev); //erase it
            //            }
            //        }
            //            for (unsigned int i=0; i<aVoxelList.size(); i++ ){
            //                for (unsigned int j=aVoxelList.size()-1; j != i; j--){
            //                    if (aVoxelList[i] == aVoxelList[j]) //if repeat
            //                        aVoxelList.erase(aVoxelList.begin() + j); //erase it
            //                }
            //            }

            //push the element into children voxels
            for (unsigned int i=0; i<aVoxelList.size(); i++ )
                aVoxelList[i]->pushinElement(aElmNo);
            //keep the element which is partly contain in children voxels in mother voxel
            //therefore to prevent missing in tracking
            if (aVoxelList.size() > 1) //an element is contained in more than one voxel
                aPartlyElementList.push_back(aElmNo);
        }

        //keep only partly contained elements in mother voxel
        if (VerboseLevel >= 3)
            G4cout << "                  "<<aPartlyElementList.size() << " elements are partly contained in children voxels."<<G4endl;
        setElementList(aPartlyElementList);
        //ask children voxels to sort it them self
        for (vector <G4UMeshfitVoxel* >::iterator itVoxel=m_ChildrenVoxels.begin();
             itVoxel !=m_ChildrenVoxels.end(); ++itVoxel ){
            G4UMeshfitVoxel* aChildVoxel = *itVoxel;
            if (aChildVoxel->hasElements())
                aChildVoxel->sortToVoxels();
        }
    }
    return true;
}

/*!
 * \brief G4UMeshfitVoxel::isAnEmptyVoxel
 *  to check if this voxel contains no elements at all
 * \return return \a true if it is an empty voxel
 */
G4bool G4UMeshfitVoxel::isAnEmptyVoxel()
{
    if (hasElements() || hasChildrenVoxels()) return false;
    else return true;
}


/*!
 * \brief G4UMeshfitVoxel::findChildVoxelwithPoint
 * find the children voxel which contains this point
 * BE CAREFUL the point on the surface.
 * for example, if a point with X = Xmax of the voxel boundary,
 * then it belongs to the next voxel because of the integer round approach.
 * \param aPoint a G4ThreeVector
 * \return the pointer to the child voxel
 */
G4UMeshfitVoxel* G4UMeshfitVoxel::findChildVoxelwithPoint(G4ThreeVector & aPoint)
{
    //for example, aPoint.x = 5.38, Xmin = 5.12, Multipler= 10
    // (5.38 - 2.12) * 10 = 2.6 ~ 2 , therefore we get the aIdxX  = 2
    G4int aIdxX = int ((aPoint.x() - m_RealBoundary->getXmin()) * m_MultiplierX);
    G4int aIdxY = int ((aPoint.y() - m_RealBoundary->getYmin()) * m_MultiplierY);
    G4int aIdxZ = int ((aPoint.z() - m_RealBoundary->getZmin()) * m_MultiplierZ);
    //check the index
    if (aIdxX < 0 || aIdxY < 0 || aIdxZ < 0 ||
           aIdxX >= m_IntervalX || aIdxY >= m_IntervalY || aIdxZ >= m_IntervalZ ) {
//        G4Exception("G4UMeshfitVoxel::getElement",
//                    "umesh", FatalException, "The point is outside the voxels!");
        return 0;
    }
    return getChildVoxel(aIdxX, aIdxY, aIdxZ);

}

/*!
 * \brief G4UMeshfitVoxel::LevelLocate
 *Search positioned volumes in mother at current top level of history
 *for volume containing globalPoint. Do not test the blocked volume.
 *If a containing volume is found, `stack' the new volume and return
 *true, else return false (the point lying in the mother but not any
 *of the daughters). localPoint = global point in local system on entry,
 *point in new system on exit.
 *###Addition###
 *This method is for finding a element with a point (and if necessary, a direction)
 *this element will be put in the history stack if found.
 *It will RECURSIVELY search its children voxels and its own element list
 * \param history : to get the Volume Navigation  history
 * \param globalPoint: a Point in the world coordinates
 * \param globalDirection: a Direction in world coordinates
 * \param pLocatedOnEdge: is the Point on the Edge?
 * \param localPoint: to return the samplePoint if found the volume
 * \return if found return \a true
 */
G4bool  G4UMeshfitVoxel::LevelLocate( G4NavigationHistory& history,
                                         const G4ThreeVector& globalPoint,
                                         const G4ThreeVector* globalDirection,
                                         const G4bool  pLocatedOnEdge,
                                               G4ThreeVector &localPoint )
{
    //transform the globalPoint in Mother LV coordiantes
    G4ThreeVector motherPoint = history.GetTopTransform().TransformPoint(globalPoint);
    //check if inside boundary box, if not return
    //perhaps not necessary because it is checked before calling this method.
    //    if (!m_RealBoundary->isInside(motherPoint)) {
    //        G4cout << "G4UMeshfitVoxel::LevelLocate: Point not inside boundary box." <<G4endl;
    //        return false;
    //    }

    //if has children voxels, searching the children voxels
    if (hasChildrenVoxels())
    {
        G4UMeshfitVoxel* foundChildVoxel =  findChildVoxelwithPoint(motherPoint);
        if (foundChildVoxel != 0)
            //call the child voxel to locate this point
            if (foundChildVoxel->LevelLocate(history, globalPoint,
                                             globalDirection, pLocatedOnEdge, localPoint))
                return true;
    }

    //if no child, or not found in inside the child, search the for its own list
    //similar as G4NormalNavigation::LevelLocate
    //initial the Checklist if necessary
    if (CheckList.empty())    {
        CheckList.assign(m_MotherLogicVol->GetNoDaughters(), false);
    }
    for (vector <G4int>::iterator itElm= m_ElementList.begin(); itElm!=m_ElementList.end(); ++itElm)
    {
        //only calculate those elements which are not calculated, saving time
        if (CheckList.at(*itElm) == false)
        {
            G4VPhysicalVolume * samplePhysical = m_MotherLogicVol->GetDaughter(*itElm);
            G4VSolid * sampleSolid = samplePhysical->GetLogicalVolume()->GetSolid();
            history.NewLevel(samplePhysical, kNormal, samplePhysical->GetCopyNo());
            G4ThreeVector samplePoint = history.GetTopTransform().TransformPoint(globalPoint);
            if( G4AuxiliaryNavServices::
                    CheckPointOnSurface(sampleSolid, samplePoint, globalDirection,
                                        history.GetTopTransform(), pLocatedOnEdge) ){
                //if found
//                G4cout << "Found in Element:\t\t"<<*itElm <<"\t";
                localPoint = samplePoint;
                return true;
            }
            else {
                //if not, pop-back a level
                history.BackLevel();
            }
            CheckList.at(*itElm) = true;//mark this element as conculated.
        }
        //for testing
//        else {
//            G4cout<<".";
//        }
    }
    return false; //not found
}

// for sorting the safety value
// This function returns true if the first pair is "less"
// than the second if the second element of the first pair
// is less than the second element of the second pair
// see http://stackoverflow.com/questions/18112773/sorting-a-vector-of-pairs
bool pairCompare(const std::pair<G4int, G4double>& firstElem, const std::pair<G4int, G4double>& secondElem) {
  return firstElem.second < secondElem.second;
}
//bool G4doubleCompare(const G4double& firstElem, const G4double& secondElem) {
//  return firstElem < secondElem;
//}
/*!
 * \brief G4UMeshfitVoxel::ComputeSafety
 *  compute the safety from a point to the mesh elements, and return this safety distance;
    this safety distance is calcualted isotropically
    recursive way to calculate the safety in children voxels if require accurate value
 * \param localPoint: a point in mother logical volume coordinates
 * \param isAccurate: if \a false, return a underestimate safety which calculates only the first layer voxel
 *                    if \a true, calculate the accurate safety to the mesh elements
 * \param aSafetyMethod: method to calcualte the safety:
 *  rigid bounding sphere: no need to calculate partly contained element safety;
    loose bounding sphere: need to calculate the partly contain element safety;
    boundary box: need to calculate the partly contain element safety, more expensive
 * \return the smallest safety, if error or no elements return the kInifinity
 */
G4double   G4UMeshfitVoxel::ComputeSafety(const G4ThreeVector &localPoint,
                                          const G4bool isAccurate,
                                          const G4UMeshSafetyMethod aSafetyMethod)
{
    //if has not element, return a Infinite large value
    if (isAnEmptyVoxel()) return kInfinity;

    //if has children voxels. remember that the top voxel always has children voxels
    G4double smallestSafety = kInfinity; //initial
    if (hasChildrenVoxels())
    {
        //first, only calculate the safety to the voxels
        vector <pair < G4int, G4double> > aMapVoxelToSafety;
        for (unsigned int i=0; i<m_ChildrenVoxels.size(); i++)
        {
            G4double tmpVoxelSafety = kInfinity;
            switch (aSafetyMethod)
            {
            case Safety_RigidSphere:
                if (m_ChildrenVoxels[i]->getRigidBoundingSphere() != 0)
                    tmpVoxelSafety = distanceToBoundingSphere(m_ChildrenVoxels[i]->getRigidBoundingSphere(), localPoint);
                break;
            case Safety_LooseSphere:
                if (m_ChildrenVoxels[i]->getLooseBoundingSphere() != 0)
                    tmpVoxelSafety = distanceToBoundingSphere(m_ChildrenVoxels[i]->getLooseBoundingSphere(), localPoint);
                break;
            case Safety_BoundaryBox:
                if (m_ChildrenVoxels[i]->getRealBoundaryBox() != 0)
                    tmpVoxelSafety = distanceToBoundaryBox(m_ChildrenVoxels[i]->getRealBoundaryBox(),localPoint);
                break;
            default:
                G4Exception("G4UMeshfitVoxel::ComputeSafety",
                            "umesh", FatalException, "Unknown safety calcualtion option!");
                return kInfinity;
            }
            if (tmpVoxelSafety != kInfinity)
                aMapVoxelToSafety.push_back(make_pair(i, tmpVoxelSafety)); // make a pair with voxel id and its safety value
        }
        if (aMapVoxelToSafety.empty()) {
            G4Exception("G4UMeshfitVoxel::ComputeSafety",
                        "umesh", FatalException, "Strange, no safety value is calcualte for all voxels!");
            return kInfinity;
        }
        //sort the array according to the safety value,
        sort(aMapVoxelToSafety.begin(), aMapVoxelToSafety.end(), pairCompare);

        if (!isAccurate) { //if no need to be accurate
            //record the smallest safety to the children voxels
            smallestSafety = aMapVoxelToSafety[0].second;
        }
        else
        {
            vector <G4double > aAccurateSafetyInVoxels; // vector to contain accurate safety calculated in each voxel
            unsigned int i;
            for (i=0; i<aMapVoxelToSafety.size() -1; i++) // -1 because exist of i+1
            {
                //obtain the voxel id of the sorted underestimated voxel safety array
                //and compute the accurate safety
                G4double tmpSafety = m_ChildrenVoxels[aMapVoxelToSafety[i].first]->ComputeSafety(localPoint, isAccurate, aSafetyMethod);
                aAccurateSafetyInVoxels.push_back(tmpSafety);
                //get the minimum accurate safety
                tmpSafety = *min_element(aAccurateSafetyInVoxels.begin(), aAccurateSafetyInVoxels.end()); //!!"*" because min_element return the iterator, not the value
                if (tmpSafety <= aMapVoxelToSafety[i+1].second) {
                   smallestSafety = tmpSafety;
                   break;
                }
                //else loop
            }
            //if reach the last voxel, get the smallest accurate safety of in all voxels
            if (i == aMapVoxelToSafety.size() -1)
            {
                G4double tmpSafety = m_ChildrenVoxels[aMapVoxelToSafety[i].first]->ComputeSafety(localPoint, isAccurate, aSafetyMethod);
                aAccurateSafetyInVoxels.push_back(tmpSafety);
                smallestSafety = *min_element(aAccurateSafetyInVoxels.begin(), aAccurateSafetyInVoxels.end()); //!!"*" because min_element return the iterator, not the value
            }
        } //end of isAccurate
    }// end of hasChildrenVoxel

    //in case of end voxels which has no children voxels, or using BoundaryBox or LooseSphere method to calculate underestimate safety
    // we have to calculate the safety to all elements in the element list, no matter isAccurate or not
    if (!hasChildrenVoxels() || aSafetyMethod == Safety_BoundaryBox || aSafetyMethod == Safety_LooseSphere )
    {
        //initial the Checklist if necessary
        if (CheckList.empty())    {
            CheckList.assign(m_MotherLogicVol->GetNoDaughters(), false);
        }
        for (unsigned int i=0; i<m_ElementList.size(); i++)
        {
            //only calculate those elements which are not calculated, saving time
            if (CheckList.at(m_ElementList[i]) == false)
            {
                //get the Physcial volume of the element
                G4VPhysicalVolume * samplePhysical = m_MotherLogicVol->GetDaughter(m_ElementList[i]);
                //get the transformation
                G4AffineTransform sampleTf(samplePhysical->GetRotation(), samplePhysical->GetTranslation());
                sampleTf.Invert(); //????? invert?
                //transform the point into local coordinate of the element
                const G4ThreeVector samplePoint = sampleTf.TransformPoint(localPoint);
                const G4VSolid *sampleSolid = samplePhysical->GetLogicalVolume()->GetSolid();
                const G4double sampleSafety = sampleSolid->DistanceToIn(samplePoint);
                G4double tmpSafety = smallestSafety;//for testing
                smallestSafety = min(smallestSafety, sampleSafety);
                if (tmpSafety != smallestSafety)
//                    G4cout << "New safety "<< smallestSafety << " \t in element \t"<<m_ElementList[i]<<G4endl;
                CheckList.at(m_ElementList[i]) = true; //mark this element as conculated.
//                G4cout<<"."; // for test
            }
            //for testing
//            else {
//                G4cout<<".";
//            }
        }
    }

    return smallestSafety;
}


/*!
 * \brief G4UMeshfitVoxel::ComputeStep
 *  compute the step a ray to elements, and return the step length
 *  this function is essential for ray tracking, helping the navigator
 *  this is a RECURSIVE function
 * \param localPoint
 * \param aVector
 * \return
 */
G4double G4UMeshfitVoxel::ComputeStep (const G4ThreeVector &localPoint,
                                         const G4ThreeVector & localDirection)
{
    //if has not element, return a Infinite large value
    if (isAnEmptyVoxel()) return kInfinity;
    //if has children voxels.
    G4double smallestStep = kInfinity; //initial
    if (hasChildrenVoxels())
    {
        //calculate the intersect points with voxel slice planes
        vector <G4ThreeVector>  aIntersectPointList = calIntersectPointWithSlicePlanes(localPoint, localDirection);

        //if has intersect points with the voxel
        if (!aIntersectPointList.empty())
        {
            /*
            //push the starting point into the list because it might ON the slice plane surface
            //        aIntersectPointList.push_back(localPoint);
            //calculate the middle points and rank them according to the distance to the starting point
            vector <G4ThreeVector> aMidPointList;
            vector <pair <G4int, G4double> > aMapMidpointDistance ;
            for (unsigned int i=0; i<aIntersectPointList.size(); i++)
            {
                G4ThreeVector aMidPoint = G4ThreeVector ((localPoint.x() + aIntersectPointList[i].x())/2,
                                                         (localPoint.y() + aIntersectPointList[i].y())/2,
                                                         (localPoint.z() + aIntersectPointList[i].z())/2);
                aMidPointList.push_back(aMidPoint);
                //distance to the start point
                G4double  aDistance = G4ThreeVector(localPoint.x() - aMidPoint.x(),
                                                    localPoint.y() - aMidPoint.y(),
                                                    localPoint.z() - aMidPoint.z()).mag();
                aMapMidpointDistance.push_back(make_pair(i, aDistance));
            }

            //sort the midpoint using the distnace to the starting point
            sort(aMapMidpointDistance.begin(), aMapMidpointDistance.end(), pairCompare);

            //calculate the steps for the voxels which the ray goes though
            for (unsigned int i=0; i<aMapMidpointDistance.size(); i++)
            {
                G4UMeshfitVoxel * aVoxel = findChildVoxelwithPoint (aMidPointList [aMapMidpointDistance[i].first]);
                if (aVoxel == 0) continue; //if outside all voxels, check next
                //calculate the step inside this voxel
                G4double aStep = aVoxel->ComputeStep(localPoint, localDirection);
                if (aStep != kInfinity){
                    smallestStep = aStep;  // the smallest step is the first found element along the ray
                    break; //quit the searching
                }//else check next
            }
            */
            //calculate the distance of intersect points with the starting point,
            vector <pair <G4int, G4double> > aMapIntersectPointDistance ;
            for (unsigned int i=0; i<aIntersectPointList.size(); i++) {
                G4double  aDistance = G4ThreeVector(localPoint.x() - aIntersectPointList[i].x(),
                                                    localPoint.y() - aIntersectPointList[i].y(),
                                                    localPoint.z() - aIntersectPointList[i].z()).mag();
                aMapIntersectPointDistance.push_back(make_pair(i, aDistance));
            }
            //sort the these intersect point using distance from small to big
            sort (aMapIntersectPointDistance.begin(), aMapIntersectPointDistance.end(), pairCompare);
            //calculate the mid point along the ray, for each segment
            //the first mid point is calculated by the starting point and the first closest intersect point
            vector <G4ThreeVector> aMidPointList;
            G4ThreeVector aMidPoint = G4ThreeVector ((localPoint.x() + aIntersectPointList[aMapIntersectPointDistance[0].first].x())/2,  //aMapIntersectPointDistance[n].first is a index for aIntersectPointList
                                                     (localPoint.y() + aIntersectPointList[aMapIntersectPointDistance[0].first].y())/2,
                                                     (localPoint.z() + aIntersectPointList[aMapIntersectPointDistance[0].first].z())/2);
            aMidPointList.push_back(aMidPoint);
            for (int i=1; i<aMapIntersectPointDistance.size(); i++) //starting from 1 because the first one is already calculated
            {
                aMidPoint = G4ThreeVector ((aIntersectPointList[aMapIntersectPointDistance[i-1].first].x() + aIntersectPointList[aMapIntersectPointDistance[i].first].x())/2,
                        (aIntersectPointList[aMapIntersectPointDistance[i-1].first].y() + aIntersectPointList[aMapIntersectPointDistance[i].first].y())/2,
                        (aIntersectPointList[aMapIntersectPointDistance[i-1].first].z() + aIntersectPointList[aMapIntersectPointDistance[i].first].z())/2);
                aMidPointList.push_back(aMidPoint);
            }
            //calculate the steps for the voxels which the ray goes though
            //the mid point is already sorted because the intersect points are already sorted
            for (unsigned int i=0; i<aMidPointList.size(); i++)
            {
                G4UMeshfitVoxel * aVoxel = findChildVoxelwithPoint (aMidPointList[i]);
                if (aVoxel == 0) continue; //if outside all voxels, check next
                //calculate the step inside this voxel
                G4double aStep = aVoxel->ComputeStep(localPoint, localDirection);
                if (aStep != kInfinity){
                    smallestStep = aStep;  // the smallest step is the first found element along the ray
                    break; //quit the searching
                }//else check next
            }


        }//else no intersect with any elements
    }//else check it own list of elements

    //initial the Checklist if necessary
    if (CheckList.empty())    {
        CheckList.assign(m_MotherLogicVol->GetNoDaughters(), false);
    }
    //the element which is partly contained inside the children voxels should be check
    for (unsigned int i=0; i<m_ElementList.size(); i++)
    {
        //only calculate those elements which are not calculated, saving time
        if (CheckList.at(m_ElementList[i]) == false)
        {
            //get the Physcial volume of the element
            G4VPhysicalVolume * samplePhysical = m_MotherLogicVol->GetDaughter(m_ElementList[i]);
            //get the transformation
            G4AffineTransform sampleTf(samplePhysical->GetRotation(), samplePhysical->GetTranslation());
            sampleTf.Invert(); //????? invert?
            //transform the point into local coordinate of the element
            const G4ThreeVector samplePoint = sampleTf.TransformPoint(localPoint);
            const G4ThreeVector sampleDirection = sampleTf.TransformAxis(localDirection);
            const G4VSolid *sampleSolid = samplePhysical->GetLogicalVolume()->GetSolid();
            const G4double sampleStep = sampleSolid->DistanceToIn(samplePoint, sampleDirection);
            smallestStep = min(smallestStep, sampleStep);  //get the smallest
            CheckList.at(m_ElementList[i]) = true; //mark this element as conculated.
//            if (sampleStep != kInfinity)
//                G4cout << "Hit element No\t" <<m_ElementList[i] << "  with Step\t"<<sampleStep<< G4endl;
        }
    }

    return smallestStep;
}


/*!
 * \brief G4UMeshfitVoxel::Init
 * Initiate data members
 */
void G4UMeshfitVoxel::Init()
{
    m_MotherLogicVol = 0;
    m_RealBoundary = 0;
    m_VoxelBoundary = 0;
    m_RigidBoundingSphere = 0;
    m_LooseBoundingSphere = 0;
    m_Level = 0; //defaut the top
}

//to initiate data members

/*!
 * \brief G4UMeshfitVoxel::getElement
 * get the Tet element using the Physcial Volume index,
 * in future, the G4Tet might be replace with a new template class, e.g. G4UMeshElement
 * \param thePVIndex the index for the element
 * \return the element,
 */
G4Tet * G4UMeshfitVoxel::getElement (const G4int & theElmNo)
{
    //massive called function, only neccessary to check for debugging
    if (m_MotherLogicVol == 0) {
        G4Exception("G4UMeshfitVoxel::getElement",
                    "umesh", FatalException, "No mother volume!");
        return 0;
    }
    if (theElmNo >= m_MotherLogicVol->GetNoDaughters()) {
        G4Exception("G4UMeshfitVoxel::getElement",
                    "umesh", FatalException, "Exceed array bound in mother logical volume!");
        return 0;
    }
    //get and check if supported element type
    G4VSolid * aSolid = m_MotherLogicVol->GetDaughter(theElmNo)->GetLogicalVolume()->GetSolid();
    G4String aType = G4String("G4Tet");
    if (aSolid->GetEntityType() != aType) {
        G4Exception("G4UMeshfitVoxel::getElement",
                    "umesh", FatalException, "Not a Tetrahedral element!");
        return 0;
    }

    return castType(aSolid );
}


/*!
 * \brief G4UMeshfitVoxel::cast_Type
 * cast from G4VSolid* to G4Tet*. it is a explicit cast because we know what it is.
 * \param aTetSolid a G4VSolid pointer
 * \return a G4Tet pointer
 */
G4Tet * G4UMeshfitVoxel::castType( G4VSolid * aTetSolid)
{
    return static_cast <G4Tet *> (aTetSolid);
}

/*!
 * \brief G4UMeshfitVoxel::calRealBoundaryBox
 *  calculate the boundary box which contain all the elements
 * \return Pointer to the boundary box, the boundary box pointer should be detele outside
 */
G4UMeshfitBoundaryBox *     G4UMeshfitVoxel::calRealBoundaryBox()
{
    if (!hasElements()) {
        G4Exception("G4UMeshfitVoxel::calRealBoundaryBox",
                    "umesh", FatalException, "No element in the list!");
        return 0;
    }
    //obtain the list of nodes;
    vector <G4ThreeVector > aNodeList;
    for (vector<G4int>::iterator iter = m_ElementList.begin(); iter != m_ElementList.end(); ++iter) {
        G4int aTmpNo = *iter; //get a element No
        vector <G4ThreeVector > aTmpNodeList = getElement(aTmpNo)->GetVertices(); //get nodes
        if (VerboseLevel >= 4) {
            for (unsigned int i=0; i<aTmpNodeList.size(); i++) {
                G4ThreeVector aPoint =aTmpNodeList[i];
                G4cout << "Node: " << i<< "\t X=" <<aPoint.x() << "\t Y="<< aPoint.y() << "\t Z=" << aPoint.z() <<G4endl;
            }
        }
        aNodeList.insert(aNodeList.end(), aTmpNodeList.begin(), aTmpNodeList.end()); //append the nodes
    }
    if (aNodeList.empty()) {
        G4Exception("G4UMeshfitVoxel::calRealBoundaryBox",
                    "umesh", FatalException, "Node list should not be empty!");
        return 0;
    }
    if (VerboseLevel >= 4) {
        for (unsigned int i=0; i<aNodeList.size(); i++) {
            G4ThreeVector aPoint = aNodeList[i];
            G4cout << "Node: " << i<< "\t X=" <<aPoint.x() << "\t Y="<< aPoint.y() << "\t Z=" << aPoint.z() <<G4endl;
        }
    }

    //calcualte the boundary box
    G4UMeshfitBoundaryBox * aBdrBox = new G4UMeshfitBoundaryBox(aNodeList);
    //adjust the real boundary box to voxel boundary
    if (m_VoxelBoundary != 0) {
        if (aBdrBox->getXmin() < m_VoxelBoundary->getXmin())
            aBdrBox->setXmin(m_VoxelBoundary->getXmin());
        if (aBdrBox->getYmin() < m_VoxelBoundary->getYmin())
            aBdrBox->setYmin(m_VoxelBoundary->getYmin());
        if (aBdrBox->getZmin() < m_VoxelBoundary->getZmin())
            aBdrBox->setZmin(m_VoxelBoundary->getZmin());
        if (aBdrBox->getXmax() > m_VoxelBoundary->getXmax())
            aBdrBox->setXmax(m_VoxelBoundary->getXmax());
        if (aBdrBox->getYmax() > m_VoxelBoundary->getYmax())
            aBdrBox->setYmax(m_VoxelBoundary->getYmax());
        if (aBdrBox->getZmax() > m_VoxelBoundary->getZmax())
            aBdrBox->setZmax(m_VoxelBoundary->getZmax());
    }
    return aBdrBox;
}

//calculate the real boundary box using the elements in the element list;

/*!
 * \brief G4UMeshfitVoxel::setRealBoundaryBox
 * method to set the boundary box
 * \param aBdrBox a Pointer to a boundary box
 */
void G4UMeshfitVoxel::setRealBoundaryBox(G4UMeshfitBoundaryBox * aBdrBox)
{
    //first detele the old one
    if (m_RealBoundary != 0)
        delete m_RealBoundary;
    m_RealBoundary = aBdrBox;
}
/*!
 * \brief G4UMeshfitVoxel::calBoundingSphere
 *  method to calculate the bounding sphere
 * the m_RigidBoundingSphere and m_LooseBoundingSphere will be changed
 *  call this method BEFORE sorting
 */
void  G4UMeshfitVoxel::calBoundingSphere()
{
    if (!hasElements()) {
        G4Exception("G4UMeshfitVoxel::calBoundingSphere",
                    "umesh", FatalException, "No element in the list!");
        return;
    }
    //obtain the list of nodes;
    vector <G4ThreeVector > aNodeList;
    vector <G4ThreeVector > aNodeListWithinBoundaryBox;
    for (vector<G4int>::iterator iter = m_ElementList.begin(); iter != m_ElementList.end(); ++iter)
    {
        vector <G4ThreeVector > aTmpNodeList = getElement(*iter)->GetVertices(); //get nodes
        aNodeList.insert(aNodeList.end(), aTmpNodeList.begin(), aTmpNodeList.end()); //append the nodes
        //check if nodes within the boundary box
        for (unsigned int i=0; i< aTmpNodeList.size(); i++) {
            if (m_RealBoundary->isInside(aTmpNodeList[i]))
                aNodeListWithinBoundaryBox.push_back(aTmpNodeList[i]);
        }
    }
    //calculate the Bounding sphere
    if (m_RigidBoundingSphere != 0)
        delete m_RigidBoundingSphere;
    if (m_LooseBoundingSphere != 0)
        delete m_LooseBoundingSphere;
    m_RigidBoundingSphere = new G4UMeshfitBoundingSphere (aNodeList);
    m_LooseBoundingSphere = new G4UMeshfitBoundingSphere (aNodeListWithinBoundaryBox);
}

/*!
 * \brief G4UMeshfitVoxel::distanceToBoundaryBox
 *  calculate the distane to a boundary box. if inside:0.0
 * \param aBdrBox a boundary box
 * \param aPoint a point which is to be meansure
 * \return the distance
 */
G4double G4UMeshfitVoxel::distanceToBoundaryBox(const G4UMeshfitBoundaryBox *aBdrBox, const G4ThreeVector &aPoint)
{
    //if null box
    if (aBdrBox == 0) {
        G4Exception("G4UMeshfitVoxel::distanceToBoundaryBox",
                    "umesh", FatalException, "a null boundary box!");
        return kInfinity;
    }

    //if inside
    if (aBdrBox->isInside(aPoint)) return 0.0;

    //region where closed to the 8 vertices
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() <= aBdrBox->getYmin() && aPoint.z() <= aBdrBox->getZmin()) // pos 1
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() <= aBdrBox->getYmin() && aPoint.z() <= aBdrBox->getZmin()) // pos 2
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() >= aBdrBox->getYmax() && aPoint.z() <= aBdrBox->getZmin()) // pos 3
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() >= aBdrBox->getYmax() && aPoint.z() <= aBdrBox->getZmin()) // pos 4
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() <= aBdrBox->getYmin() && aPoint.z() >= aBdrBox->getZmax()) // pos 5
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() <= aBdrBox->getYmin() && aPoint.z() >= aBdrBox->getZmax()) // pos 6
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() >= aBdrBox->getYmax() && aPoint.z() >= aBdrBox->getZmax()) // pos 7
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() >= aBdrBox->getYmax() && aPoint.z() >= aBdrBox->getZmax()) // pos 8
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();

    //region where close to the 12 edge
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() <= aBdrBox->getYmin() && aPoint.z() <= aBdrBox->getZmin()) // pos 1
        return G4ThreeVector(aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() >= aBdrBox->getYmin() &&
            aPoint.y() <= aBdrBox->getYmax() && aPoint.z() <= aBdrBox->getZmin())  // pos 2
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmax() && aPoint.z() <= aBdrBox->getZmin()) // pos 3
        return G4ThreeVector( aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() >= aBdrBox->getYmin() &&
            aPoint.y() <= aBdrBox->getYmax() && aPoint.z() <= aBdrBox->getZmin()) // pos 4
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(),  aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() <= aBdrBox->getYmin() && aPoint.z() >= aBdrBox->getZmax()) // pos 5
        return G4ThreeVector(aPoint.y()-aBdrBox->getYmin(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() >= aBdrBox->getYmin() &&
            aPoint.y() <= aBdrBox->getYmax() && aPoint.z() >= aBdrBox->getZmax())  // pos 6
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmax() && aPoint.z() >= aBdrBox->getZmax()) // pos 7
        return G4ThreeVector( aPoint.y()-aBdrBox->getYmax(), aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() >= aBdrBox->getYmin() &&
            aPoint.y() <= aBdrBox->getYmax() && aPoint.z() >= aBdrBox->getZmax()) // pos 8
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(),  aPoint.z()-aBdrBox->getZmin()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() <= aBdrBox->getYmin() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax()) // pos 9
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() <= aBdrBox->getYmin() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax()) // pos 10
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmin()).mag();
    if (aPoint.x() >= aBdrBox->getXmax() && aPoint.y() >= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax()) // pos 11
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmax(), aPoint.y()-aBdrBox->getYmax()).mag();
    if (aPoint.x() <= aBdrBox->getXmin() && aPoint.y() >= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax()) // pos 12
        return G4ThreeVector(aPoint.x()-aBdrBox->getXmin(), aPoint.y()-aBdrBox->getYmax()).mag();

    //region where close to the 6 face
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() <= aBdrBox->getYmin() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax())   // pos 1
        return aBdrBox->getYmin() - aPoint.y();
    if (aPoint.x() >= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmin() && aPoint.y() <= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax())   // pos 2
        return aPoint.x() - aBdrBox->getXmax();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax())   // pos 3
        return aPoint.y() - aBdrBox->getYmax();
    if (aPoint.x() <= aBdrBox->getXmin() &&
            aPoint.y() >= aBdrBox->getYmin() && aPoint.y() <= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmin() && aPoint.z() <= aBdrBox->getZmax())   // pos 4
        return aBdrBox->getXmin() - aPoint.x();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmin() && aPoint.y() <= aBdrBox->getYmax() &&
            aPoint.z() <= aBdrBox->getZmin() )                                      // pos 5
        return aBdrBox->getZmin() - aPoint.z();
    if (aPoint.x() >= aBdrBox->getXmin() && aPoint.x() <= aBdrBox->getXmax() &&
            aPoint.y() >= aBdrBox->getYmin() && aPoint.y() <= aBdrBox->getYmax() &&
            aPoint.z() >= aBdrBox->getZmax())                                       // pos 6
        return aPoint.z() - aBdrBox->getZmax();

    //if not found at the end: exception!
    G4Exception("G4UMeshfitVoxel::distanceToBoundaryBox",
                "umesh", FatalException, "Strange, cannot calcualte the distance!");
    return kInfinity;
}

///*!
// * \brief G4UMeshfitVoxel::distanceToRealBoundaryBox
// *  calculate the distane to \a m_RealBoundary
// * \param aPoint a point to calculat the distance
// * \return the distance, if inside or on the surface: 0.0; if no real boundary box: kInifinty
// */
//G4double  G4UMeshfitVoxel::distanceToRealBoundaryBox(const G4ThreeVector &aPoint)
//{
//    if (m_RealBoundary == 0) {
//        G4Exception("G4UMeshfitVoxel::distanceToRealBoundaryBox",
//                    "umesh", FatalException, "No real bounadry box!");
//        return kInfinity;
//    }
//    return distanceToBoundaryBox(m_RealBoundary, aPoint);
//}



/*!
 * \brief G4UMeshfitVoxel::distanceToBoundingSphere
 * calculate the distance to a bounding sphere, if inside: 0.0
 * MUCH CHEAP than the distanceToBoundaryBox.
 * \param aBdSph a bounding sphere
 * \param aPoint a point to be meansure
 * \return the distance, if inside: 0.0
 */
G4double G4UMeshfitVoxel::distanceToBoundingSphere(const G4UMeshfitBoundingSphere * aBdSph, const G4ThreeVector &aPoint)
{
    //if null
    if (aBdSph == 0) {
        G4Exception("G4UMeshfitVoxel::distanceToBoundingSphere",
                    "umesh", FatalException, "a null bounding sphere!");
        return kInfinity;
    }

    //get the modulus
    G4double aModule = aBdSph->DistanceToCeneter(aPoint);
    if (aModule <= aBdSph->getRadius()) return 0.0; //if inside or on the surface, return 0.0
    else return aModule - aBdSph->getRadius(); // else return the distance to the sphere surface
}

/*!
 * \brief G4UMeshfitVoxel::isIntersectBoundingSphere
 *  calculate if the ray(a point with a vector) intersect with the bounding sphere
 *  this function help the ComputeStep to save time calculating the intersect with boundary box
 *  if the point is inside, return true;
 *  if the point is on the surface, calulcate the consine of the ray vector and sphere normal, if negative true;
 *  if outside, and has one or two intersect point, return true.
 * \param aBdSph a bounding sphere
 * \param aPoint a point where the ray start
 * \param aVector direction of the ray
 * \return if intersect, true
 * \see http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
 */
G4bool G4UMeshfitVoxel::isIntersectBoundingSphere(const G4UMeshfitBoundingSphere * aBdSph,const  G4ThreeVector & aPoint,
                                                      const G4ThreeVector & aVector)
{
    //calculate the distance to the sphere center
    G4double aDistanceToCenter = aBdSph->DistanceToCeneter(aPoint);
    G4ThreeVector aVectorToCenter = G4ThreeVector(aBdSph->getCenterPoint().x() - aPoint.x(),
                                                  aBdSph->getCenterPoint().y() - aPoint.y(),
                                                  aBdSph->getCenterPoint().z() - aPoint.z() );
    G4double aConsine = aVectorToCenter.dot(aVector); // actually this is not exactly the cosine, but the project
    //if inside , return true
    if (aDistanceToCenter < aBdSph->getRadius())
        return true;
    //if on the surface, calculat the cosine therefore check the angle
    else if (aDistanceToCenter == aBdSph->getRadius())
    {
        if (aConsine >= 0) return true;  //if the ray goes inside or perpendicular
        else return false; //goes outside
    }
    //if outside, calculate the solution
    else
    {
        if (aConsine <= 0) return false; //if the ray not goes toward the sphere
        else // calculate the solution of a line with a sphere
        {
            //see http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection for details
            G4double atmpDot = aVector.dot(aVectorToCenter);
            G4double atmpMod2 = aVectorToCenter.mag2();
            G4double aSqure = atmpDot * atmpDot - atmpMod2 + aBdSph->getRadius() * aBdSph->getRadius();
            if (aSqure >= 0) return true; // if has one or two solutions
            else return false;  //if no solution
        }
    }
}


///*!
// * \brief G4UMeshfitVoxel::distanceToRigidBoundingSphere
// *      calculate the distance to \a m_RigidBoundingSphere
// * \param aPoint a point to calculat the distance
// * \return the distance, if inside or on the surface: 0.0; if no rigid bounding sphere: kInifinty
// */
//G4double   G4UMeshfitVoxel::distanceToRigidBoundingSphere(const G4ThreeVector &aPoint)
//{
//    if (m_RigidBoundingSphere == 0) {
//        G4Exception("G4UMeshfitVoxel::distanceToRigidBoundingSphere",
//                    "umesh", FatalException, "No rigid bounding sphere!");
//        return kInfinity;
//    }
//    return distanceToBoundingSphere(m_RigidBoundingSphere, aPoint);
//}

///*!
// * \brief G4UMeshfitVoxel::distanceToLooseBoundingSphere
// *  calculate the distance to \a m_LooseBoundingSphere
// * \param aPoint a point to calculat the distance
// * \return the distance, if inside or on the surface: 0.0; if no loose bounding sphere: kInifinty
// */
//G4double   G4UMeshfitVoxel::distanceToLooseBoundingSphere(const G4ThreeVector &aPoint)
//{
//    if (m_LooseBoundingSphere == 0) {
//        G4Exception("G4UMeshfitVoxel::distanceToLooseBoundingSphere",
//                    "umesh", FatalException, "No Loose bounding sphere!");
//        return kInfinity;
//    }
//    return distanceToBoundingSphere(m_LooseBoundingSphere, aPoint);
//}



/*!
 * \brief G4UMeshfitVoxel::createChildrenVoxels
 * moethod to create children voxels. the m_ChilrenVoxels
 * the real boundary box should be calcualted first
 */
void    G4UMeshfitVoxel::createChildrenVoxels()
{
    if (m_RealBoundary == 0) {
        G4Exception("G4UMeshfitVoxel::createChildrenVoxels",
                    "umesh", FatalException, "The real boundary is not calculated!");
        return ;
    }

    //get the range of the boundary box
    G4double aXRange = m_RealBoundary->getDeltaX();
    G4double aYRange = m_RealBoundary->getDeltaY();
    G4double aZRange = m_RealBoundary->getDeltaZ();

    //#####calculate appropriate multiplier and intervals####
    //see General Design document for clear explaination
    //this method require more optimization in the future
    G4double tmpMultiplierX = 1.0, tmpMultiplierY = 1.0, tmpMultiplierZ = 1.0;
    //we currently get make the intervals as 1.0~10.0;
    while (tmpMultiplierX*aXRange <= 1.5) //should avoid to = 1 interval
        tmpMultiplierX *= 10.0;
    while (tmpMultiplierX*aXRange > 15.0) //can be = 10.0
        tmpMultiplierX /= 10.0;
    while (tmpMultiplierY*aYRange <= 1.5)
        tmpMultiplierY *= 10.0;
    while (tmpMultiplierY*aYRange > 15.0)
        tmpMultiplierY /= 10.0;
    while (tmpMultiplierZ*aZRange <= 1.5)
        tmpMultiplierZ *= 10.0;
    while (tmpMultiplierZ*aZRange > 15.0)
        tmpMultiplierZ /= 10.0;
    //calculate the sutiable intervals
    G4int tmpIntervalX = G4int (ceil(tmpMultiplierX*aXRange));  //using ceil() for rounding to upper integer
    G4int tmpIntervalY = G4int (ceil(tmpMultiplierY*aYRange));
    G4int tmpIntervalZ = G4int (ceil(tmpMultiplierZ*aZRange));
    //create children voxles, in the order or X->Y->Z
    G4double aStepX = 1/ tmpMultiplierX; //a step in X direction
    G4double aStepY = 1/ tmpMultiplierY; //a step in Y direction
    G4double aStepZ = 1/ tmpMultiplierZ; //a step in Z direction
    for (int k=0; k<tmpIntervalZ; k++) {
        for (int j=0; j<tmpIntervalY; j++) {
            for (int i=0; i<tmpIntervalX; i++) {
                //create voxel boundaries
                G4double tmpXmin = m_RealBoundary->getXmin() + aStepX * i;  //voxel boundary Xmin
                G4double tmpXmax = m_RealBoundary->getXmin() + aStepX * (i+1);//voxel boundary Xmax
                G4double tmpYmin = m_RealBoundary->getYmin() + aStepY * j;  //voxel boundary Ymin
                G4double tmpYmax = m_RealBoundary->getYmin() + aStepY * (j+1);//voxel boundary Ymax
                G4double tmpZmin = m_RealBoundary->getZmin() + aStepZ * k;  //voxel boundary Zmin
                G4double tmpZmax = m_RealBoundary->getZmin() + aStepZ * (k+1);//voxel boundary Zmax
                G4UMeshfitBoundaryBox* aVoxelBdrBox = new G4UMeshfitBoundaryBox (tmpXmin, tmpYmin, tmpZmin,tmpXmax,tmpYmax, tmpZmax);
                G4UMeshfitVoxel * aChildrenVoxel = new G4UMeshfitVoxel(m_MotherLogicVol, m_Level +1 , aVoxelBdrBox);
                m_ChildrenVoxels.push_back(aChildrenVoxel);
            }
        }
    }
    //assign to members
    m_MultiplierX = tmpMultiplierX;
    m_MultiplierY = tmpMultiplierY;
    m_MultiplierZ = tmpMultiplierZ;
    m_IntervalX = tmpIntervalX;
    m_IntervalY = tmpIntervalY;
    m_IntervalZ = tmpIntervalZ;

    return;
}

/*!
 * \brief G4UMeshfitVoxel::getChildVoxel
 * //get a voxel from the list using the index
 * \param IdxZ index (k) in Z direction
 * \param IdxY index (j) in Y direction
 * \param IdxY index (i) in X direction
 * \return the pointer to the children voxel
 */
G4UMeshfitVoxel * G4UMeshfitVoxel::getChildVoxel(G4int IdxX, G4int IdxY, G4int IdxZ)
{
    //the checking code is not always necessary in realse version
    if (m_ChildrenVoxels.empty()) {
        G4Exception("G4UMeshfitVoxel::getChildrenVoxel",
                    "umesh", FatalException, "No children voxel!");
        return 0;
    }

    G4int idx = IdxZ * (m_IntervalY*m_IntervalX) + IdxY * m_IntervalX + IdxX;
    if (idx >= m_ChildrenVoxels.size()) {
        G4Exception("G4UMeshfitVoxel::getChildrenVoxel",
                    "umesh", FatalException, "Exceed childrenVoxel list bound!");
        return 0;
    }
    return m_ChildrenVoxels.at(idx);
}

/*!
 * \brief G4UMeshfitVoxel::calIntersectPointWithSlicePlanes
 *  calculate the intersect points with voxel slice planes,
 *  the voxel slice plane is obtain from the private memebers:
 *  m_VoxelBoundary, m_IntervalX/Y/Z, m_MultiplierX/Y/Z
 *  only those point long the direction is returned
 * \param aPoint    start point of the ray
 * \param aVector  direction of the ray
 * \return intersect points long the ray
 */
vector <G4ThreeVector>      G4UMeshfitVoxel::calIntersectPointWithSlicePlanes (const G4ThreeVector & aPoint,
                                                              const G4ThreeVector & aVector)
{
    //get the tolerance for parallel checking
    G4double aTolerance = m_MotherLogicVol->GetSolid()->GetTolerance();
    vector <G4ThreeVector> aIntersectPointList;
    //##first X planes
    //check if parallel
    G4double aSign = aVector.x() < 0 ? -1.0 : 1.0; // get the sign of the X component of the vector
    G4ThreeVector aVec = G4ThreeVector(aSign * m_RealBoundary->getDeltaX(), 0, 0); // a Vector to X direction, or opposit of X direction
    G4double aAngle = aVec.angle(aVector);
    G4double aTinyEdge = tan(aAngle) *  m_RealBoundary->getDeltaX(); // calculat how small is the two direction, using radian
    if (aTinyEdge > aTolerance /2) //only if not parallel, conduct the calculation
    {
        G4double aStepX = 1/ m_MultiplierX; //a step in X direction
        for (int i=0; i <= m_IntervalX ; i++) //use <= because number of slices = m_IntervalX +1
        {
            //method to solve the intersect equation
            // the equation for the line is x=x0 + vx * t; y=y0 + vy * t; z=z0 + vz * t;
            // now we know a x, then we get the t = (x - x0) / vx;
            // then get the y  and the z using the t;
            G4double tmpXSlice = m_RealBoundary->getXmin() + aStepX * i;
            G4double theT = ( tmpXSlice - aPoint.x() ) / aVector.x();
            G4ThreeVector theIntersectPoint = G4ThreeVector(tmpXSlice,
                                                            aPoint.y() + aVector.y() * theT,
                                                            aPoint.z() + aVector.z() * theT);
            //check if along the right direction, if yes store it
            //special situation: when aPoint = theIntersectPoint, the cosine calculated here is 0, this point is not useful
//            if (aVector.cosTheta(G4ThreeVector(aPoint.x() - theIntersectPoint.x(),
//                                               aPoint.y() - theIntersectPoint.y(),
//                                               aPoint.z() - theIntersectPoint.z())) > 0) //if cosine is 1
            //above is NOT CORRECT!!, should be theIntersectPoint - aPoint
            if (aVector.cosTheta(G4ThreeVector( theIntersectPoint.x() - aPoint.x() ,
                                                theIntersectPoint.y() - aPoint.y(),
                                                theIntersectPoint.z() - aPoint.z() ) ) > 0) //if cosine is 1
                aIntersectPointList.push_back(theIntersectPoint);
            //else throw it away
        }
    } //else not calculate the intersect points

    //##Y planes
    //check if parallel
    aSign = aVector.y() < 0 ? -1.0 : 1.0; // get the sign of the Y component of the vector
    aVec = G4ThreeVector( 0, aSign * m_RealBoundary->getDeltaY(), 0); // a Vector to Y direction, or opposit of Y direction
    aAngle = aVec.angle(aVector);
    aTinyEdge = tan(aAngle) *  m_RealBoundary->getDeltaY(); // calculat how small is the two direction, using radian
    if (aTinyEdge > aTolerance /2) //only if not parallel, conduct the calculation
    {
        G4double aStepY = 1/ m_MultiplierY; //a step in Y direction
        for (int i=0; i <= m_IntervalY ; i++) //use <= because number of slices = m_IntervalY +1
        {
            G4double tmpYSlice = m_RealBoundary->getYmin() + aStepY * i;
            G4double theT = ( tmpYSlice - aPoint.y() ) / aVector.y();
            G4ThreeVector theIntersectPoint = G4ThreeVector(aPoint.x() + aVector.x() * theT,
                                                            tmpYSlice,
                                                            aPoint.z() + aVector.z() * theT);

            //check if along the right direction, if yes store it
            //            if (aVector.cosTheta(G4ThreeVector(aPoint.x() - theIntersectPoint.x(),
            //                                               aPoint.y() - theIntersectPoint.y(),
            //                                               aPoint.z() - theIntersectPoint.z())) > 0) //if cosine is 1
            if (aVector.cosTheta(G4ThreeVector( theIntersectPoint.x() - aPoint.x() ,
                                                theIntersectPoint.y() - aPoint.y(),
                                                theIntersectPoint.z() - aPoint.z() ) ) > 0) //if cosine is 1
                aIntersectPointList.push_back(theIntersectPoint);
            //else throw it away
        }
    } //else not calculate the intersect points

    //##Z planes
    //check if parallel
    aSign = aVector.z() < 0 ? -1.0 : 1.0; // get the sign of the Z component of the vector
    aVec = G4ThreeVector( 0, 0, aSign * m_RealBoundary->getDeltaZ()); // a Vector to Z direction, or opposit of Z direction
    aAngle = aVec.angle(aVector);
    aTinyEdge = tan(aAngle) *  m_RealBoundary->getDeltaZ(); // calculat how small is the two direction, using radian
    if (aTinyEdge > aTolerance /2) //only if not parallel, conduct the calculation
    {
        G4double aStepZ = 1/ m_MultiplierZ; //a step in Z direction
        for (int i=0; i <= m_IntervalZ ; i++) //use <= because number of slices = m_IntervalZ +1
        {
            G4double tmpZSlice = m_RealBoundary->getZmin() + aStepZ * i;
            G4double theT = ( tmpZSlice - aPoint.z() ) / aVector.z();
            G4ThreeVector theIntersectPoint = G4ThreeVector(aPoint.x() + aVector.x() * theT,
                                                            aPoint.y() + aVector.y() * theT,
                                                            tmpZSlice);

            //check if along the right direction, if yes store it
            //            if (aVector.cosTheta(G4ThreeVector(aPoint.x() - theIntersectPoint.x(),
            //                                               aPoint.y() - theIntersectPoint.y(),
            //                                               aPoint.z() - theIntersectPoint.z())) > 0) //if cosine is 1
            if (aVector.cosTheta(G4ThreeVector( theIntersectPoint.x() - aPoint.x() ,
                                                theIntersectPoint.y() - aPoint.y(),
                                                theIntersectPoint.z() - aPoint.z() ) ) > 0) //if cosine is 1
                aIntersectPointList.push_back(theIntersectPoint);
            //else throw it away
        }
    } //else not calculate the intersect points

    return aIntersectPointList;
}






