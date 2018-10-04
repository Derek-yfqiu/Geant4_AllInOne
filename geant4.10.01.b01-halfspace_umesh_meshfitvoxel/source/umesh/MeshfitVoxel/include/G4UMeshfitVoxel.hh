#ifndef G4UMESHFITVOXEL_HH
#define G4UMESHFITVOXEL_HH

#include "G4UMeshfitUtilities.hh"

#include <vector>

#include <G4Tet.hh>
#include <globals.hh>
#include <G4ThreeVector.hh>
#include <G4Tet.hh>
#include <G4LogicalVolume.hh>
#include <G4VSolid.hh>
#include <G4NavigationHistory.hh>
#include <G4AuxiliaryNavServices.hh>    // Needed for inline methods



using namespace std;

//class G4UMeshfitBoundaryBox;
//class G4VSolid;
//class G4Tet;
//class G4LogicalVolume;

//method to calculate the safety:
// rigid bounding sphere: no need to calculate partly contained element safety;
// loose bounding sphere: need to calculate the partly contain element safety;
// boundary box: need to calculate the partly contain element safety, more expensive
enum G4UMeshSafetyMethod {
    Safety_RigidSphere,
    Safety_LooseSphere,
    Safety_BoundaryBox
};

class G4UMeshfitVoxel
{
public:
    G4UMeshfitVoxel();
    G4UMeshfitVoxel(G4LogicalVolume * theMotherLogicVolume);
    G4UMeshfitVoxel(G4LogicalVolume * theMotherLogicVolume,
                    G4int theVoxelLevel,
                    G4UMeshfitBoundaryBox * theVoxelBoundaryBox);
    ~G4UMeshfitVoxel();

    void                        setVerboseLevel(G4int aVerboseLevel);

    void                        setElementList (vector<G4int> & aElementList);
    //set the element list, with the Physcial Volume index inside
    //the array of the mother logical volume

    void                        setSubdivideCriteria(G4int exceedNoOfElm);
    //set a criteria, if number of elements exceed this criteria,
    //children voxels should be create
    //SHOULD BE SET for the top level

    void                        setCheckList (vector <G4bool> & aCheckList);
    vector <G4bool>             getCheckList ();

    bool                        sortToVoxels (const G4bool isConservative = false);
    //IMPORTANT METHOD
    //sort the elements to children voxles
    //if isConservative, random point will be generated on the solid surface, and test
    //them also for which voxel contains them

    inline void                 pushinElement(G4int & theElmNo)
                                {m_ElementList.push_back(theElmNo);};
    //push back an element into the element list

    G4bool                      isAnEmptyVoxel();
    //return true if this voxel contains no elements at all

    G4UMeshfitVoxel*            findChildVoxelwithPoint(G4ThreeVector & aPoint);
    //find the children voxel which contains this point
    //BE CAREFULL the point on the surface. for example, if X = Xmax of the voxel boundary,
    //then it belongs to the next voxel because of the integer round approach.

    G4bool                      LevelLocate( G4NavigationHistory& history,
                                             const G4ThreeVector& globalPoint,
                                             const G4ThreeVector* globalDirection,
                                             const G4bool  pLocatedOnEdge,
                                                   G4ThreeVector &localPoint );
    // Search positioned volumes in mother at current top level of history
    // for volume containing globalPoint. Do not test the blocked volume.
    // If a containing volume is found, `stack' the new volume and return
    // true, else return false (the point lying in the mother but not any
    // of the daughters). localPoint = global point in local system on entry,
    // point in new system on exit.
    //### Addition###
    // This method is for finding a element with a point (and if necessary, a direction)
    // this element will be put in the history stack if found.
    // It will RECURSIVELY search its children voxels and its own element list

    G4double                    ComputeSafety(const G4ThreeVector &localPoint,
                                              const G4bool isAccurate = false,
                                              const G4UMeshSafetyMethod aSafetyMethod = Safety_RigidSphere );
    //compute the safety from a point to the mesh elements, and return this safety distance;
    //this safety distance is calcualted isotropically
    //this is a RECURSIVELY function

    G4double                    ComputeStep (const G4ThreeVector &localPoint,
                                             const G4ThreeVector & localDirection);
    //compute the step a ray to elements, and return the step length
    //this is a RECURSIVELY function




private:
    void                        Init();
    //to initiate data members

    inline G4Tet *              getElement (const G4int & theElmNo);
    //get the Tet element using the Physcial Volume index in mother logical volume
    //in future, the G4Tet might be replace with a new template class, e.g. G4UMeshElement
    inline G4Tet *              castType(G4VSolid *aTetSolid);
    //cast from G4VSolid* to G4Tet*. it is a explicit cast because we know what it is

    G4UMeshfitBoundaryBox *     calRealBoundaryBox();
    //calculate the real boundary box using the elements in the element list;

    void                        setRealBoundaryBox(G4UMeshfitBoundaryBox * aBdrBox);
    //method to set the boundary box

    inline const G4UMeshfitBoundaryBox *            getRealBoundaryBox()
                                {return m_RealBoundary;};
    //get the m_RealBoundary
    inline const G4UMeshfitBoundingSphere *        getRigidBoundingSphere()
                                {return m_RigidBoundingSphere;};
    //get the m_RigidBoundingSphere
    inline const G4UMeshfitBoundingSphere *        getLooseBoundingSphere()
                                {return m_LooseBoundingSphere;};
    //get the m_LooseBoundingSphere

    void                        calBoundingSphere();
    //method to calculate the bounding sphere
    //call this method BEFORE sorting

    G4double                    distanceToBoundaryBox(const G4UMeshfitBoundaryBox * aBdrBox, const G4ThreeVector & aPoint);
    //calculate the isotropic distance to a boundary box. if inside:0.0

    G4double                    distanceToBoundaryBox(const G4UMeshfitBoundaryBox * aBdrBox, const G4ThreeVector & aPoint,
                                                      const G4ThreeVector & aVector);
    //calculate the distance to a boundary box along the vector. if inside: 0.0


    G4double                    distanceToBoundingSphere(const G4UMeshfitBoundingSphere *aBdSph, const  G4ThreeVector & aPoint);
    //calculate the isotropic distance to a bounding sphere, if inside: 0.0

    G4double                    distanceToBoundingSphere(const G4UMeshfitBoundingSphere * aBdSph,const  G4ThreeVector & aPoint,
                                                         const G4ThreeVector & aVector);
    //calculate the distance to a bounding spherealong the vector , if inside: 0.0

//    G4double                    distanceToRealBoundaryBox(const G4ThreeVector & aPoint);
    //calculate the isotropic distance to \a m_RealBoundary

//    G4double                    distanceToRigidBoundingSphere(const G4ThreeVector & aPoint);
    //calculate the isotropic distance to \a m_RigidBoundingSphere

//    G4double                    distanceToLooseBoundingSphere(const  G4ThreeVector & aPoint);
    //calculate the isotropic distance to \a m_LooseBoundingSphere

    G4bool                      isIntersectBoundingSphere(const G4UMeshfitBoundingSphere * aBdSph,const  G4ThreeVector & aPoint,
                                                          const G4ThreeVector & aVector);
    //calculate if the ray(a point with a vector) intersect with the bounding sphere

    void                        createChildrenVoxels();
    //method to create children voxels. the real boundary box should be calcualted first

    G4UMeshfitVoxel *           getChildVoxel(G4int IdxX, G4int IdxY, G4int IdxZ);
    //get a voxel from the list using the index

    inline G4bool               hasChildrenVoxels() {return !m_ChildrenVoxels.empty();};
    //return true if has children voxels

    inline G4bool               hasElements() {return !m_ElementList.empty();};
    //return true if has elements

    vector <G4ThreeVector>      calIntersectPointWithSlicePlanes (const G4ThreeVector & aPoint,
                                                                  const G4ThreeVector & aVector);
    //calculate the intersect points with voxel slice planes

private:
    vector <G4int>              m_ElementList;
    //Vector to contain all elements in this voxel,here use
    //"Index inside mother logic volume array" to refer to the elmenet
    //currently support Tetrahedron
    //after sorting, the elements full contained inside children voxel will be excluded;
    vector< G4UMeshfitVoxel * > m_ChildrenVoxels;
    //a 3D dynamic of children voxel which is made as 1D array,
    //and using m_IntervalX/Y/Z to index the correct one
    G4LogicalVolume *           m_MotherLogicVol;
    //mother logical volume which contains the mesh


    G4double                    m_MultiplierX;
    G4double                    m_MultiplierY;
    G4double                    m_MultiplierZ;
    //Multiplier for X, Y, Z coordiantes, to make the coordinates of a point
    //able to index the voxel

    G4UMeshfitBoundaryBox *     m_RealBoundary;
    //real boundary of the mesh with all element contained in this voxel
    G4UMeshfitBoundaryBox *     m_VoxelBoundary;
    //the voxel boundary which is larger than the real boundary
    G4UMeshfitBoundingSphere *  m_RigidBoundingSphere;
    //a "rigid" bounding sphere which contain all elements in it
    //when using it to calcualte safety, no need to calcualted safety to
    //partly contained elements
    G4UMeshfitBoundingSphere *  m_LooseBoundingSphere;
    //a "Loose" bounding sphere which contain only nodes of elements inside m_RealBoundary
    //when using it to calcualte safety,it needs to calcualted safety to
    //partly contained elements to garantee the correctness.

    G4int                       m_IntervalX;
    G4int                       m_IntervalY;
    G4int                       m_IntervalZ;
    //divided interval in X, Y, Z directions
    G4double                    m_Tolerance;
    //Tolerance to build the voxels. using the tolerance in mother logical volume
    G4int                       m_Level;
    //current level of voxel. The Top is 0, then goes to 1, 2, 3,...

    //###Static members, should be accessible by all voxel objects
    static G4int                SubdivideCriteria;
    //if Number of elements inside this voxel exceed this criteria,
    //children voxels will be generated
    //static member which is share by all objects
    static G4int                VerboseLevel;
    // verbose level for output messages

    static vector <G4bool>      CheckList;
    //list of boolean, marks for all daugther PV in mother LV if they are
    //calculated (for safetyDistance or Step),
    //the size() of this container SHOULD identical with GetNoDaughter










};

#endif // G4UMESHFITVOXEL_HH
