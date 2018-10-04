#ifndef G4UMESHFITUTILITIES_HH
#define G4UMESHFITUTILITIES_HH

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <vector>


//a class to make a boundary box
//ATTENTION!! this class performs no checkings!
class G4UMeshfitBoundaryBox
{
public:
    G4UMeshfitBoundaryBox
    (G4double Xmin, G4double Ymin,  G4double Zmin, G4double Xmax, G4double Ymax, G4double Zmax);
    //generate the boundary box by specified limites

    G4UMeshfitBoundaryBox
    (G4ThreeVector lowerPoint, G4ThreeVector upperPoint );
    //generate the boundary box by two points

//    G4UMeshfitBoundaryBox
//    (std::vector <G4double> XYZMinMax);
    G4UMeshfitBoundaryBox
    (std::vector <G4ThreeVector> &aNodeList );
    //generate the boundary box by a list of nodes;


    inline G4double getXmin() const {return m_Xmin;};
    inline G4double getYmin() const {return m_Ymin;};
    inline G4double getZmin() const {return m_Zmin;};
    inline G4double getXmax() const {return m_Xmax;};
    inline G4double getYmax() const {return m_Ymax;};
    inline G4double getZmax() const {return m_Zmax;};

    inline void     setXmin(G4double  aDouble) { m_Xmin = aDouble;};
    inline void     setYmin(G4double  aDouble) { m_Ymin = aDouble;};
    inline void     setZmin(G4double  aDouble) { m_Zmin = aDouble;};
    inline void     setXmax(G4double  aDouble) { m_Xmax = aDouble;};
    inline void     setYmax(G4double  aDouble) { m_Ymax = aDouble;};
    inline void     setZmax(G4double  aDouble) { m_Zmax = aDouble;};

    inline G4double getDeltaX() const {return m_Xmax - m_Xmin;};
    inline G4double getDeltaY() const {return m_Ymax - m_Ymin;};
    inline G4double getDeltaZ() const {return m_Zmax - m_Zmin;};
    //get the range of the boundary box

    G4bool          isInside (const G4ThreeVector &aPoint ) const;
    //return if the point is inside the boundary box
    //if a point is lying on the slice surface Xmax, Ymax, Zmax,
    //it is consider NOT inside because normaly this point will be
    //search in next voxel

    void            expandSize(const G4double  aMargin);
    //expand the boundary box, the min=min-aMargin and max= max+aMargin
    //it need to expand the boundary box then adding a tolerance

    void            printMyself ();
    //print the boundary;


private:
    G4double m_Xmin;
    G4double m_Ymin;
    G4double m_Zmin;
    G4double m_Xmax;
    G4double m_Ymax;
    G4double m_Zmax;

};


class G4UMeshfitBoundingSphere
{
public:
    G4UMeshfitBoundingSphere(G4double Radius, G4ThreeVector & Center);
    G4UMeshfitBoundingSphere(std::vector <G4ThreeVector> &aPointList);
    //calculate the Bounding sphere of a list of points

    inline G4double         getRadius() const {return m_Radius;};
    inline G4ThreeVector    getCenterPoint() const {return m_Center;};
    G4bool                  isValid() const {return m_Radius <= 1e-9? false : true;};
    G4double                DistanceToCeneter (const G4ThreeVector &aPoint) const ;
    void                    expandSize(const G4double  aMargin);
    //expand the radius, it need to expand the boundary box then adding a tolerance
private:
    G4double                m_Radius;
    G4ThreeVector           m_Center;
};

#endif // G4UMESHFITUTILITIES_HH
