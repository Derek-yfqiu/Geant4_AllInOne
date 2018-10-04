#ifndef G4UMESHPARSER_HH
#define G4UMESHPARSER_HH

#include "G4UMeshVTKLegacyReader.hh"
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>

#include <string>

class G4UMeshParser
{
public:
    G4UMeshParser();
    ~G4UMeshParser();




    bool                readMesh (string & FileName , G4UM::MeshFormat theFormat );
    //method to read the mesh, for outside calling
    bool                readVTKLegacyMesh (string & FileName);
    //read the mesh in VTK legacy format
    bool                checkMesh ();
    //method to check the mesh if satisfied for the limites
    bool                buildMesh();
    //to build the mesh and place into the logical volume

    inline void         setMaterial (G4Material * aMaterial)
                        {m_Material = aMaterial;};
    //set the material for the whole mesh

    G4Material *        getMaterial () {return m_Material;};
    //get the material of the mesh

    inline void         setLogicalEnvelop (G4LogicalVolume * aLogicalEnvelop)
                        {m_LogicalEnvelop = aLogicalEnvelop;};
    //set the enveloped for putting the mesh element into it

    G4LogicalVolume *   getLogicalEnvelop (){return m_LogicalEnvelop;};
    //get the logical enveloped (when the mesh is set)

    inline void         setVerboseLevel(G4int  verboseLevel)
                        {m_verbose = verboseLevel;};
    //set the verbose level

    inline void         setCheckOverlap(G4bool  isCheckOverlap)
                        {m_checkOverlap = isCheckOverlap;};
    //set the verbose level

    void                setDimFactor (double aFactor );
//    set the dimension factor to convert to mm
//    factor might be set before or after reading mesh;
//    in case before, set the value to m_DimensionFactor;
//    in case after, set and change the dimension

    void                calBoundaryBox(G4ThreeVector & LowerPoint, G4ThreeVector & UpperPoint );

    vector <G4VPhysicalVolume *> getElementPVList () {return m_ElementPVList;};

    void                setTranslation(const G4ThreeVector aTranslation)
                        {m_Translation = aTranslation;};

    void                setRotationMatrix(G4RotationMatrix * aRotationMatrix )
                        {m_RotateMatrix = aRotationMatrix;}


    //get the data
    inline vector <G4ThreeVector> *     getPointList()
                        {return m_PointList;};
//    inline vector <G4ThreeVector>      getPointList()
//                        {return m_PointList;};
    inline vector <int > *              getConnectivities()
                        {return m_Connectivities;};
    inline vector <int > *              getConnectivityIndex()
                        {return m_ConnectivityIndex;};
    inline vector <int > *              getPolyhedronFaces()
                        {return m_PolyhedronFaces;};
    inline vector <int > *              getPolyhedronElements()
                        {return m_PolyhedronElements;};


private:
    G4Material *        m_Material;
    G4LogicalVolume *   m_LogicalEnvelop;
    G4int               m_verbose;
    G4bool              m_checkOverlap;
    //readers
    G4UMeshVTKLegacyReader * m_VTKLegacyReader;

    vector <G4VPhysicalVolume *>  m_ElementPVList;

    //container to get the data :
    vector              <G4ThreeVector> * m_PointList;
    //list of point to construct the mesh
    vector              <int> * m_Connectivities;
    //connectivities for the elements, for each element: element type follows a list of node id(start from 0)
    vector              <int> * m_ConnectivityIndex;
    //starting location m_Connectivities for each element, start from 0
    vector              <int> * m_PolyhedronFaces;
    //A vector to contain Arbitary polyhedron faces: using CGNS format – start with a “number of node of this face” then follow the node id (start from 0)
    vector              <int> * m_PolyhedronElements;
    //A vector to contain arbitrary polyhedron element: using CGNS format – start with “number of faces of this element” then follow the face id (start from 0)

    map <G4UM::ElementType, bool>  supportTypes;
    //to check if the Element type is supported or not

    double              m_DimensionFactor;
    //dimension factor for converting to mm

    G4RotationMatrix *  m_RotateMatrix;
    G4ThreeVector       m_Translation;
};

#endif // G4UMESHPARSER_HH
