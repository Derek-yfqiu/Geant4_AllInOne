#ifndef G4UMESHWRITER_HH
#define G4UMESHWRITER_HH

#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4THitsMap.hh"
#include "G4UMeshGeneral.hh"

#include <vector>
#include <map>
#include <fstream>



using namespace std;

class G4UMeshParser;

class G4UMeshWriter
{
public:
    G4UMeshWriter();
    ~G4UMeshWriter();

    virtual bool        dumpToFile(G4String  aPrimitiveScoreName = "") = 0;
    //method to dump mesh and result into the file

    void                setFileName(const G4String & aFileName)
                        {m_FileName = aFileName;};
    void                setMeshParser(G4UMeshParser * aUMeshParser);
    void                setFieldData(std::map<G4String,G4THitsMap<G4double>* > & aFieldMap)
                        { m_FieldOnMeshElement  = aFieldMap ;};
    void                setFieldUnits(std::map<G4String,G4double > & aFieldUnits)
                        { m_FieldUnits  = aFieldUnits ;};
protected:

    G4String            m_FileName;
    //file name of the mesh file to be read
    ofstream            m_FileStream;


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

    std::map<G4String,G4THitsMap<G4double>* > m_FieldOnMeshElement; //serval tally data for one mesh, G4THitsMap maps mesh element id to value
    std::map<G4String,G4double > m_FieldUnits; //Unite for the field in m_FieldOnMeshElement




    std::map <G4UM::ElementType, int > mapNodesPerElement;
        //map for match element type with number of nodes

};

#endif // G4UMESHWRITER_HH
