#ifndef G4UMESHREADER_HH
#define G4UMESHREADER_HH

#include "G4UMeshGeneral.hh"

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <fstream>
#include <vector>

class G4UMeshReader
{
public:
     G4UMeshReader();
     ~G4UMeshReader();

public:
     virtual G4UM::ReadStatus        ProcessFile () =0;
     //method to be implement in children class



    inline void         SetFileName (string & FileName )
                        {Init(); m_FileName = FileName;};


    void                SetDimFactor (double  aFactor )
                        {m_DimensionFactor = aFactor;};
    //set the dimension factor to convert to mm

    inline void         SetVerboseLevel (G4int & vbs)
                        {m_verbose = vbs;};
    //set the verbose level

    //get the data
    inline vector <G4ThreeVector> *     getPointList()
                        {return &m_PointList;};
//    inline vector <G4ThreeVector>      getPointList()
//                        {return m_PointList;};
    inline vector <int > *              getConnectivities()
                        {return &m_Connectivities;};
    inline vector <int > *              getConnectivityIndex()
                        {return &m_ConnectivityIndex;};
    inline vector <int > *              getPolyhedronFaces()
                        {return &m_PolyhedronFaces;};
    inline vector <int > *              getPolyhedronElements()
                        {return &m_PolyhedronElements;};

protected:

    virtual void       Init() {};

    vector              <G4ThreeVector> m_PointList;
    //list of point to construct the mesh
    vector              <int> m_Connectivities;
    //connectivities for the elements, for each element: element type follows a list of node id(start from 0)
    vector              <int> m_ConnectivityIndex;
    //starting location m_Connectivities for each element, start from 0
    vector              <int> m_PolyhedronFaces;
    //A vector to contain Arbitary polyhedron faces: using CGNS format – start with a “number of node of this face” then follow the node id (start from 0)
    vector              <int> m_PolyhedronElements;
    //A vector to contain arbitrary polyhedron element: using CGNS format – start with “number of faces of this element” then follow the face id (start from 0)

    string              m_FileName;
    //file name of the mesh file to be read
    ifstream            m_FileStream;

//    static  map <G4UM::ElementType, int > mapNodesPerElement;
    //map for match element type with number of nodes

    double              m_DimensionFactor;
    //dimension factor for converting to mm

    G4int               m_verbose;
    //verbose level for output infomation
    //0-General output, 1-section status, 2-output file content, 3/4-output even debug data
};

#endif // G4UMESHREADER_HH
