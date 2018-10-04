#include "G4UMeshWriter.hh"
#include "G4UMeshParser.hh"


G4UMeshWriter::G4UMeshWriter()
{
    m_PointList = 0 ;
    m_Connectivities = 0;
    m_ConnectivityIndex = 0;
    m_PolyhedronElements = 0;
    m_PolyhedronFaces = 0;

    m_FileName = "Untitled";

    //initiating the number of nodes  for each element type
     mapNodesPerElement[G4UM::TETRA_4]  = 4;
     mapNodesPerElement[G4UM::TETRA_10] = 10;
     mapNodesPerElement[G4UM::PYRA_5]   =5;
     mapNodesPerElement[G4UM::PYRA_14]  =14;
     mapNodesPerElement[G4UM::PENTA_6]  =6;
     mapNodesPerElement[G4UM::PENTA_15] =15;
     mapNodesPerElement[G4UM::PENTA_18] =18;
     mapNodesPerElement[G4UM::HEXA_8]   =8;
     mapNodesPerElement[G4UM::HEXA_20]  =20;
     mapNodesPerElement[G4UM::HEXA_27]  =27;
     mapNodesPerElement[G4UM::PYRA_13]  =13;
     mapNodesPerElement[G4UM::POLYHED_n] =0;

}

G4UMeshWriter::~G4UMeshWriter()
{

}

void G4UMeshWriter::setMeshParser(G4UMeshParser * aUMeshParser)
{
    if (aUMeshParser != 0) {
        m_PointList = aUMeshParser->getPointList() ;
        m_Connectivities = aUMeshParser->getConnectivities();
        m_ConnectivityIndex = aUMeshParser->getConnectivityIndex();
        m_PolyhedronElements = aUMeshParser->getPolyhedronElements();
        m_PolyhedronFaces = aUMeshParser->getPolyhedronFaces();
    }
}

