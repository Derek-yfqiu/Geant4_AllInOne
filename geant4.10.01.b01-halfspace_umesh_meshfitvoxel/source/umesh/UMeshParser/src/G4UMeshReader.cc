#include "G4UMeshReader.hh"



 G4UMeshReader::G4UMeshReader()
{
     m_DimensionFactor = 1; //default using mm
     m_verbose = 1;

     //initiating the number of nodes  for each element type
//     mapNodesPerElement[G4UM::TETRA_4]  = 4;
//     mapNodesPerElement[G4UM::TETRA_10] = 10;
//     mapNodesPerElement[G4UM::PYRA_5]   =5;
//     mapNodesPerElement[G4UM::PYRA_14]  =14;
//     mapNodesPerElement[G4UM::PENTA_6]  =6;
//     mapNodesPerElement[G4UM::PENTA_15] =15;
//     mapNodesPerElement[G4UM::PENTA_18] =18;
//     mapNodesPerElement[G4UM::HEXA_8]   =8;
//     mapNodesPerElement[G4UM::HEXA_20]  =20;
//     mapNodesPerElement[G4UM::HEXA_27]  =27;
//     mapNodesPerElement[G4UM::PYRA_13]  =13;
//     mapNodesPerElement[G4UM::POLYHED_n] =0;

}

 G4UMeshReader::~G4UMeshReader()
{

}
