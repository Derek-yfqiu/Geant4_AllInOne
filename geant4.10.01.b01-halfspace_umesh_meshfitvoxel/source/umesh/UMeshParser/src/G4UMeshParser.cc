#include "G4UMeshParser.hh"

#include <G4Tet.hh>
#include <G4UMeshHex.hh>
#include <G4UMeshPyrm.hh>
#include <G4UMeshPent.hh>
//#include <G4PVPlacement.hh>

#include <stdio.h>
#include <stdlib.h>

G4UMeshParser::G4UMeshParser()
{
    m_Material = 0;
    m_LogicalEnvelop = 0;
    m_VTKLegacyReader = 0;
    m_verbose = 1; //default verbose level
    m_checkOverlap = false;
    m_DimensionFactor = 1.0;

    m_RotateMatrix = 0;
    m_Translation = G4ThreeVector(0);

    m_PointList = 0 ;
    m_Connectivities = 0;
    m_ConnectivityIndex = 0;
    m_PolyhedronElements = 0;
    m_PolyhedronFaces = 0;

    supportTypes.clear();
    supportTypes[G4UM::TETRA_4]  = true;
    supportTypes[G4UM::TETRA_10] = false;
    supportTypes[G4UM::PYRA_5]   = true;
    supportTypes[G4UM::PYRA_14]  = false;
    supportTypes[G4UM::PENTA_6]  = true;
    supportTypes[G4UM::PENTA_15] = false;
    supportTypes[G4UM::PENTA_18] = false;
    supportTypes[G4UM::HEXA_8]   = true;
    supportTypes[G4UM::HEXA_20]  = false;
    supportTypes[G4UM::HEXA_27]  = false;
    supportTypes[G4UM::PYRA_13]  = false;
    supportTypes[G4UM::POLYHED_n]= false;
    supportTypes[G4UM::UNKNOWN_TYPE]= false;
}

G4UMeshParser::~G4UMeshParser()
{
    if (m_VTKLegacyReader) delete m_VTKLegacyReader;
    if (m_RotateMatrix) delete m_RotateMatrix;
}


/*!
 * \brief readMesh
 * method to read the mesh, for outside calling
 * this method is the main method for process mesh,
 * it carried out the mesh reading, checking and building
 * by calling other method
 * \param FileName
 * \param theFormat Mesh format
 * \return  false if error
 */
bool  G4UMeshParser::readMesh (string & FileName , G4UM::MeshFormat theFormat )
{
    if (FileName.empty()) {
        G4Exception("G4UMeshParser::readMesh",
                    "umesh", FatalErrorInArgument, "The file name is empty!");
        return false;
    }

    //switch to the corresponding reader
    bool isOk = true;
    switch (theFormat)
    {
    case G4UM::Mesh_VTKLegacy:
        isOk = readVTKLegacyMesh(FileName);
        break;
    default:
        G4Exception("G4UMeshParser::readMesh",
                    "umesh", FatalErrorInArgument, "Unrecognized format!");
        return false;
    }
    if (!isOk) return false;

    //build the mesh into the logical volume
//    if (!buildMesh()) return false;// comment out to execute it in after reading mesh

    return true;
}

/*!
 * \brief readVTKLegacyMesh
 *read the mesh in VTK legacy format
 * \param FileName
 * \return  false if error
 */
bool  G4UMeshParser::readVTKLegacyMesh (string & FileName)
{
    //read the file
    if (!m_VTKLegacyReader)
        m_VTKLegacyReader = new G4UMeshVTKLegacyReader();
    m_VTKLegacyReader->SetFileName(FileName);
    m_VTKLegacyReader->SetVerboseLevel(m_verbose);
    m_VTKLegacyReader->SetDimFactor(m_DimensionFactor);
    if (m_VTKLegacyReader->ProcessFile() == G4UM::Read_failed)
        return false;

    //get the data
    m_PointList = m_VTKLegacyReader->getPointList();
    m_Connectivities = m_VTKLegacyReader->getConnectivities();
    m_ConnectivityIndex = m_VTKLegacyReader->getConnectivityIndex();
    return true;
}

/*!
 * \brief checkMesh
 *method to check the mesh if satisfied for the limites
 * check the element type if supported
 * \return false if error
 */
bool   G4UMeshParser::checkMesh ()
{
    //check the array first
    if (m_PointList->empty() || m_Connectivities->empty()||
            m_ConnectivityIndex->empty())
        return false;

    if (m_verbose >= 4) {
        for (unsigned int i=0; i<m_PointList->size(); i++) {
            G4ThreeVector aPoint = m_PointList->at(i);
            G4cout << "Point "<< i << "\t " << aPoint.x() << "\t "<< aPoint.y() << "\t "<< aPoint.z()<<G4endl;
        }
    }
    //check the element type
    for (unsigned int i=0; i<m_ConnectivityIndex->size(); i++) {
        //check bound
        if (m_ConnectivityIndex->at(i) >= m_Connectivities->size()) {
            G4Exception("G4UMeshParser::checkMesh",
                        "umesh", FatalErrorInArgument, "Connectivity array does not match with index array!");
            return false;
        }
        G4UM::ElementType tmpType =G4UM::ElementType(m_Connectivities->at(m_ConnectivityIndex->at(i))) ;
        if (!supportTypes[tmpType]) {
            G4Exception("G4UMeshParser::checkMesh",
                        "umesh", FatalErrorInArgument, "mesh type not supported!");
            return false;
        }
    }

    return true;
}

/*!
 * \brief buildMesh
 *to build the mesh and place into the logical volume
 * \return
 */
bool G4UMeshParser::buildMesh()
{
    //fisrt check the mesh and the logical volume
    if (!checkMesh()) return false;
    if (!m_LogicalEnvelop ) {
        G4Exception("G4UMeshParser::buildMesh",
                    "umesh", FatalErrorInArgument, "Logical volume mesh envelop not set!");
        return false;
    }
//    if (!m_Material ) {
//        G4Exception("G4UMeshParser::buildMesh",
//                    "umesh", FatalErrorInArgument, "Material not set!");
//        return false;
//    }

    //build the mesh
    for (int i=0; i<m_ConnectivityIndex->size(); i++) {
        G4UM::ElementType tmpType =G4UM::ElementType(m_Connectivities->at(m_ConnectivityIndex->at(i))) ;

        //switch the element type
        switch(tmpType)
        {
        case G4UM::TETRA_4:
        {
            char temp[64];
            sprintf(temp, "%d", i);
            G4String aName = "Tet";
            aName += temp;
            G4ThreeVector p1 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+1));
            G4ThreeVector p2 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+2));
            G4ThreeVector p3 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+3));
            G4ThreeVector p4 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+4));
            G4Tet* aTet = new G4Tet(aName, p1, p2, p3, p4);
            G4LogicalVolume* aLogTet = new G4LogicalVolume (aTet, m_Material, aName);
            //##!!here is a trick, we set the copyNo as i, therefore
            //we can identified them in MultifunctionalDetector using copyNo as index##
            G4VPhysicalVolume * aPhyTet = new G4PVPlacement(m_RotateMatrix, m_Translation, aLogTet,aName, m_LogicalEnvelop, false, i, m_checkOverlap);
            m_ElementPVList.push_back(aPhyTet);
            break;
        }
        case G4UM::PYRA_5:
        {
            char temp[64];
            sprintf(temp, "%d", i);
            G4String aName = "Pyrm";
            aName += temp;
            G4ThreeVector p1 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+1));
            G4ThreeVector p2 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+2));
            G4ThreeVector p3 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+3));
            G4ThreeVector p4 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+4));
            G4ThreeVector p5 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+5));
            G4UMeshPyrm* aPyrm = new G4UMeshPyrm(aName, p1, p2, p3, p4,p5);
            G4LogicalVolume* aLogPyrm = new G4LogicalVolume (aPyrm, m_Material, aName);
            //##!!here is a trick, we set the copyNo as i, therefore
            //we can identified them in MultifunctionalDetector using copyNo as index##
            G4VPhysicalVolume * aPhyPyrm = new G4PVPlacement(m_RotateMatrix, m_Translation, aLogPyrm,aName, m_LogicalEnvelop, false, i, m_checkOverlap);
            m_ElementPVList.push_back(aPhyPyrm);
            break;
        }
        case G4UM::PENTA_6:
        {
            char temp[64];
            sprintf(temp, "%d", i);
            G4String aName = "Pent";
            aName += temp;
            G4ThreeVector p1 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+1));
            G4ThreeVector p2 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+2));
            G4ThreeVector p3 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+3));
            G4ThreeVector p4 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+4));
            G4ThreeVector p5 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+5));
            G4ThreeVector p6 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+6));
            G4UMeshPent* aPent = new G4UMeshPent(aName, p1, p2, p3, p4, p5, p6);
            G4LogicalVolume* aLogPent = new G4LogicalVolume ( aPent, m_Material, aName);
            //##!!here is a trick, we set the copyNo as i, therefore
            //we can identified them in MultifunctionalDetector using copyNo as index##
            G4VPhysicalVolume * aPhyPent = new G4PVPlacement(m_RotateMatrix, m_Translation, aLogPent,aName, m_LogicalEnvelop, false, i, m_checkOverlap);
            m_ElementPVList.push_back(aPhyPent);
            break;
        }
        case G4UM::HEXA_8:
        {
            char temp[64];
            sprintf(temp, "%d", i);
            G4String aName = "Hex";
            aName += temp;
            G4ThreeVector p1 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+1));
            G4ThreeVector p2 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+2));
            G4ThreeVector p3 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+3));
            G4ThreeVector p4 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+4));
            G4ThreeVector p5 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+5));
            G4ThreeVector p6 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+6));
            G4ThreeVector p7 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+7));
            G4ThreeVector p8 = m_PointList->at(m_Connectivities->at(m_ConnectivityIndex->at(i)+8));
            G4UMeshHex* aHex = new G4UMeshHex(aName, p1, p2, p3, p4, p5, p6, p7, p8);
            G4LogicalVolume* aLogHex = new G4LogicalVolume (aHex, m_Material, aName);
            //##!!here is a trick, we set the copyNo as i, therefore
            //we can identified them in MultifunctionalDetector using copyNo as index##
            G4VPhysicalVolume * aPhyHex = new G4PVPlacement(m_RotateMatrix, m_Translation, aLogHex,aName, m_LogicalEnvelop, false, i, m_checkOverlap);
            m_ElementPVList.push_back(aPhyHex);
            break;
        }
        default:
            G4Exception("G4UMeshParser::buildMesh",
                        "umesh", FatalErrorInArgument, "not supported element type!");
            return false;
        }

    }

    return true;
}

/*!
 * \brief G4UMeshParser::calBoundaryBox
 * calculate the boundary box of the mesh,
 * using Cartesian coordinates, and return
 * the LowerPoint and the UpperPoint of the
 * Box
 * \param LowerPoint (out) Lower point of the boundary box,
 * \param UpperPoint (out) Upper point of the boundary box,
 */
void  G4UMeshParser::calBoundaryBox(G4ThreeVector & LowerPoint, G4ThreeVector & UpperPoint )
{
    G4double tmpXmin = 1e99,  tmpYmin = 1e99, tmpZmin = 1e99;
    G4double tmpXmax = -1e99,  tmpYmax = -1e99, tmpZmax = -1e99;

    for (unsigned int i=0; i< m_PointList->size(); i++) {
        G4ThreeVector aPoint =  m_PointList->at(i);
        if (aPoint.x() < tmpXmin) tmpXmin = aPoint.x();
        if (aPoint.y() < tmpYmin) tmpYmin = aPoint.y();
        if (aPoint.z() < tmpZmin) tmpZmin = aPoint.z();
        if (aPoint.x() > tmpXmax) tmpXmax = aPoint.x();
        if (aPoint.y() > tmpYmax) tmpYmax = aPoint.y();
        if (aPoint.z() > tmpZmax) tmpZmax = aPoint.z();
    }

    LowerPoint.set(tmpXmin, tmpYmin, tmpZmin);
    UpperPoint.set(tmpXmax, tmpYmax, tmpZmax);
}

/*!
 * \brief G4UMeshParser::setDimFactor
 *   set the dimension factor to convert to mm
 *   factor might be set before or after reading mesh;
 *   in case before, set the value to m_DimensionFactor;
 *   in case after, set and change the dimension
 * \param aFactor
 */
void G4UMeshParser::setDimFactor (double aFactor )
{
    //in case before reading mesh
    if (m_PointList == 0)
        m_DimensionFactor = aFactor;
    else if (m_PointList->empty())
        m_DimensionFactor = aFactor;
    else {
        for (unsigned int i=0; i<m_PointList->size(); i++) {
            G4ThreeVector aPoint = G4ThreeVector(m_PointList->at(i).x()* aFactor,
                                  m_PointList->at(i).y()* aFactor,
                                  m_PointList->at(i).z()* aFactor);
            m_PointList->at(i) = aPoint;
        }
        m_DimensionFactor = 1; //avoid setting factor again
    }
}


