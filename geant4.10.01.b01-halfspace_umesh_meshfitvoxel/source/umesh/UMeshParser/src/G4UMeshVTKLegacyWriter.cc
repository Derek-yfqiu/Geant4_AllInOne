#include "G4UMeshVTKLegacyWriter.hh"
#include "G4UMeshWriter.hh"

G4UMeshVTKLegacyWriter::G4UMeshVTKLegacyWriter()
{
    m_FileName = "Untitled.vtk";

    ElementType2VTK.clear();
    ElementType2VTK[G4UM::TETRA_4] = 10;
    ElementType2VTK[G4UM::HEXA_8] = 12;
    ElementType2VTK[G4UM::PENTA_6] =13 ;
    ElementType2VTK[G4UM::PYRA_5] = 14;
    ElementType2VTK[G4UM::TETRA_10] =24 ;
    ElementType2VTK[G4UM::HEXA_20] =25 ;
}

bool        G4UMeshVTKLegacyWriter::dumpToFile(G4String aPrimitiveScoreName)
{
    //check the mesh first
    if (m_PointList == 0 || m_PointList->empty()) {
        G4cerr <<"Error: dump file failed because no Point data! "<<G4endl;
        return false;
    }
    if (m_Connectivities == 0 || m_Connectivities->empty()) {
        G4cerr <<"Error: dump file failed because no Connectivity data! "<<G4endl;
        return false;
    }
    if (m_ConnectivityIndex == 0 || m_ConnectivityIndex->empty()) {
        G4cerr <<"Error: dump file failed because no Connectivity data! "<<G4endl;
        return false;
    }

    //Open the file
    m_FileStream.open( m_FileName, ios_base::trunc); //to overide the
    if(m_FileStream.is_open())
        cout << "File " << m_FileName << " is open to export vtk mesh.\n";
    else    {
        cout << "Error opening " << m_FileName << ".\n";
        return false;
    }

    //dump the heading info
    m_FileStream << "# vtk DataFile Version 3.0" <<endl;
    m_FileStream << "Geant4 unstructured mesh tally output, programer Yuefeng Qiu" <<endl;
    m_FileStream << "ASCII" <<endl;
    m_FileStream << "DATASET UNSTRUCTURED_GRID" <<endl;

    //dump point data
    m_FileStream << "POINTS " <<m_PointList->size() << " float"<<endl;
    for (unsigned int i=0; i<m_PointList->size(); i++)
        m_FileStream << m_PointList->at(i).x() << " "
                     << m_PointList->at(i).y() << " "
                     << m_PointList->at(i).z() << endl;


    //dump cell data
    //the mechanism of the m_ConnectivityIndex and m_Connectivities please check the SALOME MED for more information
    m_FileStream << "CELLS " <<m_ConnectivityIndex->size()<< " "<< m_Connectivities->size()<<endl;
    for (unsigned int i=0; i<m_ConnectivityIndex->size(); i++) {
        G4UM::ElementType aElementType = (G4UM::ElementType) m_Connectivities->at(m_ConnectivityIndex->at(i)); //get element type
        int aNodeCnt = mapNodesPerElement[aElementType];
        if (aNodeCnt == 0) {
            G4cerr <<"Error: Found an element without node! "<<G4endl;
            return false;
        }
        m_FileStream << aNodeCnt;
        for ( int j=1; j<=aNodeCnt; j++) {  //start from 1
            m_FileStream << " " << m_Connectivities->at(m_ConnectivityIndex->at(i) + j);
        }
        m_FileStream << endl;
    }

    //dump cell type
    m_FileStream << "CELL_TYPES " <<m_ConnectivityIndex->size()<<endl;
    for (unsigned int i=0; i<m_ConnectivityIndex->size(); i++) {
        G4UM::ElementType aElementType = (G4UM::ElementType) m_Connectivities->at(m_ConnectivityIndex->at(i)); //get element type
        int aVTKType = ElementType2VTK[aElementType];
        if (aVTKType == 0) {  //if not found, the map return 0
            G4cerr <<"Error: Found an unknown or unsupported cell type! "<<G4endl;
            return false;
        }
        m_FileStream << aVTKType <<endl;
    }

    //see if has mesh data, if not return
    if (m_FieldOnMeshElement.empty()) {
        G4cout<< "File " << m_FileName << " Export done.\n";
        return true;
    }
    else {
        m_FileStream << "CELL_DATA " <<m_ConnectivityIndex->size()<<endl;

        //!!here perform without checking the bound!!
        //please take care to assign the data!
        std::map<G4int, G4double*> * aScore; //pointer to get the data
        if (aPrimitiveScoreName.isNull())
        {  //output all the result into one file
            for (std::map<G4String,G4THitsMap<G4double>* >::iterator itr = m_FieldOnMeshElement.begin();
                 itr != m_FieldOnMeshElement.end(); itr ++)
            {
                G4String aName = itr->first;
                aScore = itr->second->GetMap();
                //output the data heading
                m_FileStream << "SCALARS " <<aName<< " float 1 "<<endl; //SCALARS dataName dataType numComp
                m_FileStream << "LOOKUP_TABLE default"<<endl;
                //find the unit
                G4double aUnit = 1.0; //default unit
                std::map<G4String, G4double>::iterator aUnitItr = m_FieldUnits.find(aName);
                if (aUnitItr != m_FieldUnits.end())  //if found, use the unit, else use default
                    aUnit = aUnitItr->second;
                //output the results
                for (unsigned int i=0; i<m_ConnectivityIndex->size(); i++) {
                    std::map<G4int, G4double*>::iterator aValue = aScore->find(i);
                    if(aValue == aScore->end()) {
                        m_FileStream  << "0.000000000000000" <<endl; //use this strange 0 for result checking
                    } else {
                        m_FileStream << *(aValue->second)/aUnit <<endl ;
                    }
                }
            }
        }
        else {
            //find the score and only ouptut it
            std::map <G4String,G4THitsMap<G4double>* >::iterator itr = m_FieldOnMeshElement.begin();
            for (;itr != m_FieldOnMeshElement.end(); itr ++)
            {
                G4String aName = itr->first;
                aScore = itr->second->GetMap();
                if (aName = aPrimitiveScoreName)
                {
                    //output the data heading
                    m_FileStream << "SCALARS " <<aName<< " float 1 "<<endl; //SCALARS dataName dataType numComp
                    m_FileStream << "LOOKUP_TABLE default"<<endl;
                    //find the unit
                    G4double aUnit = 1.0; // default unit
                    std::map<G4String, G4double>::iterator aUnitItr = m_FieldUnits.find(aName);
                    if (aUnitItr != m_FieldUnits.end())
                        aUnit = aUnitItr->second;
                    //for testing
//                    std::map<G4int, G4double*>::iterator aValueItr = aScore->begin();
//                    for (; aValueItr != aScore->end(); ++aValueItr) {
//                        G4cout<<"Key\t"<<aValueItr->first<<"\t\tData\t"<<aValueItr->second<<G4endl;
//                    }

                    //output the results
                    for (unsigned int i=0; i<m_ConnectivityIndex->size(); i++) {
                        std::map<G4int, G4double*>::iterator aValue = aScore->find(i);
                        if(aValue == aScore->end()) {
                            m_FileStream  << "0.000000000000000" <<endl;
                        } else {
                            m_FileStream << *(aValue->second)/aUnit <<endl ;
                        }
                    }
                    break;
                }
            }
            if (itr == m_FieldOnMeshElement.end())
                G4cerr<< "Warning: cannot find the Primitive Score "<<aPrimitiveScoreName <<" to output!"<<G4endl;
        }
    }
    G4cout<< "File " << m_FileName << " Export done.\n";
    return true;
}
