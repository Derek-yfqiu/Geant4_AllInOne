#ifndef G4UMESHVTKLEGACYREADER_HH
#define G4UMESHVTKLEGACYREADER_HH
#include <stdio.h>
#include <iostream>
#include <string>

#include "G4UMeshReader.hh"



enum VTKLegacyKey {
    Key_heading,
    Key_points,
    Key_cells,
    Key_cell_types,
    Key_unknow
};

class G4UMeshVTKLegacyReader : public G4UMeshReader
{
public:
    G4UMeshVTKLegacyReader();
    G4UMeshVTKLegacyReader(string &FileName);
    virtual G4UM::ReadStatus            ProcessFile();
    //process the file

private:
    virtual void                        Init();
    void                                getLine();
    G4UM::LineType                      getLineAndCheck();
    //get a line from the file and judge the line type
    VTKLegacyKey                        checkKeys(string & aString);
    //check the key word
    bool                                readHeading();
    bool                                readPoints();
    bool                                readCells();
    bool                                readCellTypes();
    bool                                processConnectivities();


private:
    string  currentLine;  //a Line get from the file
    vector <int> m_CellConnectivityRaw;
    vector <int> m_CellTypeRaw;
    map <int, G4UM::ElementType> VTK2ElementType;



};

#endif // G4UMESHVTKLEGACYREADER_HH
