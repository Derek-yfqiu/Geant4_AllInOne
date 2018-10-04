#include "G4UMeshVTKLegacyReader.hh"

#include <strstream>



string trim(string s)
{
    if (s.empty()) return s;
    //remove the space a the beginning
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    //remvoe the space at the end
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

G4UMeshVTKLegacyReader::G4UMeshVTKLegacyReader()
{
    Init();
}

G4UMeshVTKLegacyReader::G4UMeshVTKLegacyReader(string & FileName)
{
    SetFileName(FileName);
}

void    G4UMeshVTKLegacyReader::Init()
{
    m_FileName.clear();
    m_DimensionFactor = 1;
    m_PointList.clear();
    m_Connectivities.clear();
    m_ConnectivityIndex.clear();
    m_PolyhedronFaces.clear();
    m_PolyhedronElements.clear();
    m_verbose = 0;

    currentLine.clear();
    m_CellConnectivityRaw.clear();
    m_CellTypeRaw.clear();

    //matching VTK element type to G4UM::ElementType
    VTK2ElementType.clear();
    VTK2ElementType[10] = G4UM::TETRA_4;
    VTK2ElementType[12] = G4UM::HEXA_8;
    VTK2ElementType[13] = G4UM::PENTA_6;
    VTK2ElementType[14] = G4UM::PYRA_5;
    VTK2ElementType[24] = G4UM::TETRA_10;
    VTK2ElementType[25] = G4UM::HEXA_20;
}


/*!
 * \brief G4UMeshVTKLegacyReader::ProcessFile
 *  method to precess the mesh file
 * \return \a G4UMReadStatus Ok, warning or failed
 */
G4UM::ReadStatus G4UMeshVTKLegacyReader::ProcessFile()
{
    //check the file name
    if (m_FileName.empty()) return G4UM::Read_failed;

    //open the file
    m_FileStream.open(m_FileName.c_str());
    if(m_FileStream.is_open())
        cout << "File " << m_FileName << " is open.\n";
    else    {
        G4ExceptionDescription msg;
        msg << "File " << m_FileName << " Open Failed!\n";
        G4Exception("G4UMeshVTKLegacyReader::ProcessFile()",
                    "umesh", FatalException, msg);
        return G4UM::Read_failed;
    }

    //######### start to read the file#########

    G4UM::LineType aLineType = getLineAndCheck() ;
    VTKLegacyKey   aKey;
    bool isOk = true;

    while (aLineType != G4UM::Line_eof)
    {
        switch (aLineType)
        {
        case G4UM::Line_keyword:
            aKey = checkKeys(currentLine);
            switch (aKey)
            {
            case Key_heading :
                isOk = readHeading();
                if (m_verbose >=1 )
                    cout << "...Reading Heading End."<<G4endl;
                break;
            case Key_points:
                isOk = readPoints();
                if (m_verbose >=1 )
                    cout << "...Reading Points End."<<G4endl;
                break;
            case Key_cells:
                isOk = readCells();
                if (m_verbose >=1 )
                    cout << "...Reading Cells End."<<G4endl;
                break;
            case Key_cell_types:
                isOk = readCellTypes();
                if (m_verbose >=1 )
                    cout << "...Reading Cell types End."<<G4endl;
                break;
            default: //skip useless lines
                if (m_verbose >=3 )
                    cout << "#Skip Line: "<<currentLine<<G4endl;
                break ;
            }
            if (!isOk) return G4UM::Read_failed;
            aLineType = getLineAndCheck() ; //read a new line
            break;
        case G4UM::Line_blank:
            if (m_verbose >=2 )
                cout << "#Skip a blanck Line. "<<G4endl;
            aLineType = getLineAndCheck();
            break;
        default:
            aLineType = getLineAndCheck();
            break;
        }
    }

    //process the connectivities
    if (!processConnectivities()) return G4UM::Read_failed;

    return G4UM::Read_ok;
}
/*!
 * \brief G4UMeshVTKLegacyReader::getLine
 *  get a line
 * \return
 */
void  G4UMeshVTKLegacyReader::getLine()
{
    currentLine.clear();
    getline(m_FileStream, currentLine);
    if (m_verbose >=2 )
        cout << currentLine <<endl;
}


/*!
 * \brief G4UMeshVTKLegacyReader::getLine
 *  get a line and check the Line Type
 * \return Line type : eof/blank/keyword
 */
G4UM::LineType  G4UMeshVTKLegacyReader::getLineAndCheck()
{
    getLine();
    if (m_FileStream.eof())
        return G4UM::Line_eof;
    if (trim(currentLine).empty())  //if trimmed is empty
        return G4UM::Line_blank;
    else return G4UM::Line_keyword;
}

/*!
 * \brief G4UMeshVTKLegacyReader::checkKeys
 *  check the key word of the line for swtiching the reading functions
 * \param aString a line
 * \return the key word type
 */
VTKLegacyKey G4UMeshVTKLegacyReader::checkKeys(string & aString)
{
    if (aString.find("vtk DataFile") != string::npos)
        return Key_heading;
    if (aString.find("POINTS") != string::npos)
        return Key_points;
    if (aString.find("CELLS") != string::npos)
        return Key_cells;
    if (aString.find("CELL_TYPES") != string::npos)
        return Key_cell_types;
    return Key_unknow;
}


bool  G4UMeshVTKLegacyReader::readHeading()
{
    getLine();  //skip the file discription
    getLine();  //get the BINARY/ASCII
    if (currentLine.find("ASCII") == string::npos){ //if not found
        G4Exception("G4UMeshVTKLegacyReader::readHeading",
                    "umesh", FatalException, "The File format should be ASCII!");
        return false;
    }
    getLine();  //get the data-set format: UNSTRUCTURED_GRID
    if (currentLine.find("UNSTRUCTURED_GRID") == string::npos){ //if not found
        G4Exception("G4UMeshVTKLegacyReader::readHeading",
                    "umesh", FatalException, "The data_set should be UNSTRUCTURED_GRID!");
        return false;
    }
    return true;
}

bool G4UMeshVTKLegacyReader::readPoints()
{
    //get the point number in the current line: e.g. POINTS 925 float
    int NumberOfPoints = 0; // total point account
    stringstream tmpsstr (stringstream::in | stringstream::out);
    string dumy;

    tmpsstr << currentLine;
    tmpsstr >> dumy ; //dump the "POINTS"
    tmpsstr >> NumberOfPoints; //save the point acount
    tmpsstr.str("") ; //clear

    //read the points
    G4cout <<"Dimension Factor"<< m_DimensionFactor <<G4endl;

    for (int i=0; i<NumberOfPoints; i++) {
        G4ThreeVector aPoint;
        double x=0, y=0, z=0;
        m_FileStream >> x >> y >> z ;
        if (m_FileStream.eof()) {
            G4Exception("G4UMeshVTKLegacyReader::readPoints",
                        "umesh", FatalException, "The mesh file ended unexpectedly!");
            return false;
        }
        aPoint.set( x*m_DimensionFactor , y*m_DimensionFactor , z*m_DimensionFactor );
        m_PointList.push_back(aPoint); //push to the point list
        if (m_verbose >=3 )
            cout << "#Read Point:  "<<x<<"\t"<<y<<"\t"<<z<<G4endl;
    }
    return true;
}

bool  G4UMeshVTKLegacyReader::readCells()
{
    //get the number of cells and number of data
    int NumberOfCells = 0, NumberOfData = 0;
    stringstream tmpsstr (stringstream::in | stringstream::out);
    string dumy;

    tmpsstr << currentLine;
    tmpsstr >> dumy ; //dump the "CELLS"
    tmpsstr >> NumberOfCells >> NumberOfData;
    tmpsstr.str("") ; //clear

    //read the unprocessed cell connectivities
    for (int i=0; i<NumberOfData; i++) {
        int tmpInt = 0;
        m_FileStream >> tmpInt ;
        if (m_FileStream.eof()) {
            G4Exception("G4UMeshVTKLegacyReader::readCells",
                        "umesh", FatalException, "The mesh file ended unexpectedly!");
            return false;
        }
        m_CellConnectivityRaw.push_back(tmpInt);
        if (m_verbose >=3 )
            cout << "#Cell Conn:  "<<tmpInt<<G4endl;
    }

    return true;
}

bool  G4UMeshVTKLegacyReader::readCellTypes()
{
    //get the number of cells
    int NumberOfCells = 0; // total cell account
    stringstream tmpsstr (stringstream::in | stringstream::out);
    string dumy;

    tmpsstr << currentLine;
    tmpsstr >> dumy ; //dump the "CELL_TYPES"
    tmpsstr >> NumberOfCells; //save the point acount
    tmpsstr.str("") ; //clear

    //read the unprocessed cell connectivities
    for (int i=0; i<NumberOfCells; i++) {
        int tmpInt = 0;
        m_FileStream >> tmpInt ;
        if (m_FileStream.eof()) {
            G4Exception("G4UMeshVTKLegacyReader::readCellTypes",
                        "umesh", FatalException, "The mesh file ended unexpectedly!");
            return false;
        }
        m_CellTypeRaw.push_back(tmpInt);
        if (m_verbose >=3 )
            cout << "#Cell Type:  "<<tmpInt<<G4endl;
    }
    return true;
}

bool  G4UMeshVTKLegacyReader::processConnectivities()
{
    //see if data exits
    if (m_CellConnectivityRaw.empty() || m_CellTypeRaw.empty() ) {
        G4Exception("G4UMeshVTKLegacyReader::processConnectivities",
                    "umesh", FatalException, "The raw connectvities or cell type data should not be empty!");
        return false;
    }

    //start to process the connectvities data
    unsigned int  RawConnIndex =0;  //index for m_CellConnectivityRaw
    unsigned int CellIndex = 0;    //index for cells, and for m_CellTypeRaw
    while (RawConnIndex < m_CellConnectivityRaw.size())
    {

        //first check the element type
        //if cell type cannot been found, VTK2ElementType[] returns a 0, which match Element Type UNKNOWN
        G4UM::ElementType aElmType = VTK2ElementType[m_CellTypeRaw [CellIndex]] ;
        if (aElmType == G4UM::UNKNOWN_TYPE) {
            G4Exception("G4UMeshVTKLegacyReader::processConnectivities",
                        "umesh", FatalException, "The non-3D cell type is not supported!");
            return false;
        }
        //push the element type to the m_Connectivities
        m_ConnectivityIndex.push_back(m_Connectivities.size()); //push the conn index
        m_Connectivities.push_back(int(aElmType));
        if (m_verbose >=4 )
            cout << "#Cell "<<CellIndex<<"\t Type: "<<int(aElmType) <<"\t Conn: ";

        int tmpNbofNodes = m_CellConnectivityRaw[RawConnIndex]; //get the Nb of nodes in this element

        //push the connectivities into the m_Connectivities
        for (int i=0; i<tmpNbofNodes; i++) {
            RawConnIndex ++ ;
            m_Connectivities.push_back(m_CellConnectivityRaw[RawConnIndex]);
            if (m_verbose >=4 )
                cout << m_CellConnectivityRaw[RawConnIndex]<<"\t";
        }
        if (m_verbose >=4 )
            cout << G4endl;
        //move
        RawConnIndex ++ ; //index to VTK element type of next element
        CellIndex ++;     //index to the next cell
    }

    //check if correct data are read
    if (m_Connectivities.size() != m_CellConnectivityRaw.size()) {
        G4Exception("G4UMeshVTKLegacyReader::processConnectivities",
                    "umesh", FatalException, "Data not correctly processed !");
        return false;
    }
    return true;
}





