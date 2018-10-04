#ifndef G4UMESHGENERAL_HH
#define G4UMESHGENERAL_HH

#include <map>
//qiu #include <initializer_list>
using namespace std;

//G4UMesh name space
namespace G4UM {

enum ElementType {
    TETRA_4  =10,
    TETRA_10 =11,
    PYRA_5 =12,
    PYRA_14 =13,
    PENTA_6 =14,
    PENTA_15 =15,
    PENTA_18 =16,
    HEXA_8 =17,
    HEXA_20 =18,
    HEXA_27 =19,
    PYRA_13 =21,
    POLYHED_n =22,
    UNKNOWN_TYPE = 0
};



enum ReadStatus {
    Read_ok ,
    Read_warning,
    Read_failed
};

enum LineType {
    Line_keyword,
    Line_comment,
    Line_data,
    Line_blank,
    Line_eof
};

enum MeshFormat {
    Mesh_UnKnown,
    Mesh_VTKLegacy
};

} //END G4UM


#endif // G4UMESHGENERAL_HH

