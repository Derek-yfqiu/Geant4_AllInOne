#ifndef G4UMESHVTKLEGACYWRITER_HH
#define G4UMESHVTKLEGACYWRITER_HH


#include "G4UMeshWriter.hh"
#include "G4UMeshGeneral.hh"

class G4UMeshVTKLegacyWriter : public G4UMeshWriter
{
public:
    G4UMeshVTKLegacyWriter();

    bool        dumpToFile(G4String  aPrimitiveScoreName = "");


private:
    map <G4UM::ElementType, int> ElementType2VTK;

};

#endif // G4UMESHVTKLEGACYWRITER_HH
