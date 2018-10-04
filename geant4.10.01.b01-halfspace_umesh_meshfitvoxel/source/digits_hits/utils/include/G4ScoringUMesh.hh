#ifndef G4SCORINGUMESH_H
#define G4SCORINGUMESH_H


#include "globals.hh"
#include "G4VScoringMesh.hh"
#include "G4RotationMatrix.hh"
//#include "G4UMeshGeneral.hh"



class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VPrimitiveScorer;
class G4VScoreColorMap;
class G4UMeshParser;
#include <vector>

class G4ScoringUMesh: public G4VScoringMesh
{
public:
    G4ScoringUMesh(G4String wName);
    ~G4ScoringUMesh();

public:
    virtual void  Construct(G4VPhysicalVolume* fWorldPhys);

    void          List() const;
    void          Draw(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap, G4int axflg=111);

    virtual void          openUMeshFile();
    void          DrawColumn(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap,
                             G4int idxProj, G4int idxColumn) ;


//    void          setFileName(G4String & aFileName) {fMeshFileName = aFileName;};
//    void          setMeshFormat (G4UM::MeshFormat & aFormat) {fMeshFormat = aFormat;};
//    void          setDimensionFactor (const G4double aFactor) {fDimensionFactor = aFactor;};
    inline G4UMeshParser * getUMeshParser()
                    {return fUMeshParser;};

private:

    G4UMeshParser *  fUMeshParser;
    //pointer to unstructured mesh parser


//    G4String fMeshFileName;
//    G4UM::MeshFormat fMeshFormat;

//    G4double      fDimensionFactor;
//    vector <G4VPhysicalVolume * > fElementPVList;



};

#endif // G4SCORINGUMESH_H
