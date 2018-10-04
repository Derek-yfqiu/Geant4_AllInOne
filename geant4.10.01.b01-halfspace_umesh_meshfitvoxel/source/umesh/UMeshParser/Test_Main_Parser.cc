#include <G4UMeshVTKLegacyReader.hh>
#include <G4UMeshParser.hh>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Box.hh>
#include <G4SystemOfUnits.hh>
#include <globals.hh>
#include <string>
#include <G4String.hh>


int main (int argc, char **argv) {

    G4String FileName = "/mnt/shared/Geant4/Sphere_tet_vtkLegacy.vtk";
    G4int    verbose = 4;

//    G4UMeshVTKLegacyReader aVTKReader;
//    aVTKReader.SetFileName(FileName);
//    aVTKReader.SetVerboseLevel(verbose);
//    aVTKReader.ProcessFile();

    G4double world_sizeXY = 100*cm;
    G4double world_sizeZ  = 100*cm;
    G4double A = 1.01*g/mole;
    G4double Z;
    G4Material* world_mat = new G4Material("H1", Z = 1.0, A, 0.0001*g/cm3);

    G4Material * Tet_mat = new G4Material("H1", Z = 1.0, A, 1.0*g/cm3);
    G4Box* solidWorld =
      new G4Box("World",                       //its name
         0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

    G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,          //its solid
                          world_mat,           //its material
                          "World");            //its name

    G4UMeshParser aParser;
    aParser.setCheckOverlap(G4bool(false));
    aParser.setLogicalEnvelop(logicWorld);
    aParser.setMaterial(Tet_mat);
    aParser.setVerboseLevel(G4int(4));
    aParser.readMesh(FileName, G4UM::Mesh_VTKLegacy);




    return 0;
}
