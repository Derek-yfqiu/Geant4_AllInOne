#include <G4UMeshVTKLegacyReader.hh>
#include <G4UMeshParser.hh>
#include <G4UMeshfitVoxel.hh>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Box.hh>
#include <G4SystemOfUnits.hh>
#include <globals.hh>
#include <string>
#include <G4String.hh>
#include <G4NavigationHistory.hh>
#include <G4NormalNavigation.hh>
#include <G4UMeshHex.hh>

#include <math.h>
#include <time.h>


int main (int argc, char **argv) {

    G4String FileName = "/mnt/shared/Geant4/Sphere_tet_vtkLegacy.vtk";
    G4int    verbose = 2;

//    G4UMeshVTKLegacyReader aVTKReader;
//    aVTKReader.SetFileName(FileName);
//    aVTKReader.SetVerboseLevel(verbose);
//    aVTKReader.ProcessFile();

    G4double world_sizeXY = 1000*cm;
    G4double world_sizeZ  = 1000*cm;
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
    G4VPhysicalVolume* physWorld =
      new G4PVPlacement(0,                     //no rotation
                        G4ThreeVector(),       //at (0,0,0)
                        logicWorld,            //its logical volume
                        "World",               //its name
                        0,                     //its mother  volume
                        false,                 //no boolean operation
                        0,                     //copy number
                        false);        //overlaps checking





//    // ### test of mesh parser###
//    G4UMeshParser aParser;
//    aParser.setCheckOverlap(G4bool(false));
//    aParser.setLogicalEnvelop(logicWorld);
//    aParser.setMaterial(Tet_mat);
//    aParser.setVerboseLevel(verbose);
//    G4double aFactor = 1000.0;
//    aParser.setDimFactor(aFactor);
//    aParser.readMesh(FileName, G4UM::Mesh_VTKLegacy);

//    // ### test the Meshfit voxel###
//    G4UMeshfitVoxel aMeshfitVoxel (logicWorld);  // set the world as mother logical
//    vector <G4int > aElmList;
//    for (int i=0; i< logicWorld->GetNoDaughters(); i++)
//        aElmList.push_back(i); //all the daughters are mesh elements
//    aMeshfitVoxel.setElementList(aElmList);
//    aMeshfitVoxel.setSubdivideCriteria(100);
//    aMeshfitVoxel.setVerboseLevel(verbose);
//    aMeshfitVoxel.sortToVoxels();  //sort the elements to voxels,
//    // testing the Level locate function
//    G4NavigationHistory aHistory;
//    aHistory.SetFirstEntry(physWorld);
//    G4NormalNavigation aNormalNavigator;
//    double start, finish;
//    double MeshfitTime = 0.0, NormalNavigatorTime = 0.0;
//    for (int i=0; i<10000; i++) {
//        //the mesh is from -1 to 1 in three direction if DimFactor = 1;
//        G4double aRandX = (drand48() * 2 - 1) * aFactor;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);
//        G4ThreeVector aDumyPoint;
//        vector <G4bool> aCheckList ;
//        aCheckList.assign(logicWorld->GetNoDaughters(), false);
//        aMeshfitVoxel.setCheckList(aCheckList);
//        start = (double)clock() / CLOCKS_PER_SEC;
//        G4bool isFound = aMeshfitVoxel.LevelLocate(aHistory, aTargetPoint,
//                                                   /*no direction */0,
//                                                   /*not on edge */false,
//                                                   /*return localpoint*/ aDumyPoint);
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        MeshfitTime += finish - start;

//        if (isFound)
//            aHistory.BackLevel();
//        else
//            G4cout << "Not Found!\t";
//        start = (double)clock() / CLOCKS_PER_SEC;

//        isFound = aNormalNavigator.LevelLocate(aHistory, 0, 0, aTargetPoint,
//                                     /*no direction */0,
//                                     /*not on edge */false,
//                                     /*return localpoint*/ aDumyPoint);
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        NormalNavigatorTime += finish - start;
//        if (isFound)
//        {
//            aHistory.BackLevel();
//        }
//        else
//            G4cout << "Not Found too!"<<G4endl;
//    }

//    G4cout <<"CPU Time: Meshfit-" << MeshfitTime << "\t NormalNavigator-"<<NormalNavigatorTime <<G4endl;

//    //testing the ComputeSafety
//    MeshfitTime = 0.0; NormalNavigatorTime = 0.0;
//    G4cout << "########### Testing safety #########"<<G4endl;
//    for (int i=0; i<10000; i++) {
//        //the mesh is from -1 to 1 in three direction if DimFactor = 1;
//        G4double aRandX = (drand48() * 2 - 1) * aFactor*10; //make the random point in larger range
//        G4double aRandY = (drand48() * 2 - 1) * aFactor*10;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor*10;
//        G4cout.precision(12);
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);
//        vector <G4bool> aCheckList ;
//        aCheckList.assign(logicWorld->GetNoDaughters(), false);
//        aMeshfitVoxel.setCheckList(aCheckList);
//        start = (double)clock() / CLOCKS_PER_SEC;
//        G4double aSafety = aMeshfitVoxel.ComputeSafety(aTargetPoint, /*is accurate*/ true, Safety_LooseSphere  );
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        MeshfitTime += finish - start;
//        G4cout <<"Has Safety Distance:\t"<<aSafety<<"\t";

//        start = (double)clock() / CLOCKS_PER_SEC;
//        G4double bSafety = aNormalNavigator.ComputeSafety(aTargetPoint, aHistory);
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        NormalNavigatorTime += finish - start;
//        G4cout <<bSafety<<G4endl ;
//        continue;
//    }

//    G4cout <<"CPU Time: Meshfit-" << MeshfitTime << "\t NormalNavigator-"<<NormalNavigatorTime <<G4endl;


//    //testing the Compute Step
//    MeshfitTime = 0.0; NormalNavigatorTime = 0.0;
//    G4cout << "########### Testing step #########"<<G4endl;
//    for (int i=0; i<10000; i++) {
//        //the mesh is from -1 to 1 in three direction if DimFactor = 1;
//        G4double aRandX = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandXVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandYVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4cout <<"Vector\t"<<aRandXVec<<"\t\t"<<aRandYVec<<"\t\t"<<aRandZVec<<"\t\t" ;

//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);
//        G4ThreeVector aTargetVec = G4ThreeVector(aRandXVec, aRandYVec, aRandZVec);
//        vector <G4bool> aCheckList ;
//        aCheckList.assign(logicWorld->GetNoDaughters(), false);
//        aMeshfitVoxel.setCheckList(aCheckList);
//        start = (double)clock() / CLOCKS_PER_SEC;
//        G4double aStep = aMeshfitVoxel.ComputeStep(aTargetPoint,aTargetVec);
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        MeshfitTime += finish - start;
//        G4cout <<"Has step:\t"<<aStep<<"\t";

//        G4double aDumySafety= kInfinity;
//        G4bool   aDumyBool = false, bDumyBool = false, cDumyBool = false;
//        G4int    aDumyInt = -1;
//        G4ThreeVector aDumyVec ;
//        G4VPhysicalVolume * aDumyPV;

//        start = (double)clock() / CLOCKS_PER_SEC;
//        G4double bStep = aNormalNavigator.ComputeStep(aTargetPoint,aTargetVec,
//                                                        /*propose step */kInfinity,
//                                                        aDumySafety, aHistory,aDumyBool, aDumyVec, bDumyBool,cDumyBool ,&aDumyPV, aDumyInt );
//        finish = (double)clock() / CLOCKS_PER_SEC;
//        NormalNavigatorTime += finish - start;
//        G4cout <<bStep<<G4endl ;
//        continue;
//    }
//    G4cout <<"CPU Time: Meshfit-" << MeshfitTime << "\t NormalNavigator-"<<NormalNavigatorTime <<G4endl;

    // ### test of G4UMeshHex###
    //create a G4UMeshHex as a cubic with a=200cm
    G4ThreeVector p0 (-100*cm, -100*cm, -100*cm);
    G4ThreeVector p1 ( 100*cm, -100*cm, -100*cm);
    G4ThreeVector p2 ( 100*cm,  100*cm, -100*cm);
    G4ThreeVector p3 (-100*cm,  100*cm, -100*cm);
    G4ThreeVector p4 (-100*cm, -100*cm,  100*cm);
    G4ThreeVector p5 ( 100*cm, -100*cm,  100*cm);
    G4ThreeVector p6 ( 100*cm,  100*cm,  100*cm);
    G4ThreeVector p7 (-100*cm,  100*cm,  100*cm);
    G4UMeshHex * aHex = new G4UMeshHex("Hex", p0, p1, p2, p3, p4, p5, p6, p7);
    //create a counterpart: G4Box with the same size
    G4double aSize = 200*cm;
    G4Box * aBox = new G4Box("Box", aSize/2, aSize/2, aSize/2);
//    //1. test inside()
//    G4double aFactor = 1000.0; //mm
//    G4cout <<"####Testing inside()..."<<G4endl;
//    for (int i=0; i<100000; i++) {

//        G4double aRandX = (drand48() * 2 - 1) * aFactor *2;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor *2;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor *2;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);

//        EInside Position1 = aHex->Inside(aTargetPoint);
//        EInside Position2 = aBox->Inside(aTargetPoint);
//        if (Position1 == Position2)
//            G4cout <<"is Pass." <<G4endl;
//        else
//            G4cout <<"NOT Pass!!" <<G4endl;
//    }
//    //2. test surface normal
//    G4cout <<"####Testing SurfaceNormal()..."<<G4endl;
//    for (int i=0; i<100000; i++) {

//        G4double aRandX = (drand48() * 2 - 1) * aFactor *2;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor *2;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor *2;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);

//        G4ThreeVector Normal1 = aHex->SurfaceNormal(aTargetPoint);
//        G4ThreeVector Normal2 = aBox->SurfaceNormal(aTargetPoint);
//        if (Normal1 == Normal2)
//            G4cout <<"is Pass." <<G4endl;
//        else
//            G4cout <<"NOT Pass!!" <<G4endl;
//    }
//    //3. test DistanceToIn(p,v)
//    G4cout <<"####Testing DistanceToIn(p,v)..."<<G4endl;
//    for (int i=0; i<100000; i++) {
//        G4double aRandX = (drand48() * 2 - 1) * aFactor*1.2;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor*1.2;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor*1.2;
//        G4double aRandXVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandYVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZVec = (drand48() * 2 - 1) * aFactor*1.;
//        //exclude the points inside
//        if (aRandX>(-aFactor) && aRandX < aFactor && aRandY>(-aFactor) && aRandY < aFactor && aRandZ>(-aFactor) && aRandZ < aFactor )
//            continue;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4cout <<"Vector\t"<<aRandXVec<<"\t\t"<<aRandYVec<<"\t\t"<<aRandZVec<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);
//        G4ThreeVector aTargetVec = G4ThreeVector(aRandXVec, aRandYVec, aRandZVec).unit();

//        G4double aDistance1 = aHex->DistanceToIn(aTargetPoint, aTargetVec);
//        G4double aDistance2 = aBox->DistanceToIn(aTargetPoint, aTargetVec);
//        if (std::fabs(aDistance1 - aDistance2) < 1e-9)
//            G4cout <<"is Pass." <<aDistance1<<"\t\t"<<aDistance2<<G4endl;
//        else
//            G4cout <<"NOT Pass!!\t\t" <<aDistance1<<"\t\t"<<aDistance2 <<G4endl;
//    }

        //4. test DistanceToIn(p)
//        G4cout <<"####Testing DistanceToIn(p)..., The safety distance might be different."<<G4endl;
//        for (int i=0; i<100; i++) {
//            G4double aRandX = (drand48() * 2 - 1) * aFactor*1.2;
//            G4double aRandY = (drand48() * 2 - 1) * aFactor*1.2;
//            G4double aRandZ = (drand48() * 2 - 1) * aFactor*1.2;
//            //exclude the points inside
//            if (aRandX>(-aFactor) && aRandX < aFactor && aRandY>(-aFactor) && aRandY < aFactor && aRandZ>(-aFactor) && aRandZ < aFactor )
//                continue;
//            G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//            G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);

//            G4double aDistance1 = aHex->DistanceToIn(aTargetPoint);
//            G4double aDistance2 = aBox->DistanceToIn(aTargetPoint);
//            //because the calculation method is different, therefore the safety is also different
//            if (aDistance1 <= aDistance2) //safety distance
//                G4cout <<"is Pass." <<aDistance1<<"\t\t"<<aDistance2<<G4endl;
//            else
//                G4cout <<"NOT Pass!!\t\t" <<aDistance1<<"\t\t"<<aDistance2 <<G4endl;
//        }

//    //5. test DistanceToOut(p,v)
//    G4cout <<"####Testing DistanceToOut(p,v)..."<<G4endl;
//    for (int i=0; i<100000; i++) {
//        G4double aRandX = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandXVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandYVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZVec = (drand48() * 2 - 1) * aFactor*1.;
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4cout <<"Vector\t"<<aRandXVec<<"\t\t"<<aRandYVec<<"\t\t"<<aRandZVec<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);
//        G4ThreeVector aTargetVec = G4ThreeVector(aRandXVec, aRandYVec, aRandZVec).unit();

//        G4double aDistance1 = aHex->DistanceToOut(aTargetPoint, aTargetVec);
//        G4double aDistance2 = aBox->DistanceToOut(aTargetPoint, aTargetVec);
//        if (std::fabs(aDistance1 - aDistance2) < 1e-9)
//            G4cout <<"is Pass." <<aDistance1<<"\t\t"<<aDistance2<<G4endl;
//        else
//            G4cout <<"NOT Pass!!\t\t" <<aDistance1<<"\t\t"<<aDistance2 <<G4endl;
//    }

//    //6. test DistanceToOut(p)
//    G4cout <<"####Testing DistanceToOut(p)..."<<G4endl;
//    for (int i=0; i<100000; i++) {
//        G4double aRandX = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandY = (drand48() * 2 - 1) * aFactor*1.;
//        G4double aRandZ = (drand48() * 2 - 1) * aFactor*1.;
//        //exclude the points inside
//        G4cout <<"Point\t"<<aRandX<<"\t\t"<<aRandY<<"\t\t"<<aRandZ<<"\t\t" ;
//        G4ThreeVector aTargetPoint = G4ThreeVector(aRandX, aRandY, aRandZ);

//        G4double aDistance1 = aHex->DistanceToOut(aTargetPoint);
//        G4double aDistance2 = aBox->DistanceToOut(aTargetPoint);
//        //because the calculation method is different, therefore the safety is also different
//        if (std::fabs(aDistance1 - aDistance2) < 1e-9) //safety distance
//            G4cout <<"is Pass." <<aDistance1<<"\t\t"<<aDistance2<<G4endl;
//        else
//            G4cout <<"NOT Pass!!\t\t" <<aDistance1<<"\t\t"<<aDistance2 <<G4endl;
//    }

//    //7. test GetPointOnSurface
//    G4cout <<"####Testing GetPointOnSurface()..."<<G4endl;
//    for (int i=0; i<100000; i++) {
//       G4ThreeVector aPoint = aHex->GetPointOnSurface();
//       G4cout <<aPoint.x()<<","<<aPoint.y()<<","<<aPoint.z()<<G4endl ;

//    }



    return 0;
}
