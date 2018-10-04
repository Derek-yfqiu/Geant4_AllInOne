#include "G4HalfSpaceSolid.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceSphere.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceCylinder.hh"
#include "G4HalfSpaceCylinderOnAxis.hh"
#include "G4HalfSpaceCone.hh"
#include "G4HalfSpaceConeOnAxis.hh"
#include "G4HalfSpaceTorus.hh"
#include "G4HalfSpaceBoolean.hh"


#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4CutTubs.hh"
#include "G4Ellipsoid.hh"

#include <G4SystemOfUnits.hh>
#include <globals.hh>
#include <string>
#include <G4String.hh>
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include <math.h>
#include <time.h>
#include <fstream>
using namespace std;
const G4double EPS = 1e-10;
//const G4double PI = 3.141592626535;
ofstream outlog;
G4ThreeVector randPoint(G4ThreeVector aLower, G4ThreeVector aUpper);
G4ThreeVector randVec();
void    testInside(G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample );
void    testSurfaceNormal (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample );
void    testDistanceToIn (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample);
void    testSafetyToIn (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample);
void    testDistanceToOut (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample);
void    testSafetyToOut (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample);





G4ThreeVector randPoint(G4ThreeVector aLower, G4ThreeVector aUpper)
{
    G4ThreeVector aPoint = G4ThreeVector(G4RandFlat::shoot(aLower.x(), aUpper.x()),
                                         G4RandFlat::shoot(aLower.y(), aUpper.y()),
                                         G4RandFlat::shoot(aLower.z(), aUpper.z()));
    return aPoint;

}
G4ThreeVector randVec()
{
    G4ThreeVector aVector = G4ThreeVector(G4RandFlat::shoot(-10.0, 10.0),
                                          G4RandFlat::shoot(-10.0, 10.0),
                                          G4RandFlat::shoot(-10.0, 10.0));
    return aVector.unit();
}

void    testVolume(G4HalfSpaceSolid * mySolid, G4VSolid * refSolid )
{
    outlog << "######## Volume Test #########"<<endl;

//    G4double RefVolume = refSolid->GetCubicVolume();
    G4double myVolume = mySolid->GetCubicVolume();
//    G4double aTmpDouble = RefVolume - myVolume;
//    outlog<<RefVolume<<"\t"<<myVolume<<"\t"<</*(RefVolume-myVolume)/RefVolume*/ aTmpDouble<<"\t";
//    if (aTmpDouble <= EPS)
//        outlog <<"is Pass." <<endl;
//    else
//        outlog <<"NOT Pass!!  : " <<aTmpDouble<<endl;

}

void    testInside(G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample )
{
    outlog << "######## Inside Test #########"<<endl;
    //get the boundary box
    G4BoundingBox3D aBox = *mySolid->getBoundaryBox();
    aBox.Margin(aBox.GetSize()/2); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
        EInside Position1 = mySolid->Inside(aRandPoint);
//        EInside Position2 = refSolid->Inside(aRandPoint);
//        if (Position1 == Position2)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!" <<endl;
    }

}

void    testSurfaceNormal (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample )
{
    outlog << "######## SurfaceNormal Test #########"<<endl;
    for (int i=0; i<Sample; i++){
        G4ThreeVector aRandPoint = refSolid->GetPointOnSurface();
        G4ThreeVector myNormal = mySolid->SurfaceNormal(aRandPoint);
//        G4ThreeVector refNormal = refSolid ->SurfaceNormal(aRandPoint);
//        G4double aTmpDouble = (refNormal - myNormal).mag2();
//        outlog<< aTmpDouble<<"\t";
//        if (aTmpDouble <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!  " <<endl;
    }
}

void    testDistanceToIn (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample)
{
    outlog << "######## DistanceToIn Test #########"<<endl;
    G4BoundingBox3D aBox = *mySolid->getBoundaryBox();
    aBox.Margin(aBox.GetSize()); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kInside) {
//            i--;
//            continue;
//        }
        G4ThreeVector aRandVec = randVec();
//        G4double refDist = refSolid->DistanceToIn(aRandPoint, aRandVec);
//        if (refDist == kInfinity) {
//            i--;
//            continue;
//        }
        G4double myDist = mySolid->DistanceToIn(aRandPoint, aRandVec);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//           outlog <<"is Pass.\t" <<myDist<<"\t:\t"<<refDist<<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;

    }

}
void    testSafetyToIn (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample)
{
    outlog << "######## SafetyToIn Test #########"<<endl;
    G4BoundingBox3D aBox = *mySolid->getBoundaryBox();
    aBox.Margin(aBox.GetSize() * 4); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kInside) {
//            i--;
//            continue;
//        }
        G4double myDist = mySolid->DistanceToIn(aRandPoint);
//        G4double refDist = refSolid->DistanceToIn(aRandPoint);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        outlog <<"myDist\t" <<myDist<<"\t refDist \t"<<refDist<<"\t Diff \t"<<diff<<endl;
    }
}

void    testDistanceToOut (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample)
{
    outlog << "######## DistanceToOut Test #########"<<endl;
    G4BoundingBox3D aBox = *mySolid->getBoundaryBox();
//    aBox.Margin(aBox.size); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
        if (refSolid->Inside(aRandPoint) == kOutside) {
            i--;
            continue;
        }
        G4ThreeVector aRandVec = randVec();

//        G4double refDist = refSolid->DistanceToOut(aRandPoint, aRandVec);
//        if (refDist == kInfinity) {
//            i--;
//            continue;
//        }
        G4double myDist = mySolid->DistanceToOut(aRandPoint, aRandVec);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;

    }

}
void    testSafetyToOut (G4HalfSpaceSolid * mySolid, G4VSolid * refSolid, G4int Sample)
{
    outlog << "######## SafetyToOut Test #########"<<endl;
    G4BoundingBox3D aBox = *mySolid->getBoundaryBox();
//    aBox.Margin(aBox.size * 4); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kOutside) {
//           i--;
//            continue;
//        }
        G4double myDist = mySolid->DistanceToOut(aRandPoint);
//        G4double refDist = refSolid->DistanceToOut(aRandPoint);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;
    }
}

void    testBox()
{
    outlog << "//////////////////// Box Test /////////////////////"<<endl;

    //reference
    G4Box * refSolid = new G4Box("ref", 30, 40, 60);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisX, -30, 1));
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisX, 30, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisY, -40, 1));
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisY, 40, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisZ, -60, 1));
    aSurfaceList.push_back(new G4HalfSpacePlane(AxisZ, 60, -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-30,-40,-60),
                                                      G4ThreeVector(30,40,60),
                                                     /*576000*/0,43200);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testSphere()
{
    outlog << "//////////////////// Sphere Test /////////////////////"<<endl;
    //reference
    G4Orb * refSolid = new G4Orb("ref", 100);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceSphere(0,0,0, 100, -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-100,-100,-100),
                                                      G4ThreeVector(100,100,100),
                                                     /*4188790.204786*/0,125662.715188);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testCylinder()
{
    outlog << "//////////////////// Cylinder Test /////////////////////"<<endl;
    G4Tubs * refSolid = new G4Tubs("ref",0,15,20,0, 2*CLHEP::pi);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, -20, 1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 20, -1));
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 15, -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-15,-15,-20),
                                                      G4ThreeVector(15,15,20),
                                                     /*28274.333882*/ 0,5183.627878);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

//    aSurfaceList.pop_back();
//    aSurfaceList.push_back(new G4HalfSpaceCylinder (AxisZ, 15,0 , 0, -1));
//     mySolid = new G4HalfSpaceSolid("mysolid",
//                                                      aSurfaceList,
//                                                      G4ThreeVector(-15,-15,-20),
//                                                      G4ThreeVector(15,15,20),
//                                                     28274.333882,5183.627878);
//    testInside(mySolid, refSolid, 100000);
//    testSurfaceNormal(mySolid, refSolid, 100000);
//    testDistanceToIn(mySolid, refSolid, 100000);
//    testSafetyToIn(mySolid, refSolid, 100000);
//    testDistanceToOut(mySolid, refSolid, 100000);
//    testSafetyToOut(mySolid, refSolid, 100000);

}

void    testCone()
{
    outlog << "//////////////////// Cone Test /////////////////////"<<endl;
    G4Cons * refSolid = new G4Cons("ref", 0,0,0,25,40,0,2*CLHEP::pi);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 40, -1));
    aSurfaceList.push_back(new G4HalfSpaceConeOnAxis(AxisZ, 0.09765625, -40, -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-25,-25,-40),
                                                      G4ThreeVector(25,25,40),
                                                     /*52359.87756*/0,8546.331562);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

    aSurfaceList.pop_back();
//    aSurfaceList.push_back(new G4HalfSpaceCone(AxisZ, 0.09765625, 0,0,-40, -1));
//    mySolid = new G4HalfSpaceSolid("mysolid",
//                                   aSurfaceList,
//                                   G4ThreeVector(-25,-25,-40),
//                                   G4ThreeVector(25,25,40),
//                                   52359.87756, 8546.331562);
//    testInside(mySolid, refSolid, 100000);
//    testSurfaceNormal(mySolid, refSolid, 100000);
//    testDistanceToIn(mySolid, refSolid, 100000);
//    testSafetyToIn(mySolid, refSolid, 100000);
//    testDistanceToOut(mySolid, refSolid, 100000);
//    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testEllipsoid()
{
    outlog << "//////////////////// Ellipsoid Test /////////////////////"<<endl;
    G4Ellipsoid * refSolid = new G4Ellipsoid("ref", 10, 20, 50);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceQuadric(0.01,
                                                  0,
                                                  0,
                                                  0,
                                                  0.0025,
                                                  0,
                                                  0,
                                                  0.0004,
                                                  0,
                                                  -1,
                                                  -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-11,-21,-51),
                                                      G4ThreeVector(11,21,51),
                                                     /*41870*/0,7868);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

    aSurfaceList.pop_back();
//    aSurfaceList.push_back(new G4HalfSpaceCone(AxisZ, 0.09765625, 0,0,-40, -1));
//    mySolid = new G4HalfSpaceSolid("mysolid",
//                                   aSurfaceList,
//                                   G4ThreeVector(-25,-25,-40),
//                                   G4ThreeVector(25,25,40),
//                                   52359.87756, 8546.331562);
//    testInside(mySolid, refSolid, 100000);
//    testSurfaceNormal(mySolid, refSolid, 100000);
//    testDistanceToIn(mySolid, refSolid, 100000);
//    testSafetyToIn(mySolid, refSolid, 100000);
//    testDistanceToOut(mySolid, refSolid, 100000);
//    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testTorus()
{
    outlog << "//////////////////// Torus Test /////////////////////"<<endl;
    G4Torus * refSolid = new G4Torus ("ref", 0, 60, 200, 0, 2*CLHEP::pi);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceTorus(G4ThreeVector(0,0,0), G4ThreeVector(0,0,1), 200, 60, -1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-260,-260,-60),
                                                      G4ThreeVector(260,260,60),
                                                     /*14212230.337569*/0,473892.774465);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

}

void    testTrapezoid()
{
    outlog << "//////////////////// Trapezoid Test /////////////////////"<<endl;
    G4Trd * refSolid = new G4Trd("ref", 30, 10, 40, 15, 60);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, -60, 1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 60, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0.0000000, 0.978980419738,-0.203954254112,-26.921961543,1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0.0000000,-0.978980419738,-0.203954254112,-26.921961543,1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0.986393923832,-0.0000000,0.164398987305, 19.727878477,-1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0.986393923832,0.0000000,-0.164398987305, -19.727878477 ,1));
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-30,-40,-60),
                                                      G4ThreeVector(30,40,60),
                                                     /*284000*/0,28588.198104);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testCylinderSection()
{
    outlog << "//////////////////// CylinderSection Test /////////////////////"<<endl;
    G4Tubs * refSolid = new G4Tubs("ref", 10, 15, 20, 0, CLHEP::pi/2 );
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 10,  1));
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 15, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, -20, 1)); //pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 20, -1));//pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisX, 0,  1));//px
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-0,-0,-20),
                                                      G4ThreeVector(15,15,20),
                                                     /*3926.990817*/0,1767.14586775);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

}

void    testCutTube()
{
    outlog << "//////////////////// Cut-Tube Test /////////////////////"<<endl;
    G4CutTubs * refSolid = new G4CutTubs("ref", 12, 20, 30, 0,CLHEP::pi /*1.5*/,
                                         G4ThreeVector(0, -0.7020742, -0.7121038),
                                         G4ThreeVector(0.7020742, 0, 0.7121038) );
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 12,  1));
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 20, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0, -0.7020742, -0.7121038,21.363114 , -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (0.7020742, 0, 0.7121038, 21.363114  , -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-20,-0,-50),
                                                      G4ThreeVector(20,20,50),
                                                     /*16186.157104*/0,5865.697735);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);



}

void    testConeSection()
{
    outlog << "//////////////////// Cone Section Test /////////////////////"<<endl;
    G4Cons * refSolid = new G4Cons("ref", 5, 10, 20, 25, 40, 0 , CLHEP::pi);
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceConeOnAxis(AxisZ,0.03515625,  -93.33333333333333333, -1));
    aSurfaceList.push_back(new G4HalfSpaceConeOnAxis(AxisZ,0.03515625,  -66.66666666666666667,   1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, -40, 1)); //pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 40, -1));//pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisY, 0, 1));//py
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-25,-25,-40),
                                                      G4ThreeVector(25,25,40),
                                                     /*18849.555922*/0,8942.452377);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);

}

void    testTorusSection()
{
    outlog << "//////////////////// Torus Section Test /////////////////////"<<endl;
    G4Torus * refSolid = new G4Torus ("ref", 40, 60, 200, 0, 0.5*CLHEP::pi );
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceTorus(G4ThreeVector(0,0,0), G4ThreeVector(0,0,1), 200, 60, -1));
    aSurfaceList.push_back(new G4HalfSpaceTorus(G4ThreeVector(0,0,0), G4ThreeVector(0,0,1), 200, 40,  1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisX, 0,  1));//px
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(0,0,-60),
                                                      G4ThreeVector(260,260,60),
                                                     /*14212230.337569*/0,473892.774465);
    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000);
    testSafetyToIn(mySolid, refSolid, 100000);
    testDistanceToOut(mySolid, refSolid, 100000);
    testSafetyToOut(mySolid, refSolid, 100000);
}

void    testBooleanCylinderSection()
{
    outlog << "//////////////////// CylinderSection Test /////////////////////"<<endl;
    G4Tubs * refSolid = new G4Tubs("ref", 10, 15, 20, 0, CLHEP::pi/2 );
    //mysolid
    vector < G4HalfSpaceSurface * > aSurfaceList;
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 10,  1));
    aSurfaceList.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 15, -1));
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, -20, 1)); //pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisZ, 20, -1));//pz
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisX, 0,  1));//px
    aSurfaceList.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    G4HalfSpaceSolid * mySolid = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList,
                                                      G4ThreeVector(-0,-0,-20),
                                                      G4ThreeVector(15,15,20),
                                                     /*3926.990817*/0,1767.14586775);
    //my boolean solid
    vector < G4HalfSpaceSurface * > aSurfaceList1,  aSurfaceList2;
    aSurfaceList1.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 10, -1));
    aSurfaceList2.push_back(new G4HalfSpaceCylinderOnAxis (AxisZ, 15, -1));
    aSurfaceList1.push_back(new G4HalfSpacePlane (AxisZ, -20, 1)); //pz
    aSurfaceList1.push_back(new G4HalfSpacePlane (AxisZ, 20, -1));//pz
    aSurfaceList1.push_back(new G4HalfSpacePlane (AxisX, 0,  1));//px
    aSurfaceList1.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    aSurfaceList2.push_back(new G4HalfSpacePlane (AxisZ, -20, 1)); //pz
    aSurfaceList2.push_back(new G4HalfSpacePlane (AxisZ, 20, -1));//pz
    aSurfaceList2.push_back(new G4HalfSpacePlane (AxisX, 0,  1));//px
    aSurfaceList2.push_back(new G4HalfSpacePlane (AxisY, 0,  1));//py
    G4HalfSpaceSolid * mySolid1 = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList1,
                                                      G4ThreeVector(-0,-0,-20),
                                                      G4ThreeVector(15,15,20),
                                                     /*3926.990817*/0,1767.14586775);
    G4HalfSpaceSolid * mySolid2 = new G4HalfSpaceSolid("mysolid",
                                                      aSurfaceList2,
                                                      G4ThreeVector(-0,-0,-20),
                                                      G4ThreeVector(15,15,20),
                                                     /*3926.990817*/0,1767.14586775);
    G4VSolid * boolSolid = G4HalfSpaceBoolean::Subtract(mySolid2, mySolid1);

    testVolume(mySolid, boolSolid);
    testInside(mySolid, boolSolid, 100000);
    testSurfaceNormal(mySolid, boolSolid, 100000);
    testDistanceToIn(mySolid, boolSolid, 100000);
    testSafetyToIn(mySolid, boolSolid, 100000);
    testDistanceToOut(mySolid, boolSolid, 100000);
    testSafetyToOut(mySolid, boolSolid, 100000);

}

int main (int argc, char **argv) {

    outlog.open("out.log");
    double start, finish;
    start = (double)clock() / CLOCKS_PER_SEC;
    testBox();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"box: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testSphere();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Sphere: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testCylinder();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Cylinder: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testCone();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Cone: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testTorus();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Torus: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testTrapezoid();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Trapezoid: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testCylinderSection();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"CylinderSection: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testCutTube();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"CutTube: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testConeSection();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"ConeSection: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testTorusSection();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"TorusSection: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testBooleanCylinderSection();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"BooleanCylinderSection: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testEllipsoid();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Ellipsoid: "<<finish - start <<endl;


    return 0;
}
