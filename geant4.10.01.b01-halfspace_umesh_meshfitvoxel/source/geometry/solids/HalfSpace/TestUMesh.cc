#include "G4BoundingBox3D.hh"
#include "G4UMeshElement1st.hh"
#include "G4UMeshHex.hh"
#include "G4UMeshTet.hh"
#include "G4UMeshPyrm.hh"
#include "G4UMeshPent.hh"


#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4CutTubs.hh"
#include "G4Tet.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"

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
void    testInside(G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample );
void    testSurfaceNormal (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample );
void    testDistanceToIn (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample);
void    testSafetyToIn (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample);
void    testDistanceToOut (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample);
void    testSafetyToOut (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample);





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

void    testVolume(G4UMeshElement1st * mySolid, G4VSolid * refSolid )
{
    outlog << "######## Volume Test #########"<<endl;

    G4double RefVolume = refSolid->GetCubicVolume();
//    G4double myVolume = mySolid->GetCubicVolume();
//    G4double aTmpDouble = RefVolume - myVolume;
//    outlog<<RefVolume<<"\t"<<myVolume<<"\t"<</*(RefVolume-myVolume)/RefVolume*/ aTmpDouble<<"\t";
//    if (aTmpDouble <= EPS)
//        outlog <<"is Pass." <<endl;
//    else
//        outlog <<"NOT Pass!!  : " <<aTmpDouble<<endl;

}

void    testInside(G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample, G4BoundingBox3D & BBox )
{
    outlog << "######## Inside Test #########"<<endl;
    //get the boundary box
    G4BoundingBox3D aBox = BBox;
    aBox.Margin(aBox.GetSize()/2); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
 //       EInside Position1 = mySolid->Inside(aRandPoint);
        EInside Position2 = refSolid->Inside(aRandPoint);
//        if (Position1 == Position2)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!" <<endl;
    }

}

void    testSurfaceNormal (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample )
{
    outlog << "######## SurfaceNormal Test #########"<<endl;
    for (int i=0; i<Sample; i++){
        G4ThreeVector aRandPoint = refSolid->GetPointOnSurface();
 //       G4ThreeVector myNormal = mySolid->SurfaceNormal(aRandPoint);
        G4ThreeVector refNormal = refSolid ->SurfaceNormal(aRandPoint);
//        G4double aTmpDouble = (refNormal - myNormal).mag2();
//        outlog<< aTmpDouble<<"\t";
//        if (aTmpDouble <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!  " <<endl;
    }
}

void    testDistanceToIn (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample, G4BoundingBox3D & BBox)
{
    outlog << "######## DistanceToIn Test #########"<<endl;
    G4BoundingBox3D aBox = BBox;
    aBox.Margin(aBox.GetSize()); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kInside) {
//            i--;
//            continue;
//        }
        G4ThreeVector aRandVec = randVec();
        G4double refDist = refSolid->DistanceToIn(aRandPoint, aRandVec);
//        if (refDist == kInfinity) {
//            i--;
//            continue;
//        }
//        G4double myDist = mySolid->DistanceToIn(aRandPoint, aRandVec);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//            outlog <<"is Pass.\t" <<myDist<<"\t:\t"<<refDist<<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;

    }

}
void    testSafetyToIn (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample, G4BoundingBox3D & BBox)
{
    outlog << "######## SafetyToIn Test #########"<<endl;
    G4BoundingBox3D aBox = BBox;
    aBox.Margin(aBox.GetSize() * 4); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kInside) {
//            i--;
//            continue;
//        }
 //       G4double myDist = mySolid->DistanceToIn(aRandPoint);
        G4double refDist = refSolid->DistanceToIn(aRandPoint);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        outlog <<"myDist\t" <<myDist<<"\t refDist \t"<<refDist<<"\t Diff \t"<<diff<<endl;
    }
}

void    testDistanceToOut (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample,G4BoundingBox3D & BBox)
{
    outlog << "######## DistanceToOut Test #########"<<endl;
    G4BoundingBox3D aBox = BBox;
//    aBox.Margin(aBox.size); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
        if (refSolid->Inside(aRandPoint) == kOutside) {
            i--;
            continue;
        }
        G4ThreeVector aRandVec = randVec();

        G4double refDist = refSolid->DistanceToOut(aRandPoint, aRandVec);
//        if (refDist == kInfinity) {
//            i--;
//            continue;
//        }
 //       G4double myDist = mySolid->DistanceToOut(aRandPoint, aRandVec);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;

    }

}
void    testSafetyToOut (G4UMeshElement1st * mySolid, G4VSolid * refSolid, G4int Sample,G4BoundingBox3D & BBox)
{
    outlog << "######## SafetyToOut Test #########"<<endl;
    G4BoundingBox3D aBox = BBox;
//    aBox.Margin(aBox.size * 4); //enlarge the boundary box
    for (int i=0; i<Sample; i++)
    {
        G4ThreeVector aRandPoint = randPoint(aBox.GetBoxMin(), aBox.GetBoxMax());
//        if (refSolid->Inside(aRandPoint) == kOutside) {
//            i--;
//            continue;
//        }
//        G4double myDist = mySolid->DistanceToOut(aRandPoint);
        G4double refDist = refSolid->DistanceToOut(aRandPoint);
//        G4double diff = /*fabs*/(refDist-myDist);
//        outlog<<diff<<"\t";
//        if (fabs(myDist - refDist) <= EPS)
//            outlog <<"is Pass." <<endl;
//        else
//            outlog <<"NOT Pass!!   " <<myDist<<"\t:\t"<<refDist<<"\t Diff \t"<<diff<<endl;
    }
}


void testUMeshBox()
{
    G4ThreeVector p0 (-100*mm, -100*mm, -100*mm);
    G4ThreeVector p1 ( 100*mm, -100*mm, -100*mm);
    G4ThreeVector p2 ( 100*mm,  100*mm, -100*mm);
    G4ThreeVector p3 (-100*mm,  100*mm, -100*mm);
    G4ThreeVector p4 (-100*mm, -100*mm,  100*mm);
    G4ThreeVector p5 ( 100*mm, -100*mm,  100*mm);
    G4ThreeVector p6 ( 100*mm,  100*mm,  100*mm);
    G4ThreeVector p7 (-100*mm,  100*mm,  100*mm);
    G4UMeshHex * mySolid = new G4UMeshHex("Hex", p0, p1, p2, p3, p4, p5, p6, p7);
    //create a counterpart: G4Box with the same size
    G4double aSize = 200*mm;
    G4Box * refSolid = new G4Box("Box", aSize/2, aSize/2, aSize/2);
    G4BoundingBox3D aBox  = G4BoundingBox3D(G4ThreeVector(-101,-101,-101),G4ThreeVector(101,101,101));

    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000,aBox);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000,aBox);
    testSafetyToIn(mySolid, refSolid, 100000,aBox);
    testDistanceToOut(mySolid, refSolid, 100000,aBox);
    testSafetyToOut(mySolid, refSolid, 100000,aBox);

}

void testUMeshTrapozoid()
{
    outlog << "//////////////////// Trapezoid Test /////////////////////"<<endl;
    G4Trd * refSolid = new G4Trd("ref", 30, 10, 40, 15, 60);
    //mysolid

    G4ThreeVector p0 (-30*mm, -40*mm, -60*mm);
    G4ThreeVector p1 ( 30*mm, -40*mm, -60*mm);
    G4ThreeVector p2 ( 30*mm,  40*mm, -60*mm);
    G4ThreeVector p3 (-30*mm,  40*mm, -60*mm);

    G4ThreeVector p4 (-10*mm, -15*mm,  60*mm);
    G4ThreeVector p5 ( 10*mm, -15*mm,  60*mm);
    G4ThreeVector p6 ( 10*mm,  15*mm,  60*mm);
    G4ThreeVector p7 (-10*mm,  15*mm,  60*mm);

    G4UMeshHex * mySolid = new G4UMeshHex("Hex", p0, p1, p2, p3, p4, p5, p6, p7);
    G4BoundingBox3D aBox  = G4BoundingBox3D(G4ThreeVector(-31,-41,-61),G4ThreeVector(31, 41, 61));

    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000,aBox);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000,aBox);
    testSafetyToIn(mySolid, refSolid, 100000,aBox);
    testDistanceToOut(mySolid, refSolid, 100000,aBox);
    testSafetyToOut(mySolid, refSolid, 100000,aBox);

}

void testUMeshPyramid()
{
    outlog << "//////////////////// Pyramid Test /////////////////////"<<endl;
    G4Trd * refSolid = new G4Trd("ref", 30, 0, 40, 0, 60);
    //mysolid

    G4ThreeVector p0 (-30*mm, -40*mm, -60*mm);
    G4ThreeVector p1 ( 30*mm, -40*mm, -60*mm);
    G4ThreeVector p2 ( 30*mm,  40*mm, -60*mm);
    G4ThreeVector p3 (-30*mm,  40*mm, -60*mm);

    G4ThreeVector p4 (0*mm, 0*mm,  60*mm);
    G4UMeshPyrm * mySolid = new G4UMeshPyrm("pyr", p0, p1, p2, p3, p4);
    G4BoundingBox3D aBox  = G4BoundingBox3D(G4ThreeVector(-31,-41,-61),G4ThreeVector(31, 41, 61));

    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000,aBox);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000,aBox);
    testSafetyToIn(mySolid, refSolid, 100000,aBox);
    testDistanceToOut(mySolid, refSolid, 100000,aBox);
    testSafetyToOut(mySolid, refSolid, 100000,aBox);

}

void testUMeshTet()
{
    outlog << "//////////////////// tetrahedron Test /////////////////////"<<endl;
    //mysolid

    G4ThreeVector p0 (-30*mm, -40*mm, -60*mm);
    G4ThreeVector p1 ( 30*mm, -40*mm, -60*mm);
    G4ThreeVector p2 ( 30*mm,  40*mm, -60*mm);

    G4ThreeVector p3 (0*mm, 0*mm,  60*mm);
    G4Tet* refSolid = new G4Tet("ref", p0, p1, p2, p3);

    G4UMeshTet * mySolid = new G4UMeshTet("tet", p0, p1, p2, p3);
    G4BoundingBox3D aBox  = G4BoundingBox3D(G4ThreeVector(-31,-41,-61),G4ThreeVector(31, 41, 61));

    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000,aBox);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000,aBox);
    testSafetyToIn(mySolid, refSolid, 100000,aBox);
    testDistanceToOut(mySolid, refSolid, 100000,aBox);
    testSafetyToOut(mySolid, refSolid, 100000,aBox);

}

void testUMeshPent()
{
    outlog << "//////////////////// pentahedron Test /////////////////////"<<endl;
    //mysolid

    G4ThreeVector p0 (-30*mm, -40*mm, -60*mm);
    G4ThreeVector p1 ( 30*mm, -40*mm, -60*mm);
    G4ThreeVector p2 ( 30*mm,  40*mm, -60*mm);
    G4ThreeVector p3 (-30*mm, -40*mm,  60*mm);
    G4ThreeVector p4 ( 30*mm, -40*mm,  60*mm);
    G4ThreeVector p5 ( 30*mm,  40*mm,  60*mm);


    G4UMeshPent * mySolid = new G4UMeshPent("tet", p0, p1, p2, p3, p4, p5);
    G4BoundingBox3D aBox  = G4BoundingBox3D(G4ThreeVector(-31,-41,-61),G4ThreeVector(31, 41, 61));

    std::vector <G4TwoVector> aPointList;
    aPointList.push_back(G4TwoVector  (-30*mm, -40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm,  40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm, -40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm, -40*mm));

    aPointList.push_back(G4TwoVector  (-30*mm, -40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm,  40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm, -40*mm));
    aPointList.push_back(G4TwoVector  ( 30*mm, -40*mm));

    G4GenericTrap* refSolid  = new G4GenericTrap("GenTra",60,aPointList);

    testVolume(mySolid,refSolid);
    testInside(mySolid, refSolid, 100000,aBox);
    testSurfaceNormal(mySolid, refSolid, 100000);
    testDistanceToIn(mySolid, refSolid, 100000,aBox);
    testSafetyToIn(mySolid, refSolid, 100000,aBox);
    testDistanceToOut(mySolid, refSolid, 100000,aBox);
    testSafetyToOut(mySolid, refSolid, 100000,aBox);

}

int main (int argc, char **argv) {

    outlog.open("out.log");
    double start, finish;
    start = (double)clock() / CLOCKS_PER_SEC;
    testUMeshBox();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Box: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testUMeshTrapozoid();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"rapozoid: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testUMeshPyramid();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Pyramid: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testUMeshTet();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Tet: "<<finish - start <<endl;
    start = (double)clock() / CLOCKS_PER_SEC;
    testUMeshPent();
    finish = (double)clock() / CLOCKS_PER_SEC;
    outlog <<"Pent: "<<finish - start <<endl;

    return 0;
}
