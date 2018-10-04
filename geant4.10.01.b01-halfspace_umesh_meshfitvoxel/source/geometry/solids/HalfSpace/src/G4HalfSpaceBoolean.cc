#include "G4HalfSpaceBoolean.hh"
namespace G4HalfSpaceBoolean
{

G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1)
{
    return new G4IntersectionSolid("Unnamed",A,B1);
}

G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2)
{
    return Intersect( Intersect(A, B1), B2);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3)
{
    return Intersect( Intersect(A, B1, B2), B3);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4)
{
    return Intersect( Intersect(A, B1, B2, B3), B4);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5)
{
    return Intersect( Intersect(A, B1, B2, B3, B4), B5);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6)
{
    return Intersect( Intersect(A, B1, B2, B3, B4, B5), B6);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7)
{
    return Intersect( Intersect(A, B1, B2, B3, B4, B5, B6), B7);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8)
{
    return Intersect( Intersect(A, B1, B2, B3, B4, B5, B6, B7), B8);
}
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8,
                      G4VSolid * B9)
{
    return Intersect( Intersect(A, B1, B2, B3, B4, B5, B6, B7, B8), B9);
}
G4VSolid * Intersect( std::vector <G4VSolid*> aSolidList)
{
    if (aSolidList.size() < 2) {
        G4Exception("G4HalfSpaceBoolean::Intersect","GeomSolids1003",
                    FatalException, "Should have at lease to two solids!" );
    }
    G4VSolid * aSolid = Intersect(aSolidList[0], aSolidList[1]);
    for (unsigned int i=2; i<aSolidList.size(); i++){
        aSolid = Intersect(aSolid, aSolidList[i]);
    }
    return aSolid;
}


/////////////////////////////////
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1)
{
    return new G4UnionSolid("Unnamed",A,B1);
}

G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2)
{
    return Union( Union(A, B1), B2);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3)
{
    return Union( Union(A, B1, B2), B3);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4)
{
    return Union( Union(A, B1, B2, B3), B4);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4,
                  G4VSolid * B5)
{
    return Union( Union(A, B1, B2, B3, B4), B5);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4,
                  G4VSolid * B5,
                  G4VSolid * B6)
{
    return Union( Union(A, B1, B2, B3, B4, B5), B6);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4,
                  G4VSolid * B5,
                  G4VSolid * B6,
                  G4VSolid * B7)
{
    return Union( Union(A, B1, B2, B3, B4, B5, B6), B7);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4,
                  G4VSolid * B5,
                  G4VSolid * B6,
                  G4VSolid * B7,
                  G4VSolid * B8)
{
    return Union( Union(A, B1, B2, B3, B4, B5, B6, B7), B8);
}
G4VSolid * Union( G4VSolid * A,
                  G4VSolid * B1,
                  G4VSolid * B2,
                  G4VSolid * B3,
                  G4VSolid * B4,
                  G4VSolid * B5,
                  G4VSolid * B6,
                  G4VSolid * B7,
                  G4VSolid * B8,
                  G4VSolid * B9)
{
    return Union( Union(A, B1, B2, B3, B4, B5, B6, B7, B8), B9);
}
G4VSolid * Union( std::vector <G4VSolid*> aSolidList)
{
    if (aSolidList.size() < 2) {
        G4Exception("G4HalfSpaceBoolean::Union","GeomSolids1003",
                    FatalException, "Should have at lease to two solids!" );
    }
    G4VSolid * aSolid = Union(aSolidList[0], aSolidList[1]);
    for (unsigned int i=2; i<aSolidList.size(); i++){
        aSolid = Union(aSolid, aSolidList[i]);
    }
    return aSolid;
}

////////////////////////////////////////

G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1)
{
    return new G4SubtractionSolid("Unnamed",A,B1);
}

G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2)
{
    return Subtract( Subtract(A, B1), B2);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3)
{
    return Subtract( Subtract(A, B1, B2), B3);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4)
{
    return Subtract( Subtract(A, B1, B2, B3), B4);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4,
                     G4VSolid * B5)
{
    return Subtract( Subtract(A, B1, B2, B3, B4), B5);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4,
                     G4VSolid * B5,
                     G4VSolid * B6)
{
    return Subtract( Subtract(A, B1, B2, B3, B4, B5), B6);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4,
                     G4VSolid * B5,
                     G4VSolid * B6,
                     G4VSolid * B7)
{
    return Subtract( Subtract(A, B1, B2, B3, B4, B5, B6), B7);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4,
                     G4VSolid * B5,
                     G4VSolid * B6,
                     G4VSolid * B7,
                     G4VSolid * B8)
{
    return Subtract( Subtract(A, B1, B2, B3, B4, B5, B6, B7), B8);
}
G4VSolid * Subtract( G4VSolid * A,
                     G4VSolid * B1,
                     G4VSolid * B2,
                     G4VSolid * B3,
                     G4VSolid * B4,
                     G4VSolid * B5,
                     G4VSolid * B6,
                     G4VSolid * B7,
                     G4VSolid * B8,
                     G4VSolid * B9)
{
    return Subtract( Subtract(A, B1, B2, B3, B4, B5, B6, B7, B8), B9);
}
G4VSolid * Subtract( std::vector <G4VSolid*> aSolidList)
{
    if (aSolidList.size() < 2) {
        G4Exception("G4HalfSpaceBoolean::Subtract","GeomSolids1003",
                    FatalException, "Should have at lease to two solids!" );
    }
    G4VSolid * aSolid = Subtract(aSolidList[0], aSolidList[1]);
    for (unsigned int i=2; i<aSolidList.size(); i++){
        aSolid = Subtract(aSolid, aSolidList[i]);
    }
    return aSolid;
}

}//end of namespace

