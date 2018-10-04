#ifndef G4HALFSPACEBOOLEAN_HH
#define G4HALFSPACEBOOLEAN_HH
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include <vector>

/*!
 *  Implementation of methods for Boolean operation
 *  Intersect: A and B
 *  Union: A or B
 *  Subtract: A not B
 */
namespace G4HalfSpaceBoolean
{

//implement in this way so that the method is easy for use
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8);
G4VSolid * Intersect( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8,
                      G4VSolid * B9);
G4VSolid * Intersect( std::vector <G4VSolid*> aSolidList);


G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8);
G4VSolid * Union( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8,
                      G4VSolid * B9);
G4VSolid * Union( std::vector <G4VSolid*> aSolidList);


G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8);
G4VSolid * Subtract( G4VSolid * A,
                      G4VSolid * B1,
                      G4VSolid * B2,
                      G4VSolid * B3,
                      G4VSolid * B4,
                      G4VSolid * B5,
                      G4VSolid * B6,
                      G4VSolid * B7,
                      G4VSolid * B8,
                      G4VSolid * B9);
G4VSolid * Subtract( std::vector <G4VSolid*> aSolidList);
}


#endif // G4HALFSPACEBOOLEAN_HH
