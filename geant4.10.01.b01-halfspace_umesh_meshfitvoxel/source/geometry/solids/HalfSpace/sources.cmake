#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomHalfSpace
# Package: Geant4.src.G4geometry..G4geomHalfSpace
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 76607 2013-11-13 08:46:40Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/usolids/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomHalfSpace
    HEADERS
        G4HalfSpaceSurface.hh
        G4HalfSpacePlane.hh
        G4HalfSpaceSolid.hh
        G4BoundingBox3D.hh
        G4BoundingBox3D.icc
        G4HalfSpaceSphere.hh
        G4HalfSpaceQuadric.hh
        G4HalfSpaceCylinder.hh
        G4HalfSpaceCylinderOnAxis.hh
        G4HalfSpaceCone.hh
        G4HalfSpaceConeOnAxis.hh
        G4HalfSpaceTorus.hh
        G4HalfSpaceBoolean.hh
#        G4HalfSpaceEllipticalTorus.hh
#        G4HalfSpaceSpecialQuadric.hh
#        G4Plane.hh
#        G4Ray.hh
#        G4Ray.icc
#        G4Sort.hh
    SOURCES
        G4HalfSpaceSurface.cc
        G4HalfSpacePlane.cc
        G4HalfSpaceSolid.cc
        G4BoundingBox3D.cc
        G4HalfSpaceSphere.cc
        G4HalfSpaceQuadric.cc
        G4HalfSpaceCylinder.cc
        G4HalfSpaceCylinderOnAxis.cc
        G4HalfSpaceCone.cc
        G4HalfSpaceConeOnAxis.cc
        G4HalfSpaceTorus.cc
        G4HalfSpaceBoolean.cc
#        G4HalfSpaceEllipticalTorus.cc
#        G4HalfSpaceSpecialQuadric.cc
#        G4Ray.cc
#        G4Sort.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4graphics_reps
        G4volumes
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
    LINK_LIBRARIES
)

# List any source specific properties here

