#------------------------------------------------------------------------------
# sources.cmake
# Module : G4UMeshParser
# Package: Geant4.src.G4UMesh.G4UMeshParser
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 15/07/2013
#
# $Id: sources.cmake 72956 2013-08-14 14:24:35Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/usolids/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)



#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4UMeshParser
    HEADERS
        G4UMeshReader.hh
        G4UMeshVTKLegacyReader.hh
        G4UMeshGeneral.hh
        G4UMeshParserMessenger.hh
        G4UMeshParser.hh
        G4UMeshWriter.hh
        G4UMeshVTKLegacyWriter.hh
    SOURCES
        G4UMeshReader.cc
        G4UMeshVTKLegacyReader.cc
        G4UMeshParserMessenger.cc
        G4UMeshParser.cc
        G4UMeshWriter.cc
        G4UMeshVTKLegacyWriter.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4geometrymng
        G4volumes
    GLOBAL_DEPENDENCIES
        G4global
        G4geometry
    LINK_LIBRARIES
)

# List any source specific properties here


