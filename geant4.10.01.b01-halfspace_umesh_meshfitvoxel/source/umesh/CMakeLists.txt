#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4UMesh
# Package: Geant4.src.G4UMesh
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 75578 2013-11-04 12:03:33Z gcosmo $
#
#------------------------------------------------------------------------------
# SET(CMAKE_BUILD_TYPE Debug)
add_subdirectory(UMeshParser)
add_subdirectory(MeshfitVoxel)
#add_subdirectory(Navigator)
#add_subdirectory(PhysicalVolume)
#add_subdirectory(Scoring)



if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4UMesh
        COMPONENTS
            UMeshParser/sources.cmake
            MeshfitVoxel/sources.cmake
#            Navigator/sources.cmake
#            PhysicalVolume/sources.cmake
#            Scoring/sources.cmake
)
endif()


# executable for testing
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

set(TestMain_SRCS
    Test_Main.cc
)
#SET(CMAKE_BUILD_TYPE Debug)
#include_directories(./UMeshParser/include)
#include_directories(./MeshfitVoxel/include)
#LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/outputs/library/Linux-g++)
#add_executable(TestMain ${TestMain_SRCS})
#target_link_libraries(TestMain G4UMesh )

