#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geomHalfSpace
# Package: Geant4.src.G4geometry..G4geomHalfSpace
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
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
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/HalfSpace/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/usolids/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)

set(TestHalfSpace_SRCS
    TestHalfSpace.cc
)
# SET(CMAKE_BUILD_TYPE Debug)
#LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/outputs/library/Linux-g++)
#add_executable(TestHalfSpace ${TestHalfSpace_SRCS})
#target_link_libraries(TestHalfSpace G4geometry )

#set(TestUMesh_SRCS
#    TestUMesh.cc
#)
#add_executable(TestUMesh ${TestUMesh_SRCS})
#target_link_libraries(TestUMesh G4geometry )
