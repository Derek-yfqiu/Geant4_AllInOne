#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4UMesh
# Package: Geant4.src.G4UMesh.G4Parser
#
# CMakeLists.txt for building a single granular library.
#
# Created on : 24.02.2014
#
# 
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()


# executable for testing
# List external includes needed.
#include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
#include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/usolids/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

#set(TestMain_Parser_SRCS
#    Test_Main_Parser.cc
#)
#SET(CMAKE_BUILD_TYPE Debug)
#include_directories(./include)
#LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/outputs/library/Linux-g++)
#add_executable(TestMain_Parser ${TestMain_Parser_SRCS})
#target_link_libraries(TestMain_Parser G4UMesh )
