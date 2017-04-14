project(SDK_CommandPluginExample)

 #cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 2.8)

###############################################################################
# Options for building the plugin
###############################################################################

option(BUILD_DAGMC_EXPORTER  "Include DagMC export capability."     ON )
option(BUILD_MCNP_IMPORTER   "Include MCNP import capability."      ON )
option(BUILD_IGEOM_TESTS     "Include tests for iGeom."             OFF)


###############################################################################
# Compilation flags
###############################################################################


set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

# These are used in MyPlugin.cpp to keep the plugin working 
# if not all parts are built.
if(${BUILD_DAGMC_EXPORTER})
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_DAGMC_EXPORT" )
endif()
if(${BUILD_MCNP_IMPORTER})
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_MCNP_IMPORT" )
endif()
if(${BUILD_IGEOM_TESTS})
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_IGEOM_TESTS" )
endif()

# adjust compiler setting for Linux using gcc version 5.0 and higher
if(CMAKE_SYSTEM_NAME MATCHES Linux AND CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
    if(NOT CMAKE_CXX_FLAGS MATCHES _GLIBCXX_USE_CXX11_ABI)
      set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 ${CMAKE_CXX_FLAGS}")
    endif()
  endif()
endif()




###############################################################################
# Find necessary packages, libraries, and include directories
###############################################################################

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIRS})


# Find the other libraries and directories
include_directories("${DAGMC_DIR}/../include/")
include_directories("${CMAKE_SOURCE_DIR}/igeom")

if(Cubit_DIR)
    find_library(CUBITI19 NAMES cubiti19 PATHS ${Cubit_DIR})
else()
    find_library(CUBITI19 NAMES cubiti19 PATHS ${CMAKE_PREFIX_PATH})
endif()

if(${BUILD_MCNP_IMPORTER})
    if(MCNP2CAD_DIR)
        find_library(MCNP2CAD NAMES mcnp2cad PATHS ${MCNP2CAD_DIR})
    elseif(Cubit_DIR)
        find_library(MCNP2CAD NAMES mcnp2cad PATHS ${Cubit_DIR})
    else()
        find_library(MCNP2CAD NAMES mcnp2cad PATHS ${CMAKE_PREFIX_PATH})
    endif()
endif()

include_directories(${MCNP2CAD_DIR}/..)

if(${BUILD_DAGMC_EXPORTER})
    find_library(MAKE_WATERTIGHT NAMES makeWatertight PATHS ${DAGMC_DIR})
endif()

###############################################################################
# Source variables
###############################################################################

set(PLUGIN_SRC
    MyPlugin.cpp
    MyPlugin.hpp)
    
if(${BUILD_DAGMC_EXPORTER})
    set(DAGMC_SRC
        DAGMCExportCommand.cpp
        DAGMCExportCommand.hpp)
endif()

if(${BUILD_MCNP_IMPORTER})
    set(MCNP_SRC
        MCNPImp.cpp
        MCNPImp.hpp)
endif()

if(${BUILD_IGEOM_TESTS})
    set(IGEOMTEST_SRC
        iGeom_test.cpp
        iGeom_test.hpp)
endif()

SET(IGEOM_SRC
    igeom/iGeom.cpp
    igeom/iGeom.h
    igeom/iGeomError.cc)

###############################################################################
# Build commands
###############################################################################

add_library(iGeom SHARED ${IGEOM_SRC})
target_link_libraries(iGeom cubit_geom ${CUBITI19})

add_library(dagmc_export_plugin MODULE ${PLUGIN_SRC} ${DAGMC_SRC} ${MCNP_SRC} ${IGEOMTEST_SRC})
target_link_libraries(dagmc_export_plugin cubiti cubit_util cubit_geom ${MAKE_WATERTIGHT} ${MOAB_LIBRARIES} ${MCNP2CAD})
