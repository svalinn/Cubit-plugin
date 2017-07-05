project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 2.8)

###############################################################################
# Find necessary packages, libraries, and include directories
###############################################################################

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

# adjust compiler setting for Linux using gcc version 5.0 and higher
if(CMAKE_SYSTEM_NAME MATCHES Linux AND CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
    if(NOT CMAKE_CXX_FLAGS MATCHES _GLIBCXX_USE_CXX11_ABI)
      set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 ${CMAKE_CXX_FLAGS}")
    endif()
  endif()
endif()

if(Cubit_DIR)
    find_library(CUBITI19 NAMES cubiti19 PATHS ${Cubit_DIR})
else()
    find_library(CUBITI19 NAMES cubiti19 PATHS ${CMAKE_PREFIX_PATH})
endif()

###############################################################################
# Options for building the plugin
###############################################################################

option(BUILD_DAGMC_EXPORTER  "Include DagMC export capability."     ON )
option(BUILD_MCNP_IMPORTER   "Include MCNP import capability."      ON )
option(BUILD_IGEOM           "Include iGeom."                       OFF)
option(BUILD_IGEOM_TESTS     "Include tests for iGeom."             OFF)

###############################################################################
# Compilation flags
###############################################################################

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

###############################################################################
# Top-level plugin files
###############################################################################

set(PLUGIN_SRC
    SvalinnPlugin.cpp
    SvalinnPlugin.hpp)

  
###############################################################################
# DAGMC Exporter
###############################################################################
if(${BUILD_DAGMC_EXPORTER})
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_DAGMC_EXPORT" )

    find_package(MOAB REQUIRED)
    include_directories(${MOAB_INCLUDE_DIRS})

    # Find the other libraries and directories
    include_directories("${DAGMC_DIR}/include")

    find_library(MAKE_WATERTIGHT NAMES makeWatertight PATHS "${DAGMC_DIR}/lib" )

    set(DAGMC_SRC
        export_dagmc_cmd/DAGMCExportCommand.cpp
        export_dagmc_cmd/DAGMCExportCommand.hpp)

      
endif()

###############################################################################
# MCNP Importer
###############################################################################
if(${BUILD_MCNP_IMPORTER})
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_MCNP_IMPORT" )
    set( BUILD_IGEOM ON)

    include_directories("${CMAKE_SOURCE_DIR}/iGeom")

    add_subdirectory(mcnp2cad)
    include_directories(mcnp2cad)

    set( MCNP2CAD_LIB "mcnp2cad")
    set(MCNP_SRC
        import_mcnp_cmd/MCNPImp.cpp
        import_mcnp_cmd/MCNPImp.hpp)


endif()


###############################################################################
# iGeom Tests
###############################################################################

if(${BUILD_IGEOM_TESTS})
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_IGEOM_TESTS" )
  set( BUILD_IGEOM ON)

    include_directories("${CMAKE_SOURCE_DIR}/iGeom")

    set(IGEOMTEST_SRC
      iGeom/tests/iGeom_test.cpp
      iGeom/tests/iGeom_test.hpp)

endif()

###############################################################################
# iGeom
###############################################################################


if(${BUILD_IGEOM})
  SET(IGEOM_SRC
      iGeom/iGeom.cpp
      iGeom/iGeom.h
      iGeom/iGeomError.cc)

  add_library(iGeom SHARED ${IGEOM_SRC})
  target_link_libraries(iGeom cubit_geom ${CUBITI19})



endif()

###############################################################################
# Build commands
###############################################################################

add_library(svalinn_plugin MODULE ${PLUGIN_SRC} ${DAGMC_SRC} ${MCNP_SRC} ${IGEOMTEST_SRC})
target_link_libraries(svalinn_plugin cubiti cubit_util cubit_geom ${MAKE_WATERTIGHT} ${MOAB_LIBRARIES} ${MCNP2CAD_LIB})

if(${BUILD_IGEOM})
    add_dependencies(svalinn_plugin iGeom)
endif()

if(${BUILD_MCNP_IMPORTER})
  add_dependencies(svalinn_plugin mcnp2cad)
  add_dependencies(mcnp2cad iGeom)
endif()
