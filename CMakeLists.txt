project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 3.1)

# Build options
option(BUILD_DAGMC_EXPORTER "Include DagMC export capability" ON)
option(BUILD_MCNP_IMPORTER  "Include MCNP import capability"  ON)
option(BUILD_IGEOM          "Build iGeom library"             ON)
option(BUILD_IGEOM_TESTS    "Include tests for iGeom"         ON)

# Default to a release build
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Release")
  set(CMAKE_BUILD_TYPE Release)
endif ()

# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Adjust compiler setting for Linux using gcc version 5.0 and higher
if (CMAKE_SYSTEM_NAME MATCHES Linux AND CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
    if (NOT CMAKE_CXX_FLAGS MATCHES _GLIBCXX_USE_CXX11_ABI)
      set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 ${CMAKE_CXX_FLAGS}")
    endif ()
  endif ()
endif ()

# Find Cubit/Trelis
list(APPEND CMAKE_PREFIX_PATH ${CUBIT_ROOT}/bin)
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})
link_directories(${Cubit_DIR})
set(CUBIT_LIBS cubiti19 cubit_geom cubit_util)
message(STATUS "Cubit_DIR: ${Cubit_DIR}")
message(STATUS "CUBIT_INCLUDE_DIRS: ${CUBIT_INCLUDE_DIRS}")

# Setup RPATH
if (NOT CMAKE_INSTALL_RPATH)
  set(CMAKE_INSTALL_RPATH ${Cubit_DIR} ${Cubit_DIR}/plugins/svalinn)
endif ()
message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")

# Add source files and required libraries for the plugin to these variables
# based on which build options are enabled
set(PLUGIN_SRC_FILES   SvalinnPlugin.cpp)
set(PLUGIN_PUB_HEADERS SvalinnPlugin.hpp)
set(PLUGIN_LINK_LIBS ${CUBIT_LIBS})

# DAGMC exporter
if (BUILD_DAGMC_EXPORTER)
  message(STATUS "Building DAGMC exporter")
  add_definitions(-DBUILD_DAGMC_EXPORT=ON)

  # Find DAGMC (DAGMCCMakeConfig.cmake already knows where MOAB is, so we don't
  # need to look for it here separately)
  find_package(DAGMC REQUIRED CONFIG PATHS ${DAGMC_DIR}/lib/cmake)
  link_directories(${DAGMC_LIBRARY_DIRS})
  include_directories(${DAGMC_INCLUDE_DIRS})
  message(STATUS "DAGMC_LIBRARIES: ${DAGMC_LIBRARIES}")
  message(STATUS "DAGMC_LIBRARY_DIRS: ${DAGMC_LIBRARY_DIRS}")
  message(STATUS "DAGMC_INCLUDE_DIRS: ${DAGMC_INCLUDE_DIRS}")

  list(APPEND PLUGIN_SRC_FILES   export_dagmc_cmd/DAGMCExportCommand.cpp)
  list(APPEND PLUGIN_PUB_HEADERS export_dagmc_cmd/DAGMCExportCommand.hpp)
  list(APPEND PLUGIN_LINK_LIBS ${DAGMC_LIBRARIES})
endif ()

if (BUILD_MCNP_IMPORTER OR BUILD_IGEOM_TESTS)
  set(BUILD_IGEOM ON)
endif ()

# iGeom
if (BUILD_IGEOM)
  message(STATUS "Building iGeom")
  set(IGEOM_SRC_FILES iGeom/iGeom.cpp iGeom/iGeomError.cc)
  include_directories(iGeom)
  add_library(iGeom SHARED ${IGEOM_SRC_FILES})
  target_link_libraries(iGeom ${CUBIT_LIBS})
  install(TARGETS iGeom LIBRARY DESTINATION lib)
  set(IGEOM_PUB_HEADERS iGeom/CATag.hpp
                        iGeom/iGeomError.h
                        iGeom/iGeom.h)
  list(APPEND PLUGIN_PUB_HEADERS ${IGEOM_PUB_HEADERS})
  list(APPEND PLUGIN_LINK_LIBS iGeom)
endif ()

# iGeom tests
if (BUILD_IGEOM_TESTS)
  message(STATUS "Building iGeom tests")
  add_definitions(-DBUILD_IGEOM_TESTS=ON)
  list(APPEND PLUGIN_SRC_FILES   iGeom/tests/iGeom_test.cpp)
  list(APPEND PLUGIN_PUB_HEADERS iGeom/tests/iGeom_test.hpp)
endif ()

# MCNP importer
if (BUILD_MCNP_IMPORTER)
  message(STATUS "Building MCNP importer")
  add_definitions(-DBUILD_MCNP_IMPORT=ON)
  set(BUILD_MCNP2CAD_TESTS OFF CACHE BOOL "")
  add_subdirectory(mcnp2cad)
  include_directories(mcnp2cad)
  list(APPEND PLUGIN_SRC_FILES   import_mcnp_cmd/MCNPImp.cpp)
  list(APPEND PLUGIN_PUB_HEADERS import_mcnp_cmd/MCNPImp.hpp)
  list(APPEND PLUGIN_LINK_LIBS mcnp2cad)
endif ()

# Build the plugin
message(STATUS "Building Svalinn plugin")
message(STATUS "PLUGIN_SRC_FILES: ${PLUGIN_SRC_FILES}")
message(STATUS "PLUGIN_PUB_HEADERS: ${PLUGIN_PUB_HEADERS}")
message(STATUS "PLUGIN_LINK_LIBS: ${PLUGIN_LINK_LIBS}")
add_library(svalinn_plugin MODULE ${PLUGIN_SRC_FILES})
set_target_properties(svalinn_plugin PROPERTIES PUBLIC_HEADER
                      "${PLUGIN_PUB_HEADERS}")
target_link_libraries(svalinn_plugin ${PLUGIN_LINK_LIBS})
install(TARGETS svalinn_plugin LIBRARY DESTINATION lib)
