project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 3.1)

# Build options
option(BUILD_DAGMC_EXPORTER "Include DagMC export capability" ON)
option(BUILD_MCNP_IMPORTER  "Include MCNP import capability"  ON)
option(BUILD_IGEOM          "Build iGeom library"             ON)
option(BUILD_IGEOM_TESTS    "Include tests for iGeom"         ON)

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
set(INSTALL_RPATH_DIRS ${Cubit_DIR} ${Cubit_DIR}/plugins/svalinn)
message(STATUS "INSTALL_RPATH_DIRS: ${INSTALL_RPATH_DIRS}")

# Add source files and required libraries for the plugin to these variables
# based on which build options are enabled
set(PLUGIN_SRC_FILES SvalinnPlugin.cpp)
set(PLUGIN_LINK_LIBS ${CUBIT_LIBS})

# DAGMC exporter
if (BUILD_DAGMC_EXPORTER)
  message(STATUS "Building DAGMC exporter")
  add_definitions(-DBUILD_DAGMC_EXPORT=ON)
  list(APPEND PLUGIN_SRC_FILES export_dagmc_cmd/DAGMCExportCommand.cpp)

  # Find DAGMC (DAGMCCMakeConfig.cmake already knows where MOAB is, so we don't
  # need to look for it here separately)
  find_package(DAGMC REQUIRED CONFIG PATHS ${DAGMC_DIR}/lib/cmake)
  link_directories(${DAGMC_LIBRARY_DIRS})
  include_directories(${DAGMC_INCLUDE_DIRS})
  message(STATUS "DAGMC_LIBRARIES: ${DAGMC_LIBRARIES}")
  message(STATUS "DAGMC_LIBRARY_DIRS: ${DAGMC_LIBRARY_DIRS}")
  message(STATUS "DAGMC_INCLUDE_DIRS: ${DAGMC_INCLUDE_DIRS}")

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
  set_target_properties(iGeom PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}")
  target_link_libraries(iGeom ${CUBIT_LIBS})
  install(TARGETS iGeom LIBRARY DESTINATION lib)
  list(APPEND PLUGIN_LINK_LIBS iGeom)
endif ()

# iGeom tests
if (BUILD_IGEOM_TESTS)
  message(STATUS "Building iGeom tests")
  add_definitions(-DBUILD_IGEOM_TESTS=ON)
  list(APPEND PLUGIN_SRC_FILES iGeom/tests/iGeom_test.cpp)
endif ()

# MCNP importer
if (BUILD_MCNP_IMPORTER)
  message(STATUS "Building MCNP importer")
  add_definitions(-DHAVE_IGEOM_CONE=ON)
  add_definitions(-DBUILD_MCNP_IMPORT=ON)
  file(GLOB MCNP2CAD_SRC_FILES "mcnp2cad/*.cpp")
  include_directories(mcnp2cad)
  add_library(mcnp2cad SHARED ${MCNP2CAD_SRC_FILES})
  set_target_properties(mcnp2cad PROPERTIES INSTALL_RPATH
                        "${INSTALL_RPATH_DIRS}")
  target_link_libraries(mcnp2cad ${CUBIT_LIBS} armadillo iGeom)
  install(TARGETS mcnp2cad LIBRARY DESTINATION lib)
  list(APPEND PLUGIN_SRC_FILES import_mcnp_cmd/MCNPImp.cpp)
  list(APPEND PLUGIN_LINK_LIBS mcnp2cad)
endif ()

# Build the plugin
message(STATUS "Building Svalinn plugin")
message(STATUS "PLUGIN_LINK_LIBS: ${PLUGIN_LINK_LIBS}")
add_library(svalinn_plugin MODULE ${PLUGIN_SRC_FILES})
set_target_properties(svalinn_plugin PROPERTIES INSTALL_RPATH
                      "${INSTALL_RPATH_DIRS}")
target_link_libraries(svalinn_plugin ${PLUGIN_LINK_LIBS})
install(TARGETS svalinn_plugin LIBRARY DESTINATION lib)
