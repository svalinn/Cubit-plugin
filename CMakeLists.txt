project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 3.1)

# Build options
option(BUILD_DAGMC_EXPORTER "Include DagMC export capability" ON)
option(BUILD_MCNP_IMPORTER  "Include MCNP import capability"  ON)
option(BUILD_IGEOM          "Build iGeom library"             ON)
option(BUILD_IGEOM_TESTS    "Include tests for iGeom"         ON)

# Use C++11
set(CMAKE_CXX_STANDARD 11)

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
  find_library(MOAB_LIBRARY NAMES MOAB PATHS ${MOAB_DIR}/lib)
  find_library(DAGMC_LIBRARY NAMES dagmc PATHS ${DAGMC_DIR}/lib)
  find_library(PYNE_LIBRARY NAMES pyne_dagmc PATHS ${DAGMC_DIR}/lib)
  find_library(MW_LIBRARY NAMES makeWatertight PATHS ${DAGMC_DIR}/lib)
  include_directories(${MOAB_DIR}/include)
  include_directories(${DAGMC_DIR}/include)
  set(DAGMC_LIBRARIES ${MOAB_LIBRARY} ${DAGMC_LIBRARY} ${PYNE_LIBRARY}
      ${MW_LIBRARY})
  message(STATUS "DAGMC_LIBRARIES: ${DAGMC_LIBRARIES}")
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
