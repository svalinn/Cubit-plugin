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
find_library(CUBITI19 NAMES cubiti19 PATHS ${Cubit_DIR} NO DEFAULT_PATH)
message(STATUS "Cubit_DIR: ${Cubit_DIR}")
message(STATUS "CUBIT_INCLUDE_DIRS: ${CUBIT_INCLUDE_DIRS}")
message(STATUS "CUBITI19: ${CUBITI19}")

# Setup RPATH
set(INSTALL_RPATH_DIRS ${CMAKE_INSTALL_PREFIX}/lib ${Cubit_DIR})
message(STATUS "INSTALL_RPATH_DIRS: ${INSTALL_RPATH_DIRS}")

# Add source files and required libraries for the plugin to these variables based on which build options are enabled
set(PLUGIN_SRC_FILES SvalinnPlugin.cpp)
set(PLUGIN_LINK_LIBS)

# DAGMC exporter
if (BUILD_DAGMC_EXPORTER)
  list(APPEND PLUGIN_SRC_FILES export_dagmc_cmd/DAGMCExportCommand.cpp)
  add_definitions(-DBUILD_DAGMC_EXPORT=ON)

  # Find MOAB
  list(APPEND CMAKE_PREFIX_PATH ${MOAB_ROOT}/lib/cmake/MOAB)
  find_package(MOAB REQUIRED)
  include_directories(${MOAB_INCLUDE_DIRS})
  find_library(MOAB_LIBRARY NAMES MOAB PATHS ${MOAB_DIR}/../.. NO_DEFAULT_PATH)
  list(APPEND PLUGIN_LINK_LIBS ${MOAB_LIBRARY})
  message(STATUS "MOAB_DIR: ${MOAB_DIR}")
  message(STATUS "MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
  message(STATUS "MOAB_LIBRARY: ${MOAB_LIBRARY}")

  # Copy external MOAB library with appropriate symlinks; e.g. libMOAB.so.0 and libMOAB.so.0.0.0
  get_filename_component(MOAB_LIB_DIR ${MOAB_LIBRARY} DIRECTORY)
  get_filename_component(MOAB_LIB_NAME ${MOAB_LIBRARY} NAME)
  file(GLOB MOAB_LIB_INCL_SYMLINKS "${MOAB_LIB_DIR}/${MOAB_LIB_NAME}*")
  install(FILES ${MOAB_LIB_INCL_SYMLINKS} DESTINATION lib)

  # Build make_watertight
  set(MW_SRC_FILES "DAGMC/src/make_watertight/Arc.cpp"
                   "DAGMC/src/make_watertight/check_watertight.cpp"
                   "DAGMC/src/make_watertight/CheckWatertight.cpp"
                   "DAGMC/src/make_watertight/Cleanup.cpp"
                   "DAGMC/src/make_watertight/Gen.cpp"
                   "DAGMC/src/make_watertight/MakeWatertight.cpp"
                   "DAGMC/src/make_watertight/Zip.cpp")
  include_directories(DAGMC/src)
  add_library(makeWatertight SHARED ${MW_SRC_FILES})
  set_target_properties(makeWatertight PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}")
  target_link_libraries(makeWatertight ${MOAB_LIBRARY})
  install(TARGETS makeWatertight LIBRARY DESTINATION lib)
endif ()

if (BUILD_MCNP_IMPORTER OR BUILD_IGEOM_TESTS)
  set(BUILD_IGEOM ON)
endif ()

# iGeom
if (BUILD_IGEOM)
  set(IGEOM_SRC_FILES iGeom/iGeom.cpp iGeom/iGeomError.cc)
  include_directories(iGeom)
  add_library(iGeom SHARED ${IGEOM_SRC_FILES})
  set_target_properties(iGeom PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}")
  target_link_libraries(iGeom cubit_geom ${CUBITI19})
  install(TARGETS iGeom LIBRARY DESTINATION lib)
endif ()

# iGeom tests
if (BUILD_IGEOM_TESTS)
  list(APPEND PLUGIN_SRC_FILES iGeom/tests/iGeom_test.cpp)
  add_definitions(-DBUILD_IGEOM_TESTS=ON)
endif ()

# MCNP importer
if (BUILD_MCNP_IMPORTER)
  add_definitions(-DHAVE_IGEOM_CONE=ON)
  add_definitions(-DBUILD_MCNP_IMPORT=ON)
  file(GLOB MCNP2CAD_SRC_FILES "mcnp2cad/*.cpp")
  include_directories(mcnp2cad)
  add_library(mcnp2cad SHARED ${MCNP2CAD_SRC_FILES})
  set_target_properties(mcnp2cad PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}")
  target_link_libraries(mcnp2cad armadillo iGeom)
  install(TARGETS mcnp2cad LIBRARY DESTINATION lib)
endif ()

# Build the plugin
add_library(svalinn_plugin MODULE ${PLUGIN_SRC_FILES})
set_target_properties(svalinn_plugin PROPERTIES INSTALL_RPATH "${INSTALL_RPATH_DIRS}")
target_link_libraries(svalinn_plugin ${PLUGIN_LINK_LIBS})
install(TARGETS svalinn_plugin LIBRARY DESTINATION lib)
