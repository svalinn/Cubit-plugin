project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 2.8)

# Find Cubit/Trelis
list(APPEND CMAKE_PREFIX_PATH ${CUBIT_DIR} ${CUBIT_DIR}/bin)
find_package(Cubit REQUIRED CONFIG)
message(STATUS "Cubit_DIR: ${Cubit_DIR}")
message(STATUS "CUBIT_INCLUDE_DIRS: ${CUBIT_INCLUDE_DIRS}")
include_directories(${CUBIT_INCLUDE_DIRS})
find_library(CUBITI19 NAMES cubiti19 PATHS ${Cubit_DIR} NO DEFAULT_PATH)
message(STATUS "CUBITI19: ${CUBITI19}")

# Options
option(BUILD_DAGMC_EXPORTER  "Include DagMC export capability"  ON)
option(BUILD_MCNP_IMPORTER   "Include MCNP import capability"   ON)
option(BUILD_IGEOM           "Include iGeom"                    ON)
option(BUILD_IGEOM_TESTS     "Include tests for iGeom"          ON)

# C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Top-level plugin files
set(PLUGIN_SRC SvalinnPlugin.cpp)

# DAGMC Exporter
if (BUILD_DAGMC_EXPORTER)
  set(DAGMC_SRC export_dagmc_cmd/DAGMCExportCommand.cpp)
  add_definitions(-DBUILD_DAGMC_EXPORT=ON)

  # Find MOAB
  list(APPEND CMAKE_PREFIX_PATH ${MOAB_DIR} ${MOAB_DIR}/lib {MOAB_DIR}/lib/cmake
                                ${MOAB_DIR}/lib/cmake/MOAB)
  find_package(MOAB REQUIRED)
  message(STATUS "MOAB_DIR: ${MOAB_DIR}")
  message(STATUS "MOAB_INCLUDE_DIRS: ${MOAB_INCLUDE_DIRS}")
  include_directories(${MOAB_INCLUDE_DIRS})
  set(MOAB_LIBRARIES)
  find_library(MOAB_LIBRARIES NAMES MOAB PATHS ${MOAB_DIR}/../.. NO_DEFAULT_PATH)
  message(STATUS "MOAB_LIBRARIES: ${MOAB_LIBRARIES}")

  # Find DAGMC
  set(DAGMC_INCLUDE_DIRS ${DAGMC_DIR}/include)
  message(STATUS "DAGMC_INCLUDE_DIRS: ${DAGMC_INCLUDE_DIRS}")
  include_directories(${DAGMC_INCLUDE_DIRS})
  find_library(MAKE_WATERTIGHT_LIBRARIES NAMES makeWatertight PATHS ${DAGMC_DIR}/lib NO_DEFAULT_PATH)
  message(STATUS "MAKE_WATERTIGHT_LIBRARIES: ${MAKE_WATERTIGHT_LIBRARIES}")
endif ()

if (BUILD_MCNP_IMPORTER OR BUILD_IGEOM_TESTS)
  set(BUILD_IGEOM ON)
endif ()

# iGeom
if (BUILD_IGEOM)
  set(IGEOM_SRC iGeom/iGeom.cpp iGeom/iGeomError.cc)
  include_directories(iGeom)
  add_library(iGeom SHARED ${IGEOM_SRC})
  target_link_libraries(iGeom cubit_geom ${CUBITI19})
  install(TARGETS iGeom LIBRARY DESTINATION lib)
  set(IGEOM_LIB iGeom)
endif ()

# iGeom Tests
if (BUILD_IGEOM_TESTS)
  set(IGEOMTEST_SRC iGeom/tests/iGeom_test.cpp)
  add_definitions(-DBUILD_IGEOM_TESTS=ON)
endif ()

# MCNP Importer
if (BUILD_MCNP_IMPORTER)
  add_subdirectory(mcnp2cad)
  set(MCNP_SRC import_mcnp_cmd/MCNPImp.cpp)
  include_directories(mcnp2cad)
  add_definitions(-DBUILD_MCNP_IMPORT=ON)
endif ()

# Build commands
add_library(svalinn_plugin MODULE ${PLUGIN_SRC} ${DAGMC_SRC} ${MCNP_SRC} ${IGEOMTEST_SRC})
target_link_libraries(svalinn_plugin cubiti cubit_util cubit_geom ${MAKE_WATERTIGHT_LIBRARIES} ${MOAB_LIBRARIES} mcnp2cad)
install(TARGETS svalinn_plugin LIBRARY DESTINATION lib)

# Copy external libraries
get_filename_component(MOABLIB ${MOAB_LIBRARIES} REALPATH)
get_filename_component(MWLIB ${MAKE_WATERTIGHT_LIBRARIES} REALPATH)
install(FILES ${MOABLIB} ${MOAB_LIBRARIES} DESTINATION lib)
install(FILES ${MWLIB} ${MAKE_WATERTIGHT_LIBRARIES} DESTINATION lib)
