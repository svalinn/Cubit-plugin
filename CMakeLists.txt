project(SvalinnCommandPlugin)

cmake_minimum_required(VERSION 3.0)

# Build options
option(UPDATE_SUBMODULES    "Automatically update git submodules" ON)
option(BUILD_DAGMC_EXPORTER "Include DagMC export capability"     ON)
option(BUILD_MCNP_IMPORTER  "Include MCNP import capability"      ON)
option(BUILD_IGEOM          "Build iGeom library"                 ON)
option(BUILD_IGEOM_TESTS    "Include tests for iGeom"             ON)

# MCNP importer and iGeom tests require iGeom
if (BUILD_MCNP_IMPORTER OR BUILD_IGEOM_TESTS)
  set(BUILD_IGEOM ON)
endif ()

# Default to a release build
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE not specified, defaulting to Release")
  set(CMAKE_BUILD_TYPE Release)
endif ()

#enable windows export flags
if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set(BUILD_SHARED_LIBS TRUE)
endif()

# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Find Cubit
if (NOT CUBIT_ROOT)
  message(FATAL_ERROR "CUBIT_ROOT not specified")
endif ()
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
set(PLUGIN_SRC_FILES SvalinnPlugin.cpp)
set(PLUGIN_LINK_LIBS ${CUBIT_LIBS})

# HDF5 for mac
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(HDF5_USE_STATIC_LIBRARIES ON)
  find_package(HDF5 REQUIRED)
  include_directories(${HDF5_INCLUDE_DIRS})
  list(APPEND PLUGIN_LINK_LIBS ${HDF5_LIBRARIES})
endif()

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

  list(APPEND PLUGIN_SRC_FILES export_dagmc_cmd/DAGMCExportCommand.cpp)
  list(APPEND PLUGIN_LINK_LIBS ${DAGMC_LIBRARIES})
endif ()

# iGeom
if (BUILD_IGEOM)
  message(STATUS "Building iGeom")
  set(IGEOM_SRC_FILES iGeom/iGeom.cpp iGeom/iGeomError.cc)
  include_directories(iGeom)
  add_library(iGeom SHARED ${IGEOM_SRC_FILES})
  target_link_libraries(iGeom ${CUBIT_LIBS})
   target_include_directories(iGeom PUBLIC
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
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
  # Update mcnp2cad submodule as neeeded
  if (UPDATE_SUBMODULES)
    find_package(Git)
    if (GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
      message(STATUS "Updating mcnp2cad submodule")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if (NOT GIT_SUBMOD_RESULT EQUAL 0)
        message(FATAL_ERROR "git submodule update --init failed with \
            ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif ()
    endif ()
  endif ()

  message(STATUS "Building MCNP importer")
  add_definitions(-DBUILD_MCNP_IMPORT=ON)
  set(STANDALONE_MCNP2CAD  OFF CACHE BOOL "")
  set(BUILD_MCNP2CAD_TESTS OFF CACHE BOOL "")
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mcnp2cad/CMakeLists.txt)
    add_subdirectory(mcnp2cad)
  else ()
    message(FATAL_ERROR "mcnp2cad submodule not found.\n"
                        "Run \"git submodule update --init --recursive\" "
                        "to pull the mcnp2cad submodule.")
  endif ()
  include_directories(mcnp2cad)
  list(APPEND PLUGIN_SRC_FILES import_mcnp_cmd/MCNPImp.cpp)
  list(APPEND PLUGIN_LINK_LIBS mcnp2cad)
endif ()

# Build the plugin
message(STATUS "Building Svalinn plugin")
message(STATUS "PLUGIN_SRC_FILES: ${PLUGIN_SRC_FILES}")
message(STATUS "PLUGIN_LINK_LIBS: ${PLUGIN_LINK_LIBS}")
add_library(svalinn_plugin MODULE ${PLUGIN_SRC_FILES})
target_include_directories(svalinn_plugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(svalinn_plugin ${PLUGIN_LINK_LIBS})
install(TARGETS svalinn_plugin LIBRARY DESTINATION lib)
