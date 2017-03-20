project(SDK_CommandPluginExample)

# cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DHAVE_IGEOM_CONE" )

# adjust compiler setting for Linux using gcc version 5.0 and higher
if(CMAKE_SYSTEM_NAME MATCHES Linux AND CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
    if(NOT CMAKE_CXX_FLAGS MATCHES _GLIBCXX_USE_CXX11_ABI)
      set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 ${CMAKE_CXX_FLAGS}")
    endif()
  endif()
endif()

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIRS})

include_directories("${DAGMC_DIR}/../include/")

include_directories("${CMAKE_SOURCE_DIR}/mcnp2cad")

set(PLUGIN_SRC
    MyPlugin.cpp
    MyPlugin.hpp
    DAGMCExportCommand.cpp
    DAGMCExportCommand.hpp
    MCNPImp.cpp
    MCNPImp.hpp
    iGeom_test.cpp
    iGeom_test.hpp)

SET(IGEOM_SRC
    mcnp2cad/iGeom.cpp
    mcnp2cad/iGeom.h
    mcnp2cad/iGeomError.cc)

SET(MCNP2CAD_SRC
    mcnp2cad.cpp
    mcnp2cad.hpp
    mcnp2cad/MCNPInput.cpp
    mcnp2cad/MCNPInput.hpp
    mcnp2cad/geometry.cpp
    mcnp2cad/geometry.hpp
    mcnp2cad/volumes.cpp
    mcnp2cad/volumes.hpp)

add_library(igeom SHARED ${IGEOM_SRC})
target_link_libraries(igeom cubit_geom)

add_library(mcnp2cad SHARED ${MCNP2CAD_SRC})
target_link_libraries(mcnp2cad igeom)

add_library(dagmc_export_plugin MODULE ${PLUGIN_SRC})
target_link_libraries(dagmc_export_plugin cubiti cubit_util cubit_geom ${DAGMC_DIR}/libmakeWatertight.so ${MOAB_LIBRARIES} mcnp2cad)
