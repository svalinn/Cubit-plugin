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

set(SRC
    iGeom_test.cpp
    iGeom_test.hpp
    MyPlugin.cpp
    MyPlugin.hpp
    DAGMCExportCommand.cpp
    DAGMCExportCommand.hpp
    iGeom_funcs.cpp
    iGeom_funcs.hpp
    mcnp2cad.cpp
    mcnp2cad.hpp
    MCNPImp.cpp
    MCNPImp.hpp
    mcnp2cad/MCNPInput.cpp
    mcnp2cad/MCNPInput.hpp
    mcnp2cad/geometry.cpp
    mcnp2cad/geometry.hpp
    mcnp2cad/volumes.cpp
    mcnp2cad/volumes.hpp
    mcnp2cad/iGeomError.cc)

add_library(dagmc_export_plugin MODULE ${SRC})
target_link_libraries(dagmc_export_plugin cubiti cubit_util cubit_geom ${DAGMC_DIR}/libmakeWatertight.so ${MOAB_LIBRARIES})
