project(SDK_CommandPluginExample)

# cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIRS})

set(SRC
    MyPlugin.cpp
    MyPlugin.hpp
    DAGMCExportCommand.cpp
    DAGMCExportCommand.hpp)

add_library(dagmc_export_plugin MODULE ${SRC})
target_link_libraries(dagmc_export_plugin cubiti cubit_util ${MOAB_LIBRARIES})
