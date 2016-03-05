project(SDK_CommandPluginExample)

cmake_minimum_required(VERSION 3.1)

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

set(SRC
    MyPlugin.cpp
    MyPlugin.hpp
    DAGMCExportCommand.cpp
    DAGMCExportCommand.hpp)

add_library(my_plugin MODULE ${SRC})
target_link_libraries(my_plugin cubiti cubit_util)
