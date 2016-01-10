project(SDK_CommandPluginExample)

cmake_minimum_required(VERSION 2.8)

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

set(SRC
    MyPlugin.cpp
    MyPlugin.hpp
    MyVersionCommand.cpp
    MyVersionCommand.hpp
    MyExportCommand.cpp
    MyExportCommand.hpp)

add_library(my_plugin MODULE ${SRC})
target_link_libraries(my_plugin cubiti cubit_util)
