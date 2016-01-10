project(ComponentExample)

cmake_minimum_required(VERSION 2.8)

# Get the packages needed to use the SDK
find_package(Cubit REQUIRED CONFIG)
include_directories(${CUBIT_INCLUDE_DIRS})

find_package(Qt4 REQUIRED QtCore QtGui)
include(${QT_USE_FILE})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR})

# Include hdf5 files for the .cub file example
include_directories(${HDF5_SOURCE_DIR}/src
                    ${HDF5_SOURCE_DIR}/hl/src
                    ${HDF5_SOURCE_DIR}/hl/c++/src
                    ${HDF5_BINARY_DIR})


# Include our GUI code
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GUI)
add_subdirectory(GUI)

# Include our EventObserver code
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/EventObservers)
add_subdirectory(EventObservers)

# Files that need to run through Qt's MOC
set(MYCOMP_MOCS
    ${GUI_MOC}
    MenuManager.hpp
    ToolbarManager.hpp
    MyWidgetFactory.hpp
    ExportManager.hpp)

qt4_wrap_cpp(moc_src ${MYCOMP_MOCS})
qt4_wrap_ui(ui_src ${GUI_UI})

set(MYCOMP_SRCS
    ${moc_src}
    ${ui_src}
    ${GUI_SRC}
    ${EO_SRC}
    MyComp.cpp
    MyComp.hpp
    MenuManager.cpp
    ToolbarManager.cpp
    PanelManager.cpp
    PanelManager.hpp
    MyWidgetFactory.cpp
    ExportManager.cpp)

# Setup swig. For more information see
set(DOT_I_DEPENDS ${MYCOMP_SRCS})
set(CMAKE_SWIG_FLAGS "")
set(CMAKE_SWIG_OUTDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set_source_files_properties(mycomp.i PROPERTIES CPLUSPLUS ON OBJECT_DEPENDS "${DOT_I_DEPENDS}")
set(SWIG_MODULE_mycomp_EXTRA_DEPS ${DOT_I_DEPENDS})

# Create the python library using swig
swig_add_module(mycomp python mycomp.i ${MYCOMP_SRCS})
swig_link_libraries(mycomp
                    component clarogui navigation pickwidget cubiti clarofw
                    hdf5_hl_cpp hdf5_hl hdf5
                    ${QT_LIBRARIES}
                    ${PYTHON_LIBRARIES})
