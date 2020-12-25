#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "rlm_activate" for configuration "Release"
set_property(TARGET rlm_activate APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(rlm_activate PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/MacOS/rlm_activate"
  )

list(APPEND _IMPORT_CHECK_TARGETS rlm_activate )
list(APPEND _IMPORT_CHECK_FILES_FOR_rlm_activate "${_IMPORT_PREFIX}/MacOS/rlm_activate" )

# Import target "trelis_lm" for configuration "Release"
set_property(TARGET trelis_lm APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(trelis_lm PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "Qt5::Core;Qt5::Network;Qt5::Widgets"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libtrelis_lm.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libtrelis_lm.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS trelis_lm )
list(APPEND _IMPORT_CHECK_FILES_FOR_trelis_lm "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libtrelis_lm.dylib" )

# Import target "clarofw" for configuration "Release"
set_property(TARGET clarofw APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(clarofw PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "trelis_lm;Qt5::Gui;Qt5::Core;cubit_python"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libclarofw.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libclarofw.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS clarofw )
list(APPEND _IMPORT_CHECK_FILES_FOR_clarofw "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libclarofw.dylib" )

# Import target "navigation" for configuration "Release"
set_property(TARGET navigation APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(navigation PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "Qt5::Widgets"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libnavigation.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libnavigation.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS navigation )
list(APPEND _IMPORT_CHECK_FILES_FOR_navigation "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libnavigation.dylib" )

# Import target "clarogui" for configuration "Release"
set_property(TARGET clarogui APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(clarogui PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "clarofw;navigation;Qt5::Core;Qt5::Widgets;cubit_python"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libclarogui.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libclarogui.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS clarogui )
list(APPEND _IMPORT_CHECK_FILES_FOR_clarogui "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libclarogui.dylib" )

# Import target "gtcAttrib" for configuration "Release"
set_property(TARGET gtcAttrib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gtcAttrib PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "SpaACIS"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libgtcAttrib.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libgtcAttrib.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS gtcAttrib )
list(APPEND _IMPORT_CHECK_FILES_FOR_gtcAttrib "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libgtcAttrib.dylib" )

# Import target "cubit_sizing_source" for configuration "Release"
set_property(TARGET cubit_sizing_source APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_sizing_source PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "cubit_util"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_sizing_source.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libcubit_sizing_source.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_sizing_source )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_sizing_source "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_sizing_source.dylib" )

# Import target "cubit_mdb" for configuration "Release"
set_property(TARGET cubit_mdb APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_mdb PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_mdb.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_mdb )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_mdb "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_mdb.a" )

# Import target "cubit_sim" for configuration "Release"
set_property(TARGET cubit_sim APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_sim PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_sim.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_sim )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_sim "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_sim.a" )

# Import target "material_commands" for configuration "Release"
set_property(TARGET material_commands APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(material_commands PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/plugins/libmaterial_commands.so"
  IMPORTED_NO_SONAME_RELEASE "TRUE"
  )

list(APPEND _IMPORT_CHECK_TARGETS material_commands )
list(APPEND _IMPORT_CHECK_FILES_FOR_material_commands "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/plugins/libmaterial_commands.so" )

# Import target "cubit_parsing_core" for configuration "Release"
set_property(TARGET cubit_parsing_core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_parsing_core PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_parsing_core.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_parsing_core )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_parsing_core "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_parsing_core.a" )

# Import target "cubiti" for configuration "Release"
set_property(TARGET cubiti APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubiti PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "trelis_lm;cubit_sizing_source;showviz_cubit;showviz_sim;showviz_mesh;showviz_base;cubit_smd;gtcAttrib;SpaACIS;cubit_geom;mg_tetra_hpc;mg_adapt;mg_tetra;mg_surfopt;mg_cadsurf;meshgems;meshgems_stubs;Qt5::Core"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubiti19.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libcubiti19.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubiti )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubiti "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubiti19.dylib" )

# Import target "cubit_python" for configuration "Release"
set_property(TARGET cubit_python APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_python PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_python.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libcubit_python.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_python )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_python "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libcubit_python.dylib" )

# Import target "cubit_python3" for configuration "Release"
set_property(TARGET cubit_python3 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_python3 PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/python3/libcubit_python3.so"
  IMPORTED_NO_SONAME_RELEASE "TRUE"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_python3 )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_python3 "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/python3/libcubit_python3.so" )

# Import target "cubit_python2" for configuration "Release"
set_property(TARGET cubit_python2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cubit_python2 PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/python2/libcubit_python2.so"
  IMPORTED_NO_SONAME_RELEASE "TRUE"
  )

list(APPEND _IMPORT_CHECK_TARGETS cubit_python2 )
list(APPEND _IMPORT_CHECK_FILES_FOR_cubit_python2 "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/python2/libcubit_python2.so" )

# Import target "treemodel" for configuration "Release"
set_property(TARGET treemodel APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(treemodel PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "cubiti;Qt5::Core"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libtreemodel.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libtreemodel.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS treemodel )
list(APPEND _IMPORT_CHECK_FILES_FOR_treemodel "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libtreemodel.dylib" )

# Import target "pickwidget" for configuration "Release"
set_property(TARGET pickwidget APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(pickwidget PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "Qt5::Widgets"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libpickwidget.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libpickwidget.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS pickwidget )
list(APPEND _IMPORT_CHECK_FILES_FOR_pickwidget "${_IMPORT_PREFIX}/Trelis-17.1.app/Contents/MacOS/libpickwidget.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
