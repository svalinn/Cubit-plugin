
  # MOAB Library
  set ( MOABIO_LIB_SRCS
    ExoIIUtil.cpp
    FileTokenizer.cpp
    GmshUtil.cpp
    ReadABAQUS.cpp
    ReadGmsh.cpp
    ReadIDEAS.cpp
    ReadMCNP5.cpp
    ReadNASTRAN.cpp
    ReadSmf.cpp
    ReadSms.cpp
    ReadSTL.cpp
    ReadTemplate.cpp
    ReadTetGen.cpp
    ReadTxt.cpp
    ReadVtk.cpp
    SMF_State.cpp
    Tqdcfr.cpp
    VtkUtil.cpp
    WriteAns.cpp
    WriteGMV.cpp
    WriteGmsh.cpp
    WriteSTL.cpp
    WriteSmf.cpp 
    WriteTemplate.cpp
    WriteVtk.cpp
  )
  include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_BINARY_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${MOAB_BINARY_DIR}/src/parallel
  )

  if ( MOAB_USE_NETCDF )
    find_package( NetCDF )
    if ( NetCDF_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DNETCDF_FILE" )
      set ( MOABIO_LIB_SRCS
        ${MOABIO_LIB_SRCS}
        ReadNCDF.cpp
        WriteNCDF.cpp
        WriteSLAC.cpp
      )
      include_directories(
        ${NetCDF_INCLUDE_DIRECTORIES}
      )
    endif ( NetCDF_FOUND )
  endif ( MOAB_USE_NETCDF )

  if ( MOAB_USE_HDF )
    # HDF5
    find_package( HDF5 )
    if ( HDF5_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DHDF5_FILE" )
      check_function_exists( H5Pset_fapl_mpio MOAB_HDF_HAVE_PARALLEL )
      set ( MOABIO_LIB_SRCS
        ${MOABIO_LIB_SRCS}
        ReadHDF5.cpp
        WriteHDF5.cpp
      )
      include_directories(
        ${HDF5_INCLUDE_DIRECTORIES}
        ${MOAB_SOURCE_DIR}/mhdf/include
      )
      add_subdirectory( mhdf )
    endif ( HDF5_FOUND )
  endif ( MOAB_USE_HDF )

  set_source_files_properties( ${MOABIO_LIB_SRCS}
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
  )
  add_library( moabio
    ${MOABIO_LIB_SRCS}
  )

  if ( MOAB_USE_HDF AND HDF5_FOUND )
    target_link_libraries( moabio
      mhdf
      ${HDF5_LIBRARIES}
    )
  endif ( MOAB_USE_HDF AND HDF5_FOUND )

  if ( MOAB_USE_NETCDF AND NetCDF_FOUND )
    target_link_libraries( moabio
      ${NetCDF_LIBRARIES}
    )
  endif ( MOAB_USE_NETCDF AND NetCDF_FOUND )
