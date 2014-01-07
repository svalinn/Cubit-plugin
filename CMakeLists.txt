include_directories( ${MOAB_SOURCE_DIR}/src/io
                     ${MOAB_SOURCE_DIR}/src/mhdf )
set( MOABIO_LIBS )
set( MOABIO_DEFINES "-DIS_BUILDING_MB" )
set( MOABIO_INSTALL_HEADERS moab/ExoIIInterface.hpp )
set( MOABIO_LIB_SRCS
        exodus_order.h
        patran_order.h
        IODebugTrack.hpp  IODebugTrack.cpp
        ExoIIUtil.hpp     ExoIIUtil.cpp
        FileTokenizer.hpp FileTokenizer.cpp
        GmshUtil.hpp      GmshUtil.cpp
        ReadABAQUS.hpp    ReadABAQUS.cpp
        ReadGmsh.hpp      ReadGmsh.cpp
        ReadIDEAS.hpp     ReadIDEAS.cpp
        ReadMCNP5.hpp     ReadMCNP5.cpp
        ReadNASTRAN.hpp   ReadNASTRAN.cpp
        ReadSmf.hpp       ReadSmf.cpp
        ReadSms.hpp       ReadSms.cpp
        ReadSTL.hpp       ReadSTL.cpp
        ReadTemplate.hpp  ReadTemplate.cpp
        ReadTetGen.hpp    ReadTetGen.cpp
        ReadTxt.hpp       ReadTxt.cpp
        ReadVtk.hpp       ReadVtk.cpp
        SMF_State.hpp     SMF_State.cpp
        Tqdcfr.hpp        Tqdcfr.cpp
        VtkUtil.hpp       VtkUtil.cpp
        WriteAns.hpp      WriteAns.cpp
        WriteGMV.hpp      WriteGMV.cpp
        WriteGmsh.hpp     WriteGmsh.cpp
        WriteSTL.hpp      WriteSTL.cpp
        WriteSmf.hpp      WriteSmf.cpp
        WriteTemplate.hpp WriteTemplate.cpp
        WriteVtk.hpp      WriteVtk.cpp )

if( MOAB_USE_NETCDF )
  set( MOAB_NETCDF_SRCS
          ReadNCDF.hpp      ReadNCDF.cpp
          WriteNCDF.hpp     WriteNCDF.cpp
          WriteSLAC.hpp     WriteSLAC.cpp
          ReadNC.hpp        ReadNC.cpp
          NCHelper.hpp      NCHelper.cpp
          NCHelperEuler.hpp NCHelperEuler.cpp
          NCHelperFV.hpp    NCHelperFV.cpp
          NCHelperHOMME.hpp NCHelperHOMME.cpp
          NCHelperMPAS.hpp  NCHelperMPAS.cpp
          ReadGCRM.hpp      ReadGCRM.cpp )
  set( MOABIO_LIBS ${MOABIO_LIBS} ${NetCDF_LIBRARIES} )
elseif( MOAB_USE_PNETCDF )
  set( MOAB_NETCDF_SRCS
          ReadNC.hpp        ReadNC.cpp
          NCHelper.hpp      NCHelper.cpp
          NCHelperEuler.hpp NCHelperEuler.cpp
          NCHelperFV.hpp    NCHelperFV.cpp
          NCHelperHOMME.hpp NCHelperHOMME.cpp
          NCHelperMPAS.hpp  NCHelperMPAS.cpp
          ReadGCRM.hpp      ReadGCRM.cpp )
endif()

if ( MOAB_USE_CGNS )
  set( MOAB_CGNS_SRCS ReadCGNS.hpp ReadCGNS.cpp
                      WriteCGNS.hpp WriteCGNS.cpp)
endif()

if( MOAB_USE_HDF )
  add_subdirectory( mhdf )
  set( MOAB_HDF5_SRCS
          HDF5Common.hpp HDF5Common.cpp
          ReadHDF5.hpp ReadHDF5.cpp
          ReadHDF5Dataset.hpp ReadHDF5Dataset.cpp
          ReadHDF5VarLen.hpp ReadHDF5VarLen.cpp
          WriteHDF5.cpp WriteHDF5.hpp )
  set( MDHF_OBJ $<TARGET_OBJECTS:mhdf> PARENT_SCOPE )
  set( MOABIO_LIBS ${MOABIO_LIBS} ${HDF5_LIBRARIES} )
  check_function_exists( H5Pset_fapl_mpio MOAB_HDF_HAVE_PARALLEL )
endif()

if( MOAB_USE_CCMIO )
  set( MOAB_CCMIO_SRCS WriteCCMIO.hpp WriteCCMIO.cpp 
                       ReadCCMIO.hpp ReadCCMIO.cpp )
endif()

if( MOAD_USE_DAMSEL )
  set( MOAB_DAMSEL_SRCS
          WriteDamsel.hpp WriteDamsel.cpp 
          ReadDamsel.hpp  ReadDamsel.cpp
          DamselUtil.hpp  DamselUtil.cpp )
endif()

if( CGM_FOUND )
  set( MOAB_CGM_SRCS ReadCGM.hpp ReadCGM.cpp )
endif()

add_library( moabio OBJECT ${MOABIO_LIB_SRCS} ${MOAB_NETCDF_SRCS}
                           ${MOAB_CGNS_SRCS}  ${MOAB_HDF5_SRCS}
                           ${MOAB_CCMIO_SRCS} ${MOAB_DAMSEL_SRCS}
                           ${MOAB_CGM_SRCS}   ${MOABIO_INSTALL_HEADERS} )
set_target_properties( moabio PROPERTIES COMPILE_FLAGS "${MOABIO_DEFINES} ${MOAB_DEFINES}")
include(AutoconfHeader)
moab_install_headers(${MOABIO_INSTALL_HEADERS})
