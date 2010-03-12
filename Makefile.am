# Subdirectories to build
if HDF5_FILE
  SUBDIRS = mhdf .
else
  SUBDIRS = .
endif

# Things to build
noinst_LTLIBRARIES = libmoabio.la 
libmoabio_la_LIBADD = $(NETCDF_LIBS) $(HDF5_LIBS) 

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I..  -I$(srcdir)/.. -I$(srcdir)/../parallel

if NETCDF_FILE
  MOAB_NETCDF_SRCS = ReadNCDF.cpp  ReadNCDF.hpp \
                     WriteNCDF.cpp WriteNCDF.hpp \
                     WriteSLAC.cpp WriteSLAC.hpp 
else
  MOAB_NETCDF_SRCS =
endif

if HDF5_FILE
  libmoabio_la_LIBADD += mhdf/libmhdf.la
  MOAB_HDF5_SRCS = ReadHDF5.cpp ReadHDF5.hpp WriteHDF5.cpp WriteHDF5.hpp
  AM_CPPFLAGS += -I$(srcdir)/mhdf/include
else
  MOAB_HDF5_SRCS = 
endif
if CCMIO_FILE
  MOAB_CCMIO_SRCS = WriteCCMIO.cpp WriteCCMIO.hpp ReadCCMIO.cpp ReadCCMIO.hpp
  libmoabio_la_LIBADD += $(CCMIO_LIBS)
else
  MOAB_CCMIO_SRCS = 
endif

if HAVE_CGM
  AM_CPPFLAGS += @CGM_CPPFLAGS@ -DCGM @MOAB_CGM_DEFINES@
  MOAB_CGM_SRCS =  ReadCGM.cpp  ReadCGM.hpp
  libmoabio_la_LIBADD += @CGM_LIBS@
else
  MOAB_CGM_SRCS = 
endif

# The list of source files, and any header files that do not need to be installed
libmoabio_la_SOURCES = \
  exodus_order.h \
  patran_order.h \
  ExoIIUtil.cpp \
  ExoIIUtil.hpp \
  FileTokenizer.cpp \
  FileTokenizer.hpp \
  GmshUtil.cpp \
  GmshUtil.hpp \
  ReadABAQUS.cpp \
  ReadABAQUS.hpp \
  ReadGmsh.cpp \
  ReadGmsh.hpp \
  ReadIDEAS.cpp \
  ReadIDEAS.hpp \
  ReadMCNP5.cpp \
  ReadMCNP5.hpp \
  ReadNASTRAN.cpp \
  ReadNASTRAN.hpp \
  ReadSmf.hpp \
  ReadSmf.cpp \
  ReadSms.hpp \
  ReadSms.cpp \
  ReadSTL.cpp \
  ReadSTL.hpp \
  ReadTetGen.cpp \
  ReadTetGen.hpp \
  ReadVtk.cpp \
  ReadVtk.hpp \
  SMF_State.cpp \
  SMF_State.hpp \
  Tqdcfr.cpp \
  Tqdcfr.hpp \
  VtkUtil.cpp \
  VtkUtil.hpp \
  WriteAns.cpp \
  WriteAns.hpp \
  WriteGMV.cpp \
  WriteGMV.hpp \
  WriteGmsh.cpp \
  WriteGmsh.hpp \
  WriteSTL.cpp \
  WriteSTL.hpp \
  WriteSmf.cpp \
  WriteSmf.hpp \
  WriteTemplate.cpp \
  WriteTemplate.hpp \
  WriteVtk.cpp \
  WriteVtk.hpp \
  $(MOAB_NETCDF_SRCS) \
  $(MOAB_HDF5_SRCS) \
  $(MOAB_CCMIO_SRCS) \
  $(MOAB_CGM_SRCS)


libmoabio_la_includedir = $(includedir)
nobase_libmoabio_la_include_HEADERS = moab/ExoIIInterface.hpp

