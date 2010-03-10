# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign
LDADD = $(top_builddir)/libMOAB.la

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

# Subdirectories to build
SUBDIRS = 
if HDF5_FILE
  SUBDIRS += mhdf
endif
if PARALLEL
  SUBDIRS += parallel
endif
SUBDIRS += . test tools doc doxygen 

BUILT_SOURCES = MBCN_protos.h

# Things to build
lib_LTLIBRARIES = libMOAB.la 

TESTS = moab_test \
        homxform_test \
        scdseq_test \
        tqdcfr \
        test_adj \
        seq_man_test \
        geom_util_test \
        xform_test \
        obb_test \
	vtk_test \
	adaptive_kd_tree_tests \
	bsp_tree_test \
	file_options_test \
        kd_tree_test \
        var_len_test var_len_test_no_template \
	tag_test \
	mesh_set_test \
	cub_file_test \
	exodus_test \
        gmsh_test \
	nastran_test \
        stl_test \
        mbcn_test \
        bsp_tree_poly_test
#                 merge_test \         # input files no longer exist?
#                 test_tag_server \    # fails

check_PROGRAMS = $(TESTS) \
		 kd_tree_tool \
		 kd_tree_time

#noinst_PROGRAMS = scdseq_timing


# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_builddir) 

# The directory in which to install headers
libMOAB_la_includedir = $(includedir)



# Conditional sources
MOAB_EXTRA_SRCS =
MOAB_EXTRA_HDRS =

if NETCDF_FILE
  MOAB_EXTRA_SRCS +=  ReadNCDF.cpp  ReadNCDF.hpp \
                     WriteNCDF.cpp WriteNCDF.hpp \
                     WriteSLAC.cpp WriteSLAC.hpp 
endif
libMOAB_la_LIBADD = $(NETCDF_LIBS) $(HDF5_LIBS) 
if HDF5_FILE
  libMOAB_la_LIBADD += $(top_builddir)/mhdf/libmhdf.la
  INCLUDES += -I$(srcdir)/mhdf/include
  MOAB_EXTRA_SRCS += ReadHDF5.cpp ReadHDF5.hpp WriteHDF5.cpp WriteHDF5.hpp
endif
if CCMIO_FILE
  MOAB_EXTRA_SRCS += WriteCCMIO.cpp WriteCCMIO.hpp ReadCCMIO.cpp ReadCCMIO.hpp
  libMOAB_la_LIBADD += $(CCMIO_LIBS)
endif

if PARALLEL
  libMOAB_la_LIBADD += $(top_builddir)/parallel/libMOABpar.la
  INCLUDES += -I$(srcdir)/parallel
endif

if HAVE_CGM
  CPPFLAGS += @CGM_CPPFLAGS@ -DCGM @MOAB_CGM_DEFINES@
  MOAB_EXTRA_SRCS +=  ReadCGM.cpp  ReadCGM.hpp
  libMOAB_la_LIBADD += @CGM_LIBS@
endif


# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = moab.make moab.config

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/moab.make
	echo "MOAB_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/moab.make


# Generate a file to be installed in $libdir containing the configuration
# options used for this MOAB build.
# Note: If you modify this, verify that it works correctly for option
#       values containing spaces.
# First line: get unformatted data from config.status
# Second line: extract just the list of options
# Third line: put individual options on separate lines
moab.config: config.status Makefile.am
	./config.status --version | \
	  sed -e 's/.*options "\(.*\)"/\1/p' -e 'd' | \
	  sed -e "s/'[[:space:]]'/'\n'/g" >$@

# The list of source files, and any header files that do not need to be installed
libMOAB_la_SOURCES = \
  AEntityFactory.cpp \
  AEntityFactory.hpp \
  BitTagServer.cpp \
  BitTagServer.hpp \
  DualTool.cpp \
  ElementSequence.hpp \
  EntitySequence.cpp \
  EntitySequence.hpp \
  exodus_order.h \
  ExoIIUtil.cpp \
  ExoIIUtil.hpp \
  FileOptions.cpp \
  FileOptions.hpp \
  FileTokenizer.cpp \
  FileTokenizer.hpp \
  GeomTopoTool.cpp \
  GmshUtil.cpp \
  GmshUtil.hpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  MBAdaptiveKDTree.cpp \
  MBAffineXform.cpp \
  MBAffineXform.hpp \
  MBAlloc.cpp \
  MBAlloc.hpp \
  MBAxisBox.cpp \
  MBAxisBox.hpp \
  MBBSPTree.cpp \
  MBCN.cpp \
  MBCNArrays.hpp \
  MBCartVect.cpp \
  MBMatrix3.cpp \
  MBCore.cpp \
  MBFactory.cpp \
  MBGeomUtil.cpp \
  MBInternals.hpp \
  MBMem.hpp \
  MBMeshSet.cpp \
  MBMeshSet.hpp \
  MBOrientedBox.cpp \
  MBOrientedBoxTreeTool.cpp \
  BSPTreePoly.cpp \
  MBRange.cpp \
  MBRangeSeqIntersectIter.cpp \
  MBRangeSeqIntersectIter.hpp \
  MBReadUtil.cpp \
  MBReadUtil.hpp \
  MBReaderWriterSet.cpp \
  MBSkinner.cpp \
  MBSysUtil.cpp \
  MBSysUtil.hpp \
  MBUtil.cpp \
  MBWriteUtil.cpp \
  MBWriteUtil.hpp \
  MeshSetSequence.cpp \
  MeshSetSequence.hpp \
  MeshTopoUtil.cpp \
  patran_order.h \
  PolyElementSeq.cpp \
  PolyElementSeq.hpp \
  RangeMap.hpp \
  ReadGmsh.cpp \
  ReadGmsh.hpp \
  ReadSms.hpp \
  ReadSms.cpp \
  ReadSTL.cpp \
  ReadSTL.hpp \
  ReadVtk.cpp \
  ReadVtk.hpp \
  ReadSmf.hpp \
  ReadSmf.cpp \
  SMF_State.cpp \
  SMF_State.hpp \
  ReadIDEAS.cpp \
  ReadIDEAS.hpp \
  ReadMCNP5.cpp \
  ReadMCNP5.hpp \
  ReadNASTRAN.cpp \
  ReadNASTRAN.hpp \
  ReadABAQUS.cpp \
  ReadABAQUS.hpp \
  ReadTetGen.cpp \
  ReadTetGen.hpp \
  ScdElementData.cpp \
  ScdElementData.hpp \
  ScdVertexData.cpp \
  ScdVertexData.hpp \
  SequenceData.hpp \
  SequenceData.cpp \
  SequenceManager.cpp \
  SequenceManager.hpp \
  SparseTagSuperCollection.cpp \
  SparseTagSuperCollection.hpp \
  SparseTagCollection.cpp \
  SparseTagCollection.hpp \
  StructuredElementSeq.cpp \
  StructuredElementSeq.hpp \
  SweptElementData.cpp \
  SweptElementData.hpp \
  SweptElementSeq.cpp \
  SweptElementSeq.hpp \
  SweptVertexData.cpp \
  SweptVertexData.hpp \
  TagCompare.hpp \
  TagInfo.cpp \
  TagInfo.hpp \
  TagServer.cpp \
  TagServer.hpp \
  Tqdcfr.cpp \
  Tqdcfr.hpp \
  TypeSequenceManager.cpp \
  TypeSequenceManager.hpp \
  UnstructuredElemSeq.cpp \
  UnstructuredElemSeq.hpp \
  VarLenTag.hpp \
  VertexSequence.hpp \
  VertexSequence.cpp \
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
  WriteTemplate.cpp \
  WriteVtk.cpp \
  WriteVtk.hpp \
  WriteSmf.cpp \
  WriteSmf.hpp \
  $(MOAB_EXTRA_SRCS)

#  MBMem.cpp

# The list of header files which are to be installed
libMOAB_la_include_HEADERS = \
  DualTool.hpp \
  ExoIIInterface.hpp \
  GeomTopoTool.hpp \
  HigherOrderFactory.hpp \
  HomXform.hpp \
  MBAdaptiveKDTree.hpp \
  MBBSPTree.hpp \
  MBCN.hpp \
  MBCN.h \
  MBCN_FCDefs.h \
  MBCN_protos.h \
  MBCartVect.hpp \
  MBCore.hpp \
  MBEntityType.h \
  MBEntityHandle.h \
  MBError.hpp \
  MBForward.hpp \
  MBGeomUtil.hpp \
  MBInterface.hpp \
  MBMatrix3.hpp \
  MBOrientedBox.hpp \
  MBOrientedBoxTreeTool.hpp \
  MBParallelConventions.h \
  BSPTreePoly.hpp \
  MBRange.hpp \
  MBReadUtilIface.hpp \
  MBReaderIface.hpp \
  MBReaderWriterSet.hpp \
  MBSkinner.hpp \
  MBTagConventions.hpp \
  MBTypes.h \
  MBUnknownInterface.hpp \
  MBUtil.hpp \
  MBWriteUtilIface.hpp \
  MBWriterIface.hpp \
  MeshTopoUtil.hpp \
  WriteTemplate.hpp \
  $(MOAB_EXTRA_HDRS)

MBCN_protos.h: MBCN.h tools/iMesh/mkprotos.sh
	$(srcdir)/tools/iMesh/mkprotos.sh MBCN MBCN $< $@ MBCN_FCDefs.h

MBCN_FCDefs.h: iBase_FCDefs.h
	./config.status MBCN_FCDefs.h

CPPFLAGS += -DSRCDIR=$(srcdir)

# Tests and such

moab_test_SOURCES = MBTest.cpp
moab_test_DEPENDENCIES = test/mb_big_test.g test/cell1.gen test/cell2.gen $(top_builddir)/libMOAB.la

homxform_test_SOURCES = HomXform.cpp
homxform_test_CXXFLAGS = -DTEST $(CXXFLAGS)
homxform_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

scdseq_test_SOURCES = scdseq_test.cpp
scdseq_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

#scdseq_timing_SOURCES = scdseq_timing.cpp
#scdseq_timing_DEPENDENCIES = $(top_builddir)/libMOAB.la

tqdcfr_SOURCES = Tqdcfr.cpp
tqdcfr_CXXFLAGS = -DTEST_TQDCFR $(CXXFLAGS)
tqdcfr_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_rms_SOURCES = test_rms.cpp
#test_rms_DEPENDENCIES = $(top_builddir)/libMOAB.la

#merge_test_SOURCES = merge_test.cpp
#merge_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_exo_SOURCES = test_exo.cpp
#test_exo_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_tag_server_SOURCES = TagServer.cpp
#test_tag_server_CXXFLAGS = -DTEST $(CXXFLAGS)
#test_tag_server_DEPENDENCIES = $(top_builddir)/libMOAB.la

test_adj_SOURCES = test_adj.cpp
test_adj_DEPENDENCIES = $(top_builddir)/libMOAB.la

seq_man_test_SOURCES = TestUtil.hpp TestTypeSequenceManager.cpp 

mesh_set_test_SOURCES = Test_MBMeshSet.cpp TestUtil.hpp

geom_util_test_SOURCES = TestUtil.hpp GeomUtilTests.cpp 
geom_util_test_DEPENDENCIES = $(geom_util_test_LDADD)

xform_test_SOURCES = MBAffineXform.cpp
xform_test_CXXFLAGS = -DTEST $(CXXFLAGS)
xform_test_DEPENDENCIES = $(xform_test_LDADD)

obb_test_SOURCES = OBBTest.cpp
obb_test_DEPENDENCIES = $(obb_test_LDADD)

vtk_test_SOURCES = TestUtil.hpp VtkTest.cpp
vtk_test_DEPENDENCIES = $(vtk_test_LDADD)

adaptive_kd_tree_tests_SOURCES = adaptive_kd_tree_tests.cpp
adaptive_kd_tree_tests_DEPENDENCIES = $(adaptive_kd_tree_tests_LDADD)

kd_tree_tool_SOURCES = kd_tree_tool.cpp
kd_tree_tool_DEPENDENCIES = $(kd_tree_tool_LDADD)

kd_tree_time_SOURCES = kd_tree_time.cpp
kd_tree_time_DEPENDENCIES = $(kd_tree_time_LDADD)

kd_tree_test_SOURCES = kd_tree_test.cpp

bsp_tree_test_SOURCES = bsp_tree_test.cpp

file_options_test_SOURCES = TestUtil.hpp FileOptions.cpp 
file_options_test_CPPFLAGS = -DTEST
file_options_test_LDADD = 

var_len_test_SOURCES = TestUtil.hpp VarLenTagTest.cpp
var_len_test_no_template_SOURCES = $(var_len_test_SOURCES)
var_len_test_no_template_CPPFLAGS = -UTEMPLATE_SPECIALIZATION

tag_test_SOURCES = TestUtil.hpp TagTest.cpp
tag_test_DEPENDENCIES = $(tag_test_LDADD)

cub_file_test_SOURCES = TestUtil.hpp cub_file_test.cc
cub_file_test_DEPENDENCIES = $(tag_test_LDADD)

exodus_test_SOURCES = TestUtil.hpp exodus_test.cc
exodus_test_DEPENDENCIES = $(tag_test_LDADD)

gmsh_test_SOURCES = TestUtil.hpp gmsh_test.cc
gmsh_test_DEPENDENCIES = $(tag_test_LDADD)

nastran_test_SOURCES = TestUtil.hpp nastran_test.cc
nastran_test_DEPENDENCIES = $(tag_test_LDADD)

stl_test_SOURCES = TestUtil.hpp stl_test.cc
stl_test_DEPENDENCIES = $(tag_test_LDADD)

mbcn_test_SOURCES = MBCN.hpp MBCN.cpp mbcn_test.cc
mbcn_test_CPPFLAGS = -DSRCDIR=$(srcdir)  # define anything to work around build issue
mbcn_test_LDADD = 

bsp_tree_poly_test_SOURCES = bsp_tree_poly_test.cpp

# Other files to clean up (e.g. output from tests)
DISTCLEANFILES = MBCN_FCDefs.h
MOSTLYCLEANFILES = dumped_acis.sat tree.h5m moab.config

test/mb_big_test.g: $(top_srcdir)/test/mb_big_test.g.gz
	$(ZCAT) $(top_srcdir)/test/mb_big_test.g.gz > test/mb_big_test.g

test/cell1.gen: $(top_srcdir)/test/cell1.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell1.gen.gz > test/cell1.gen

test/cell2.gen: $(top_srcdir)/test/cell2.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell2.gen.gz > test/cell2.gen


# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
