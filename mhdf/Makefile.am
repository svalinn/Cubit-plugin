
noinst_LTLIBRARIES = libmhdf.la
libmhdf_la_LIBADD = $(HDF5_LIBS)
bin_PROGRAMS = h5minfo

libmhdf_la_SOURCES = include/mhdf.h \
                     src/adjacency.c \
                     src/connectivity.c \
                     src/file-desc.c \
                     src/file-handle.c \
                     src/file-handle.h \
                     src/file.c \
                     src/names-and-paths.h \
                     src/nodes.c \
                     src/sets.c \
                     src/status.c \
                     src/status.h \
                     src/tags.c \
                     src/util.c \
                     src/util.h

h5minfo_SOURCES = include/mhdf.h example/info.c
h5minfo_LDADD = libmhdf.la $(HDF5_LIBS)

DEFS = -I$(srcdir)/include 

EXTRA_DIST = doc/Doxyfile example/hexes_to_gmsh.c

