project(mhdf)

include_directories(
  ${HDF5_INCLUDE_DIRECTORIES}
)

set ( mhdf_LIB_SRCS
  src/adjacency.c
  src/connectivity.c
  src/file-handle.c
  src/file-desc.c
  src/file.c
  src/nodes.c
  src/sets.c
  src/status.c
  src/status.h
  src/tags.c
  src/util.c
)

add_library( mhdf
  ${mhdf_LIB_SRCS}
)

target_link_libraries(
  mhdf
  ${HDF5_LIBRARIES}
)

#include/mhdf.h \
#src/file-handle.h \
#src/names-and-paths.h \
#src/util.h
