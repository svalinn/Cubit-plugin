project(mhdf)

include_directories( ${HDF5_INCLUDE_DIR} )

set ( mhdf_LIB_SRCS
        include/mhdf.h
        src/adjacency.c
        src/connectivity.c
        src/file-desc.c
        src/file-handle.h  src/file-handle.c
        src/file.c
        src/names-and-paths.h
        src/nodes.c
        src/sets.c
        src/status.h src/status.c
        src/tags.c
        src/util.h src/util.c
)

add_library( mhdf STATIC ${mhdf_LIB_SRCS} )
target_link_libraries( mhdf  ${HDF5_LIBRARIES} )
install(TARGETS mhdf EXPORT MOABTargets
                     RUNTIME DESTINATION bin
                     LIBRARY DESTINATION lib
                     ARCHIVE DESTINATION lib )

add_executable( h5minfo include/mhdf.h example/info.c )
target_link_libraries( h5minfo mhdf )
install(TARGETS h5minfo RUNTIME DESTINATION bin
                        LIBRARY DESTINATION lib
                        ARCHIVE DESTINATION lib )

add_executable( h5mvalidate example/validate.c )
target_link_libraries( h5mvalidate mhdf )
install(TARGETS h5mvalidate RUNTIME DESTINATION bin
                            LIBRARY DESTINATION lib
                            ARCHIVE DESTINATION lib )
