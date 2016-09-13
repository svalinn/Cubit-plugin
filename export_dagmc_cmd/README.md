```
mkdir bld
cd bld
cmake .. -DCMAKE_PREFIX_PATH=/path/to/Trelis-16.0/bin -DCMAKE_INSTALL_BINARY_DIR=bin -DMOAB_DIR=/path/to/MOAB/lib
make
```

Install
=======

Installation of the additional command is as simple as copying the generated
shared object library (`my_plugin.so`) to the Trelis plugin driectory:
`/path/to/Trelis-16.0/bin/plugins`.

Installation of a component requires that it be manually loaded in the Components dialog accessible from the Trelis Tools menu.  When loading the component, you will asked for the case-sensitive C++ class name `MyComp`.

It doesn't matter what directory the component is in, but proper installation is probably best in a subirectory of `/path/to/Trelis16.0/bin/plugins`.

Notes & Limitations
====================

This does not currently have a wise/intelligent configuration of the rpath in the plugin and thus requires the correct versions of MOAB and HDF5 to be available in the LD_LIBRARY_PATH.  (This may not be enough??)

To date, PPHW has successfully built this by building & installing a special instance of MOAB into the plugins directory: .../plugins/dagmc/lib, .../plugins/dagmc/include, etc.  It should be possible to only place libMOAB.so and libhdf5.so in that location and have the plugin work.
