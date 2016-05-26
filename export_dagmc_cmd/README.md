
THIS PLUGIN IS BROKEN ON AN HDF5 CONFLICT.

Purpose
==========

This adds a command to Trelis to export a geometry in the MOAB-based facet format necessary for DAGMC.

Dependencies
=============

This plugin/extension requires a licensed installation of Trelis 16.x and
the Trelis 16.x SDK.  The SDK should be untarred in the Trelis install
directory.

This plugin also depends on:
* swig
* python2.7-dev
* MOAB (and thus HDF5)

Build
======

```
mkdir bld
cd blb
cmake .. -DCMAKE_PREFIX_PATH=/path/to/Trelis-16.0/bin -DCMAKE_INSTALL_BINARY_DIR=bin -DMOAB_DIR=/path/to/MOAB/lib
make
```

Install
=======

Installation of the additional command is as simple as copying the generated
shared object library (`my_plugin.so`) to the Trelis plugin driectory:
`/path/to/Trelis-15.1/bin/plugins`.

Installation of a component requires that it be manually loaded in the Components dialog accessible from the Trelis Tools menu.  When loading the component, you will asked for the case-sensitive C++ class name `MyComp`.

It doesn't matter what directory the component is in, but proper installation is probably best in a subirectory of `/path/to/Trelis15.1/bin/plugins`.

