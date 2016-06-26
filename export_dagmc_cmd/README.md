bld
cd bld
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

