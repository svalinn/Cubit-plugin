Svalinn plugins and command extensions for Trelis
=================================================

**Beta:** This software is currently under early development. It has been
demonstrated to work on a wide range of problems, but the build system is not
finalized.

The plugin has been tested and is confirmed to work with Trelis versions 16.5
and 17.1. These build instructions assume that Trelis 17.1 is being used, but
other versions can be used by replacing instances of "17.1" with the version
being used.

Prerequisites
=============

In order to build the plugin, you must have access to Trelis and the Trelis SDK.
Additionally, the following system packages must be present on
your computer:

* EIGEN3
* HDF5

On Ubuntu, these packages can be obtained by running

```
sudo apt install libeigen3-dev libhdf5-dev
```

The following packages are not available from the package manager and must be
built yourself:

* MOAB 5.1.0
* DAGMC

Install Trelis
==============

Trelis can be installed by obtaining the Trelis `.deb` package and installing it
with the package manager; i.e.

```
sudo dpkg -i Trelis-17.1.0-Lin64.deb
```

This installs Trelis to `/opt/Trelis-17.1`.

The Trelis 17 SDK can be installed with these commands:

```
cd /opt
sudo tar -xzvf /path/to/Trelis-SDK-17.1.0-Lin64.tar.gz
```

### Note for Trelis 16

The SDK for Trelis 16 should be unpacked from `/opt/Trelis-16.5` instead of
just `/opt`; i.e.

```
cd /opt/Trelis-16.5
sudo tar -xzvf /path/to/Trelis-SDK-16.5.4-Lin64.tar.gz
```

There is also a bug (or some other unknown issue) in Trelis 16 which requires a
file in the Trelis SDK to be modified. The following commands show how to make
this change. This issue is not present in Trelis 17, so these commands do not
need to be run for Trelis 17.

```
cd /opt/Trelis-16.5/bin
sudo cp -pv CubitExport.cmake CubitExport.cmake.orig
sudo sed -i "s/\"cubit_util\" \"showviz_base\"/\"cubit_util\"/" CubitExport.cmake
```

Notes on Build Instructions
===========================

A non-source directory build is recommended. These build instructions assume
that the plugin build will take place in the `${HOME}/plugin-build` directory,
and they assume that the Trelis-plugin repo has been cloned into
`${HOME}/plugin-build/Trelis-plugin`.

**Before building anything, ensure that the `LD_LIBRARY_PATH` environment
variable is empty**. This variable is not needed to build MOAB, DAGMC, or the
plugin, and if it is not empty it can only cause problems. Ensure that it
remains empty when running Trelis as well.

Build MOAB
==========

MOAB must be built with HDF5 enabled. On Ubuntu 18.04, HDF5 is located in the
`/usr/lib/x86_64-linux-gnu/hdf5/serial` directory, but it may be located
somewhere else on other flavors or versions of Linux. MOAB should be built with
the Eigen matrix algebra library instead of LAPACK. The
`_GLIBCXX_USE_CXX11_ABI=0` flag is required for compatibility with Trelis.

The following commands show how to correctly build the MOAB dependency. If HDF5
is located somewhere other than `/usr/lib/x86_64-linux-gnu/hdf5/serial`, then
replace the directory with the correct one.

```
cd ${HOME}/plugin-build
mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
cd moab
autoreconf -fi
cd ../bld
../moab/configure CXXFLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
                  --enable-shared \
                  --enable-optimize \
                  --disable-debug \
                  --disable-blaslapack \
                  --with-eigen3=/usr/include/eigen3 \
                  --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
                  --prefix=${HOME}/plugin-build/moab
make -j`grep -c processor /proc/cpuinfo`
make install
```

Build DAGMC
===========

The following commands show how to build the DAGMC dependency. The `uwuw` and
`make_watertight` features should be turned on, while other features should be
turned off. The `MOAB_DIR` variable should point to the location of the
previously-built MOAB library. The `_GLIBCXX_USE_CXX11_ABI=0` flag is once again
required.

The following commands show how to correctly build the DAGMC dependency.

```
cd ${HOME}/plugin-build
mkdir -pv DAGMC/bld
cd DAGMC
git clone https://github.com/svalinn/DAGMC -b develop
cd bld
cmake ../DAGMC -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
               -DMOAB_DIR=${HOME}/plugin-build/moab \
               -DBUILD_UWUW=ON \
               -DBUILD_TALLY=OFF \
               -DBUILD_BUILD_OBB=OFF \
               -DBUILD_MAKE_WATERTIGHT=ON \
               -DBUILD_SHARED_LIBS=ON \
               -DBUILD_STATIC_LIBS=OFF \
               -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_INSTALL_PREFIX=${HOME}/plugin-build/DAGMC
make -j`grep -c processor /proc/cpuinfo`
make install
```

Build the Plugin
================

The following commands show how to build the plugin itself. The `CUBIT_ROOT`
variable should point to the location of Trelis. The `DAGMC_DIR` variable should
point to the location of the previously-built DAGMC library.

```
cd ${HOME}/plugin-build
mkdir -pv bld
cd bld
cmake ../Trelis-plugin -DCUBIT_ROOT=/opt/Trelis-17.1 \
                       -DDAGMC_DIR=${HOME}/plugin-build/DAGMC \
                       -DCMAKE_BUILD_TYPE=Release \
                       -DCMAKE_INSTALL_PREFIX=${HOME}/plugin-build
make -j`grep -c processor /proc/cpuinfo`
make install
```

### Submodules

The plugin depends on another external repository called mcnp2cad. mcnp2cad is
available in this repo as a git submodule. It is pulled by default during the
`cmake` configuration step above.

If a custom version of mcnp2cad is needed, this behavior pulling can be disabled
by adding `-DUPDATE_SUBMODULES=OFF` to the `cmake` configuration. mcnp2cad can
then be manually updated with the following commands:

```
cd ${HOME}/plugin-build/Trelis-plugin
git submodule update --init
```

Create the Tarball
==================

The following commands show how to create the tarall for the plugin. Once again,
the location of HDF5 might be different than what is presented here depending on
what flavor or version of Linux is being used.

```
# Set up the directory which will contain the libraries
cd ${HOME}/plugin-build
mkdir -p pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn

# Copy all needed libraries into current directory
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
cp -pPv ${HOME}/plugin-build/moab/lib/libMOAB.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libdagmc.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libmakeWatertight.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libpyne_dagmc.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libuwuw.so* .
cp -pPv ${HOME}/plugin-build/lib/* .
chmod 644 *

# Set the RPATH to be the current directory for the DAGMC libraries
patchelf --set-rpath /opt/Trelis-17.1/bin/plugins/svalinn libMOAB.so
patchelf --set-rpath /opt/Trelis-17.1/bin/plugins/svalinn libdagmc.so
patchelf --set-rpath /opt/Trelis-17.1/bin/plugins/svalinn libmakeWatertight.so
patchelf --set-rpath /opt/Trelis-17.1/bin/plugins/svalinn libpyne_dagmc.so
patchelf --set-rpath /opt/Trelis-17.1/bin/plugins/svalinn libuwuw.so

# Create the Svalinn plugin tarball
cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ../..
tar --sort=name -czvf svalinn-plugin-17.1.tgz bin
mv -v svalinn-plugin-17.1.tgz ..
cd ..
rm -rf pack
```

The Svalinn plugin tarball should now be located at
`${HOME}/plugin-build/svalinn-plugin-17.1.tgz`.

Install the Plugin
==================

To install the plugin, simply run

```
cd /opt/Trelis-17.1
sudo tar -xzvf ${HOME}/plugin-build/svalinn-plugin-17.1.tgz
```

Test the Plugin
===============

Run `trelis`. If the plugin was installed correctly, after the Trelis GUI
finishes loading, the following output should appear in the Trelis command line:

```
Loaded Svalinn plugin.
-- DAGMC export command available.
-- iGeom_test command available.
-- MCNP import command available.
Journaled Command: undo on

Trelis>
```

If this output does not appear, then the plugin was not installed correctly.

To view the available command line options for the MCNP importer, type
`help mcnp` in the Trelis command line. Similarly for the DAGMC exporter, type
`help dagmc` in the Trelis command line.

Some sample files have been included in the `test_plugin` directory of this
repository. Navigate to that directory, then run `trelis`. In the Trelis command
line, type `import mcnp test.i`. This should result in the MCNP geometry and
material mapping contained within the MCNP input file `test.i` being imported
into Trelis.

Next, run `export acis test.sat overwrite attributes_on`. This will save the
geometry in ACIS format to `test.sat`.

Lastly, run
`export dagmc "test.h5m" faceting_tolerance 1e-3 make_watertight verbose`. This
will facet the geometry and save it in a format that can be used by DAGMC.

Test the Plugin (Command Line Mode)
===================================

The plugin can also be run in command line mode, without needing to load the
Trelis GUI. The file `test.jou` in the `test_plugin` directory contains the five
commands mentioned in the previous section of this readme. To execute these
commands in command line mode, run
`trelis -batch -nographics -nojournal test.jou` from the regular command line.

It is through this command line interface that one would replicate the workflows
of years past which involved the now-defunct `mcnp2cad` and `dagmc_preproc`
executables.
