Svalinn plugins and command extensions for Trelis
=================================================

**Beta:** This software is currently under early development.  It has been
demonstrated to work on a wide range of problems, but the build system is not
well developed.

Prerequisites
=============

In order to build the plugin, you must have access to Trelis-16.5 and the
Trelis-16.5 SDK. Additionally, the following system packages must be present on
your computer:

* BLAS/LAPACK
* Armadillo
* HDF5

On Ubuntu, these packages can be obtained by running

```
sudo apt install libblas-dev liblapack-dev libarmadillo-dev libhdf5-dev
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
sudo dpkg -i Trelis-16.5.3-Lin64.deb
```

This installs Trelis to `/opt/Trelis-16.5`.

The Trelis SDK can be installed with these commands:

```
cd /opt/Trelis-16.5
sudo tar -xzvf /path/to/Trelis-SDK-16.5.3-Lin64.tar.gz
```

There is currently a bug (or some other unknown issue) which requires a file in
the Trelis SDK to be modified. The following commands show how to make this
change.

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

The following commands show how to build the MOAB dependency.

```
cd ${HOME}/plugin-build
mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
cd moab
autoreconf -fi
cd ..
ln -sv moab src
cd bld
../src/configure --disable-blaslapack \
                 --enable-shared \
                 --enable-optimize \
                 --disable-debug \
                 --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
                 --prefix=${HOME}/plugin-build/moab
make -j`grep -c processor /proc/cpuinfo`
make install
```

This reslts in the MOAB library being built against the system HDF5 libraries.

Build DAGMC
===========

The following commands show how to build the DAGMC dependency.

```
cd ${HOME}/plugin-build
mkdir -pv DAGMC/bld
cd DAGMC
git clone https://github.com/svalinn/DAGMC -b develop
ln -sv DAGMC src
cd bld
cmake ../src -DMOAB_DIR=${HOME}/plugin-build/moab \
             -DBUILD_UWUW=OFF \
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

This reslts in the DAGMC library being built against the previously-built MOAB
library.

Build the Plugin
================

Before building the plugin, the external mcnp2cad repository must first be
cloned.

```
cd ${HOME}/plugin-build/Trelis-plugin
git clone https://github.com/svalinn/mcnp2cad -b master
```

The following commands show how to build the plugin itself.

```
cd ${HOME}/plugin-build
ln -sv Trelis-plugin src
mkdir -pv bld
cd bld
cmake ../src -DCUBIT_ROOT=/opt/Trelis-16.5 \
             -DMOAB_DIR=${HOME}/plugin-build/moab \
             -DDAGMC_DIR=${HOME}/plugin-build/DAGMC \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=${HOME}/plugin-build \
             -DCMAKE_INSTALL_RPATH=/opt/Trelis-16.5:/opt/Trelis-16.5/bin/plugins/svalinn
make -j`grep -c processor /proc/cpuinfo`
make install
```

Create the Tarball
==================

The following commands show how to create the tarall for the plugin. These
commands have only been tested on Ubuntu 18.04.

```
# Set up the directory which will contain the libraries
cd ${HOME}/plugin-build
mkdir -p pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn

# Copy all needed libraries into current directory
cp -pPv ${HOME}/plugin-build/lib/* .
cp -pPv ${HOME}/plugin-build/moab/lib/libMOAB.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libdagmc.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libmakeWatertight.so* .
cp -pPv /usr/lib/libarmadillo.so.8* .
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
chmod 644 *

# Set the RPATH to be the current directory for the DAGMC libraries
patchelf --set-rpath /opt/Trelis-16.5/bin/plugins/svalinn libdagmc.so
patchelf --set-rpath /opt/Trelis-16.5/bin/plugins/svalinn libmakeWatertight.so

# Create the Svalinn plugin tarball
cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ../..
tar --sort=name -czvf svalinn-plugin.tgz bin
mv -v svalinn-plugin.tgz ..
cd ..
rm -rf pack
```

The Svalinn plugin tarball should now be located at
`${HOME}/plugin-build/svalinn-plugin.tgz`.

Install the Plugin
==================

To install the plugin, simply run

```
cd /opt/Trelis-16.5
sudo tar -xzvf ${HOME}/plugin-build/svalinn-plugin.tgz
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
