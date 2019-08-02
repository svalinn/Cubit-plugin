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

Install Trelis
==============

Trelis can be installed by obtaining the Trelis `.deb` package and installing it
with the package manager; i.e.

```
sudo dpkg -i Trelis-16.5.4-Lin64.deb
```

This installs Trelis to `/opt/Trelis-16.5`.

The following commands show how to install the SDK.

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
that the plugin build will take place in the `~/plugin-build` directory, and
they assume that the DAGMC-Trelis repo has been cloned into
`~/plugin-build/DAGMC-Trelis`.

Build MOAB
==========

The following commands show how to build the MOAB dependency.

```
cd ~/plugin-build
mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
ln -sv moab src
cd bld
cmake ../src -DBUILD_SHARED_LIBS=ON \
             -DENABLE_HDF5=ON \
             -DENABLE_BLASLAPACK=ON \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=~/plugin-build/moab \
             -DCMAKE_INSTALL_RPATH=/opt/Trelis-16.5/bin/plugins/svalinn
make -j`grep -c processor /proc/cpuinfo`
make install
```

This reslts in the shared MOAB library being built against the system HDF5
libraries. The MOAB library is now located at
`~/plugin-build/moab/lib/libMOAB.so.5.1.0`.

Build the Plugin
================

Before building the plugin, some external repositories must first be cloned.

```
cd ~/plugin-build/DAGMC-Trelis
git clone https://github.com/svalinn/DAGMC
git clone https://github.com/svalinn/mcnp2cad
```

The following commands show how to build the plugin itself.

```
cd ~/plugin-build
ln -sv DAGMC-Trelis src
mkdir -pv bld
cd bld
cmake ../src -DCUBIT_ROOT=/opt/Trelis-16.5 \
             -DMOAB_ROOT=~/plugin-build/moab \
             -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=~/plugin-build \
             -DCMAKE_INSTALL_RPATH=/opt/Trelis-16.5:/opt/Trelis-16.5/bin/plugins/svalinn
make -j`grep -c processor /proc/cpuinfo`
make install
```

Create the Tarball
==================

The following commands show how to create the tarall for the plugin. These
commands have only been tested on Ubuntu 18.04.

```
cd ~/plugin-build
mkdir -p pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn
cp -pPv ~/plugin-build/lib/* .
cp -pPv /usr/lib/libarmadillo.so.8* .
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
chmod 644 *
cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ..
ln -sv plugins/svalinn/libiGeom.so .
ln -sv plugins/svalinn/libmakeWatertight.so .
ln -sv plugins/svalinn/libmcnp2cad.so .
ln -sv plugins/svalinn/libMOAB.so.5 .
ln -sv plugins/svalinn/libsvalinn_plugin.so .
cd ..
tar --sort=name -czvf svalinn-plugin.tgz bin
mv -v svalinn-plugin.tgz ..
```

The Svalinn pluigin tarball should now be located at
`~/plugin-build/svalinn-plugin.tgz`.

Install the Plugin
==================

To install the plugin, simply run

```
cd /opt/Trelis-16.5
sudo tar -xzvf /path/to/svalinn-plugin.tgz
```
