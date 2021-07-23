Building the DAGMC Cubit plugin from source
===========================================

This guide is intended for developers and maintainers of the Cubit plugin.
For installation instructions see the [installation guide](README_dev.md).


Prerequisites
=============

In order to build the plugin, you must have access to Cubit and the Cubit SDK.
Additionally, the following system packages must be present on
your computer:

* EIGEN3
* HDF5

On Ubuntu, these packages can be obtained by running

```bash
sudo apt install libeigen3-dev libhdf5-dev
```

The following packages are not available from the package manager and must be
built yourself:

* [MOAB 5.3.0](https://bitbucket.org/fathomteam/moab/src/master/)
* [DAGMC 3.2](https://github.com/svalinn/DAGMC)


Notes on Build Instructions
===========================

A non-source directory build is recommended. These build instructions assume
that the plugin build will take place in the `${HOME}/plugin-build` directory,
and they assume that the Cubit-plugin repo has been cloned into
`${HOME}/plugin-build/Cubit-plugin`.

**Before building anything, ensure that the `LD_LIBRARY_PATH` environment
variable is empty**. Ensure that it remains empty when running Cubit as well.

```bash
unset LD_LIBRARY_PATH
```

Build MOAB
==========

MOAB must be built with HDF5 enabled. On Ubuntu 18.04, HDF5 is located in the
`/usr/lib/x86_64-linux-gnu/hdf5/serial` directory, but it may be located
somewhere else on other flavors or versions of Linux. MOAB should be built with
the Eigen matrix algebra library instead of LAPACK. The
`_GLIBCXX_USE_CXX11_ABI=0` flag is required for compatibility with Cubit.

The following commands show how to correctly build the MOAB dependency. If HDF5
is located somewhere other than `/usr/lib/x86_64-linux-gnu/hdf5/serial`, then
replace the directory with the correct one.

```bash
cd ${HOME}/plugin-build
mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b 5.3.0
cd moab
# patching MOAB CMakeLists.txt to use default find(HDF5)
$SED -i "s/HDF5_MOAB/HDF5/" CMakeLists.txt
cd ..
#end of patch
cd bld
cmake ../moab -DENABLE_HDF5=ON \
        -DHDF5_ROOT=$HDF5_PATH \
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
        -DENABLE_BLASLAPACK=OFF \
        -DENABLE_FORTRAN=OFF \
        $CMAKE_ADDITIONAL_FLAGS \
        -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/moab
make
make install
cd ../..
rm -rf moab/moab moab/bld
```

Build DAGMC
===========

The following commands show how to build the DAGMC dependency. The `uwuw` and
`make_watertight` features should be turned on, while other features should be
turned off. The `MOAB_DIR` variable should point to the location of the
previously-built MOAB library. The `_GLIBCXX_USE_CXX11_ABI=0` flag is once again
required.

The following commands show how to correctly build the DAGMC dependency.

```bash
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
variable should point to the location of Cubit. The `DAGMC_DIR` variable should
point to the location of the previously-built DAGMC library.

```bash
cd ${HOME}/plugin-build
mkdir -pv bld
cd bld
cmake ../Cubit-plugin -DCUBIT_ROOT=PATH_TO_CUBIT \
                       -DDAGMC_DIR=${HOME}/plugin-build/DAGMC \
                       -DCMAKE_BUILD_TYPE=Release \
                       -DCMAKE_INSTALL_PREFIX=${HOME}/plugin-build
make -j`grep -c processor /proc/cpuinfo`
make install
```

Submodules
==========

The plugin depends on another external repository called mcnp2cad. mcnp2cad is
available in this repo as a git submodule. It is pulled by default during the
`cmake` configuration step above.

If a custom version of mcnp2cad is needed, this behavior pulling can be disabled
by adding `-DUPDATE_SUBMODULES=OFF` to the `cmake` configuration. mcnp2cad can
then be manually updated with the following commands:

```bash
cd ${HOME}/plugin-build/Cubit-plugin
git submodule update --init
```

Create the Tarball
==================

The following commands show how to create the tarall for the plugin. Once again,
the location of HDF5 might be different than what is presented here depending on
what flavor or version of Linux is being used.


Set up the directory which will contain the libraries
```bash
cd ${HOME}/plugin-build
mkdir -p pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn
```

Copy all needed libraries into current directory
```bash
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
cp -pPv ${HOME}/plugin-build/moab/lib/libMOAB.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libdagmc.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libmakeWatertight.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libpyne_dagmc.so* .
cp -pPv ${HOME}/plugin-build/DAGMC/lib/libuwuw.so* .
cp -pPv ${HOME}/plugin-build/lib/* .
chmod 644 *
```

The resulting shared library objects require an update to their RPATH (runtime path) 
attribute to ensure that the correct set of libraries is discovered when starting 
Cubit and loading the plugin library.
```bash
patchelf --set-rpath PATH_TO_CUBIT/bin/plugins/svalinn libMOAB.so
patchelf --set-rpath PATH_TO_CUBIT/bin/plugins/svalinn libdagmc.so
patchelf --set-rpath PATH_TO_CUBIT/bin/plugins/svalinn libmakeWatertight.so
patchelf --set-rpath PATH_TO_CUBIT/bin/plugins/svalinn libpyne_dagmc.so
patchelf --set-rpath PATH_TO_CUBIT/bin/plugins/svalinn libuwuw.so
```

Create the Svalinn plugin tarball
```bash
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
```bash
cd PATH_TO_CUBIT
sudo tar -xzvf ${HOME}/plugin-build/svalinn-plugin.tgz
```
