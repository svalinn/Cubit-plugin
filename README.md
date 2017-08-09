Svalinn plugins and command extensions for Trelis
===================================================

**Beta:** This software is currently under early development.  It has been
demonstrated to work on a wide range of problems, but the build system is not
well developed.

Requirements
============

- Trelis 16.x SDK
- Armadillo
- MOAB (and therefore HDF5 & Lapack/blas)
- DAGMC

Build
======

If you are building the MCNP importer, you should first clone the mcnp2cad library repo into this repo.
```
git clone -b modular --single-branch https://github.com/svalinn/mcnp2cad
```

From here building should work as follows:
```
mkdir bld
cd bld
cmake .. -DCMAKE_PREFIX_PATH=/path/to/Trelis-16.x/bin \
         -DMOAB_DIR=/path/to/MOAB/lib \
         -DDAGMC_DIR=/path/to/DAGMC/
make
```

If not building DAGMC exporter, replace `-DDAGMC_DIR` portion with `-DBUILD_DAGMC_EXPORTER=OFF`.

If not building MCNP importer, replace `-DMCNP2CAD_DIR` portion with `-DBUILD_MCNP_IMPORTER=OFF`.

If testing iGeom functions, add `-DBUILD_IGEOM_TESTS=ON`.

Install
=======

(some of these many need to be completed as `sudo`)
```
PLUGINDIR=/path/to/Trelis-16.x/bin/plugins/svalinn
mkdir $PLUGINDIR
cp libsvalinn_plugin.so $PLUGINDIR
cp libiGeom.so $PLUGINDIR
cp libmcnp2cad.so $PLUGINDIR
cp /path/to/MOAB/lib/libMOAB.so.0 $PLUGINDIR
cp /path/to/DAGMC/lib/libmakeWatertight.so $PLUGINDIR
cp install.sh $PLUGINDIR
cd $PLUGINDIR/../..
bash plugins/svalinn/install.sh
```

You may also need to find and "install" a copy of your HDF5 library in a
fashion similar to the MOAB library above.

# Windows install
Find the Trelis folder, probably "C:\Program Files\Trelis 16.x\"

Copy MOAB.dll to "path\to\Trelis #\bin\" and svalinn_plugin.dll to "path\to\Trelis #\bin\plugins\".

Ensure that you have a copy of hdf5 installed and in your path.  If you do not, download from https://support.hdfgroup.org/HDF5/release/obtain5.html the distribution built for Windows 64-bit using VS 2013. Copy hdf5.dll into "path\to\Trelis #\bin\". 

Distribution
============

The simplest way to make a tarball for distribution is the following command
from the Trelis bin directory on a system with a complete/valid installation:

```
tar czhf ~/tmp/svalinn-plugin.tgz plugins/svalinn
```

This can then be deployed with the following commands from the same directory
on another system:

```
tar xzf ~/Downloads/svalinn-plugin.tgz
bash plugins/svalinn/install.sh
```

Notes & Limitations
====================

This does not currently have a wise/intelligent configuration of the rpath in the plugin and thus requires the correct versions of MOAB and HDF5 to be available in the LD_LIBRARY_PATH.  (This may not be enough??)

