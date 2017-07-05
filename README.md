Svalinn plugins and command extensions for Trelis
===================================================

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
cmake .. -DCMAKE_PREFIX_PATH=/path/to/Trelis-16.x/bin -DCMAKE_INSTALL_BINARY_DIR=bin -DMOAB_DIR=/path/to/MOAB/lib -DDAGMC_DIR=/path/to/DAGMC/ -DMCNP2CAD_DIR=/path/to/libmcnp2cad
make
```

If not building DAGMC exporter, replace `-DDAGMC_DIR` portion with `-DBUILD_DAGMC_EXPORTER=OFF`.

If not building MCNP importer, replace `-DMCNP2CAD_DIR` portion with `-DBUILD_MCNP_IMPORTER=OFF`.

If testing iGeom functions, add `-DBUILD_IGEOM_TESTS=ON`.

Install
=======

(some of these many need to be completed as `sudo`)
```
mkdir /path/to/Trelis-16.x/bin/plugins/svalinn
cp libsvalinn_plugin.so /path/to/Trelis-16.x/bin/plugins/svalinn
cp libiGeom.so /path/to/Trelis-16.x/bin/plugins/svalinn
cp libmcnp2cad.so /path/to/Trelis-16.x/bin/plugins/svalinn
cp /path/to/MOAB/lib/libMOAB.so.0 /path/to/Trelis-16.x/bin/plugins/svalinn
cp /path/to/DAGMC/lib/libmakeWatertight.so /path/to/Trelis-16.x/bin/plugins/svalinn
cd /path/to/Trelis-16.x/bin
ln -s plugins/svalinn/libiGeom.so .
ln -s plugins/svalinn/libmcnp2cad.so .
ln -s plugins/svalinn/libMOAB.so .
ln -s plugins/svalinn/libmakeWatertight.so .
cd plugins
ln -s svalinn/libsvalinn_plugin.so .
```

You may also need to find and "install" a copy of your HDF5 library in a
fashion similar to the MOAB library above.

# Windows install
Find the Trelis folder, probably "C:\Program Files\Trelis 16.x\"

Copy MOAB.dll to "path\to\Trelis #\bin\" and svalinn_plugin.dll to "path\to\Trelis #\bin\plugins\".

Ensure that you have a copy of hdf5 installed and in your path.  If you do not, download from https://support.hdfgroup.org/HDF5/release/obtain5.html the distribution built for Windows 64-bit using VS 2013. Copy hdf5.dll into "path\to\Trelis #\bin\". 


Notes & Limitations
====================

This does not currently have a wise/intelligent configuration of the rpath in the plugin and thus requires the correct versions of MOAB and HDF5 to be available in the LD_LIBRARY_PATH.  (This may not be enough??)

