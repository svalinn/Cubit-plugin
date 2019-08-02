#!/bin/bash

set -e

cd ..

jobs=`grep -c processor /proc/cpuinfo`
base_dir=${PWD}
trelis_dir=/opt/Trelis-16.5
plugin_dir=${trelis_dir}/bin/plugins/svalinn
tarball=svalinn-plugin.tgz

# ------------------------------------------------------------------------------
# Build MOAB
# ------------------------------------------------------------------------------

cd ${base_dir}
rm -rf moab
mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
ln -sv moab src
cd bld

cmake_string=""
cmake_string+=" -DBUILD_SHARED_LIBS=ON"
cmake_string+=" -DENABLE_HDF5=ON"
cmake_string+=" -DENABLE_BLASLAPACK=ON"
cmake_string+=" -DCMAKE_BUILD_TYPE=Release"
cmake_string+=" -DCMAKE_INSTALL_PREFIX=${base_dir}/moab"
cmake_string+=" -DCMAKE_INSTALL_RPATH=${plugin_dir}"

cmake ../src ${cmake_string}
make -j${jobs}
make install

# ------------------------------------------------------------------------------
# Build Plugin
# ------------------------------------------------------------------------------

cd ${base_dir}
git clone https://github.com/svalinn/DAGMC
git clone https://github.com/svalinn/mcnp2cad
rm -rf bld
mkdir -pv bld
cd bld

cmake_string=""
cmake_string+=" -DCUBIT_ROOT=${trelis_dir}"
cmake_string+=" -DMOAB_ROOT=${base_dir}/moab"
cmake_string+=" -DCMAKE_BUILD_TYPE=Release"
cmake_string+=" -DCMAKE_INSTALL_PREFIX=${base_dir}"
cmake_string+=" -DCMAKE_INSTALL_RPATH=${trelis_dir}:${plugin_dir}"

cmake .. ${cmake_string}
make -j${jobs}
make install

# ------------------------------------------------------------------------------
# Create tarball
# ------------------------------------------------------------------------------

cd ${base_dir}
rm -rf pack
mkdir -pv pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn

cp -pPv ${base_dir}/lib/* .
cp -pPv /usr/lib/libarmadillo.so.8* .
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
chmod 644 *

cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ..
ln -sv plugins/svalinn/libarmadillo.so.8 .
ln -sv plugins/svalinn/libarmadillo.so.8.400.0 .
ln -sv plugins/svalinn/libhdf5_serial.so.100 .
ln -sv plugins/svalinn/libhdf5_serial.so.100.0.1 .
ln -sv plugins/svalinn/libiGeom.so .
ln -sv plugins/svalinn/libmakeWatertight.so .
ln -sv plugins/svalinn/libmcnp2cad.so .
ln -sv plugins/svalinn/libMOAB.so .
ln -sv plugins/svalinn/libMOAB.so.5 .
ln -sv plugins/svalinn/libMOAB.so.5.1.0 .
ln -sv plugins/svalinn/libsvalinn_plugin.so .

cd ..
tar --sort=name -czvf ${tarball} bin
mv -v ${tarball} ..

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------

cd ${base_dir}
#rm -rf bld DAGMC lib mcnp2cad moab pack
