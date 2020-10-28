#!/bin/bash

cd /Trelis-sdk 
dpkg -i Trelis-$1-Lin64.deb

cd /opt
tar -xzvf /Trelis-sdk/Trelis-SDK-$1-Lin64.tar.gz
cd /opt/Trelis-16.5
tar -xzvf /Trelis-sdk/Trelis-SDK-$1-Lin64.tar.gz

apt-get update -y
apt-get install -y autogen autoconf libtool libeigen3-dev libhdf5-dev patchelf gfortran git cmake

cd 

# Setup
CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`

PLUGIN_DIR="plugin-build"

mkdir ${PLUGIN_DIR}
PLUGIN_ABS_PATH=${CURRENT}/${PLUGIN_DIR}

echo "Building the Trelis plugin in ${CURRENT}\\${PLUGIN_DIR}"

unset LD_LIBRARY_PATH

cd ${PLUGIN_ABS_PATH}
ln -s $SCRIPTPATH/ ./

mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
cd moab
autoreconf -fi
cd ../bld
../moab/configure CXXFLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
                  --disable-blaslapack \
                  --enable-shared \
                  --enable-optimize \
                  --disable-debug \
                  --disable-blaslapack \
                  --with-eigen3=/usr/include/eigen3 \
                  --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
                  --prefix=${PLUGIN_ABS_PATH}/moab
make -j`grep -c processor /proc/cpuinfo`
make install

cd ${PLUGIN_ABS_PATH}
mkdir -pv DAGMC/bld
cd DAGMC
git clone https://github.com/bam241/DAGMC -b preproc_plugin
cd bld
cmake ../DAGMC -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
               -DMOAB_DIR=${PLUGIN_ABS_PATH}/moab \
               -DBUILD_UWUW=ON \
               -DBUILD_TALLY=OFF \
               -DBUILD_BUILD_OBB=OFF \
               -DBUILD_MAKE_WATERTIGHT=ON \
               -DBUILD_SHARED_LIBS=ON \
               -DBUILD_STATIC_LIBS=OFF \
               -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/DAGMC
make -j`grep -c processor /proc/cpuinfo`
make install


cd ${PLUGIN_ABS_PATH}/Trelis-plugin
git submodule update --init

cd ${PLUGIN_ABS_PATH}
mkdir -pv bld
cd bld
ls /opt/Trelis-${1::4}
ls /opt/Trelis-*
ls /opt
cmake ../Trelis-plugin -DCUBIT_ROOT=/opt/Trelis-${1::4} \
                       -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
                       -DCMAKE_BUILD_TYPE=Release \
                       -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
make -j`grep -c processor /proc/cpuinfo`
make install
echo " cmake ../Trelis-plugin -DCubit_DIR=/opt/Trelis-${1::4} \
                       -DCUBIT_ROOT=/opt/Trelis-${1::4} \
                       -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
                       -DCMAKE_BUILD_TYPE=Release \
                       -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH} "

cd ${PLUGIN_ABS_PATH}
mkdir -p pack/bin/plugins/svalinn
cd pack/bin/plugins/svalinn

# Copy all needed libraries into current directory
cp -pPv ${PLUGIN_ABS_PATH}/lib/* .
cp -pPv ${PLUGIN_ABS_PATH}/moab/lib/libMOAB.so* .
cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libdagmc.so* .
cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libmakeWatertight.so* .
cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libpyne_dagmc.so .
cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libuwuw.so .
cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100* .
chmod 644 *

# Set the RPATH to be the current directory for the DAGMC libraries
patchelf --set-rpath /opt/Trelis-${1::4}/bin/plugins/svalinn libMOAB.so
patchelf --set-rpath /opt/Trelis-${1::4}/bin/plugins/svalinn libdagmc.so
patchelf --set-rpath /opt/Trelis-${1::4}/bin/plugins/svalinn libmakeWatertight.so
patchelf --set-rpath /opt/Trelis-${1::4}/bin/plugins/svalinn libpyne_dagmc.so
patchelf --set-rpath /opt/Trelis-${1::4}/bin/plugins/svalinn libuwuw.so

# Create the Svalinn plugin tarball
cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ../..
tar --sort=name -czvf svalinn-plugin_linux_$1.tgz bin
mv -v svalinn-plugin_linux_$1.tgz /Trelis-sdk
cd ..
rm -rf pack bld DAGMC lib moab
rm Trelis-plugin
