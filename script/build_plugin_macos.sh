#!/bin/bash

#Getting Macport
#curl https://distfiles.macports.org/MacPorts/MacPorts-2.6.3-10.15-Catalina.pkg --output MacPorts-2.6.3-10.15-Catalina.pkg
#sudo installer -pkg MacPorts-2.6.3-10.15-Catalina.pkg -target /
#rm -rf  MacPorts-2.6.3-10.15-Catalina.pkg

export PATH=/opt/local/bin:/opt/local/sbin:$PATH
export MANPATH=/opt/local/share/man:$MANPATH
export LD_LIBRARY_PATH=/opt/local/lib:$LD_LIBRARY_PATH


#sudo port -N selfupdate
#sudo port -N install libtool eigen3 hdf5 cmake gcc6 wget realpath

#wget https://github.com/fxcoudert/gfortran-for-macOS/releases/download/10.2/gfortran-10.2-Catalina.dmg
#hdiutil attach gfortran-10.2-Catalina.dmg
#sudo installer -pkg /Volumes/gfortran-10.2-Catalina/gfortran.pkg -target /
#hdiutil detach /Volumes/gfortran-10.2-Catalina


cd 

# Setup
CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`

PLUGIN_DIR="plugin-build"
mkdir ${PLUGIN_DIR}
hdiutil attach -quiet -nobrowse -noverify -noautoopen SDK/Trelis-17.1.0-Mac64.dmg
cp -rf /Volumes/Trelis-17.1.0-Mac64/Trelis-17.1.app ${PLUGIN_DIR}/
hdiutil detach /Volumes/Trelis-17.1.0-Mac64

PLUGIN_ABS_PATH=${CURRENT}/${PLUGIN_DIR}

echo "Building the Trelis plugin in ${CURRENT}\\${PLUGIN_DIR}"

unset LD_LIBRARY_PATH

cd ${PLUGIN_ABS_PATH}
ln -s $SCRIPTPATH/ ./

mkdir -pv moab/bld
cd moab
git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
cd bld
cmake ../moab -DENABLE_HDF5=ON \
          -DCMAKE_PREFIX_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial \
          -DBUILD_SHARED_LIBS=OFF \
          -DENABLE_BLASLAPACK=OFF -DENABLE_FORTRAN=OFF \
          -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
          -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/moab
make -j
make install

cd ${PLUGIN_ABS_PATH}
mkdir -pv DAGMC/bld
cd DAGMC
git clone https://github.com/bam241/DAGMC -b build_exe
cd bld
cmake ../DAGMC -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 \
               -DMOAB_DIR=${PLUGIN_ABS_PATH}/moab \
               -DBUILD_UWUW=ON \
               -DBUILD_TALLY=OFF \
               -DBUILD_BUILD_OBB=OFF \
               -DBUILD_MAKE_WATERTIGHT=ON \
               -DBUILD_SHARED_LIBS=OFF \
               -DBUILD_STATIC_LIBS=ON \
               -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
               -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/DAGMC
make -j
make install

cd
TRELIS_INSTALL_LOC="${PLUGIN_ABS_PATH}/Trelis-17.1.app/Contents"
cd $TRELIS_INSTALL_LOC 
tar -xzf /Users/mouginot/SDK/Trelis-SDK-17.1.0-Mac64.tar
mv Trelis-17.1/* ./
cp -f Trelis-17.1.app/Contents/MacOS/* MacOS/
cp -f bin/* MacOS/
rm -rf bin Trelis-17.1.app
ln -s MacOS bin
#cd bin
#sudo cp -pv CubitExport-Release.cmake CubitExport-Release.cmake.orig
#sudo port install gsed
#sudo gsed -i "s/\"Trelis-17.1.app\/Contents/\MacOS\"/\"bin\"/" CubitExport-Release.cmake
cd 

cd ${PLUGIN_ABS_PATH}/Trelis-plugin
git submodule update --init

cd ${PLUGIN_ABS_PATH}/Trelis-plugin
cp script/CubitExport-release.cmake ${TRELIS_INSTALL_LOC}/MacOS/
rm -rf bld
mkdir -pv bld
cd bld
cmake .. -DCubit_DIR=${TRELIS_INSTALL_LOC}/MacOS \
         -DCUBIT_ROOT=${TRELIS_INSTALL_LOC}/MacOS \
         -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
make -j
make install


#                       -DCMAKE_BUILD_TYPE=Release \



cd ${PLUGIN_ABS_PATH}
mkdir -p pack/MacOS/plugins/svalinn
cd pack/MacOS/plugins/svalinn

# Copy all needed libraries into current directory
cp -pPv ${PLUGIN_ABS_PATH}/lib/* .

# Create the Svalinn plugin tarball
cd ..
ln -sv svalinn/libsvalinn_plugin.so .
cd ../..
tar -czvf svalinn-plugin_mac.tgz MacOS
mv -v svalinn-plugin_mac.tgz ~/
cd ..
# rm -rf pack bld DAGMC lib moab
# rm Trelis-plugin
