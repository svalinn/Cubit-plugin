#!/bin/bash
cd /c/Users/mouginot/Desktop
mkdir plugin_build
cd plugin_build

git clone https://bitbucket.org/bam241/moab -b windows
mkdir moab_build moab_install
cd moab_build
cmake ../moab -DENABLE_BLASLAPACK=OFF -DENABLE_FORTRAN=OFF -DENABLE_IMESH=OFF -DENABLE_TESTING=OFF -DENABLE_HDF5=ON -DEIGEN3_DIR=/c/Users/mouginot/Desktop/eigen-3.3.8/  -G"Visual Studio 16 2019"  -DCMAKE_INSTALL_PREFIX=../moab_install/ -DHDF5_hdf5_LIBRARY_RELEASE="C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5_hl.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libzlib.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libszip.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5_cpp.lib" -DCMAKE_EXE_LINKER_FLAGS="" -DCMAKE_MODULE_LINKER_FLAGS="" -DCMAKE_SHARED_LINKER_FLAGS="" -DCMAKE_STATIC_LINKER_FLAGS="" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.27.29110/bin/Hostx64/x64/cl.exe" -DCMAKE_CXX_COMPILER="C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.27.29110/bin/Hostx64/x64/cl.exe" -DBUILD_SHARED_LIBS=ON
 
#-DMOAB_BUILD_HEXMODOPS=OFF -DMOAB_BUILD_MBCONVERT=OFF -DMOAB_BUILD_MBCOUPLER=OFF -DMOAB_BUILD_MBDEPTH=OFF -DMOAB_BUILD_MBGSETS=OFF -DMOAB_BUILD_MBHONODES=OFF -DMOAB_BUILD_MBMEM=OFF -DMOAB_BUILD_MBPART=OFF -DMOAB_BUILD_MBQUALITY=OFF -DMOAB_BUILD_MBSIZE=OFF -DMOAB_BUILD_MBSKIN=OFF -DMOAB_BUILD_MBSLAVEPART=OFF -DMOAB_BUILD_MBSURFPLOT=OFF -DMOAB_BUILD_MBTAGPROP=OFF -DMOAB_BUILD_MBTEMPEST=OFF -DMOAB_BUILD_MBUMR=OFF -DMOAB_BUILD_SPHEREDECOMP=OFF 

cmake --build . --config Release
cmake --install . --config Release


# add \"\" in the cmakeMOAB.config 
# copy MOAB.lib into MOAB.a ...
cd /c/Users/mouginot/Desktop/plugin_build
mkdir dagmc_build dagmc_install

git clone https://github.com/bam241/dagmc -b windows
cd dagmc_build
cmake ../dagmc -G"Visual Studio 16 2019" -DBUILD_EXE=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DBUILD_TALLY=OFF -DBUILD_BUILD_OBB=OFF -DBUILD_UWUW=ON -DBUILD_MAKE_WATERTIGHT=ON -DBUILD_TESTS=OFF -DMOAB_DIR=../moab_install -DCMAKE_PREFIX_PATH="C:/Users/mouginot/Desktop/eigen-3.3.8/cmake"  -DCMAKE_INSTALL_PREFIX=../dagmc_install/ -DCMAKE_EXE_LINKER_FLAGS="" -DCMAKE_MODULE_LINKER_FLAGS="" -DCMAKE_SHARED_LINKER_FLAGS="" -DCMAKE_STATIC_LINKER_FLAGS="" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
cmake --install . --config Release

cd /c/Users/mouginot/Desktop/plugin_build
cp -r ../Trelis-plugin ./
cd Trelis-plugin
#git submodule update --init
cd ..
mkdir plugin_build plugin_install
cd plugin_build

cmake ../Trelis-plugin/ -G"Visual Studio 16 2019" -DCubit_DIR="C:/Program Files/Trelis 17.1/bin" -DCUBIT_ROOT="C:/Program Files/Trelis 17.1/bin" -DDAGMC_DIR="../dagmc_install" -DEIGEN3_DIR=C:/Users/mouginot/Desktop/eigen-3.3.8 -DEigen3_DIR=C:/Users/mouginot/Desktop/eigen-3.3.8  -DCMAKE_PREFIX_PATH="C:/Users/mouginot/Desktop/eigen-3.3.8/cmake"   -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../plugin_install  -DCMAKE_EXE_LINKER_FLAGS="" -DCMAKE_MODULE_LINKER_FLAGS="" -DCMAKE_SHARED_LINKER_FLAGS="" -DCMAKE_STATIC_LINKER_FLAGS="" 
cmake --build . --config Release
cmake --install . --config Release

cd /c/Users/mouginot/Desktop/plugin_build
mkdir -p bin/plugins
cd bin/plugins
cp ../../moab_install/bin/MOAB.dll ./
cp ../../dagmc_install/lib/*.lib ./
cp ../../plugin_install/bin/* ./
cp ../../plugin_install/lib/* ./
cd /c/Users/mouginot/Desktop/plugin_build
tar -cvf svalin_plugin_windows.tar  bin

mv svalin_plugin_windows.tar ../
cd ..
#rm -rf plugin_build