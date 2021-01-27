#!/bin/bash

cd /c/Users/mouginot/Desktop/moab_build
rm -rf *
cmake ../moab -DENABLE_BLASLAPACK=OFF -DENABLE_FORTRAN=OFF -DENABLE_IMESH=OFF -DENABLE_TESTING=OFF -DENABLE_HDF5=ON -DEIGEN3_DIR=/c/Users/mouginot/Desktop/eigen-3.3.8/  -G"Visual Studio 16 2019"  -DCMAKE_INSTALL_PREFIX=../moab_install/ -DHDF5_hdf5_LIBRARY_RELEASE="C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5_hl.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libzlib.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libszip.lib;C:/Program Files/HDF_Group/HDF5/1.8.21/lib/libhdf5_cpp.lib" -DCMAKE_EXE_LINKER_FLAGS="" -DCMAKE_MODULE_LINKER_FLAGS="" -DCMAKE_SHARED_LINKER_FLAGS="" -DCMAKE_STATIC_LINKER_FLAGS="" -DMOAB_BUILD_MBCONVERT=OFF -DMOAB_BUILD_HEXMODOPS=OFF -DMOAB_BUILD_MBSIZE=OFF -DMOAB_BUILD_MBMEM=OFF -DMOAB_BUILD_MBSKIN=OFF -DMOAB_BUILD_MBDEPTH=OFF -DMOAB_BUILD_MBTAGPROP=OFF -DMOAB_BUILD_MBGSETS=OFF -DMOAB_BUILD_SPHEREDECOMP=OFF -DMOAB_BUILD_MBSURFPLOT=OFF -DMOAB_BUILD_MBPART=OFF -DMOAB_BUILD_MBSLAVEPART=OFF -DMOAB_BUILD_MBCOUPLER=OFF -DMOAB_BUILD_MBHONODES=OFF -DMOAB_BUILD_MBUMR=OFF -DMOAB_BUILD_MBQUALITY=OFF -DMOAB_BUILD_MBTEMPEST=OFF -DCMAKE_BUILD_TYPE=Release
 
cmake --build . --config Release

cmake --install .


# add \"\" in the cmakeMOAB.config 
# copy MOAB.lib into MOAB.a ...
cd /c/Users/mouginot/Desktop/dagmc_build
rm -rf *

cmake ../dagmc -G"Visual Studio 16 2019" -DBUILD_EXE=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DBUILD_TALLY=OFF -DBUILD_BUILD_OBB=OFF -DBUILD_UWUW=ON -DBUILD_MAKE_WATERTIGHT=ON -DBUILD_TESTS=OFF -DMOAB_DIR=../moab_install -DCMAKE_PREFIX_PATH=../eigen-3.3.8/cmake  -DCMAKE_INSTALL_PREFIX=../dagmc_install/ -DCMAKE_EXE_LINKER_FLAGS="" -DCMAKE_MODULE_LINKER_FLAGS="" -DCMAKE_SHARED_LINKER_FLAGS="" -DCMAKE_STATIC_LINKER_FLAGS="" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
cmake --install . --config Release

cd /c/Users/mouginot/Desktop/plugin_build
rm -rf *

cmake ../Trelis-plugin/ -G"Visual Studio 16 2019" -DCubit_DIR="C:/Program Files/Trelis 17.1/bin" -DCUBIT_ROOT="C:/Program Files/Trelis 17.1/bin" -DDAGMC_DIR="../dagmc_install" -DEIGEN3_DIR=C:/Users/mouginot/Desktop/eigen-3.3.8 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../plugin_install
cmake --build . --config Release
cmake --install . --config Release
