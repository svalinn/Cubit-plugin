#!/bin/bash

#DAGMC Installation
INSTALL_ROOT=$HOME
cd $INSTALL_ROOT
mkdir dagmc
cd dagmc
git clone https://github.com/svalinn/DAGMC
cd DAGMC
git checkout develop
cd ..
mkdir bld

DAGMC_INSTALL_DIR=$INSTALL_ROOT/dagmc/install
MOAB_CONFIG_DIR=$INSTALL_ROOT/MOAB/install/lib/cmake/MOAB
mkdir $DAGMC_INSTALL_DIR
cd bld
cmake ../DAGMC -DCMAKE_INSTALL_PREFIX=$DAGMC_INSTALL_DIR -DMOAB_CMAKE_CONFIG=$MOAB_CONFIG_DIR
make 
make install

printf '\nConsider running the following lines and adding them to .bashrc:\n'
printf "export PATH=$DAGMC_INSTALL_DIR/bin"':$PATH\n'
printf "export LD_LIBRARY_PATH=$DAGMC_INSTALL_DIR/lib"':$LD_LIBRARY_PATH\n'

cd $INSTALL_ROOT
