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
mkdir $DAGMC_INSTALL_DIR
cd bld
cmake ../DAGMC -DCMAKE_INSTALL_PREFIX=$DAGMC_INSTALL_DIR
make 
make install
printf '\nConsider adding the following lines to .bashrc:\n'
printf 'export PATH=$DAGMC_INSTALL_DIR/bin:$PATH\n'
export PATH=$DAGMC_INSTALL_DIR/bin:$PATH
printf 'export LD_LIBRARY_PATH=$DAGMC_INSTALL_DIR/lib:$LD_LIBRARY_PATH\n'
export LD_LIBRARY_PATH=$DAGMC_INSTALL_DIR/lib:$LD_LIBRARY_PATH
