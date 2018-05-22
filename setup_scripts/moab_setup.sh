#!/bin/bash
 
#MOAB Installation
INSTALL_ROUTE=$HOME
cd $INSTALL_ROUTE

apt install -y autoconf libtool make mpich libblas-dev liblapack-dev libhdf5-dev cmake

mkdir MOAB
cd MOAB
mkdir bld
mkdir install
git clone https://bitbucket.org/fathomteam/moab
cd moab
git checkout Version5.0
autoreconf -fi
cd ../bld
../moab/configure --prefix=$INSTALL_ROUTE/MOAB/install --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial --enable-shared;
make -j4
make check
make install
printf 'export PATH=$PATH:%s/MOAB/install/bin' "$INSTALL_ROUTE">>$HOME/.bashrc
export PATH=$PATH:$INSTALL_ROUTE/MOAB/install/bin
printf 'export LD_LIBRARY_PATH=%s/MOAB/install/lib' "$INSTALL_ROUTE">>$HOME/.bashrc
export LD_LIBRARY_PATH=$INSTALL_ROUTE/MOAB/install/lib
