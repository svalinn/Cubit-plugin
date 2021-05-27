#!/bin/bash
set -ex

# need to clear the LD_LIBRARY_PATH to avoid lib conflict
unset LD_LIBRARY_PATH

CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`
FOLDER_PKG="$2"

PLUGIN_DIR="plugin-build"

PLUGIN_ABS_PATH=""
CUBIT_PATH=""
CUBIT_PKG="$3"
CUBIT_SDK_PKG="$4"   
CMAKE_ADDITIONAL_FLAGS=""
UBUNTU_VERSION=""
HDF5_PATH=""
BUILD_SHARED_LIBS=""
SUDO=""
OS="UBUNTU"

source ${SCRIPTPATH}/scripts/unix_share_build.sh

if [ -f /.dockerenv ]; then
  install_prerequisites   
fi

setup 
setup_var $1

build_hdf5

build_moab
build_dagmc

# $1 is the version of Trelis/Cubit one are trying to compile against i.e. 17.1.0
linux_setup_cubit_sdk $1
build_plugin
linux_build_plugin_pkg $1

mv -v svalinn-plugin_linux_cubit_$1.tgz ${FOLDER_PKG}
cd ..
rm -rf pack bld DAGMC lib moab
rm Trelis-plugin
