#!/bin/bash
set -ex

# need to clear the LD_LIBRARY_PATH to avoid lib conflict
unset LD_LIBRARY_PATH

CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`
FOLDER_PKG="${2}"

PLUGIN_DIR="plugin-build"

PLUGIN_ABS_PATH=""
TRELIS_PATH=""
TRELIS_PKG="$3"
TRELIS_SDK_PKG="$4"   
CMAKE_ADDITIONAL_FLAGS="" 

source ${SCRIPTPATH}/scripts/linux_share_build.sh

install_prerequisites

setup 
setup_var $1

build_hdf5

build_moab
build_dagmc

# $1 is the version of Trelis/Cubit one are trying to compile against i.e. 17.1.0
setup_Trelis_sdk $1
build_plugin $1
build_plugin_pkg $1

mv -v svalinn-plugin_linux_cubit_$1.tgz ${FOLDER_PKG}
cd ..
rm -rf pack bld DAGMC lib moab
rm Trelis-plugin
