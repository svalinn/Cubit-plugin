#!/bin/bash
CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`

PLUGIN_DIR="plugin-build"
PLUGIN_ABS_PATH=""

FOLDER_PKG="/trelis-sdk"
TRELIS_PATH=""
TRELIS_PKG=""
TRELIS_SDK_PKG=""    

source ${SCRIPTPATH}/script/linux_share_build.sh

install_prerequise
setup_Trelis_sdk $1
setup_folder


build_moab
build_dagmc

build_plugin
build_plugin_pkg $1



mv -v svalinn-plugin_linux_$1.tgz ${FOLDER_PKG}
cd ..
rm -rf pack bld DAGMC lib moab
rm Trelis-plugin
