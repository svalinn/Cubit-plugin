#!/bin/bash
CURRENT=$(pwd)
SCRIPTPATH=`dirname $(dirname $(realpath $0))`

PLUGIN_DIR="plugin-build"
PLUGIN_ABS_PATH=""

source ${SCRIPTPATH}/script/linux_share_build.sh

install_prerequise

setup_folder

build_moab
build_dagmc

setup_Trelis_sdk $1
build_plugin $1
build_plugin_pkg $1



mv -v svalinn-plugin_linux_$1.tgz /Trelis-sdk
cd ..
rm -rf pack bld DAGMC lib moab
rm Trelis-plugin
