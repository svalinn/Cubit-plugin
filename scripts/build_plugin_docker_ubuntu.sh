#!/bin/bash


# Usage: 
# $1: docker image to use (tested with ubuntu:18.04 and ubuntu:20.04)
# $2 path to local pkg folder containing both trelis/cubit install package and sdk
# $3: Trelis/Cubit version to install, 17.1.0 or 2020.2

SCRIPTPATH=`dirname $(dirname $(realpath $0))`
docker run -v "$SCRIPTPATH:/Trelis-plugin" -v "$2:/trelis-sdk" -it $1 bash -c "/Trelis-plugin/scripts/build_plugin_linux.sh $3 /trelis-sdk; bash"
