#!/bin/bash


# Usage: 
# $1: docker image to use (tested with ubuntu:18.04 and ubuntu:20.04)
# $2: Trelis/Cubit version to install, 17.1.0 or 2020.2
# $3 path to local pkg folder containing both trelis/cubit install package and sdk
# $4 name of the Trelis deb package
# $5 name of the Trelis sdk tarball

SCRIPTPATH=`dirname $(dirname $(realpath $0))`
docker run -v "$SCRIPTPATH:/Trelis-plugin" -v "$3:/trelis-sdk" -it $1 bash -c "/Trelis-plugin/scripts/build_plugin_linux.sh $2 /trelis-sdk $4 $5; bash"
