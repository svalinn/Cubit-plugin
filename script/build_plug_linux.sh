#!/bin/bash

SCRIPTPATH=`dirname $(dirname $(realpath $0))`
docker run -v "$SCRIPTPATH:/Trelis-plugin" -v "$2:/Trelis-sdk" -it $1 bash -c "/Trelis-plugin/script/build_plugin.sh $3"