SVALINN PLUGIN SCRIPTS
======================

This folder contains scripts to build the Svalinn plugin on multiple platforms.

Linux (tested on Ubuntu 18.04 and 20.04):
-----------------------------------------

- `linux_share_build.sh`: contains multiple methods allowing to setup and build the required dependencies as well as the plugin
    - `linux/mac_install_prerequisites()`: apt-get install system dependencies as function of the operating system
    - `setup_folder()`: build the folder layout required to plugin build
    - `linux/mac__build_hdf5`: install/build hdf5 as function of the operating system
    - `build_moab()`: build and install MOAB
    - `build_dagmc()`: build and install DAGMC
    - `linux/mac_setup_cubit()`: install Trelis/CUbit and the corresponding SDK as function of the operating system
    - `build_plugin()`: build the plugin
    - `linux/mac_build_plugin_pkg()`: build a tarball containing the plugin as function of the operating system

- `linux/mac/build_plugin_linux.sh`: build the plugin on a linux computer (relies on apt-get to install required dependencies)
    - arguments: 
        - 1 Cubit version to build the plugin for (only supported yet `17.1.0`and Cubit `2020.2`)
        - 2 path to Cubit package folder (containing Cubit install package and SDK)
        - 3 name of the Cubit deb package
        - 4 name of the Cubit sdk tarball
- `build_plugin_docker_ubuntu.sh`: build the plugin in a docker container 
    - arguments: Arguments have to be provided in order
        - 1 docker image (tested on `ubuntu:18.04` and `ubuntu:20.04`)
        - 2 version of Cubit Cubit to build the plugin for (only tested with `17.1.0` and Cubit `2020.2`)
        - 3 local path to the folder containing the Cubit install deb and SDK
        - 4 name of the Cubit deb package
        - 5 name of the Cubit sdk tarball
