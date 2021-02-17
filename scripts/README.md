SVALINN PLUGIN SCRIPTS
======================

This folder contains scripts to build the Svalinn plugin for multiple platform.

Linux (tested on Ubuntu 18.04 and 20.04):
-----------------------------------------

- `linux_share_build.sh`: contains multiple methods allowing to setup and build the required dependencies as well as the plugin
    - `install_prerequisites()`: apt-get install system dependencies
    - `setup_folder()`: build the fodler layout required to plugin build
    - `build_moab()`: build and install MOAB
    - `build_dagmc()`: build and install DAGMC
    - `setup_Trelis_sdk()`: install Trelis/CUbit and the corresponding SDK
    - `build_plugin()`: build the plugin
    - `build_plugin_pkg()`:  build a tarball containing the plugin

- `build_plugin_linux.sh`: build the plugin on a linux computer (relies on apt-get to install required dependencies)
    - arguments: Trelis/CUbit version to build the plugin for (only supported yet `17.1.0`)
- `build_plugin_docker_ubuntu.sh`: build the plugin in a docker container 
    - arguments: Arguments have to be provided in order
        - 1 docker image (tested on `ubuntu:18.04` and `ubuntu:20.04`)
        - 2 local path to the folder containing the Trelis/Cubit install deb and SDK
        - 3 version of CUbit Trelis to build the plugin for (only supported yet `17.1.0`)
