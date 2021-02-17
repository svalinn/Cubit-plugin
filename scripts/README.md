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