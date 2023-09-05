#!/bin/bash

# THIS SET OF FUNCTION EXPECTS THE FOLLOWING VARIABLE TO BE DEFINED/AVAILABLE:
# - PLUGIN_ABS_PATH : absolute path to Cubit-plugin folder
# - CUBIT_PKG : name of the cubit installation package
# - CUBIT_SDK_PKG : name of the cubit SDK installation package (only used for version <= 1.7.0)
# - CUBIT_BASE_NAME : root name of the cubit installation
# - CUBIT_PATH : absolute path to the cubit installation folder
# - CMAKE_ADDITIONAL_FLAGS : additional flag required for the cmake setup (shared between plugin and dependencies)
# - BUILD_SHARED_LIBS : ON/OFF specifiies if shared libs will be built for MOAB and DAGMC
# - BUILD_STATIC_LIBS : ON/OFF specifiies if static libs will be built for MOAB and DAGMC
# - SUDO : specify the command prefix to run as root (if any)
# - SED : command to run sed (usually sed for linux, gsed for mac)
# - OS : reference to the Operating system, used to name the plugin tarball (and when using setup_var() function)
# - OS_VERSION : reference to the Operating system version, used to name the plugin tarball (and when using setup_var() function)

set -ex


function mac_install_brew() {
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
}

function mac_install_prerequisites() {
    brew install eigen gcc@6 gsed
}

function unix_version() {
    export LINUX_VERSION=$(lsb_release -rs |cut -d"." -f1)
    echo "Linux Version: " $LINUX_VERSION
}

function linux_install_prerequisites() {
    TZ=America/Chicago
    $SUDO ln -snf /usr/share/zoneinfo/$TZ /etc/localtime
    $SUDO sh -c 'echo $TZ > /etc/timezone'
    $SUDO apt-get update -y
    if [ "$OS" == "ubuntu" ]; then
    $SUDO apt-get install -y g++ libeigen3-dev patchelf git cmake curl lsb-release python3 lsb-core
    fi
    if [ "$OS" == "debian" ]; then
    $SUDO apt-get install -y g++ libeigen3-dev patchelf git cmake curl lsb-release python3
    fi
    $SUDO update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
    if [ $LINUX_VERSION -lt 21 ] && [ "$OS" == "ubuntu" ]; then
    $SUDO update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    fi
}

function setup() {
    unset LD_LIBRARY_PATH

    echo "Building the Cubit plugin in ${CURRENT}\\${PLUGIN_DIR}"
    cd ${CURRENT}
    mkdir ${PLUGIN_DIR}
    cd ${PLUGIN_DIR}
    PLUGIN_ABS_PATH=$(pwd)
    ln -s ${SCRIPTPATH}/ ./
    mkdir $SCRIPTPATH/release
}



function setup_var() {
    # Setup the variables
    CMAKE_ADDITIONAL_FLAGS="-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0"

    if [ "$1" == "2020.2" ]; then
        CUBIT_PATH="/opt/Coreform-Cubit-2020.2"
        unset CMAKE_ADDITIONAL_FLAGS
    elif [ "$1" == "17.1.0" ]; then
        CUBIT_PATH="/opt/Trelis-17.1"
    elif [ "$1" == "2021.3" ] ; then
        CUBIT_PATH="/opt/Coreform-Cubit-2021.3"
    elif [ "$1" == "2021.4" ] ; then
        CUBIT_PATH="/opt/Coreform-Cubit-2021.4"
    elif [ "$1" == "2022.4" ] ; then
        CUBIT_PATH="/opt/Coreform-Cubit-2022.4"
    elif [ "$1" == "2023.6-dev" ] ; then
        CUBIT_PATH="/opt/Coreform-Cubit-2023.6"
    else
        echo "unknown Cubit version"
        return 1
    fi

    if [ "$OS" == "MAC" ]; then
        BUILD_SHARED_LIBS="OFF"
        BUILD_STATIC_LIBS="ON"
    elif [ "$OS" == "ubuntu"]; then
        BUILD_SHARED_LIBS="ON"
        BUILD_STATIC_LIBS="OFF"
    elif [ "$OS" == "debian"]; then
        BUILD_SHARED_LIBS="ON"
        BUILD_STATIC_LIBS="OFF"
    else
        echo "OS ENV variable needs to be defined to either UBUNTU, DEBIAN or MAC"
        return 1
    fi
}

function linux_build_hdf5() {
    # if ubuntu 18.04 or lower rely on apt-get hdf5
    unix_version
    if [ $LINUX_VERSION -lt 20 ] && [ "$OS" == "ubuntu" ]; then
        $SUDO apt-get install -y libhdf5-serial-dev
    else
        cd ${PLUGIN_ABS_PATH}
        mkdir -p hdf5/bld
        cd hdf5
        git clone https://github.com/HDFGroup/hdf5.git -b hdf5-1_12_0 --depth 1 --shallow-submodules
        cd bld
        cmake ../hdf5 -DBUILD_SHARED_LIBS:BOOL=ON
        make
        $SUDO make install
    fi
}

function mac_build_hdf5() {
    brew install hdf5
    brew link hdf5
}

function build_moab() {
    cd ${PLUGIN_ABS_PATH}
    mkdir -pv moab/bld
    cd moab
    git clone https://bitbucket.org/fathomteam/moab -b 5.3.0 --depth 1 --shallow-submodules
    cd moab
    # patching MOAB CMakeLists.txt to use default find(HDF5)
    $SED -i "s/HDF5_MOAB/HDF5/" CMakeLists.txt
    cd ..
    #end of patch
    cd bld
    cmake ../moab -DENABLE_HDF5=ON \
            -DHDF5_ROOT=$HDF5_PATH \
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
            -DENABLE_BLASLAPACK=OFF \
            -DENABLE_FORTRAN=OFF \
            $CMAKE_ADDITIONAL_FLAGS \
            -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/moab
    make
    make install
    cd ../..
    rm -rf moab/moab moab/bld
}

function build_dagmc(){
    cd ${PLUGIN_ABS_PATH}
    mkdir -pv DAGMC/bld
    cd DAGMC
    git clone https://github.com/svalinn/DAGMC -b develop --depth 1 --shallow-submodules
    cd bld
    cmake ../DAGMC -DMOAB_DIR=${PLUGIN_ABS_PATH}/moab \
                -DBUILD_UWUW=ON \
                -DBUILD_TALLY=OFF \
                -DBUILD_BUILD_OBB=OFF \
                -DBUILD_MAKE_WATERTIGHT=ON \
                -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
                -DBUILD_STATIC_LIBS=${BUILD_STATIC_LIBS} \
                -DBUILD_EXE=OFF \
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
                -DCMAKE_BUILD_TYPE=Release \
                $CMAKE_ADDITIONAL_FLAGS \
                -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/DAGMC

    make
    make install
    cd ../..
    rm -rf DAGMC/DAGMC DAGMC/bld
}

function remove_app_logger() {
    cd ${CUBIT_PATH}/bin
    $SUDO cp -pv CubitExport.cmake CubitExport.cmake.orig
    $SUDO $SED -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
    $SUDO cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
    $SUDO $SED -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
}

function mac_setup_cubit () {
    if [ "${CUBIT_PKG##*.}" == "pkg" ]; then
        mac_pkg_setup_cubit $1
    else
        mac_dmg_setup_cubit $1
    fi
}
function mac_pkg_setup_cubit() {
    cd ${FOLDER_PKG}
    sudo installer -pkg ${CUBIT_PKG} -target /
}

function mac_dmg_setup_cubit() {
    cd ${FOLDER_PKG}
    hdiutil convert ${CUBIT_PKG} -format UDTO -o cubit_eula.dmg.cdr
    hdiutil attach cubit_eula.dmg.cdr -mountpoint /Volumes/Cubit
    mv /Volumes/Cubit/*.app /Applications/
    rm -rf cubit.dmg

    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "${1}" = "2020.2" ]; then #|| [ "$1" == "2021.3" ] || [ "$1" == "2021.4" ]; then
        remove_app_logger
    fi

    cd /Applications

    # 17.1.0 comes with a separate SDK. It is unclear yet on how it is supposed to be installed and used.
    # this is a way to have it working...
    if [ "$1" == "17.1.0" ] ; then
        cd ${CUBIT_PATH}
        $SUDO tar -xzf ${FOLDER_PKG}/${CUBIT_SDK_PKG}
        $SUDO rsync -a  ${CUBIT_BASE_NAME}/* ./
        $SUDO rsync -a  ${CUBIT_BASE_NAME}.app/Contents/MacOS/* MacOS/
        $SUDO rsync -a bin/* MacOS/
        $SUDO rm -rf bin ${CUBIT_BASE_NAME}.app
        $SUDO ln -s MacOS bin
        $SUDO ln -s ${CUBIT_PATH}/include /Applications/include

        #  # fixing the path to Contents/Include
        $SUDO cp -pv ${CUBIT_PATH}/MacOS/CubitExport-release.cmake ${CUBIT_PATH}/MacOS/CubitExport-release.cmake.orig
        $SUDO $SED -i "s/\${_IMPORT_PREFIX}/\/Applications/" ${CUBIT_PATH}/MacOS/CubitExport-release.cmake
    fi

    hdiutil detach /Volumes/Cubit
}

function linux_setup_cubit() {

    cd ${FOLDER_PKG}
    $SUDO apt-get install -y ./${CUBIT_PKG}

    if [ "$1" == "2023.4" ]; then
        cd ${CUBIT_PATH}
        mkdir license_server
        cd license_server
        ln -sf ../bin/libcf_license_server.so .
        ln -sf ../bin/libcf_license_renewals.so .
    fi

    if [ "$1" == "2021.3" ] || [ "$1" == "2021.4" ] || [ "$1" == "2021.5" ] || [ "$1" == "2021.11" ] || [ "$1" == "2022.4" ] || [ "$1" == "2023.6-dev" ] ; then
	    return
    fi

    cd /opt
    $SUDO tar -xzf ${FOLDER_PKG}/${CUBIT_SDK_PKG}
    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "$1" = "2020.2" ]; then
        remove_app_logger
    fi
}

function build_plugin(){
    cd ${PLUGIN_ABS_PATH}
    cd Cubit-plugin
    git config --global --add safe.directory ${GITHUB_WORKSPACE}
    git submodule update --init
    cd ../
    mkdir -pv bld
    cd bld
    cmake ../Cubit-plugin -DCMAKE_PREFIX_PATH=${CUBIT_PATH} \
                           -DCUBIT_ROOT=${CUBIT_PATH} \
                           -DCubit_DIR=${CUBIT_PATH} \
                           -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
                           -DCMAKE_BUILD_TYPE=Release \
                            $CMAKE_ADDITIONAL_FLAGS \
                           -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
    make -j$PROC
    make install
}

function linux_build_plugin_pkg(){
    cd ${PLUGIN_ABS_PATH}
    mkdir -p pack/bin/plugins/svalinn
    cd pack/bin/plugins/svalinn

    # Copy all needed libraries into current directory
    cp -pPv ${PLUGIN_ABS_PATH}/lib/* .
    cp -pPv ${PLUGIN_ABS_PATH}/moab/lib/libMOAB.so* .
    cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libdagmc.so* .
    cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libmakeWatertight.so* .
    cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libpyne_dagmc.so* .
    cp -pPv ${PLUGIN_ABS_PATH}/DAGMC/lib/libuwuw.so* .
    cp -vL $HDF5_PATH/lib/libhdf5.so* .
    chmod 644 *

    # Set the RPATH to be the current directory for the DAGMC libraries
    patchelf --set-rpath ${CUBIT_PATH}/bin/plugins/svalinn libMOAB.so
    patchelf --set-rpath ${CUBIT_PATH}/bin/plugins/svalinn libdagmc.so
    patchelf --set-rpath ${CUBIT_PATH}/bin/plugins/svalinn libmakeWatertight.so
    patchelf --set-rpath ${CUBIT_PATH}/bin/plugins/svalinn libpyne_dagmc.so
    patchelf --set-rpath ${CUBIT_PATH}/bin/plugins/svalinn libuwuw.so

    # Create the Svalinn plugin tarball
    cd ..
    ln -sv svalinn/libsvalinn_plugin.so .
    cd ../..
    PLUGIN_FILENAME=svalinn-plugin_${OS}-${OS_VERSION}_cubit_$1.tgz
    tar --sort=name -czvf ${PLUGIN_FILENAME} bin
    chmod 666 ${PLUGIN_FILENAME}
    cp ${PLUGIN_FILENAME} $SCRIPTPATH/release/
}

function mac_build_plugin_pkg(){
    cd ${PLUGIN_ABS_PATH}
    mkdir -p pack/MacOS/plugins/svalinn
    cd pack/MacOS/plugins/svalinn

    # Copy all needed libraries into current directory
    cp -pPv ${PLUGIN_ABS_PATH}/lib/* .
    cp /usr/local/opt/libaec/lib/libsz.dylib .
    install_name_tool -change /usr/local/opt/libaec/lib/libsz.dylib @rpath/libsz.dylib libsvalinn_plugin.so

    libsz=libsz.dylib
    if [ "$1" == "2022.4" ] || [ "$1" == "2023.6-dev" ] ; then
        libsz=libsz.2.dylib
    fi

    cp /usr/local/opt/libaec/lib/$libsz .
    install_name_tool -change /usr/local/opt/libaec/lib/$libsz @rpath/$libsz libsvalinn_plugin.so

    # restoring correct RPATH and BIN for 17.1 (bin does not exist as it is not shipped with SDK)
    if [ "$1" == "17.1.0" ] ; then
        # Correcting the RPATH for the svalinn and dependent libs
        rpath_old=${CUBIT_PATH}/bin/plugins/svalinn
        rpath_fix=${CUBIT_PATH}/MacOS/plugins/svalinn
        install_name_tool -rpath ${rpath_old} ${rpath_fix} libsvalinn_plugin.so
        install_name_tool -rpath ${rpath_old} ${rpath_fix} libiGeom.dylib
        install_name_tool -rpath ${rpath_old} ${rpath_fix} libmcnp2cad.dylib

        # Correcting the BIN for the svalinm and dependent libs
        bin_old=${CUBIT_PATH}/bin
        bin_fix=${CUBIT_PATH}/MacOS
        install_name_tool -rpath ${bin_old} ${bin_fix} libmcnp2cad.dylib
        install_name_tool -rpath ${bin_old} ${bin_fix} libiGeom.dylib
        install_name_tool -rpath ${bin_old} ${bin_fix} libsvalinn_plugin.so
    fi

    # Create the Svalinn plugin tarball
    cd ..
    ln -sv svalinn/libsvalinn_plugin.so .
    cd ../..
    PLUGIN_FILENAME=svalinn-plugin_${OS}_cubit_$1.tgz
    tar -czvf ${PLUGIN_FILENAME} MacOS
    chmod 666 ${PLUGIN_FILENAME}
    cp ${PLUGIN_FILENAME} $SCRIPTPATH/release/
}
