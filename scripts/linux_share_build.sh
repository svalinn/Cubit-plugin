#!/bin/bash
PROC=$((`grep -c processor /proc/cpuinfo`))

function install_prerequisites() {
    TZ=America/Chicago
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
    apt-get update -y
    apt-get install -y g++ libeigen3-dev libhdf5-dev patchelf git cmake
}


function setup() {
    unset LD_LIBRARY_PATH
    
    echo "Building the Trelis plugin in ${CURRENT}\\${PLUGIN_DIR}"
    cd ${CURRENT}
    mkdir ${PLUGIN_DIR}
    cd ${PLUGIN_DIR}
    PLUGIN_ABS_PATH=$(pwd)
    ln -s ${SCRIPTPATH}/ ./

    # upload the variable in GithubAction
    if [ "$GHA" = "ON" ]; then
        echo "PLUGIN_ABS_PATH=$PLUGIN_ABS_PATH" >> $GITHUB_ENV
    fi
}

function setup_var() {
    # Setup the variables
    if [ "$1" = "2020.2" ]; then
        TRELIS_PATH="/opt/Coreform-Cubit-2020.2"
        TRELIS_PKG="Coreform-Cubit-2020.2-Lin64.deb"
        TRELIS_SDK_PKG="Coreform-Cubit-2020.2-Lin64-SDK.tar.gz"
    elif [ "$1" = "17.1.0" ]; then
        TRELIS_PATH="/opt/Trelis-17.1"
        TRELIS_PKG="Trelis-17.1.0-Lin64.deb"
        TRELIS_SDK_PKG="Trelis-SDK-17.1.0-Lin64.tar.gz"
        CMAKE_ADDITIONAL_FLAGS="-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0"
    else
        echo "unknown Trelis/Cubit version, use: \"17.1\" or \"2020.2\""
    fi
    echo  "GHA" $GHA

    # upload the variable in GithubAction
    if [ "$GHA" = "ON" ]; then
        echo  "TRELIS_PATH" $TRELIS_PATH
        echo "TRELIS_PATH=$TRELIS_PATH" >> $GITHUB_ENV
        echo "TRELIS_PKG=$TRELIS_PKG" >> $GITHUB_ENV
        echo "TRELIS_SDK_PKG=$TRELIS_SDK_PKG" >> $GITHUB_ENV
        echo "CMAKE_ADDITIONAL_FLAGS=$CMAKE_ADDITIONAL_FLAGS" >> $GITHUB_ENV
    fi
}


function build_moab() {
    cd ${PLUGIN_ABS_PATH}
    mkdir -pv moab/bld
    cd moab
    git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
    cd bld
    cmake ../moab -DENABLE_HDF5=ON \
            -DCMAKE_PREFIX_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial \
            -DBUILD_SHARED_LIBS=ON \
            -DENABLE_BLASLAPACK=OFF \
            -DENABLE_FORTRAN=OFF \
            $CMAKE_ADDITIONAL_FLAGS \
            -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/moab

    make -j`grep -c processor /proc/cpuinfo`
    make install
    cd ../..
    rm -rf moab/moab moab/bld
}


function build_dagmc(){
    cd ${PLUGIN_ABS_PATH}
    mkdir -pv DAGMC/bld
    cd DAGMC
    git clone https://github.com/svalinn/DAGMC -b develop
    cd bld
    cmake ../DAGMC -DMOAB_DIR=${PLUGIN_ABS_PATH}/moab \
                -DBUILD_UWUW=ON \
                -DBUILD_TALLY=OFF \
                -DBUILD_BUILD_OBB=OFF \
                -DBUILD_MAKE_WATERTIGHT=ON \
                -DBUILD_SHARED_LIBS=ON \
                -DBUILD_STATIC_LIBS=OFF \
                -DBUILD_EXE=OFF \
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
                -DCMAKE_BUILD_TYPE=Release \
                $CMAKE_ADDITIONAL_FLAGS \
                -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/DAGMC
                
    
    make -j`grep -c processor /proc/cpuinfo`
    make install
    cd ../..
    rm -rf DAGMC/DAGMC DAGCM/bld
}


function setup_Trelis_sdk() {

    cd ${FOLDER_PKG}
    apt install -y ./${TRELIS_PKG}
    cd /opt
    tar -xzvf ${FOLDER_PKG}/${TRELIS_SDK_PKG}
    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "$1" = "2020.2" ]; then
        cd ${TRELIS_PATH}/bin
        cp -pv CubitExport.cmake CubitExport.cmake.orig
        sed -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
        cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
        sed -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
    fi
}


function build_plugin(){
    cd ${PLUGIN_ABS_PATH}
    cd Trelis-plugin
    git submodule update --init
    cd ../
    mkdir -pv bld
    cd bld
    cmake ../Trelis-plugin -DCUBIT_ROOT=${TRELIS_PATH} \
                           -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
                           -DCMAKE_BUILD_TYPE=Release \
                            $CMAKE_ADDITIONAL_FLAGS \
                           -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
    make -j$PROC
    make install
}


function build_plugin_pkg(){
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
    cp -pPv /usr/lib/x86_64-linux-gnu/libhdf5_serial.so* .
    chmod 644 *

    # Set the RPATH to be the current directory for the DAGMC libraries
    patchelf --set-rpath ${TRELIS_PATH}/bin/plugins/svalinn libMOAB.so
    patchelf --set-rpath ${TRELIS_PATH}/bin/plugins/svalinn libdagmc.so
    patchelf --set-rpath ${TRELIS_PATH}/bin/plugins/svalinn libmakeWatertight.so
    patchelf --set-rpath ${TRELIS_PATH}/bin/plugins/svalinn libpyne_dagmc.so
    patchelf --set-rpath ${TRELIS_PATH}/bin/plugins/svalinn libuwuw.so

    # Create the Svalinn plugin tarball
    cd ..
    ln -sv svalinn/libsvalinn_plugin.so .
    cd ../..
    tar --sort=name -czvf svalinn-plugin_linux_$1.tgz bin
}
