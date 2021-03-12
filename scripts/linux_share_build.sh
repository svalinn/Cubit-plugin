#!/bin/bash
PROC=$((`grep -c processor /proc/cpuinfo`))

function install_prerequisites() {
    TZ=America/Chicago
    $SUDO ln -snf /usr/share/zoneinfo/$TZ /etc/localtime
    $SUDO sh -c 'echo $TZ > /etc/timezone'
    $SUDO apt-get update -y
    $SUDO apt-get install -y g++ libeigen3-dev patchelf git cmake curl
}

function setup() {
    unset LD_LIBRARY_PATH
    
    echo "Building the Trelis plugin in ${CURRENT}\\${PLUGIN_DIR}"
    cd ${CURRENT}
    mkdir ${PLUGIN_DIR}
    cd ${PLUGIN_DIR}
    PLUGIN_ABS_PATH=$(pwd)
    ln -s ${SCRIPTPATH}/ ./
}

function setup_var() {
    # Setup the variables
    if [ "$1" = "2020.2" ]; then
        TRELIS_PATH="/opt/Coreform-Cubit-2020.2"
    elif [ "$1" = "17.1.0" ]; then
        TRELIS_PATH="/opt/Trelis-17.1"
        CMAKE_ADDITIONAL_FLAGS="-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0"
    else
        echo "unknown Trelis/Cubit version, use: \"17.1.0\" or \"2020.2\""
        return 1
        
    fi

}

function build_hdf5() {
    cd ${PLUGIN_ABS_PATH}
    mkdir -p hdf5/bld
    cd hdf5
    git clone https://github.com/HDFGroup/hdf5.git -b hdf5-1_12_0
    cd bld
    cmake ../hdf5 -DBUILD_SHARED_LIBS:BOOL=ON
    make -j$PROC
    $SUDO make install
}

function build_moab() {
    cd ${PLUGIN_ABS_PATH}
    mkdir -pv moab/bld
    cd moab
    git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
    cd moab
    # patching MOAB CMakeLists.txt to use default find(HDF5)
    sed -i "s/HDF5_MOAB/HDF5/" CMakeLists.txt
    cd ..
    #end of patch
    cd bld
    cmake ../moab -DENABLE_HDF5=ON \
            -DHDF5_ROOT=/usr/local/HDF_Group/HDF5/1.12.0 \
            -DBUILD_SHARED_LIBS=ON \
            -DENABLE_BLASLAPACK=OFF \
            -DENABLE_FORTRAN=OFF \
            $CMAKE_ADDITIONAL_FLAGS \
            -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}/moab

    make -j$PROC
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
    
    make -j$PROC
    make install
    cd ../..
    rm -rf DAGMC/DAGMC DAGMC/bld
}

function setup_Trelis_sdk() {

    cd ${FOLDER_PKG}
    $SUDO apt-get install -y ./${TRELIS_PKG}
    cd /opt
    $SUDO tar -xzvf ${FOLDER_PKG}/${TRELIS_SDK_PKG}
    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "$1" = "2020.2" ]; then
        cd ${TRELIS_PATH}/bin
        $SUDO cp -pv CubitExport.cmake CubitExport.cmake.orig
        $SUDO sed -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
        $SUDO cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
        $SUDO sed -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
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
    cp -pPv /usr/local/HDF_Group/HDF5/1.12.0/lib/libhdf5.so* .
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
    tar --sort=name -czvf svalinn-plugin_linux_cubit_$1.tgz bin
}
