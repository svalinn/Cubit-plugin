#!/bin/bash
set -ex


function mac_install_brew() {
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
}

function mac_install_prerequisites() {
    brew install eigen gcc@6 gsed
}

function ubuntu_version() {
    export UBUNTU_VERSION=$(lsb_release -rs |cut -d"." -f1)
    echo "Ubuntu Version: " $UBUNTU_VERSION
}

function linux_install_prerequisites() {
    TZ=America/Chicago
    $SUDO ln -snf /usr/share/zoneinfo/$TZ /etc/localtime
    $SUDO sh -c 'echo $TZ > /etc/timezone'
    $SUDO apt-get update -y
    $SUDO apt-get install -y g++ libeigen3-dev patchelf git cmake curl lsb-release python3  lsb-core
    $SUDO update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
    $SUDO update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
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
    CMAKE_ADDITIONAL_FLAGS="-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0"

    if [ "$1" == "2020.2" ]; then
        TRELIS_PATH="/opt/Coreform-Cubit-2020.2"
        unset CMAKE_ADDITIONAL_FLAGS
    elif [ "$1" == "17.1.0" ]; then
        TRELIS_PATH="/opt/Trelis-17.1"
    elif [ "$1" == "2021.3" ] ; then
        TRELIS_PATH="/opt/Coreform-Cubit-2021.3"
    elif [ "$1" == "2021.4" ] ; then
        TRELIS_PATH="/opt/Coreform-Cubit-2021.4"
    else
        echo "unknown Trelis/Cubit version, use: \"17.1.0\" or \"2020.2\""
        return 1
    fi

    if [ "OS" == "MAC" ]; then
        BUILD_SHARED_LIBS="OFF"
        BUILD_STATIC_LIBS="ON"
    elif [ "OS" == "UBUNTU"]; then
        BUILD_SHARED_LIBS="ON"
        BUILD_STATIC_LIBS="OFF"
    else
        echo "OS ENV variable needs to be defined to either UBUNTU or MAC"
        return 1
    fi
}

function linux_build_hdf5() {
    # if ubuntu 18.04 or lower rely of apt-get hdf5
    ubuntu_version
    if [ $UBUNTU_VERSION -lt 20 ]; then
        $SUDO apt-get install -y libhdf5-serial-dev
    else
        cd ${PLUGIN_ABS_PATH}
        mkdir -p hdf5/bld
        cd hdf5
        git clone https://github.com/HDFGroup/hdf5.git -b hdf5-1_12_0
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
    git clone https://bitbucket.org/fathomteam/moab -b Version5.1.0
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
    git clone https://github.com/svalinn/DAGMC -b develop
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

function mac_setup_cubit_sdk() {
    


    cd ${FOLDER_PKG}
    if [ "${1}" == "17.1.0" ]; then
        hdiutil convert ${CUBIT_PKG} -format UDTO -o trelis_eula.dmg.cdr
        hdiutil attach trelis_eula.dmg.cdr -mountpoint /Volumes/Cubit
        mv /Volumes/Cubit/*.app /Applications/
        hdiutil detach /Volumes/Cubit
        rm -rf trelis.dmg
    elif [ "${1}" == "2020.2" ]; then
        sudo installer -pkg ${CUBIT_PKG} -target /
        rm -rf cubit.pkg
    fi
    if [ "$1" == "2021.3" ] || [ "$1" == "2021.4" ] ; then
	    return
    fi
    cd ${CUBIT_PATH}

    sudo tar -xzf ${FOLDER_PKG}/${CUBIT_SDK_PKG}
    sudo mv ${CUBIT_BASE_NAME}/* ./
    sudo mv ${CUBIT_BASE_NAME}.app/Contents/MacOS/* MacOS/
    sudo mv bin/* MacOS/
    sudo rm -rf bin ${CUBIT_BASE_NAME}.app
    sudo ln -s MacOS bin
    sudo ln -s ${CUBIT_PATH}/include /Applications/include


    # fixing the path to Contents/Include
    sudo cp -pv ${CUBIT_PATH}/MacOS/CubitExport-release.cmake ${CUBIT_PATH}/MacOS/CubitExport-release.cmake.orig
    sudo $SED -i "s/\/${CUBIT_BASE_NAME}.app\/Contents//" ${CUBIT_PATH}/MacOS/CubitExport-release.cmake

    sudo cp -pv ${CUBIT_PATH}/MacOS/CubitGeomConfig.cmake ${CUBIT_PATH}/MacOS/CubitGeomConfig.cmake.orig
    sudo $SED -i "s/\${_IMPORT_PREFIX}\/include/\${_IMPORT_PREFIX}\/${CUBIT_BASE_NAME}.app\/include/" ${CUBIT_PATH}/MacOS/CubitGeomConfig.cmake


    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "${1}" = "2020.2" ]; then
        cd ${CUBIT_PATH}/bin
        sudo cp -pv CubitExport.cmake CubitExport.cmake.orig
        sudo $SED -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
        sudo cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
        sudo $SED -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
    fi

}

function linux_setup_cubit_sdk() {

    if [ "$1" == "2021.3" ] || [ "$1" == "2021.4" ] ; then
	    return
    fi
    
    cd ${FOLDER_PKG}
    $SUDO apt-get install -y ./${CUBIT_PKG}
    cd /opt
    $SUDO tar -xzf ${FOLDER_PKG}/${CUBIT_SDK_PKG}
    # removing app_loger that seems to not be present in Cubit 2020.2
    if [ "$1" = "2020.2" ]; then
        cd ${CUBIT_PATH}/bin
        $SUDO cp -pv CubitExport.cmake CubitExport.cmake.orig
        $SUDO $SED -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
        $SUDO cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
        $SUDO $SED -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
    fi
}

function build_plugin(){
    cd ${PLUGIN_ABS_PATH}
    cd Trelis-plugin
    git submodule update --init
    cd ../
    mkdir -pv bld
    cd bld
    ls ${CUBIT_PATH}
    cmake ../Trelis-plugin -DCUBIT_ROOT=${CUBIT_PATH} \
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
    tar --sort=name -czvf svalinn-plugin_linux_cubit_$1.tgz bin
    chmod 666 svalinn-plugin_linux_cubit_$1.tgz
}

function mac_build_plugin_pkg(){
    cd ${PLUGIN_ABS_PATH}
    mkdir -p pack/MacOS/plugins/svalinn
    cd pack/MacOS/plugins/svalinn

    # Copy all needed libraries into current directory
    cp -pPv ${PLUGIN_ABS_PATH}/lib/* .
    cp /usr/local/opt/szip/lib/libsz.2.dylib .
    install_name_tool -change /usr/local/opt/szip/lib/libsz.2.dylib @rpath/libsz.2.dylib libsvalinn_plugin.so

    # Create the Svalinn plugin tarball
    cd ..
    ln -sv svalinn/libsvalinn_plugin.so .
    cd ../..
    tar -czvf svalinn-plugin_mac_cubit_${1}.tgz MacOS
    chmod 666 svalinn-plugin_mac_cubit_$1.tgz
}