#!/bin/bash

function install_brew() {
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
}

function install_prerequisites() {
    brew install eigen hdf5 gcc@6 gsed
    brew link hdf5
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
        TRELIS_PATH="/Applications/Coreform-Cubit-2020.2/Contents"
    elif [ "$1" = "17.1.0" ]; then
        TRELIS_PATH="/Applications/Trelis-17.1.app/Contents"
        CMAKE_ADDITIONAL_FLAGS="-DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0"
    else
        echo "unknown Trelis/Cubit version, use: \"17.1.0\" or \"2020.2\""
        return 1
    fi

}

function build_hdf5() {
    # if ubuntu 18.04 or lower rely of apt-get hdf5
    if [[ $UBUNTU_VERSION < 20 ]]; then
        $SUDO apt-get install -y libhdf5-serial-dev
        HDF5_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial"
    else
        cd ${PLUGIN_ABS_PATH}
        mkdir -p hdf5/bld
        cd hdf5
        git clone https://github.com/HDFGroup/hdf5.git -b hdf5-1_12_0
        cd bld
        cmake ../hdf5 -DBUILD_SHARED_LIBS:BOOL=ON
        make
        $SUDO make install
        HDF5_PATH="/usr/local/HDF_Group/HDF5/1.12.0"
    fi
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
            -DHDF5_ROOT=$HDF5_PATH \
            -DBUILD_SHARED_LIBS=OFF \
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
                -DBUILD_SHARED_LIBS=OFF \
                -DBUILD_STATIC_LIBS=ON \
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

function setup_Trelis {

    cd ${PLUGIN_ABS_PATH}
    if [ "${1}" = "17_1_0" ]; then
        hdiutil convert trelis.dmg -format UDTO -o ${TRELIS_PKG}
        hdiutil attach trelis_eula.dmg.cdr -mountpoint /Volumes/Cubit
        mv /Volumes/Cubit/*.app /Applications/
        hdiutil detach /Volumes/Cubit
        rm -rf trelis.dmg
    elif [ "${1}" = "2020_2" ]; then
        sudo installer -pkg ${TRELIS_PKG} -target /
        rm -rf cubit.pkg
    fi

}

function setup_trelis_sdk() {
   
    cd ${CUBIT_PATH}
    if [ "${1}" = "2020.2" ]; then
        CUBIT_BASE_NAME="Coreform-Cubit-2020.2"
    elif [ "${1}" = "17.1.0" ]; then
        CUBIT_BASE_NAME="Trelis-17.1"
    fi
    ls -al
    sudo tar -xzf ${FOLDER_PKG}/${TRELIS_SDK_PKG}
    echo "ARG 1: ${1}"
    echo "CUBIT_BASE_NAME: ${CUBIT_BASE_NAME}"
    sudo mv ${CUBIT_BASE_NAME}/* ./
    sudo mv ${CUBIT_BASE_NAME}.app/Contents/MacOS/* MacOS/
    sudo mv bin/* MacOS/
    sudo rm -rf bin ${CUBIT_BASE_NAME}.app
    sudo ln -s MacOS bin
    sudo ln -s ${TRELIS_PATH}/include /Applications/include

    cd ${PLUGIN_ABS_PATH}/Trelis-plugin
    sudo cp scripts/*.cmake ${TRELIS_PATH}/MacOS/
    if [ "${1}" = "2020.2" ]; then
        cd ${TRELIS_PATH}/bin
        sudo cp -pv CubitExport.cmake CubitExport.cmake.orig
        sudo gsed -i "s/\"\/\.\.\/app_logger\"/\"\"/" CubitExport.cmake
        sudo gsed -i "s/Trelis-17.1.app/${CUBIT_BASE_NAME}.app/" CubitExport.cmake
        sudo cp -pv CubitUtilConfig.cmake CubitUtilConfig.cmake.orig
        sudo gsed -i "s/\/\.\.\/app_logger\;//" CubitUtilConfig.cmake
        sudo gsed -i "s/Trelis-17.1.app/${CUBIT_BASE_NAME}.app/" CubitGeomConfig.cmake
    fi

}

function build_plugin(){
    cd ${PLUGIN_ABS_PATH}
    cd Trelis-plugin
    rm -rf mcnp2cad
    git clone  https://github.com/bam241/mcnp2cad -b mac_build
    cd ../
    mkdir -pv bld
    cd bld
    cmake ../Trelis-plugin $CMAKE_ADDITIONAL_FLAG \
            -DCubit_DIR=${TRELIS_PATH}/MacOS \
            -DCUBIT_ROOT=${TRELIS_PATH}/MacOS \
            -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
    make -j${PROC}
    make install
}

function build_plugin_pkg(){
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
    chmod 666 svalinn-plugin_linux_cubit_$1.tgz
}
