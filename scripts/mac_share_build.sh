#!/bin/bash


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
            -DCubit_DIR=${CUBIT_PATH}/MacOS \
            -DCUBIT_ROOT=${CUBIT_PATH}/MacOS \
            -DDAGMC_DIR=${PLUGIN_ABS_PATH}/DAGMC \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX=${PLUGIN_ABS_PATH}
    make -j${PROC}
    make install
}


