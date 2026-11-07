#!/bin/sh

# usage: % /bin/sh ./build_tpl.sh [install directory path]

# To use this, you need to be in a directory, $BASE_DIR, with clones of RAJA
# and umpire 2022.03.1 versions, with both having loaded submodules using
# "git submodule init; git submodule update".  The camp submodule in RAJA is
# used for the camp build. The *.cmake host-config files must also be in the
# current directory

# Steps before using
#
# % git clone git@github.com:LLNL/raja.git
# % cd raja
# % git checkout v2022.03.1
# % git submodule init; git submodule update
# % cd ..
# % git clone git@github.com:LLNL/umpire.git
# % cd umpire
# % git checkout v2022.03.1
# % git submodule init; git submodule update
# % cd ..

# Make sure you have an appropriate node allocation according to the
# policies of the platform you are using.

# This can be ported to other compilers and/or platforms by changing
# the hard-coded host-configs to others that define different configurations.
# Check environment variables

sys_type=${SYS_TYPE:-""}
if [[ -z ${sys_type} ]]
then
    sys_type=${OSTYPE:-""}
    if [[ -z ${sys_type} ]]
    then
        echo "System type not found (both SYS_TYPE and OSTYPE are undefined)"
        exit 1
    fi
fi

if [ $# -ne 2 ]; then
   >&2 echo "usage: /bin/sh ./build_tpl.sh [install directory path] [compiler-name]"
   exit 1
fi

if [ ! -d "$1" ]; then 
   mkdir $1
fi

if [ ! -d "$1" ]; then
   >&2 echo "build_tpl.sh:  Unable to create install directory."
   exit 1
fi

INSTALL_DIR="$(cd "$1" && pwd -P)"
echo $INSTALL_DIR

BASE_DIR=$(pwd)
RAJA_SRC=$BASE_DIR/raja
UMPIRE_SRC=$BASE_DIR/umpire
CAMP_SRC=$BASE_DIR/raja/tpl/camp
RAJA_CONFIG=$BASE_DIR/$2-raja.cmake
UMPIRE_CONFIG=$BASE_DIR/$2-umpire.cmake
CAMP_CONFIG=$BASE_DIR/$2-camp.cmake
cmake_cmd="/usr/tce/packages/cmake/cmake-3.25.2/bin/cmake"
if [[ "${sys_type}" == "toss_4_x86_64_ib_cray" ]]
then
    cmake_cmd="/usr/tce/packages/cmake/cmake-3.24.2/bin/cmake"
fi

if [ ! -f "$RAJA_CONFIG" ]; then
   >&2 echo "build_tpl.sh:  Unable to find file $RAJA_CONFIG"
   exit 1
fi
if [ ! -f "$UMPIRE_CONFIG" ]; then
   >&2 echo "build_tpl.sh:  Unable to find file $UMPIRE_CONFIG"
   exit 1
fi
if [ ! -f "$CAMP_CONFIG" ]; then
   >&2 echo "build_tpl.sh:  Unable to find file $CAMP_CONFIG"
   exit 1
fi

if [ ! -d "$BASE_DIR/raja-build" ]; then
  mkdir $BASE_DIR/raja-build
fi
cd $BASE_DIR/raja-build
${cmake_cmd} $RAJA_SRC -C $RAJA_CONFIG -DBLT_SOURCE_DIR=$BASE_DIR/../blt -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/raja
make -j
make install
cd $BASE_DIR

if [ ! -d "$BASE_DIR/umpire-build" ]; then
  mkdir $BASE_DIR/umpire-build
fi
cd $BASE_DIR/umpire-build
${cmake_cmd} $UMPIRE_SRC -C $UMPIRE_CONFIG -DBLT_SOURCE_DIR=$BASE_DIR/../blt -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/umpire
make -j
make install
cd $BASE_DIR

if [ ! -d "$BASE_DIR/camp-build" ]; then
  mkdir $BASE_DIR/camp-build
fi
cd $BASE_DIR/camp-build
${cmake_cmd} $CAMP_SRC -C $CAMP_CONFIG -DBLT_SOURCE_DIR=$BASE_DIR/../blt -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/camp
make -j
make install
cd $BASE_DIR
