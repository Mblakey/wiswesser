#!/bin/bash

set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if ! command -v git &> /dev/null;then
    echo "git could not be found"
    exit 1
fi

if ! command -v cmake &> /dev/null;then
    echo "cmake could not be found"
    exit 1
fi

if ! command -v make &> /dev/null;then
    echo "make could not be found"
    exit 1
fi

if [ ! -d "$SCRIPT_DIR/external/" ]
then
  mkdir "$SCRIPT_DIR/external/"
fi

if [ -d "$SCRIPT_DIR/external/openbabel" ] 
then
    mkdir -p "$SCRIPT_DIR/external/openbabel/build"
    cd "$SCRIPT_DIR/external/openbabel/build"
    cmake .. -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/external"
    make -j 10
    make install
else
    cd external
    git clone https://github.com/openbabel/openbabel.git
    mkdir "$SCRIPT_DIR/external/openbabel/build"
    cd "$SCRIPT_DIR/external/openbabel/build"
    cmake .. -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/external"
    make -j 10
    make install 
fi

mkdir -p "$SCRIPT_DIR/bin"
cd "$SCRIPT_DIR/bin"
cmake ..
make 

echo ""
echo "Build successful"
exit 0
