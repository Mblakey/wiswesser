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


if [ -d "$SCRIPT_DIR/src/openbabel" ] 
then
    mkdir -p "$SCRIPT_DIR/src/openbabel/build"
    cd "$SCRIPT_DIR/src/openbabel/build"
    cmake ..
    make -j 10
else
    cd src
    git clone https://github.com/openbabel/openbabel.git
    mkdir "$SCRIPT_DIR/src/openbabel/build"
    cd "$SCRIPT_DIR/src/openbabel/build"
    cmake ..
    make -j 10
fi

mkdir -p "$SCRIPT_DIR/src/parser/build"
cd "$SCRIPT_DIR/src/parser/build"
cmake ..
make 

echo ""
echo "Build successful - please run:" 
echo "  source ./link.sh"
echo "before using parser"
exit 0