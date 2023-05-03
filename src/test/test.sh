#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARSE="${SCRIPT_DIR}/../parser/build/wln-writer3"
WLN="$1"
SMILES=$($PARSE "$WLN")
echo -ne "$WLN\t$SMILES\n"
exit 0
