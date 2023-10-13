#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

FILE=$1
PARSE="${SCRIPT_DIR}/../../src/parser/build/readwln"

if [ -z ${FILE} ]; then
  echo "Error: no input file detected! should be first arg!"
  exit 1
fi;

LINE=0
while read ENTRY; do
  ((LINE++));
  echo -ne "$LINE: "
	SMILES=$($PARSE -osmi -s '${ENTRY^^}' 2> /dev/null)
	if [ -n "$SMILES" ]; then
		echo -ne "${ENTRY^^} = $SMILES\n"
	else
    echo -ne "\r"
  fi
done <$FILE

echo "unit test complete"
exit 0