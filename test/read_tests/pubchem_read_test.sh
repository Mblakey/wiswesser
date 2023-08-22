#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PUB="${SCRIPT_DIR}/../../data/pubchem.tsv"
PARSE="${SCRIPT_DIR}/../../src/parser/build/readwln"
CANONICAL="${SCRIPT_DIR}/../../src/parser/build/obabel_strip"

COUNT=0
TOTAL=$(wc -l < $PUB)

LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  
  CAN_SMILES=$($CANONICAL "$SMILES" 2> /dev/null)
  NEW_SMILES=$($PARSE -ocan -s "${WLN}" 2> /dev/null) # chembl is canonical smiles

  if [ -z $NEW_SMILES ]; then
    echo "$WLN != anything"
    continue
  fi;

  NEW_SMILES=${NEW_SMILES:0:${#NEW_SMILES}-1}
  CAN_SMILES=${CAN_SMILES:0:${#CAN_SMILES}-1}

  if [[ "$CAN_SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
    echo "$WLN != $CAN_SMILES    $NEW_SMILES"
  fi;

done <$PUB

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"