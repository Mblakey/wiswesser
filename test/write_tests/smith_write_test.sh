#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SMITH="${SCRIPT_DIR}/../../data/smith.tsv"
WRITER="${SCRIPT_DIR}/../../src/parser/build/writewln"
READER="${SCRIPT_DIR}/../../src/parser/build/readwln"
CANONICAL="${SCRIPT_DIR}/../../src/parser/build/obabel_strip"

COUNT=0
TOTAL=$(wc -l < $SMITH)

LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
  WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  CAN_SMILES=$($CANONICAL "$SMILES" 2> /dev/null)

  NEW_WLN=$($WRITER -ismi -s "${CAN_SMILES}" 2> /dev/null) # chembl is canonical smiles

  if [ -z "$NEW_WLN" ]; then
    echo "$SMILES != anything"
    continue
  fi;

  NEW_SMILES=$($READER -ocan -s "${NEW_WLN}" 2> /dev/null) 
  
  if [ -z "$NEW_SMILES" ]; then
    echo "$NEW_WLN != anything $CAN_SMILES"
    continue
  fi;

  if [[ "$CAN_SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
    echo "$WLN != $NEW_WLN   $CAN_SMILES"
  fi;

done <$SMITH

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"