#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CHEMBL="${SCRIPT_DIR}/../../data/chembl24.tsv"
WRITER="${SCRIPT_DIR}/../../src/parser/build/writewln"
READER="${SCRIPT_DIR}/../../src/parser/build/readwln"
CANONICAL="${SCRIPT_DIR}/../../src/parser/build/obabel_strip"

COUNT=0
TOTAL=$(wc -l < $CHEMBL)

LINE=0
while read p; do
  ((LINE++));

  WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f3)
  CAN_SMILES=$($CANONICAL "$SMILES" 2> /dev/null)

  NEW_WLN=$($WRITER -ismi -s "${CAN_SMILES}" 2> /dev/null) # chembl is canonical smiles

  if [ -z "$NEW_WLN" ]; then
    echo "$LINE: $SMILES != anything"
    continue
  fi;

  NEW_SMILES=$($READER -ocan -s "${NEW_WLN}" 2> /dev/null) 
  
  if [ -z "$NEW_SMILES" ]; then
    echo "$LINE: $NEW_WLN != anything"
    continue
  fi;

  if [[ "$CAN_SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
  else
    echo "$LINE: $WLN != $NEW_WLN   $CAN_SMILES"
  fi;

done <$CHEMBL

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"