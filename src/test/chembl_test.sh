#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CHEMBL="${SCRIPT_DIR}/../../data/chembl24.tsv"
PARSE="${SCRIPT_DIR}/../parser/build/readwln"
CANONICAL="${SCRIPT_DIR}/../parser/build/obabel_strip"

rm chembl_false.tsv 2> /dev/null

COUNT=0
TOTAL=$(wc -l < $CHEMBL)

while read p; do
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f3)
  

  CAN_SMILES=$($CANONICAL "$SMILES" 2> /dev/null)

  NEW_SMILES=$($PARSE -c -s "${WLN}" 2> /dev/null) # chembl is canonical smiles

  if [ -z $NEW_SMILES ]; then
    echo "$WLN != $SMILES"
    continue
  fi;

  NEW_SMILES=${NEW_SMILES:0:${#NEW_SMILES}-1}
  CAN_SMILES=${CAN_SMILES:0:${#CAN_SMILES}-1}

  if [[ "$CAN_SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
  else
  	echo "$WLN != $CAN_SMILES"
    echo "$WLN  $CAN_SMILES $NEW_SMILES" >> chembl_false.tsv
  fi;

done <$CHEMBL

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"