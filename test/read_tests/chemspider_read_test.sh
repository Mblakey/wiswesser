#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SPIDER="${SCRIPT_DIR}/../../data/chemspider.tsv "
PARSE="${SCRIPT_DIR}/../../src/parser/build/readwln"
CANONICAL="${SCRIPT_DIR}/../../src/parser/build/obabel_strip"

COUNT=0
TOTAL=$(wc -l < $SPIDER)
LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  
  NEW_SMILES=$($PARSE -c -ocan -s "${WLN}" 2> /dev/null)
  CAN_SMILES=$($CANONICAL "$SMILES" 2> /dev/null)

  if [ -z $NEW_SMILES ]; then
    echo "$WLN != $CAN_SMILES"
    continue
  fi;

  #NEW_SMILES=${NEW_SMILES::-1}
  NEW_SMILES=${NEW_SMILES:0:${#NEW_SMILES}-1}
  CAN_SMILES=${CAN_SMILES:0:${#CAN_SMILES}-1}

  if [[ "$CAN_SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
  	echo "$WLN != $CAN_SMILES    $NEW_SMILES"
  fi;

done <$SPIDER

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"