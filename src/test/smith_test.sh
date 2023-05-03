#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SMITH="${SCRIPT_DIR}/../../data/smith_wln.tsv"
PARSE="${SCRIPT_DIR}/../parser/build/wln-writer3"

COUNT=0
TOTAL=$(wc -l < $SMITH)

while read p; do
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  
  NEW_SMILES=$($PARSE "${WLN}" 2> /dev/null)

  if [ -z $NEW_SMILES ]; then
    echo "$WLN != $SMILES"
    continue
  fi;

  NEW_SMILES=${NEW_SMILES::-1}

  if [[ "$SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
  else
  	echo "$WLN != $SMILES"
  fi;

done <$SMITH

echo -ne "\r$COUNT/$TOTAL correct\n"