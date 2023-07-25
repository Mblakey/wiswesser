#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SMITH="${SCRIPT_DIR}/../../data/chembl24.tsv"
PARSE="${SCRIPT_DIR}/../parser/build/readwln"

COUNT=0
TOTAL=$(wc -l < $SMITH)

while read p; do
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f3)
  
  NEW_SMILES=$($PARSE -c -s "${WLN}" 2> /dev/null) # chembl is canonical smiles

  if [ -z $NEW_SMILES ]; then
    echo "$WLN != $SMILES"
    continue
  fi;

  #NEW_SMILES=${NEW_SMILES::-1}
  NEW_SMILES=${NEW_SMILES:0:${#NEW_SMILES}-1}

  if [[ "$SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
  else
  	echo "$WLN != $SMILES"
  fi;

done <$SMITH

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"