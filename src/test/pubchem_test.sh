#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PUBCHEM="${SCRIPT_DIR}/../../data/pubchem.txt"
PARSE="${SCRIPT_DIR}/../parser/build/readwln"

COUNT=0
TOTAL=$(wc -l < $PUBCHEM)

LINE=0
while read p; do
  ((LINE++));

	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  #SMILES=$(echo -n "$p" | cut -d $'\t' -f3)
  
  NEW_SMILES=$($PARSE -c -s "${WLN}" 2> /dev/null) # chembl is canonical smiles

  if [ -z $NEW_SMILES ]; then
    echo "$LINE: $WLN != anything"
  else
    ((COUNT++));
  fi;

done <$PUBCHEM

echo -ne "\r$COUNT/$TOTAL correct\n"
echo "unit test complete"