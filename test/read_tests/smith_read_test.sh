#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SMITH="${SCRIPT_DIR}/../../data/smith.tsv"
PARSE="${SCRIPT_DIR}/../../src/parser/build/readwln"

COUNT=0
MISSED=0
WRONG=0
TOTAL=$(wc -l < $SMITH)
LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  
  NEW_SMILES=$($PARSE -osmi -s "${WLN}" 2> /dev/null)

  if [ -z $NEW_SMILES ]; then
    ((MISSED++));
    echo "$WLN != $SMILES"
    continue
  fi;

  #NEW_SMILES=${NEW_SMILES::-1}
  NEW_SMILES=${NEW_SMILES:0:${#NEW_SMILES}-1}

  if [[ "$SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
  	echo "$WLN != $SMILES   $NEW_SMILES"
    ((WRONG++));
  fi;

done <$SMITH

echo -ne "\r$COUNT/$TOTAL correct\n"
echo -ne "$MISSED completely missed\n"
echo -ne "$WRONG wrong output\n"
echo "unit test complete"