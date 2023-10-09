#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SPIDER="${SCRIPT_DIR}/../../data/chemspider.tsv "
PARSE="${SCRIPT_DIR}/../../bin/readwln"
COMP="${SCRIPT_DIR}/../../bin/obcomp"

COUNT=0
MISSED=0
WRONG=0
TOTAL=$(wc -l < $SPIDER)

LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
  NEW_SMILES=$($PARSE -osmi -s "${WLN}" 2> /dev/null) 

  if [ -z $NEW_SMILES ]; then
    ((MISSED++));
    echo "$WLN != anything - $SMILES"
    continue
  fi;

  SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"
  NEW_SMILES="$(sed -e 's/[[:space:]]*$//' <<<${NEW_SMILES})"
  
  SAME=$($COMP "$SMILES" "$NEW_SMILES" 2> /dev/null)

  if [[ "$SAME" == 1 ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
  	echo "$WLN != $SMILES   $NEW_SMILES"
    ((WRONG++));
  fi;

done <$SPIDER

echo -ne "\r$COUNT/$TOTAL correct\n"
echo -ne "$MISSED completely missed\n"
echo -ne "$WRONG wrong output\n"
echo "unit test complete"
