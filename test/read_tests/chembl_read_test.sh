#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CHEMBL="${SCRIPT_DIR}/../../data/chembl24.tsv"
PARSE="${SCRIPT_DIR}/../../bin/readwln"
COMP="${SCRIPT_DIR}/../../bin/obcomp"

COUNT=0
MISSED=0
WRONG=0
TOTAL=$(wc -l < $CHEMBL)

LINE=0
while read p; do
  ((LINE++));
  echo -ne "$LINE: "
	WLN=$(echo -n "$p" | cut -d $'\t' -f1)
  SMILES=$(echo -n "$p" | cut -d $'\t' -f3)
  NEW_SMILES=$($PARSE -ocan -s "${WLN}" 2> /dev/null) 

  if [ -z $NEW_SMILES ]; then
    ((MISSED++));
    echo "$WLN != anything - $SMILES"
    continue
  fi;

  SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"
  NEW_SMILES="$(sed -e 's/[[:space:]]*$//' <<<${NEW_SMILES})"

  SAME=$($COMP "$SMILES" "$NEW_SMILES")

  if [[ "$SAME" == "1" ]]; then
  	((COUNT++));
    echo -ne "\r"
  else
  	echo -ne  "$WLN != $SMILES\t$NEW_SMILES\n"
    ((WRONG++));
  fi;

done <$CHEMBL

echo -ne "\r$COUNT/$TOTAL correct\n"
echo -ne "$MISSED completely missed\n"
echo -ne "$WRONG wrong output\n"
echo "unit test complete"
