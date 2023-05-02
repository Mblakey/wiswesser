#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SMITH="${SCRIPT_DIR}/../../data/smith_wln.tsv"
PARSE="${SCRIPT_DIR}/../parser/build/wln-writer3"

export BABEL_LIBDIR="${SCRIPT_DIR}/../openbabel/build/lib/"

COUNT=0
TOTAL=$(wc -l < $SMITH)

while read p; do
	WLN=$(echo "$p" | cut -d $'\t' -f1)
  SMILES=$(echo "$p" | cut -d $'\t' -f2)
  
  NEW_SMILES=$($PARSE $WLN 2> /dev/null)
  NEW_SMILES=${NEW_SMILES::-1}

  if [[ "$SMILES" == "$NEW_SMILES" ]]; then
  	((COUNT++));
  else
  	echo "$SMILES - $NEW_SMILES"
  fi;

done <$SMITH

echo -ne "\r$COUNT/$TOTAL correct\n"