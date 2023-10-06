#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

ENG="${SCRIPT_DIR}/../../data/english_words.txt"
PARSE="${SCRIPT_DIR}/../../bin/readwln"

echo "Checking whole english language for WLN notation ..."
echo "Performing large seg fault unit test, this will take multiple minutes ..."


LINE=0
while read ENTRY; do
  ((LINE++));
  echo -ne "$LINE: "
	SMILES=$($PARSE -osmi -s '${ENTRY^^}' 2> /dev/null)
	if [ -n "$SMILES" ]; then
		echo -ne "${ENTRY^^} = $SMILES\n"
	else
    echo -ne "\r"
  fi
done <$ENG

echo "unit test complete - english_wln.tsv is now in the data directory"
exit 0
