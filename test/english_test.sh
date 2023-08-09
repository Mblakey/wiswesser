#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

ENG="${SCRIPT_DIR}/../data/english_words.txt"
PARSE="${SCRIPT_DIR}/../src/parser/build/readwln"

echo -n "" > "${SCRIPT_DIR}/../data/english_wln.tsv"

echo "Checking whole english language for WLN notation ..."
echo "Performing large seg fault unit test, this will take multiple minutes ..."

while read ENTRY; do
	SMILES=$($PARSE -osmi -s "${ENTRY^^}" 2> /dev/null)
	if [ -n "$SMILES" ]; then
		echo -ne "${ENTRY^^}\t$SMILES\n"
		echo -ne "${ENTRY^^}\t$SMILES\n" >> "${SCRIPT_DIR}/../../data/english_wln.tsv"
	fi

done <$ENG

echo "unit test complete - english_wln.tsv is now in the data directory"
exit 0