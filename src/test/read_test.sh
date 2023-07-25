#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

FILE=$1
PARSE="${SCRIPT_DIR}/../parser/build/readwln"

RED='\033[0;41m'
NC='\033[0m' # No Color


if [[ -z "$FILE" ]]; then
  echo "Error: no file path provided"
  exit 1;
fi; 

if [ ! -f $FILE ]; then
  echo "Error: no file found at path"
  exit 1;
fi 

COUNT=0
TOTAL=$(wc -l < $FILE)

while read p; do
	WLN=$(echo -n "$p")
  NEW_SMILES=$($PARSE -s "${WLN}" 2> /dev/null)

  echo -ne "${WLN} :"

  if [ -z $NEW_SMILES ]; then
    echo -e "${RED}FAIL${NC}"
  else
    echo "PASS"
    ((COUNT++));
  fi;

done <$FILE

echo -ne "\r$COUNT/$TOTAL passed read\n"
echo "unit test complete"