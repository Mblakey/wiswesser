#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


PARSE="${SCRIPT_DIR}/../build/readwln"
COMP="${SCRIPT_DIR}/../build/obcomp"

COUNT=0
MISSED=0
WRONG=0

MODE=""
FILE=""

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: reading.sh <options> <mode>"
        echo "Modes:"
        echo "  chemspider"
        echo "  pubchem"
        echo "  smith"
        echo "  chembl"
        echo "Options:"
        echo "  -h, --help          Show this help message"
        exit 0;
        ;;
      *)
        MODE=$arg
        ;;
    esac
    shift # Shift to the next argument
  done
}

main(){

  if [ -z $MODE ]; then 
    echo "No mode selected!, for help use -h flag";
    exit 1;
  elif [ "$MODE" == "chemspider" ]; then
    FILE="${SCRIPT_DIR}/../data/unit_test/chemspider.tsv"
  elif [ "$MODE" == "pubchem" ]; then
    FILE="${SCRIPT_DIR}/../data/unit_test/pubchem.tsv"
  elif [ "$MODE" == "smith" ]; then
    FILE="${SCRIPT_DIR}/../data/unit_test/smith.tsv"
  elif [ "$MODE" == "chembl" ]; then
    FILE="${SCRIPT_DIR}/../data/unit_test/chembl24.tsv"
  else
    echo "mode choice invalid!, for help use -h flag";
    exit 1;
  fi;

  TOTAL=$(wc -l < $FILE)

  LINE=0
  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"


    NEW_SMILES=$($PARSE -osmi "${WLN}" 2> /dev/null) 

    if [ -z $NEW_SMILES ]; then
      ((MISSED++));
      echo "$WLN != anything - $SMILES"
      continue
    fi;

    NEW_SMILES="$(sed -e 's/[[:space:]]*$//' <<<${NEW_SMILES})"

    SAME=$($COMP "$SMILES" "$NEW_SMILES")

    if [[ "$SAME" == "1" ]]; then
      ((COUNT++));
      echo -ne "\r"
    else
      echo -ne  "$WLN != $SMILES\t$NEW_SMILES\n"
      ((WRONG++));
    fi;

  done <$FILE

  echo -ne "\r$COUNT/$TOTAL correct\n"
  echo -ne "$MISSED completely missed\n"
  echo -ne "$WRONG wrong output\n"
}

process_arguments "$@"
main
exit 0
