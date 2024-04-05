#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

INCORRECT="${SCRIPT_DIR}/incorrect.txt"
CORRECT="${SCRIPT_DIR}/correct.txt"
WRITER="${SCRIPT_DIR}/../build/writewln"
READER="${SCRIPT_DIR}/../build/readwln"
COMP="${SCRIPT_DIR}/../build/obcomp"
MODE=""
FILE=""


ARG_COUNT=0
process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -w|--wrong)
        WRONG=1
        echo "Failed SMILES" > $INCORRECT
        ;;
       -c|--correct)
        RIGHT=1
        echo "Correct SMILES" > $CORRECT
        ;;
      -h|--help)
        # Display a help message
        echo "Usage: writing.sh <options> <mode> <file>"
        echo "modes:"
        echo "  chemspider"
        echo "  pubchem"
        echo "  smith"
        echo "  chembl"
        echo "  external"
        echo "options:"
        echo "  -h, --help          show this help message"
        exit 0;
        ;;
      *)
        if [ $ARG_COUNT -eq 0 ]; then
          MODE=$arg;
        elif [ $ARG_COUNT -eq 1 ]; then
          FILE=$arg
        fi
        ((ARG_COUNT++))
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
  elif [ "$MODE" == "external" ]; then
    echo "performing unit test on $FILE"
  else
    echo "mode choice invalid!, for help use -h flag";
    exit 1;
  fi;

  COUNT=0
  TOTAL=$(wc -l < $FILE)

  LINE=0
  while read p; do
    ((LINE++));
    if [ -t 1 ]; then
      echo -ne "$LINE: "
    fi;
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"

    NEW_WLN=$($WRITER -ismi "${SMILES}" 2> /dev/null) # chembl is canonical smiles

    if [ -z "$NEW_WLN" ]; then
      echo -ne "$WLN\tno write\t${SMILES}\n"
      continue
    fi;

    NEW_WLN="$(sed -e 's/[[:space:]]*$//' <<<${NEW_WLN})"
    NEW_SMILES=$($READER -osmi "${NEW_WLN}" 2> /dev/null) 
    
    if [ -z "$NEW_SMILES" ]; then
      echo -ne "$NEW_WLN\tre-read fail\t$SMILES\n"

      continue
    fi;

    SAME=$($COMP "$SMILES" "$NEW_SMILES")
    if [[ "$SAME" == "1" ]]; then
      ((COUNT++));
      if [ -t 1 ]; then
        echo -ne "\r"
      fi; 
    else
      echo -ne  "$NEW_WLN\tnot equal\t$SMILES"
    fi;

  done <$FILE

  if [ -t 1 ]; then
    echo -ne "\r$COUNT/$TOTAL correct\n"
  fi; 
}

process_arguments "$@"
main
exit 0
