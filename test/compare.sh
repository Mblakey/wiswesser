#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

NEW="${SCRIPT_DIR}/../bin/readwln"
OLD="${SCRIPT_DIR}/../bin/OLD_readwln"
EXTRACT="${SCRIPT_DIR}/../bin/wlngrep"
COMP="${SCRIPT_DIR}/../bin/obcomp"
MODE=""


OLD_COUNT=0
NEW_COUNT=0

OPT_EXACT=0
OPT_GREEDY=0

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: compare.sh <options> <mode>"
        echo "Modes:"
        echo "  chemspider"
        echo "  pubchem"
        echo "  smith"
        echo "  chembl"
        echo "Options:"
        echo "  -h, --help          Show this help message"
        echo "  -x, --exact-match   Use wlngrep for exact matching"
        echo "  -g, --greedy-match   Use wlngrep for greedy matching"
        exit 0;
        ;;
      
      -x|--exact-match)
        shift;
        OPT_EXACT=1;
      ;;
      -g|--greedy-match)
        shift;
        OPT_GREEDY=1;
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

  LINE=0
  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"

    if [ $OPT_EXACT -eq 1 ]; then
      pWLN=$($EXTRACT -x -s "${WLN}"  2> /dev/null)
      if [ -z "$pWLN" ]; then
        echo -ne "\r"
        continue
      fi;
    fi;

    NEW_SMILES=$($NEW -osmi -c -s "${WLN}" 2> /dev/null)
    OLD_SMILES=$($OLD -osmi -s "${WLN}" 2> /dev/null) 
    
    SAME_NEW=""
    SAME_OLD=""

    if [ -n "$NEW_SMILES" ];then 
      SAME_NEW=$($COMP "$SMILES" "$NEW_SMILES" 2> /dev/null)
    fi;

    if [ -n "$OLD_SMILES" ];then 
      SAME_OLD=$($COMP "$SMILES" "$OLD_SMILES" 2> /dev/null)
    fi;

    if [[ "$SAME_NEW" == "1" ]]; then
      ((NEW_COUNT++));
    fi;

    if [[ "$SAME_OLD" == "1" ]]; then
      ((OLD_COUNT++));
      if [ -z $SAME_NEW ] || [[ "$SAME_NEW" == "0" ]]; then
        echo -ne "old got it: $WLN\t$SMILES\n";
      fi;
    fi;

    echo -ne "\r"

  done <$FILE

  echo "new parser: $NEW_COUNT"
  echo "old parser: $OLD_COUNT"
}


process_arguments "$@"
main
exit 0
