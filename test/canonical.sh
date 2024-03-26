#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PARSE="${SCRIPT_DIR}/../build/readwln"
COMP="${SCRIPT_DIR}/../build/obcomp"
MODE=""

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: canonical.sh <options> <mode>"
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

  LINE=0
  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)

    SMILES=$($PARSE -ocan "${WLN}" 2> /dev/null)
    CAN_WLN=$($PARSE -owln "${WLN}" 2> /dev/null) 
    RE_SMILES=$($PARSE -ocan "${CAN_WLN}" 2> /dev/null) 
          
    SAME_NEW=$($COMP "$SMILES" "$RE_SMILES" 2> /dev/null)
    if [ -z $SAME_NEW ] || [[ "$SAME_NEW" == "0" ]]; then
      echo -ne "$WLN failed canonicalisation - $CAN_WLN\n";
    elif [ ${#WLN} -gt ${#CAN_WLN} ]; then
      echo -ne "$CAN_WLN passes but is longer than $WLN\n"
    fi;

    echo -ne "\r"
  done <$FILE
}


process_arguments "$@"
main
