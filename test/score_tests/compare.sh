#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

NEW="${SCRIPT_DIR}/../../bin/readwln"
OLD="${SCRIPT_DIR}/../../bin/OLD_readwln"
EXTRACT="${SCRIPT_DIR}/../../bin/wlngrep"
COMP="${SCRIPT_DIR}/../../bin/obcomp"
FILE=""


OLD_COUNT=0
NEW_COUNT=0

OPT_LEGACY=0

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: pubchem_score_test.sh [OPTIONS] <FILE.tsv>"
        echo "Options:"
        echo "  -l, --legacy        Show this help message"
        echo "  -h, --help          Show this help message"
        ;;
      -l|--legacy)
        shift
        OPT_LEGACY=1
        ;;
      *)
        FILE=$arg
        ;;
    esac
    shift # Shift to the next argument
  done
}

main(){

  if [ -z $FILE ]; then 
    echo "No file inputted!";
    exit 0;
  fi; 

  LINE=0
  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"

    pWLN=$($EXTRACT -x -s "${WLN}"  2> /dev/null)
    if [ -z "$pWLN" ]; then
      echo -ne "\r"
      continue
    fi;

    if [ $OPT_LEGACY -eq 1 ]; then
      NEW_SMILES=$($NEW -osmi -c -l -s "${WLN}" 2> /dev/null)  
    else
      NEW_SMILES=$($NEW -osmi -c -s "${WLN}" 2> /dev/null)
    fi; 
    
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
