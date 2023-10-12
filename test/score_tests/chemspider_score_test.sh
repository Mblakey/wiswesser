#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SPIDER="${SCRIPT_DIR}/../../data/chemspider.tsv"
PARSE="${SCRIPT_DIR}/../../bin/readwln"
EXTRACT="${SCRIPT_DIR}/../../bin/wlngrep"
COMP="${SCRIPT_DIR}/../../bin/obcomp"

MISSED_FILE="${SCRIPT_DIR}/missed.tsv"
WRONG_FILE="${SCRIPT_DIR}/wrong.tsv"

COUNT=0
MISSED=0
WRONG=0
TOTAL=0
LINE=0

OPT_EXACT=0
OPT_SCORRECT=0
OPT_GREEDY=0

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: pubchem_score_test.sh [OPTIONS] ARGUMENTS"
        echo "Options:"
        echo "  -h, --help          Show this help message"
        echo "  -x, --fsm-exact     Only match wln if whole line is valid"
        echo "  -c, --correct       use run-time spelling correction"
        echo "  -g, --greedy-parse  compare all valid portions of the wln string"
        exit 0;
        ;;

      -x|--fsm-exact)
        shift # Shift to the next argument
        OPT_EXACT=1
        ;;
      -c|--correct)
        # Handle the -d or --directory option
        shift # Shift to the next argument
        OPT_SCORRECT=1
        ;;
      -g|--greedy-parse)
        shift # Shift to the next argument
        OPT_GREEDY=1
        ;;
      *)
        ;;
    esac
    shift # Shift to the next argument
  done
}

main_a(){

  echo -ne "WLN\tSMILES\n" > "$MISSED_FILE"
  echo -ne "WLN\tSMILES\tATTEMPT\n" > "$WRONG_FILE"

  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)

    if [ "$OPT_EXACT" -eq 1 ]; then
      pWLN=$($EXTRACT -x -s "${WLN}"  2> /dev/null)
      if [ -z "$pWLN" ]; then
        echo "$WLN not valid string"
        continue
      fi;
    fi;

    ((TOTAL++));

    if [ "$OPT_SCORRECT" -eq 1 ]; then
      NEW_SMILES=$($PARSE -osmi -c -s "${WLN}" 2> /dev/null) 
    else
      NEW_SMILES=$($PARSE -osmi -s "${WLN}" 2> /dev/null) 
    fi;


    if [ -z $NEW_SMILES ]; then
      ((MISSED++));
      echo -ne "$WLN\t$SMILES\n" >> "$MISSED_FILE"
      echo "$WLN != anything - $SMILES"
      continue
    fi;

    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"
    NEW_SMILES="$(sed -e 's/[[:space:]]*$//' <<<${NEW_SMILES})"
    
    SAME=$($COMP "$SMILES" "$NEW_SMILES" 2> /dev/null)

    if [[ "$SAME" == 1 ]]; then
      ((COUNT++));
      echo -ne "\r"
    else
      echo "$WLN != $SMILES   $NEW_SMILES"
      echo -ne "$WLN\t$SMILES\t$NEW_SMILES\n" >> "$WRONG_FILE"
      ((WRONG++));
    fi;

  done <$SPIDER

  echo -ne "\r$COUNT/$TOTAL correct\n"
  echo -ne "$MISSED completely missed\n"
  echo -ne "$WRONG wrong output\n"
  echo "unit test complete"
}

main_b(){
  while read p; do
    ((LINE++));
    echo -ne "$LINE: "
    WLN=$(echo -n "$p" | cut -d $'\t' -f1)
    SMILES=$(echo -n "$p" | cut -d $'\t' -f2)
    SMILES="$(sed -e 's/[[:space:]]*$//' <<<${SMILES})"

    WLNarray=$($EXTRACT -o -s "${WLN}"  2> /dev/null)

    while IFS= read -r line; do
      if [ "$OPT_SCORRECT" -eq 1 ]; then
        NEW_SMILES=$($PARSE -osmi -c -s "${line}" 2> /dev/null) 
      else
        NEW_SMILES=$($PARSE -osmi -s "${line}" 2> /dev/null) 
      fi;

      if [ -z $NEW_SMILES ]; then
        continue
      fi;

      NEW_SMILES="$(sed -e 's/[[:space:]]*$//' <<<${NEW_SMILES})"
      SAME=$($COMP "$SMILES" "$NEW_SMILES" 2> /dev/null)

      if [[ "$SAME" == 1 ]]; then
        ((COUNT++));
        break;
      fi;

    done <<< "$WLNarray"
    echo -ne "\r"
  done <$SPIDER

  echo -ne "\r$COUNT/$LINE correctly found\n"
  echo "unit test complete"
}


process_arguments "$@"


if [ "$OPT_GREEDY" -eq 0 ]; then
  main_a
else
  main_b 
fi;

