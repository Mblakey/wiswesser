#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
READ="${SCRIPT_DIR}/../bin/readwln"
WRITE="${SCRIPT_DIR}/../bin/writewln"
FILE=""
MODE=0

process_arguments() {
  # Loop through all the arguments
  for arg in "$@"; do
    case "$arg" in
      -h|--help)
        # Display a help message
        echo "Usage: file.sh <mode> <file>"
        echo "mode"
        echo "  -r, --read     read wln to smiles"
        echo "  -w, --write    write smiles to wln"
        exit 0;
        ;;
      -r|--read)
        MODE=1
        ;;
      -w|--write)
        MODE=2
        ;;
      *)
        FILE=$arg
        ;;
    esac
    shift # Shift to the next argument
  done

  if [ -z $FILE ]; then
    echo "Error: no input file detected! should be first arg!"
    exit 1
  fi;

  if [ $MODE -eq 0 ]; then
    echo "Error: mode is not set!"
    exit 1; 
  fi;
}

main(){
  LINE=0
  
  while read ENTRY; do
    ((LINE++));
    OUT=""
    if [ $MODE -eq 1 ]; then
      OUT=$($READ -osmi -s "${ENTRY}" 2> /dev/null)
    elif [ $MODE -eq 2 ]; then 
      OUT=$($WRITE -ismi -s "${ENTRY}" 2> /dev/null)
    fi; 

    if [ -n "$OUT" ]; then
      echo -ne "${ENTRY}\t$OUT\n"
    else
      echo -ne "${ENTRY}\tnull"
    fi
  done <$FILE
  exit 0
}

process_arguments "$@"
main