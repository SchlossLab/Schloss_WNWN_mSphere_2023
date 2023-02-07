#!/usr/bin/env bash

# Concatenates N-1 files and outputs them to the Nth filename. Uses the header
# row from the first file

N_ARGS=${#@}
INPUT_FILES=${@:1:($N_ARGS-1)}
FIRST_FILE=$1
OUTPUT_FILE=${@:$N_ARGS:1}

head -n1 $FIRST_FILE > $OUTPUT_FILE
tail -n +2 -q $INPUT_FILES | sort -n -k 3 -t _ >> $OUTPUT_FILE
