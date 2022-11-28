#!/usr/bin/env bash

PM=$1
PDS=$2

COMPARE_FILE=$(echo $PM | sed -E 's/html/compare.txt/')
echo "$COMPARE_FILE"

sdiff $PM $PDS > "$COMPARE_FILE"

echo "done."
