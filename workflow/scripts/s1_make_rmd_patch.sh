#!/usr/bin/env bash

UNMODIFIED=$1
MODIFIED=$2
PATCH=$3

# echo $UNMODIFIED
# echo $MODIFIED
# echo $PATCH

diff -u $UNMODIFIED $MODIFIED > $PATCH

exit 0