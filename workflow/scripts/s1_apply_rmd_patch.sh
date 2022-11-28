#!/usr/bin/env bash

ORIG_RMD=$1
PATCHED_RMD=$2
PATCH=$3

patch -u $ORIG_RMD -i $PATCH -o $PATCHED_RMD 
