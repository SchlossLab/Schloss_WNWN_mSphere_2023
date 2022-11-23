#!/usr/bin/env bash

wget -N https://doi.org/10.1371/journal.pcbi.1003531.s001 -O norarefy-source.zip

unzip -o norarefy-source.zip

mv norarefy-source.zip norarefy-source/

rm -rf __MACOSX