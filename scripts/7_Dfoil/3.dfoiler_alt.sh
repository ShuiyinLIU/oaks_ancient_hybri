#!/bin/bash
mkdir -p dfoil
mkdir -p precheck
cpu=$1

find ./counts/ -type f -name '*.counts' | parallel -j $cpu dfoil.py --mode dfoilalt --infile {} --out dfoil/{/.}.dfoil_alt ">" precheck/{/.}.precheck_alt
