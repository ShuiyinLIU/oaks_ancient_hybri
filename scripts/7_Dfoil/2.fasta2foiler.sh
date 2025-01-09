#!/bin/bash
mkdir -p counts/
names=$1
fasta=$2
outgroup=$3
cpu=$4

cat $names | parallel -j $cpu --env outgroup --env fasta 'fasta2dfoil.py <(selectSeqs.pl -f <(echo {} | tr " " $"\n"; cat '$outgroup') '$fasta') --out counts/$(echo {} | tr " " ".").counts --names $(echo {} | tr " " ","),$(cat '$outgroup')'


