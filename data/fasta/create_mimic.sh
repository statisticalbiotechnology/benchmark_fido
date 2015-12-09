#!/bin/bash

cp ups.with_contaminants.non_redundant.fasta ups.with_contaminants.non_redundant.with_mimic.fasta

for i in {1..4}; do 
  mimic mimic${i}_ ups.with_contaminants.non_redundant.with_mimic.fasta >> ups.with_contaminants.non_redundant.with_mimic.fasta
  sleep 2
done

cd ../../bin

python reverse_fasta.py
