#!/bin/bash

cd ../../bin

base="simulated_pi0_75.frac1_75"
outputdir="perc_validation"

#python simulate_experiment_psms.py ../data/fasta/swissprot_human.fasta --outputPath ../data/sim_percolator/${base}.pin.tab --fracAbsentProteins 0.75 --frac1 0.75

cd ../data/sim_percolator
mkdir -p ${outputdir}

#percolator ${base}.pin.tab -A fisher -f /home/matthew/mergespec/data/db/SP_hum_20110921_clean.fa -X ${outputdir}/${base}.pout.xml -l ${outputdir}/${base}.pout.proteins -L ${outputdir}/${base}.pout.decoy.proteins -r ${outputdir}/${base}.pout.peptides -B ${outputdir}/${base}.pout.decoy_peptides -q -y -P decoy > ${outputdir}/${base}.pout.log 2>&1

#percolator ${base}.pin.tab -X ${outputdir}/${base}.pout.xml -r ${outputdir}/${base}.pout.peptides -B ${outputdir}/${base}.pout.decoy_peptides -q -P decoy > ${outputdir}/${base}.pout.log 2>&1

percolator ${base}.pin.tab -X ${outputdir}/${base}.no_tdc.pout.xml -r ${outputdir}/${base}.no_tdc.pout.peptides -B ${outputdir}/${base}.no_tdc.pout.decoy_peptides -q -y -P decoy > ${outputdir}/${base}.no_tdc.pout.log 2>&1
