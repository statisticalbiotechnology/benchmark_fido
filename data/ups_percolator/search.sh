#/bin/bash

outputdir="151009_fido_protein_inference"
base="20080315_CPTAC6_22_6QC1"
#base="merged"
db="ups.with_contaminants.non_redundant.fasta"

mkdir -p ${outputdir}

#percolator ${base}.pin.tab -r ${outputdir}/${base}.pout.peptides > ${outputdir}/${base}.pout.log 2>&1

percolator -A ${base}.pin.tab -X ${outputdir}/${base}.pout.xml -l ${outputdir}/${base}.pout.proteins -L ${outputdir}/${base}.pout.decoy.proteins -r ${outputdir}/${base}.pout.peptides -q -Z -d 4 > ${outputdir}/${base}.pout.log 2>&1

#percolator ${base}.pin.tab -f ../fasta/${db} -c -g -z 6,40 -X ${outputdir}/${base}.pout.xml -l ${outputdir}/${base}.pout.proteins -L ${outputdir}/${base}.pout.decoy.proteins -r ${outputdir}/${base}.pout.peptides -q -Z > ${outputdir}/${base}.pout.log 2>&1
