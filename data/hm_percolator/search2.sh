#/bin/bash

outputdir="151009_fido_protein_inference"
base="103111-Yeast-2hr"
db="Saccharomyces_cerevisiae.EF3.64.pep.all.fa"
#base="112111-Human-2h-01"
#db="Homo_sapiens.GRCh38.pep.all.fa"

mkdir -p ${outputdir}

#percolator ${base}.pin.tab -r ${outputdir}/${base}.pout.peptides > ${outputdir}/${base}.pout.log 2>&1

rep_base=${base}-01
percolator -A ${rep_base}.tab -X ${outputdir}/${rep_base}.pout.xml -l ${outputdir}/${rep_base}.pout.proteins -L ${outputdir}/${rep_base}.pout.decoy.proteins -r ${outputdir}/${rep_base}.pout.peptides -q -Z -d 3 > ${outputdir}/${rep_base}.pout.log 2>&1

rep_base=${base}-02
percolator -A ${rep_base}.tab -X ${outputdir}/${rep_base}.pout.xml -l ${outputdir}/${rep_base}.pout.proteins -L ${outputdir}/${rep_base}.pout.decoy.proteins -r ${outputdir}/${rep_base}.pout.peptides -q -Z -d 3 > ${outputdir}/${rep_base}.pout.log 2>&1

rep_base=${base}-03
percolator -A ${rep_base}.tab -X ${outputdir}/${rep_base}.pout.xml -l ${outputdir}/${rep_base}.pout.proteins -L ${outputdir}/${rep_base}.pout.decoy.proteins -r ${outputdir}/${rep_base}.pout.peptides -q -Z -d 3 > ${outputdir}/${rep_base}.pout.log 2>&1

#percolator -A ${base}.tab -X ${outputdir}/${base}.pout.xml -l ${outputdir}/${base}.pout.proteins -L ${outputdir}/${base}.pout.decoy.proteins -r ${outputdir}/${base}.pout.peptides -q -Z -d 3 > ${outputdir}/${base}.pout.log 2>&1

#percolator ${base}.tab -f ../fasta/${db} -c -g -z 7,50 -X ${outputdir}/${base}.pout.xml -l ${outputdir}/${base}.pout.proteins -L ${outputdir}/${base}.pout.decoy.proteins -r ${outputdir}/${base}.pout.peptides -q -Z > ${outputdir}/${base}.pout.log 2>&1
