#/bin/bash

outputdir="test"
inputbase="simulated_pi0_0.75"
outputbase="simulated_pi0_0.75_fisher"

mkdir -p ${outputdir}

#percolator ${outputbase}.pin.tab -r ${outputdir}/${outputbase}.pout.peptides -y > ${outputdir}/${outputbase}.pout.log 2>&1

#percolator ${inputbase}.pin.tab -A fido -X ${outputdir}/${outputbase}.pout.xml -l ${outputdir}/${outputbase}.pout.proteins -L ${outputdir}/${outputbase}.pout.decoy.proteins -r ${outputdir}/${outputbase}.pout.peptides -q -d 4 -P decoy -I 0.25 > ${outputdir}/${outputbase}.pout.log 2>&1

percolator ${inputbase}.pin.tab -A fisher -f /home/matthew/mergespec/data/db/SP_hum_20110921_clean.fa -X ${outputdir}/${outputbase}.pout.xml -l ${outputdir}/${outputbase}.pout.proteins -L ${outputdir}/${outputbase}.pout.decoy.proteins -r ${outputdir}/${outputbase}.pout.peptides -B ${outputdir}/${outputbase}.pout.decoy_peptides -q -P decoy > ${outputdir}/${outputbase}.pout.log 2>&1

