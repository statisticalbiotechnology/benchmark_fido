#!/bin/bash

pi0="0.25"
frac1="0.75"
numspec="20000"
#database="Homo_sapiens.GRCh38.pep.all.with.contaminants.fasta"
database="swissprot_human.fasta"

if [[ $database == swiss* ]]; then
  base="swissprot"
else
  base="ensembl"
fi
inputdir="${base}_pin_rep_pi0_${pi0}_frac1_${frac1}"
outputdir="${base}_fisher_rep_pi0_${pi0}_frac1_${frac1}"
inputbase="${base}_simulated_pi0_${pi0}_frac1_${frac1}"

outputbase=${inputbase}
mkdir -p ${inputdir}
mkdir -p ${outputdir}

cd ../../bin

for i in {1..10}
do
  echo "Generating exp ${i}"
  if [ ! -f ../data/sim_percolator/${inputdir}/${inputbase}_rep${i}.pin.tab ]; then
    python simulate_experiment_psms.py ../data/fasta/${database} --outputPath ../data/sim_percolator/${inputdir}/${inputbase}_rep${i}.pin.tab --present ../data/sim_percolator/${inputdir}/${inputbase}_rep${i}.present_proteins.fasta --fracAbsentProteins ${pi0} --frac1 ${frac1} --numSpectra ${numspec} --seed ${i} >> ../data/sim_percolator/${inputdir}/${outputbase}.generate.log 2>&1
  fi
done

cd ../data/sim_percolator/${outputdir}

for i in {1..10}
do
  echo "Running percolator for exp ${i}"
  outputs="-X ${outputbase}_rep${i}.pout.xml -l ${outputbase}_rep${i}.percolator.tab.proteins -L ${outputbase}_rep${i}.percolator.decoys.tab.proteins -r ${outputbase}_rep${i}.percolator.tab.peptides -B ${outputbase}_rep${i}.percolator.decoys.tab.peptides"
  percolator ../${inputdir}/${inputbase}_rep${i}.pin.tab ${outputs} -f ../../fasta/${database} -c -g -P decoy >> ${outputbase}.percolator.log 2>&1
  
  #outputs="-X ${outputbase}_rep${i}.pout.xml -r ${outputbase}_rep${i}.percolator.tab.peptides -B ${outputbase}_rep${i}.percolator.decoys.tab.peptides"
  #percolator ../${inputdir}/${inputbase}_rep${i}.pin.tab ${outputs} -y -P decoy >> ${outputbase}.percolator.log 2>&1
  
  #percolator -A ../${inputdir}/${inputbase}_rep${i}.pin.tab ${outputs} -q -d 3 -P decoy >> ${outputbase}.percolator.log 2>&1
done
