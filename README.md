# benchmark_fido
Alicia Menéndez Hurtado's scripts for benchmarking Fido

**Definitions.py:**

This script will give as an output two plots. 
The first one is the expected q value (calculated as qv = entrapment/target in each iteration) 
versus the observed FDR (calculated as FDR = ent * Pi_0 / (ent + target) in each iteration). 
The second plot gives back the number of target proteins found vs the FDR. 
The input file for this script is a protein output file from Percolator (-l), using as an input for it 
the ups database (with ~50 proteins) plus a contaminants database with a mimic entrapment of ~450 proteins.

**Different peptides nr.py:**

This script takes as an input three different Percolator protein outputs (-l). The input for Percolator is a simulated data file created with the simulation script, changing the number of peptides in each file (--numSpectra). It also groups proteins with the same peptides together. The output is a plot with the expected q value (taken from the file) and the observed FDR (calculated as above) of all three input files together. 

**Standard deviation pept nr.py:**

In this case we have 10 different simulations for each number of peptides, and the script calculates and plots the average of all of them with the correspondant standard deviation. The input files are the same as before, only a bigger number of them. In this case, they all have the same name but changin a number from 1 to 10, which makes opening easier. As every file has a different length, the script takes the shortest length for the ten simulations. 

**Expectedqv-fdr.py:**

This script will give as an output a plot similar to the first one in *Definitions.py*. It will plot the expected q value against the observed FDR, but it won't group the proteins with the same peptides together. Furthermore, this script also takes as an input the protein output file from Percolator, but the input for Percolator should be a simulated one. Otherwise, the words in c to calculate p and fp should be changed. 

**Exp qv vs fdr.py:**

Another way to do what the above script does.

**Grouping qv vs fdr.py:**

It will do the same as the previous script, only grouping the proteins that have the same peptides, as *Definitions.py* does. 

**FDR-vs-pep.py:**

In this case, the script will plot the observed FDR (calculated as in the other scripts) against the pep value taken from the file. It takes as an input the protein output file from Percolator. Again, if the data is not simulated, the words in c must be changed appropriately.

**pv histogram2.py:**

This script takes as an input the .xml file from Percolator, and plots a histogram of all the p values found there. The range of the i was calculated on the back of the envelope (sorry about that), it should be changed for each file, most likely. 
