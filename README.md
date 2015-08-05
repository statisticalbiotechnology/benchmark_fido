# benchmark_fido
Alicia Menéndez Hurtado's scripts for benchmarking Fido

**Definitions.py:**

This script will give as an output two plots. 
The first one is the expected q value (calculated as qv = entrapment/target in each iteration) 
versus the observed FDR (calculated as FDR = ent * Pi_0 / (ent + target) in each iteration). 
The second plot gives back the number of target proteins found vs the FDR. 
The input file for this script is a protein output file from Percolator (-l), using as an input for it 
the ups database (with ~50 proteins) plus a contaminants database with a mimic entrapment of ~450 proteins.

`$ python Definitions.py File.tab`

**Different peptides nr.py:**

This script takes as an input three different Percolator protein outputs (-l) and the names for the tags on the plot.

`$ python Different\peptides\nr.py File1.tab Tag1 File2.tab Tag2 File3.tab Tag3`

The input for Percolator is a simulated data file created with the simulation script, changing the number of peptides in each file (--numSpectra). It also groups proteins with the same peptides together. The output is a plot with the expected q value (taken from the file) and the observed FDR (calculated as above) of all three input files together. 

**Standard deviation pept nr.py:**

In this case we have X different simulations for each number of peptides, and the script calculates and plots the average of all of them with the correspondant standard deviation. The input files are the same as before, only a bigger number of them. In this case, they all have the same name but changin a number from 1 to X, which makes opening easier. As every file has a different length, the script takes the shortest length for the ten simulations. Here's an example of how to run it with input files First1.tab, First2.tab,..., FirstX.tab, Second1.tab, Second2.tab, etc:

`$ python Standard\deviation\pept\nr.py First%d.tab Tag1 Second%d.tab Tag2 Third%d.tab Tag3 X`

**Expectedqv-fdr.py:**

This script will give as an output a plot similar to the first one in *Definitions.py*. It will plot the expected q value against the observed FDR, but it won't group the proteins with the same peptides together. Furthermore, this script also takes as an input the protein output file from Percolator, but the input for Percolator should be a simulated one. Otherwise, the words in c to calculate p and fp should be changed. 

**Exp qv vs fdr.py:**

Another way to do what the above script does.

**Grouping qv vs fdr.py:**

It will do the same as the previous script, only grouping the proteins that have the same peptides, as *Definitions.py* does. 

**FDR-vs-pep.py:**

In this case, the script will plot the observed FDR (calculated as in the other scripts) against the pep value taken from the file. It takes as an input the protein output file from Percolator. Again, if the data is not simulated, the words in c must be changed appropriately.

**pv histogram2.py:**

This script takes as an input the .xml file from Percolator, and plots a histogram (with as many bricks as in  Nr_of_bins) of all the p values found there. The range of the i (Nr_of_rows) was calculated on the back of the envelope (sorry about that), it should be changed for each file, most likely. It has three input arguments:

`Filename.xml Nr_of_rows Nr_of_bins`

###**Example:**

Let's have a look at the whole path to get to *Different peptides nr*'s output. 

First, we get the script [SimulateExperiment](https://github.com/statisticalbiotechnology/inferrensim/blob/master/scripts/simulateExperiment.py) and we run it as:

`$ python SimulateExperiment.py --outputPath 15000peptides.tab --numSpectra 15000 swissprot_human.fasta`

So we have our first file. We want to do this with three different numbers of peptides, so we run it a couple of times more.

`$ python SimulateExperiment.py --outputPath 30000peptides.tab --numSpectra 30000 swissprot_human.fasta`

`$ python SimulateExperiment.py --outputPath 60000peptides.tab --numSpectra 60000 swissprot_human.fasta`

Now we need to run this files in Percolator and get the protein output file. 

`$ percolator -X 15000peps.xml -A -q -d 2 -l 15000peps.prot.tab 15000peptides.tab > 15000peps.psms`

We ask Percolator to make a fine grid search with -d 2, and to give us an extra file with the proteins list with -l name.tab. This is the file we are going to use as an input for our script. We run Percolator another two times to get the other two protein files, and we can run the script.

`$ python Different\peptides\nr.py 15000peps.prot.tab 15.000\peptides 30000peps.prot.tab 30.000\peptides 60000peps.prot.tab 60.000\peptides`


