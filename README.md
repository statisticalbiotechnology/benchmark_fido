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

**Example:**

Let's have a look at the whole path to get to *Different peptides nr*'s output. 

First, we get the script [SimulateExperiment](https://github.com/statisticalbiotechnology/inferrensim/blob/master/scripts/simulateExperiment.py) and we run it as:

`$ python SimulateExperiment.py --outputPath 10000peptides.tab --numSpectra 10000 swissprot_human.fasta`

So we have our first file. We want to do this with three different numbers of peptides, so we run it a couple of times more.

`$ python SimulateExperiment.py --outputPath 50000peptides.tab --numSpectra 50000 swissprot_human.fasta`

`$ python SimulateExperiment.py --outputPath 100000peptides.tab --numSpectra 100000 swissprot_human.fasta`

Now we need to run this files in Percolator and get the protein output file. 

`$ percolator -X 10000peps.xml -A -q -d 2 -l 10000peps.prot.tab 10000peptides.tab > 10000peps.psms`

We ask Percolator to make a fine grid search with -d 2, and to give us an extra file with the proteins list with -l name.tab. This is the file we are going to use as an input for our script. We run Percolator another two times to get the other two protein files, and we go to the script.

In lines from 6 to 14, we have:

`data15 = csv.reader(open('Filename1.tab', 'rb'), delimiter='\t')`

`table15 = [row for row in data15]`

So we substitute the Filenames with our own.

`data15 = csv.reader(open('10000peps.prot.tab', 'rb'), delimiter='\t')`

`table15 = [row for row in data15]`

And we write the number of peptides for each line in the plot in lines 129, 130, 131. If Python gives back an error such as "division by 0 in line 122", check in the tab file that the names in the first column are actually "absent" and "present". If not, change them. 
