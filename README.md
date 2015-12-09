# benchmark_fido
Alicia Menéndez Hurtado's scripts for benchmarking Fido

Definitions.py:

This script will give as an output two plots. 
The first one is the expected q value (calculated as qv = entrapment/target in each iteration) 
versus the observed FDR (calculated as FDR = ent * Pi_0 / (ent + target) in each iteration). 
The second plot gives back the number of target proteins found vs the FDR. 

The input file for this script is a protein output file from Percolator (-l), using as an input for it 
the ups database (with ~50 proteins) plus a contaminants database with a mimic entrapment of ~450 proteins.

