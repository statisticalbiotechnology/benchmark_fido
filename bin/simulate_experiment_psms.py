#!/usr/bin/python
#
# Copyright 2014 Lukas Kall, KTH
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import operator as op
import argparse
import numpy as np
#from numpy import random as nprnd
import networkx as nx
import itertools as it
import time
from scipy.stats import gumbel_l

# method for reading the given fasta database
def readFasta(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, "".join(seq))
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield (name, "".join(seq))

# method for generating the list of peptides
def getPeptideByTrypticDigest(seq):
  lenS, start = len(seq), 0
  for i in range(lenS):
    if (seq[i] == 'K' or seq[i] == 'R') and seq[min(lenS-1,i+1)] != 'P':
      lenP = i - start + 1 
      if lenP >= 6 and lenP <= 50 : yield (seq[start : i + 1])
      start = i + 1
  lenP = lenS - start
  if lenP >= 6 and lenP <= 50 : yield (seq[start : ])

def getDecoyPeptide(peptides):
  for peptide in peptides:
    yield peptide

def insertToGraph(g, allProteins, peptideSet, proteins, prefix, peptideLabel, lookUp, getPeptide = getPeptideByTrypticDigest, decoyStrategy = ""):

  # get peptide lists for each protein and build the graph
  for protein in proteins:
    proteinName = prefix + protein
    allProteins.add(proteinName)  
    g.add_node(proteinName, bipartite = 0, pep = None)
    
    dPeptides = set()
    proteinSeq = lookUp[protein]
    if decoyStrategy == "reverse-protein":
      proteinSeq = proteinSeq[::-1]
    for peptide in getPeptide(proteinSeq):
      if decoyStrategy == "reverse-peptide": 
        peptide = peptide[::-1]
      
      peptideSet.add(peptide)
      if not peptide in g: 
        g.add_node(peptide, bipartite = 1, pep = None, label = peptideLabel) # t for target, d for decoy
      g.add_edge(proteinName, peptide)
 
# method for generating random sample of size n, without repeatation
def mySample(population, size, localRand):
  
  assert size <= len(population)
  population = list(population)  
  localRand.shuffle(population)
  return set(population[:size])
  '''
  lenP = len(population)
  #print [lenP, size]
  randIndices = localRand.choice(lenP, size, replace = False)
  #print randIndices
  sample = set(population[i] for i in randIndices[0: size])
  return sample
  '''


#################################################### main method 
def main(**kwargs): 
  
  start_time = time.clock()  
  argList = kwargs.get('args', None)
  #print argList
  
  docu = '''A script for simulating the outcome of a mass spec experiment.  
  The script takes a sequence database as input and randomly assignes a 
  fraction of them as present and the rest of the proteins as
  absent. The proteins constituent peptides are subsequently simulated
  as being matched against a set of (virtual) spectra. In the process
  some of the peptides are simulated as being correctly and some
  incorrectly matched. peptides from absent proteins are simulated as
  allways incorrectly matched, while correctly matched peptides are all
  simulated as steming from present poteins.'''

  parser = argparse.ArgumentParser(description = docu)
  parser.add_argument('fastaDatabase',
                   help = 'The sequence database which we have matched with our simulated search engine')
  parser.add_argument('--fracAbsentProteins', default = 0.75, type = float, metavar = 'FRAC',
                   help = 'The fraction of all proteins in the file that will be simulated as absent')
  parser.add_argument('--numSpectra', type = int, default = 20000, metavar = 'N',
                   help = 'The number of peptides that will be matced in the experiment')
  parser.add_argument('--frac1', type = float, default = 0.5,  
                    help = 'The number of peptides that will be simulated as correct')
  parser.add_argument('--seed', type = int, default = 1,  metavar = 'S',
                   help = 'The random seed for different reproducable behavior of the simulation')
  parser.add_argument('--present', default = None,  metavar = 'PATH',
                   help = 'Output the present protein database as a fasta file named by PATH')
  parser.add_argument('--decoy', default = None,  metavar = 'PATH',
                   help = 'Output the manifactured decoy database as a fasta file named by PATH')
  parser.add_argument('--outputPath', default = 'simulatorOutput.txt',  metavar = 'OUTPATH',
                   help = 'File path where the simulator output will be saved')
  
  if argList: args = parser.parse_args(argList) 
  #  print args
  else: args = parser.parse_args() # parse sys.argv, if not otherwise specified

  localRand = np.random #make a local instance of random to make it thread-safe
  localRand.seed(args.seed) # Setting random seed to assure reproducible behavior

  ####################### skip this step if a dict is received as an optional argument 

  allProteinSeqs = kwargs.get('allProteinSeqs', None)
  if allProteinSeqs is None:
    #print "reading fasta database"
    allProteinSeqs = {}
    with open(args.fastaDatabase) as f:
      for proteinName, proteinSeq in readFasta(f):
        allProteinSeqs[proteinName] = proteinSeq  
    #print "file read",
    #print time.clock() - start_time, "seconds"

  allTargetProteins = set(allProteinSeqs)
  numProteins = len(allTargetProteins)
  #print numProteins, 
  #print args.fracAbsentProteins
  numAbsentProteins = int(numProteins * args.fracAbsentProteins) 
  numPresentProteins = numProteins - numAbsentProteins
  #numDecoys = numProteins

  # split the proteins into two sets: present/absent, generate the smaller one first
  presentProteins, absentProteins = (), ()
  if numPresentProteins < numAbsentProteins:
    presentProteins = mySample(allTargetProteins, numPresentProteins, localRand)
    absentProteins = allTargetProteins - presentProteins  
  else:
    absentProteins = mySample(allTargetProteins, numAbsentProteins, localRand)
    presentProteins = allTargetProteins - absentProteins
  
  #print "present/absent sampled", 
  #print time.clock() - start_time, "seconds"

  #initialize the protein-peptide bipartite graph and other sets  
  graph, allProteins, decoyProteins = nx.Graph(), set(), {}
  presentPeptides, absentPeptides, decoyPeptides = set(), set(), set()

  # insert the proteins and their respective peptides into the graph
  insertToGraph(graph, allProteins, presentPeptides, presentProteins, "", 't', lookUp = allProteinSeqs)
  insertToGraph(graph, allProteins, absentPeptides, absentProteins, "", 't', lookUp = allProteinSeqs) 
  insertToGraph(graph, allProteins, decoyPeptides, allTargetProteins, "decoy_", 'd', lookUp = allProteinSeqs, decoyStrategy = "reverse-protein")
  
  print "graph is built", 
  print time.clock() - start_time, "seconds"
  
  # print the present protein database
  if args.present:
    with open(args.present,"w") as ofp:
      # reverse the protein sequences
      for protein in presentProteins:
        print >> ofp, ">%s" % protein
  
  # list all target Peptides
  targetPeptides = presentPeptides.union(absentPeptides)  
  # make the three sets exclusive
  absentPeptides = absentPeptides - presentPeptides 
  decoyPeptides = decoyPeptides - targetPeptides

  num1 = int(args.numSpectra * args.frac1) # correct PSMs
  num0 = args.numSpectra - num1 # incorrect PSMs
  
  print "num1 = %i num0 = %i, total = %i" % (num1, num0, len(targetPeptides) )

  print >> sys.stderr, "Simulating an experiment with %i spectra, and %i of the %i proteins in the database present" % (args.numSpectra, numPresentProteins, numProteins)
  print >> sys.stderr, "Approximately %i of the (unique) PSMs are correct, %i are incorrect" % (num1, num0)

  # Output stats
  sizeSet = len(presentPeptides)
  sizeSample = num1
  print >> sys.stderr,"Will sample %i peptides out of %i present, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)
  sizeSet = len(targetPeptides)
  sizeSample = num0
  print >> sys.stderr,"Will sample %i peptides out of %i target, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)
  sizeSet = len(decoyPeptides)
  sizeSample = num0 + num1
  print >> sys.stderr,"Will sample %i peptides out of %i decoys, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)

  # randomly select peptides from each category for the above counts
  presentSample = mySample(presentPeptides, num1, localRand) 
  targetPeptides = targetPeptides - presentSample
  targetSample = mySample(targetPeptides, num0, localRand) # sample this from allPeptides instead -LK 
  decoySample = mySample(decoyPeptides, num0 + num1, localRand)
  
  presentProbs = gumbel_l.rvs(0.4, 0.7, size = num1)
  # assign probabilities to peptides in the graph
  for peptide, prob in zip(presentSample, presentProbs): 
    graph.node[peptide]['pep'] = prob
    graph.node[peptide]['prefix'] = 'c'
    #graph.node[peptide]['prefix'] = 'i'
  
  targetProbs = gumbel_l.rvs(-0.75, 0.3, size = num0)
  for peptide, prob in zip(targetSample, targetProbs): 
    graph.node[peptide]['pep'] = prob
    graph.node[peptide]['prefix'] = 'i'
  
  decoyProbs = gumbel_l.rvs(-0.75, 0.3, size = num0 + num1)
  for peptide, prob in zip(decoySample, decoyProbs): 
    graph.node[peptide]['pep'] = prob
    graph.node[peptide]['prefix'] = 'i'
    
  # get all the peptides as node list and shuffle
  outList = list((n, d) for n, d in graph.nodes(data = True) if d['bipartite'] == 1 and d['pep'] is not None)
  np.random.shuffle(outList)

  print "after shuffling"
  print time.clock() - start_time, "seconds"

  #print outList[0 : 10]

  output = open(args.outputPath, "w")

  print >> output, "SpecID\tLabel\tScanNr\tfeature1\tpeptide\tProteins"
  fp = 0
  
  numTarget,numDecoy = 0, 0
  for i, item in enumerate(outList):
    peptide = item[0]

    linkedProteins = graph.neighbors(peptide)
    proteins = '\t'.join(linkedProteins)
    
    label = -1
    if item[1]['label'] == 't': 
      label = 1
      numTarget += 1
      scannr = numTarget
    else:
      label = -1
      numDecoy += 1
      scannr = numDecoy
    
    peptide = "-.[" + item[1]['prefix'] + "]" + peptide + ".-"
    print >> output, "PSM_%i\t%i\t%i\t%f\t%s\t%s"%(scannr, label, scannr, item[1]['pep'], peptide, proteins)
  
  # print the decoy database
  if args.decoy:
    with open(args.decoy,"w") as ofp:
      # reverse the protein sequences
      for protein in allProteins:
        if protein.startswith("decoy"): continue
        proteinId = protein.split("_")[1]
        print >> ofp, ">%s" % protein
        decoySeq = "".join(decoyProteins[proteinId])
        print >> ofp, "%s" % decoySeq

  
  print "done",
  print time.clock() - start_time, "seconds"

  #print dLen, "decoys"
  

#############################################################################################################
if __name__ == "__main__":
    main()
