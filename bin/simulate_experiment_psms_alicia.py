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
    if seq[i] == 'K' or seq[i] == 'R':
      lenP = i - start + 1 
      if lenP >= 10 and lenP <= 30 : yield (seq[start : i + 1])
      start = i + 1
  lenP = lenS - start
  if lenP >= 10 and lenP <= 30 : yield (seq[start : ])


def getDecoyPeptide(peptides):

  for peptide in peptides:
    yield peptide


def insertToGraph(g, allProteins, peptideSet, proteins, prefix, peptideLabel, lookUp, getPeptide = getPeptideByTrypticDigest, decoy = None):

  # get peptide lists for each protein and build the graph
  for protein in proteins:
    proteinName = prefix + protein
    allProteins.add(proteinName)  
    g.add_node(proteinName, bipartite = 0, pep = None)
    
    dPeptides = set()
    for peptide in getPeptide(lookUp[protein]):
      if not peptide in g: 
        g.add_node(peptide, bipartite = 1, pep = None, label = peptideLabel) # t for target, d for decoy
      g.add_edge(proteinName, peptide)
      
      peptideSet.add(peptide)
      if decoy is not None: dPeptides.add(peptide[::-1])
      
    if decoy is not None: decoy[protein] = dPeptides 
 
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
  parser.add_argument('--fracPresentProteins', default = 0.224, type = float, metavar = 'FRAC',
                   help = 'The fraction of all proteins in the file that will be simulated as present')
  parser.add_argument('--numSpectra', type = int, default = 15355,  metavar = 'N',
                   help = 'The number of peptides tha will be matced in the experiment')
  parser.add_argument('--frac0', type = float, default = 0.598,  
                    help = 'The number of peptides that will be simulated as correctly matched with pep=0')
  parser.add_argument('--frac1', type = float, default = 0.195, 
                   help = 'The number of peptides that will be simulated as incorrectly matched with pep=1')
  parser.add_argument('--seed', type = int, default = 1,  metavar = 'S',
                   help = 'The random seed for different reproducable behavior of the simulation')

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

  allProteins = set(allProteinSeqs)
  numProteins = len(allProteins)
  #print numProteins, 
  #print args.fracPresentProteins
  numPresentProteins = int(numProteins * args.fracPresentProteins) 
  numAbsentProteins = numProteins - numPresentProteins
  #numDecoys = numProteins

  # split the proteins into two sets: present/absent, generate the smaller one first
  presentProteins, absentProteins = (), ()
  if numPresentProteins < numAbsentProteins:
    presentProteins = mySample(allProteins, numPresentProteins, localRand)
    absentProteins = allProteins - presentProteins  
  else:
    absentProteins = mySample(allProteins, numAbsentProteins, localRand)
    presentProteins = allProteins - absentProteins
  
  #print "present/absent sampled", 
  #print time.clock() - start_time, "seconds"

  #initialize the protein-peptide bipartite graph and other sets  
  graph, allProteins, decoyProteins = nx.Graph(), set(), {}
  presentPeptides, absentPeptides, decoyPeptides = set(), set(), set()

  # insert the proteins and their respective peptides into the graph
  insertToGraph(graph, allProteins, presentPeptides, presentProteins, "present_", 't', lookUp = allProteinSeqs, decoy = decoyProteins)
  insertToGraph(graph, allProteins, absentPeptides, absentProteins, "absent_", 't', lookUp = allProteinSeqs, decoy = decoyProteins) 
  insertToGraph(graph, allProteins, decoyPeptides, decoyProteins, "decoy_", 'd', lookUp = decoyProteins, getPeptide = getDecoyPeptide, decoy = None)
  
  print "graph is built", 
  print time.clock() - start_time, "seconds"

  # list all target Peptides
  targetPeptides = presentPeptides.union(absentPeptides)  
  # make the three sets exclusive
  absentPeptides = absentPeptides - presentPeptides 
  decoyPeptides = decoyPeptides - targetPeptides

  num1 = int(args.numSpectra * args.frac1) # PEP 1
  num0 = int(args.numSpectra * args.frac0) # PEP 0
  numVague = args.numSpectra - num0 - num1
  vagueStep = 1.0/float(numVague)
  
  print "num1 = %i num0 = %i numvague = %i, total = %i" % (num1, num0, numVague, len(targetPeptides) )

  print >> sys.stderr, "Simulating an experiment with %i spectra, and %i of the %i proteins in the database present" % (args.numSpectra, numPresentProteins, numProteins)
  print >> sys.stderr, "Approximately %i of the (unique) PSMs have a pep=0, %i a pep=1, and %i have a probability somewhere inbetween" % (num0, num1, numVague)
  
  #print >> sys.stderr, "Assigning probabilities to peptides"

  #print >> sys.stderr, "First the peptides with 0<pep<1"

  probabilities = np.arange(0, 1.0, vagueStep) # a ramp function
  us = localRand.uniform(0, 1, numVague) # instead of generating a random number each time, generate them all together

  presentProbs, targetProbs, decoyProbs = [], [], []
  pCount, tCount, dCount = 0, 0, 0
  
  # determine the case and store the respective peps
  for u, peptideProb in zip(us, probabilities):
    totalProb = 2.0 - peptideProb    
    if u < peptideProb/totalProb: 
      pCount += 1
      presentProbs.append(peptideProb)      
    elif u < 1.0/totalProb: 
      tCount += 1
      targetProbs.append(peptideProb)
    else: 
      dCount += 1
      decoyProbs.append(peptideProb)

  #print presentProbs, "presentProb"
  #print targetProbs, "targetProb"
  #print decoyProbs, "decoyProb"

  # Output stats
  sizeSet = len(presentPeptides)
  sizeSample = pCount + num0
  print >> sys.stderr,"Will sample %i peptides out of %i present, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)
  sizeSet = len(targetPeptides)
  sizeSample = pCount+tCount + num1/2 + num0
  print >> sys.stderr,"Will sample %i peptides out of %i target, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)
  sizeSet = len(decoyPeptides)
  sizeSample = dCount + num1/2
  print >> sys.stderr,"Will sample %i peptides out of %i decoys, i.e. %f%%"%(sizeSample, sizeSet ,sizeSample/float(sizeSet)*100)



  # randomly select peptides from each category for the above counts
  presentSample = mySample(presentPeptides, pCount, localRand) 
  targetPeptides = targetPeptides - presentSample
  targetSample = mySample(targetPeptides, tCount, localRand) # sample this from allPeptides instead -LK 
  decoySample = mySample(decoyPeptides, dCount, localRand)

  # exclude the current selection from the three sets
  targetPeptides = targetPeptides - targetSample
  presentPeptides = presentPeptides - presentSample - targetSample
  absentPeptides = absentPeptides - targetSample
  decoyPeptides = decoyPeptides - decoySample
  
  # assign probabilities to peptides in the graph
  for peptide, prob in zip(presentSample, presentProbs): 
    graph.node[peptide]['pep'] = 1.0 - prob
    graph.node[peptide]['prefix'] = 'c'
    #graph.node[peptide]['prefix'] = 'i'
    
  for peptide, prob in zip(targetSample, targetProbs): 
    graph.node[peptide]['pep'] = 1.0 - prob
    graph.node[peptide]['prefix'] = 'i'

  for peptide, prob in zip(decoySample, decoyProbs): 
    graph.node[peptide]['pep'] = 1.0 - prob
    graph.node[peptide]['prefix'] = 'i'
  
  #print >> sys.stderr, "Then the peptides with pep=0"
  
  targetSample = mySample(presentPeptides, num0, localRand)
  for peptide in targetSample: 
    graph.node[peptide]['pep'] = 0.0
    graph.node[peptide]['prefix'] = 'c'

  #remove them to prevent repeated selection
  presentPeptides = presentPeptides - targetSample
  targetPeptides = targetPeptides - targetSample

  print >> sys.stderr, "Last, the peptides with pep=1.0"
  #Adding a small delta to the probability to make these peptides sortable
  sampleSize = num1/2
  #maxD = vagueStep # for handling int, actually should have divided it with 100 

  # generate random deltas
  deltas = localRand.uniform(0, vagueStep, num1)
 
  # sample this from allPeptides instead -LK
  targetSample = mySample(targetPeptides, sampleSize, localRand)
  for peptide, delta in zip(targetSample, deltas[: sampleSize]): 
    graph.node[peptide]['pep'] = 1.0 - delta
    graph.node[peptide]['prefix'] = 'i'
  '''
  pTargetShare = sampleSize * args.fracPresentProteins
  aTargetShare = sampleSize - pTargetShare
  pSample = mySample(presentPeptides, pTargetShare, localRand)
  aSample = mySample(absentPeptides, aTargetShare, localRand)
  for peptide, delta in zip(pSample, deltas[: pTargetShare]): graph.node[peptide]['pep'] = 1.0 - delta
  for peptide, delta in zip(aSample, deltas[pTargetShare: sampleSize]): graph.node[peptide]['pep'] = 1.0 - delta
  '''

  #deltas = localRand.uniform(0, maxD, sampleSize) # sample again for the other half

  decoySample = mySample(decoyPeptides, sampleSize, localRand)
  for peptide, delta in zip(decoySample, deltas[sampleSize:]): 
    graph.node[peptide]['pep'] = 1.0 - delta 
    graph.node[peptide]['prefix'] = 'i'

  print >> sys.stderr, "Sorting all peptides according to their peps"
  
  # get all the peptides as node list and sort by pep
  allTargetDecoyPeptides = list((n, d) for n, d in graph.nodes(data = True) if d['bipartite'] == 1 and d['pep'] is not None)
  outList = sorted(allTargetDecoyPeptides, key = lambda x : x[1]['pep'])

  '''
  debug = open("simDebug.txt", "w")
  prevPep = 0
  prevPeptide = ""
  for item in outList:
    peptide = item[0]
    pep = item[1]['pep']
    if pep > 0:
      if abs(pep - prevPep) < 1e-30:
        print >> debug, peptide,
        print >> debug, pep, " : ",
        print >> debug, prevPep
        print >> debug, prevPeptide
    prevPep = pep 
    prevPeptide = peptide
    
  '''
  print "after sorting"
  print time.clock() - start_time, "seconds"

  #print outList[0 : 10]

  output = open(args.outputPath, "w")
  '''
  print >> output, "PSM\tscore\tqvalue\tpep\tpeptide\tProteins"

  numTarget,numDecoy = 0, 0
  for item in outList:
    peptide = item[0]

    if item[1]['label'] == 't': numTarget += 1
    else: numDecoy += 1
    num = numTarget + numDecoy
  
    linkedProteins = graph.neighbors(peptide)
    proteins = '\t'.join(linkedProteins)
    #print >> output, "%s PSM_%i\t%f\t%f\t%f\t%s\t%s"%(item[1]['label'], num, num, numDecoy/float(numTarget), item[1]['pep'], peptide,proteins)
    print >> output, "PSM_%i\t%f\t%f\t%.16f\t[%s]%s\t%s"%(num, num, numDecoy/float(numTarget), item[1]['pep'], item[1]['prefix'], peptide,proteins)
  '''
  print >> output, "SpecID\tLabel\tScanNr\tfeature1\tpeptide\tProteins"
  fp = 0
  
  numTarget,numDecoy = 0, 0
  for item in outList:
    peptide = item[0]

    if item[1]['label'] == 't': numTarget += 1
    else: numDecoy += 1
    num = numTarget + numDecoy
  
    linkedProteins = graph.neighbors(peptide)
    proteins = '\t'.join(linkedProteins)
    
    
    if 'decoy' in proteins:
      fp = fp + 1
    
    if 'decoy' in proteins:
      lb = -1
    else: 
      lb = 1
        
    peptide = "-." + peptide + ".-"
    
    #print >> output, "%s PSM_%i\t%f\t%f\t%f\t%s\t%s"%(item[1]['label'], num, num, numDecoy/float(numTarget), item[1]['pep'], peptide,proteins)
    #feature1 = numDecoy/float(numTarget)
    #feature2 = item[1]['pep']
    feature1 = num
    print >> output, "PSM_%i\t%i\t%i\t%f\t%s\t%s"%(num, lb, num, feature1, peptide, proteins)
    
  
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
