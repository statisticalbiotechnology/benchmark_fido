import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm
import bisect
import collections
import scipy.stats as stats
import subprocess

def main():
  dataSet = "pandey"
  percPin = ""
  
  percTabBases = list()
  if dataSet == "prest":
    percTabFolder = "/media/storage/mergespec/data/FE_Kall_vial_%d/percolator_tdc_with_randoms_full_digest/tab_no_pept_cutoff/"
    percTabBase = os.path.join(percTabFolder, "FE_Kall_vial_%d.percolator")
    for i in range(1,4):
      percTabBases.append(percTabBase % (i,i))
  elif dataSet == "sim":
    pi0 = "0.25"
    frac1 = "0.25"
    percTabBase = "/home/matthew/benchmark_fido/data/sim_percolator/swissprot_fisher_rep_pi0_%s_frac1_%s/swissprot_simulated_pi0_%s_frac1_%s_rep%s.percolator" % (pi0, frac1, pi0, frac1,"%d")
    for i in range(1,11):
      percTabBases.append(percTabBase % i)
      break
  elif dataSet == "pandey":
    percPin = "/media/storage/mergespec/data/Pandey/percolator_tdc/pin/Pandey.tab"
    percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_1M/Pandey.percolator")
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator")
      
  for percTabBase in percTabBases:
    writeProteinFeatures(percTabBase, percPin)

def writeProteinFeatures(percTabBase, percPin):
  targetFN = percTabBase + ".tab.peptides"
  decoyFN = percTabBase + ".decoys.tab.peptides"
  targetPsmFN = percTabBase + ".tab.psms"
  decoyPsmFN = percTabBase + ".decoys.tab.psms"
  outputFN = percTabBase + ".pin.peptides"
  
  print targetFN
  targetPsmIds = getPsmIds(targetFN)
  decoyPsmIds = getPsmIds(decoyFN)
  
  #peptideCounts = collections.defaultdict(int)
  #peptideCounts = getPeptideCountsPSM(targetPsmFN, peptideCounts)
  #peptideCounts = getPeptideCountsPSM(decoyPsmFN, peptideCounts)
  
  #peptideCounts = getPeptideCounts(percPin)
  
  reader = csv.reader(open(percPin, 'r'), delimiter = '\t')
  header = reader.next()
  defaultDirection = reader.next()
  #header.insert(24,"NumPSMs")
  #defaultDirection.insert(24,0)
  with open(outputFN,'w') as f:
    writer = csv.writer(f, delimiter = '\t')
    writer.writerow(header)
    writer.writerow(defaultDirection)
    print outputFN
    
    for row in reader:
      if (int(row[1]) == 1 and row[0] in targetPsmIds) or (int(row[1]) == -1 and row[0] in decoyPsmIds):
        #row.insert(24, peptideCounts[row[24]])
        writer.writerow(row)
  
  print "Running percolator"
  cmd = "percolator -r %s.peptides -B %s.decoys.peptides %s -N 100000 > %s.log 2>&1" % (outputFN, outputFN, outputFN, outputFN)
  print cmd
  subprocess.call(cmd, shell = True)

def getPeptideCounts(percPin):
  reader = csv.reader(open(percPin, 'r'), delimiter = '\t')
  header = reader.next()
  defaultDirection = reader.next()
  
  peptideCounts = collections.defaultdict(int)
  for row in reader:
    peptideCounts[row[24]] += 1
  return peptideCounts

def getPeptideCountsPSM(percPsms, peptideCounts):
  reader = csv.reader(open(percPsms, 'r'), delimiter = '\t')
  header = reader.next()
  
  for row in reader:
    #if float(row[2]) < 0.01:
      peptideCounts[row[4]] += 1
    #else:
    #  break
  return peptideCounts
  
def getPsmIds(peptFile):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  psmIds = list()
  for row in reader:
    psmIds.append(row[0])
  return set(psmIds)

if __name__ == "__main__":
  main()

