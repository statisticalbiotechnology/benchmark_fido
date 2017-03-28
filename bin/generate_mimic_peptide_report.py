import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import bisect

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main():
  #percPeptFile = "/media/storage/mergespec/data/103111-Yeast-2hr/percolator_tandem_tdc/tab/103111-Yeast-2hr.percolator.tab.peptides"
  #percPeptFile = "/media/storage/mergespec/data/Percolator_system_tests/percolator_tdc/tab/Percolator_system_tests.percolator.tab.peptides"
  #percPeptFile = "/media/storage/mergespec/data/Percolator_system_tests/percolator_tdc_concat/tab/Percolator_system_tests.percolator.tab.peptides"
  percPeptFile = sys.argv[1]
  fastaFile = "/home/matthew/data/db/swissprot_yeast_160315.with_mimic_s0.04.fasta"
  samplePeptides = getSamplePeptides(open(fastaFile, 'rb'))
  print percPeptFile
  print "NumSamplePeptides = ", len(samplePeptides)
  qvalsReported, qvalsObserved = getQvalues(percPeptFile, samplePeptides)
  plotQvalues(qvalsReported, qvalsObserved)
  plt.tight_layout()
  plt.show()

def getPeptideByTrypticDigest(seq):
  lenS, start = len(seq), 0
  for i in range(lenS):
    if seq[i] == 'K' or seq[i] == 'R':
      lenP = i - start + 1 
      if lenP >= 7 and lenP <= 50 : yield (seq[start : i + 1])
      start = i + 1
  lenP = lenS - start
  if lenP >= 7 and lenP <= 50 : yield (seq[start : ])

# method for reading the given fasta database
def getSamplePeptides(fp):
  name, seq, samplePeptides = None, [], []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name and not "mimic" in name: samplePeptides.extend(list(getPeptideByTrypticDigest("".join(seq))))
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name and not "mimic" in name: samplePeptides.extend(list(getPeptideByTrypticDigest("".join(seq))))
  return set(samplePeptides)

def cleanPeptide(peptide):
  cleanPeptide = ""
  mods = list()
  i = 0
  while i < len(peptide):
    aa = peptide[i]
    cleanPeptide += aa
    if i+1 < len(peptide) and peptide[i+1] == "[":
      i += peptide[i+1:].find("]")+1
    i += 1  
  return cleanPeptide
  
def getQvalues(percPeptFile, samplePeptides):
  file = open(percPeptFile, 'rb') # The input is the peptide output file (.tab) from Percolator (-l)
  reader = csv.reader(file, delimiter='\t')
  reader.next()
  fp = 1
  tp = 0

  fdrs, qvalsReported, numIds = list(), list(), list()
  for row in reader:
    reportedQvalue = float(row[2])
    
    if cleanPeptide(row[4][2:-2]) in samplePeptides or "M" + cleanPeptide(row[4][2:-2]) in samplePeptides:
      tp = tp + 1
    else:
      fp = fp + 1 # Actual false positives
      if fp < 10:
        print cleanPeptide(row[4][2:-2]), row[4], row[5]
    fdr = float(fp) / (tp + fp)
    
    numIds.append(tp+fp)
    qvalsReported.append(reportedQvalue)
    fdrs.append(fdr)
  
  qvalsObserved = fdrsToQvals(fdrs)
  print "Number of identifications at reported q-val < 0.01:", numIds[bisect.bisect_left(qvalsReported, 0.01)]
  print "Number of identifications at observed q-val < 0.01:", numIds[bisect.bisect_left(qvalsObserved, 0.01)]
  
  return qvalsReported, qvalsObserved

def fdrsToQvals(fdrs):
  qvals = [0] * len(fdrs)
  qvals[len(fdrs)-1] = fdrs[-1]
  for i in range(len(fdrs)-2, -1, -1):
    qvals[i] = min(qvals[i+1], fdrs[i])
  return qvals
    
def plotQvalues(qvalsReported, qvalsObserved):
  global plotIdx_
  
  xlabel = "Reported q value"
  ylabel = "Observed FDR"

  plotIdx_ += 1
  
  plt.subplot(1, 2, plotIdx_)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvalsReported, qvalsObserved)
  plt.plot(x, x, 'k--')
  plt.plot(x, x * 0.9, 'k:') # on PSM level we expect 1 extra false positive from the sample database, for each 9 false positives to the entrapment database
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(1, 2, plotIdx_)
  plt.plot(qvalsReported, qvalsObserved)
  plt.plot(x, x, 'k--')
  plt.plot(x, x * 0.9, 'k:') # on PSM level we expect 1 extra false positive from the sample database, for each 9 false positives to the entrapment database
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
   
if __name__ == "__main__":
  main()

