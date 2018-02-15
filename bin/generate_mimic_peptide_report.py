import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import bisect

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main(argv):
  #percPeptFile = "/media/storage/mergespec/data/103111-Yeast-2hr/percolator_tandem_tdc/tab/103111-Yeast-2hr.percolator.tab.peptides"
  #percPeptFile = "/media/storage/mergespec/data/Percolator_system_tests/percolator_tdc/tab/Percolator_system_tests.percolator.tab.peptides"
  #percPeptFile = "/media/storage/mergespec/data/Percolator_system_tests/percolator_tdc_concat/tab/Percolator_system_tests.percolator.tab.peptides"
  fastaFile = "/home/matthew/data/db/swissprot_yeast_160315.with_mimic_s0.04.fasta"
  samplePeptides = getSamplePeptides(open(fastaFile, 'rb'))
  print "NumSamplePeptides = ", len(samplePeptides)
  
  maxNumIds = 0
  for percPeptFile in argv:
    print percPeptFile
    qvalsReported, qvalsObserved, numIds = getQvaluesAndNumIds(percPeptFile, samplePeptides)
    label = percPeptFile.split('/')[-2].replace('tab_','').replace('_', '\_')
    plotQvalues(qvalsReported, qvalsObserved, label = label)
    plotPerformance(qvalsReported, qvalsObserved, numIds, label = label)
    m = (max(numIds[:bisect.bisect_right(qvalsReported, 0.05)]) / 1000 + 1) * 1000
    maxNumIds = max([m, maxNumIds])
  
  xlabel = "Reported q value"
  ylabel = "Observed FDR"
  
  plt.figure(1)
  plt.subplot(1, 2, 1)
  x = np.linspace(0, 1, num=1000)
  plt.plot(x, x, 'k--')
  plt.plot(x, x / 1.116, 'k:') # on PSM level we expect 1 extra false positive from the sample database, for each 9 false positives to the entrapment database
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  plt.legend(loc = 'upper left')
  
  plt.subplot(1, 2, 2)
  plt.plot(x, x, 'k--')
  plt.plot(x, x / 1.116, 'k:') # on PSM level we expect 1 extra false positive from the sample database, for each 9 false positives to the entrapment database
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  plt.legend(loc = 'upper left')
  
  plt.tight_layout()
  
  xlabel = "Reported q value"
  ylabel = "Num significant PSMs"

  plt.figure(2)
  plt.subplot(1, 2, 1)
  x = np.linspace(0, 1, num=1000)
  plt.plot([0.01, 0.01], [0, maxNumIds], 'k--')
  plt.axis([0, 0.05, 0, maxNumIds])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  plt.legend(loc = 'lower right')
  
  xlabel = "Observed q value"
  
  plt.subplot(1, 2, 2)
  plt.plot([0.01, 0.01], [0, maxNumIds], 'k--')
  plt.axis([0, 0.05, 0, maxNumIds])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  plt.legend(loc = 'lower right')
  
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
  
def getQvaluesAndNumIds(percPeptFile, samplePeptides):
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
  
  return qvalsReported, qvalsObserved, numIds

def fdrsToQvals(fdrs):
  qvals = [0] * len(fdrs)
  qvals[len(fdrs)-1] = fdrs[-1]
  for i in range(len(fdrs)-2, -1, -1):
    qvals[i] = min(qvals[i+1], fdrs[i])
  return qvals
    
def plotQvalues(qvalsReported, qvalsObserved, label = ""):  
  plt.figure(1)
  plt.subplot(1, 2, 1)
  plt.plot(qvalsReported, qvalsObserved, label = label)
  
  plt.subplot(1, 2, 2)
  plt.plot(qvalsReported, qvalsObserved, label = label)

def plotPerformance(qvalsReported, qvalsObserved, numIds, label = ""):
  plt.figure(2)
  plt.subplot(1, 2, 1)
  plt.plot(qvalsReported, numIds, label = label)
  
  plt.subplot(1, 2, 2)
  plt.plot(qvalsObserved, numIds, label = label)
  
if __name__ == "__main__":
  main(sys.argv[1:])

