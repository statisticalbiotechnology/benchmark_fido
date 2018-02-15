import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main(argv):
  #percTargetFile = "/media/storage/mergespec/data/FE_Kall_vial_v2_%d/percolator_tdc/tab/FE_Kall_vial_v2_%d.percolator.tab.peptides"
  #vials = range(1,4)
  
  #percTargetFile = "/media/storage/mergespec/data/FE_Kall_vial_%s/percolator_tdc/tab/FE_Kall_vial_%s.percolator.tab.peptides"
  #percTargetFile = "/media/storage/mergespec/data/FE_Kall_vial_%s/percolator_tide_concat/tab_tdc/FE_Kall_vial_%s.percolator.tab.peptides"
  #vials = ["A","B","AB"]
  
  percTargetFile = argv[0]
  vials = ["AB"]
  
  poolAFastaFile = "/media/storage/mergespec/data/db/prest_pool_a.fasta"
  poolBFastaFile = "/media/storage/mergespec/data/db/prest_pool_b.fasta"
  
  poolAProteins = list(getProteinIds(open(poolAFastaFile,'r')))
  poolBProteins = list(getProteinIds(open(poolBFastaFile,'r')))
  
  for i in vials:
    #percPeptFile = percTargetFile % (i,i)
    percPeptFile = percTargetFile
    print percPeptFile
    if i == "A":
      presentProteins = poolAProteins
    elif i == "B":
      presentProteins = poolBProteins
    elif i == "AB":
      presentProteins = poolAProteins + poolBProteins
    qvals, fdrs = getQvalues(percPeptFile, presentProteins)
    plotQvalues(qvals, fdrs)
  plt.show()

def getProteinIds(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield name
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield name 

def getQvalues(percPeptFile, presentProteins):
  file = open(percPeptFile, 'rb') # The input is the peptide output file (.tab) from Percolator (-l)
  reader = csv.reader(file, delimiter='\t')
  reader.next()
  fp = 1
  tp = 0

  qvals = []
  fdrs = []
  
  dbCorrection = 1.0 - (len(presentProteins) / 1383.0) # for every (1383 - presentProteins) false positives, we expect 1383 false positives
  
  print "\t".join(["PSMId", "PercolatorScore", "reportedQvalue", "observedFDR", "peptide", "protein"])
  for row in reader:
    reportedQvalue = float(row[2])
    
    fpSeen, pSeen = False, False
    #if len(row) > 6:
    #  continue
    for proteinId in row[5:]:
      proteinAbsent = not proteinId in presentProteins
      if proteinAbsent:
        fpSeen = True
      else: 
        pSeen = True
    
    if pSeen: 
      tp = tp + 1
    elif fpSeen:
      fp = fp + 1 # Actual false positives
      if fp <= 20:
        print "\t".join(row[:2] + [str(reportedQvalue), str(round(float(fp) / (tp + fp),4))] + row[4:])
    fdr = float(fp) / dbCorrection / (tp + fp)
    
    qvals.append(reportedQvalue)
    fdrs.append(fdr)
  return qvals, fdrs

def plotQvalues(qvals, fdrs):
  global plotIdx_
  
  xlabel = "Reported q value"
  ylabel = "Observed $FDR_A$"

  plotIdx_ += 1
  
  plt.subplot(3,2, plotIdx_)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(3,2, plotIdx_)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
   
if __name__ == "__main__":
  main(sys.argv[1:])

