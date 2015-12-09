import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main():
  #percProtFile = "/media/storage/mergespec/data/FE_Kall_vial_%d/percolator_tdc/tab/FE_Kall_vial_%d.percolator.tab.proteins"
  #percProtFile = "/media/storage/mergespec/data/FE_Kall_vial_%d/percolator_tdc/tab/FE_Kall_vial_%d.percolator.tab.picked_proteins"
  percProtFile = "/media/storage/mergespec/data/FE_Kall_vial_%d/percolator_tdc_with_randoms_full_digest/tab_no_pept_cutoff/FE_Kall_vial_%d.percolator.tab.picked_proteins"
  poolAFastaFile = "/home/matthew/mergespec/data/db/prest_pool_a.fasta"
  poolBFastaFile = "/home/matthew/mergespec/data/db/prest_pool_b.fasta"
  
  poolAProteins = list(getProteinIds(open(poolAFastaFile,'r')))
  poolBProteins = list(getProteinIds(open(poolBFastaFile,'r')))
  
  for i in range(1,4):
    fileName = percProtFile % (i,i)
    print fileName
    if i == 1:
      presentProteins = poolAProteins + poolBProteins
    elif i == 2:
      presentProteins = poolBProteins
    elif i == 3:
      presentProteins = poolAProteins
    qvals, fdrs, fomrs = getQvalues(fileName, presentProteins, 1383)
    print "#significant proteins:", sum(1 if qval < 0.01 else 0 for qval in qvals)
    plotQvalues(qvals, fdrs, fomrs)
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

def getQvalues(fileName, presentProteins, totalProteins):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  table = [row for row in data]
  x = np.linspace(0, 1, num=1000)
  fp = 0
  tp = 0

  qvals = []
  fdrs = []
  fomrs = []

  fpSeen = False
  pSeen = False
  nextProteinGroup = -1
  
  numAbsentProteins = totalProteins - len(presentProteins)

  for i in range (1, len(table)):
    reportedQvalue = float(table[i][2])
    
    proteinName = table[i][0] # Protein name
    curProteinGroup = table[i][1] # Protein group ID
    if i < len(table) - 1:
      nextProteinGroup = table[i+1][1]
    else:
      nextProteinGroup = -1
    
    for proteinId in proteinName.split(","):
      proteinAbsent = not proteinId in presentProteins
      if proteinAbsent:
        fpSeen = True
      else: 
        pSeen = True
    
    if curProteinGroup == nextProteinGroup: 
      continue 
    else:
      if pSeen: 
        tp = tp + 1
      elif fpSeen:
        fp = fp + 1 # Actual false positives
        if fp < 10:
          print proteinId, float(fp) / (tp + fp), reportedQvalue
        #print float(fp) / (p + fp), p, fp, table[i][0], table[i][2], table[i][4]
      fpSeen = False
      pSeen = False
    
    tn = numAbsentProteins - fp
    fn = len(presentProteins) - tp
    fdr = float(fp) / (tp + fp)
    fomr = float(fn) / (tn + fn)
    # Q-value * Positives = Expected false positives

    qvals.append(reportedQvalue)
    fdrs.append(fdr)
    fomrs.append(fomr)
  return qvals, fdrs, fomrs

def plotQvalues(qvals, fdrs, fomrs):
  global plotIdx_
  
  xlabel = "Reported q value"
  ylabel = "Observed $FDR_A$"

  plotIdx_ += 1
  
  plt.subplot(3,3, plotIdx_)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(3,3, plotIdx_)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(3,3, plotIdx_)
  plt.plot(fomrs, fdrs)
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel("FOR", fontsize = 18)
  plt.ylabel("FDR", fontsize = 18)
   
if __name__ == "__main__":
  main()

