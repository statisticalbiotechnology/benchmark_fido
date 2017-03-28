import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main():  
  percProtFile = "/media/storage/mergespec/data/FE_Kall_vial_%s/percolator_tdc/tab/FE_Kall_vial_%s.percolator.tab.proteins"
  vials = ["A","B","AB"]
  
  poolAFastaFile = "/media/storage/mergespec/data/db/prest_pool_a.fasta"
  poolBFastaFile = "/media/storage/mergespec/data/db/prest_pool_b.fasta"
  totalDbProteins = 1383
  
  poolAProteins = list(getProteinIds(open(poolAFastaFile,'r')))
  poolBProteins = list(getProteinIds(open(poolBFastaFile,'r')))
  
  for i in vials:
    fileName = percProtFile % (i,i)
    print fileName
    if i == vials[0]:
      presentProteins = poolAProteins
    elif i == vials[1]:
      presentProteins = poolBProteins
    elif i == vials[2]:
      presentProteins = poolAProteins + poolBProteins
    qvals, fdrs, fomrs, tpfps = getQvalues(fileName, presentProteins, totalDbProteins)
    print "#significant proteins:", sum(1 if qval < 0.01 else 0 for qval in qvals)
    plotQvalues(qvals, fdrs, fomrs, tpfps, i)
  plt.subplots_adjust(hspace = 0.4)
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
  fp = 1
  tp = 0

  qvals, fdrs, fomrs, tpfps = list(), list(), list(), list()
    
  numAbsentProteins = totalProteins - len(presentProteins)

  for i in range (1, len(table)):
    reportedQvalue = float(table[i][2])
    
    proteinName = table[i][0]    
    
    fpSeen = False
    pSeen = False
    for proteinId in proteinName.split(","):
      proteinAbsent = not proteinId in presentProteins
      if proteinAbsent:
        fpSeen = True
      else: 
        pSeen = True
    
    if pSeen: 
      tp = tp + 1
    elif fpSeen:
      fp = fp + 1 # Actual false positives
    
    tn = numAbsentProteins - fp
    fn = len(presentProteins) - tp
    fdr = float(fp) / (tp + fp)
    fomr = float(fn) / (tn + fn)

    tpfps.append([tp,fp])
    qvals.append(reportedQvalue)
    fdrs.append(fdr)
    fomrs.append(fomr)
  return qvals, fdrs, fomrs, tpfps

def plotQvalues(qvals, fdrs, fomrs, tpfps, vial):
  global plotIdx_
  numRows = 3
  numCols = 3
  
  xlabel = "Reported q value"
  ylabel = "Observed $FDR_A$"
  labelFontSize = 18
  
  plotIdx_ += 1
  
  plt.subplot(numRows,numCols, plotIdx_)
  
  numIds = sum(1 for fdr in qvals if fdr < 0.01)
  plt.title("Vial " + str(vial) + " (" + str(numIds) + " proteins at 1\% FDR)", fontsize = 24, fontweight = 'bold')
  
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  
  plotIdx_ += 1
  
  plt.subplot(numRows,numCols, plotIdx_)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  
  if False:
    plotIdx_ += 1
    
    plt.subplot(numRows,numCols, plotIdx_)
    plt.plot(fomrs, fdrs)
    plt.axis([0, 0.2, 0, 0.2])
    plt.xlabel("False Ommission Rate", fontsize = labelFontSize)
    plt.ylabel(ylabel, fontsize = labelFontSize)
  
  plotIdx_ += 1
  
  plt.subplot(numRows,numCols, plotIdx_)
  plt.plot([x[1]/float(x[1]+x[0]) for x in tpfps], [x[0]+x[1] for x in tpfps])
  plt.plot([0.01, 0.01], [0, 2000], 'k', linestyle = 'dotted')
  plt.axis([0, 0.1, 0, 400])
  plt.xlabel(ylabel, fontsize = labelFontSize)
  plt.ylabel("Number of protein groups", fontsize = labelFontSize)
   
if __name__ == "__main__":
  main()

