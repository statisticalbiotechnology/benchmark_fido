import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main():
  pi0 = sys.argv[1]
  frac1 = sys.argv[2]
  percProtFile = "/home/matthew/benchmark_fido/data/sim_percolator/swissprot_fisher_rep_pi0_%s_frac1_%s/swissprot_simulated_pi0_%s_frac1_%s_rep%s.percolator.tab.classic_proteins" % (pi0, frac1, pi0, frac1,"%d")
  presentProteinFile = "/home/matthew/benchmark_fido/data/sim_percolator/swissprot_pin_rep_pi0_%s_frac1_%s/swissprot_simulated_pi0_%s_frac1_%s_rep%s.present_proteins.fasta" % (pi0, frac1, pi0, frac1,"%d")
  
  reportedQvalList, fdrAsList, fdrIsList = [], [], []
  for i in range(1,11):
    print "Protein file:", percProtFile % i
    presentProteins = list(getProteinIds(open(presentProteinFile % i,'r')))
    fileName = percProtFile % i
    qvals, fdrAs, fdrIs, fomrs = getQvalues(fileName, presentProteins, 20201)
    reportedQvalList.append(qvals)
    fdrAsList.append(fdrAs)
    fdrIsList.append(fdrIs)
    
  plotQvaluesErrorBar(reportedQvalList, fdrAsList, fdrIsList)
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
  fp, tp, fpI, tpI = 0, 0, 0, 0

  qvals = []
  fdrAs, fdrIs = [], []
  fomrs = []

  fpSeen, pSeen, fpISeen, pISeen = False, False, False, False
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
      bestPeptideIncorrect = table[i][4][1] == 'i'
      noCorrectPeptide = sum(1 for p in table[i][4].split() if p[1] == 'c') == 0
      proteinAbsent = not proteinId in presentProteins
      if proteinAbsent:
        fpSeen = True
      else: 
        pSeen = True
      if bestPeptideIncorrect:
        fpISeen = True
      else: 
        pISeen = True
      
    if curProteinGroup == nextProteinGroup: 
      continue 
    else:
      if pSeen: 
        tp = tp + 1
      elif fpSeen:
        fp = fp + 1 # Actual false positives
        #if fp < 10:
        #  print proteinId, float(fp) / (tp + fp), reportedQvalue
        #print float(fp) / (p + fp), p, fp, table[i][0], table[i][2], table[i][4]
      if pISeen: 
        tpI = tpI + 1
      elif fpISeen:
        fpI = fpI + 1
        if fpI < 10:
          print proteinId, float(fp) / (tp + fp), reportedQvalue
      fpSeen, pSeen, fpISeen, pISeen = False, False, False, False
    
    tn = numAbsentProteins - fp
    fn = len(presentProteins) - tp
    fdrA = float(fp) / (tp + fp)
    fdrI = float(fpI) / (tpI + fpI)
    fomr = float(fn) / (tn + fn)
    # Q-value * Positives = Expected false positives

    qvals.append(reportedQvalue)
    fdrAs.append(fdrA)
    fdrIs.append(fdrI)
    fomrs.append(fomr)
  print "False positives FDR_A", fp
  print "False positives FDR_I", fpI
  return qvals, fdrAs, fdrIs, fomrs

def plotQvalues(qvals, fdrs, fomrs):
  global plotIdx_
  
  xlabel = "Reported q value"
  ylabel = "Observed $FDR_A$"

  plotIdx_ += 1
  
  plt.subplot(3,1, plotIdx_)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(3,1, plotIdx_)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(3,1, plotIdx_)
  plt.plot(fomrs, fdrs)
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel("FOR", fontsize = 18)
  plt.ylabel("FDR", fontsize = 18)

def plotQvaluesErrorBar(reportedQvalList, fdrAsList, fdrIsList):
  plotQvaluesErrorBarLocal(reportedQvalList, fdrAsList, "Reported q value", "Observed $FDR_A$", 1)
  plotQvaluesErrorBarLocal(reportedQvalList, fdrIsList, "Reported q value", "Observed $FDR_I$", 3)

def plotQvaluesErrorBarLocal(reportedQvalList, fdrsList, xlabel, ylabel, plotIdx):
  x = list(np.linspace(0, 1, 1001))
  thresholdFdrs = [] 
  ind = len(x)
  for xi, fdr in zip(reportedQvalList, fdrsList): 
    newFdr = fdrByThresholding(x, xi, fdr)
    ind = min([ind, newFdr.index(fdr[-1])])
    thresholdFdrs.append(newFdr)

  thresholdFdrs = np.array(thresholdFdrs)
  means = list(thresholdFdrs.mean(axis = 0))
  stdevs = list(thresholdFdrs.std(axis = 0))
  
  # draw error bars
  plt.subplot(2,2,plotIdx)
  plt.errorbar(x[: ind], means[: ind], yerr = stdevs[: ind], marker = 'o', c = 'r')
  plt.plot(x, x, c = 'k')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

  plt.subplot(2,2,plotIdx+1)
  plt.errorbar(x[: ind], means[: ind], yerr = stdevs[: ind], marker = 'o', c = 'r')
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)


def fdrByThresholding(newFoAs, foAs, fdrs):

  qFoAs = qValues(foAs)
  qFdrs = qValues(fdrs)
  
  #print newFoAs, "cutoffs"
  #plt.figure()
  #plt.plot(qFoAs, qFdrs)
  #plt.plot(newFoAs)
  #plt.show()   

  newFdrs = []
  for threshold in newFoAs:
    closestFoAPosition = bisect_right(qFoAs, threshold)
    if closestFoAPosition == 0: newFdr = 0
    else: newFdr = qFdrs[closestFoAPosition - 1]   
    newFdrs.append(newFdr)

  return newFdrs

#methods for qValue, fdr thresholding etc
def qValues(inputL):
  qValues = []
  # calculate q-values for each item in the inputL
  nextVal = 1
  for curVal in inputL[:: -1]:
    qVal = min(curVal, nextVal)
    qValues.append(qVal)
    nextVal = qVal
  qValues.reverse()
  return qValues
   
if __name__ == "__main__":
  main()

