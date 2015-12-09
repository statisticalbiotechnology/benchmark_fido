import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

xlabel = "Protein q value"
ylabel = "Frac proteins with no correct peptides"
#ylabel = "Frac mimic proteins"

def main():
  if sys.argv[1].endswith("proteins"):
    fileName = sys.argv[1]
    fidoQvalues, fdrs = getQvalues(fileName)
    plotQvalues(fidoQvalues, fdrs)
  else:
    fileNameBase = sys.argv[1]
    fidoQvaluesList, fdrsList = [], []
    for i in range(1,11):
      fileName = fileNameBase + "_rep" + str(i) + ".pout.proteins"
      fidoQvalues, fdrs = getQvalues(fileName)
      fidoQvaluesList.append(fidoQvalues)
      fdrsList.append(fdrs)
    plotQvaluesErrorBar(fidoQvaluesList, fdrsList)
  plt.show()
    
def getQvalues(fileName):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  table = [row for row in data]
  x = np.linspace(0, 1, num=1000)
  fp = 0
  p = 0

  fidoQvalues = []
  fdrs = []
  cl = []

  fpSeen = False
  pSeen = False
  nextProteinGroup = -1

  for i in range (1, len(table)):
    fidoQvalue = float(table[i][2])
    
    proteinName = table[i][0] # Protein name
    curProteinGroup = table[i][1] # Protein group ID
    if i < len(table) - 1:
      nextProteinGroup = table[i+1][1]
    else:
      nextProteinGroup = -1
    
    bestPeptideIncorrect = table[i][4][1] == 'i'
    noCorrectPeptide = sum(1 for p in table[i][4].split() if p[1] == 'c') == 0
    proteinAbsent = not 'present' in proteinName
    if noCorrectPeptide or proteinName.startswith('mimic'):
      fpSeen = True
    else: 
      pSeen = True
    
    if curProteinGroup == nextProteinGroup: 
      continue 
    else:
      if pSeen: 
        p = p + 1
      elif fpSeen:
        fp = fp + 1 # Actual false positives
        #print float(fp) / (p + fp), p, fp, table[i][0], table[i][2], table[i][4]
      fpSeen = False
      pSeen = False
    
    fdr = float(fp) / (p + fp)
    # Q-value * Positives = Expected false positives

    fidoQvalues.append(fidoQvalue)
    fdrs.append(fdr)
  return fidoQvalues, fdrs

def plotQvalues(fidoQvalues, fdrs):
  x = np.linspace(0, 1, num=1000)
  plt.plot(fidoQvalues, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

  plt.figure()
  plt.plot(fidoQvalues, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

def plotQvaluesErrorBar(fidoQvaluesList, fdrsList):
  x = list(np.linspace(0, 1, 1001))
  thresholdFdrs = [] 
  ind = len(x)
  for xi, fdr in zip(fidoQvaluesList, fdrsList): 
    newFdr = fdrByThresholding(x, xi, fdr)
    ind = min([ind, newFdr.index(fdr[-1])])
    thresholdFdrs.append(newFdr)

  thresholdFdrs = np.array(thresholdFdrs)
  means = list(thresholdFdrs.mean(axis = 0))
  stdevs = list(thresholdFdrs.std(axis = 0))
  
  # draw error bars
  plt.figure()
  plt.errorbar(x[: ind], means[: ind], yerr = stdevs[: ind], marker = 'o', c = 'r')
  plt.plot(x, x, c = 'k')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

  plt.figure()
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

