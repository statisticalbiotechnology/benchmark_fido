import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

xlabel = "Picked q value"
ylabel = "Frac proteins with incorrect best peptide"

def main():
  if sys.argv[1].endswith("proteins"):
    fileName = sys.argv[1]
    decoyFileName = fileName.replace(".proteins", ".decoy.proteins")
    pickedQvalues, fdrs = getQvalues(fileName, decoyFileName)
    plotQvalues(pickedQvalues, fdrs)
  else:
    fileNameBase = sys.argv[1]
    pickedQvaluesList, fdrsList = [], []
    for i in range(1,11):
      fileName = fileNameBase + "_rep" + str(i) + ".pout.proteins"
      decoyFileName = fileNameBase + "_rep" + str(i) + ".pout.decoy.proteins"
      pickedQvalues, fdrs = getQvalues(fileName, decoyFileName)
      pickedQvaluesList.append(pickedQvalues)
      fdrsList.append(fdrs)
    plotQvaluesErrorBar(pickedQvaluesList, fdrsList)
  plt.show()
    
def getQvalues(fileName, decoyFileName):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  
  decoyFile = open(decoyFileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  decoyData = csv.reader(decoyFile, delimiter='\t')
  
  table = [row for row in data]
  decoyTable = [row for row in decoyData]
  x = np.linspace(0, 1, num=1000)
  fp = 0
  p = 0

  pickedQvalues = []
  fdrs = []
  cl = []

  fpSeen = False
  pSeen = False
  nextProteinGroup = -1
  j = 1
  d = 0
  t = 0
  
  seenProteins = list()
  for i in range (1, len(table)):
    while j < len(decoyTable) and float(decoyTable[j][2]) <= float(table[i][2]):
      decoyProteins = map(lambda x : x[6:], decoyTable[j][0].split(","))
      seen = False
      for decoyProtein in decoyProteins:
        if decoyProtein in seenProteins:
          seen = True
          break
      if not seen:
        d += 1
      seenProteins.extend(decoyProteins)
      j += 1
    
    targetProteins = map(lambda x : "_".join(x.split("_")[1:]), table[i][0].split(","))
    seen = False
    for targetProtein in targetProteins:
      if targetProtein in seenProteins:
        seen = True
        break
    if not seen:
      t += 1
    #else:
      #print targetProteins[0]
      #print float(table[i][2]), pickedQvalue, fdr, targetProteins[0], 
    seenProteins.extend(targetProteins)
    pickedQvalue = float(d)/t
    
    proteinName = table[i][0] # Protein name
    curProteinGroup = table[i][1] # Protein group ID
    if i < len(table) - 1:
      nextProteinGroup = table[i+1][1]
    else:
      nextProteinGroup = -1
    
    #if 'absent' in proteinName or proteinName.startswith('mimic'):
    if table[i][4][1] == 'i':
      fpSeen = True
    else:
#    if 'present' in proteinName:
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

    pickedQvalues.append(pickedQvalue)
    fdrs.append(fdr)
    #print float(table[i][2]), pickedQvalue, fdr, targetProteins[0]
  return pickedQvalues, fdrs

def plotQvalues(pickedQvalues, fdrs):
  x = np.linspace(0, 1, num=1000)
  plt.plot(pickedQvalues, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

  plt.figure()
  plt.plot(pickedQvalues, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

def plotQvaluesErrorBar(pickedQvaluesList, fdrsList):
  x = list(np.linspace(0, 1, 1001))
  thresholdFdrs = [] 
  ind = len(x)
  for xi, fdr in zip(pickedQvaluesList, fdrsList): 
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

