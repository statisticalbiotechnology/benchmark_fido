import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm
import bisect

def main():
  picked = False
  dataSet = "pandey"
  
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
    percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_1M/Pandey.percolator")
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator")
      
  for percTabBase in percTabBases:
    plotProteinFeatures(percTabBase, picked)
  
  plt.show()

def plotProteinFeatures(percTabBase, picked):
  targetFN = percTabBase + ".tab.peptides"
  decoyFN = percTabBase + ".decoys.tab.peptides"
  
  print targetFN
  decoyScores = getPercScores(decoyFN)
  targetScores, targetPvalues = getPvalues(targetFN, decoyScores)
  
  logTargetPvalues = [np.log10(x) for x in targetPvalues]
  plt.plot(targetScores, logTargetPvalues)
  if False:
    minLogPvalue = np.log10(float(10) / len(decoyScores))
    print minLogPvalue, logTargetPvalues[0]
    polyfitTargets = [x for x in logTargetPvalues if x > minLogPvalue]
    pfit = np.polyfit(targetScores[:len(polyfitTargets)], polyfitTargets, 5)
  else:
    polyfitTargets, scoreTargets = zip(*[(x,y) for x,y in zip(logTargetPvalues,targetScores) if x > -4 and x < -3])
    pfit = np.polyfit(scoreTargets, polyfitTargets, 1)
  plt.plot(targetScores, np.polyval(pfit,targetScores),'y', linestyle = 'dashed', linewidth = 3, markeredgecolor = 'none')
  
def getPercScores(peptFile):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  scores = list()
  for row in reader:
    scores.append(float(row[1]))
  return list(reversed(scores))

def getPvalues(peptFile, decoyScores):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  pvalues, scores = list(), list()
  for row in reader:
    scores.append(float(row[1]))
    pvalues.append(1.0 - (bisect.bisect_left(decoyScores, float(row[1])) + 0.5) / (len(decoyScores)+1))
  return list(reversed(scores)), list(reversed(pvalues))
  
if __name__ == "__main__":
  main()

