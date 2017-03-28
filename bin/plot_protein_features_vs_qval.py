import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm

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
  elif dataSet == "pandey":
    percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_15M/Pandey.percolator")
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator")
      
  for percTabBase in percTabBases:
    plotProteinFeatures(percTabBase, picked)
  
  plt.show()

def plotProteinFeatures(percTabBase, picked):
  if picked:
    targetFN = percTabBase + ".tab.picked_proteins"
    decoyFN = percTabBase + ".decoys.tab.picked_proteins"
  else:
    targetFN = percTabBase + ".tab.proteins"
    decoyFN = percTabBase + ".decoys.tab.proteins"
    
  print targetFN
  targetQvalues, targetFeatures = getFeatures(targetFN, picked)
  decoyQvalues, decoyFeatures = getFeatures(decoyFN, picked)
  
  plt.subplot(3,1,1)
  plt.plot(targetQvalues, targetFeatures, 'b.', label = "target")
  plt.plot(decoyQvalues, decoyFeatures, 'r.', label = "decoy")
  
  plt.ylabel("log(number of peptides)")
  plt.xlabel("log(q-value)")
  plt.legend()
  
  plt.subplot(3,1,2)
  plt.plot(targetQvalues, targetFeatures, 'b.', label = "target")
  
  plt.ylabel("log(number of peptides)")
  plt.xlabel("log(q-value)")
  plt.legend()
  
  plt.subplot(3,1,3)
  plt.plot(decoyQvalues, decoyFeatures, 'r.', label = "decoy")
  
  plt.ylabel("log(number of peptides)")
  plt.xlabel("log(q-value)")
  plt.ylim([0,3.5])
  plt.legend()
  
  #plotDensity(targetQvalues, targetFeatures)
  #plotDensity(decoyQvalues, decoyFeatures)

def plotDensity(x, y):
  plt.figure()
  plt.hexbin(x,y)
  plt.colorbar()
  plt.xlim([-4, 0])
  plt.ylim([0, 3.5])

def getFeatures(protFile, picked):
  csv.field_size_limit(sys.maxsize)
  reader = csv.reader(open(protFile, 'r'), delimiter = '\t')
  reader.next()
  features, qvalues = list(), list()
  for row in reader:
    if float(row[2]) > 0.0001:
      if picked:
        features.append(np.log10(len(row[4:])))
        qvalues.append(np.log10(float(row[2])+1e-4))
      else:
        features.append(np.log10(len(row[4].split())))
        qvalues.append(np.log10(float(row[2])))
  return qvalues, features
  
if __name__ == "__main__":
  main()

