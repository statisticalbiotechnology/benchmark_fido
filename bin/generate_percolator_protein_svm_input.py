import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.cm as cm
import bisect
import collections
import scipy.stats as stats
import subprocess

def main():
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
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_subset_500k/Pandey.percolator")
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_15M/Pandey.percolator")
    percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator")
      
  for percTabBase in percTabBases:
    outputFN = writeProteinFeatures(percTabBase)
    runPercolator(outputFN)

def writeProteinFeatures(percTabBase):
  targetFN = percTabBase + ".tab.peptides"
  decoyFN = percTabBase + ".decoys.tab.peptides"
  
  targetProtFN = percTabBase + ".tab.proteins"
  decoyProtFN = percTabBase + ".decoys.tab.proteins"
  
  outputFN = percTabBase + ".pin.proteins"
  
  headers = ["SpecId","Label","ScanNr"]
    
  #headers.append("BestScore")
  #headers.append("ScoreRatio")
  #headers.append("SecondScore")
  #headers.append("ThirdScore")
  #headers.append("AvgScore")
  #headers.append("logFisher")
  #headers.append("logFisherCut")
  #headers.append("BinomTest")
  headers.append("MultPEP")
  #headers.append("FracBelowThresh")
  #headers.append("NumBelowThresh")
  headers.append("logNumPepts")
  
  headers += ["Peptide","Proteins"]
  
  print targetFN
  if "logFisher" in headers or "logFisherCut" in headers:
    decoyScores = getPercScores(decoyFN)
  else:
    decoyScores = list()
  
  decoyProteinGroups = getProteinGroups(decoyProtFN)
  targetProteinGroups = getProteinGroups(targetProtFN)
  
  decoyProteins = getProteinFeatures(decoyFN, decoyScores, decoyProteinGroups)
  targetProteins = getProteinFeatures(targetFN, decoyScores, targetProteinGroups)
  
  with open(outputFN,'w') as f:
    writer = csv.writer(f, delimiter = '\t')
    writer.writerow(headers)
    #writer.writerow(["DefaultDirection","-","-",0.5,-0.5,0])
    
    print outputFN
    
    idx = 1
    for protein in targetProteins:
      decoyProtein = ",".join(["decoy_" + p for p in protein.split(",")])
      if decoyProtein in decoyProteins:
        writer.writerow(getFeatureRow(protein, targetProteins[protein], idx, 1, headers))
        writer.writerow(getFeatureRow(decoyProtein, decoyProteins[decoyProtein], idx, -1, headers))
        decoyProteins.pop(decoyProtein, None)
      else:
        writer.writerow(getFeatureRow(protein, targetProteins[protein], idx, 1, headers))
      idx += 1
    
    for decoyProtein in decoyProteins:
      writer.writerow(getFeatureRow(decoyProtein, decoyProteins[decoyProtein], idx, -1, headers))
      idx += 1
  return outputFN

def runPercolator(outputFN):  
  print "Running percolator"
  print outputFN + ".log"
  cmd = "percolator -r %s.peptides -B %s.decoys.peptides %s -S 5 > %s.log 2>&1" % (outputFN, outputFN, outputFN, outputFN)
  subprocess.call(cmd, shell = True)
  
def getFeatureRow(protein, features, idx, label, headers):
  outRow = [protein, label, idx]
  
  pvalCutoff = 0.001
  minScore = -3
  sortedScores = sorted([max([minScore, x[0]]) for x in features])
  if len(sortedScores) > 1:
    secondScore = sortedScores[-2]
  else:
    secondScore = minScore
    
  sortedPvalues = sorted([x[1] for x in features])
  
  numPvalsBelowThresh = sum([1 for x in features if x[1] < pvalCutoff])
  numPepts = float(len(features))
  #ks, p = stats.kstest(sortedPvalues, 'uniform', alternative = 'greater')
  
  if "BestScore" in headers:
    outRow += [sortedScores[-1]]
  if "SecondScore" in headers:
    outRow += [secondScore]
  if "ScoreRatio" in headers:
    outRow += [(secondScore-sortedScores[-1])/max([sortedScores[-1], 1.0])]
    #outRow += [(secondScore-sortedScores[-1])/(secondScore+sortedScores[-1]+2*minScore)]
    #outRow += [(secondScore-minScore)/(sortedScores[-1]-minScore)]
  if "ThirdScore" in headers:  
    if len(sortedScores) > 2:
      thirdScore = sortedScores[-3]
    else:
      thirdScore = -3
    outRow += [thirdScore]
  if "AvgScore" in headers:
    outRow += [np.mean(sortedScores)]
  if "logFisher" in headers:
    outRow += [fisherLogP(sortedPvalues)]
  if "logFisherCut" in headers:
    sortedPvaluesThresh = [p/pvalCutoff for p in sortedPvalues if p < pvalCutoff]
    outRow += [fisherLogP(sortedPvaluesThresh)]
  if "BinomTest" in headers:
    if numPvalsBelowThresh > 0:
      binomTest = 1.0 - stats.binom.cdf(numPvalsBelowThresh-1, numPepts, pvalCutoff)
    else:
      binomTest = 1.0
    binomTest = np.log10(binomTest + 1.0)
    outRow += [binomTest]
  if "MultPEP" in headers:
    multPEP = 1.0
    for x in features:
      multPEP *= x[2]
    if multPEP < 1e-307:
      multPEP = 1e-307
    outRow += [np.log10(multPEP)]
  if "NumBelowThresh" in headers:
    outRow += [numPvalsBelowThresh]
  if "FracBelowThresh" in headers:
    outRow += [numPvalsBelowThresh / numPepts]
  if "logNumPepts" in headers:
    outRow += [np.log10(numPepts)]
  
  outRow += ["-.[" + protein + "].-"] + [x[3] + "[" + str(x[0]) + "]" for x in features]
  return outRow
  #return [protein, label, idx, sortedPvalues[0], "-.[" + protein + "].-"] + [x[3] + "[" + x[4] + "]" for x in features]

def fisherLogP(pvalues):
  if len(pvalues) > 0:
    c = -2* sum(map(np.log10, pvalues))
    fisherP = 1 - stats.chi2.cdf(c, len(pvalues) * 2)
  else:
    fisherP = 1.0
  return np.log10(fisherP + 1)
  
def getPercScores(peptFile):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  scores = list()
  for i, row in enumerate(reader):
    score = float(row[1])
    scores.append(score)
  return list(reversed(scores))

def getProteinFeatures(peptFile, decoyScores, proteinGroups):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  proteinFeatures = collections.defaultdict(list)
  for row in reader:
    #if len(row) == 6:
    proteinId = getProteinGroup(row[5:])
    if proteinId in proteinGroups:
      score = float(row[1])
      if len(decoyScores) > 0:
        pvalue = 1.0 - (bisect.bisect_left(decoyScores, float(row[1])) + 0.5) / (len(decoyScores)+1)
      else:
        pvalue = 1.0
      pep = float(row[3])
      proteinFeatures[proteinId].append([score, pvalue, pep, row[4][2:-2], row[2]])
  return proteinFeatures

def getProteinGroup(proteinIds):
  return ",".join(sorted(list(set(proteinIds))))
  
def getProteinGroups(protFile):
  csv.field_size_limit(sys.maxsize)
  reader = csv.reader(open(protFile, 'r'), delimiter = '\t')
  reader.next()
  proteinGroups = set()
  for row in reader:
    proteinGroups.add(",".join(sorted(row[0].split(","))))
  return proteinGroups
    
if __name__ == "__main__":
  main()

