import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import os

def main():
  picked = True
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
    percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc/tab_15M_pval_0.001/Pandey.percolator")
    #percTabBases.append("/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator")
      
  for percTabBase in percTabBases:
    writeProteinFdrs(percTabBase, picked)
    writeProteinFdrsFromProts(percTabBase, picked)

def writeProteinFdrs(percTabBase, picked):
  targetFN = percTabBase + ".tab.peptides"
  decoyFN = percTabBase + ".decoys.tab.peptides"
  if picked:
    targetOutFN = percTabBase + ".tab.picked_proteins"
    decoyOutFN = percTabBase + ".decoys.tab.picked_proteins"
  else:
    targetOutFN = percTabBase + ".tab.classic_proteins"
    decoyOutFN = percTabBase + ".decoys.tab.classic_proteins"
  
  outHeader = ["ProteinId", "ProteinGroupId", "q-value", "posterior_error_prob", "peptideIds"]
  
  print targetFN
  print decoyFN
  targetProteins = getProteinQvalueMap(targetFN)
  decoyProteins = getProteinQvalueMap(decoyFN)
  
  pickedProteins = pickedProteinCompetition(targetProteins, decoyProteins, picked)
  
  targetWriter = csv.writer(open(targetOutFN, 'w'), delimiter = '\t')
  decoyWriter = csv.writer(open(decoyOutFN, 'w'), delimiter = '\t')
  
  targetWriter.writerow(outHeader)
  decoyWriter.writerow(outHeader)
  significantProteins = 0
  for i, proteinQvaluePair in enumerate(pickedProteins):
    outRow = [proteinQvaluePair[0], i+1, proteinQvaluePair[1][0], 1] + proteinQvaluePair[1][1]
    if proteinQvaluePair[0].startswith("decoy"):
      decoyWriter.writerow(outRow)
    else:
      targetWriter.writerow(outRow)
      if proteinQvaluePair[1][0] < 0.01:
        significantProteins += 1
  print "Proteins with q < 0.01:", significantProteins
    
def getProteinQvalueMap(peptFile):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  proteinQvalueMap = dict()
  for row in reader:
    if len(row) == 6:
      if row[5] not in proteinQvalueMap:
        proteinQvalueMap[row[5]] = [float(row[3]), [row[4][2:-2]]]
      else:
        proteinQvalueMap[row[5]][1].append(row[4][2:-2])
  return proteinQvalueMap
  
def pickedProteinCompetition(targetProteins, decoyProteins, picked):
  pickedProteins = list()
  for protein in targetProteins:
    decoyProtein = "decoy_" + protein
    if picked and decoyProtein in decoyProteins:
      if targetProteins[protein][0] <= decoyProteins[decoyProtein][0]:
        pickedProteins.append([protein, targetProteins[protein]])
      else:
        pickedProteins.append([decoyProtein, decoyProteins[decoyProtein]])
      decoyProteins.pop(decoyProtein, None)
    else:
      pickedProteins.append([protein, targetProteins[protein]])
  
  for decoyProtein in decoyProteins:
    pickedProteins.append([decoyProtein, decoyProteins[decoyProtein]])
  
  pickedProteins = sorted(pickedProteins, key = lambda x : x[1][0])
  
  decoys, targets = 0, 0
  for proteinQvaluePair in pickedProteins:
    if proteinQvaluePair[0].startswith("decoy"):
      decoys += 1
    else:
      targets += 1
    proteinQvaluePair[1][0] = float(decoys) / max([1,targets])
  
  lastMin = 1.0
  for proteinQvaluePair in reversed(pickedProteins):
    proteinQvaluePair[1][0] = min([lastMin, proteinQvaluePair[1][0]])
    lastMin = proteinQvaluePair[1][0]
    
  return pickedProteins

def writeProteinFdrsFromProts(percTabBase, picked):
  targetFN = percTabBase + ".tab.proteins"
  decoyFN = percTabBase + ".decoys.tab.proteins"
  if picked:
    targetOutFN = percTabBase + ".tab.picked_proteins_from_prots"
    decoyOutFN = percTabBase + ".decoys.tab.picked_proteins_from_prots"
  else:
    targetOutFN = percTabBase + ".tab.classic_proteins_from_prots"
    decoyOutFN = percTabBase + ".decoys.tab.classic_proteins_from_prots"
  
  outHeader = ["ProteinId", "ProteinGroupId", "q-value", "posterior_error_prob", "peptideIds"]
  
  print targetFN
  print decoyFN
  targetProteins = getProteinQvalueMapFromProts(targetFN)
  decoyProteins = getProteinQvalueMapFromProts(decoyFN)
  
  pickedProteins = pickedProteinCompetitionFromProts(targetProteins, decoyProteins, picked)
  
  targetWriter = csv.writer(open(targetOutFN, 'w'), delimiter = '\t')
  decoyWriter = csv.writer(open(decoyOutFN, 'w'), delimiter = '\t')
  
  targetWriter.writerow(outHeader)
  decoyWriter.writerow(outHeader)
  significantProteins = 0
  for i, outRow in enumerate(pickedProteins):
    outRow[1] = i+1
    if outRow[0].startswith("decoy"):
      decoyWriter.writerow(outRow)
    else:
      targetWriter.writerow(outRow)
      if outRow[2] < 0.01:
        significantProteins += 1
  print "Proteins with q < 0.01:", significantProteins

def getProteinQvalueMapFromProts(protFile):
  csv.field_size_limit(sys.maxsize)
  reader = csv.reader(open(protFile, 'r'), delimiter = '\t')
  reader.next()
  proteinQvalueMap = dict()
  for row in reader:
    for protein in row[0].split(","):
      proteinQvalueMap[protein] = row
  return proteinQvalueMap
    
def pickedProteinCompetitionFromProts(targetProteins, decoyProteins, picked):
  pickedProteins = list()
  for protein in targetProteins:
    decoyProtein = "decoy_" + protein
    if picked and decoyProtein in decoyProteins:
      if float(targetProteins[protein][2]) <= float(decoyProteins[decoyProtein][2]):
        pickedProteins.append(targetProteins[protein])
      else:
        pickedProteins.append(decoyProteins[decoyProtein])
      decoyProteins.pop(decoyProtein, None)
    else:
      pickedProteins.append(targetProteins[protein])
  
  for decoyProtein in decoyProteins:
    pickedProteins.append(decoyProteins[decoyProtein])
  
  pickedProteins = sorted(pickedProteins, key = lambda x : float(x[2]))
  
  decoys, targets = 0, 0
  seenProteins = set()
  newPickedProteins = list()
  for row in pickedProteins:
    if not row[0] in seenProteins:
      seenProteins.add(row[0])
      if row[0].startswith("decoy"):
        decoys += 1
      else:
        targets += 1
      row[2] = float(decoys) / max([1,targets])
      newPickedProteins.append(row)
  
  pickedProteins = newPickedProteins
  lastMin = 1.0
  for proteinQvaluePair in reversed(pickedProteins):
    proteinQvaluePair[2] = min([lastMin, proteinQvaluePair[2]])
    lastMin = proteinQvaluePair[2]
    
  return pickedProteins
  
if __name__ == "__main__":
  main()

