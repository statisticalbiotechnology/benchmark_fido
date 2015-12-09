import sys
import os
import csv 
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

xlabel = "Protein q value"
ylabel = "Frac proteins with no correct peptides"
#ylabel = "Frac mimic proteins"

def main():
  folder = "/home/matthew/benchmark_fido/data/hm_percolator/151009_fido_protein_inference"
  #fileNames = ["103111-Yeast-2hr-01.pout.proteins", "103111-Yeast-2hr-02.pout.proteins", "103111-Yeast-2hr-03.pout.proteins", "103111-Yeast-2hr.pout.proteins"]
  fileNames = ["112111-Human-2h-01.pout.proteins", "112111-Human-2h-02.pout.proteins", "112111-Human-2h-03.pout.proteins", "112111-Human-2h.pout.proteins"]
  proteinLists = dict()
  for fileName in fileNames:
    proteinLists[fileName] = getProteinList(os.path.join(folder,fileName))
    print fileName, len(proteinLists[fileName])
  
  print ""
  
  atLeastOnce = proteinLists[fileNames[0]] | proteinLists[fileNames[1]] | proteinLists[fileNames[2]]
  print "At least once:", len(atLeastOnce)
  print "At least once and in merged:", len(atLeastOnce.intersection(proteinLists[fileNames[3]]))
  
  atLeastTwice = proteinLists[fileNames[0]].intersection(proteinLists[fileNames[1]]) | proteinLists[fileNames[1]].intersection(proteinLists[fileNames[2]]) | proteinLists[fileNames[0]].intersection(proteinLists[fileNames[2]])
  print "At least twice:", len(atLeastTwice)
  print "At least twice and in merged:", len(atLeastTwice.intersection(proteinLists[fileNames[3]]))
  
  atLeastThrice = proteinLists[fileNames[0]].intersection(proteinLists[fileNames[1]]).intersection(proteinLists[fileNames[2]])
  print "At least thrice:", len(atLeastThrice)
  print "At least thrice and in merged:", len(atLeastThrice.intersection(proteinLists[fileNames[3]]))
  
def getProteinList(fileName):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')  
  data.next() # skip header line
  table = [row for row in data]
  proteinList = []

  fpSeen = False
  pSeen = False
  nextProteinGroup = -1
  i = 0
  groupProteins = set()
  
  for row in table:
    proteinQvalue = float(row[2])
    if proteinQvalue > 0.01:
      break
    
    proteinName = row[0] # Protein name
    curProteinGroup = row[1] # Protein group ID
    if i < len(table) - 1:
      nextProteinGroup = table[i+1][1]
    else:
      nextProteinGroup = -1
    
    i += 1
    
    groupProteins.add(proteinName)
    if curProteinGroup == nextProteinGroup: 
      continue 
    else:
      proteinList.append(",".join(sorted(list(groupProteins))))
      groupProteins = set()
    
  return set(proteinList)
   
if __name__ == "__main__":
  main()

