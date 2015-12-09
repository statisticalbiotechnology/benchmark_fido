import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np

file = open(sys.argv[1], 'rb') # The input is the protein output file (.tab) from Percolator (-l)
data = csv.reader(file, delimiter='\t')
table = [row for row in data]
x = np.linspace(0, 1, num=1000)
fp = 0
p = 0

proteinName = table[1][0]
if proteinName.startswith('mimic') or proteinName.startswith('sp') or proteinName.startswith('ups'):
  pi0 = 15.0/16
else:
  pi0 = 1.0

px = []
py = []
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
  
  bestPeptideIncorrect = table[i][4] == 'i'
  noCorrectPeptide = sum(1 if p[1] == 'c' for p in table[i][4:] else 0) == 0
  if noCorrectPeptide:
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
    fpSeen = False
    pSeen = False
  
  fdr = pi0 * float(fp) / (p + fp)
  # Q-value * Positives = Expected false positives

  px.append(fidoQvalue)
  py.append(fdr)

plt.plot(px, py)
plt.plot(x, x, 'k--')
#plt.axis([0, 0.1, 0, 0.3])
plt.xlabel("Fido q value")
plt.ylabel("Observed FDR")

plt.figure()
plt.plot(px, py)
plt.plot(x, x, 'k--')
plt.axis([0, 0.1, 0, 0.1])
plt.xlabel("Fido q value")
plt.ylabel("Observed FDR")

plt.show()

