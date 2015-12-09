import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np

file = open(sys.argv[1], 'rb') # The input is the peptide output file (.tab) from Percolator (-r)
data = csv.reader(file, delimiter='\t')
table = [row for row in data]
x = np.linspace(0, 1, num=1000)
fp = 0
p = 0

px = []
py = []
cl = []

for i in range (1, len(table)):
  peptide = table[i][4]
  tag = peptide[3]
  if peptide[3] == 'c':
    p = p + 1
  else:
    fp = fp + 1
  
  fdr = float(fp) / (p + fp)
  # Q-value * Positives = Expected false positives

  px.append(float(table[i][2]))
  py.append(fdr)

plt.plot(px, py)
plt.plot(x, x, 'k--')
#plt.axis([0, 0.1, 0, 0.3])
plt.xlabel("Percolator peptide q value")
plt.ylabel("Observed peptide FDR")

plt.figure()
plt.plot(px, py)
plt.plot(x, x, 'k--')
plt.axis([0, 0.1, 0, 0.1])
plt.xlabel("Percolator peptide q value")
plt.ylabel("Observed peptide FDR")
plt.show()

