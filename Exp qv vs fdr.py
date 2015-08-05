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
pl0 = float(450) / 500

px = []
py = []
cl = []

fpSeen = False
pSeen = False
    
for i in range (1, len(table)):
    a = table [i][2]
    a = float(a)
    
    c = table [i][0] # Protein name

    if 'absent' in c:
        fpSeen = True
        if fpSeen:
            fp = fp + 1 # Actual false positives
    if 'present' in c: 
        pSeen = True
        if pSeen: 
            p = p + 1
    
    fpSeen = False
    pSeen = False
    fpp = float(fp) / (p + fp)
    # Q-value * Positives = Expected false positives

    px.append(a)
    py.append(fpp)
    
plt.plot(px, py)
plt.plot(x, x, 'k--')
plt.axis([0, 0.1, 0, 0.3])
plt.xlabel("Expected q value")
plt.ylabel("Observed FDR")
plt.show()

