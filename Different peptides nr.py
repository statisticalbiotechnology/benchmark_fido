import csv 
import matplotlib.pyplot as plt
import pylab
import numpy as np

file = open('Filename.tab', 'rb') # The input is the protein output file (.tab) from Percolator (-l). In this case, each file has a different peptide number. 
data15 = csv.reader(file, delimiter='\t') # Before Percolator, the file was created with the simulation script, each with different --numSpectra value.
table15 = [row for row in data15]

data30 = csv.reader(open('Filename2.tab', 'rb'), delimiter='\t')
table30 = [row for row in data30]

data60 = csv.reader(open('Filename3.tab', 'rb'), delimiter='\t')
table60 = [row for row in data60]

x = np.linspace(0, 1, num=1000)

fp1 = 0
p1 = 0
fp3 = 0
p3 = 0
fp6 = 0
p6 = 0

pl0 = float(450) / 500

px = []
py = []
px3 = []
py3 = []
px6 = []
py6 = []

fpSeen1 = False
pSeen1 = False
fpSeen3 = False
pSeen3 = False
fpSeen6 = False
pSeen6 = False
    
for i in range (1, 5723):
    a = table15 [i][2]
    a = float(a)
    
    c = table15 [i][0] # Protein name
    d = table15 [i][1] # Protein group ID
    j = i + 1
    e = table15 [j][1]
    if 'absent' in c:
        fpSeen1 = True
    if 'present' in c: 
        pSeen1 = True
    if d == e: 
        continue
    
    else:
        if fpSeen1:
            fp1 = fp1 + 1 # Actual false positives
        if pSeen1: 
            p1 = p1 + 1
    
    fpSeen1 = False
    pSeen1 = False
    fpp = float(fp1) / (p1 + fp1)
    # Q-value * Positives = Expected false positives

    px.append(a)
    py.append(fpp)
    
for i in range (1, 8090):
    f = table30 [i][2]
    f = float(f)
    
    g = table30 [i][0] # Protein name
    h = table30 [i][1] # Protein group ID
    j = i + 1
    k = table30 [j][1]
    if 'absent' in g:
        fpSeen3 = True
    if 'present' in g: 
        pSeen3 = True
    if h == k: 
        continue
    
    else:
        if fpSeen3:
            fp3 = fp3 + 1 # Actual false positives
        if pSeen3: 
            p3 = p3 + 1
    
    fpSeen3 = False
    pSeen3 = False
    fpp3 = float(fp3) / (p3 + fp3)
    # Q-value * Positives = Expected false positives

    px3.append(f)
    py3.append(fpp3)
    
for i in range (1, 10876):
    l = table60 [i][2]
    l = float(l)
    
    m = table60 [i][0] # Protein name
    n = table60 [i][1] # Protein group ID
    j = i + 1
    o = table60 [j][1]
    if 'absent' in m:
        fpSeen6 = True
    if 'present' in m: 
        pSeen6 = True
    if n == o: 
        continue
    
    else:
        if fpSeen6:
            fp6 = fp6 + 1 # Actual false positives
        if pSeen3: 
            p6 = p6 + 1
    
    fpSeen6 = False
    pSeen6 = False
    fpp6 = float(fp6) / (p6 + fp6)
    # Q-value * Positives = Expected false positives

    px6.append(l)
    py6.append(fpp6)
    
    
first = plt.plot(px, py, 'b', label='15000 peptides')
second = plt.plot(px3, py3, 'g', label='30000 peptides')
third = plt.plot(px6, py6, 'r', label='60000 peptides')
pylab.legend(loc='upper left')
plt.plot(x, x, 'k--')
#plt.axis([0, 0.4, 0, 0.4])
plt.xlabel("Expected q value")
plt.ylabel("Observed FDR")
plt.show()

