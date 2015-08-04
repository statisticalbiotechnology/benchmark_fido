import csv 
import matplotlib.pyplot as plt
import numpy as np

file = open('Filename', 'rb') # The input is the protein output file (.tab) from Percolator (-l)
data = csv.reader(file, delimiter='\t')
table = [row for row in data]
x = np.linspace(0, 1, num=1000)
ent = 0
tg = 0

Pi0 = float(450) / 500

px = []
py = []
pt = []

#fpSeen = False
#pSeen = False


    
for i in range (3, len(table)-1):
#for c in table[3:][0]:
    
    c = table [i][0] # Protein name
    
    if c[0] == 'e':
        continue
       
    if 'mimic' in c:    
        ent = ent + 1   # Entrapment
        
    if 'sp' in c or 'ups' in c: 
        tg = tg + 1     # Target
       
    a = float (ent) / tg # Q value

    fpp = float(ent) * Pi0 / (tg + ent)
    
    pt.append(tg)
    px.append(a)
    py.append(fpp)
    

plt.figure(1)
plt.plot(px, py)
plt.plot(x, x, 'k--')
plt.axis([0, 0.1, 0, 0.1])
plt.xlabel("Expected q value")
plt.ylabel("Observed FDR")
plt.show()

plt.figure(2)
plt.plot(py, pt)
plt.xlabel('Actual FDR')
plt.ylabel('Number of target proteins')
plt.show()
