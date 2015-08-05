import csv 
import matplotlib.pyplot as plt
import numpy as np

file = open('Filename.tab', 'rb') # The input is the protein output file (.tab) from Percolator (-l)
data = csv.reader(file, delimiter='\t')
table = [row for row in data]
x = np.linspace(0, 1, num=1000)
fp = 0
p = 0

px = []
py = []
for i in range (1, len(table)):
    a = table [i][3]
    a = float(a)
    # a is the column of q-values 
    
    c = table [i][0]
    if 'absent' in c:
        fp = fp + 1 # Actual false positives
    if  'present' in c:
        p = p + 1
     
    fpp = float(fp) / (fp + p)
    # q-value * positives = expected false positives 
        
    px.append(a)
    py.append(fpp)


plt.plot(px, py)
plt.plot(x, x, 'k--')
plt.axis([0, 0.1, 0, 0.1])
plt.xlabel("Pep")
plt.ylabel("FDR")
plt.show()



