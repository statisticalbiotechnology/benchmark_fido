import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import pylab

fp = 0
p = 0
fp3 = 0
p3 = 0
fp6 = 0
p6 = 0
x = np.linspace(0, 1, num=1000)

Pi0 = float(450) / 500
qvalues = []
fdrs = []
px = []
py = []
ex = []
qvalues3 = []
fdrs3 = []
px3 = []
py3 = []
ex3 = []
qvalues6 = []
fdrs6 = []
px6 = []
py6 = []
ex6 = []

for i in range(1,10):
    t = csv.reader(open(sys.argv[1] % (i),'rb'), delimiter = '\t')
    t.next()
    localQvalues = [float(row[2]) for row in t]
    qvalues.append(localQvalues)
    q = csv.reader(open(sys.argv[1] % (i),'rb'), delimiter = '\t')
    q.next()
    falsep = [row[0] for row in q]
    fdrs.append(falsep)
        
minLength = min([len(lq) for lq in qvalues])

for i in range(minLength):
    qvalRow = [qvalues[j][i] for j in range(int(sys.argv[7])-1)]
    avqv = np.mean(qvalRow)
    sdqv = np.std(qvalRow) 
    
    c = fdrs[0][i]
    
    if 'absent' in c:
        fp = fp + 1 
    
    if  'present' in c: 
        p = p + 1
    
    
    fpp = float(fp * Pi0) / (p + fp)
    
    px.append(avqv)
    py.append(fpp)
    ex.append(sdqv)
   
for i in range(1,10):
    t3 = csv.reader(open(sys.argv[3] % (i),'rb'), delimiter = '\t')
    t3.next()
    localQvalues3 = [float(row[2]) for row in t3]
    qvalues3.append(localQvalues3)
    q3 = csv.reader(open(sys.argv[3] % (i),'rb'), delimiter = '\t')
    q3.next()
    falsep3 = [row[0] for row in q3]
    fdrs3.append(falsep3)
        
minLength = min([len(lq) for lq in qvalues3])

for i in range(minLength):
    qvalRow3 = [qvalues3[j][i] for j in range(int(sys.argv[7])-1)]
    avqv3 = np.mean(qvalRow3)
    sdqv3 = np.std(qvalRow3) 
    
    d = fdrs3[0][i]
    
    if 'absent' in d:
        fp3 = fp3 + 1 
    
    if  'present' in d: 
        p3 = p3 + 1
    
    
    fpp3 = float(fp3 * Pi0) / (p3 + fp3)
    
    px3.append(avqv3)
    py3.append(fpp3)
    ex3.append(sdqv3)
       
for i in range(1,10):
    t6 = csv.reader(open(sys.argv[5] % (i),'rb'), delimiter = '\t')
    t6.next()
    localQvalues6 = [float(row[2]) for row in t6]
    qvalues6.append(localQvalues6)
    q6 = csv.reader(open(sys.argv[5] % (i),'rb'), delimiter = '\t')
    q6.next()
    falsep6 = [row[0] for row in q6]
    fdrs6.append(falsep6)
        
minLength = min([len(lq) for lq in qvalues6])

for i in range(minLength):
    qvalRow6 = [qvalues6[j][i] for j in range(int(sys.argv[7])-1)]
    avqv6 = np.mean(qvalRow6)
    sdqv6 = np.std(qvalRow6) 
    
    e = fdrs6[0][i]
    
    if 'absent' in e:
        fp6 = fp6 + 1 
    
    if  'present' in e: 
        p6 = p6 + 1
    
    
    fpp6 = float(fp6 * Pi0) / (p6 + fp6)
    
    px6.append(avqv6)
    py6.append(fpp6)
    ex6.append(sdqv6)
    

plt.plot(x, x, 'k--')
plt.errorbar(px, py, xerr=ex, yerr=None, color='r', ecolor = 'maroon', errorevery = 60, label=sys.argv[2])
plt.errorbar(px3, py3, xerr=ex3, yerr=None, color='g', ecolor = 'darkgreen', errorevery = 60, label=sys.argv[4])
plt.errorbar(px6, py6, xerr=ex6, yerr=None, color='royalblue', ecolor = 'steelblue', errorevery = 80, label=sys.argv[6])
pylab.legend(loc='lower right')
plt.xlabel("Expected q value")
plt.ylabel("Observed FDR")
plt.show()
