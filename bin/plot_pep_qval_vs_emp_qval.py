import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right

xlabel = "Q-values from PEPs"
ylabel = "Q-values from decoys"

def main():
  poutFile = sys.argv[1]
  reader = csv.reader(open(poutFile, 'rb'), delimiter = '\t')
  reader.next()
  pepSum = 0.0
  pepQvals, decoyQvals = list(), list()
  for i, row in enumerate(reader):
    pepSum += float(row[3])
    pepQvals.append(pepSum / (i+1))
    decoyQvals.append(float(row[2]))
  
  plotQvalues(pepQvals, decoyQvals)
  plt.show()
 
def plotQvalues(qv1, qv2):
  x = np.linspace(0, 1, num=1000)
  plt.plot(qv1, qv2)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

  plt.figure()
  plt.plot(qv1, qv2)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)

if __name__ == "__main__":
  main()

