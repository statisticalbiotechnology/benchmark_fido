#!/usr/bin/python

import csv
import sys
import math

reader = csv.reader(open(sys.argv[1],'r'), delimiter = '\t')
headers = reader.next()
reader.next()
colIdx = headers.index("lnNumSP")
count = 0
numSp = 0
for row in reader:
  if int(row[1]) == 1:
    numSp += math.exp(float(row[colIdx]))
    count += 1

print float(numSp) / count
  
