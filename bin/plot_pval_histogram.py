import csv
import matplotlib.pyplot as plt
import sys
import numpy as np
import bisect

# python plot_score_histogram.py ~/GPRT/data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.IPI_db/percolator_peptides.no_tdc.pout.tsv ~/GPRT/data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.IPI_db/percolator_peptides_decoy.no_tdc.pout.tsv

targetFile = sys.argv[1]
targetScores = []
targetReader = csv.reader(open(targetFile, 'rb'), delimiter = '\t')
targetReader.next()
for row in targetReader:
  targetScores.append(float(row[1]))

if len(sys.argv) > 2:
  decoyFile = sys.argv[2]
else:
  decoyFile = targetFile.replace(".percolator.", ".percolator.decoys.")
decoyScores = []
decoyReader = csv.reader(open(decoyFile, 'rb'), delimiter = '\t')
decoyReader.next()
for row in decoyReader:
  decoyScores.append(float(row[1]))

reverse = False
if targetScores[0] > targetScores[-1]:
  reverse = True
  decoyScores = decoyScores[::-1]
  targetScores = targetScores[::-1]
  
pvals = list()
for score in targetScores:
  pval = (bisect.bisect_left(decoyScores, score) + 0.5) / (len(decoyScores)+1)
  if reverse:
    pvals.append(1.0 - pval)
  else:
    pvals.append(pval)

plt.hist(pvals, bins = np.arange(0, 1.01, 0.025), normed = True)
plt.ylim([0, 1.5])
plt.show()
