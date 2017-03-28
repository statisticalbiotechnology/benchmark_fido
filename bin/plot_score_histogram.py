import csv
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.stats import gumbel_l
from scipy.stats import gumbel_r

# python plot_score_histogram.py ~/GPRT/data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.IPI_db/percolator_peptides.no_tdc.pout.tsv ~/GPRT/data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.IPI_db/percolator_peptides_decoy.no_tdc.pout.tsv
pepsTarget = []
targetReader = csv.reader(open(sys.argv[1], 'rb'), delimiter = '\t')
targetReader.next()
for row in targetReader:
  pepsTarget.append(float(row[1]))

pepsDecoy = []
if len(sys.argv) > 2:
  decoyReader = csv.reader(open(sys.argv[2], 'rb'), delimiter = '\t')
  decoyReader.next()
  for row in decoyReader:
    pepsDecoy.append(float(row[1]))

numBins = 50
maxScore = max(pepsTarget + pepsDecoy)
minScore = min(pepsTarget + pepsDecoy)
bins = np.arange(minScore, maxScore + 1e-10, (maxScore - minScore) / numBins)
plt.hist(pepsTarget, bins = bins, label = 'Target', color = 'b', normed = False, alpha = 0.5)
plt.hist(pepsDecoy, bins = bins, label = 'Decoy', color = 'r', normed = False, alpha = 0.5)

if False:
  loc, scale = gumbel_l.fit(pepsDecoy)
  print loc, scale
  x = np.linspace(-3, 3, num=1000)

  #plt.plot(x, gumbel_l.pdf(x, loc, scale))
  #plt.plot(x, 0.7*gumbel_l.pdf(x, -0.83, 0.45) + 0.3*gumbel_l.pdf(x, 0.5, 0.7))
  plt.plot(x, gumbel_r.pdf(x, -0.5, 0.3))
  plt.plot(x, gumbel_l.pdf(x, -0.75, 0.3))

plt.show()
