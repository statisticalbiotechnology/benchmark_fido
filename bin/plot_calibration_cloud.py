import csv
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy import stats

import scatter

# python plot_calibration_cloud.py /media/storage/mergespec/data/112111-Human-2h/percolator_tdc_shuffle1/tab/112111-Human-2h.percolator.decoys.tab.psms /media/storage/mergespec/data/112111-Human-2h/percolator_tdc_shuffle2/tab/112111-Human-2h.percolator.decoys.tab.psms

# Using the original sqt files from the Percolator publication:
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_publication/tab/yeast-01.shuffled.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_publication/tab/yeast-01.shuffled2.percolator.decoys.tab.psms

# Comet non-specific digestion:
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_shuffle1/tab/Percolator_original_publication.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_shuffle2/tab/Percolator_original_publication.percolator.decoys.tab.psms

# Comet non-specific digestion; tryptic features removed:
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_shuffle1_no_tryptic_feats/tab/Percolator_original_publication.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_shuffle2_no_tryptic_feats/tab/Percolator_original_publication.percolator.decoys.tab.psms

# Comet tryptic
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_tryptic_shuffle1/tab/Percolator_original_publication.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_comet_tdc_tryptic_shuffle2/tab/Percolator_original_publication.percolator.decoys.tab.psms

# Tide tryptic
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_shuffle_seed2/tab/Percolator_original_publication.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_shuffle_seed3/tab/Percolator_original_publication.percolator.decoys.tab.psms

# Tide semi-tryptic
# python plot_calibration_cloud.py /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_semi_tryptic_shuffle_seed2/tab/Percolator_original_publication.percolator.decoys.tab.psms /media/storage/mergespec/data/Percolator_original_publication/percolator_tdc_semi_tryptic_shuffle_seed3/tab/Percolator_original_publication.percolator.decoys.tab.psms

svmScoresTarget1 = dict()
targetReader = csv.reader(open(sys.argv[1], 'rb'), delimiter = '\t')
targetReader.next()
for row in targetReader:
  svmScoresTarget1[row[0]] = float(row[1])

svmScoresTarget2 = dict()
if len(sys.argv) > 2:
  targetReader2 = csv.reader(open(sys.argv[2], 'rb'), delimiter = '\t')
  targetReader2.next()
  for row in targetReader2:
    svmScoresTarget2[row[0]] = float(row[1])

svmScorePairs = list()
for psmId in svmScoresTarget1:
  if psmId in svmScoresTarget2:
    svmScorePairs.append((svmScoresTarget1[psmId], svmScoresTarget2[psmId]))

scores1, scores2 = zip(*svmScorePairs)

print "Pearson correlation:", np.corrcoef(scores1, scores2)
print "Kruskal-Wallis test (groups = 2):", stats.kruskal(scores1, scores2)
print "Kruskal-Wallis test (groups = numSpectra):", stats.kruskal(*svmScorePairs)

scatter.plotDensityScatterHist(scores1, scores2)
plt.xlabel("Decoy SVM score shuffle1", fontsize = 16)
plt.ylabel("Decoy SVM score shuffle2", fontsize = 16)
plt.axis('equal')
plt.xlim([-3, 1])
plt.ylim([-3, 1])
plt.show()
