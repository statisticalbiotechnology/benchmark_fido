import numpy as np
import scipy.stats as stats

pvalues = [0.11, 0.30]
c = -2* sum(map(np.log, pvalues))
pvalue = 1 - stats.chi2.cdf(c, len(pvalues) * 2)
print pvalue
      
