import xml.etree.ElementTree as ET
import sys
import matplotlib.pyplot as plt
import numpy as np

ns = {'pout': 'http://per-colator.com/percolator_out/15'}
tree = ET.parse(sys.argv[1])
root = tree.getroot()

proteins = root.find('pout:proteins', ns)
fidoQvs, empQvs, trueQvs = [], [], []
prevProteinGroup = -1
p, fp = 0, 0
pi0 = 1.0
for protein in proteins:
  fidoQvs.append(float(protein.find('pout:q_value', ns).text))
  empQvs.append(float(protein.find('pout:q_value_emp', ns).text))
  proteinName = protein.get('{http://per-colator.com/percolator_out/15}protein_id')
  
  curProteinGroup = fidoQvs[-1] # Protein group ID
  
  if proteinName.startswith('absent') or proteinName.startswith('mimic'):
    fpSeen = True
  if proteinName.startswith('present') or proteinName.startswith('sp') or proteinName.startswith('ups'): 
    pSeen = True
  
  if curProteinGroup != prevProteinGroup:
    prevProteinGroup = curProteinGroup
    if pSeen:
      p = p + 1
    elif fpSeen:
      fp = fp + 1 # Actual false positives
    fpSeen = False
    pSeen = False
  
  fdr = pi0 * float(fp) / (p + fp)
  
  trueQvs.append(fdr)
  
x = np.linspace(0, 1, num=1000)

plt.plot(trueQvs, empQvs)
plt.plot(x, x, 'k--')

plt.xlabel("Observed FDR")
plt.ylabel("Empirical q value")
plt.show()
  
