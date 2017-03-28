import xml.etree.ElementTree as ET
import sys
import matplotlib.pyplot as plt
import numpy as np

ns = {'pout': 'http://per-colator.com/percolator_out/15'}
tree = ET.parse(sys.argv[1])
root = tree.getroot()

proteins = root.find('pout:proteins', ns)
fidoQvs, empQvs = [], []
for protein in proteins:
  #isDecoy = protein.get('{http://per-colator.com/percolator_out/15}decoy') != "false"
  #if not isDecoy:
    fidoQvs.append(float(protein.find('pout:q_value', ns).text))
    empQvs.append(float(protein.find('pout:q_value_emp', ns).text))

x = np.linspace(0, 1, num=1000)

plt.plot(fidoQvs, empQvs)
plt.plot(x, x, 'k--')

plt.xlabel("Fido q value")
plt.ylabel("Empirical q value")

if False:
  plt.figure()
  plt.plot(fidoQvs)

  plt.figure()
  plt.plot(empQvs)
plt.show()
  
