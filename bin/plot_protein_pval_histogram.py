import xml.etree.ElementTree as ET
import sys
import matplotlib.pyplot as plt
import numpy as np

ns = {'pout': 'http://per-colator.com/percolator_out/15'}
tree = ET.parse(sys.argv[1])
root = tree.getroot()

proteins = root.find('pout:proteins', ns)
incorrectPvalues, correctPvalues, decoyPvalues = list(), list(), list()
for protein in proteins:
  isDecoy = protein.get('{http://per-colator.com/percolator_out/15}decoy') != "false"
  pvalue = float(protein.find('pout:p_value', ns).text)
  if not isDecoy:
    for peptide_element in protein.iter("{http://per-colator.com/percolator_out/15}peptide_seq"):
      seq = peptide_element.get('seq')
      if seq[1] == 'i':
        incorrectPvalues.append(pvalue)
      else:
        correctPvalues.append(pvalue)
      break
  else:
    decoyPvalues.append(pvalue)

plt.hist(incorrectPvalues + correctPvalues, bins = 20)
plt.title("Target p-values")
plt.xlabel("Fisher p-value")
plt.ylabel("counts")

plt.figure()

plt.hist(incorrectPvalues, bins = 20)
plt.title("Incorrect target p-values")
plt.xlabel("Fisher p-value")
plt.ylabel("counts")

plt.figure()

plt.hist(correctPvalues, bins = 20)
plt.title("Correct target p-values")
plt.xlabel("Fisher p-value")
plt.ylabel("counts")

plt.figure()

plt.hist(decoyPvalues, bins = 20)
plt.title("Decoy p-values")
plt.xlabel("Fisher p-value")
plt.ylabel("counts")

plt.show()
  