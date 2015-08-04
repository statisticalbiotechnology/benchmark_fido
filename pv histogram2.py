import xmltodict
import numpy as np
import matplotlib.pyplot as plt

docfile = open('Filename.xml', 'r') # The input is the xml output file from Percolator 
oridoc = docfile.read()
doc = xmltodict.parse(oridoc)
perc = doc['percolator_output']['proteins']['protein']

px = []

for i in range(0, 656):
    trial = perc[i]['@p:protein_id']
    
    if 'Random' in trial:
        pv = perc[i]['p_value']
        pv = float(pv)
        px.append(pv)

plt.hist(px, bins=60)
plt.show()
