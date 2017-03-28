import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy

def plotLogIntensitiesKDE(intx, inty, lims):
  x = np.log10(intx)
  y = np.log10(inty)
  
  _,_,z = getKDE(x,y)
  
  plt.scatter(intx, inty, c = z, edgecolor = '')
  plt.plot(lims,lims,'k--')
  plt.gca().set_xscale("log", nonposx='clip')
  plt.gca().set_yscale("log", nonposy='clip')
  plt.xlim(lims)
  plt.ylim(lims)

def plotDensityScatterKDE(x, y):
  x,y,z = getKDE(x,y)

  plt.scatter(x, y, c=z, s=50, edgecolor='')
  #plt.colorbar()

def getKDE(x,y, sort = False):
  # Calculate the point density
  xy = np.vstack([x,y])
  z = gaussian_kde(xy)(xy)

  # Sort the points by density, so that the densest points are plotted last
  if sort:
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
  return x,y,z

def plotLogDensityScatterHist(xdat, ydat, bins = [100,100]):
  xdat = np.log10(xdat)
  ydat = np.log10(ydat)
  plotDensityScatterHist(xdat, ydat, bins)
  
def plotDensityScatterHist(xdat, ydat, bins = [100,100]):
  xdat = np.array(xdat)
  ydat = np.array(ydat)
  #histogram definition
  xyrange = [[min(xdat),max(xdat)],[min(ydat), max(ydat)]] # data range
  thresh = 1  #density threshold

  hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
  posx = np.digitize(xdat, locx)
  posy = np.digitize(ydat, locy)

  #select points within the histogram
  ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
  hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
  xdat1 = xdat[ind][hhsub < thresh] # low density points
  ydat1 = ydat[ind][hhsub < thresh]
  hh[hh < thresh] = np.nan # fill the areas with low density by NaNs
  
  #plt.figure()
  plt.imshow(np.flipud(hh.T),cmap='jet',extent=np.array(xyrange).flatten(), interpolation='none', origin='upper', aspect='auto')
  plt.colorbar()   
  plt.plot(xdat1, ydat1, '.',color='darkblue')

def prepareSubplots():
  fig = plt.figure()
  fig.subplots_adjust(bottom = 0.05)
  fig.subplots_adjust(top = 0.95)
  fig.subplots_adjust(right = 0.95)
  fig.subplots_adjust(left = 0.05)
