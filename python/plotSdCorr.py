#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from  matplotlib import cm
from itertools import accumulate

parser = argparse.ArgumentParser(description="Plot SD correlation results from spatialDecompose")
parser.add_argument('sdcorrData', help="SD correlation data file <sdcorr.h5>")
parser.add_argument('-o', '--out', default='corr', help="output figure base filename, default = 'sdcorr'")
args = parser.parse_args()

def zipIndexPair(idx_r, idx_c, size):
  """
  Returns the single index based the row index and column index
  """
  return idx_r * size + idx_c

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(accumulate(range(4))))[idx_r] + idx_c - idx_r

with h5py.File(args.sdcorrData, 'r') as f:
  timeLags = f['timeLags'][...]
  rBins = f['rBins'][...]
  sdCorr = f['sdCorr'][...]  # nm^2 / ps^2
  numMol = f.attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2;

sdCorr2 = np.empty([numIonTypePairs, rBins.size, timeLags.size])
for i in range(numIonTypes):
  for j in range(i, numIonTypes):
    if (i == j):
      sdCorr2[zipIndexPair2(i, j, numIonTypes)] = sdCorr[zipIndexPair(i, j, numIonTypes)]
    else:
      sdCorr2[zipIndexPair2(i, j, numIonTypes)] = (sdCorr[zipIndexPair(i, j, numIonTypes)] +
                                                     sdCorr[zipIndexPair(j, i, numIonTypes)]) / 2

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.42, 1]
        return np.ma.masked_array(np.interp(value, x, y))

nm2AA = 10

fig, axs = plt.subplots(1, numIonTypePairs, sharex=True, sharey=True, figsize=(18, 5))
tmin, tmax = 0, 70
rmin, rmax = 25, 120
vmin, vmax = (np.nanmin(sdCorr2[:, rmin:rmax, tmin:tmax]) * nm2AA**2,
              np.nanmax(sdCorr2[:, rmin:rmax, tmin:tmax]) * nm2AA**2)
T, R = np.meshgrid(timeLags[tmin:tmax], rBins[rmin:rmax] * nm2AA)
#bounds = np.linspace(vmin, vmax, endpoint=True)
bounds = np.arange(-0.05, 0.301, 0.025)
#cmap = cm.get_cmap('RdBu_r')
cmap = cm.get_cmap('RdYlBu_r', 28)
norm = MidpointNormalize(midpoint=0, vmin=-0.3, vmax=0.3)

for i, (ax, sd) in enumerate(zip(axs.flat, sdCorr2)):
  if (i < numIonTypePairs):
    c = ax.contourf(T, R, sd[rmin:rmax, tmin:tmax] * nm2AA**2, bounds, norm=norm, cmap=cmap)
    ax.set_xlabel(r'$t$  (ps)')
    ax.set_title('{}'.format(i))
    if (i == 0):
      ax.set_ylabel(r'$r$  ($\AA$)')

plt.tight_layout()
cb = plt.colorbar(c, ax=axs.ravel().tolist(), ticks=np.arange(-0.05, 0.301, 0.05))
cb.set_label(r'$c_{IL}^{(2)}(t;r)$  ($\AA^2$ ps$^{-2}$)')

plt.ion()
plt.show()
