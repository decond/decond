#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from itertools import accumulate

parser = argparse.ArgumentParser(description="Plot and examine the results from fitSdDCesaro.py and fitNDCesaro.py")
parser.add_argument('-ND', '--NDCesaroFit', help="fitted ND results data file <NDCesaro.fit.h5>")
parser.add_argument('sdDCesaroFit', help="fitted sdD results data file <sdDCesaro.fit.h5>")
parser.add_argument('-o', '--out', default='sdDCesaro.fit', help="output figure base filename, default = 'sdDCesaro.fit'")
parser.add_argument('-T', '--temp', type=float, required=True, help="temperature in K")
args = parser.parse_args()

class Const:
  """
  Defines some constants
  """
  kB = 1.3806488E-23 #(J K-1)
  beta = 1 / (kB * args.temp) #(J-1)
  basicCharge = 1.60217646E-19 #(Coulomb)
  ps = 1.0E-12 #(s)
  nm = 1.0E-9 #(m)
  nm2AA = 10

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(accumulate(range(4))))[idx_r] + idx_c - idx_r

def loadDictFromH5(h5g):
  dict = {}
  def func(k, v):
    dict[k] = v[...]
  h5g.visititems(func)
  return dict

with h5py.File(args.sdDCesaroFit, 'r') as fSdD:
  charge = fSdD.attrs['charge'][...]
  numMol = fSdD.attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
  rBins = fSdD['rBins'][...]
  zzCross = fSdD['zzCross'][...]
  volumeSdD = fSdD['volume'][...]
  sdDCesaro = fSdD['sdDCesaro'][...]
  sdDCesaro_err = fSdD['sdDCesaro_err'][...]
  rho2 = fSdD['rho2'][...]
  rho2_err = fSdD['rho2_err'][...]
  sdD = loadDictFromH5(fSdD['sdD'])
  sdD_err = loadDictFromH5(fSdD['sdD_err'])

if (args.NDCesaroFit != None):
  with h5py.File(args.NDCesaroFit, 'r') as fND:
    zz = fND['zz'][...]
    ND = loadDictFromH5(fND['ND'])
    volumeND = fND['volume'][...]
    Const.ND2ecSI = Const.beta * Const.basicCharge**2 / (volumeND*(Const.nm**3)) * Const.nm**2 / Const.ps

  DI = {}
  sigAutoI = {}
  for fit in ND:
    DI[fit] = ND[fit][:numMol.size] / numMol * Const.nm2AA**2  # AA^2 / ps
    sigAutoI[fit] = ND[fit][:numMol.size] * zz[:numMol.size] * Const.ND2ecSI

density = numMol / volumeSdD
dr = rBins[1] - rBins[0]
dv = 4 * np.pi * rBins**2 * dr
vol = volumeSdD
rho_Vdv = rho2
rho_V = rho_Vdv / dv
rho_dv = rho_Vdv / vol
rho = rho_V / vol
g = rho / np.array([d1 * d2 for (e1, d1) in enumerate(density)
                            for (e2, d2) in enumerate(density) if e2 >= e1]
                  )[:, np.newaxis]


rBins *= Const.nm2AA
figs = []
numPlots = 3 if (args.NDCesaroFit != None) else 2

if (args.NDCesaroFit != None):
  sigIL = {}
  sigI = {}
  for i, fit in enumerate(sdD):
    sigIL[fit] = rho_dv / Const.nm**3 * sdD[fit] * Const.nm**2 / Const.ps * \
                 zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
    sigIL[fit][np.isnan(sigIL[fit])] = 0
    sigIL[fit] = integrate.cumtrapz(sigIL[fit], rBins, initial=0)
    sigI[fit] = sigAutoI[fit][:, np.newaxis] * np.ones_like(rBins)
    for r in range(numMol.size):
      for c in range(r, numMol.size):
        sigI[fit][r] += sigIL[fit][zipIndexPair2(r,c, numMol.size)]
        if (r != c):
          sigI[fit][c] += sigIL[fit][zipIndexPair2(r,c, numMol.size)]

threshold = 0.1
smallRegion = [list(range(35)), list(range(25)), list(range(39))]
for fitKey in sorted(sdD, key=lambda x:x.split(sep='-')[0]):
  fig, axs = plt.subplots(numPlots, 1, sharex=True)
  figs.append(fig)

  sdD[fitKey] *= Const.nm2AA**2

  for i, rdf in enumerate(g):
    axs[0].plot(rBins, rdf, label='{}'.format(i))
  axs[0].legend()
  axs[0].set_title("Fit {} ps".format(fitKey))
  axs[0].set_ylabel(r"$g_{IL}(r)$")

  i = -1
  if (args.NDCesaroFit != None):
    for i, D in enumerate(DI[fitKey]):
      axs[1].plot(rBins, np.ones_like(rBins)*D, label='{}'.format(i), linestyle='--')
      axs[1].legend()
    axs[1].set_ylabel(r"$D_I, D_{IL}(r)$  ($\AA^2$ ps$^{-1}$)")
  else:
    axs[1].set_ylabel(r"$D_{IL}(r)$  ($\AA^2$ ps$^{-1}$)")
  for j, D in enumerate(sdD[fitKey]):
    D[g[j][smallRegion[j]] < threshold] = np.nan
    axs[1].plot(rBins, D, label='{}'.format(i+j+1))
    axs[1].legend()
    axs[1].set_title("threshold {}".format(threshold))

  if (args.NDCesaroFit != None):
    for i, sig in enumerate(sigI[fitKey]):
      axs[2].plot(rBins, sig, label='{}'.format(i))
      axs[2].legend()
    axs[2].set_ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)")
    axs[2].set_xlabel(r"$r$  ($\AA$)")

  plt.xlim(xmax=rBins[rBins.size / 2])

plt.ion()
plt.show()
