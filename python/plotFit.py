#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate
import itertools as it

parser = argparse.ArgumentParser(description="Plot and examine the results from fitCesaro.py")
parser.add_argument('cesaroFit', help="fitted data file <cesaro.fit.h5>")
parser.add_argument('-o', '--out', default='cesaro.fit', help="output figure base filename, default = 'cesaro.fit'")
parser.add_argument('-T', '--temp', type=float, required=True, help="temperature in K")
parser.add_argument('--threshold', type=float, default=0, help="RDF threshold for D_IL(r) figure, default = 0")
parser.add_argument('--color', nargs='*', help="manually assign line color for each auto and cross terms. "
                                    "<auto1>...<autoN> <cross11>...<cross1N> <cross22>...<cross2N> .. <crossNN>")
parser.add_argument('--label', nargs='*', help="manually assign label for each component. <mol1>...<molN>")
args = parser.parse_args()

threshold = args.threshold

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
  nm2cm = 1.0E-7
  ps2s = 1.0E-12

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(it.accumulate(range(4))))[idx_r] + idx_c - idx_r

def loadDictFromH5(h5g):
  dict = {}
  def func(k, v):
    dict[k] = v[...]
  h5g.visititems(func)
  return dict

with h5py.File(args.cesaroFit, 'r') as f:
  numMol = f.attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2;
  charge = f.attrs['charge'][...]
  timeLags = f['timeLags'][...]
  rBins = f['rBins'][...]
  zzCross = f['zzCross'][...]
  zz = f['zz'][...]
  ww = f['ww'][...]
  volume = f['volume'][...]
  volume_err = f['volume_err'][...]
  nDCesaro = f['nDCesaro'][...]
  nDCesaro_err = f['nDCesaro_err'][...]
  nDCesaroTotal = f['nDCesaroTotal'][...]
  nDCesaroTotal_err = f['nDCesaroTotal_err'][...]
  nD = loadDictFromH5(f['nD'])
  nD_err = loadDictFromH5(f['nD_err'])
  nDTotal = loadDictFromH5(f['nDTotal'])
  nDTotal_err = loadDictFromH5(f['nDTotal_err'])
  sdDCesaro = f['sdDCesaro'][...]
  sdDCesaro_err = f['sdDCesaro_err'][...]
  rho2 = f['rho2'][...]
  rho2_err = f['rho2_err'][...]
  sdD = loadDictFromH5(f['sdD'])
  sdD_err = loadDictFromH5(f['sdD_err'])
  Const.nD2ecSI = Const.beta * Const.basicCharge**2 / (volume*(Const.nm**3)) * Const.nm**2 / Const.ps
  Const.D2AA2_ps = Const.nm2AA**2
  Const.D2cm2_s = Const.nm2cm**2 / Const.ps2s

# validate arguments
if (args.color is not None):
  assert(len(args.color) == numIonTypes + numIonTypePairs )
  mpl.rcParams['axes.color_cycle'] = args.color

def connectLabel(label):
  return label[0] + '-' + label[1]

if (args.label is not None):
  assert(len(args.label) == numIonTypes)
  label = args.label
else:
  label = ['{}'.format(i+1) for i in range(numIonTypes)]
label += [connectLabel(l) for l in it.combinations_with_replacement(label, 2)]

DI = {}
DI_err = {}
sigAutoI = {}
for fit in nD:
  DI[fit] = nD[fit][:numMol.size] / numMol # nm^2 / ps
  DI_err[fit] = nD_err[fit][:numIonTypes] / numMol  # nm^2 / ps
  sigAutoI[fit] = nD[fit][:numMol.size] * zz[:numMol.size] * Const.nD2ecSI

ec = {}
ec_err = {}
ecTotal = {}
ecTotal_err = {}
sortedKeys = sorted(nD.keys(), key=lambda x:x.split(sep='-')[0])
for k in sortedKeys:
  print(k + ':')
  print("========================")
  print("Electrical conductivity in S / m:")
  ec[k] = nD[k] * Const.nD2ecSI * zz
  ec_err[k] = nD_err[k] * Const.nD2ecSI
  print(ec[k])
  print("+/-\n", ec_err[k], sep="")
  ecTotal[k] = Const.nD2ecSI * sum(nD[k]*zz*ww)
  ecTotal_err[k] = Const.nD2ecSI * sum(abs(nD_err[k]))
  print("Total: ", ecTotal[k], " +/- ", ecTotal_err[k], '\n', sep="")
  print("Diffusion constant in 10^-5 cm^2 / s:")
  print(DI[k] * Const.D2cm2_s * 1E5)
  print("+/-\n", DI_err[k] * Const.D2cm2_s * 1E5, '\n', sep="")

print("Total")
for k in sorted(nDTotal.keys(), key=lambda x:x.split(sep='-')[0]):
  print(k + ':')
  print(nDTotal[k] * Const.nD2ecSI, " +/- ", nDTotal_err[k] * Const.nD2ecSI, '\n', sep="")

lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs
# plot fitting results of nCorr
plt.figure()
for i in range(numIonTypes + numIonTypePairs):
  plt.plot(range(len(ecTotal)), [ec[k][i] for k in sortedKeys], linestyle=lineStyle[i], label=label[i])
  if (args.color is None and i == numIonTypes - 1): plt.gca().set_color_cycle(None)

plt.plot(range(len(ecTotal)), [ecTotal[k] for k in sortedKeys], linestyle=':', label='total')
plt.xticks(range(len(ecTotal)), sortedKeys)
plt.legend()
plt.xlabel("fit range  (ps)")
plt.ylabel(r"$\sigma$  (S m$^{-1}$)")

plt.figure()
numErrBars = 5
for i, (nDC, nDC_err)  in enumerate(zip(nDCesaro, nDCesaro_err)):
    plt.errorbar(timeLags, nDC, yerr=nDC_err, errorevery=timeLags.size//numErrBars, linestyle=lineStyle[i], label=label[i])
    if (args.color is None and i == numIonTypes - 1): plt.gca().set_color_cycle(None)
plt.legend(loc='upper left')
plt.xlabel("$\Lambda$  (ps)")
plt.ylabel(r"$\tilde D_I(\Lambda)$, $\tilde D_{IL}(\Lambda)$  ($\AA^2$)")


# plot g-D-sig
density = numMol / volume
cellLengthHalf = volume**(1/3) / 2
dr = rBins[1] - rBins[0]
dv = 4 * np.pi * dr * rBins**2
dvsim = dv.copy()
filter = (rBins > cellLengthHalf) & (rBins < np.sqrt(2) * cellLengthHalf)
dvsim[filter] = 4 * np.pi * dr * rBins[filter] * (3 * cellLengthHalf - 2 * rBins[filter])
filter = (rBins >= np.sqrt(2) * cellLengthHalf)
dvsim[filter] = np.nan
rho_Vdvsim = rho2
rho_V = rho_Vdvsim / dvsim
rho_dvsim = rho_Vdvsim / volume
rho = rho_V / volume

density2 = np.array([d1 * d2 for (e1, d1) in enumerate(density)
                             for (e2, d2) in enumerate(density) if e2 >= e1]
                    )[:, np.newaxis]
g = rho / density2

sigIL = {}
for i, fit in enumerate(sdD):
  sigIL[fit] = rho_dvsim / Const.nm**3 * sdD[fit] * Const.nm**2 / Const.ps * \
               zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
  sigIL[fit][np.isnan(sigIL[fit])] = 0
  sigIL[fit] = integrate.cumtrapz(sigIL[fit], initial=0)

sigI = {}
for i, fit in enumerate(sdD):
  sigI[fit] = sigAutoI[fit][:, np.newaxis] * np.ones_like(rBins)
  for r in range(numIonTypes):
    for c in range(r, numIonTypes):
      sigI[fit][r] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)]
      if (r != c):
        sigI[fit][c] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)]

rBins *= Const.nm2AA
numPlots = 3

smallRegion = []
for rdf in g:
  smallRegion.append(next(i for i, v in enumerate(rdf) if v >= 1))
print("small-rdf region =", smallRegion)

for fitKey in sorted(sdD, key=lambda x:x.split(sep='-')[0]):
  fig, axs = plt.subplots(numPlots, 1, sharex=True)

  # plot rdf
  if (args.color is not None):
    axs[0].set_color_cycle(args.color[numIonTypes:])
  axs[0].axhline(1, linestyle=':', color='black', linewidth=1.0)
  for i, rdf in enumerate(g):
    axs[0].plot(rBins, rdf, label=label[numIonTypes + i])
  axs[0].legend(loc='upper right')
  axs[0].set_title("Fit {} ps".format(fitKey))
  axs[0].set_ylabel(r"$\mathsf{g}_{IL}(r)$")

  # plot D
  DI[fitKey] *= Const.D2AA2_ps
  for i, D in enumerate(DI[fitKey]):
    axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i], linestyle=lineStyle[i])
  axs[1].set_ylabel(r"$D_I$, $D_{IL}(r)$  ($\AA^2$ ps$^{-1}$)")

  sdD[fitKey] *= Const.D2AA2_ps
  for i, D in enumerate(sdD[fitKey]):
    g_masked = np.where(np.isnan(g[i]), -1, g[i])
    D_masked = np.ma.masked_where([c if j <= smallRegion[i] else False
                                   for j, c in enumerate(g_masked < threshold)], D)
    axs[1].plot(rBins, D_masked, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i])
    axs[1].legend(loc='upper right')
    axs[1].set_title("threshold {}".format(threshold))

  # plot sig
  for i, sig in enumerate(sigI[fitKey]):
    axs[2].plot(rBins, sig, label=label[i])
    axs[2].legend(loc='upper right')
  axs[2].set_ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)")
  axs[2].set_xlabel(r"$r$  ($\AA$)")

  plt.xlim(xmax=rBins[rBins.size / 2])

plt.ion()
plt.show()
