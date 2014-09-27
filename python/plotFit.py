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
parser.add_argument('-p', '--plugin', nargs='*', help="plugin files which will be executed at the end")
parser.add_argument('--nosd', action='store_true', help="no-SD mode, i.e. one-two only mode")
args = parser.parse_args()

threshold = args.threshold

if (not args.nosd):
  with h5py.File(args.cesaroFit, 'r') as f:
    try:
      rBins = f['rBins'][...]
    except KeyError as e:
      print("Warning: no 'rBins' dataset is found in", args.cesaroFit)
      print("Automatically change to --nosd mode")
      args.nosd = True

if (args.out is None):
  if (args.nosd):
    outFilename = 'cesaro-nosd.fit'
  else:
    outFilename = 'cesaro.fit'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

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
  if (not args.nosd):
    rBins = f['rBins'][...]
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

if (args.label is not None):
  assert(len(args.label) == numIonTypes)
  label = args.label
else:
  label = ['{}'.format(i+1) for i in range(numIonTypes)]
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

DI = {}
DI_err = {}
if (not args.nosd):
  sigAutoI = {}
for fit in nD:
  DI[fit] = nD[fit][:numMol.size] / numMol # nm^2 / ps
  DI_err[fit] = nD_err[fit][:numIonTypes] / numMol  # nm^2 / ps
  if (not args.nosd):
    sigAutoI[fit] = nD[fit][:numMol.size] * zz[:numMol.size] * Const.nD2ecSI

ec = {}
ec_err = {}
ecTotal = {}
ecTotal_err = {}
sortedKeys = sorted(nD.keys(), key=lambda x:x.split(sep='-')[0])
for k in sortedKeys:
  print("Fit", k + ' ps:')
  print("========================")
  print("Electrical conductivity in S / m:")
  ec[k] = nD[k] * Const.nD2ecSI * zz
  ec_err[k] = nD_err[k] * Const.nD2ecSI
  print(ec[k])
  print("+/-\n", ec_err[k], sep="")
  ecTotal[k] = Const.nD2ecSI * sum(nD[k]*zz*ww)
  ecTotal_err[k] = Const.nD2ecSI * sum(abs(nD_err[k]))
  print("Fit then sum: ", ecTotal[k], " +/- ", ecTotal_err[k], sep="")
  print("Sum then fit: ", nDTotal[k] * Const.nD2ecSI, " +/- ", nDTotal_err[k] * Const.nD2ecSI, '\n', sep="")
  print("Diffusion constant in 10^-5 cm^2 / s:")
  print(DI[k] * Const.D2cm2_s * 1E5)
  print("+/-\n", DI_err[k] * Const.D2cm2_s * 1E5, '\n', sep="")

rc = {'font': {'size': 34,
               'family': 'serif',
               'serif': 'Times New Roman'},
      'text': {'usetex': True},
      'legend': {'fontsize': 34},
      'axes': {'labelsize': 34},
      'xtick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1},
      'ytick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1},
      'lines': {'linewidth': 3}
     }

for key in rc:
  mpl.rc(key, **rc[key])

labelpad = 10
spineLineWidth = 1.1

figsize1 = (12, 10)
figsize3 = (11, 24)
format='eps'

lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs
# plot fitting results of nCorr
plt.figure(figsize=figsize1)
for i in range(numIonTypes + numIonTypePairs):
  plt.plot(range(len(ecTotal)), [ec[k][i] for k in sortedKeys],
           linestyle=lineStyle[i], label=label[i])
  if (args.color is None and i == numIonTypes - 1): plt.gca().set_color_cycle(None)

plt.plot(range(len(ecTotal)), [ecTotal[k] for k in sortedKeys], linestyle=':', label='total')
plt.xticks(range(len(ecTotal)), sortedKeys)
plt.legend()
plt.xlabel("fit range\ \ (ps)", labelpad=labelpad)
plt.ylabel(r"$\sigma$\ \ (S m$^{-1}$)", labelpad=labelpad)
for sp in plt.gca().spines.values():
  sp.set_linewidth(spineLineWidth)
#plt.tight_layout()
plt.savefig('ec-fitrange.' + format, bbox_inches="tight")

plt.figure(figsize=figsize1)
numErrBars = 5
for i, (nDC, nDC_err)  in enumerate(zip(nDCesaro, nDCesaro_err)):
    plt.errorbar(timeLags, nDC, yerr=nDC_err, errorevery=timeLags.size//numErrBars,
                 linestyle=lineStyle[i], label=label[i])
    if (args.color is None and i == numIonTypes - 1): plt.gca().set_color_cycle(None)

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
plt.legend(handles, labels, loc='upper left', fontsize=31, labelspacing=0.2,
           borderpad=0.2, handletextpad=0.4)
plt.xlabel("$\Lambda$\ \ (ps)", labelpad=labelpad)
plt.ylabel(r"$\tilde D^{(1)}_I(\Lambda)$, $\tilde D^{(2)}_{IL}(\Lambda)$\ \ (\AA$^2$)", labelpad=labelpad)
for sp in plt.gca().spines.values():
  sp.set_linewidth(spineLineWidth)
plt.xlim([0, 1000])
ax.xaxis.labelpad = 6
ax.yaxis.labelpad = 1
plt.savefig('cesaro.' + format, bbox_inches="tight")


if (not args.nosd):
  # plot g-D-sig
  # TODO: when dealing with ion pair with unequal charges, extra care is needed
  # for the numMolPair, since N_I * (N_L - 1) != N_L * (N_I - 1) in that case
  numMolPair = np.array([n1 * (n2-1) if e1 == e2 else n1*n2
                         for (e1, n1) in enumerate(numMol)
                         for (e2, n2) in enumerate(numMol) if e2 >= e1])
  density2 = numMolPair / volume**2
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

  g = rho / density2[:, np.newaxis]

  # determine D_IL(\infty)
  sdDInfty = {}
  sigNonLocal = {}
  halfCellIndex = rBins.size / np.sqrt(3)
  aveWidth = 0.25 / dr
  for fit in sdD:
    sdDInfty[fit] = np.mean(sdD[fit][:, halfCellIndex - aveWidth: halfCellIndex + aveWidth], axis=1)
#    sdDInfty[fit] = np.zeros_like(sdDInfty[fit])
    sigNonLocal[fit] = numMolPair / (volume * Const.nm**3) * sdDInfty[fit] * Const.nm**2 / Const.ps *\
                       zzCross * Const.beta * Const.basicCharge**2

  sigIL = {}
  for i, fit in enumerate(sdD):
#    sigIL[fit] = rho_dvsim / Const.nm**3 * sdD[fit] * Const.nm**2 / Const.ps * \
#                 zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
    sigIL[fit] = rho_dvsim / Const.nm**3 * (sdD[fit] - sdDInfty[fit][:, np.newaxis]) *\
                 Const.nm**2 / Const.ps * zzCross[:, np.newaxis] * Const.beta * Const.basicCharge**2
    sigIL[fit][np.isnan(sigIL[fit])] = 0
    sigIL[fit] = integrate.cumtrapz(sigIL[fit], initial=0)

  sigI = {}
  for i, fit in enumerate(sdD):
    sigI[fit] = sigAutoI[fit][:, np.newaxis] * np.ones_like(rBins)
    for r in range(numIonTypes):
      for c in range(r, numIonTypes):
        sigI[fit][r] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)] +\
                        sigNonLocal[fit][zipIndexPair2(r,c, numIonTypes)]
        if (r != c):
          sigI[fit][c] += sigIL[fit][zipIndexPair2(r,c, numIonTypes)] +\
                          sigNonLocal[fit][zipIndexPair2(r,c, numIonTypes)]

  rBins *= Const.nm2AA
  numPlots = 3

  smallRegion = []
  for rdf in g:
    smallRegion.append(next(i for i, v in enumerate(rdf) if v >= 1))
  print("small-rdf region =", smallRegion)

  for fitKey in sorted(sdD, key=lambda x:x.split(sep='-')[0]):
    fig, axs = plt.subplots(numPlots, 1, sharex=True, figsize=figsize3)

    # plot rdf
    if (args.color is not None):
      axs[0].set_color_cycle(args.color[numIonTypes:])
    axs[0].axhline(1, linestyle=':', color='black', linewidth=1.0)
    for i, rdf in enumerate(g):
      axs[0].plot(rBins, rdf, label=label[numIonTypes + i])
    axs[0].legend(loc='upper right')
#    axs[0].set_title("Fit {} ps".format(fitKey))
    axs[0].set_ylabel(r"$\mathsf{g}_{IL}(r)$", labelpad=labelpad)

    # plot D
    DI[fitKey] *= Const.D2AA2_ps
    axs[1].axhline(0, linestyle=':', color='black', linewidth=1.0)
    for i, D in enumerate(DI[fitKey]):
      axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i], linestyle=lineStyle[i])

    if (args.color is None): axs[1].set_color_cycle(None)
    sdD[fitKey] *= Const.D2AA2_ps
    for i, D in enumerate(sdD[fitKey]):
      g_masked = np.where(np.isnan(g[i]), -1, g[i])
      D_masked = np.ma.masked_where([c if j <= smallRegion[i] else False
                                     for j, c in enumerate(g_masked < threshold)], D)
      axs[1].plot(rBins, D_masked, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i])

    axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(r)$  ($\AA^2$ ps$^{-1}$)", labelpad=labelpad)
#    axs[1].legend(loc='center right')
    axs[1].legend(loc=(0.48, 0.2))
#    axs[1].set_title("threshold {}".format(threshold))

    # plot sig
    for i, sig in enumerate(sigI[fitKey]):
      axs[2].plot(rBins, sig, label=label[i])
      axs[2].legend(loc='upper right')
    axs[2].set_ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)", labelpad=labelpad)
    axs[2].set_xlabel(r"$r$  ($\AA$)", labelpad=labelpad)

    plt.xlim(xmax=cellLengthHalf*Const.nm2AA)

    axs[1].set_ylim(bottom=-0.0005, top=0.0040)
    axs[1].set_yticks(np.arange(-0.000,0.0041,0.001))

    for ax in axs:
      for sp in ax.spines.values():
        sp.set_linewidth(spineLineWidth)

#    plt.tight_layout()
#    plt.savefig('g-D-sig.' + fitKey + '.' + format)

plt.ion()
#plt.show()

# execute plugin scripts
if (args.plugin is not None):
  for plug in args.plugin:
    with open(plug) as f:
      code = compile(f.read(), plug, 'exec')
      exec(code)
