import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import itertools as it

rc = {'font': {'size': 34,
               'family': 'serif',
               'serif': 'Times New Roman'},
      'pdf': {'fonttype': 42},
      'text': {'usetex': True},
      'legend': {'fontsize': 34},
      'axes': {'labelsize': 34},
      'xtick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1,
                'minor.size': 4,
                'minor.width': 1},
      'ytick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1,
                'minor.size': 4,
                'minor.width': 1},
      'lines': {'linewidth': 3}
     }



for key in rc:
  mpl.rc(key, **rc[key])

labelpad = 10
spineLineWidth = 1.1

figsize3 = (10, 28)
format='eps'

numIonTypes = 2
numIonTypePairs = numIonTypes * (numIonTypes+1) // 2
lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs

label = ['cation', 'anion']
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

color = ['b', 'g', 'b', 'r', 'g']
mpl.rcParams['axes.color_cycle'] = color

class Const:
  """
  Defines some constants
  """
  kB = 1.3806488E-23 #(J K-1)
  beta = 1 / (kB * 303) #(J-1)
  basicCharge = 1.60217646E-19 #(Coulomb)
  ps = 1.0E-12 #(s)
  nm = 1.0E-9 #(m)
  nm2AA = 10
  nm2cm = 1.0E-7
  ps2s = 1.0E-12
  D2AA2_ps = nm2AA**2
  D2cm2_s = nm2cm**2 / ps2s

fitKey = '400.0-600.0'

w1Data = np.load('g.npz')
rBins_w1 = w1Data['rBins']
g_w1 = w1Data['g']

w3Data = np.load('g-D-sigI.w3.npz')
rBins = w3Data['rBins']
g = w3Data['g']
DI = w3Data['DI'].item()
sdD = w3Data['sdD'].item()
sigI = w3Data['sigI'].item()

numPlots = 3
threshold = 0.1

#halfCellIndex_w1 = rBins_w1.size / np.sqrt(3)
#halfCellLength_w1 = rBins_w1[halfCellIndex_w1]

halfCellIndex = rBins.size / np.sqrt(3)
halfCellLength = rBins[halfCellIndex]

smallRegion = []
for rdf in g:
  smallRegion.append(next(i for i, v in enumerate(rdf) if v >= 1))

fig, axs = plt.subplots(numPlots, 1, sharex=False, figsize=figsize3)

abcPos = (0.03, 0.965)
xticks = np.arange(0, 21, 5)

# plot rdf
axs[0].set_color_cycle(color[numIonTypes:])
axs[0].axhline(1, linestyle=':', color='black', linewidth=1.0)
for i, rdf in enumerate(g_w1):
  axs[0].plot(rBins_w1, rdf, label=label[numIonTypes + i])
axs[0].legend(loc='upper right')
#    axs[0].set_title("Fit {} ps".format(fitKey))
axs[0].set_xlabel(r"$r$  (\AA)", labelpad=labelpad)
axs[0].set_ylabel(r"$\textsl{\textrm{g}}_{IL}(r)$", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(a)', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot D
axs[1].axhline(0, linestyle=':', color='black', linewidth=1.0)
for i, D in enumerate(DI[fitKey]):
  axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i], linestyle=lineStyle[i])

sdD[fitKey]
for i, D in enumerate(sdD[fitKey]):
  g_masked = np.where(np.isnan(g[i]), -1, g[i])
  D_masked = np.ma.masked_where([c if j <= smallRegion[i] else False
                                 for j, c in enumerate(g_masked < threshold)], D)
  axs[1].plot(rBins, D_masked, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i])

axs[1].set_xlabel(r"$r$  (\AA)", labelpad=labelpad)
axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(r)$  (\AA$^2$ ps$^{-1}$)", labelpad=labelpad)
#    axs[1].legend(loc='center right')
axs[1].legend(loc=(0.437, 0.18), labelspacing=0.2)
#    axs[1].set_title("threshold {}".format(threshold))
plt.text(abcPos[0], abcPos[1], '(b)', transform=axs[1].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot sig
for i, sig in enumerate(sigI[fitKey]):
  axs[2].plot(rBins, sig, label=label[i])
  axs[2].legend(loc='upper right')
axs[2].set_xlabel(r"$\lambda$  (\AA)", labelpad=labelpad)
axs[2].set_ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(c)', transform=axs[2].transAxes,
         horizontalalignment='left', verticalalignment='top')

axs[1].set_ylim(bottom=-0.0003, top=0.0043)
axs[1].set_yticks(np.arange(-0.000,0.0041,0.001))
axs[2].set_ylim(top=0.75)

for ax in axs:
  ax.set_xticks(xticks)
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
  ax.set_xlim(xmax=halfCellLength)
  ax.xaxis.labelpad = 1 
  ax.yaxis.set_label_coords(-0.165, 0.5)
  for sp in ax.spines.values():
    sp.set_linewidth(spineLineWidth)

#plt.tight_layout()
#  plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.subplots_adjust(hspace=0.25)
plt.savefig('g-D-sig.' + format, bbox_inches="tight")

plt.ion()
