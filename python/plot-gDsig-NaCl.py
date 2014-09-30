import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import itertools as it

rc = {'font': {'size': 28,
               'family': 'serif',
               'serif': 'Times New Roman'},
      'pdf': {'fonttype': 42},
      'text': {'usetex': True},
      'legend': {'fontsize': 28},
      'axes': {'labelsize': 28},
      'xtick': {'labelsize': 28,
                'major.pad': 8,
                'major.size': 8,
                'major.width': 1,
                'minor.size': 4,
                'minor.width': 1},
      'ytick': {'labelsize': 28,
                'major.pad': 8,
                'major.size': 8,
                'major.width': 1,
                'minor.size': 4,
                'minor.width': 1},
      'lines': {'linewidth': 2,
                'markeredgewidth': 2,
                'markersize': 10}
     }


for key in rc:
  mpl.rc(key, **rc[key])

labelpad = 10
spineLineWidth = 1.1

figsize3 = (10, 12)
format='eps'

numIonTypes = 2
numIonTypePairs = numIonTypes * (numIonTypes+1) // 2
#lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs
lineStyle = ['-', '--', '-', '-.', '--']
#marker = ['s', '^', 's', 'x', '^']
marker = [None, None, None, None, None]
#markevery = [(start, 10) for start in np.arange(numIonTypes + numIonTypePairs)*2] 

label = ['Na$^+$', 'Cl$^-$']
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

color = ['k', 'k', 'k', 'k', 'k']
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

fitKey = '80.0-100.0'

data = np.load('g-D-sigI.npz')
rBins = data['rBins']
g = data['g']
DI = data['DI'].item()
sdD = data['sdD'].item()
sigI = data['sigI'].item()

numPlots = 2
threshold = 0.1

#halfCellIndex_w1 = rBins.size / np.sqrt(3)
#halfCellLength_w1 = rBins[halfCellIndex_w1]

halfCellIndex = rBins.size / np.sqrt(3)
halfCellLength = rBins[halfCellIndex]

smallRegion = []
for rdf in g:
  smallRegion.append(next(i for i, v in enumerate(rdf) if v >= 1))

fig, axs = plt.subplots(numPlots, 1, sharex=False, figsize=figsize3)

abcPos = (0.03, 0.05)
xticks = np.arange(0, 21, 5)

# plot rdf
markevery = [(2, 15)]*(numIonTypes + numIonTypePairs)
axs[0].set_color_cycle(color[numIonTypes:])
axs[0].axhline(1, linestyle=':', color='black', linewidth=1.0)
for i, rdf in enumerate(g):
  axs[0].plot(rBins, rdf, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i],
              marker=marker[numIonTypes + i], markevery=markevery[numIonTypes + i], fillstyle='none')
axs[0].legend(loc='upper right', numpoints=1, labelspacing=0.2)
#    axs[0].set_title("Fit {} ps".format(fitKey))
axs[0].set_xlabel(r"$r$\ \ (\AA)", labelpad=labelpad)
axs[0].set_ylabel(r"$\textsl{\textrm{g}}_{IL}(r)$", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(a)', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='bottom')
axs[0].set_ylim(top=6)
plt.text(0.22, 0.94, 'up to 43', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot D
#markevery = [(3, 9)]*(numIonTypes + numIonTypePairs)
#axs[1].axhline(0, linestyle=':', color='black', linewidth=1.0)
#for i, D in enumerate(DI[fitKey]):
#  axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i], linestyle=lineStyle[i],
#              marker=marker[i], markevery=markevery[i], fillstyle='none')
#
#sdD[fitKey]
#for i, D in enumerate(sdD[fitKey]):
#  g_masked = np.where(np.isnan(g[i]), -1, g[i])
#  D_masked = np.ma.masked_where([c if j <= smallRegion[i] else False
#                                 for j, c in enumerate(g_masked < threshold)], D)
#  axs[1].plot(rBins, D_masked, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i],
#              marker=marker[numIonTypes + i], markevery=markevery[numIonTypes + i], fillstyle='none')
#
#axs[1].set_xlabel(r"$r$\ \ (\AA)", labelpad=labelpad)
#axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(r)$\ \ (\AA$^2$ ps$^{-1}$)", labelpad=labelpad)
##    axs[1].legend(loc='center right')
#axs[1].legend(loc=(0.56, 0.16), labelspacing=0.05, numpoints=1, borderpad=0.15)
##    axs[1].set_title("threshold {}".format(threshold))
#plt.text(abcPos[0], abcPos[1], '(b)', transform=axs[1].transAxes,
#         horizontalalignment='left', verticalalignment='bottom')

# plot sig
markevery = [(2, 15)]*(numIonTypes + numIonTypePairs)
for i, sig in enumerate(sigI[fitKey]):
  axs[1].plot(rBins, sig, label=label[i], linestyle=lineStyle[i],
              marker=marker[i], markevery=markevery[i], fillstyle='none')
  axs[1].legend(loc='upper right', numpoints=1, labelspacing=0.2)
axs[1].set_xlabel(r"$\lambda$\ \ (\AA)", labelpad=labelpad)
axs[1].set_ylabel(r"$\sigma_I(\lambda)$\ \ (S m$^{-1}$)", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(b)', transform=axs[1].transAxes,
         horizontalalignment='left', verticalalignment='bottom')
#axs[1].set_ylim(top=0.7)

for ax in axs:
  ax.set_xticks(xticks)
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
  ax.set_xlim(xmax=halfCellLength)
  ax.xaxis.labelpad = 0.6 
#  ax.yaxis.set_label_coords(-0.15, 0.5)
  ax.yaxis.labelpad = 8
  for sp in ax.spines.values():
    sp.set_linewidth(spineLineWidth)

#plt.tight_layout()
#  plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.subplots_adjust(hspace=0.26)
#plt.savefig('g-D-sig.' + format, bbox_inches="tight")
plt.savefig('g-sig.' + format, bbox_inches="tight")

plt.ion()
