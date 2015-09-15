import argparse
import numpy as np
import itertools as it
import decond.analyze as da
import h5py
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

default_outbasename = "g-D-ecdec"
parser = argparse.ArgumentParser(description="Plot rdf-D-ecdec")
parser.add_argument('decond', help="decond analysis file. <decond.d5>")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
args = parser.parse_args()

rc = {'font': {'size': 36,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': True},
      'legend': {'fontsize': 34},
      'axes': {'labelsize': 36},
      'xtick': {'labelsize': 36,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5,
                'minor.size': 4,
                'minor.width': 1.5},
      'ytick': {'labelsize': 36,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5,
                'minor.size': 4,
                'minor.width': 1.5},
      'lines': {'linewidth': 3}
     }

for key in rc:
  mpl.rc(key, **rc[key])

labelpad = 10
spineLineWidth = 1.6
reflinewidth = 1.5

figsize3 = (10, 28)
format='eps'

numIonTypes = 2
numIonTypePairs = numIonTypes * (numIonTypes+1) // 2
lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs

label = ['Na$^+$', 'Cl$^-$']
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

color = ['b', 'g', 'b', 'r', 'g']
mpl.rcParams['axes.color_cycle'] = color

fitKey = 0

g, rBins = da.get_rdf(args.decond)[0:2]
DI, _, _, fit = da.get_diffusion(args.decond)[0:4]
sdD = da.get_decD(args.decond, da.DecType.spatial)[0]
sigI = da.get_ec_dec(args.decond, da.DecType.spatial)[0]

rBins /= da.const.angstrom
DI /= da.const.angstrom**2 / da.const.pico
sdD /= da.const.angstrom**2 / da.const.pico

numPlots = 3
threshold = 0.1

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
axs[0].axhline(1, linestyle=':', color='black', linewidth=reflinewidth)
for i, rdf in enumerate(g):
  axs[0].plot(rBins, rdf, label=label[numIonTypes + i])
axs[0].legend(loc='upper right')
#    axs[0].set_title("Fit {} ps".format(fitKey))
axs[0].set_xlabel(r"$r$\ \ (\AA)", labelpad=labelpad)
axs[0].set_ylabel(r"$\textsl{\textrm{g}}_{IL}(r)$", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(a)', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot D
axs[1].axhline(0, linestyle=':', color='black', linewidth=reflinewidth)
for i, D in enumerate(DI[fitKey]):
  axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i], linestyle=lineStyle[i])

for i, D in enumerate(sdD[fitKey]):
  g_masked = np.where(np.isnan(g[i]), -1, g[i])
  D_masked = np.ma.masked_where([c if j <= smallRegion[i] else False
                                 for j, c in enumerate(g_masked < threshold)], D)
  axs[1].plot(rBins, D_masked, label=label[numIonTypes + i], linestyle=lineStyle[numIonTypes + i])

axs[1].set_xlabel(r"$r$\ \ (\AA)", labelpad=labelpad)
axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(r)$\ \ (\AA$^2$ ps$^{-1}$)", labelpad=labelpad)
#    axs[1].legend(loc='center right')
axs[1].legend(loc=(0.515, 0.245), labelspacing=0.2)
#    axs[1].set_title("threshold {}".format(threshold))
plt.text(abcPos[0], abcPos[1], '(b)', transform=axs[1].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot sig
for i, sig in enumerate(sigI[fitKey]):
  axs[2].plot(rBins, sig, label=label[i])
  axs[2].legend(loc='upper right')
axs[2].set_xlabel(r"$\lambda$\ \ (\AA)", labelpad=labelpad)
axs[2].set_ylabel(r"$\sigma_I(\lambda)$\ \ (S m$^{-1}$)", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(c)', transform=axs[2].transAxes,
         horizontalalignment='left', verticalalignment='top')

axs[0].set_ylim(top=6)
axs[1].set_ylim(bottom=-0.035, top=0.225)
axs[1].set_yticks(np.arange(0.0,0.3,0.1))
axs[2].set_ylim(bottom=2.5, top=8.6)

for ax in axs:
  ax.set_xticks(xticks)
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
  ax.set_xlim(xmax=halfCellLength)
  ax.xaxis.labelpad = 1 
  ax.yaxis.set_label_coords(-0.12, 0.5)
  for sp in ax.spines.values():
    sp.set_linewidth(spineLineWidth)

#plt.tight_layout()
#  plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.subplots_adjust(hspace=0.25)
plt.savefig(args.out + '.' + format, bbox_inches="tight")
