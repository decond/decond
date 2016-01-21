#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import itertools as it
import decond.analyze as da
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import interpolate

default_outbasename = "rdf-D-ecdec"
parser = argparse.ArgumentParser(description="Plot rdf-D-ecdec")
parser.add_argument('decond', help="decond analysis file. <decond.d5>")
parser.add_argument('--decond_D', metavar='DECOND',
                    help="decond analysis file for plotting D. <decond.d5>")
parser.add_argument('--decond_ecdec', metavar='DECOND',
                    help="decond analysis file for plotting ecdec."
                         " <decond.d5>")
parser.add_argument('--smooth', action='store_true',
                    help="smooth the data")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
args = parser.parse_args()

# ===================== customization =======================
# set usetex to False
# if UnicodeDecodeError occurs or the output eps is blank
usetex = True

# set which fitting results to plot
# only meaningful when multiple fit ranges are included in decond.d5
fitkey = 0

# e.g. label = ['cation', 'anion']
label = None

# the oder of color is
# [ auto-1, ..., auto-N,
#   cross-11, cross-12, cross-13, ..., cross-1N,
#             cross-22, cross-23, ..., cross-2N,
#                       ... ... ... ... ... ...
#                                      cross-NN ]
#
# e.g. color = ['b', 'g', 'b', 'r', 'g']
#
# if the provided number of colors is not enough,
# the pattern will be repeated for the rest of terms
#
# set to None for default color list (may be ugly)
# see available colors: http://matplotlib.org/api/colors_api.html
color = None

# data(r) will not be plotted if g(r) < threshold
threshold = 0  # e.g. threshold = 0.1

# the plotting range of x-axis, None for auto
xmin = None  # e.g. xmin = 0
xmax = None  # e.g. xmax = 3

rdf_top = None     # rdf_top = 2.5
D_top = None       # D_top = 0.004
D_bottom = None    # D_bottom = -0.001
sig_top = None     # sig_top = 0.75
sig_bottom = None  # sig_bottom = 0

# ticks for x-axis
xticks = None  # xticks = np.arange(0, 21, 5)
xticks_minor = None

# ticks for y-axis
rdf_yticks = None
D_yticks = None  # D_yticks = np.arange(0, 2.5, 0.5)
sig_yticks = None

rdf_legend_loc = None  # rdf_legend_loc = 'upper right'
D_legend_loc = None    # D_legend_loc = 'upper right'
sig_legend_loc = None  # sig_legend_loc = 'upper right'

# set to None to plot all components
# or set to a list to select certain indexes
# such as: rdf_plot_list = [0, 2]
# which plots the 0th and 2nd compondent of rdf
rdf_plot_list = None
DI_plot_list = None
sdD_plot_list = None
sig_plot_list = None

xlabelpad = 1  # controls the distance between x-axis and x-axis label
ylabel_coord = (-0.18, 0.5)  # relative position of ylabel

spineLineWidth = 1.6  # line widith of bouding box
reflinewidth = 1.5  # line width of zero-reference line

figsize3 = (10, 28)  # figure size (width, height)
format = 'eps'

# relative position of (a) (b) (c) labels within each sub-figure
abc_pos = (0.03, 0.965)

# smoothing method
# http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.interp1d.html
# ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’
smooth = 'cubic'
num_smooth_point = 500

# other adjustment
rc = {'font': {'size': 36,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': usetex},
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
      'lines': {'linewidth': 3},
      'savefig': {'transparent': True}}
# ===========================================================

for key in rc:
    mpl.rc(key, **rc[key])

with h5py.File(args.decond, 'r') as f:
    numMol = f['numMol'][...]
    numIonTypes = numMol.size

numIonTypePairs = numIonTypes * (numIonTypes+1) // 2

if color is None:
    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

while len(color) < numIonTypes + numIonTypePairs:
    color += color

assert(len(color) >= numIonTypes + numIonTypePairs)

if label is None:
    label = ['{}'.format(i+1) for i in range(numIonTypes)]

assert(len(label) == numIonTypes)

label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs

if (args.decond_D is None):
    decond_D = args.decond
else:
    decond_D = args.decond_D

if (args.decond_ecdec is None):
    decond_ecdec = args.decond
else:
    decond_ecdec = args.decond_ecdec

g, rBins = da.get_rdf(args.decond)[0:2]
DI, _, _, fit = da.get_diffusion(decond_D)[0:4]
sdD, _, _, rBins_sdD = da.get_decD(decond_D, da.DecType.spatial)[0:4]
g_sdD = da.get_rdf(decond_D)[0]
sigI, _, rBins_sigI = da.get_ec_dec(decond_ecdec, da.DecType.spatial)[0:3]

rBins /= da.const.angstrom
rBins_sdD /= da.const.angstrom
rBins_sigI /= da.const.angstrom
DI /= da.const.angstrom**2 / da.const.pico
sdD /= da.const.angstrom**2 / da.const.pico

numPlots = 3

halfCellIndex = rBins.size / np.sqrt(3)
halfCellLength = rBins[halfCellIndex]

fig, axs = plt.subplots(numPlots, 1, sharex=False, figsize=figsize3)

# plot rdf
if rdf_plot_list is None:
    rdf_plot_list = list(range(numIonTypePairs))

axs[0].axhline(1, linestyle=':', color='black', linewidth=reflinewidth)

for i, rdf in enumerate(g):
    if i in rdf_plot_list:
        axs[0].plot(rBins, rdf, label=label[numIonTypes + i],
                    color=color[numIonTypes + i])

axs[0].legend(loc=rdf_legend_loc)
#    axs[0].set_title("Fit {} ps".format(fitkey))
axs[0].set_xlabel(r"$r$\ \ (\AA)")
axs[0].set_ylabel(r"$\textsl{\textrm{g}}_{IL}(r)$")
plt.text(abc_pos[0], abc_pos[1], '(a)', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot D
axs[1].axhline(0, linestyle=':', color='black', linewidth=reflinewidth)

if DI_plot_list is None:
    DI_plot_list = list(range(numIonTypes))

for i, D in enumerate(DI[fitkey]):
    if i in DI_plot_list:
        axs[1].plot(rBins, np.ones_like(rBins)*D, label=label[i],
                    linestyle=lineStyle[i], color=color[i])

if sdD_plot_list is None:
    sdD_plot_list = list(range(numIonTypePairs))

for i, D in enumerate(sdD[fitkey]):
    if i in sdD_plot_list:
        g_masked = np.where(np.isnan(g_sdD[i]), -1, g_sdD[i])
        idx_threshold = next(i for i, g in enumerate(g_masked) if g >= threshold)

        _rBins_sdD = rBins_sdD[idx_threshold:]
        D = D[idx_threshold:]

        not_nan_D = np.logical_not(np.isnan(D))
        _rBins_sdD = _rBins_sdD[not_nan_D]
        D = D[not_nan_D]

        if args.smooth:
            D_interp = interpolate.interp1d(_rBins_sdD, D, kind=smooth)
            _rBins_sdD = np.linspace(_rBins_sdD[0], _rBins_sdD[-1], num_smooth_point)
            D = D_interp(_rBins_sdD)

        axs[1].plot(_rBins_sdD, D, label=label[numIonTypes + i],
                    linestyle=lineStyle[numIonTypes + i],
                    color=color[numIonTypes + i])

axs[1].set_xlabel(r"$r$\ \ (\AA)")
axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(r)$\ \ (\AA$^2$ ps$^{-1}$)")
axs[1].legend(loc=D_legend_loc)
# axs[1].legend(loc=(0.515, 0.245), labelspacing=0.2)
# axs[1].set_title("threshold {}".format(threshold))
plt.text(abc_pos[0], abc_pos[1], '(b)', transform=axs[1].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot sig
if sig_plot_list is None:
    sig_plot_list = list(range(numIonTypes))

for i, sig in enumerate(sigI[fitkey]):
    if i in sig_plot_list:
        axs[2].plot(rBins_sigI, sig, label=label[i], color=color[i])
        axs[2].legend(loc=sig_legend_loc)
axs[2].set_xlabel(r"$\lambda$\ \ (\AA)")
axs[2].set_ylabel(r"$\sigma_I(\lambda)$\ \ (S m$^{-1}$)")
plt.text(abc_pos[0], abc_pos[1], '(c)', transform=axs[2].transAxes,
         horizontalalignment='left', verticalalignment='top')

axs[0].set_ylim(top=rdf_top)
if rdf_yticks is not None:
    axs[0].set_yticks(rdf_yticks)

axs[1].set_ylim(bottom=D_bottom, top=D_top)
if D_yticks is not None:
    axs[1].set_yticks(D_yticks)

axs[2].set_ylim(bottom=sig_bottom, top=sig_top)
if sig_yticks is not None:
    axs[2].set_yticks(sig_yticks)

if xmax is None:
    xmax = halfCellLength
for ax in axs:
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.set_xlim(left=xmin, right=xmax)
    if xticks_minor is not None:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(xticks_minor))
    ax.xaxis.labelpad = xlabelpad
    ax.yaxis.set_label_coords(ylabel_coord[0], ylabel_coord[1])
    for sp in ax.spines.values():
        sp.set_linewidth(spineLineWidth)

# plt.tight_layout()
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
# wspace=None, hspace=None)
plt.subplots_adjust(hspace=0.25)
plt.savefig(args.out + '.' + format, bbox_inches="tight")
