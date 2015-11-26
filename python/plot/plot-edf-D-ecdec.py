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

default_outbasename = "edf-D-ecdec"
parser = argparse.ArgumentParser(description="Plot edf-D-ecdec")
parser.add_argument('decond', help="decond analysis file. <decond.d5>")
parser.add_argument('--decond_D', metavar='DECOND',
                    help="decond analysis file for plotting D. <decond.d5>")
parser.add_argument('--decond_ecdec', metavar='DECOND',
                    help="decond analysis file for plotting ecdec. <decond.d5>")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
parser.add_argument('-c', '--custom', action='store_true',
                    help="Read the customized parameters in the script")
args = parser.parse_args()

# ======= basic customization ==========
if args.custom:
    label = ['cation', 'anion']
    color = ['b', 'g', 'b', 'r', 'g']

    threshold_D = 1e-5
    edf_top = 0.20
    D_top = 0.004
    D_bottom = -0.0005
    sig_top = 3.0
    sig_bottom = -3.0

    # set to None for auto-ticks
    xticks = [-40, -20, 0, 20, 40]
    xticks_minor = 2
    yticks_edf = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2]
    yticks_D = np.arange(0, 0.0041, 0.001)
    edf_legend_loc = 'lower left'
    D_legend_loc = 'upper center'
    sig_legend_loc = 'upper right'
# ======================================
else:
    threshold_D = 0
    edf_legend_loc = 'upper right'
    D_legend_loc = 'upper right'
    sig_legend_loc = 'upper right'

rc = {'font': {'size': 36,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': True},
      'legend': {'fontsize': 24},
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
format = 'eps'

with h5py.File(args.decond, 'r') as f:
    numMol = f['numMol'][...]
    numIonTypes = numMol.size

numIonTypePairs = numIonTypes * (numIonTypes+1) // 2

if (not args.custom):
    label = ['{}'.format(i+1) for i in range(numIonTypes)]
lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

if (args.custom):
    mpl.rcParams['axes.color_cycle'] = color

fitKey = 0

if (args.decond_D is None):
    decond_D = args.decond
else:
    decond_D = args.decond_D

if (args.decond_ecdec is None):
    decond_ecdec = args.decond
else:
    decond_ecdec = args.decond_ecdec

edf, _, eBins = da.get_edf(args.decond)[0:3]
DI, _, _, fit = da.get_diffusion(decond_D)[0:4]
edD, _, _, eBins_edD = da.get_decD(decond_D, da.DecType.energy)[0:4]
edf_edD = da.get_edf(decond_D)[0]
sig_IL, _, eBins_sig = da.get_ec_dec_energy(decond_ecdec)[0:3]

eBins /= da.const.calorie
eBins_edD /= da.const.calorie
eBins_sig /= da.const.calorie
edf *= da.const.angstrom**3 * da.const.calorie
edf_edD *= da.const.angstrom**3 * da.const.calorie
DI /= da.const.angstrom**2 / da.const.pico
edD /= da.const.angstrom**2 / da.const.pico

numPlots = 3

fig, axs = plt.subplots(numPlots, 1, sharex=False, figsize=figsize3)

abcPos = (0.03, 0.965)

# plot edf
if args.custom:
    axs[0].set_color_cycle(color[numIonTypes:])
for i, iedf in enumerate(edf):
    # edf_masked = np.ma.masked_where([np.isnan(e) or e < threshold_edf for e in iedf], iedf)
    edf_masked = np.ma.masked_where(np.isnan(iedf), iedf)
    axs[0].plot(eBins, edf_masked, label=label[numIonTypes + i])
axs[0].legend(loc=edf_legend_loc)
#    axs[0].set_title("Fit {} ps".format(fitKey))
axs[0].set_xlabel(r"$\epsilon$\ \ (kcal mol$^{-1}$)", labelpad=labelpad)
# axs[0].set_ylabel(r"$\rho_{IL}(\epsilon)$\ \ (nm$^{-3}$ J$^{-1}$ mol$^{1}$)", labelpad=labelpad)
axs[0].set_ylabel(r"$\rho_{IL}^{(2)}(\epsilon)$\ \ ($\AA^{-3}$ kcal$^{-1}$ mol)", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(a)', transform=axs[0].transAxes,
         horizontalalignment='left', verticalalignment='top')
axs[0].set_yscale('log')

# plot D
axs[1].axhline(0, linestyle=':', color='black', linewidth=reflinewidth)
for i, D in enumerate(DI[fitKey]):
    axs[1].plot(eBins, np.ones_like(eBins)*D, label=label[i],
                linestyle=lineStyle[i])

for i, D in enumerate(edD[fitKey]):
    D_masked = np.ma.masked_where([np.isnan(e) or e < threshold_D for e in edf_edD[i]], D)
    axs[1].plot(eBins_edD, D_masked, label=label[numIonTypes + i],
                linestyle=lineStyle[numIonTypes + i])

axs[1].set_xlabel(r"$\epsilon$\ \ (kcal mol$^{-1}$)", labelpad=labelpad)
axs[1].set_ylabel(r"$D^{(1)}_I$, $D^{(2)}_{IL}(\epsilon)$\ \ (\AA$^2$ ps$^{-1}$)", labelpad=labelpad)
axs[1].legend(loc=D_legend_loc)
# axs[1].legend(loc=(0.515, 0.245), labelspacing=0.2)
plt.text(abcPos[0], abcPos[1], '(b)', transform=axs[1].transAxes,
         horizontalalignment='left', verticalalignment='top')

# plot sig
# for i, sig in enumerate(sigI[fitKey]):
    # axs[2].plot(eBins_sigI, sig, label=label[i])
    # axs[2].legend(loc=sig_legend_loc)
# axs[2].set_xlabel(r"$\lambda$\ \ (kcal mol$^{-1}$)", labelpad=labelpad)
# axs[2].set_ylabel(r"$\sigma_I(\lambda)$\ \ (S m$^{-1}$)", labelpad=labelpad)
# plt.text(abcPos[0], abcPos[1], '(c)', transform=axs[2].transAxes,
         # horizontalalignment='left', verticalalignment='top')

# sig_IL
if args.custom:
    axs[2].set_color_cycle(color[numIonTypes:])
for i, sig in enumerate(sig_IL[fitKey]):
    axs[2].plot(eBins_sig, sig, label=label[numIonTypes+i])
    axs[2].legend(loc=sig_legend_loc)
axs[2].set_xlabel(r"$\lambda$\ \ (kcal mol$^{-1}$)", labelpad=labelpad)
axs[2].set_ylabel(r"$\sigma_{IL}^{(2)}(\lambda)$\ \ (S m$^{-1}$)", labelpad=labelpad)
plt.text(abcPos[0], abcPos[1], '(c)', transform=axs[2].transAxes,
         horizontalalignment='left', verticalalignment='top')

if args.custom:
    axs[0].set_ylim(top=edf_top)
    if yticks_edf is not None:
        axs[0].set_yticks(yticks_edf)
    axs[1].set_ylim(bottom=D_bottom, top=D_top)
    if yticks_D is not None:
        axs[1].set_yticks(yticks_D)
    axs[2].set_ylim(bottom=sig_bottom, top=sig_top)

for i, ax in enumerate(axs):
    if (args.custom and xticks is not None):
        ax.set_xticks(xticks)
        if xticks_minor is not None:
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(xticks_minor))
    # if i == 2:
        # ax.set_xlim(xmin=eBins_sig[0], xmax=eBins_sig[-1])
    # else:
        # ax.set_xlim(xmin=eBins[0], xmax=eBins[-1])
    ax.set_xlim(xmin=eBins[0], xmax=eBins[-1])
    ax.xaxis.labelpad = 1
    ax.yaxis.set_label_coords(-0.18, 0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(spineLineWidth)

# plt.tight_layout()
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
# wspace=None, hspace=None)
plt.subplots_adjust(hspace=0.25)
plt.savefig(args.out + '.' + format, bbox_inches="tight")
