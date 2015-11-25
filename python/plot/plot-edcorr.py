#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import itertools as it
import decond.analyze as da
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

default_outbasename = "edcorr"
parser = argparse.ArgumentParser(description="Plot decomposed correlation")
parser.add_argument('corrData', help="correlation data file <c5 or d5>")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
parser.add_argument('-c', '--custom', action='store_true',
                    help="Read the customized parameters in the script")
args = parser.parse_args()

# ======= basic customization ==========
if (args.custom):
    label = ['cation', 'anion']
    threshold = 0  # 1e-6
    tmin, tmax, tstep = 0, 125, 1
    emin, emax, estep = 0, 400, 1
    window = 1
    cbounds = np.arange(-0.05, 0.301, 0.01)
    colorbar_ticks = np.arange(-0.05, 0.301, 0.05)
    cmap = cm.get_cmap('RdYlBu_r')
    cmin, cmax = (-0.3, 0.3)
    # set to None for auto-yticks
    xticks = None
    yticks = None  
# ======================================
    emin //= window
    emax //= window
    estep //= window
    if estep < 1:
        estep = 1
else:
    tmin, tmax, tstep = 0, 400, 2
    emin, emax, estep = 0, 400, 2
    colorbar_ticks = None
    cmap = cm.get_cmap('RdYlBu_r')
    threshold = 0
    cnum = 31

with h5py.File(args.corrData, 'r') as f:
    timeLags = f['timeLags'][...]
    volume = f['volume'][...]
    numMol = f['numMol'][...]
    numIonTypes = numMol.size
    numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2
    decgrp = f[da.DecType.energy.value]
    eBins = decgrp['decBins'][...]  # kcal / mol
    edCorr = decgrp['decCorr'][...]  # nm^2 / ps^2
    edCorr *= (da.const.nano / da.const.angstrom)**2  # AA^2 / ps^2

edf = da.get_edf(args.corrData)[0]
edf *= da.const.angstrom**3 * da.const.calorie  # AA^-3 kcal^-1 mol

# validate arguments
if (args.custom):
    assert(len(label) == numIonTypes)
else:
    label = ['{}'.format(i+1) for i in range(numIonTypes)]

label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

# plot edCorr
rc = {'font': {'size': 46,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': True},
      'legend': {'fontsize': 46},
      'axes': {'labelsize': 46,
               'titlesize': 50},
      'xtick': {'labelsize': 46,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5},
      'ytick': {'labelsize': 46,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5},
      'lines': {'linewidth': 3}}

for key in rc:
    mpl.rc(key, **rc[key])

xlabelpad = 5
ylabelpad = 0.5
clabelpad = 20
spineLineWidth = 1.6
figsize3 = (30, 10)
format = 'eps'


class CustomNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, vanchor=None,
                 clip=False, canchor=0.5):
        self.vanchor = vanchor
        self.canchor = canchor
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vanchor, self.vmax], [0, self.canchor, 1]
        return np.ma.masked_array(np.interp(value, x, y))


edCorr_masked = []
for i in range(numIonTypePairs):
    edf_masked = np.ma.masked_where(np.isnan(edf[i]), edf[i])
    edCorr_masked.append(np.ma.masked_array(
            edCorr[i], np.ones_like(edCorr[i]) *
            np.array(edf_masked < threshold)[:, np.newaxis]))
edCorr_masked = np.ma.masked_array(edCorr_masked)


# ====== Determine the overall data range and color map ========
T, R = np.meshgrid(timeLags[tmin:tmax:tstep], eBins[emin:emax:estep])

if (not args.custom):
    cmin, cmax = (
            np.nanmin(edCorr_masked[:, emin:emax:estep, tmin:tmax:tstep]),
            np.nanmax(edCorr_masked[:, emin:emax:estep, tmin:tmax:tstep]))
    cbounds = np.arange(cmin, cmax, (cmax-cmin)/(cnum+1))

norm = CustomNormalize(vanchor=0, canchor=0.42, vmin=cmin, vmax=cmax)
fig, axs = plt.subplots(1, numIonTypePairs, sharex=True, sharey=True,
                        figsize=figsize3)
for i, (ax, ed) in enumerate(zip(axs.flat, edCorr_masked)):
    c = ax.contourf(T, R, ed[emin:emax:estep, tmin:tmax:tstep],
                    cbounds, norm=norm, cmap=cmap)
    #  ax.contour(T, R, ed[emin:emax:estep, tmin:tmax:tstep], [0],
    #  colors='black')
    ax.set_xlabel(r'$t$\ \ (ps)', labelpad=xlabelpad)
    #    ax.set_title(label[numIonTypes + i])
    plt.sca(ax)
    plt.title(label[numIonTypes + i], y=1.02)
    if (args.custom and xticks is not None):
        plt.xticks(xticks)
    if (i == 0):
        ax.set_ylabel(r'$\epsilon$\ \ (kcal mol$^{-1}$)', labelpad=ylabelpad)
        if (args.custom and yticks is not None):
            ax.set_yticks(yticks)

#  plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
#  wspace=None, hspace=None)
plt.subplots_adjust(left=0.05, bottom=0.15, right=1.05, wspace=0.07)
cb = plt.colorbar(c, ax=axs.ravel().tolist(), ticks=colorbar_ticks, pad=0.01)
cb.set_label(r'$c_{IL}^{(2)}(t;r)$\ \ (\AA$^2$ ps$^{-2}$)', labelpad=clabelpad)
#  plt.tight_layout()

for sp in cb.ax.spines.values():
    sp.set_linewidth(spineLineWidth)
for ax in axs:
    for sp in ax.spines.values():
        sp.set_linewidth(spineLineWidth)

plt.savefig(args.out + '.' + format, bbox_inches="tight", pad_inches=0.15)


# ====== Draw figures one by one ========
# to help determining the appropriate cmin and cmax
if (not args.custom):
    for i, ed in enumerate(edCorr_masked):
        cmin, cmax = (
                np.nanmin(edCorr_masked[i, emin:emax:estep, tmin:tmax:tstep]),
                np.nanmax(edCorr_masked[i, emin:emax:estep, tmin:tmax:tstep]))
        norm = CustomNormalize(vanchor=0, canchor=0.42, vmin=cmin, vmax=cmax)
        plt.figure()
        c = plt.contourf(T, R, ed[emin:emax:estep, tmin:tmax:tstep],
                         32, norm=norm, cmap=cmap)
        # plt.contour(T, R, ed[emin:emax:estep, tmin:tmax:tstep], [0],
        #       colors='black')
        ax = plt.gca()
        ax.set_xlabel(r'$t$\ \ (ps)', labelpad=xlabelpad)
        ax.set_ylabel(r'$\epsilon$\ \ (kcal mol$^{-1}$)', labelpad=ylabelpad)
        ax.set_title(label[numIonTypes + i])
        cb = plt.colorbar(c)
        cb.set_label(r'$c_{IL}^{(2)}(t;r)$\ \ (\AA$^2$ ps$^{-2}$)')
        #    plt.tight_layout()
        plt.savefig(args.out + '-' + str(i) + '.' + format,
                    bbox_inches="tight", pad_inches=0.15)
