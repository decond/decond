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

default_outbasename = "sdcorr"
parser = argparse.ArgumentParser(description="Plot decomposed correlation")
parser.add_argument('corrData', help="correlation data file <c5 or d5>")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
parser.add_argument('-s', '--split', action='store_true',
                    help="plot contour maps separately for each component")
args = parser.parse_args()

# ===================== customization =======================
# set usetex to False
# if UnicodeDecodeError occurs or the output eps is blank
usetex = True

# e.g. label = ['cation', 'anion']
label = None

# data(r) will not be plotted if g(r) < threshold
threshold = 0  # e.g. threshold = 0.1

# set to None to plot the whole data (CAUTION: may take a lot of time!)
# unit: ps
# e.g. tmin, tmax = 0, 2.5
tmin, tmax = None, None

# timestep in the unit of frame
tstep_idx = 1

# set to None to plot the whole data
# unit: Angstrom
# e.g. rmin, rmax, rstep = 3.5, 10
rmin, rmax = None, None

# rstep in the unit of frame
rstep_idx = 1

# ticks for axes, set to None for auto-ticks
xticks = None  # xticks = np.arange(0, 2.5, 0.5)
yticks = None

# e.g. cbounds = np.arange(-0.05, 0.301, 0.01)
cbounds = None

# e.g. colorbar_ticks = np.arange(-0.05, 0.301, 0.05)
colorbar_ticks = None

# e.g. crange = 0.3
crange = None

cmap = cm.get_cmap('RdYlBu_r')
cnum = 31

# distance between axis and axis label
xlabelpad = 5
ylabelpad = 0.5
clabelpad = 20

spineLineWidth = 1.6  # line widith of bouding box
figsize3 = (30, 10)  # figure size (width, height)
format = 'eps'

# other adjustment
rc = {'font': {'size': 46,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': usetex},
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
      'lines': {'linewidth': 3},
      'savefig': {'transparent': True}}
# ===========================================================

if crange is None:
    cmax = None
    cmin = None
else:
    cmax = abs(crange)
    cmin = -cmax

with h5py.File(args.corrData, 'r') as f:
    timeLags = f['timeLags'][...]
    volume = f['volume'][...]
    numMol = f['numMol'][...]
    numIonTypes = numMol.size
    numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2
    decgrp = f[da.DecType.spatial.value]
    rBins = decgrp['decBins'][...]  # nm
    rBins *= da.const.nano / da.const.angstrom  # AA
    sdCorr = decgrp['decCorr'][...]  # nm^2 / ps^2
    sdCorr *= (da.const.nano / da.const.angstrom)**2  # AA^2 / ps^2

g = da.get_rdf(args.corrData)[0]

if label is None:
    label = ['{}'.format(i+1) for i in range(numIonTypes)]

assert(len(label) == numIonTypes)

label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

for key in rc:
    mpl.rc(key, **rc[key])


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

smallRegion = []
for rdf in g:
    smallRegion.append(next(i for i, v in enumerate(rdf) if v >= 1))
print("smallRegion =", smallRegion)

sdCorr_masked = []
for i in range(numIonTypePairs):
    g_masked = np.where(np.isnan(g[i]), -1, g[i])
    sdCorr_masked.append(np.ma.masked_array(
            sdCorr[i], np.ones_like(sdCorr[i]) * np.array(
                [c if j <= smallRegion[i] else False
                 for j, c in enumerate(g_masked < threshold)])[:, np.newaxis]))
sdCorr_masked = np.ma.masked_array(sdCorr_masked)

# Determine the overall data range
if tmin is None:
    tmin = timeLags[0]

if tmax is None:
    tmax = timeLags[-1]

tmin_idx = next(
        (i for i, t in reversed(list(enumerate(timeLags))) if t <= tmin), 0)

tmax_idx = next(
        (i for i, t in enumerate(timeLags) if t >= tmax), len(timeLags)-1) + 1

if rmin is None:
    rmin = rBins[0]

if rmax is None:
    rmax = rBins[-1]

rmin_idx = next(
        (i for i, r in reversed(list(enumerate(rBins))) if r <= rmin), 0)

rmax_idx = next(
        (i for i, r in enumerate(rBins) if r >= rmax), len(rBins)-1) + 1

T, R = np.meshgrid(timeLags[tmin_idx:tmax_idx:tstep_idx],
                   rBins[rmin_idx:rmax_idx:rstep_idx])

if (not args.split):
    if cmin is None:
        _cmin = np.nanmin(sdCorr_masked[
                :, rmin_idx:rmax_idx:rstep_idx, tmin_idx:tmax_idx:tstep_idx])
        print("sdcorr min: {}".format(_cmin))
    else:
        _cmin = cmin

    if cmax is None:
        _cmax = np.nanmax(sdCorr_masked[
                :, rmin_idx:rmax_idx:rstep_idx, tmin_idx:tmax_idx:tstep_idx])
        print("sdcorr max: {}".format(_cmax))
    else:
        _cmax = cmax

    if cbounds is None:
        _cbounds = np.linspace(_cmin, _cmax, cnum+1)
    else:
        _cbounds = cbounds

    _cmax = max(abs(_cmin), abs(_cmax))
    _cmin = -_cmax

    norm = CustomNormalize(vanchor=0, canchor=0.42, vmin=_cmin, vmax=_cmax)
    fig, axs = plt.subplots(1, numIonTypePairs, sharex=True, sharey=True,
                            figsize=figsize3)
    for i, (ax, sd) in enumerate(zip(axs.flat, sdCorr_masked)):
        c = ax.contourf(T, R, sd[rmin_idx:rmax_idx:rstep_idx,
                                 tmin_idx:tmax_idx:tstep_idx],
                        _cbounds, norm=norm, cmap=cmap)
        #  ax.contour(T, R, sd[rmin_idx:rmax_idx:rstep_idx,
        #                      tmin_idx:tmax_idx:tstep_idx], [0],
        #             colors='black')
        ax.set_xlabel(r'$t$\ \ (ps)', labelpad=xlabelpad)
        #    ax.set_title(label[numIonTypes + i])
        plt.sca(ax)
        plt.title(label[numIonTypes + i], y=1.02)
        if xticks is not None:
            ax.set_xticks(xticks)
        if (i == 0):
            ax.set_ylabel(r'$r$\ \ (\AA)', labelpad=ylabelpad)
            if yticks is not None:
                ax.set_yticks(yticks)

    #  plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
    #  wspace=None, hspace=None)
    plt.subplots_adjust(left=0.05, bottom=0.15, right=1.05, wspace=0.07)
    cb = plt.colorbar(c, ax=axs.ravel().tolist(), ticks=colorbar_ticks,
                      pad=0.01)
    cb.set_label(r'$c_{IL}^{(2)}(t;r)$\ \ (\AA$^2$ ps$^{-2}$)',
                 labelpad=clabelpad)
    #  plt.tight_layout()

    for sp in cb.ax.spines.values():
        sp.set_linewidth(spineLineWidth)
    for ax in axs:
        for sp in ax.spines.values():
            sp.set_linewidth(spineLineWidth)

    plt.savefig(args.out + '.' + format, bbox_inches="tight", pad_inches=0.15)

else:
    # ====== Draw figures one by one ========
    for i, sd in enumerate(sdCorr_masked):
        if cmin is None:
            _cmin = np.nanmin(sdCorr_masked[
                i, rmin_idx:rmax_idx:rstep_idx, tmin_idx:tmax_idx:tstep_idx])
        else:
            _cmin = cmin

        if cmax is None:
            _cmax = np.nanmax(sdCorr_masked[
                i, rmin_idx:rmax_idx:rstep_idx, tmin_idx:tmax_idx:tstep_idx])
        else:
            _cmax = cmax

        if cbounds is None:
            _cbounds = np.linspace(_cmin, _cmax, cnum+1)
        else:
            _cbounds = cbounds

        _cmax = max(abs(_cmin), abs(_cmax))
        _cmin = -_cmax

        norm = CustomNormalize(vanchor=0, canchor=0.42, vmin=_cmin, vmax=_cmax)
        plt.figure()
        c = plt.contourf(T, R, sd[rmin_idx:rmax_idx:rstep_idx,
                                  tmin_idx:tmax_idx:tstep_idx],
                         _cbounds, norm=norm, cmap=cmap)
        # plt.contour(T, R, sd[rmin_idx:rmax_idx:rstep_idx,
        #                      tmin_idx:tmax_idx:tstep_idx], [0],
        #             colors='black')
        ax = plt.gca()
        ax.set_xlabel(r'$t$\ \ (ps)', labelpad=xlabelpad)
        ax.set_ylabel(r'$r$\ \ (\AA)', labelpad=ylabelpad)
        ax.set_title(label[numIonTypes + i])
        cb = plt.colorbar(c, ticks=colorbar_ticks)
        cb.set_label(r'$c_{IL}^{(2)}(t;r)$\ \ (\AA$^2$ ps$^{-2}$)')
        # plt.tight_layout()
        plt.savefig(args.out + '-' + str(i) + '.' + format,
                    bbox_inches="tight", pad_inches=0.15)
