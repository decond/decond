#!/usr/bin/env python3
import argparse
import h5py
import itertools as it
import decond.analyze as da
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

default_outbasename = "corr"
parser = argparse.ArgumentParser(description="Plot correlation")
parser.add_argument('corrData', help="correlation data file <c5 or d5>")
parser.add_argument('-o', '--out', default=default_outbasename,
                    help="output plot file, default <{0}>".format(
                        default_outbasename))
args = parser.parse_args()

# ===================== customization =======================
# set usetex to False
# if UnicodeDecodeError occurs or the output eps is blank
usetex = True

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

xmax = None

# distance between axis and axis label
xlabelpad = 5
ylabelpad = 0.5

spineLineWidth = 1.6  # line widith of bouding box
reflinewidth = 1.0  # line width of zero-reference line

figsize1 = (10, 8.3)  # figure size (width, height)
format = 'eps'

# other adjustment
rc = {'font': {'size': 34,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': usetex},
      'legend': {'fontsize': 34},
      'axes': {'labelsize': 34,
               'titlesize': 34},
      'xtick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5},
      'ytick': {'labelsize': 34,
                'major.pad': 10,
                'major.size': 8,
                'major.width': 1.5},
      'lines': {'linewidth': 3},
      'savefig': {'transparent': True}}
# ===========================================================

with h5py.File(args.corrData, 'r') as f:
    timeLags = f['timeLags'][...]
    nCorr = f['nCorr'][...]  # nm^2 / ps^2
    nCorr *= (da.const.nano / da.const.angstrom)**2  # AA^2 / ps^2
    numMol = f['numMol'][...]
    numIonTypes = numMol.size
    numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2

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

for key in rc:
    mpl.rc(key, **rc[key])

plt.figure(figsize=figsize1)
plt.gca().axhline(0, linestyle=':', color='black', linewidth=reflinewidth)
for i, corr in enumerate(nCorr):
    plt.plot(timeLags, corr, label=label[i], linestyle=lineStyle[i],
             color=color[i])

leg = plt.legend()
plt.xlim(xmax=xmax)
plt.xlabel(r'$t$\ \ (ps)', labelpad=xlabelpad)
plt.ylabel(r'$C_I^{(1)}(t)$, $C_{IL}^{(2)}(t)$\ \ (\AA$^2$ ps$^{-2}$)',
           labelpad=ylabelpad)
# plt.tight_layout()

ax = plt.gca()
for sp in ax.spines.values():
    sp.set_linewidth(spineLineWidth)

plt.savefig(args.out + '.' + format, bbox_inches="tight", pad_inches=0.20)
