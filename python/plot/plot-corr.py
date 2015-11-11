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
parser.add_argument('-c', '--custom', action='store_true',
                    help="Read the customized parameters in the script")
args = parser.parse_args()

# ======= basic customization ==========
if (args.custom):
    label = ['cation', 'anion']
    color = ['b', 'g', 'b', 'r', 'g']
    xmax = 3
# ======================================
else:
    xmax = 1

with h5py.File(args.corrData, 'r') as f:
    timeLags = f['timeLags'][...]
    nCorr = f['nCorr'][...]  # nm^2 / ps^2
    nCorr *= (da.const.nano / da.const.angstrom)**2  # AA^2 / ps^2
    numMol = f['numMol'][...]
    numIonTypes = numMol.size
    numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2

# validate arguments
if (args.custom):
    assert(len(color) == numIonTypes + numIonTypePairs)
    mpl.rcParams['axes.color_cycle'] = color
    assert(len(label) == numIonTypes)
else:
    label = ['{}'.format(i+1) for i in range(numIonTypes)]
label += ['-'.join(l) for l in it.combinations_with_replacement(label, 2)]

# plot nCorr
rc = {'font': {'size': 34,
               'family': 'serif',
               'serif': 'Times'},
      'text': {'usetex': True},
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

for key in rc:
    mpl.rc(key, **rc[key])

xlabelpad = 5
ylabelpad = 0.5
spineLineWidth = 1.6

figsize1 = (10, 8.3)
format = 'eps'

lineStyle = ['--'] * numIonTypes + ['-'] * numIonTypePairs
plt.figure(figsize=figsize1)
plt.gca().axhline(0, linestyle=':', color='black', linewidth=1.0)
for i, corr in enumerate(nCorr):
    plt.plot(timeLags, corr, label=label[i], linestyle=lineStyle[i])

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
