import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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

xlabelpad = 5
ylabelpad = 8 
spineLineWidth = 1.1

figsize = (12, 10)
format='eps'

fitKey = '80.0-100.0'
rBins = np.load('sigI-old.npz')['rBins']
sigI_old = np.load('sigI-old.npz')['sigI'].item()[fitKey]
sigI_new = np.load('sigI-new.npz')['sigI'].item()[fitKey]

labelOld = ['Na$^+$ (original)', 'Cl$^-$ (original) ']
labelNew = ['Na$^+$ (re-computed)', 'Cl$^-$ (re-computed)']

plt.figure(figsize=figsize)

for i, sig in enumerate(sigI_old):
  plt.plot(rBins, sig, label=labelOld[i], linestyle="--")

plt.gca().set_color_cycle(None)

for i, sig in enumerate(sigI_new):
  plt.plot(rBins, sig, label=labelNew[i], linestyle="-")

plt.legend(loc='upper right')
plt.ylabel(r"$\sigma_I(\lambda)$\ \ (S m$^{-1}$)", labelpad=ylabelpad)
plt.xlabel(r"$r$\ \ (\AA)", labelpad=xlabelpad)

halfCellIndex = rBins.size / np.sqrt(3)
plt.xlim(xmax=rBins[halfCellIndex])

for sp in plt.gca().spines.values():
  sp.set_linewidth(spineLineWidth)

#plt.tight_layout()
plt.savefig("sigI-NaCl." + format, bbox_inches="tight")

plt.ion()
#plt.show()
