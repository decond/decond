import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

rc = {'legend': {'fontsize': 34},
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

labelPad = 10
spineLineWidth = 1

figsize = (11.5, 10)
format='pdf'

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
plt.ylabel(r"$\sigma_I(\lambda)$  (S m$^{-1}$)", labelpad=labelPad)
plt.xlabel(r"$r$  ($\AA$)", labelpad=labelPad)

halfCellIndex = rBins.size / np.sqrt(3)
plt.xlim(xmax=rBins[halfCellIndex])

for sp in plt.gca().spines.values():
  sp.set_linewidth(spineLineWidth)

plt.tight_layout()
plt.savefig("sigI-NaCl.pdf")

plt.ion()
#plt.show()
