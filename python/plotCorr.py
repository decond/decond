#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot correlation results from oneTwoDecompose")
parser.add_argument('corrData', help="correlation data file <oneTwoDecompose.corr.h5>")
parser.add_argument('-o', '--out', default='corr', help="output figure base filename, default = 'corr'")
args = parser.parse_args()

with h5py.File(args.corrData, 'r') as f:
  timeLags = f['timeLags'][...]
  autoCorr = f['autoCorr'][...]
  crossCorr = f['crossCorr'][...]

  for i, corr in enumerate(autoCorr[:]):
    plt.plot(timeLags, corr*100, label='auto-{}'.format(i), linestyle='--')

  plt.gca().set_color_cycle(None)
  for i, corr in enumerate(crossCorr[:]):
    plt.plot(timeLags, corr*100, label='cross-{}'.format(i))

  leg = plt.legend()
  plt.xlabel(r'time  (ps)')
  plt.ylabel(r'$C_I^{(1)}(t)$, $C_{IL}^{(2)}(t)$  ($\AA^2$ ps$^{-2}$)')
  plt.ion()
  plt.show()
