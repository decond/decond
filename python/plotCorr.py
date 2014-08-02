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

  p = []

  for i, corr in enumerate(autoCorr[:]):
    p.append(plt.plot(timeLags, corr*100, label='auto{}'.format(i)))

  for i, corr in enumerate(crossCorr[:]):
    p.append(plt.plot(timeLags, corr*100, label='cross{}'.format(i)))

  p = [i for j in p for i in j]

  leg = plt.legend()
  plt.xlabel(r'$\mathrm{time\ (ps)}$')
  plt.ylabel(r'$C(t)\ (\AA^2 \mathrm{ps}^{-2}$)')
  plt.ion()
  plt.show()
