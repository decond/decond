#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot two correlation results from oneTwoDecompose")
parser.add_argument('corrData', nargs=2, help="correlation data file <oneTwoDecompose.corr.h5>")
parser.add_argument('-o', '--out', default='corr', help="output figure base filename, default = 'corr'")
args = parser.parse_args()

fs = [None] * 2
with h5py.File(args.corrData[0], 'r') as fs[0], h5py.File(args.corrData[1], 'r') as fs[1]:
  timeLagsList = []
  autoCorrList = []
  crossCorrList = []

  for f in fs:
    timeLagsList.append(f['timeLags'][...])
    autoCorrList.append(f['autoCorr'][...])
    crossCorrList.append(f['crossCorr'][...])


for n, (timeLags, autoCorr, crossCorr) in enumerate(zip(timeLagsList, autoCorrList, crossCorrList)):
  for i, corr in enumerate(autoCorr[:]):
    if (n==0):
      plt.plot(timeLags, corr*100, label='{}-auto{}'.format(n, i))
    else:
      plt.plot(timeLags, corr*100, label='{}-auto{}'.format(n, i), linestyle='--')
  for i, corr in enumerate(crossCorr[:]):
    if (n==0):
      plt.plot(timeLags, corr*100, label='{}-cross{}'.format(n, i))
    else:
      plt.plot(timeLags, corr*100, label='{}-cross{}'.format(n, i), linestyle='--')
  plt.gca().set_color_cycle(None)

leg = plt.legend()
plt.xlabel(r'$\mathrm{time\ (ps)}$')
plt.ylabel(r'$C(t)\ (\AA^2 \mathrm{ps}^{-2}$)')
plt.ion()
plt.show()
