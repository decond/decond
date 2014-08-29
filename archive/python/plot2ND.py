#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description="Plot and examine two results from fitNDCesaro.py")
parser.add_argument('NDCesaroFit', nargs=2, help="fitted ND results data file <NDCesaro.fit.h5>")
parser.add_argument('-o', '--out', default='NDCesaro.fit', help="output figure base filename, default = 'NDCesaro.fit'")
parser.add_argument('-T', '--temp', type=float, required=True, help="temperature in K")
args = parser.parse_args()

class Const:
  """
  Defines some constants
  """
  kB = 1.3806488E-23 #(J K-1)
  beta = 1 / (kB * args.temp) #(J-1)
  basicCharge = 1.60217646E-19 #(Coulomb)
  ps = 1.0E-12 #(s)
  nm = 1.0E-9 #(m)

  def __init__(self, volume):
    self.ND2ecSI = self.beta * self.basicCharge**2 / (volume*(self.nm**3)) * self.nm**2 / self.ps

def loadDictFromH5(h5g):
  dict = {}
  def func(k, v):
    dict[k] = v[...]
  h5g.visititems(func)
  return dict

f = [None] * 2
with h5py.File(args.NDCesaroFit[0], 'r') as f[0], h5py.File(args.NDCesaroFit[1], 'r') as f[1]:
  numMol = f[0].attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
  charge = f[0].attrs['charge'][...]
  zz = f[0]['zz'][...]
  ww = f[0]['ww'][...]
  volume = f[0]['volume'][...]

  NDCesaroList = []
  NDCesaro_errList = []
  NDCesaroTotalList = []
  NDCesaroTotal_errList = []
  NDList = []
  ND_errList = []
  NDTotalList = []
  NDTotal_errList = []
  timeLagsList = []
  for file in (f):
    timeLagsList.append(file['timeLags'][...])
    NDCesaroList.append(file['NDCesaro'][...])
    NDCesaro_errList.append(file['NDCesaro_err'][...])
    NDCesaroTotalList.append(file['NDCesaroTotal'][...])
    NDCesaroTotal_errList = file['NDCesaroTotal_err'][...]
    NDList = loadDictFromH5(file['ND'])
    ND_errList = loadDictFromH5(file['ND_err'])
    NDTotalList = loadDictFromH5(file['NDTotal'])
    NDTotal_errList = loadDictFromH5(file['NDTotal_err'])

for n, (timeLags, NDCesaro, NDCesaro_err) in enumerate(zip(timeLagsList, NDCesaroList, NDCesaro_errList)):
  for i, (NDC, NDC_err) in enumerate(zip(NDCesaro, NDCesaro_err)):
    plt.errorbar(timeLags, NDC, yerr=NDC_err, errorevery=10000, label='{}-{}'.format(n,i))

plt.legend()
plt.xlabel("time lag  (ps)")
plt.ylabel(r"$\tilde D_I$ $\tilde D_{IL}$  ($\AA^2$)")

plt.ion()
plt.show()
