#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot and examine the results from fitNDCesaro.py")
parser.add_argument('NDCesaroFit', help="fitted ND results data file <NDCesaro.fit.h5>")
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


with h5py.File(args.NDCesaroFit, 'r') as f:
  numMol = f.attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
  charge = f.attrs['charge'][...]
  timeLags = f['timeLags'][...]
  zz = f['zz'][...]
  ww = f['ww'][...]
  volume = f['volume'][...]
  NDCesaro = f['NDCesaro'][...]
  NDCesaro_err = f['NDCesaro_err'][...]
  NDCesaroTotal = f['NDCesaroTotal'][...]
  NDCesaroTotal_err = f['NDCesaroTotal_err'][...]
  ND = loadDictFromH5(f['ND'])
  ND_err = loadDictFromH5(f['ND_err'])
  NDTotal = loadDictFromH5(f['NDTotal'])
  NDTotal_err = loadDictFromH5(f['NDTotal_err'])
  const = Const(volume)

for k in sorted(ND.keys(), key=lambda x:x.split(sep='-')[0]):
  print(k + ':')
  print(ND[k] * const.ND2ecSI)
  print("Total:", const.ND2ecSI * sum(ND[k]*zz*ww))

print("Total")
for k in sorted(NDTotal.keys(), key=lambda x:x.split(sep='-')[0]):
  print(k + ':')
  print(NDTotal[k] * const.ND2ecSI, '\n')

for i, (NDC, NDC_err)  in enumerate(zip(NDCesaro, NDCesaro_err)):
  plt.errorbar(timeLags, NDC, yerr=NDC_err, errorevery=10000, label='{}'.format(i))
plt.legend()

plt.ion()
plt.show()
