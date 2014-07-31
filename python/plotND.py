#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser(description="Plot and examine the results from fitNDCesaro.py")
parser.add_argument('NDCesaroFit', help="fitted ND results data file <NDCesaro.fit.h5>")
parser.add_argument('-o', '--out', default='NDCesaro.fit', help="output figure base filename")
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
    self.nd2ecSI = self.beta * self.basicCharge**2 / (volume*(self.nm**3)) * self.nm**2 / self.ps

def loadDictFromH5(h5g):
  dict = {}
  def func(k, v):
    dict[k] = v[...]
  h5g.visititems(func)
  return dict


with h5py.File(args.NDCesaroFit) as f:
  numMol = f.attrs['numMol'][...]
  numIonTypes = numMol.size
  numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
  charge = f.attrs['charge'][...]
  timeLags = f['timeLags'][...]
  zz = f['zz'][...]
  volume = loadDictFromH5(f['volume'])
  NDCesaro = loadDictFromH5(f['NDCesaro'])
  NDCesaroTotal = loadDictFromH5(f['NDCesaroTotal'])
  ND = loadDictFromH5(f['ND'])
  ND_err = loadDictFromH5(f['ND_err'])
  NDTotal = loadDictFromH5(f['NDTotal'])
  NDTotal_err = loadDictFromH5(f['NDTotal_err'])
  const = Const(volume['ave'])

