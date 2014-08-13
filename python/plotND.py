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

colorCode = ['blue', 'green', 'blue', 'red', 'green']
labelMol = [r'[C$_4$mim]', r'[NTf$_2$]']
labels = labelMol + [m1 + '-' + m2 for i1, m1 in enumerate(labelMol)
                                   for i2, m2 in enumerate(labelMol)
                                   if i1 <= i2 ]
class Const:
  """
  Defines some constants
  """
  kB = 1.3806488E-23 #(J K-1)
  beta = 1 / (kB * args.temp) #(J-1)
  basicCharge = 1.60217646E-19 #(Coulomb)
  ps = 1.0E-12 #(s)
  nm = 1.0E-9 #(m)
  nm2cm = 1.0E-7
  ps2s = 1.0E-12

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
  volume_err = f['volume_err'][...]
  NDCesaro = f['NDCesaro'][...]
  NDCesaro_err = f['NDCesaro_err'][...]
  NDCesaroTotal = f['NDCesaroTotal'][...]
  NDCesaroTotal_err = f['NDCesaroTotal_err'][...]
  ND = loadDictFromH5(f['ND'])
  ND_err = loadDictFromH5(f['ND_err'])
  NDTotal = loadDictFromH5(f['NDTotal'])
  NDTotal_err = loadDictFromH5(f['NDTotal_err'])
  const = Const(volume)

DI = {}
DI_err = {}
for fit in ND:
  DI[fit] = ND[fit][:numIonTypes] / numMol * 1E5 * Const.nm2cm**2 / Const.ps2s # 10^-5 cm^2 / s
  DI_err[fit] = ND_err[fit][:numIonTypes] / numMol * 1E5 * Const.nm2cm**2 / Const.ps2s # 10^-5 cm^2 / s
ec = {}
ec_err = {}
ecTotal = {}
ecTotal_err = {}
sortedKeys = sorted(ND.keys(), key=lambda x:x.split(sep='-')[0])
for k in sortedKeys:
  print(k + ':')
  print("========================")
  print("Electrical conductivity in S / m:")
  ec[k] = ND[k] * const.ND2ecSI * zz
  ec_err[k] = ND_err[k] * const.ND2ecSI
  print(ec[k])
  print("+/-\n", ec_err[k], sep="")
  ecTotal[k] = const.ND2ecSI * sum(ND[k]*zz*ww)
  ecTotal_err[k] = const.ND2ecSI * sum(abs(ND_err[k]))
  print("Total: ", ecTotal[k], " +/- ", ecTotal_err[k], '\n', sep="")
  print("Diffusion constant in 10^-5 cm^2 / s:")
  print(DI[k])
  print("+/-\n", DI_err[k], '\n', sep="")

plt.figure(1)
for i, (color, label) in enumerate(zip(colorCode, labels)):
  if (i < numIonTypes):
    plt.plot(range(len(ecTotal)), [ec[k][i] for k in sortedKeys], linestyle='--', color=color, label=label)
  else:
    plt.plot(range(len(ecTotal)), [ec[k][i] for k in sortedKeys], color=color, label=label)
plt.plot(range(len(ecTotal)), [ecTotal[k] for k in sortedKeys], color='black', label='total')
plt.xticks(range(len(ecTotal)), sortedKeys)
plt.legend()
plt.xlabel("fit range  (ps)")
plt.ylabel(r"$\sigma$  (S m$^{-1}$)")

print("Total")
for k in sorted(NDTotal.keys(), key=lambda x:x.split(sep='-')[0]):
  print(k + ':')
  print(NDTotal[k] * const.ND2ecSI, " +/- ", NDTotal_err[k] * const.ND2ecSI, '\n', sep="")

numErrBars = 5

plt.figure(2)
for i, (NDC, NDC_err, color, label)  in enumerate(zip(NDCesaro, NDCesaro_err, colorCode, labels)):
  if (i < numIonTypes):
    plt.errorbar(timeLags, NDC, yerr=NDC_err, errorevery=timeLags.size//numErrBars,
                 linestyle='--', label=label, color=color)
  else:
    plt.errorbar(timeLags, NDC, yerr=NDC_err, errorevery=timeLags.size//numErrBars,
                 label=label, color=color)

plt.legend(loc='upper left')
plt.xlabel("time lag  (ps)")
plt.ylabel(r"$\tilde D_I$ $\tilde D_{IL}$  ($\AA^2$)")

plt.ion()
plt.show()
