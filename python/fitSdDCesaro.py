#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from scipy import integrate
from itertools import accumulate

parser = argparse.ArgumentParser(description="Fit the no-average Cesaro of sdcorr calculated by calSdDCesaro.py")
parser.add_argument('sdDCesaroData', nargs='+',
                    help="no-average Cesaro data to be averaged and fit. <data.sdDCesaro.h5>")
parser.add_argument('-o', '--out', default='sdDCesaro.fit.h5',
                    help="output file, default = 'sdDCesaro.fit.h5'")
parser.add_argument('-r', '--fitRange', nargs=2, type=float, metavar=('BEGIN', 'END'),
                    action='append', required=True, help="fitting range in ps. Multiple ranges are allowed")
parser.add_argument('-m', '--memoryFriendly', action='store_true',
                    help="turn on memory-friendly mode. Files are loaded and processed one by one,"
                         "thus slower but more memory friendly")
args = parser.parse_args()

outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(accumulate(range(4))))[idx_r] + idx_c - idx_r

numMD = len(args.sdDCesaroData)

isTimeLagsChanged = False

if (args.memoryFriendly):
  print("proceed with memory-friendly mode")
  print("determine the averages... ")
  for n, data in enumerate(args.sdDCesaroData):
    print("reading " + data)
    with h5py.File(data, 'r') as f:
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2;
        charge = f.attrs['charge']
        zzCross = np.ones([numIonTypePairs])
        for i in range(numIonTypes):
          for j in range(i,numIonTypes):
            zzCross[zipIndexPair2(i, j, numIonTypes)] = charge[i] * charge[j]
        timeLags = f['timeLags'][...]
        rBins = f['rBins'][...]
        sdDCesaro = np.zeros([numIonTypePairs, rBins.size, timeLags.size])
        rho2 = np.zeros([numIonTypePairs, rBins.size])
        volume = 0.

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          sdDCesaro = sdDCesaro[:, :, :timeLags.size]

      if (f['rBins'].size < rBins.size):
        rBins = f['rBins'][...]
        sdDCesaro = sdDCesaro[:, :rBins.size, :]
        rho2 = rho2[:, :rBins.size]

      sdDCesaro += f['sdDCesaro'][:, :rBins.size, :timeLags.size]
      rho2 += f['rho2'][:, :rBins.size]
      volume += f.attrs['cell'].prod()

  sdDCesaro /= numMD
  rho2 /= numMD
  volume /= numMD

  print("determine the errors... ")
  sdDCesaro_std = np.zeros_like(sdDCesaro)
  rho2_std = np.zeros_like(rho2)
  volume_std = 0.

  for n, data in enumerate(args.sdDCesaroData):
    print("reading " + data)
    with h5py.File(data, 'r') as f:
      sdDCesaro_std += (f['sdDCesaro'][:, :rBins.size, :timeLags.size] - sdDCesaro)**2
      rho2_std += (f['rho2'][:, :rBins.size] - rho2)**2
      volume_std += (f.attrs['cell'].prod() - volume)**2

  sdDCesaro_std = np.sqrt(sdDCesaro_std / (numMD - 1)) 
  rho2_std = np.sqrt(rho2_std / (numMD - 1))
  volume_std = np.sqrt(volume_std / (numMD - 1))

  sdDCesaro_err = sdDCesaro_std / np.sqrt(numMD)
  rho2_err = rho2_std / np.sqrt(numMD)
  volume_err = volume_std / np.sqrt(numMD)

else:
  for n, data in enumerate(args.sdDCesaroData):
    print("reading " + data)
    with h5py.File(data, 'r') as f:
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) // 2;
        charge = f.attrs['charge']
        zzCross = np.ones([numIonTypePairs])
        for i in range(numIonTypes):
          for j in range(i,numIonTypes):
            zzCross[zipIndexPair2(i, j, numIonTypes)] = charge[i] * charge[j]
        timeLags = f['timeLags'][...]
        rBins = f['rBins'][...]
        sdDCesaroRaw = np.empty([numMD, numIonTypePairs, rBins.size, timeLags.size])
        rho2Raw = np.empty([numMD, numIonTypePairs, rBins.size])
        volumeRaw = np.empty([numMD])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          sdDCesaroRaw = sdDCesaroRaw[:, :, :, :timeLags.size]

      if (f['rBins'].size < rBins.size):
        rBins = f['rBins'][...]
        sdDCesaroRaw = sdDCesaroRaw[:, :, :rBins.size, :]
        rho2Raw = rho2Raw[:, :, :rBins.size]

      sdDCesaroRaw[n] = f['sdDCesaro'][:, :rBins.size, :timeLags.size]
      rho2Raw[n] = f['rho2'][:, :rBins.size]
      volumeRaw[n] = f.attrs['cell'].prod()

  if (isTimeLagsChanged):
    print("Note: the maximum timeLags are different among the sdcorr files\n"
          "      it is now set to {} ps".format(timeLags[-1]))

  sdDCesaro = np.mean(sdDCesaroRaw, axis=0)
  sdDCesaro_std = np.std(sdDCesaroRaw, axis=0)
  sdDCesaro_err = sdDCesaro_std / np.sqrt(numMD)

  rho2 = np.mean(rho2Raw, axis=0)
  rho2_std = np.std(rho2Raw, axis=0)
  rho2_err = rho2_std / np.sqrt(numMD)

  volume = np.mean(volumeRaw)
  volume_std = np.std(volumeRaw)
  volume_err = volume_std / np.sqrt(numMD)

dt = timeLags[1] - timeLags[0]
fitRangeBoundary = (np.array(args.fitRange) / dt).astype(int)
fitRange = [list(range(i, j)) for [i, j] in fitRangeBoundary]

def getKeyFromFitBoundary(fitBoundary):
  return '{}-{}'.format(*fitBoundary)

sdD = {}
for fit, fitBoundary in zip(fitRange, args.fitRange):
  sdD[getKeyFromFitBoundary(fitBoundary)] = np.empty([numIonTypePairs, rBins.size])
  for tp in range(numIonTypePairs):
    # I don't know why the array should not be transposed?
    # sdD[getKeyFromFitBoundary(fitBoundary)][tp, ...] = np.polyfit(timeLags[fit], sdDCesaro[tp, ..., fit].T, 1)[0, :]
    sdD[getKeyFromFitBoundary(fitBoundary)][tp, ...] = np.polyfit(timeLags[fit], sdDCesaro[tp, ..., fit], 1)[0, :]

sdD_err = {}
if (numMD > 1):
  for fit, fitBoundary in zip(fitRange, args.fitRange):
    rec_sig2 = 1 / sdDCesaro_std[..., fit] ** 2  # [type, rBins, timeLags]
    S = np.sum(rec_sig2, 2) # [type, rBins]
    Sx = np.sum(timeLags[fit] * rec_sig2, 2) # [type, rBins]
    Sxx = np.sum(timeLags[fit]**2 * rec_sig2, 2) # [type, rBins]
    Sy = np.sum(sdDCesaro[..., fit] * rec_sig2, 2) # [type, rBins]
    Syy = np.sum(sdDCesaro[..., fit]**2 * rec_sig2, 2) # [type, rBins]
    Sxy = np.sum(timeLags[fit] * sdDCesaro[..., fit] * rec_sig2, 2) # [type, rBins]

    delta = S * Sxx - Sx * Sx
    slope_b = (S * Sxy - Sx * Sy) / delta  # can be used for double check
    sdD_err[getKeyFromFitBoundary(fitBoundary)] = np.sqrt(S / delta)

else:
  for fitBoundary in args.fitRange:
    sdD_err[getKeyFromFitBoundary(fitBoundary)] = np.zeros_like([sdD[getKeyFromFitBoundary(fitBoundary)]])
    sdD_err[getKeyFromFitBoundary(fitBoundary)][...] = np.nan

def saveDictToH5(h5g, name, dict):
  g = h5g.create_group(name)
  for k, v in dict.items():
    g[k] = v

with h5py.File(args.sdDCesaroData[0], 'r') as f, h5py.File(args.out, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell' and name != 'timestep'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['rBins'] = rBins
  outFile['zzCross'] = zzCross
  outFile['volume'] = volume
  outFile['volume_err'] = volume_err
  outFile['sdDCesaro'] = sdDCesaro
  outFile['sdDCesaro_err'] = sdDCesaro_err
  outFile['rho2'] = rho2
  outFile['rho2_err'] = rho2_err

  saveDictToH5(outFile, 'sdD', sdD)
  saveDictToH5(outFile, 'sdD_err', sdD_err)

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['rBins'].dims.create_scale(outFile['rBins'], 'r')
  outFile['sdDCesaro'].dims[1].attach_scale(outFile['rBins'])
  outFile['sdDCesaro'].dims[2].attach_scale(outFile['timeLags'])
  outFile['sdDCesaro_err'].dims[1].attach_scale(outFile['rBins'])
  outFile['sdDCesaro_err'].dims[2].attach_scale(outFile['timeLags'])
  outFile['rho2'].dims[1].attach_scale(outFile['rBins'])
  outFile['rho2_err'].dims[1].attach_scale(outFile['rBins'])

  for dset in outFile['sdD'].values():
    dset.dims[1].attach_scale(outFile['rBins'])

  for dset in outFile['sdD_err'].values():
    dset.dims[1].attach_scale(outFile['rBins'])

print("File is output as: " + outFilename)
