#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser(description="Fit the no-average Cesaro of oneTwoDecompose corr "
                                             "calculated by calNDCesaro.py")
parser.add_argument('NDCesaroData', nargs='+', help="no-average Cesaro data to be averaged and fit. <data.NDCesaro.h5>")
parser.add_argument('-o', '--out', default='NDCesaro.fit.h5', help="output file")
parser.add_argument('-r', '--fitRange', nargs=2, type=float, metavar=('BEGIN', 'END'),
                    action='append', required=True, help="fitting range in ps. Multiple ranges are allowed")
args = parser.parse_args()

outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * (size - 1) + idx_c

numMD = len(args.NDCesaroData)

# sum the NDCesaroData
for n, data in enumerate(args.NDCesaroData):
  with h5py.File(data, 'r') as f:
    if (n == 0):
      numMol = f.attrs['numMol']
      numIonTypes = numMol.size
      numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
      charge = f.attrs['charge']
      zz = np.ones([numIonTypes + numIonTypePairs])
      for i in range(numIonTypes):
        zz[i] = charge[i] ** 2
        for j in range(i,numIonTypes):
          zz[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = charge[i] * charge[j]
      timeLags = f['timeLags'][:]
      NDCesaroRaw = np.zeros([numMD, numIonTypes + numIonTypePairs, timeLags.size])
      volumeRaw = np.zeros([numMD])

    NDCesaroRaw[n, :numIonTypes, :] += f['autoNDCesaro']
    NDCesaroRaw[n, numIonTypes:, :] += f['crossNDCesaro']
    volumeRaw[n] = f.attrs['cell'].prod()


NDCesaroTotalRaw = np.sum(NDCesaroRaw * zz[:, np.newaxis], axis=1)
NDCesaroTotal = {}
NDCesaroTotal['ave'] = np.mean(NDCesaroTotalRaw, axis=0)
NDCesaroTotal['std'] = np.std(NDCesaroTotalRaw, axis=0)
NDCesaroTotal['err'] = NDCesaroTotal['std'] / np.sqrt(numMD)

NDCesaro = {}
NDCesaro['ave'] = np.mean(NDCesaroRaw, axis=0)
NDCesaro['std'] = np.std(NDCesaroRaw, axis=0)
NDCesaro['err'] = NDCesaro['std'] / np.sqrt(numMD)

volume = {}
volume['ave'] = np.mean(volumeRaw)
volume['std'] = np.std(volumeRaw)
volume['err'] = volume['std'] / np.sqrt(numMD)

dt = timeLags[1] - timeLags[0]
fitRangeBoundary = (np.array(args.fitRange) / dt).astype(int)
fitRange = [list(range(i, j)) for [i, j] in fitRangeBoundary]

ND = {}
NDTotal = {}
for fit, fitBoundary in zip(fitRange, args.fitRange):
  ND[str(fitBoundary)] = np.polyfit(timeLags[fit], NDCesaro['ave'][:, fit].T, 1)[0, :]
  NDTotal[str(fitBoundary)] = np.polyfit(timeLags[fit], NDCesaroTotal['ave'][fit], 1)[0]

ND_err = {}
NDTotal_err = {}
if (numMD > 1):
  for fit, fitBoundary in zip(fitRange, args.fitRange):
    rec_sig2 = 1 / NDCesaro['std'][:, fit] ** 2  # [type, timeLags]
    S = np.sum(rec_sig2, 1) # [type]
    Sx = np.sum(timeLags[fit] * rec_sig2, 1) # [type]
    Sxx = np.sum(timeLags[fit]**2 * rec_sig2, 1) # [type]
    Sy = np.sum(NDCesaro['ave'][:, fit] * rec_sig2, 1) # [type]
    Syy = np.sum(NDCesaro['ave'][:, fit]**2 * rec_sig2, 1) # [type]
    Sxy = np.sum(timeLags[fit] * NDCesaro['ave'][:, fit] * rec_sig2, 1) # [type]

    delta = S * Sxx - Sx * Sx
    slope_b = (S * Sxy - Sx * Sy) / delta  # can be used for double check
    ND_err[str(fitBoundary)] = np.sqrt(S / delta)

  for fit, fitBoundary in zip(fitRange, args.fitRange):
    rec_sig2 = 1 / NDCesaroTotal['std'][fit] ** 2  # [timeLags]
    S = np.sum(rec_sig2, 1) # scalar
    Sx = np.sum(timeLags[fit] * rec_sig2, 1) # scalar
    Sxx = np.sum(timeLags[fit]**2 * rec_sig2, 1) # scalar
    Sy = np.sum(NDCesaroTotal['ave'][:, fit] * rec_sig2, 1) # scalar
    Syy = np.sum(NDCesaroTotal['ave'][:, fit]**2 * rec_sig2, 1) # scalar
    Sxy = np.sum(timeLags[fit] * NDCesaro['ave'][:, fit] * rec_sig2, 1) # scalar

    delta = S * Sxx - Sx * Sx
    slopeTotal_b = (S * Sxy - Sx * Sy) / delta  # can be used for double check
    NDTotal_err[str(fitBoundary)] = np.sqrt(S / delta)
else:
  for fitBoundary in args.fitRange:
    ND_err[str(fitBoundary)] = 0.
    NDTotal_err[str(fitBoundary)] = 0.

def saveDictToH5(h5g, name, dict):
  g = h5g.create_group(name)
  for k, v in dict.items():
    g[k] = v

with h5py.File(args.NDCesaroData[0], 'r') as f, h5py.File(args.out, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell' and name != 'timestep'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['zz'] = zz

  saveDictToH5(outFile, 'volume', volume)
  saveDictToH5(outFile, 'NDCesaro', NDCesaro)
  saveDictToH5(outFile, 'NDCesaroTotal', NDCesaroTotal)
  saveDictToH5(outFile, 'ND', ND)
  saveDictToH5(outFile, 'ND_err', ND_err)
  saveDictToH5(outFile, 'NDTotal', NDTotal)
  saveDictToH5(outFile, 'NDTotal_err', NDTotal_err)

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['NDCesaro/ave'].dims[1].attach_scale(outFile['timeLags'])
  outFile['NDCesaro/std'].dims[1].attach_scale(outFile['timeLags'])
  outFile['NDCesaro/err'].dims[1].attach_scale(outFile['timeLags'])
  outFile['NDCesaroTotal/ave'].dims[0].attach_scale(outFile['timeLags'])
  outFile['NDCesaroTotal/std'].dims[0].attach_scale(outFile['timeLags'])
  outFile['NDCesaroTotal/err'].dims[0].attach_scale(outFile['timeLags'])

print("File is output as: " + outFilename)


