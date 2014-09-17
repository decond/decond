#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from scipy import integrate
from itertools import accumulate

parser = argparse.ArgumentParser(description="Fit the no-average Cesaro of corr calculated by calCesaro.py")
parser.add_argument('cesaroData', nargs='+',
                    help="no-average Cesaro data to be averaged and fit. <cesaro.h5>")
parser.add_argument('-o', '--out', help="output file, default = 'cesaro.fit.h5'")
parser.add_argument('-r', '--fitRange', nargs=2, type=float, metavar=('BEGIN', 'END'),
                    action='append', required=True,
                    help="fitting range in ps. Multiple ranges are allowed,"
                         "ex. -r <b1> <e1> -r <b2> <e2> ...")
parser.add_argument('-m', '--memoryFriendly', action='store_true',
                    help="turn on memory-friendly mode. Files are loaded and processed one by one,"
                         "thus slower but more memory friendly")
parser.add_argument('--nosd', action='store_true', help="no-SD mode, i.e. one-two only mode")
args = parser.parse_args()

if (not args.nosd):
  print("checking if all cesaroData contain sdData ...")
  try:
    for data in args.cesaroData:
      with h5py.File(data, 'r') as f:
        rBins = f['rBins'][...]
  except KeyError as e:
    print("Warning: no 'rBins' dataset is found in (at least)", data)
    print("Automatically change to --nosd mode")
    args.nosd = True

if (args.out is None):
  if (args.nosd):
    outFilename = 'cesaro-nosd.fit.h5'
  else:
    outFilename = 'cesaro.fit.h5'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(accumulate(range(4))))[idx_r] + idx_c - idx_r

numMD = len(args.cesaroData)

isTimeLagsChanged = False

if (args.memoryFriendly):
  print("proceed with memory-friendly mode")
  print("determine the averages... ")
  for n, data in enumerate(args.cesaroData):
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
        zz = np.ones([numIonTypes + numIonTypePairs])
        ww = np.ones([numIonTypes + numIonTypePairs])
        # ww: weight of each component. ex. for NaCl, ww = [1, 1, 1, 2, 1]
        for i in range(numIonTypes):
          zz[i] = charge[i] ** 2
          for j in range(i,numIonTypes):
            zz[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = charge[i] * charge[j]
            if (i == j):
              ww[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = 1
            else:
              ww[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = 2
        timeLags = f['timeLags'][...]
        nDCesaro = np.zeros([numIonTypes + numIonTypePairs, timeLags.size])
        volume = 0.
        if (not args.nosd):
          rBins = f['rBins'][...]
          sdDCesaro = np.zeros([numIonTypePairs, rBins.size, timeLags.size])
          rho2 = np.zeros([numIonTypePairs, rBins.size])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          nDCesaro = nDCesaro[:, :timeLags.size]
          if (not args.nosd):
            sdDCesaro = sdDCesaro[:, :, :timeLags.size]

      if (not args.nosd):
        if (f['rBins'].size < rBins.size):
          rBins = f['rBins'][...]
          sdDCesaro = sdDCesaro[:, :rBins.size, :]
          rho2 = rho2[:, :rBins.size]

      nDCesaro += f['nDCesaro'][:, :timeLags.size]
      if ('cell' in f.attrs.keys()):
        volume += f.attrs['cell'].prod()
      else:
        volume += f['volume'][...]

      if (not args.nosd):
        sdDCesaro += f['sdDCesaro'][:, :rBins.size, :timeLags.size]
        rho2 += f['rho2'][:, :rBins.size]

  nDCesaro /= numMD
  nDCesaroTotal = np.sum(nDCesaro * (zz * ww)[:, np.newaxis], axis=0)
  volume /= numMD
  if (not args.nosd):
    sdDCesaro /= numMD
    rho2 /= numMD

  if (numMD > 1):
    print("determine the errors... ")
    nDCesaro_std = np.zeros_like(nDCesaro)
    nDCesaroTotal_std = np.zeros_like(nDCesaroTotal)
    volume_std = 0.
    if (not args.nosd):
      sdDCesaro_std = np.zeros_like(sdDCesaro)
      rho2_std = np.zeros_like(rho2)

    for n, data in enumerate(args.cesaroData):
      print("reading " + data)
      with h5py.File(data, 'r') as f:
        tmp_nDCesaro = f['nDCesaro'][:, :timeLags.size]
        nDCesaro_std += (tmp_nDCesaro - nDCesaro)**2
        nDCesaroTotal_std += (np.sum(tmp_nDCesaro * (zz * ww)[:, np.newaxis], axis=0) - nDCesaroTotal)**2
        if ('cell' in f.attrs.keys()):
          volume_std += (f.attrs['cell'].prod() - volume)**2
        else:
          volume_std += (f['volume'] - volume)**2

        if (not args.nosd):
          sdDCesaro_std += (f['sdDCesaro'][:, :rBins.size, :timeLags.size] - sdDCesaro)**2
          rho2_std += (f['rho2'][:, :rBins.size] - rho2)**2

    nDCesaro_std = np.sqrt(nDCesaro_std / (numMD - 1)) 
    nDCesaroTotal_std = np.sqrt(nDCesaroTotal_std / (numMD - 1)) 
    nDCesaro_err = nDCesaro_std / np.sqrt(numMD)
    nDCesaroTotal_err = nDCesaroTotal_std / np.sqrt(numMD)
    volume_std = np.sqrt(volume_std / (numMD - 1))
    volume_err = volume_std / np.sqrt(numMD)

    if (not args.nosd):
      sdDCesaro_std = np.sqrt(sdDCesaro_std / (numMD - 1)) 
      rho2_std = np.sqrt(rho2_std / (numMD - 1))
      sdDCesaro_err = sdDCesaro_std / np.sqrt(numMD)
      rho2_err = rho2_std / np.sqrt(numMD)

  else:
    # numMD = 1
    nDCesaro_err = np.empty_like(nDCesaro)
    nDCesaro_err.fill(np.nan)
    nDCesaroTotal_err = np.empty_like(nDCesaroTotal)
    nDCesaroTotal_err.fill(np.nan)
    volume_err = np.empty_like(volume)
    volume_err.fill(np.nan)
    if (not args.nosd):
      sdDCesaro_err = np.empty_like(sdDCesaro)
      sdDCesaro_err.fill(np.nan)
      rho2_err = np.empty_like(rho2)
      rho2_err.fill(np.nan)

else:
  for n, data in enumerate(args.cesaroData):
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
        zz = np.ones([numIonTypes + numIonTypePairs])
        ww = np.ones([numIonTypes + numIonTypePairs])
        # ww: weight of each term.
        # ex. for a 2-component molecule such as NaCl, ww = [1, 1, 1, 2, 1]
        for i in range(numIonTypes):
          zz[i] = charge[i] ** 2
          for j in range(i,numIonTypes):
            zz[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = charge[i] * charge[j]
            if (i == j):
              ww[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = 1
            else:
              ww[numIonTypes + zipIndexPair2(i, j, numIonTypes)] = 2
        timeLags = f['timeLags'][...]
        nDCesaroRaw = np.empty([numMD, numIonTypes + numIonTypePairs, timeLags.size])
        volumeRaw = np.empty([numMD])
        if (not args.nosd):
          rBins = f['rBins'][...]
          sdDCesaroRaw = np.empty([numMD, numIonTypePairs, rBins.size, timeLags.size])
          rho2Raw = np.empty([numMD, numIonTypePairs, rBins.size])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          nDCesaroRaw = nDCesaroRaw[:, :, :timeLags.size]
          if (not args.nosd):
            sdDCesaroRaw = sdDCesaroRaw[:, :, :, :timeLags.size]

      if (not args.nosd):
        if (f['rBins'].size < rBins.size):
          rBins = f['rBins'][...]
          sdDCesaroRaw = sdDCesaroRaw[:, :, :rBins.size, :]
          rho2Raw = rho2Raw[:, :, :rBins.size]

      nDCesaroRaw[n] = f['nDCesaro'][:, :timeLags.size]
      if ('cell' in f.attrs.keys()):
        volumeRaw[n] = f.attrs['cell'].prod()
      else:
        volumeRaw[n] = f['volume'][...]

      if (not args.nosd):
        sdDCesaroRaw[n] = f['sdDCesaro'][:, :rBins.size, :timeLags.size]
        rho2Raw[n] = f['rho2'][:, :rBins.size]

  nDCesaroTotalRaw = np.sum(nDCesaroRaw * (zz * ww)[:, np.newaxis], axis=1)

  nDCesaroTotal = np.mean(nDCesaroTotalRaw, axis=0)
  nDCesaro = np.mean(nDCesaroRaw, axis=0)
  volume = np.mean(volumeRaw)
  if (not args.nosd):
    sdDCesaro = np.mean(sdDCesaroRaw, axis=0)
    rho2 = np.mean(rho2Raw, axis=0)

  if (numMD > 1):
    nDCesaro_std = np.std(nDCesaroRaw, axis=0)
    nDCesaro_err = nDCesaro_std / np.sqrt(numMD)
    nDCesaroTotal_std = np.std(nDCesaroTotalRaw, axis=0)
    nDCesaroTotal_err = nDCesaroTotal_std / np.sqrt(numMD)
    volume_std = np.std(volumeRaw)
    volume_err = volume_std / np.sqrt(numMD)
    if (not args.nosd):
      sdDCesaro_std = np.std(sdDCesaroRaw, axis=0)
      sdDCesaro_err = sdDCesaro_std / np.sqrt(numMD)
      rho2_std = np.std(rho2Raw, axis=0)
      rho2_err = rho2_std / np.sqrt(numMD)
  else:
    nDCesaro_err = np.empty_like(nDCesaro)
    nDCesaro_err.fill(np.nan)
    nDCesaroTotal_err = np.empty_like(nDCesaroTotal)
    nDCesaroTotal_err.fill(np.nan)
    volume_err = np.empty_like(volume)
    volume_err.fill(np.nan)
    if (not args.nosd):
      sdDCesaro_err = np.empty_like(sdDCesaro)
      sdDCesaro_err.fill(np.nan)
      rho2_err = np.empty_like(rho2)
      rho2_err.fill(np.nan)

if (isTimeLagsChanged):
  print("Note: the maximum timeLags are different among the sdcorr files\n"
        "      it is now set to {} ps".format(timeLags[-1]))

dt = timeLags[1] - timeLags[0]
fitRangeBoundary = (np.array(args.fitRange) / dt).astype(int)
fitRange = [list(range(i, j)) for [i, j] in fitRangeBoundary]

def getKeyFromFitBoundary(fitBoundary):
  return '{}-{}'.format(*fitBoundary)

print("fitting data...")
nD = {}
nDTotal = {}
for fit, fitBoundary in zip(fitRange, args.fitRange):
  nD[getKeyFromFitBoundary(fitBoundary)] = np.polyfit(timeLags[fit], nDCesaro[:, fit].T, 1)[0, :]
  nDTotal[getKeyFromFitBoundary(fitBoundary)] = np.polyfit(timeLags[fit], nDCesaroTotal[fit], 1)[0]

if (not args.nosd):
  sdD = {}
  for fit, fitBoundary in zip(fitRange, args.fitRange):
    sdD[getKeyFromFitBoundary(fitBoundary)] = np.empty([numIonTypePairs, rBins.size])
    for tp in range(numIonTypePairs):
      # I don't know why the array should not be transposed?
      # sdD[getKeyFromFitBoundary(fitBoundary)][tp, ...] = np.polyfit(timeLags[fit], sdDCesaro[tp, ..., fit].T, 1)[0, :]
      sdD[getKeyFromFitBoundary(fitBoundary)][tp, ...] = np.polyfit(timeLags[fit], sdDCesaro[tp, ..., fit], 1)[0, :]

def fitErr(timeLags, cesaro, cesaro_std, fit):
  rec_sig2 = 1 / cesaro_std[..., fit] ** 2  # [type, rBins, timeLags] or [type, timeLags]
  S = np.sum(rec_sig2, -1) # [type, rBins] or [type]
  Sx = np.sum(timeLags[fit] * rec_sig2, -1) # [type, rBins] or [type]
  Sxx = np.sum(timeLags[fit]**2 * rec_sig2, -1) # [type, rBins] or [type]
  Sy = np.sum(cesaro[..., fit] * rec_sig2, -1) # [type, rBins] or [type]
  Syy = np.sum(cesaro[..., fit]**2 * rec_sig2, -1) # [type, rBins] or [type]
  Sxy = np.sum(timeLags[fit] * cesaro[..., fit] * rec_sig2, -1) # [type, rBins] or [type]

  delta = S * Sxx - Sx * Sx
  #slope_b = (S * Sxy - Sx * Sy) / delta  # can be used for double check
  return np.sqrt(S / delta)

nD_err = {}
nDTotal_err = {}
if (not args.nosd):
  sdD_err = {}
if (numMD > 1):
  for fit, fitBoundary in zip(fitRange, args.fitRange):
    nD_err[getKeyFromFitBoundary(fitBoundary)] = fitErr(timeLags, nDCesaro, nDCesaro_std, fit)
    nDTotal_err[getKeyFromFitBoundary(fitBoundary)] = fitErr(timeLags, nDCesaroTotal, nDCesaroTotal_std, fit)
    if (not args.nosd):
      sdD_err[getKeyFromFitBoundary(fitBoundary)] = fitErr(timeLags, sdDCesaro, sdDCesaro_std, fit)

else:
  for fitBoundary in args.fitRange:
    nD_err[getKeyFromFitBoundary(fitBoundary)] = np.empty_like(nD[getKeyFromFitBoundary(fitBoundary)])
    nD_err[getKeyFromFitBoundary(fitBoundary)][...] = np.nan
    nDTotal_err[getKeyFromFitBoundary(fitBoundary)] = np.empty_like(nDTotal[getKeyFromFitBoundary(fitBoundary)])
    nDTotal_err[getKeyFromFitBoundary(fitBoundary)][...] = np.nan
    if (not args.nosd):
      sdD_err[getKeyFromFitBoundary(fitBoundary)] = np.zeros_like([sdD[getKeyFromFitBoundary(fitBoundary)]])
      sdD_err[getKeyFromFitBoundary(fitBoundary)][...] = np.nan

def saveDictToH5(h5g, name, dict):
  g = h5g.create_group(name)
  for k, v in dict.items():
    g[k] = v

with h5py.File(args.cesaroData[0], 'r') as f, h5py.File(outFilename, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell' and name != 'timestep'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['zzCross'] = zzCross
  outFile['zz'] = zz
  outFile['ww'] = ww
  outFile['volume'] = volume
  outFile['volume_err'] = volume_err
  outFile['nDCesaro'] = nDCesaro
  outFile['nDCesaro_err'] = nDCesaro_err
  outFile['nDCesaroTotal'] = nDCesaroTotal
  outFile['nDCesaroTotal_err'] = nDCesaroTotal_err

  saveDictToH5(outFile, 'nD', nD)
  saveDictToH5(outFile, 'nD_err', nD_err)
  saveDictToH5(outFile, 'nDTotal', nDTotal)
  saveDictToH5(outFile, 'nDTotal_err', nDTotal_err)

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')

  outFile['nDCesaro'].dims[1].attach_scale(outFile['timeLags'])
  outFile['nDCesaro_err'].dims[1].attach_scale(outFile['timeLags'])
  outFile['nDCesaroTotal'].dims[0].attach_scale(outFile['timeLags'])
  outFile['nDCesaroTotal_err'].dims[0].attach_scale(outFile['timeLags'])

  if (not args.nosd):
    outFile['rBins'] = rBins
    outFile['sdDCesaro'] = sdDCesaro
    outFile['sdDCesaro_err'] = sdDCesaro_err
    outFile['rho2'] = rho2
    outFile['rho2_err'] = rho2_err

    saveDictToH5(outFile, 'sdD', sdD)
    saveDictToH5(outFile, 'sdD_err', sdD_err)

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
