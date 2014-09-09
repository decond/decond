#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from itertools import accumulate

parser = argparse.ArgumentParser(description="Average correlation")
parser.add_argument('corrData', nargs='+', help="correlation data files to be averaged <corr.h5>")
parser.add_argument('-o', '--out',
                    help="output file, default = 'corr.ave<num>.h5' or 'corr.ave<num>.w<window>.h5' "
                         "if <window> > 1")
parser.add_argument('-w', '--window', type=int, default=1,
                    help="combine every <window> rBins into one bin, default = 1")
parser.add_argument('-m', '--memoryFriendly', action='store_true',
                    help="turn on memory-friendly mode. Files are loaded and processed one by one."
                         "When combined with -e, files are loaded twice thus much slower.")
parser.add_argument('-e', '--error', action='store_true',
                    help="calculate and output the errors."
                         "Files are loaded twice if memory-friendly mode is on.")
parser.add_argument('--nosd', action='store_true', help="no-SD mode, i.e. one-two only mode")
args = parser.parse_args()

if (not args.nosd):
  try:
    for data in args.corrData:
      with h5py.File(data, 'r') as f:
        rBins = f['rBins'][...]
  except KeyError as e:
    print("Warning: no 'rBins' dataset is found in (at least)", data)
    print("Automatically change to --nosd mode")
    args.nosd = True

def zipIndexPair(idx_r, idx_c, size):
  """
  Returns the single index based the row index and column index
  """
  return idx_r * size + idx_c

def zipIndexPair2(idx_r, idx_c, size):
  """
  Returns the single index of a upper-half matrix based the row index and column index

  accepts only the "upper-half" index pair, because cross-correlation should
  be the same for (i,j) and (j,i)
  """
  assert(idx_r <= idx_c)
  return idx_r * size - ([0]+list(accumulate(range(4))))[idx_r] + idx_c - idx_r

numMD = len(args.corrData)
if (args.out == None):
  if (args.window == 1):
    if (args.nosd):
      outFilename = 'corr-nosd.ave' + str(numMD) + '.h5'
    else:
      outFilename = 'corr.ave' + str(numMD) + '.h5'
  else:
    if (args.nosd):
      outFilename = 'corr-nosd.ave' + str(numMD) + '.w' + str(args.window) + '.h5'
    else:
      outFilename = 'corr.ave' + str(numMD) + '.w' + str(args.window) + '.h5'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

isTimeLagsChanged = False

if (args.memoryFriendly):
  print("proceed with memory-friendly mode")
  print("determine the averages... ")
  for n, data in enumerate(args.corrData):
    print("reading " + data)
    with h5py.File(data, 'r') as f:
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
        charge = f.attrs['charge']
        timeLags = f['timeLags'][:]
        nCorr = np.zeros([numIonTypes*(numIonTypes+1), timeLags.size])
        volume = 0
        if (not args.nosd):
          rBinsTmp = f['rBins'][:]
          rBinsTmp = rBinsTmp[:rBinsTmp.size // args.window * args.window]
          rBins = np.concatenate(np.mean(np.split(rBinsTmp, rBinsTmp.size // args.window, axis=0),
                                         axis=1, keepdims=True), axis=0)
          sdCorr = np.zeros([numIonTypes**2, rBins.size, timeLags.size])
          rho = np.zeros([numIonTypes**2, rBins.size])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          nCorr = nCorr[..., :timeLags.size]
          if (not args.nosd):
            sdCorr = sdCorr[..., :timeLags.size]

      if (not args.nosd):
        if (f['rBins'].size < rBins.size*args.window):
          rBinsTmp = f['rBins'][:]
          rBinsTmp = rBinsTmp[:rBinsTmp.size // args.window * args.window]
          rBins = np.concatenate(np.mean(np.split(rBinsTmp, rBinsTmp.size // args.window, axis=0),
                                                  axis=1, keepdims=True), axis=0)
          sdCorr = sdCorr[..., :rBins.size, :]
          rho = rho[..., :rBins.size]

      nCorr += f['nCorr'][:, :timeLags.size]

      if ('cell' in f.attrs.keys()):
        volume += f.attrs['cell'].prod()
      else:
        volume += f['volume'][...]

      if (not args.nosd):
        sdCorrTmp = f['sdCorr'][:, :rBins.size*args.window, :timeLags.size] # [type, rBin*window, time]
        rhoTmp = f['rho'][:, :rBins.size*args.window] # [type, rBin*window]
        sdCorrTmp = sdCorrTmp * rhoTmp[:, :, np.newaxis] # [type, rBin*window, time]
        sdCorrTmp = np.array(np.split(sdCorrTmp, rBins.size, axis=1)) # [rBin, type, window, time]
        sdCorrTmp = np.swapaxes(np.sum(sdCorrTmp, axis=2), 0, 1) # [type, rBin, time]
        rhoTmp = np.array(np.split(rhoTmp, rBins.size, axis=1)) # [rBin, type, window]
        rhoTmp = np.sum(rhoTmp, axis=2).T # [type, rBin]

        sdCorr += sdCorrTmp / rhoTmp[:, :, np.newaxis]
        rho += rhoTmp

  nCorr /= numMD
  volume /= numMD
  if (not args.nosd):
    sdCorr /= numMD
    rho /= numMD

  if (args.error):
    print("determine the errors... ")
    if (numMD > 1):
      nCorr_std = np.zeros_like(nCorr)
      volume_std = 0.
      if (not args.nosd):
        sdCorr_std = np.zeros_like(sdCorr)
        rho_std = np.zeros_like(rho)

      for n, data in enumerate(args.corrData):
        print("reading " + data)
        with h5py.File(data, 'r') as f:
          nCorr_std += (f['nCorr'][:, :timeLags.size] - nCorr)**2
          if ('cell' in f.attrs.keys()):
            volume_std += (f.attrs['cell'].prod() - volume)**2
          else:
            volume_std += (f['volume'] - volume)**2

          if (not args.nosd):
            sdCorr_std += (f['sdCorr'][:, :rBins.size, :timeLags.size] - sdCorr)**2
            rho_std += (f['rho'][:, :rBins.size] - rho)**2

      nCorr_std = np.sqrt(nCorr_std / (numMD - 1)) 
      volume_std = np.sqrt(volume_std / (numMD - 1))
      nCorr_err = nCorr_std / np.sqrt(numMD)
      volume_err = volume_std / np.sqrt(numMD)

      if (not args.nosd):
        sdCorr_std = np.sqrt(sdCorr_std / (numMD - 1)) 
        rho_std = np.sqrt(rho_std / (numMD - 1))
        sdCorr_err = sdCorr_std / np.sqrt(numMD)
        rho_err = rho_std / np.sqrt(numMD)

    else:
      # numMD = 1
      nCorr_err = np.empty_like(nCorr)
      nCorr_err.fill(np.nan)
      volume_err = np.empty_like(volume)
      volume_err.fill(np.nan)
      if (not args.nosd):
        sdCorr_err = np.empty_like(sdCorr)
        sdCorr_err.fill(np.nan)
        rho_err = np.empty_like(rho)
        rho_err.fill(np.nan)

else:
  for n, data in enumerate(args.corrData):
    with h5py.File(data, 'r') as f:
      print("reading " + data)
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
        charge = f.attrs['charge']
        timeLags = f['timeLags'][:]
        nCorrN = np.zeros([numMD, numIonTypes*(numIonTypes+1), timeLags.size])
        volumeN = np.zeros([numMD])
        if (not args.nosd):
          rBinsTmp = f['rBins'][:]
          rBinsTmp = rBinsTmp[:rBinsTmp.size // args.window * args.window]
          rBins = np.concatenate(np.mean(np.split(rBinsTmp, rBinsTmp.size // args.window, axis=0),
                                                axis=1, keepdims=True), axis=0)
          sdCorrN = np.zeros([numMD, numIonTypes**2, rBins.size, timeLags.size])
          rhoN = np.empty([numMD, numIonTypes**2, rBins.size])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          nCorrN = nCorrN[..., :timeLags.size]
          if (not args.nosd):
            sdCorrN = sdCorrN[..., :timeLags.size]

      if (not args.nosd):
        if (f['rBins'].size < rBins.size*args.window):
          rBinsTmp = f['rBins'][:]
          rBinsTmp = rBinsTmp[:rBinsTmp.size // args.window * args.window]
          rBins = np.concatenate(np.mean(np.split(rBinsTmp, rBinsTmp.size // args.window, axis=0),
                                                  axis=1, keepdims=True), axis=0)
          sdCorrN = sdCorrN[..., :rBins.size, :]
          rhoN = rhoN[..., :rBins.size]

      nCorrN[n] = f['nCorr'][:, :timeLags.size]

      if ('cell' in f.attrs.keys()):
        volumeN[n] = f.attrs['cell'].prod()
      else:
        volumeN[n] = f['volume'][...]

      if (not args.nosd):
        sdCorrTmp = f['sdCorr'][:, :rBins.size*args.window, :timeLags.size] # [type, rBin*window, time]
        rhoTmp = f['rho'][:, :rBins.size*args.window] # [type, rBin*window]
        sdCorrTmp = sdCorrTmp * rhoTmp[:, :, np.newaxis] # [type, rBin*window, time]
        sdCorrTmp = np.array(np.split(sdCorrTmp, rBins.size, axis=1)) # [rBin, type, window, time]
        sdCorrTmp = np.swapaxes(np.sum(sdCorrTmp, axis=2), 0, 1) # [type, rBin, time]
        rhoTmp = np.array(np.split(rhoTmp, rBins.size, axis=1)) # [rBin, type, window]
        rhoTmp = np.sum(rhoTmp, axis=2).T # [type, rBin]

        sdCorrN[n] = sdCorrTmp / rhoTmp[:, :, np.newaxis]
        rhoN[n] = rhoTmp

  nCorr = np.mean(nCorrN, axis=0)
  volume = np.mean(volumeN, axis=0)
  if (not args.nosd):
    rho = np.mean(rhoN, axis=0)
    sdCorr = np.mean(sdCorrN, axis=0)

  if (args.error):
    if (numMD > 1):
      nCorr_std = np.std(nCorrN, axis=0)
      nCorr_err = nCorr_std / np.sqrt(numMD)
      volume_std = np.std(volumeN, axis=0)
      volume_err = volume_std / np.sqrt(numMD)
      if (not args.nosd):
        sdCorr_std = np.std(sdCorrN, axis=0)
        sdCorr_err = sdCorr_std / np.sqrt(numMD)
        rho_std = np.std(rhoN, axis=0)
        rho_err = rho_std / np.sqrt(numMD)
    else:
      # numMD = 1
      nCorr_err = np.empty_like(nCorr)
      nCorr_err.fill(np.nan)
      volume_err = np.empty_like(volume)
      volume_err.fill(np.nan)
      if (not args.nosd):
        sdCorr_err = np.empty_like(sdCorr)
        sdCorr_err.fill(np.nan)
        rho_err = np.empty_like(rho)
        rho_err.fill(np.nan)

if (isTimeLagsChanged):
  print("Note: the maximum timeLags are different among the corr files\n"
        "      it is now set to {} ps".format(timeLags[-1]))

with h5py.File(args.corrData[0], 'r') as f, h5py.File(outFilename, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (numMD == 1 or (numMD > 1 and name != 'cell')):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['volume'] = volume
  outFile['nCorr'] = nCorr
  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['nCorr'].dims[1].attach_scale(outFile['timeLags'])
  if (args.error):
    outFile['volume_err'] = volume_err
    outFile['nCorr_err'] = nCorr_err
    outFile['nCorr_err'].dims[1].attach_scale(outFile['timeLags'])

  if (not args.nosd):
    outFile['rBins'] = rBins
    outFile['sdCorr'] = sdCorr
    outFile['rho'] = rho
    outFile['rBins'].dims.create_scale(outFile['rBins'], 'r')
    outFile['sdCorr'].dims[1].attach_scale(outFile['rBins'])
    outFile['sdCorr'].dims[2].attach_scale(outFile['timeLags'])
    outFile['rho'].dims[1].attach_scale(outFile['rBins'])
    if (args.error):
      outFile['sdCorr_err'] = sdCorr_err
      outFile['rho_err'] = rho_err
      outFile['sdCorr_err'].dims[1].attach_scale(outFile['rBins'])
      outFile['sdCorr_err'].dims[2].attach_scale(outFile['timeLags'])
      outFile['rho_err'].dims[1].attach_scale(outFile['rBins'])

print("File is output as: " + outFilename)
