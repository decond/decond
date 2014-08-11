#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from itertools import accumulate

parser = argparse.ArgumentParser(description="Average SD correlation")
parser.add_argument('sdcorrData', nargs='+', help="SD correlation data files to be averaged <sdcorr.h5>")
parser.add_argument('-o', '--out', help="output file, default = 'sdcorr.ave<num>.h5'")
parser.add_argument('-m', '--memoryFriendly', action='store_true',
                    help="turn on memory-friendly mode. Files are loaded and processed one by one."
                         "When combined with -e, files are loaded twice thus much slower.")
parser.add_argument('-e', '--error', action='store_true',
                    help="calculate and output the errors."
                         "Files are loaded twice if memory-friendly mode is on.")
args = parser.parse_args()

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

numMD = len(args.sdcorrData)
if (args.out == None):
  outFilename = 'sdcorr.ave' + str(numMD) + '.h5'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

isTimeLagsChanged = False

if (args.memoryFriendly):
  print("proceed with memory-friendly mode")
  print("determine the averages... ")
  for n, data in enumerate(args.sdcorrData):
    print("reading " + data)
    with h5py.File(data, 'r') as f:
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
        charge = f.attrs['charge']
        rBins = f['rBins'][...]
        timeLags = f['timeLags'][:]
        sdCorr = np.zeros([numIonTypes**2, rBins.size, timeLags.size])
        rho = np.empty([numIonTypes**2, rBins.size])
        volume = 0

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          sdCorr = sdCorr[..., :timeLags.size]

      if (f['rBins'].size < rBins.size):
        rBins = f['rBins'][...]
        sdCorr = sdCorr[..., :rBins.size, :]
        rho = rho[..., :rBins.size]

      sdCorr += f['sdCorr'][:, :rBins.size, :timeLags.size]
      rho += f['rho'][:, :rBins.size]
      volume += f.attrs['cell'].prod()

  sdCorr /= numMD
  rho /= numMD
  volume /= numMD

  if (args.error):
    print("determine the errors... ")
    sdCorr_std = np.zeros_like(sdCorr)
    rho_std = np.zeros_like(rho)
    volume_std = 0.

    for n, data in enumerate(args.sdcorrData):
      print("reading " + data)
      with h5py.File(data, 'r') as f:
        sdCorr_std += (f['sdCorr'][:, :rBins.size, :timeLags.size] - sdCorr)**2
        rho_std += (f['rho'][:, :rBins.size] - rho)**2
        volume_std += (f.attrs['cell'].prod() - volume)**2

    sdCorr_std = np.sqrt(sdCorr_std / (numMD - 1)) 
    rho_std = np.sqrt(rho_std / (numMD - 1))
    volume_std = np.sqrt(volume_std / (numMD - 1))

    sdCorr_err = sdCorr_std / np.sqrt(numMD)
    rho_err = rho_std / np.sqrt(numMD)
    volume_err = volume_std / np.sqrt(numMD)

else:
  for n, data in enumerate(args.sdcorrData):
    with h5py.File(data, 'r') as f:
      print("reading " + data)
      if (n == 0):
        numMol = f.attrs['numMol']
        numIonTypes = numMol.size
        numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
        charge = f.attrs['charge']
        rBins = f['rBins'][...]
        timeLags = f['timeLags'][:]
        sdCorrN = np.zeros([numMD, numIonTypes**2, rBins.size, timeLags.size])
        rhoN = np.empty([numMD, numIonTypes**2, rBins.size])
        volumeN = np.zeros([numMD])

      if (f['timeLags'].size != timeLags.size):
        isTimeLagsChanged = True
        if (f['timeLags'].size < timeLags.size):
          timeLags = f[timeLags][...]
          sdCorrN = sdCorrN[..., :timeLags.size]

      if (f['rBins'].size < rBins.size):
        rBins = f['rBins'][...]
        sdCorrN = sdCorrN[..., :rBins.size, :]
        rhoN = rhoN[..., :rBins.size]

      sdCorrN[n] = f['sdCorr'][:, :rBins.size, :timeLags.size]
      rhoN[n] = f['rho'][:, :rBins.size]
      volumeN[n] = f.attrs['cell'].prod()

  sdCorr = np.mean(sdCorrN, axis=0)
  volume = np.mean(volumeN, axis=0)
  rho = np.mean(rhoN, axis=0)

  if (args.error):
    sdCorr_std = np.std(sdCorrN, axis=0)
    sdCorr_err = sdCorr_std / np.sqrt(numMD)
    volume_std = np.std(volumeN, axis=0)
    volume_err = volume_std / np.sqrt(numMD)
    rho_std = np.std(rhoN, axis=0)
    rho_err = rho_std / np.sqrt(numMD)

if (isTimeLagsChanged):
  print("Note: the maximum timeLags are different among the corr files\n"
        "      it is now set to {} ps".format(timeLags[-1]))

with h5py.File(args.sdcorrData[0], 'r') as f, h5py.File(outFilename, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['rBins'] = rBins
  outFile['volume'] = volume
  outFile['sdCorr'] = sdCorr
  outFile['rho'] = rho
  if (args.error):
    outFile['volume_err'] = volume_err
    outFile['sdCorr_err'] = sdCorr_err
    outFile['rho_err'] = rho_err

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['rBins'].dims.create_scale(outFile['rBins'], 'r')
  outFile['sdCorr'].dims[1].attach_scale(outFile['rBins'])
  outFile['sdCorr'].dims[2].attach_scale(outFile['timeLags'])
  outFile['rho'].dims[1].attach_scale(outFile['rBins'])
  if (args.error):
    outFile['sdCorr_err'].dims[1].attach_scale(outFile['rBins'])
    outFile['sdCorr_err'].dims[2].attach_scale(outFile['timeLags'])
    outFile['rho_err'].dims[1].attach_scale(outFile['rBins'])

print("File is output as: " + outFilename)
