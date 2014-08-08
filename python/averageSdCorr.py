#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from itertools import accumulate

parser = argparse.ArgumentParser(description="Average SD correlation")
parser.add_argument('sdcorrData', nargs='+', help="SD correlation data files to be averaged <sdcorr.h5>")
parser.add_argument('-o', '--out', help="output file, default = 'sdcorr.ave<num>.h5'")
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
# sum the NDCesaroData
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

if (isTimeLagsChanged):
  print("Note: the maximum timeLags are different among the corr files\n"
        "      it is now set to {} ps".format(timeLags[-1]))

sdCorr = np.mean(sdCorrN, axis=0)
sdCorr_std = np.std(sdCorrN, axis=0)
sdCorr_err = sdCorr_std / np.sqrt(numMD)
volume = np.mean(volumeN, axis=0)
volume_std = np.std(volumeN, axis=0)
volume_err = volume_std / np.sqrt(numMD)
rho = np.mean(rhoN, axis=0)
rho_std = np.std(rhoN, axis=0)
rho_err = rho_std / np.sqrt(numMD)

with h5py.File(args.sdcorrData[0], 'r') as f, h5py.File(outFilename, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['rBins'] = rBins
  outFile['volume'] = volume
  outFile['volume_err'] = volume_err
  outFile['sdCorr'] = sdCorr
  outFile['sdCorr_err'] = sdCorr_err
  outFile['rho'] = rho
  outFile['rho_err'] = rho_err

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['rBins'].dims.create_scale(outFile['rBins'], 'r')
  outFile['sdCorr'].dims[1].attach_scale(outFile['rBins'])
  outFile['sdCorr'].dims[2].attach_scale(outFile['timeLags'])
  outFile['sdCorr_err'].dims[1].attach_scale(outFile['rBins'])
  outFile['sdCorr_err'].dims[2].attach_scale(outFile['timeLags'])
  outFile['rho'].dims[1].attach_scale(outFile['rBins'])
  outFile['rho_err'].dims[1].attach_scale(outFile['rBins'])

print("File is output as: " + outFilename)
