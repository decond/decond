#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser(description="Average oneTwoDecompose correlation")
parser.add_argument('corrData', nargs='+', help="correlation data files to be averaged <oneTwoDecompose.corr.h5>")
parser.add_argument('-o', '--out', default='corr.ave.h5', help="output file")
args = parser.parse_args()

outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

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
  return idx_r * (size - 1) + idx_c

numMD = len(args.corrData)

# sum the NDCesaroData
for n, data in enumerate(args.corrData):
  with h5py.File(data, 'r') as f:
    if (n == 0):
      numMol = f.attrs['numMol']
      numIonTypes = numMol.size
      numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
      charge = f.attrs['charge']
      timeLags = f['timeLags'][:]
      autoCorrN = np.zeros([numMD, numIonTypes, timeLags.size])
      crossCorrN = np.zeros([numMD, numIonTypePairs, timeLags.size])
      volumeN = np.zeros([numMD])

    volumeN[n] = f.attrs['cell'].prod()
    autoCorrN[n, :, :] += f['autoCorr']
    for i in range(numIonTypes):
      for j in range(i, numIonTypes):
        if (i == j):
          crossCorrN[n, zipIndexPair2(i, j, numIonTypes), :] += \
              f['crossCorr'][zipIndexPair(i, j, numIonTypes), :]
        else:
          crossCorrN[n, zipIndexPair2(i, j, numIonTypes), :] += \
              (f['crossCorr'][zipIndexPair(i, j, numIonTypes), :] +
              f['crossCorr'][zipIndexPair(j, i, numIonTypes), :]) / 2

autoCorr = np.mean(autoCorrN, axis=0)
crossCorr = np.mean(crossCorrN, axis=0)
volume = np.mean(volumeN, axis=0)
autoCorr_std = np.std(autoCorrN, axis=0)
crossCorr_std = np.std(crossCorrN, axis=0)
volume_std = np.std(volumeN, axis=0)
autoCorr_err = autoCorr_std / np.sqrt(numMD)
crossCorr_err = crossCorr_std / np.sqrt(numMD)
volume_err = volume_std / np.sqrt(numMD)

with h5py.File(args.corrData[0], 'r') as f, h5py.File(args.out, 'w') as outFile:
  for (name, value) in f.attrs.items():
    if (name != 'cell'):
      outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
  outFile['volume'] = volume
  outFile['volume_err'] = volume_err
  outFile['autoCorr'] = autoCorr
  outFile['crossCorr'] = crossCorr
  outFile['autoCorr_err'] = autoCorr_err
  outFile['crossCorr_err'] = crossCorr_err

  outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
  outFile['autoCorr'].dims[1].attach_scale(outFile['timeLags'])
  outFile['autoCorr_err'].dims[1].attach_scale(outFile['timeLags'])
  outFile['crossCorr'].dims[1].attach_scale(outFile['timeLags'])
  outFile['crossCorr_err'].dims[1].attach_scale(outFile['timeLags'])

print("File is output as: " + outFilename)
