#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser(description="Fit the no-average Cesaro of oneTwoDecompose corr "
                                             "calculated by calNDCesaro.py")
parser.add_argument('NDCesaroData', nargs='+', help="No-average Cesaro data to be averaged and fit. <data.NDCesaro.h5>")
parser.add_argument('--dt', required=True, type=int, help="data inverval (fs) in NDCesaroData (in integer). "
                                                            "dt = timestep * integration delta step")
parser.add_argument('-o', '--out' , default='NDCesaroData.fit.h5', help="output file")
args = parser.parse_args()

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
    volumeRaw[n] += f.attrs['cell'].prod()


NDCesaroTotal = np.sum(NDCesaroRaw * zz[:, np.newaxis], axis=1)
NDCesaro = {}
NDCesaro['ave'] = np.mean(NDCesaroRaw, axis=0)
NDCesaro['std'] = np.std(NDCesaroRaw, axis=0)
NDCesaro['err'] = NDCesaro['std'] / np.sqrt(numMD)

volume = {}
volume['ave'] = np.mean(volumeRaw)
volume['std'] = np.std(volumeRaw)
volume['err'] = volume['std'] / np.sqrt(numMD)

def saveDictToH5(h5g, name, dict):
  g = h5g.create_group(name)
  for k, v in dict.items():
    g[k] = v

with h5py.File(args.NDCesaroData[0], 'r') as f, h5py.File(args.out, 'w') as outFile:
  for (name, value) in f.attrs.items():
    outFile.attrs[name] = value

  outFile['timeLags'] = timeLags
#  saveDictToH5(outFile, 'ND', ND)
  saveDictToH5(outFile, 'volume', volume)
  saveDictToH5(outFile, 'NDCesaro', NDCesaro)

print("File is output as: " + args.out)
