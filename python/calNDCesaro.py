#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from scipy import integrate

parser = argparse.ArgumentParser(description="Calcuate no-average Cesaro ND of oneTwoDecompose corr")
parser.add_argument('corrData', help="oneTwoDecomposed correlation data file. <data.corr.h5>")
parser.add_argument('--intDelta', type=int, default=1, help="integration delta step. Default = 1")
parser.add_argument('-o', '--out', default='NDCesaro.h5', help="output file, default = 'NDCesaro.h5'")
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

with h5py.File(args.corrData, 'r') as f:
  timestep = f.attrs['timestep'][0]
  numMol = f.attrs['numMol']

  timeLags = f['timeLags'][0::args.intDelta]
  autoCorr = f['autoCorr'][:, 0::args.intDelta]
  crossCorr = f['crossCorr'][:, 0::args.intDelta]
  assert(autoCorr.shape[-1] == crossCorr.shape[-1])

  print("timestep = {}".format(timestep))
  print("numMol = {}".format(numMol))

  autoNDCesaro = integrate.cumtrapz(autoCorr, timeLags, initial=0) # nm**2 / ps
  autoNDCesaro = integrate.cumtrapz(autoNDCesaro, timeLags, initial=0) # nm**2

  numIonTypes = numMol.size
  crossNDCesaro = np.empty([numIonTypes * (numIonTypes + 1) / 2, timeLags.size])
  for i in range(numIonTypes):
    for j in range(i, numIonTypes):
      if (i == j):
        crossNDCesaro[zipIndexPair2(i, j, numIonTypes), :] = crossCorr[zipIndexPair(i, j, numIonTypes), :]
      else:
        crossNDCesaro[zipIndexPair2(i, j, numIonTypes), :] = (crossCorr[zipIndexPair(i, j, numIonTypes), :] +
                                                              crossCorr[zipIndexPair(j, i, numIonTypes), :]) / 2
  crossNDCesaro = integrate.cumtrapz(crossNDCesaro, timeLags, initial = 0)
  crossNDCesaro = integrate.cumtrapz(crossNDCesaro, timeLags, initial = 0)

  with h5py.File(outFilename, 'w') as outFile:
    for (name, value) in f.attrs.items():
      outFile.attrs[name] = value

    outFile['timeLags'] = timeLags
    outFile['autoNDCesaro'] = autoNDCesaro
    outFile['crossNDCesaro'] = crossNDCesaro

    outFile['timeLags'].dims.create_scale(outFile['timeLags'], 't')
    outFile['autoNDCesaro'].dims[1].attach_scale(outFile['timeLags'])
    outFile['crossNDCesaro'].dims[1].attach_scale(outFile['timeLags'])

  print("File is output as: " + outFilename)
