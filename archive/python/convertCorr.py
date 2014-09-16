#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from itertools import accumulate

parser = argparse.ArgumentParser(
  description="Convert corr files from separated HDF5 format to newer combined HDF5 format "
              "or from old separated Octave format to newer HDF5 format")
parser.add_argument('corrData', nargs=2, metavar=('oneTwoData', 'sdData'), help="a pair of corr files to be converted")
parser.add_argument('-o', '--out', help="output file, default = 'corr-converted.h5'")
parser.add_argument('--octave', action='store_true', help="turn on Octave mode")
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

if (args.out == None):
  outFilename = 'corr-converted.h5'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

if (args.octave):
  import scipy.io as sio
  fcorr = sio.loadmat(args.corrData[0])
  fsd = sio.loadmat(args.corrData[1])
  with h5py.File(outFilename, 'w') as outFile:
    outFile.attrs['timestep'] = np.squeeze(fcorr['timestep'][...], axis=1)
    outFile.attrs['cell'] = np.squeeze(fcorr['cell'][...])
    outFile.attrs['charge'] = np.squeeze(fcorr['charge'][...]).astype(np.int32)
    outFile.attrs['numMol'] = np.squeeze(fcorr['numAtom'][...].astype(np.int32))

    outFile['timeLags'] = np.squeeze(fcorr['timeLags'][...])
    outFile['nCorr'] = np.concatenate([fcorr['autoCorr'][...].T, fcorr['crossCorr'][...].T])

    outFile['rBins'] = np.squeeze(fsd['rBins'][...])
    outFile['rho'] = fsd['rho'][...].T
    outFile['sdCorr'] = fsd['sdCorr'][...].T

else:
  with h5py.File(args.corrData[0], 'r') as fcorr,\
       h5py.File(args.corrData[1], 'r') as fsd,\
       h5py.File(outFilename, 'w') as outFile:
    for (name, value) in fcorr.attrs.items():
      outFile.attrs[name] = value

    outFile['timeLags'] = fcorr['timeLags'][...]
    outFile['rBins'] = fsd['rBins'][...]
    outFile['nCorr'] = np.concatenate([fcorr['autoCorr'][...], fcorr['crossCorr'][...]])

    if ('autoCorr_err' in fcorr):
      outFile['nCorr_err'] = np.concatenate([fcorr['autoCorr_err'][...], fcorr['crossCorr_err'][...]])

    if ('volume' in fcorr):
      outFile['volume'] = fcorr['volume'][...]
      if ('volume_err' in fcorr):
        outFile['volume_err'] = fcorr['volume_err'][...]

    if ('sdCorr_err' in fsd):
      outFile['sdCorr_err'] = fsd['sdCorr_err'][...]

    outFile['sdCorr'] = fsd['sdCorr'][...]
    outFile['rho'] = fsd['rho'][...]

print("File is output as: " + outFilename)
