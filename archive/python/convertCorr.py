#!/home/kmtu/local/anaconda3/bin/python
import argparse
import h5py
import numpy as np
from itertools import accumulate

parser = argparse.ArgumentParser(description="Convert corr files from old separated format to newer format")
parser.add_argument('corrData', nargs=2,
                    help="pair of corr files to be converted <corr.h5> and <sdcorr.h5>")
parser.add_argument('-o', '--out', help="output file, default = 'allcorr.h5'")
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
  outFilename = 'allcorr.h5'
else:
  outFilename = args.out if args.out.split('.')[-1] == 'h5' else args.out + '.h5'

with h5py.File(args.corrData[0], 'r') as fcorr,\
     h5py.File(args.corrData[1], 'r') as fsd,\
     h5py.File(outFilename, 'w') as outFile:
  for (name, value) in fcorr.attrs.items():
    outFile.attrs[name] = value

  outFile['timeLags'] = fcorr['timeLags'][...]
  outFile['rBins'] = fsd['rBins'][...]
  outFile['nCorr'] = np.concatenate([fcorr['autoCorr'][...], fcorr['crossCorr'][...]])
  outFile['sdCorr'] = fsd['sdCorr'][...]
  outFile['rho'] = fsd['rho'][...]

print("File is output as: " + outFilename)
