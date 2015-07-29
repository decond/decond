#!/usr/bin/env python
import sys
import argparse
import h5py
import numpy as np

VERSION = "0.1.0"

parser = argparse.ArgumentParser(description="Combine the energy trajectories output from ERmod")
parser.add_argument('engtrajFirst', metavar='engtraj', nargs=1)
parser.add_argument('engtrajRest', metavar='engtraj', nargs='+')
parser.add_argument('-o', '--out', default='engtraj-all.h5', help="output file, default = 'engtraj-all.h5'")
parser.add_argument('--version', '-v' , action='version', version=VERSION)
args = parser.parse_args()

outFile = args.out
engtrajFiles = args.engtrajFirst + args.engtrajRest

def getMajor(version):
  return int(version.split(sep='.')[0])

# Check versions
for file in engtrajFiles:
  with h5py.File(file, 'r') as f:
    engtrajVersion = f.attrs['version'].decode('utf-8')
    engtrajMajor = getMajor(engtrajVersion)
    major = getMajor(VERSION)
    if (engtrajMajor < major):
      print("Error: version not compatible!")
      print("This program is version ", VERSION)
      print(file + " is version " + engtrajVersion)
      sys.exit(1)

sltspecList = []
sltfirsttagList = []
numsltList = []
nummolList = []
numpairList = []
numframeList = []

class DataInconsistentError(Exception):
  pass

# === Check data consistency ===
for file in engtrajFiles:
  with h5py.File(file, 'r') as f:
    sltspecList.append(f.attrs['sltspec'])
    sltfirsttagList.append(f.attrs['sltfirsttag'])
    numsltList.append(f.attrs['numslt'])
    nummolList.append(f.attrs['nummol'])
    numpairList.append(f.attrs['numpair'])
    numframeList.append(f.attrs['numframe'])
  
# sort the files and other attribute list according to sltspec
# ref: http://stackoverflow.com/questions/9543211/sorting-a-list-in-python-using-the-result-from-sorting-another-list
sortLists = (sltspecList, engtrajFiles, sltfirsttagList, numsltList, nummolList, numpairList, numframeList)
sortLists = zip(*sorted(zip(*sortLists)))

if not sltspecList == list(range(1, len(sltspecList) + 1)):
  print(sltspecList)
  print(list(range(1, len(sltspecList) + 1)))
  raise DataInconsistentError("sltspec's are not continous integers starting from 1")

for i in range(len(sltfirsttagList) - 1):
  if not sltfirsttagList[i] < sltfirsttagList[i+1]:
    raise DataInconsistentError("The order of sltfirsttag is not consistent with that of sltspec")

if not all(n == nummolList[0] for n in nummolList):
  raise DataInconsistentError("nummol's are not all the same")

if not all(n == numframeList[0] for n in numframeList):
  raise DataInconsistentError("numframe's are not all the same")

if not sum(numsltList) == nummolList[0]:
  raise DataInconsistentError("sum of numslt's do not equal to nummol")

if not sum(numpairList) == nummolList[0] * (nummolList[0] - 1) / 2:
  raise DataInconsistentError("numpair is not consistent with nummol")
# ===================================

numpairTotal = sum(numpairList)
fout = h5py.File(outFile, 'x')
fout.attrs['version'] = VERSION
fout.attrs['sltspec'] = sltspecList
fout.attrs['sltfirsttag'] = sltfirsttagList
fout.attrs['numslt'] = numsltList
fout.attrs['nummol'] = nummolList[0]
fout.attrs['numpair'] = numpairTotal
fout.attrs['numframe'] = numframeList[0]

energyTotal = fout.create_dataset("energy", (numpairTotal, numframeList[0]))

start = 0
for i, file in enumerate(engtrajFiles):
  with h5py.File(file, 'r') as f:
    end = start + numpairList[i]
    energyTotal[start:end, :] = f['energy'][:, :] 
    start = end

print("File is output as: " + outFile)
