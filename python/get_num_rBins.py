import argparse
import math

parser = argparse.ArgumentParser(description="Get the minimum num_rBins among multile MD's")
parser.add_argument('dataFilename', nargs=1, help="sdDNoAverageCesaro filename")
parser.add_argument('rBinWidth', nargs=1, type=float, help="width of rBins")
parser.add_argument('numMD', nargs=1, type=int, help="number of MD to be examined")
args = parser.parse_args()

dataFilename = args.dataFilename[0]
rBinWidth = args.rBinWidth[0]
numMD = args.numMD[0]
num_rBins = [0] * numMD
for (idx, dataFilename) in enumerate(['./md' + str(i) + '/' + dataFilename for i in range(numMD)]):
  count = -1
  dataFile = open(dataFilename, 'r')
  for line in dataFile:
    if 'cell' in line:
      count = 3
    elif count > 0:
      count -= 1
    elif count == 0:
      num_rBins[idx] = math.ceil(float(line) / 2 / rBinWidth * math.sqrt(3))
      break
    else:
      continue

print(num_rBins)
print("Minimum num_rBins = " + str(min(num_rBins)))
