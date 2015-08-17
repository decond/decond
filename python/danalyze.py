import argparse
import decond.analyzer as da

DEFAULT_OUTFILENAME = 'decond.h5'

parser = argparse.ArgumentParser(description="Analyze corr data from Decond")
parser.add_argument('corrData', nargs='+',
                    help="correlation data file. <corr.h5>")
parser.add_argument('-o', '--out', default=DEFAULT_OUTFILENAME,
                    help="output file, default = 'decond.h5'")
parser.add_argument('-a', '--add', action='store_true',
                    help="incrementally add data samples mode")
args = parser.parse_args()

if args.add:
    outFileMode = 'r+'  # Read/write, file must exist
else:
    outFileMode = 'w-'  # Create file, fail if exists

with da.DecondFile(args.out, outFileMode) as f:
    f.addSample(args.corrData)
