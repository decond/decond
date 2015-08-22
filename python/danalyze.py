import argparse
import decond.analyzer as da

DEFAULT_OUTFILENAME = 'decond.d5'

parser = argparse.ArgumentParser(description="Analyze corr data from Decond")
parser.add_argument('corrdata', nargs='+',
                    help="correlation data file. <corr.c5>")
parser.add_argument('-o', '--out', default=DEFAULT_OUTFILENAME,
                    help="output decond file, default <{0}>".format(
                        DEFAULT_OUTFILENAME))
parser.add_argument('-a', '--add', action='store_true',
                    help="incrementally add data samples mode")
args = parser.parse_args()

if args.add:
    outfile_mode = 'r+'  # Read/write, file must exist
else:
    outfile_mode = 'w-'  # Create file, fail if exists

da.cal_decond(args.out, args.corrdata)
print("output: " + args.out)
