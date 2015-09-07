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
parser.add_argument('-f', '--fit', nargs=2, type=float, metavar=('BEGIN', 'END'),
                    action='append', required=True,
                    help="fitting range in ps. Multiple ranges are allowed,"
                         "ex. -f <b1> <e1> -f <b2> <e2> ...")
args = parser.parse_args()

if args.add:
    outfile_mode = 'r+'  # Read/write, file must exist
else:
    outfile_mode = 'w-'  # Create file, fail if exists

da.new_decond(args.out, args.corrdata, args.fit)
print("output: " + args.out)
