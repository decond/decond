#!/usr/bin/env python3
import argparse
import decond.analyze as da

DEFAULT_OUTFILENAME = 'decond.d5'


def new(args):
    da.new_decond(args.out, args.corr, args.fit)
    print("output: " + args.out)


def add(args):
    da.extend_decond(args.out, args.decond, args.corr, args.fit)
    print("output: " + args.out)


def fit(args):
    da.fit_decond(args.out, args.decond, args.fit)
    print("output: " + args.out)


# create the top-level parser
parser = argparse.ArgumentParser(
        description="Decond analysis tool, use subcommands to perform tasks")
subparsers = parser.add_subparsers(
        description="dec subcommand -h for more specific help. "
                    "Note that all the subcommands will create new output "
                    "instead of overwriting existing data")


# create the parser for the "new" subcommand
parser_new = subparsers.add_parser(
        'new',
        help="new decond analysis from corr.c5 data")

parser_new.add_argument('corr', nargs='+',
                        help="correlation data file. <corr.c5>")
parser_new.add_argument('-f', '--fit', nargs=2, type=float,
                        metavar=('BEGIN', 'END'),
                        action='append', required=True,
                        help="fitting range in ps. Multiple ranges are allowed"
                             ", ex. -f <b1> <e1> -f <b2> <e2> ...")
parser_new.add_argument('-o', '--out', default=DEFAULT_OUTFILENAME,
                        help="output decond file, default <{0}>".format(
                            DEFAULT_OUTFILENAME))

parser_new.set_defaults(func=new)


# create the parser for the "add" subcommand
parser_add = subparsers.add_parser(
        'add',
        help="add corr.c5 data to existing decond.d5")

parser_add.add_argument('decond',
                        help="decond analysis file. <decond.d5>")
parser_add.add_argument('corr', nargs='+',
                        help="correlation data file. <corr.c5>")
parser_add.add_argument('-f', '--fit', nargs=2, type=float,
                        metavar=('BEGIN', 'END'),
                        action='append',
                        help="fitting range in ps. Multiple ranges are allowed"
                             ", ex. -f <b1> <e1> -f <b2> <e2> ...")
parser_add.add_argument('-o', '--out', default=DEFAULT_OUTFILENAME,
                        help="output decond file, default <{0}>".format(
                            DEFAULT_OUTFILENAME))

parser_add.set_defaults(func=add)


# create the parser for the "fit" subcommand
parser_add = subparsers.add_parser(
        'fit',
        help="change fit ranges of existing decond.d5")

parser_add.add_argument('decond',
                        help="decond analysis file. <decond.d5>")
parser_add.add_argument('-f', '--fit', nargs=2, type=float,
                        metavar=('BEGIN', 'END'),
                        action='append', required=True,
                        help="fitting range in ps. Multiple ranges are allowed"
                             ", ex. -f <b1> <e1> -f <b2> <e2> ...")
parser_add.add_argument('-o', '--out', default=DEFAULT_OUTFILENAME,
                        help="output decond file, default <{0}>".format(
                            DEFAULT_OUTFILENAME))

parser_add.set_defaults(func=fit)


# parse the args and call whatever function was selected
args = parser.parse_args()
args.func(args)
