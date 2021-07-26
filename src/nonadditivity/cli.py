# -*- coding: utf-8 -*-

"""Command line interface for :mod:`nonadditivity`.

Why does this file exist, and why not put this in ``__main__``? You might be tempted to import things from ``__main__``
later, but that will cause problems--the code will get executed twice:

- When you run ``python3 -m nonadditivity`` python will execute``__main__.py`` as a script. That means there won't be any
  ``pybel.__main__`` in ``sys.modules``.
- When you import __main__ it will get executed again (as a module) because
  there's no ``nonadditivity.__main__`` in ``sys.modules``.

.. seealso:: https://click.palletsprojects.com/en/latest/setuptools/#setuptools-integration
"""

import argparse

from .api import run_nonadd_calculation


def main():
    parser = argparse.ArgumentParser(
        description="Nonadditivity Analysis based on matched MMPs",
        epilog="\nINPUT FORMAT:\t Identifier [sep] SMILES [sep] data[...]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-in", required=True, dest="infile",
                        help="Input file")
    parser.add_argument("-update", action="store_true",
                        help="Use fragmentation and images from previous run")
    parser.add_argument("-max_heavy", type=int, default=70,
                        help="Maximum number of heavy atoms per ligand")
    parser.add_argument("-no_chiral", action="store_true",
                        help="Skip all transformations that include chiral centers")
    parser.add_argument("-delimiter", choices=["tab", "space", "comma", "semicolon"], default="tab",
                        help="Specify delimiter in input file")
    parser.add_argument("-include_censored", action="store_true",
                        help="Include Circles where one out of four compounds has a censored value")
    parser.add_argument("-out", default=None, dest="outfile",
                        help="Output file name")
    parser.add_argument("-series_column", default=None,
                        help="Column that identifies subseries within the dataset")
    parser.add_argument("-props", nargs='*', default=None,
                        help="Property columns for which Nonadditivity should be calculated")
    parser.add_argument("-units", nargs='*', default=None,
                        help="Unit of the activity given. Need to supply either none or as many as props. \
                             Units need to be one out of M, mM, uM, nM, pM, noconv.")
    parser.add_argument("-shorts", nargs='*', default=None,
                        help="Activity Abbreviations. Need to supply either none or as many as props.")

    # Legacy functions
    # parser.add_argument("-write_images", action="store_true",
    #                    help="Write Nonadditivity circle images to ./images/ directory")

    args = parser.parse_args()
    args.write_images = False

    if args.shorts and args.props:
        if not len(args.shorts) == len(args.props):
            print("Please either supply no or exactly as many shortcut names as property columns")
            print("Exiting.")
            exit(0)

    if args.units and args.props:
        if not len(args.units) == len(args.props):
            print("Please either supply no or exactly as many units as property columns")
            print("Exiting.")
            exit(0)

    run_nonadd_calculation(args)


if __name__ == '__main__':
    main()
