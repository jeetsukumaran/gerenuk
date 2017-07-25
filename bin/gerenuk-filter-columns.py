#! /usr/bin/env python

import os
import sys
import argparse
import traceback
import time
import tempfile
from gerenuk import utility

def main():
    parser = argparse.ArgumentParser(
            description="GERENUK Simultaneous Divergence Time Analysis -- Post-Process Columns",
            )
    parser.add_argument(
            "target_data_filepath",
            nargs="+",
            help="Path to target or observed data file.")
    filter_options = parser.add_argument_group("Filter Options")
    filter_options.add_argument(
            "--master-column-filepath",
            help="If specified, then only columns with names"
                 " in this file will be retained in the target"
                 " file(s).")
    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument('--field-delimiter',
        type=str,
        default='\t',
        help="Field delimiter (default: <TAB>').")
    run_options.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Work silently.")
    args = parser.parse_args()
    if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
        open = utility.pre_py34_open
    if args.master_column_filepath is not None:
        if not args.quiet:
            sys.stderr.write("-gerenuk- Retaining only columns defined in '{}'\n".format(args.master_column_filepath))
        with open(args.master_column_filepath) as master_src:
            columns_to_retain = utility.extract_fieldnames_from_file(
                    src=master_src,
                    field_delimiter=args.field_delimiter)
        for filepath in args.target_data_filepath:
            with open(filepath) as src:
                with open(filepath + ".filtered", "w") as dest:
                # with tempfile.NamedTemporaryFile() as dest:
                    utility.filter_columns_from_file(
                            src=src,
                            dest=dest,
                            columns_to_retain=columns_to_retain,
                            field_delimiter=args.field_delimiter)

    else:
        raise NotImplementedError


if __name__ == "__main__":
    main()


