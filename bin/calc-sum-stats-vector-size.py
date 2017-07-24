#! /usr/bin/env python

###############################################################################
##
##  Copyright 2017 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

import sys
import os
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--num-lineage-pairs",
            default=10,
            type=int,
            help="Number of lineage pairs. Default: %(default)s.")
    parser.add_argument("-l", "--num-samples-loci",
            default=10,
            type=int,
            help="Number of loci sampled. Default: %(default)s.")
    parser.add_argument("-n", "--num-samples-deme",
            default=10,
            type=int,
            help="Number of individuals sampled from each deme. Default: %(default)s.")
    args = parser.parse_args()

    num_ss_per_locus = (args.num_samples_deme * 2) + (args.num_samples_deme * args.num_samples_deme)
    num_ss = num_ss_per_locus * args.num_samples_loci * args.num_lineage_pairs
    print(num_ss)


if __name__ == '__main__':
    main()


