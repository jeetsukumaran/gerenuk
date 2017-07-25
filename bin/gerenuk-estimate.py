#! /usr/bin/env python

import math
import csv
import os
import sys
import argparse
import collections
from gerenuk import utility

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    open = utility.pre_py34_open

class GerenukEstimator(object):

    def __init__(self,
            run_logger,
            stats_field_prefix="stat",
            logging_frequency=1000,
            field_delimiter="\t",
            is_output_summary_stats=False,
            is_suppress_checks=False,
            ):
        self.run_logger = run_logger
        self.stats_field_prefix = stats_field_prefix
        self.logging_frequency = logging_frequency
        self.field_delimiter = field_delimiter
        self.is_output_summary_stats = is_output_summary_stats
        self.is_suppress_checks = is_suppress_checks
        self.other_fieldnames = None
        self.stat_fieldnames = None
        self.stat_fieldnames_check = None
        self.other_fieldname_check = None
        self.stat_values = []
        self.other_values = []

    def read_simulated_data(self, filepaths):
        for filepath in filepaths:
            self.run_logger.info("Reading simulation file: '{}'".format(filepath))
            with open(filepath) as src:
                reader = csv.DictReader(
                        src,
                        delimiter=self.field_delimiter,
                        quoting=csv.QUOTE_NONE)
                for row_idx, row in enumerate(reader):
                    if self.logging_frequency and row_idx > 0 and row_idx % self.logging_frequency == 0:
                        self.run_logger.info("- Processing row {}".format(row_idx+1))
                    if self.stat_fieldnames is None:
                        self.stat_fieldnames = []
                        self.other_fieldnames = []
                        for field in reader.fieldnames:
                            if field.startswith(self.stats_field_prefix):
                                self.stat_fieldnames.append(field)
                            else:
                                self.other_fieldnames.append(field)
                        self.stat_fieldnames_check = set(self.stat_fieldnames)
                        self.other_fieldname_check = set(self.other_fieldnames)
                    row_stat_values = []
                    row_other_values = []
                    for key_idx, key in enumerate(row):
                        if not self.is_suppress_checks:
                            if key not in self.stat_fieldnames_check and key not in self.other_fieldname_check:
                                raise ValueError("File '{}', row {}, column {}: field '{}' not recognized".format(
                                    filepath, row_idx+1, key_idx+1, key))
                        if key.startswith(self.stats_field_prefix):
                            row_stat_values.append(float(row[key]))
                        else:
                            row_other_values.append( row[key] )
                    self.stat_values.append(row_stat_values)
                    self.other_values.append(row_other_values)

    def euclidean_distance(self, vector1, vector2):
        assert len(vector1) == len(vector2)
        dist = [(a - b)**2 for a, b in zip(vector1, vector2)]
        dist = math.sqrt(sum(dist))
        return dist

    def closest_values_indexes(self, target_stat_values, num_to_retain):
        assert len(target_stat_values) == len(self.stat_fieldnames), "Expecting {} values but found {}".format(
                len(self.stat_fieldnames),
                len(target_stat_values),
                )
        distances = []
        for idx, prior_value in enumerate(self.stat_values):
            d = self.euclidean_distance(prior_value, target_stat_values)
            distances.append((d, idx,))
        distances.sort(key=lambda x: x[0])
        return distances[:num_to_retain]

    def write_posterior(self,
            target_data_filepath,
            num_to_retain,
            ):
        with open(target_data_filepath) as src:
            reader = csv.DictReader(
                    src,
                    delimiter=self.field_delimiter,
                    quoting=csv.QUOTE_NONE)
            for row_idx, row in enumerate(reader):
                target_stat_values = []
                target_other_values = []
                for key_idx, key in enumerate(row):
                    if not self.is_suppress_checks:
                        if key not in self.stat_fieldnames_check and key not in self.other_fieldname_check:
                            raise ValueError("File '{}', target {}, column {}: field '{}' not recognized".format(
                                target_data_filepath, target_idx+1, key_idx+1, key))
                    if key.startswith(self.stats_field_prefix):
                        target_stat_values.append(float(row[key]))
                    else:
                        target_other_values.append( row[key] )
                posterior_indexes = self.closest_values_indexes(
                        target_stat_values=target_stat_values,
                        num_to_retain=num_to_retain,)
                dest = open(os.path.splitext(os.path.basename(target_data_filepath))[0] + ".posterior.{}.tsv".format(row_idx+1), "w")
                dest.write(self.field_delimiter.join(str(v) for v in self.other_fieldnames))
                if self.is_output_summary_stats:
                    dest.write(self.field_delimiter.join(str(v) for v in self.stat_fieldnames))
                dest.write("\n")
                for distance, index in posterior_indexes:
                    dest.write(self.field_delimiter.join(str(v) for v in self.other_values[index]))
                    if self.is_output_summary_stats:
                        dest.write(self.field_delimiter.join(str(v) for v in self.stat_values[index]))
                    dest.write("\n")

def main():
    parser = argparse.ArgumentParser(
            description="GERENUK Simultaneous Divergence Time Analysis -- Post-Process Columns",
            )
    parser.add_argument(
            "target_data_filepath",
            help="Path to target or observed data file.")
    parser.add_argument(
            "simulations_data_filepaths",
            nargs="+",
            help="Path to samples from the prior data files.")
    processing_options = parser.add_argument_group("Processing Options")
    processing_options.add_argument('--field-delimiter',
        type=str,
        default='\t',
        help="Field delimiter (default: <TAB>').")
    processing_options.add_argument('--stats-field-prefix',
        type=str,
        default='stat',
        help="Prefix identifying summary statistic fields (default: '%d(default)s').")
    output_options = parser.add_argument_group("Run Options")
    output_options.add_argument(
            "--output-summary-stats",
            action="store_true",
            help="Include summary stats in the samples from the posterior.")
    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Work silently.")
    args = parser.parse_args()

    run_logger = utility.RunLogger(
            name="gerenuk-estimate",
            stderr_logging_level="info",
            log_to_stderr=not args.quiet,
            log_to_file=False
            )
    ge = GerenukEstimator(
            run_logger=run_logger,
            stats_field_prefix=args.stats_field_prefix,
            field_delimiter=args.field_delimiter,
            is_output_summary_stats=args.output_summary_stats,
            )
    ge.read_simulated_data(args.simulations_data_filepaths)
    ge.write_posterior(
            target_data_filepath=args.target_data_filepath,
            num_to_retain=100,
            )

if __name__ == "__main__":
    main()



