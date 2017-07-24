#! /usr/bin/env python

import os
import sys
import argparse
import traceback
import time
from gerenuk import simulate
from gerenuk import utility

def main():
    parser = argparse.ArgumentParser(
            description="GERENUK Simultaneous Divergence Time Analysis -- Site Frequency Spectrum Simulator",
            )

    model_options = parser.add_argument_group("Simulation Model")
    model_options.add_argument("model_file",
            metavar="MODEL-FILE",
            help="Path to file defining the model.")

    simulator_options = parser.add_argument_group("Simulation Options")
    simulator_options.add_argument("-n", "--num-reps",
            type=int,
            default=10,
            help="Number of replicates (default: %(default)s).")
    simulator_options.add_argument("-z", "--random-seed",
            default=None,
            help="Seed for random number generator engine.")
    simulator_options.add_argument("-s", "--site-frequency-spectrum-type",
            choices=["folded", "unfolded"],
            default="folded",
            help="Type of site frequency spectrum to generate, 'folded' or 'unfolded' (default: %(default)s).")

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument('--name',
        action='store',
        type=str,
        default=None,
        metavar='SIMULATION-NAME',
        help="Identifier for this simulation run (default: timestamp).")
    output_options.add_argument('-o', '--output-prefix',
        action='store',
        dest='output_prefix',
        type=str,
        default=None,
        metavar='OUTPUT-FILE-PREFIX',
        help="Prefix for output files (default: same as simulation name)').")
    output_options.add_argument('-w', '--working-directory',
        action='store',
        type=str,
        default=None,
        help="Directory for temporary files (default: '%(default)s').")
    output_options.add_argument("-l", "--labels",
            action="append",
            help="Labels to append to output (in format <FIELD-NAME>:value;)")
    output_options.add_argument( "--append",
            action="store_true",
            default=False,
            help="Append instead of overwriting output file(s).")
    output_options.add_argument('--summary-stats-label-prefix',
        type=str,
        default='stat',
        metavar='PREFIX',
        help="Prefix for field labels for summary statistics (default: '%(default)s').")
    output_options.add_argument( "--include-model-id-field",
            action="store_true",
            default=False,
            help="Include a 'model.id' field (with same value as 'param.divTimeModel' field) in output.")
    output_options.add_argument( "--no-write-header",
            action="store_true",
            default=False,
            help="Do not writer header row.")

    run_options = parser.add_argument_group("Run Options")
    run_options.add_argument("-m", "--num-processes",
            default=1,
            type=int,
            help="Number of processes/CPU to run (default: '%(default)s').")
    run_options.add_argument("--log-frequency",
            default=None,
            type=float,
            help="Frequency that background progress messages get written to the log (0: do not log informational messages).")
    run_options.add_argument("--file-logging-level",
            default=None,
            help="Message level threshold for file logs.")
    run_options.add_argument("--stderr-logging-level",
            default=None,
            help="Message level threshold for screen logs.")
    run_options.add_argument("--debug-mode",
            action="store_true",
            default=False,
            help="Run in debugging mode.")

    fsc2_options = parser.add_argument_group("FastSimCoal2 Options")
    fsc2_options.add_argument("--fsc2-path",
            metavar="FSC2-PATH",
            default="fsc25",
            help="Path to FastsimCoal2 application (default: %(default)s).")

    args = parser.parse_args()

    config_d = {}
    utility.parse_legacy_configuration(
            filepath=args.model_file,
            config_d=config_d)
    if args.name is None:
        config_d["name"] = time.strftime("%Y%m%d%H%M%S")
    else:
        config_d["name"] = args.name
    if args.output_prefix is None:
        config_d["output_prefix"] = os.path.splitext(os.path.basename(args.model_file))[0]
    else:
        config_d["output_prefix"] = args.output_prefix
    if args.working_directory is not None:
        config_d["working_directory"] = args.working_directory
    else:
        config_d["working_directory"] = "work" # TODO! use tempfile to get this
    config_d["logging_frequency"] = args.log_frequency
    config_d["fsc2_path"] = args.fsc2_path
    config_d["file_logging_level"] = args.file_logging_level
    config_d["standard_error_logging_level"] = args.stderr_logging_level
    # config_d["log_to_file"] = args.log_to_file
    # config_d["log_to_stderr"] = args.log_to_stderr
    config_d["site_frequency_spectrum_type"] = args.site_frequency_spectrum_type
    config_d["stat_label_prefix"] = args.summary_stats_label_prefix
    config_d["supplemental_labels"] = utility.parse_fieldname_and_value(args.labels)
    config_d["is_include_model_id_field"] = args.include_model_id_field
    gs = simulate.GerenukSimulator(
            config_d=config_d,
            num_processes=args.num_processes,
            is_verbose_setup=False)
    filepath = config_d["output_prefix"] + ".sumstats.csv"
    dest = utility.open_destput_file_for_csv_writer(filepath=filepath, is_append=args.append)
    if args.append or args.no_write_header:
        is_write_header = False
    else:
        is_write_header = True
    with dest:
        writer = utility.get_csv_writer(dest=dest)
        try:
            results = gs.execute(
                    nreps=args.num_reps,
                    results_csv_writer=writer,
                    results_store=None,
                    is_write_header=is_write_header,
                    column_separator=",")
        except Exception as e:
            sys.stderr.write("Traceback (most recent call last):\n  {}{}\n".format(
                "  ".join(traceback.format_tb(sys.exc_info()[2])),
                e))
            sys.exit(1)

if __name__ == "__main__":
    main()

