#! /usr/bin/env python

##############################################################################
## Copyright (c) 2017 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import subprocess
import random
import sys
import os
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing

from gerenuk import utility

# FSC2_CONFIG_TEMPLATE = """\
# //Number of population samples (demes)
# 2
# //Population effective sizes (number of genes)
# 1000
# 1000
# //Sample sizes
# 100
# 100
# //Growth rates: negative growth implies population expansion
# 0
# 0
# //Number of migration matrices : 0 implies no migration between demes
# 0
# //historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
# 1  historical event
# 10000 0 1 1 2 0 0
# //Number of independent loci [chromosome]
# 1 0
# //Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
# 1
# //per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
# DNA 100 0 2.5e-6 OUTEXP
# """

# FSC2_CONFIG_TEMPLATE = """\
# //Number of population samples (demes)
# 2
# //Population effective sizes (number of genes)
# 1000
# 1000
# //Sample sizes
# 10
# 10
# //Growth rates: negative growth implies population expansion
# 0
# 0
# //Number of migration matrices : 0 implies no migration between demes
# 0
# //historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
# 1  historical event
# 10000 0 1 1 2 0 0
# //Number of independent loci [chromosome]
# 1 0
# //Per chromosome: Number of linkage blocks
# 1
# //per Block: data type, num loci, rec. rate and mut rate + optional parameters
# DNA 1 0.00000 0.02 0.33
# """

FSC2_CONFIG_TEMPLATE = """\
//Number of population samples (demes)
2
//Population effective sizes (number of genes)
1000
1000
//Sample sizes
100
100
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
1  historical event
10000 0 1 1 2 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 10 0.00000 0.02 0.33
"""

class SpeciesPair(object):

    index = 0

    def __init__(self):
        SpeciesPair.index += 1
        self.label = "S{:02}".format(SpeciesPair.index)
        self.num_population_genes = [1000, 1000]
        self.num_sampled_genes = [20, 20]
        self.divergence_time = 100000
        self.num_loci = 10
        self.mutation_rate = 1e-6
        self.transition_bias = 0.3

class GerenukSimulationModel(object):

    def __init__(self):
        self.num_species_pairs = 4
        self.species_pairs = [SpeciesPair() for i in range(self.num_species_pairs)]

class SimulationWorker(multiprocessing.Process):

    def __init__(self,
            name,
            work_queue,
            results_queue,
            fsc2_path,
            random_seed,
            run_logger,
            messenger_lock,
            debug_mode,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.fsc2_path = fsc2_path
        self.rng = random.Random(random_seed)
        self.run_logger = run_logger
        self.messenger_lock = messenger_lock
        self.kill_received = False
        self.num_tasks_received = 0
        self.num_tasks_completed = 0
        self.is_debug_mode = debug_mode

        self.num_executions = 0
        self._fsc2_parameter_filepath = None
        self._current_execution_id = None

    def _get_current_execution_id(self):
        if self._current_execution_id is None:
            self._current_execution_id = "".join([self.title, "-{:06d}".format(self.num_executions),
                ])
        return self._current_execution_id
    current_execution_id = property(_get_current_execution_id)

    def _get_fsc2_parameter_filepath(self):
        if self._fsc2_parameter_filepath is None:
            if self.is_debug_mode:
                self._fsc2_parameter_filepath = ".".join([self.current_execution_id, "par"])
            else:
                self._fsc2_parameter_filepath = ".".join([self.title, "par"])
        return self._fsc2_parameter_filepath
    fsc2_parameter_filepath = property(_get_fsc2_parameter_filepath)

    def _generate_fsc2_parameter_file(self):
        assert self.fsc2_parameter_filepath
        if self.is_verbose_setup:
            self.run_logger.info("Creating parameter file: '{}'".format(self.fsc2_parameter_filepath))
        with open(self.fsc2_parameter_filepath, "w") as dest:
            config = FSC2_CONFIG_TEMPLATE
            dest.write(config)

    def _new_execution_reset(self):
        self._current_execution_id = None
        self._fsc2_parameter_filepath = None

    def _setup_for_execution(self):
        self._new_execution_reset()
        self._generate_fsc2_parameter_file()

    def _post_execution_cleanup(self):
        pass

    def send_message(self, msg, level):
        if self.run_logger is None:
            return
        # if self.run_logger.messaging_level > level or self.messenger.silent:
        #     return
        msg = "{}: {}".format(self.name, msg)
        self.messenger_lock.acquire()
        try:
            self.run_logger.log(msg, level=level)
        finally:
            self.messenger_lock.release()

    def send_info(self, msg):
        self.send_message(msg, utility.RunLogger.INFO_MESSAGING_LEVEL)

    def send_warning(self, msg):
        self.send_message(msg, utility.RunLogger.WARNING_MESSAGING_LEVEL)

    def send_error(self, msg):
        self.send_message(msg, utility.RunLogger.ERROR_MESSAGING_LEVEL)

    def run(self):
        while not self.kill_received:
            try:
                rep_idx = self.work_queue.get_nowait()
            except queue.Empty:
                break
            self.num_tasks_received += 1
            # self.send_info("Received task {task_count}: '{task_name}'".format(
            self.send_info("Received task: '{task_name}'".format(
                task_count=self.num_tasks_received,
                task_name=rep_idx))
            try:
                pass
            except (KeyboardInterrupt, Exception) as e:
                e.worker_name = self.name
                self.results_queue.put(e)
                break
            if self.kill_received:
                break
            self.num_tasks_completed += 1
            # self.send_info("Completed task {task_count}: '{task_name}'".format(
            self.send_info("Completed task: '{task_name}'".format(
                task_count=self.num_tasks_received,
                task_name=rep_idx))
        if self.kill_received:
            self.send_warning("Terminating in response to kill request")
        else:
            self.results_queue.put(rep_idx)

    def execute(self):
        self._setup_for_execution()
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend(["-n", "1"]) # number of simulations to perform
        cmds.extend(["-r", str(self.rng.randint(1, 1E6))]) # seed for random number generator (positive integer <= 1E6)
        cmds.extend(["-d", "-s0", "-x", "-I", ])
        cmds.extend(["-i", self.fsc2_parameter_filepath])
        self.run_logger.info("Running FastSimCoal2")
        self.run_logger.info(cmds)
        p = subprocess.Popen(cmds,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )
        stdout, stderr = utility.communicate_process(p)
        self.num_executions += 1
        self._post_execution_cleanup()

class GerenukSimulator(object):

    def __init__(self,
            config_d,
            num_processes=None,
            is_verbose_setup=True):
        # configure
        self.elapsed_time = 0.0 # need to be here for logging
        config_d = dict(config_d) # make copy so we can pop items
        self.is_verbose_setup = is_verbose_setup
        self.configure_simulator(config_d)
        self.num_cpus = multiprocessing.cpu_count()
        if num_processes is None or num_processes <= 0:
            self.num_processes = num_cpus
        elif num_processes == 1 and self.num_cpus > 1 and self.is_verbose_setup:
            self.run_logger.info(
                    ("Multiple processors ({num_cpus}) available:"
                    " consider using the '-M' or '-m' options to"
                    " parallelize processing of trees"
                    ).format(num_cpus=num_cpus))
            self.num_processes = 1
        else:
            self.num_processes = num_processes
        if self.is_verbose_setup:
            self.run_logger.info("Will run up to {} processes in parallel".format(self.num_processes))
            self.run_logger.info("{} pairs of species in analysis:".format(self.model.num_species_pairs))
            for spidx, sp in enumerate(self.model.species_pairs):
                self.run_logger.info("    {:>3d}. '{}': {} / {}".format(spidx, sp.label, sp.num_sampled_genes[0], sp.num_sampled_genes[1]))

    def configure_simulator(self, config_d, verbose=True):
        self.name = config_d.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.title = "gerenuk-{}".format(self.name)
        self.output_prefix = config_d.pop("output_prefix", self.title)
        self.run_logger = config_d.pop("run_logger", None)
        if self.run_logger is None:
            self.run_logger = utility.RunLogger(
                    name="gerenuk",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_path=self.output_prefix + ".log",
                    file_logging_level=config_d.pop("file_logging_level", "info"),
                    )
        self.run_logger.system = self
        if self.is_verbose_setup:
            self.run_logger.info("Configuring simulation '{}'".format(self.name))
        self.fsc2_path = config_d.pop("fsc2_path", "fsc25")
        if self.is_verbose_setup:
            self.run_logger.info("FastSimCoal2 path: '{}'".format(self.fsc2_path))
        self.rng = config_d.pop("rng", None)
        if self.rng is None:
            self.random_seed = config_d.pop("random_seed", None)
            if self.random_seed is None:
                self.random_seed = random.randint(0, sys.maxsize)
            if self.is_verbose_setup:
                self.run_logger.info("Initializing with random seed {}".format(self.random_seed))
            self.rng = random.Random(self.random_seed)
        else:
            if "random_seed" in config_d:
                raise TypeError("Cannot specify both 'rng' and 'random_seed'")
            if self.is_verbose_setup:
                self.run_logger.info("Using existing random number generator")
        self.is_debug_mode = config_d.pop("debug_mode", False)
        if self.is_verbose_setup and self.is_debug_mode:
            self.run_logger.info("Running in DEBUG mode")
        self.model = GerenukSimulationModel()
        if config_d:
            raise Exception("Unrecognized configuration entries: {}".format(config_d))

    def execute(self, nreps):
        # load up queue
        self.run_logger.info("Creating work queue")
        work_queue = multiprocessing.Queue()
        for repidx in range(nreps):
            work_queue.put(repidx+1)
        self.run_logger.info("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        for pidx in range(self.num_processes):
            worker = SimulationWorker(
                    name=str(pidx+1),
                    work_queue=work_queue,
                    results_queue=results_queue,
                    fsc2_path=self.fsc2_path,
                    random_seed=self.rng.randint(1, 1E6),
                    run_logger=self.run_logger,
                    messenger_lock=messenger_lock,
                    debug_mode=self.is_debug_mode,
                    )
            worker.start()
            workers.append(worker)

        # collate results
        result_count = 0
        results_collator = []
        try:
            while result_count < self.num_processes:
                result = results_queue.get()
                if isinstance(result, Exception) or isinstance(result, KeyboardInterrupt):
                    # self.info_message("Exception raised in worker process '{}'".format(result.worker_name))
                    raise result
                results_collator.append(result)
                # self.run_logger.info("Recovered results from worker process '{}'".format(result.worker_name))
                result_count += 1
                # self.info_message("Recovered results from {} of {} worker processes".format(result_count, self.num_processes))
        except (Exception, KeyboardInterrupt) as e:
            for worker in workers:
                worker.terminate()
            raise
        self.run_logger.info("All {} worker processes terminated".format(self.num_processes))


if __name__ == "__main__":
    config_d = {"name": "test"}
    gs = GerenukSimulator(
            config_d=config_d,
            num_processes=3,
            is_verbose_setup=True)
    gs.execute(10)
