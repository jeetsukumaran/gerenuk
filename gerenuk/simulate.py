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
import time
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing
import traceback

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
# DNA 10 0.00000 0.02 0.33
# """

FSC2_CONFIG_TEMPLATE = """\
//Number of population samples (demes)
2
//Population effective sizes (number of genes)
{d0_population_size}
{d1_population_size}
//Sample sizes
{d0_sample_size}
{d0_sample_size}
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
1  historical event
{div_time} 0 1 1 2 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 10 0.00000 0.02 0.33
"""

def weighted_choice(seq, weights, rng):
    """
    Selects an element out of seq, with probabilities of each element
    given by the list `weights` (which must be at least as long as the
    length of `seq` - 1).
    """
    if weights is None:
        weights = [1.0/len(seq) for count in range(len(seq))]
    else:
        weights = list(weights)
    if len(weights) < len(seq) - 1:
        raise Exception("Insufficient number of weights specified")
    sow = sum(weights)
    if len(weights) == len(seq) - 1:
        weights.append(1 - sow)
    return seq[weighted_index_choice(weights, sow, rng)]

def weighted_index_choice(weights, sum_of_weights, rng):
    """
    (From: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    The following is a simple function to implement weighted random choice in
    Python. Given a list of weights, it returns an index randomly, according
    to these weights [1].
    For example, given [2, 3, 5] it returns 0 (the index of the first element)
    with probability 0.2, 1 with probability 0.3 and 2 with probability 0.5.
    The weights need not sum up to anything in particular, and can actually be
    arbitrary Python floating point numbers.
    If we manage to sort the weights in descending order before passing them
    to weighted_choice_sub, it will run even faster, since the random call
    returns a uniformly distributed value and larger chunks of the total
    weight will be skipped in the beginning.
    """
    rnd = rng.uniform(0, 1) * sum_of_weights
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def sample_partition(
        number_of_elements,
        scaling_parameter,
        rng,):
    groups = []
    element_ids = [i for i in range(number_of_elements)] # decouple actual element index values, from indexes used to run Dirichlet process
    rng.shuffle(element_ids)                             # ... thus allowing for us to ensure the "first" index is randomized
    a = scaling_parameter
    for i, element_id in enumerate(element_ids):
        probs = []
        n = i + 1
        if i == 0:
            groups.append([element_id])
            continue
        p_new = a/(a + n - 1.0)
        probs.append(p_new)
        for group in groups:
            p = len(group)/(a + n - 1.0)
            probs.append(p)
        assert abs(sum(probs) - 1.0) <= 1e-5
        selected_idx = weighted_index_choice(
                weights=probs,
                sum_of_weights=1.0,
                rng=rng)
        if selected_idx == 0:
            groups.append([element_id])
        else:
            groups[selected_idx-1].append(element_id)
    return groups

class Deme(object):
    def __init__(self):
        self.population_size = 1000 # in number of genes, i.e. N for haploid or 2N for diploid
        self.sample_size = 20 # in number of genes

class LineagePair(object):
    def __init__(self):
        self.demes = (Deme(), Deme())

class GerenukSimulationModel(object):
    pass

    def __init__(self,
            rng):
        self.num_lineage_pairs = 2
        self.lineage_pairs = tuple(LineagePair() for i in range(self.num_lineage_pairs))
        self.rng = rng

        # shape and scale of Gamma hyperprior on
        # concentration parameter of Dirichlet process to partition pairs
        self.prior_num_divergences = (1000.0, 0.00437)
        # shape and scale of Gamma hyperprior on
        # theta (population) parameters for daughter deme
        self.prior_theta = (4.0, 0.001)
        # shape and scale of Gamma hyperprior on
        # theta (population) parameters for ancestral deme
        self.prior_ancestral_theta = (0, 0)
        # specification of fixed/free parameters
        self.theta_constraints = "000"
        # shape and scale of Gamma hyperprior on
        # divergence times
        self.prior_tau = (1.0, 0.02)

    def sample_parameter_values_from_prior(self):
        params = {}
        concentration_v = self.rng.gammavariate(*self.prior_num_divergences)
        params["param.concentration"] = concentration_v
        groups = sample_partition(
                number_of_elements=self.num_lineage_pairs,
                scaling_parameter=concentration_v, # sample from prior
                rng=self.rng,
                )
        params["param.numDivTimes"] = len(groups)
        div_time_values = [self.rng.gammavariate(*self.prior_tau) for i in groups]
        fsc2_run_configurations = [None for i in range(self.num_lineage_pairs)]
        div_time_model_desc = [None for i in range(self.num_lineage_pairs)]
        expected_lineage_pair_idxs = set([i for i in range(self.num_lineage_pairs)])
        for group_id, group in enumerate(groups):
            for lineage_pair_idx in group:
                assert lineage_pair_idx in expected_lineage_pair_idxs
                assert lineage_pair_idx not in fsc2_run_configurations
                fsc2_config_d = {}

                for deme_idx in range(2):
                    # effective population sizes in genes, of each lineage deme
                    # TODO: actually sample from prior instead of assigning!
                    deme_param_ne = self.lineage_pairs[lineage_pair_idx].demes[deme_idx].population_size
                    fsc2_config_d["d{}_population_size".format(deme_idx)] = deme_param_ne
                    params["param.populationSize.spp{}.deme{}".format(lineage_pair_idx, deme_idx)] = deme_param_ne

                    # num genes sampled, of each lineage deme
                    fsc2_config_d["d{}_sample_size".format(deme_idx)] = self.lineage_pairs[lineage_pair_idx].demes[deme_idx].sample_size

                    # TODO: specify data structure

                # divergence time model description
                div_time_model_desc[lineage_pair_idx] = str(group_id+1)

                # divergence time
                fsc2_config_d["div_time"] = div_time_values[group_id]
                params["param.divTime.spp{}".format(lineage_pair_idx)] = div_time_values[group_id]

                fsc2_run_configurations[lineage_pair_idx] = fsc2_config_d
        # code the model as 1100112
        # vector of div times
        # for each population pair --- other

        params["param.divModel"] = "".join(div_time_model_desc)
        return params, fsc2_run_configurations

class Fsc2Handler(object):

    def __init__(self, name, fsc2_path, working_directory):
        self.name = name
        self.fsc2_path = fsc2_path
        self.working_directory = working_directory
        self._is_file_system_staged = False
        self._num_executions = 0
        self._parameter_filepath = None
        self._results_dirpath = None
        self._current_execution_id = None

    def stage_filesystem(self):
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)
        self._is_file_system_staged = True

    def _get_current_execution_id(self):
        if self._current_execution_id is None:
            self._current_execution_id = "".join([self.name, "-{:06d}".format(self._num_executions),
                ])
        return self._current_execution_id
    current_execution_id = property(_get_current_execution_id)

    def _get_parameter_filepath(self):
        if self._parameter_filepath is None:
            # self._parameter_filepath = os.path.join(self.working_directory, ".".join([self.name, "par"]))
            self._parameter_filepath = ".".join([self.name, "par"])
        return self._parameter_filepath
    parameter_filepath = property(_get_parameter_filepath)

    def _get_results_dirpath(self):
        if self._results_dirpath is None:
            self._results_dirpath = os.path.splitext(self.parameter_filepath)[1]
        return self._results_dirpath
    results_dirpath = property(_get_results_dirpath)

    def _generate_parameter_file(self, fsc2_config_d):
        assert self.parameter_filepath
        with open(os.path.join(self.working_directory, self.parameter_filepath), "w") as dest:
            config = FSC2_CONFIG_TEMPLATE.format(**fsc2_config_d)
            dest.write(config)

    def _new_execution_reset(self):
        self._current_execution_id = None
        self._parameter_filepath = None

    def _setup_for_execution(self):
        self._new_execution_reset()
        if not self._is_file_system_staged:
            self.stage_filesystem()

    def _post_execution_cleanup(self):
        pass

    def run(self,
            fsc2_config_d,
            random_seed,):
        self._setup_for_execution()
        self._generate_parameter_file(fsc2_config_d)
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend(["-n", "1"]) # number of simulations to perform
        cmds.extend(["-r", str(random_seed)]) # seed for random number generator (positive integer <= 1E6)
        cmds.extend(["-d", "-s0", "-x", "-I", ])
        cmds.extend(["-i", self.parameter_filepath])
        p = subprocess.Popen(cmds,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.working_directory,
                )
        stdout, stderr = utility.communicate_process(p)
        self._num_executions += 1
        self._post_execution_cleanup()

class SimulationWorker(multiprocessing.Process):

    def __init__(self,
            name,
            model,
            work_queue,
            results_queue,
            fsc2_path,
            working_directory,
            run_logger,
            logging_frequency,
            messenger_lock,
            debug_mode,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.fsc2_handler = Fsc2Handler(
                name=name,
                fsc2_path=fsc2_path,
                working_directory=working_directory)
        self.model = model
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.run_logger = run_logger
        self.logging_frequency = logging_frequency
        self.messenger_lock = messenger_lock
        self.kill_received = False
        self.num_tasks_received = 0
        self.num_tasks_completed = 0
        self.is_debug_mode = debug_mode

    def send_worker_message(self, msg, level):
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

    def send_worker_critical(self, msg):
        self.send_worker_message(msg, utility.RunLogger.CRITICAL_MESSAGING_LEVEL)

    def send_worker_debug(self, msg):
        self.send_worker_message(msg, utility.RunLogger.DEBUG_MESSAGING_LEVEL)

    def send_worker_info(self, msg):
        self.send_worker_message(msg, utility.RunLogger.INFO_MESSAGING_LEVEL)

    def send_worker_warning(self, msg):
        self.send_worker_message(msg, utility.RunLogger.WARNING_MESSAGING_LEVEL)

    def send_worker_error(self, msg):
        self.send_worker_message(msg, utility.RunLogger.ERROR_MESSAGING_LEVEL)

    def run(self):
        result = None
        while not self.kill_received:
            try:
                rep_idx = self.work_queue.get_nowait()
            except queue.Empty:
                break
            self.num_tasks_received += 1
            # self.send_worker_critical("Received task: '{task_name}'".format(
            #     task_count=self.num_tasks_received,
            #     task_name=rep_idx))
            # rng = random.Random(random_seed)
            try:
                result = self.simulate()
            except (KeyboardInterrupt, Exception) as e:
                # traceback.print_exc()
                e.worker_name = self.name
                e.traceback_exc = traceback.format_exc()
                self.results_queue.put(e)
                break
            if self.kill_received:
                break
            self.results_queue.put(result)
            self.num_tasks_completed += 1
            # self.send_info("Completed task {task_count}: '{task_name}'".format(
            if rep_idx and self.logging_frequency and rep_idx % self.logging_frequency == 0:
                self.run_logger.info("Completed replicate {task_name}".format(
                    task_count=self.num_tasks_received,
                    task_name=rep_idx))
        if self.kill_received:
            self.send_worker_warning("Terminating in response to kill request")

    def simulate(self):
        result = {}
        params, fsc2_run_configurations = self.model.sample_parameter_values_from_prior()
        result.update(params)
        for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
            self.fsc2_handler.run(
                    fsc2_config_d=fsc2_run_configurations[lineage_pair_idx],
                    random_seed=self.model.rng.randint(1, 1E6),
                    )
        return result

class GerenukSimulator(object):

    def __init__(self,
            config_d,
            num_processes=None,
            logging_frequency=1000,
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
            self.run_logger.info("{} pairs of demes in analysis:".format(self.model.num_lineage_pairs))
            for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
                self.run_logger.info("    {:>3d}. {} / {}".format(
                        lineage_pair_idx+1,
                        lineage_pair.demes[0].sample_size,
                        lineage_pair.demes[1].sample_size,
                        ))

        self.worker_class = SimulationWorker

    def configure_simulator(self, config_d, verbose=True):
        self.name = config_d.pop("name", None)
        if self.name is None:
            self.name = str(id(self))
        self.title = "gerenuk-{}".format(self.name)
        self.output_prefix = config_d.pop("output_prefix", self.title)
        self.working_directory = config_d.pop("working_directory", self.title)
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
        self.logging_frequency = config_d.pop("logging_frequency", 1000)
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
        self.model = GerenukSimulationModel(
                rng=self.rng,
                )
        if config_d:
            raise Exception("Unrecognized configuration entries: {}".format(config_d))

    def execute(self, nreps):
        # load up queue
        self.run_logger.info("Creating work queue")
        work_queue = multiprocessing.Queue()
        for rep_idx in range(nreps):
            work_queue.put( rep_idx )
        time.sleep(0.1) # to avoid: 'IOError: [Errno 32] Broken pipe'; https://stackoverflow.com/questions/36359528/broken-pipe-error-with-multiprocessing-queue
        self.run_logger.info("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        for pidx in range(self.num_processes):
            worker = self.worker_class(
                    name=str(pidx+1),
                    model=self.model,
                    work_queue=work_queue,
                    results_queue=results_queue,
                    fsc2_path=self.fsc2_path,
                    working_directory=self.working_directory,
                    run_logger=self.run_logger,
                    logging_frequency=self.logging_frequency,
                    messenger_lock=messenger_lock,
                    debug_mode=self.is_debug_mode,
                    )
            worker.start()
            workers.append(worker)

        # collate results
        result_count = 0
        results_collator = []
        try:
            while result_count < nreps:
                result = results_queue.get()
                if isinstance(result, Exception) or isinstance(result, KeyboardInterrupt):
                    self.run_logger.error("Exception raised in worker process '{}'"
                                          "\n>>>\n{}<<<\n".format(
                                              result.worker_name,
                                              result.traceback_exc))
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
        return results_collator

if __name__ == "__main__":
    config_d = {"name": "test"}
    gs = GerenukSimulator(
            config_d=config_d,
            num_processes=3,
            is_verbose_setup=True)
    try:
        results = gs.execute(5)
    except Exception as e:
        sys.exit(1)
    print(results)
