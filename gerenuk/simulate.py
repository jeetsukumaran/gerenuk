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
import collections
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

FSC2_CONFIG_TEMPLATE = """\
//Number of population samples (demes)
2
//Population effective sizes (number of genes)
{d0_population_size}
{d1_population_size}
//Sample sizes
{d0_sample_size}
{d1_sample_size}
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
1  historical event
{div_time} 0 1 1 2 0 0
//Number of independent loci [chromosome]; '0' => same structure for all loci
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination rate, per generation mutation rate and optional parameters
DNA {num_sites} {recombination_rate} {mutation_rate} {ti_proportional_bias}
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

def compose_lineage_pair_label(lineage_pair_idx):
    return "spp{}".format(lineage_pair_idx)

def compose_deme_label(deme_idx):
    return "deme{}".format(deme_idx)
_DEME0_LABEL = compose_deme_label(0)
_DEME1_LABEL = compose_deme_label(1)
_ANCESTOR_DEME_LABEL = compose_deme_label("A")

class LocusDefinition(object):

    def __init__(self, locus_d):
        self.configure(locus_d)

    def configure(self, locus_d):
        # Doc/comments for parameters from, and following, PyMsBayes (Jamie Oaks; https://github.com/joaks1/PyMsBayes)
        # label for this locus
        self.locus_label = locus_d.pop("locus_label")
        # The number in this column is used to scale for differences in ploidy among loci
        # or for differences in generation-times among taxa. In our example configuration
        # file 1.0 is used for loci from a diploid nuclear genome, whereas 0.25 is used
        # for a region of the mitochondrial genome (because its haploid and maternally
        # inherited). However, if a taxon "species-3" had 1/4 the generation times of
        # the other two taxa, we would specify "1.0" for the third column for its
        # mitochondrial locus, and "4.0" for the third column for its nuclear loci.
        self.ploidy_factor = float(locus_d.pop("ploidy_factor"))
        # The number in this column is used to scale for differences in mutation rates
        # among taxa and/or loci.
        self.mutation_rate_factor = float(locus_d.pop("mutation_rate_factor"))
        # Number of genes/sequences/taxa from first population
        self.num_genes_deme0 = int(locus_d.pop("num_genes_deme0"))
        # Number of genes/sequences/taxa from second population
        self.num_genes_deme1 = int(locus_d.pop("num_genes_deme1"))
        # This is the transition/transversion rate ratio ("Kappa") of the HKY85
        # model of nucleotide substitution [4] for this alignment. NOTE: This is
        # the transition/transversion rate ratio, not the "count" ratio. I.e.,
        # Kappa = 1 is equal to the Jukes-Cantor model.
        self.ti_tv_rate_ratio = float(locus_d.pop("ti_tv_rate_ratio"))
        # Number of sites
        self.num_sites = int(locus_d.pop("num_sites"))
        # Equilibrium frequency of nucleotide
        self.freq_a = float(locus_d.pop("freq_a"))
        # Equilibrium frequency of nucleotide
        self.freq_c = float(locus_d.pop("freq_c"))
        # Equilibrium frequency of nucleotide
        self.freq_g = float(locus_d.pop("freq_g"))
        # Path to alignment file (optional)
        self.alignment_filepath = locus_d.pop("alignment_filepath", None)
        # Done!
        if locus_d:
            raise Exception("Unrecognized locus definition entries: {}".format(locus_d))

class LineagePair(object):

    def __init__(self, taxon_label):
        self.taxon_label = taxon_label
        self.locus_definitions = []

    def add_locus_definition(self, locus_d):
        locus = LocusDefinition(locus_d)
        self.locus_definitions.append(locus)
        return locus

class GerenukSimulationModel(object):

    def __init__(self, params_d, locus_info,):
        self.configure_loci(locus_info)
        self.configure_params(params_d)

    def configure_loci(self, locus_info):
        self.label_to_lineage_pair_map = {}
        self.lineage_pairs = [] # list to maintain order for indexing during Dirichlet process partitioning
        self.lineage_pairs_loci_labels = {}
        for locus_d in locus_info:
            taxon_label = locus_d.pop("taxon_label")
            try:
                lineage_pair = self.label_to_lineage_pair_map[taxon_label]
            except KeyError:
                lineage_pair = LineagePair(
                        taxon_label=taxon_label)
                self.label_to_lineage_pair_map[taxon_label] = lineage_pair
                self.lineage_pairs.append(lineage_pair)
                self.lineage_pairs_loci_labels[lineage_pair] = set()
            locus_definition = lineage_pair.add_locus_definition(locus_d)
            if locus_definition.locus_label in self.lineage_pairs_loci_labels[lineage_pair]:
                raise ValueError("Lineage pair '{}': locus with label '{}' has already been defined".format(lineage_pair.taxon_label, locus_definition.locus_label))
            self.lineage_pairs_loci_labels[lineage_pair].add(locus_definition.locus_label)

    def configure_params(self, params_d):
        # Doc/comments for parameters from, and following, PyMsBayes (Jamie Oaks; https://github.com/joaks1/PyMsBayes)
        params_d = utility.CaseInsensitiveDict(params_d)
        # Shape and scale of Gamma hyperprior on
        # concentration parameter of Dirichlet process to partition pairs
        self.prior_concentration = (
                float(params_d.pop("concentrationShape")),
                float(params_d.pop("concentrationScale"))
                )

        # # Shape and scale of Gamma hyperprior on theta PyMsBayes does not
        # # independently model N and \mu, but FastsimCoal2 does.
        self.prior_theta = (
                float(params_d.pop("thetaShape")),
                float(params_d.pop("thetaScale"))
                )

        # Shape and scale of Gamma hyperprior on population size.
        # PyMsBayes does not independently model N and \mu. Here we do.
        # BEAST* uses a gamma distribution with a mean 2\psi
        # and a shape of 2, with user specifying a (hyper-)prior on \psi.
        # Here, for now, we just use Gamma directly, with default shape
        # parameter of 2.
        # self.prior_popsize = (
        #         float(params_d.pop("popsizeShape", 2))
        #         float(params_d.pop("popsizeScale"))
        #         )
        # # Shape and scale of Gamma hyperprior on mutation rate.
        # # PyMsBayes does not independently model N and \mu. Here we do.
        # self.prior_mutRate = (
        #         float(params_d.pop("mutRateShape"))
        #         float(params_d.pop("mutRateScale"))
        #         )

        # Shape and scale of Gamma hyperprior on
        # theta (population) parameters for ancestral deme
        self.prior_ancestral_theta = (
                float(params_d.pop("ancestralThetaShape")),
                float(params_d.pop("ancestralThetaScale"))
                )
        # specification of fixed/free parameters
        self.theta_constraints = str(params_d.pop("thetaParameters", "000"))
        if len(self.theta_constraints) != 3:
            raise ValueError("Incorrectly specified 'thetaParameters' constraints: '{}'".format(self.theta_constraints))
        for idx, i in enumerate(self.theta_constraints):
            if i not in ["0", "1", "2"]:
                raise ValueError("Incorrectly specified 'thetaParameters' constraints: '{}'".format(self.theta_constraints))
        # Shape and scale of Gamma hyperprior on
        # divergence times
        self.prior_tau = (
                float(params_d.pop("tauShape")),
                float(params_d.pop("tauScale"))
                )
        # Shape and scale of Gamma hyperprior on
        # divergence times
        self.prior_migration = (
                float(params_d.pop("migrationShape")),
                float(params_d.pop("migrationScale"))
                )
        if self.prior_migration[0] != 0 or self.prior_migration[1] != 0:
            raise NotImplementedError("Migration is not yet supported")
        # 1: Time units are in expected substitutions per site. For example, a
        #    divergence of 0.05 means that, on average, 5% of sites have changed
        #    since the populations diverged (so you expect 10% divergence between
        #    the populations since the population divergence). Thus, you can
        #    convert these units to the number of generations or years by dividing
        #    by the mutation rate.
        tss = int(params_d.pop("timeInSubsPerSite"))
        self.time_in_subs_per_site = bool(tss)
        if not self.time_in_subs_per_site:
            raise NotImplementedError("Time not in expected substitutions per site not supported")
        # If both are positive, these settings define a beta prior on the magnitude of a
        # post-divergence bottleneck in each of the descendant populations.
        # bottleProportionShapeA and bottleProportionShapeB correspond to the shape
        # parameters alpha and beta, respectively, of the beta prior.
        # The bottleneck magnitude is the proportion of the effective population size
        # that remains during the bottleneck. For example, a value of 0.95 would mean
        # that bottleneck reduces the effective population size by 5%.
        # If either or both are zero or less, there is no post-divergence population
        # bottleneck in the descendant populations (i.e., the bottleneck-magnitude
        # parameters, along with the timing of each bottleneck, are removed from
        # the model).
        self.prior_bottleneck_proportion = (
                float(params_d.pop("bottleProportionShapeA", 0)),
                float(params_d.pop("bottleProportionShapeB", 0)),
                )
        # If bottleProportionShared = 0, then there are two free
        # bottleneck-magnitude parameters for each population pair (one for
        # each descendant population). If bottleProportionShared = 1, then
        # there is one bottleneck-magnitude parameter for each population pair
        # (i.e., the descendant populations of each pair share the same
        # bottleneck magnitude; the bottleneck magnitude still varies among the
        # pairs).
        self.bottle_proportion_shared = bool(int(params_d.pop("bottleProportionShared")))
        # If this setting is zero (the default), the number of divergence
        # events is free to vary according to the Dirichlet process prior on
        # divergence models. If it is greater than zero, then the model is
        # constrained to numTauClasses divergence events. This is useful for
        # simulation-based power analyses, but should not be used for empirical
        # analyses.
        self.num_tau_classes = int(params_d.pop("numTauClasses"))
        # Done!
        if params_d:
            raise Exception("Unrecognized parameter configuration entries: {}".format(params_d))

    def _get_num_lineage_pairs(self):
        return len(self.lineage_pairs)
    num_lineage_pairs = property(_get_num_lineage_pairs)

    def sample_parameter_values_from_prior(self, rng):
        params = collections.OrderedDict()

        ## div time
        concentration_v = rng.gammavariate(*self.prior_concentration)
        params["param.divTimeModel"] = "NA" # initialize here, so first column
        # params["param.concentration"] = concentration_v
        groups = sample_partition(
                number_of_elements=self.num_lineage_pairs,
                scaling_parameter=concentration_v, # sample from prior
                rng=rng,
                )
        params["param.numDivTimes"] = len(groups)
        div_time_values = [rng.gammavariate(*self.prior_tau) for i in groups]
        fsc2_run_configurations = collections.OrderedDict()
        div_time_model_desc = [None for i in range(self.num_lineage_pairs)]

        # thetas
        theta = rng.gammavariate(*self.prior_theta)

        expected_lineage_pair_idxs = set([i for i in range(self.num_lineage_pairs)])
        # groups sorted by earliest occurring lineage pair index, to retain consistent div time coding
        for group_id, group in enumerate(sorted(groups, key=lambda group: min(lineage_pair_idx for lineage_pair_idx in group))):
            for lineage_pair_idx in group:
                assert lineage_pair_idx in expected_lineage_pair_idxs
                assert lineage_pair_idx not in fsc2_run_configurations
                lineage_pair = self.lineage_pairs[lineage_pair_idx]
                ## divergence time
                div_time_model_desc[lineage_pair_idx] = str(group_id+1) # divergence time model description
                div_time = div_time_values[group_id]
                params["param.divTime.{}".format(lineage_pair.taxon_label)] = div_time

                ## population parameters --- separate N and mu parameterization
                ## -- deme 0
                # deme0_pop_size = rng.gammavariate(*self.prior_popsize)
                # deme0_mu = rng.gammavariate(*self.prior_mutRate)
                # deme0_theta = 4 * deme0_pop_size * deme_0_mu
                # params["param.popSize.{}.{}".format(lineage_pair.taxon_label, _DEME0_LABEL)] = deme0_pop_size
                # params["param.mutRate.{}.{}".format(lineage_pair.taxon_label, _DEME0_LABEL)] = deme0_mu
                # params["param.theta.{}.{}".format(lineage_pair.taxon_label, _DEME0_LABEL)] = deme0_theta
                # ## -- deme 1
                # if self.theta_constraints[1] == self.theta_constraints[0]:
                #     deme1_pop_size = deme0_pop_size
                #     deme1_mu = deme0_mu
                #     deme1_theta = deme0_theta
                # else:
                #     deme1_pop_size = rng.gammavariate(*self.prior_popsize)
                #     deme1_mu = rng.gammavariate(*self.prior_mutRate)
                #     deme1_theta = 4 * deme1_pop_size * deme1_mu
                # params["param.popSize.{}.{}".format(lineage_pair.taxon_label, _DEME1_LABEL)] = deme1_pop_size
                # params["param.mutRate.{}.{}".format(lineage_pair.taxon_label, _DEME1_LABEL)] = deme1_mu
                # params["param.theta.{}.{}".format(lineage_pair.taxon_label, _DEME1_LABEL)] = deme1_theta
                # ## -- ancestor deme 2
                # if self.theta_constraints[2] == self.theta_constraints[0]:
                #     deme2_pop_size = deme0_pop_size
                #     deme2_mu = deme0_mu
                #     deme2_theta = deme0_theta
                # elif self.theta_constraints[2] == self.theta_constraints[1]:
                #     deme2_pop_size = deme1_pop_size
                #     deme2_mu = deme1_mu
                #     deme2_theta = deme1_theta
                # elif self.prior_ancestral_theta[0] != 0 and self.prior_ancestral_theta[1] != 0:
                #     raise NotImplementedError
                #     deme2_pop_size = rng.gammavariate(*self.prior_popsize)
                #     deme2_mu = rng.gammavariate(*self.prior_mutRate)
                #     deme2_theta = 4 * deme2_pop_size * deme2_mu
                # else:
                #     deme2_pop_size = rng.gammavariate(*self.prior_popsize)
                #     deme2_mu = rng.gammavariate(*self.prior_mutRate)
                #     deme2_theta = 4 * deme2_pop_size * deme2_mu
                # params["param.popSize.{}.{}".format(lineage_pair.taxon_label, _ANCESTOR_DEME_LABEL)] = deme2_pop_size
                # params["param.mutRate.{}.{}".format(lineage_pair.taxon_label, _ANCESTOR_DEME_LABEL)] = deme2_mu
                # params["param.theta.{}.{}".format(lineage_pair.taxon_label, _ANCESTOR_DEME_LABEL)] = deme2_theta

                ## population parameters --- theta parameterization
                deme0_theta = rng.gammavariate(*self.prior_theta)
                if self.theta_constraints[1] == self.theta_constraints[0]:
                    deme1_theta = deme0_theta
                else:
                    deme1_theta = rng.gammavariate(*self.prior_theta)
                if self.theta_constraints[2] == self.theta_constraints[0]:
                    deme2_theta = deme0_theta
                elif self.theta_constraints[2] == self.theta_constraints[1]:
                    deme2_theta = deme1_theta
                elif self.prior_ancestral_theta[0] != 0 and self.prior_ancestral_theta[1] != 0:
                    deme2_theta = rng.gammavariate(*self.prior_ancestral_theta)
                else:
                    deme2_theta = rng.gammavariate(*self.prior_theta)
                params["param.theta.{}.{}".format(lineage_pair.taxon_label, _DEME0_LABEL)] = deme0_theta
                params["param.theta.{}.{}".format(lineage_pair.taxon_label, _DEME1_LABEL)] = deme1_theta
                params["param.theta.{}.{}".format(lineage_pair.taxon_label, _ANCESTOR_DEME_LABEL)] = deme2_theta

                for locus_id, locus_definition in enumerate(lineage_pair.locus_definitions):
                    # Fastsimecoal2 separates pop size and mutation rate, but
                    # the msBayes/PyMsBayes model does not separate the two,
                    # using theta.
                    #
                    # We could just reparameterize the PyMsBayes model here,
                    # sampling over N and mu independently. But say we want to
                    # stick to the theta parameterization.
                    #
                    # We could simply scale everything by mutation rate --
                    # i.e., population size and div time specified in units of
                    # N * mu, and in the sequence generation assume a base
                    # mutation rate of 1.0. Problem with this is that
                    # Fastsimcoal coerces the population size to an integer,
                    # and so anything less than 1 becomes zero. So we multiply
                    # the population size time by a large number and adjust
                    # this in the actual mutation rate so N mu remains the
                    # same:
                    #
                    #   theta = 4 N mu = 4 * (N * C) * (mu/C)
                    #
                    # Of course, VERY important to also apply the adjustment
                    # factor to the divergence time, or, indeed, any other time
                    # variable!
                    adjustment_hack = 1E8
                    #
                    fsc2_config_d = {
                        "d0_population_size": deme0_theta/4.0 * locus_definition.ploidy_factor * adjustment_hack,
                        "d1_population_size": deme1_theta/4.0 * locus_definition.ploidy_factor * adjustment_hack,
                        "d0_sample_size": locus_definition.num_genes_deme0,
                        "d1_sample_size": locus_definition.num_genes_deme1,
                        "div_time": div_time * adjustment_hack, # ditto
                        "num_sites": locus_definition.num_sites,
                        "recombination_rate": 0,
                        "mutation_rate": locus_definition.mutation_rate_factor / adjustment_hack,
                        "ti_proportional_bias": (1.0 * locus_definition.ti_tv_rate_ratio)/3.0,
                        }
                    fsc2_run_configurations[locus_definition] = fsc2_config_d
        params["param.divTimeModel"] = "Model{}".format("".join(div_time_model_desc))
        return params, fsc2_run_configurations

class Fsc2RuntimeError(RuntimeError):
    def __init__(self, msg):
        RuntimeError.__init__(self, msg)

class Fsc2Handler(object):

    def __init__(self,
            name,
            fsc2_path,
            working_directory,
            is_calculate_single_population_sfs,
            is_calculate_joint_population_sfs,
            is_unfolded_site_frequency_spectrum,
            ):
        self.name = name
        self.fsc2_path = fsc2_path
        self.working_directory = working_directory
        self.is_unfolded_site_frequency_spectrum = is_unfolded_site_frequency_spectrum
        if self.is_unfolded_site_frequency_spectrum:
            self.sfs_file_prefix = "DAF"
            self.fsc2_sfs_generation_command = "-d"
        else:
            self.sfs_file_prefix = "MAF"
            self.fsc2_sfs_generation_command = "-m"
        self.is_calculate_single_population_sfs = is_calculate_single_population_sfs
        self.is_calculate_joint_population_sfs = is_calculate_joint_population_sfs
        self._is_file_system_staged = False
        self._num_executions = 0
        self._current_execution_id = None
        self._parameter_filepath = None
        self._results_dirpath = None
        self._deme0_site_frequency_filepath = None
        self._deme1_site_frequency_filepath = None
        self._joint_site_frequency_filepath = None

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
            self._results_dirpath = os.path.join(self.working_directory, os.path.splitext(self.parameter_filepath)[0])
        return self._results_dirpath
    results_dirpath = property(_get_results_dirpath)

    def _get_result_deme0_site_frequency_filepath(self):
        if self._deme0_site_frequency_filepath is None:
            self._deme0_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_{}pop0.obs".format(self.name, self.sfs_file_prefix))
        return self._deme0_site_frequency_filepath
    deme0_site_frequency_filepath = property(_get_result_deme0_site_frequency_filepath)

    def _get_result_deme1_site_frequency_filepath(self):
        if self._deme1_site_frequency_filepath is None:
            self._deme1_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_{}pop1.obs".format(self.name, self.sfs_file_prefix))
        return self._deme1_site_frequency_filepath
    deme1_site_frequency_filepath = property(_get_result_deme1_site_frequency_filepath)

    def _get_result_joint_site_frequency_filepath(self):
        if self._joint_site_frequency_filepath is None:
            self._joint_site_frequency_filepath = os.path.join(self.results_dirpath, "{}_joint{}pop1_0.obs".format(self.name, self.sfs_file_prefix))
        return self._joint_site_frequency_filepath
    joint_site_frequency_filepath = property(_get_result_joint_site_frequency_filepath)

    def _new_execution_reset(self):
        self._current_execution_id = None
        self._parameter_filepath = None

    def _setup_for_execution(self):
        self._new_execution_reset()
        if not self._is_file_system_staged:
            self._stage_filesystem()

    def _stage_filesystem(self):
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)
        self._is_file_system_staged = True

    def _generate_parameter_file(self, fsc2_config_d):
        assert self.parameter_filepath
        with open(os.path.join(self.working_directory, self.parameter_filepath), "w") as dest:
            self._write_parameter_configuration(
                    dest=dest,
                    fsc2_config_d=fsc2_config_d,
                    )

    def _write_parameter_configuration(self, dest, fsc2_config_d):
            config = FSC2_CONFIG_TEMPLATE.format(**fsc2_config_d)
            dest.write(config)

    def _parse_deme_derived_allele_frequencies(self,
            filepath,
            field_name_prefix,
            results_d):
        # results_d = collections.OrderedDict()
        with open(filepath) as src:
            lines = src.read().split("\n")
            assert len(lines) == 4 and lines[3] == ""
            header_row = lines[1].split("\t")
            results_d_row = lines[2].split("\t")
            assert len(header_row) == len(results_d_row)
            for key, val in zip(header_row, results_d_row):
                if not val:
                    continue
                results_d["{}.{}".format(field_name_prefix, key)] = float(val)
        return results_d

    def _parse_joint_derived_allele_frequencies(self,
            filepath,
            field_name_prefix,
            results_d):
        # results_d = collections.OrderedDict()
        with open(filepath) as src:
            lines = src.read().split("\n")
            col_keys = lines[1].split("\t")[1:]
            for line in lines[2:]:
                if not line:
                    continue
                cols = line.split("\t")
                assert len(cols) - 1 == len(col_keys)
                row_key = cols[0]
                for col_key, val in zip(col_keys, cols[1:]):
                    results_d["{}.{}.{}".format(field_name_prefix, row_key, col_key)] = float(val)
        return results_d

    def _harvest_run_results(self, field_name_prefix, results_d):
        if self.is_calculate_single_population_sfs:
            self._parse_deme_derived_allele_frequencies(
                    filepath=self.deme0_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(0)),
                    results_d=results_d)
            self._parse_deme_derived_allele_frequencies(
                    filepath=self.deme1_site_frequency_filepath,
                    field_name_prefix="{}.{}.sfs".format(field_name_prefix, compose_deme_label(1)),
                    results_d=results_d)
        if self.is_calculate_joint_population_sfs:
            self._parse_joint_derived_allele_frequencies(
                    filepath=self.joint_site_frequency_filepath,
                    field_name_prefix="{}.joint.sfs".format(field_name_prefix),
                    results_d=results_d)
        return results_d

    def _post_execution_cleanup(self):
        pass

    def run(self,
            field_name_prefix,
            fsc2_config_d,
            random_seed,
            results_d,):
        self._setup_for_execution()
        self._generate_parameter_file(fsc2_config_d)
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend(["-n", "1"]) # number of simulations to perform
        cmds.extend(["-r", str(random_seed)]) # seed for random number generator (positive integer <= 1E6)
        cmds.append(self.fsc2_sfs_generation_command)
        cmds.extend(["-s0", "-x", "-I", ])
        cmds.extend(["-i", self.parameter_filepath])
        p = subprocess.Popen(cmds,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.working_directory,
                )
        stdout, stderr = utility.communicate_process(p)
        if p.returncode != 0:
            raise Fsc2RuntimeError("FastSimCoal2 execution failure: {}".format(stderr))
        self._num_executions += 1
        if results_d is None:
            results_d = collections.OrderedDict()
        self._harvest_run_results(
                field_name_prefix=field_name_prefix,
                results_d=results_d)
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
            random_seed,
            is_calculate_single_population_sfs,
            is_calculate_joint_population_sfs,
            is_unfolded_site_frequency_spectrum,
            stat_label_prefix,
            is_include_model_id_field,
            supplemental_labels,
            debug_mode,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.fsc2_handler = Fsc2Handler(
                name=name,
                fsc2_path=fsc2_path,
                working_directory=working_directory,
                is_calculate_single_population_sfs=is_calculate_single_population_sfs,
                is_calculate_joint_population_sfs=is_calculate_joint_population_sfs,
                is_unfolded_site_frequency_spectrum=is_unfolded_site_frequency_spectrum)
        self.model = model
        self.rng = random.Random(random_seed)
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.run_logger = run_logger
        self.logging_frequency = logging_frequency
        self.messenger_lock = messenger_lock
        self.is_unfolded_site_frequency_spectrum = is_unfolded_site_frequency_spectrum
        self.stat_label_prefix = stat_label_prefix
        self.is_include_model_id_field = is_include_model_id_field
        self.supplemental_labels = supplemental_labels
        self.is_debug_mode = debug_mode
        self.kill_received = False
        self.num_tasks_received = 0
        self.num_tasks_completed = 0

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
        results_d = collections.OrderedDict()
        if self.is_include_model_id_field:
            results_d["model.id"] = None
        if self.supplemental_labels:
            for key in self.supplemental_labels:
                results_d[key] = self.supplemental_labels[key]
        params, fsc2_run_configurations = self.model.sample_parameter_values_from_prior(rng=self.rng)
        results_d.update(params)
        for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
            for locus_definition in lineage_pair.locus_definitions:
                self.fsc2_handler.run(
                        field_name_prefix="{}.{}.{}".format(
                                self.stat_label_prefix,
                                lineage_pair.taxon_label,
                                locus_definition.locus_label),
                        fsc2_config_d=fsc2_run_configurations[locus_definition],
                        random_seed=self.rng.randint(1, 1E6),
                        results_d=results_d,
                        )
        if self.is_include_model_id_field:
            results_d["model.id"] = results_d["param.divTimeModel"]
        return results_d

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
                    ).format(num_cpus=self.num_cpus))
            self.num_processes = 1
        else:
            self.num_processes = num_processes
        if self.is_verbose_setup:
            self.run_logger.info("Will run up to {} processes in parallel".format(self.num_processes))
            self.run_logger.info("{} lineage pairs in analysis:".format(self.model.num_lineage_pairs))
            for lineage_pair_idx, lineage_pair in enumerate(self.model.lineage_pairs):
                self.run_logger.info("  - '{}': {:>2d} loci (Samples: {})".format(
                        lineage_pair.taxon_label,
                        len(lineage_pair.locus_definitions),
                        ", ".join("{}/{}".format(locus.num_genes_deme0, locus.num_genes_deme1) for locus in lineage_pair.locus_definitions),
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
                    name="gerenuk-simulate",
                    stderr_logging_level=config_d.pop("standard_error_logging_level", "info"),
                    log_to_file=config_d.pop("log_to_file", True),
                    log_to_stderr=config_d.pop("log_to_stderr", True),
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
        if self.is_verbose_setup:
            self.run_logger.info("Working directory: '{}'".format(self.working_directory))
        self.is_debug_mode = config_d.pop("debug_mode", False)
        if self.is_verbose_setup and self.is_debug_mode:
            self.run_logger.info("Running in DEBUG mode")
        self.site_frequency_spectrum_type = config_d.pop("site_frequency_spectrum_type", "unfolded").lower()
        self.is_unfolded_site_frequency_spectrum = config_d.pop("is_unfolded_site_frequency_spectrum", False)
        self.is_calculate_single_population_sfs = config_d.pop("is_calculate_single_population_sfs", False)
        self.is_calculate_joint_population_sfs = config_d.pop("is_calculate_joint_population_sfs", True)
        if not self.is_calculate_single_population_sfs and not self.is_calculate_joint_population_sfs:
            raise ValueError("Neither single-population nor joint site frequency spectrum will be calculated!")
        self.stat_label_prefix = config_d.pop("stat_label_prefix", "stat")
        self.supplemental_labels = config_d.pop("supplemental_labels", None)
        self.is_include_model_id_field = config_d.pop("is_include_model_id_field", False)
        if "params" not in config_d:
            raise ValueError("Missing 'params' entry in configuration")
        params_d = config_d.pop("params")
        if "locus_info" not in config_d:
            raise ValueError("Missing 'locus_info' entry in configuration")
        locus_info = config_d.pop("locus_info")
        self.model = GerenukSimulationModel(params_d=params_d, locus_info=locus_info,)
        if config_d:
            raise Exception("Unrecognized configuration entries: {}".format(config_d))

    def execute(self,
            nreps,
            results_csv_writer=None,
            results_store=None,
            is_write_header=True,
            ):
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
                    # name=str(pidx+1),
                    name="{}-{}".format(self.title, pidx+1),
                    model=self.model,
                    work_queue=work_queue,
                    results_queue=results_queue,
                    fsc2_path=self.fsc2_path,
                    working_directory=self.working_directory,
                    run_logger=self.run_logger,
                    logging_frequency=self.logging_frequency,
                    messenger_lock=messenger_lock,
                    random_seed=self.rng.randint(1, sys.maxint),
                    is_calculate_single_population_sfs=self.is_calculate_single_population_sfs,
                    is_calculate_joint_population_sfs=self.is_calculate_joint_population_sfs,
                    is_unfolded_site_frequency_spectrum=self.is_unfolded_site_frequency_spectrum,
                    stat_label_prefix=self.stat_label_prefix,
                    is_include_model_id_field=self.is_include_model_id_field,
                    supplemental_labels=self.supplemental_labels,
                    debug_mode=self.is_debug_mode,
                    )
            worker.start()
            workers.append(worker)

        # collate results
        result_count = 0
        try:
            while result_count < nreps:
                result = results_queue.get()
                if isinstance(result, KeyboardInterrupt):
                    raise result
                elif isinstance(result, Exception):
                    self.run_logger.error("Exception raised in worker process '{}'"
                                          "\n>>>\n{}<<<\n".format(
                                              result.worker_name,
                                              result.traceback_exc))
                    raise result
                if results_store is not None:
                    results_store.append(result)
                if results_csv_writer is not None:
                    if result_count == 0 and is_write_header:
                        results_csv_writer.fieldnames = result.keys()
                        results_csv_writer.writeheader()
                    results_csv_writer.writerow(result)
                # self.run_logger.info("Recovered results from worker process '{}'".format(result.worker_name))
                result_count += 1
                # self.info_message("Recovered results from {} of {} worker processes".format(result_count, self.num_processes))
        except (Exception, KeyboardInterrupt) as e:
            for worker in workers:
                worker.terminate()
            raise
        self.run_logger.info("All {} worker processes terminated".format(self.num_processes))
        return results_store

