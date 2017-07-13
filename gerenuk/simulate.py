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

from gerenuk import utility

FSC2_CONFIG_TEMPLATE = """\
//Number of population samples (demes)
2
//Population effective sizes (number of genes)
1000
1000
//Sample sizes
10
10
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event
1  historical event
10 0 1 1 2 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 1 0.00000 0.00000002 0.33
"""

class GerenukSimulate(object):

    def __init__(self, **kwargs):
        self.fsc2_path = kwargs.pop("fsc2_path", "fsc25")
        if "rng" in kwargs:
            self.rng = kwargs.pop("rng")
        else:
            self.rng = random.Random()
        if kwargs:
            raise Exception("Unrecognized keywords arguments: {}".format(kwargs))

    def generate_configuration_file(self, filepath):
        with open(filepath, "w") as dest:
            config = FSC2_CONFIG_TEMPLATE
            dest.write(config)

    def execute(self):
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend("-n", "1") # number of simulations to perform
        cmds.extend("--allsites") # output the whole DNA sequence, incl. monomorphic sites
        # cmds.extend("--inf") # generates DNA mutations according to an infinite site (IS) mutation model

if __name__ == "__main__":
    gs = GerenukSimulate()

