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
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
## BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import subprocess
import random
import sys

##############################################################################
## Process Control/Handling

import locale
import codecs
import sys

try:
    ENCODING = locale.getdefaultlocale()[1]
except ValueError:
    ENCODING = None # let default value be assigned below

if ENCODING == None:
    ENCODING = 'UTF-8'

def _bytes_to_text(s):
    """
    Converts a byte string (as read from, e.g., standard input)
    to a text string.

    In Python 3, this is from type ``bytes`` to ``str``.
    In Python 2, this is, confusingly, from type ``str`` to ``unicode``.

    """
    s = codecs.decode(s, ENCODING)
    if sys.hexversion < 0x03000000:
        s = codecs.encode(s, "utf-8")
    return s

def _communicate_process(p, commands=None, timeout=None):
    if isinstance(commands, list) or isinstance(commands, tuple):
        commands = "\n".join(str(c) for c in commands)
    if commands is not None:
        commands = str.encode(commands)
    if timeout is None:
        stdout, stderr = p.communicate(commands)
    else:
        try:
            stdout, stderr = p.communicate(commands, timeout=timeout)
        except TypeError as e:
            if "unexpected keyword argument 'timeout'" in str(e):
                stdout, stderr = p.communicate(commands)
            else:
                raise
    if stdout is not None:
        stdout = _bytes_to_text(stdout)
    if stderr is not None:
        stderr = _bytes_to_text(stderr)
    return stdout, stderr

## Process Control/Handling
##############################################################################


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

    def execute(self):
        cmds = []
        cmds.append(self.fsc2_path)
        cmds.extend("-n", "1") # number of simulations to perform
        cmds.extend("--allsites") # output the whole DNA sequence, incl. monomorphic sites
        # cmds.extend("--inf") # generates DNA mutations according to an infinite site (IS) mutation model

if __name__ == "__main__":
    gs = GerenukSimulate()

