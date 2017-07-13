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

import locale
import codecs
import sys
import logging
import tempfile

##############################################################################
## Process Control/Handling

try:
    ENCODING = locale.getdefaultlocale()[1]
except ValueError:
    ENCODING = None # let default value be assigned below

if ENCODING == None:
    ENCODING = 'UTF-8'

def bytes_to_text(s):
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

def communicate_process(p, commands=None, timeout=None):
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
        stdout = bytes_to_text(stdout)
    if stderr is not None:
        stderr = bytes_to_text(stderr)
    return stdout, stderr

##############################################################################
## Logging

_LOGGING_LEVEL_ENVAR = "GERENUK_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "GERENUK_LOGGING_FORMAT"

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        self.handlers = []
        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)
        if kwargs.get("log_to_file", True):
            if "log_stream" in kwargs:
                log_stream = kwargs.get("log_stream")
            else:
                log_stream = open(kwargs.get("log_path", self.name + ".log"), "w")
            handler2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            handler2.setLevel(file_logging_level)
            handler2.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler2)
            self.handlers.append(handler2)
        self._system = None

    def _get_system(self):
        return self._system

    def _set_system(self, system):
        self._system = system
        if self._system is None:
            for handler in self.handlers:
                handler.setFormatter(self.get_default_formatter())
        else:
            for handler in self.handlers:
                handler.setFormatter(self.get_simulation_generation_formatter())

    system = property(_get_system, _set_system)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simulation_generation_formatter(self):
        # f = logging.Formatter("[%(asctime)s] t = %(elapsed_time)10.6f: %(message)s")
        f = logging.Formatter("[%(asctime)s] %(simulation_time)s%(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def supplemental_info_d(self):
        if self._system is None or self._system.elapsed_time == 0 :
            return {
                    "simulation_time" : "Setup: ",
                    }
        else:
            return {
                    "simulation_time" : "[t = {:13.6f}] ".format(self._system.elapsed_time),
                    }

    def debug(self, msg):
        self._log.debug("[DEBUG] {}".format(msg), extra=self.supplemental_info_d())

    def info(self, msg):
        self._log.info(msg, extra=self.supplemental_info_d())

    def warning(self, msg):
        self._log.warning(msg, extra=self.supplemental_info_d())

    def error(self, msg):
        self._log.error(msg, extra=self.supplemental_info_d())

    def critical(self, msg):
        self._log.critical(msg, extra=self.supplemental_info_d())

    def flush(self):
        for handler in self.handlers:
            handler.flush()
