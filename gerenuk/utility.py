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

import shutil
import errno
import collections
import csv
import locale
import codecs
import sys
import os
import logging
import tempfile
import re

##############################################################################
## StringIO
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3

##############################################################################
## Opening files

def pre_py34_open(file,
        mode='r',
        buffering=-1,
        encoding=None,
        errors=None,
        newline=None,
        closefd=True,
        opener=None):
    if encoding is not None:
        raise NotImplementedError
    if errors is not None:
        raise NotImplementedError
    if newline is None:
        if mode.startswith("r"):
            mode = mode + "U"
    else:
        raise NotImplementedError
    if closefd is not True:
        raise NotImplementedError
    if opener is not None:
        raise NotImplementedError
    return open(
            file,
            mode=mode,
            buffering=buffering)

##############################################################################
## CSV File Handling

def open_destput_file_for_csv_writer(filepath, is_append=False):
    if filepath is None or filepath == "-":
        dest = sys.stdout
    elif sys.version_info >= (3,0,0):
        dest = open(filepath, "a" if is_append else "w", newline='')
    else:
        dest = open(filepath, "ab" if is_append else "wb")
    return dest

def get_csv_writer(
        dest,
        fieldnames=None,
        delimiter=",",
        ):
    if isinstance(dest, str):
        dest = open_destput_file_for_csv_writer(filepath=dest)
    writer = csv.DictWriter(
            dest,
            fieldnames=fieldnames,
            restval="NA",
            delimiter=delimiter,
            lineterminator=os.linesep,
            )
    return writer

def write_dict_csv(list_of_dicts, filepath, fieldnames=None, is_no_header_row=False, is_append=False):
    dest = open_destput_file_for_csv_writer(
            filepath=filepath,
            is_append=is_append)
    if fieldnames is None:
        fieldnames = list_of_dicts[0].keys()
    with dest:
        writer = get_csv_writer(
                dest=dest,
                delimiter=",",
                lineterminator=os.linesep,
                )
        writer.fieldnames = fieldnames
        if not is_no_header_row:
            writer.writeheader()
        writer.writerows(list_of_dicts)

##############################################################################
## Configuration File Handling

def parse_legacy_configuration(filepath, config_d=None):
    recognized_preamble_keys = collections.OrderedDict([
            ("concentrationshape", float),
            ("concentrationscale", float),
            ("thetashape", float),
            ("thetascale", float),
            ("ancestralthetashape", float),
            ("ancestralthetascale", float),
            ("thetaparameters", str),
            ("taushape", float),
            ("tauscale", float),
            ("timeinsubspersite", float),
            ("bottleproportionshapea", float),
            ("bottleproportionshapeb", float),
            ("bottleproportionshared", float),
            ("migrationshape", float),
            ("migrationscale", float),
            ("numtauclasses", float),
            ])
    sample_table_keys = [
            ("taxon_label", str),
            ("locus_label", str),
            ("ploidy_factor", float),
            ("mutation_rate_factor", float),
            ("num_genes_deme0", int),
            ("num_genes_deme1", int),
            ("ti_tv_rate_ratio", float),
            ("num_sites", int),
            ("freq_a", float),
            ("freq_c", float),
            ("freq_g", float),
            ("alignment_filepath", str),
            ]
    sample_table_begin_pattern = re.compile(r"^\s*BEGIN\s+SAMPLE_TBL\s*$", re.I)
    sample_table_end_pattern = re.compile(r"^\s*END\s+SAMPLE_TBL\s*$", re.I)
    sample_table_splitter = re.compile("\s+", re.I)
    if config_d is None:
        config_d = CaseInsensitiveDict()
    config_d["params"] = {}
    config_d["locus_info"] = []
    section = "preamble"
    src = open(filepath)
    for row_idx, row in enumerate(src):
        row = row.strip()
        if not row:
            continue
        comment_start_idx = row.find("#")
        if section == "sample-table" and comment_start_idx > 0:
            raise ValueError("Configuration file '{}', row {}: sample table section cannot contain mid-row comment".format(
                filepath,
                row_idx+1))
        if comment_start_idx > -1:
            row = row[0:comment_start_idx]
        row = row.strip()
        if not row:
            continue
        if section == "preamble":
            if sample_table_begin_pattern.match(row):
                section = "sample-table"
                continue
            row_parts = row.split("=")
            if len(row_parts) != 2:
                raise ValueError("Configuration file '{}', row {}: sample table section cannot contain mid-row comment".format(
                    filepath,
                    row_idx+1))
            key = row_parts[0].strip()
            case_normalized_key = key.lower()
            if case_normalized_key not in recognized_preamble_keys:
                raise ValueError("Configuration file '{}', row {}: unrecognized preamble key '{}'".format(
                    filepath,
                    row_idx+1,
                    key
                    ))
            config_d["params"][key] = recognized_preamble_keys[case_normalized_key](row_parts[1].strip())
        else:
            if sample_table_end_pattern.match(row):
                continue
            cols = sample_table_splitter.split(row)
            locus_info = {}
            if len(cols) != len(sample_table_keys):
                raise ValueError("Configuration file '{}', row {}: expecting {} columns but only found {}".format(
                    filepath,
                    row_idx+1,
                    len(sample_table_keys),
                    len(cols),
                    ))
            for (key, val_type), val in zip(sample_table_keys, cols):
                locus_info[key] = val_type(val)
            config_d["locus_info"].append(locus_info)
    return config_d

##############################################################################
## Temporary Directory Handling ( [somewhat] replicates 3.2 'TemporaryDirectory')

class TemporaryDirectory(object):

    def __init__(self,
            suffix="",
            prefix="",
            parent_dir=None,
            ):
        self.suffix = suffix
        self.prefix = prefix
        self.parent_dir = parent_dir
        self.temp_dir_path = None
        self.suppress_cleanup = False

    def cleanup(self):
        try:
            shutil.rmtree(self.temp_dir_path)
        except OSError as e:
            # Reraise unless ENOENT: No such file or directory
            # (ok if directory has already been deleted)
            if e.errno != errno.ENOENT:
                raise

    def __enter__(self):
        self.temp_dir_path = tempfile.mkdtemp(
                suffix=self.suffix,
                prefix=self.prefix,
                dir=self.parent_dir)
        return self.temp_dir_path

    def __exit__(self, *args):
        if not self.suppress_cleanup:
            self.cleanup()

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
## Command line processing

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*?)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError("Cannot parse fieldname and label (format required: fieldname:value): {}".format(label))
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map

##############################################################################
## Logging

_LOGGING_LEVEL_ENVAR = "GERENUK_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "GERENUK_LOGGING_FORMAT"

class RunLogger(object):

    NOTSET_MESSAGING_LEVEL = logging.NOTSET
    DEBUG_MESSAGING_LEVEL = logging.DEBUG
    INFO_MESSAGING_LEVEL = logging.INFO
    WARNING_MESSAGING_LEVEL = logging.WARNING
    ERROR_MESSAGING_LEVEL = logging.ERROR
    CRITICAL_MESSAGING_LEVEL = logging.CRITICAL

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(RunLogger.DEBUG_MESSAGING_LEVEL)
        self.handlers = []
        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", RunLogger.INFO_MESSAGING_LEVEL))
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
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", RunLogger.DEBUG_MESSAGING_LEVEL))
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
        if level in [
                RunLogger.NOTSET_MESSAGING_LEVEL,
                RunLogger.DEBUG_MESSAGING_LEVEL,
                RunLogger.INFO_MESSAGING_LEVEL,
                RunLogger.WARNING_MESSAGING_LEVEL,
                RunLogger.ERROR_MESSAGING_LEVEL,
                RunLogger.CRITICAL_MESSAGING_LEVEL,
                ]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = RunLogger.NOTSET_MESSAGING_LEVEL
        elif level_name == "DEBUG":
            level = RunLogger.DEBUG_MESSAGING_LEVEL
        elif level_name == "INFO":
            level = RunLogger.INFO_MESSAGING_LEVEL
        elif level_name == "WARNING":
            level = RunLogger.WARNING_MESSAGING_LEVEL
        elif level_name == "ERROR":
            level = RunLogger.ERROR_MESSAGING_LEVEL
        elif level_name == "CRITICAL":
            level = RunLogger.CRITICAL_MESSAGING_LEVEL
        else:
            level = RunLogger.NOTSET_MESSAGING_LEVEL
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simulation_generation_formatter(self):
        # f = logging.Formatter("[%(asctime)s] %(simulation_time)s%(message)s")
        f = logging.Formatter("[%(asctime)s] %(message)s")
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

    def debug(self, msg):
        self._log.debug("[DEBUG] {}".format(msg))

    def info(self, msg):
        self._log.info(msg)

    def warning(self, msg):
        self._log.warning(msg)

    def error(self, msg):
        self._log.error(msg)

    def critical(self, msg):
        self._log.critical(msg)

    def log(self, msg, level):
        if level == RunLogger.DEBUG_MESSAGING_LEVEL:
            self.debug(msg)
        elif level == RunLogger.INFO_MESSAGING_LEVEL:
            self.info(msg)
        elif level == RunLogger.WARNING_MESSAGING_LEVEL:
            self.warning(msg)
        elif level == RunLogger.ERROR_MESSAGING_LEVEL:
            self.error(msg)
        elif level == RunLogger.CRITICAL_MESSAGING_LEVEL:
            self.critical(msg)
        else:
            raise Exception("Unrecognized messaging level: {}".format(level))

    def flush(self):
        for handler in self.handlers:
            handler.flush()

###############################################################################
# CaseInsensitiveDict
#
# From:
#        https://github.com/kennethreitz/requests
#
# Copyright 2014 Kenneth Reitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

class CaseInsensitiveDict(collections.MutableMapping):
    """
    A case-insensitive ``dict``-like object.

    Implements all methods and operations of
    ``collections.MutableMapping`` as well as dict's ``copy``. Also
    provides ``lower_items``.

    All keys are expected to be strings. The structure remembers the
    case of the last key to be set, and ``iter(instance)``,
    ``keys()``, ``items()``, ``iterkeys()``, and ``iteritems()``
    will contain case-sensitive keys. However, querying and contains
    testing is case insensitive:

        cid = CaseInsensitiveDict()
        cid['Accept'] = 'application/json'
        cid['aCCEPT'] == 'application/json'  # True
        list(cid) == ['Accept']  # True

    For example, ``headers['content-encoding']`` will return the
    value of a ``'Content-Encoding'`` response header, regardless
    of how the header name was originally stored.

    If the constructor, ``.update``, or equality comparison
    operations are given keys that have equal ``.lower()``s, the
    behavior is undefined.

    """
    def __init__(self, data=None, **kwargs):
        self._store = dict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def __setitem__(self, key, value):
        # Use the lowercased key for lookups, but store the actual
        # key alongside the value.
        self._store[key.lower()] = (key, value)

    def __getitem__(self, key):
        return self._store[key.lower()][1]

    def __delitem__(self, key):
        del self._store[key.lower()]

    def __iter__(self):
        return (casedkey for casedkey, mappedvalue in self._store.values())

    def __len__(self):
        return len(self._store)

    def lower_items(self):
        """Like iteritems(), but with all lowercase keys."""
        return (
            (lowerkey, keyval[1])
            for (lowerkey, keyval)
            in self._store.items()
        )

    def __eq__(self, other):
        if isinstance(other, collections.Mapping):
            other = CaseInsensitiveDict(other)
        else:
            return NotImplementedError
        # Compare insensitively
        return dict(self.lower_items()) == dict(other.lower_items())

    # Copy is required
    def copy(self):
        return CaseInsensitiveDict(self._store.values())

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, dict(self.items()))

# CaseInsensitiveDict
###############################################################################


###############################################################################
## OrderedCaselessDict
##
## From DendroPy (http://dendropy.org), under the BSD License:
##
##      Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##      All rights reserved.
##

class OrderedCaselessDict(dict):
    """
    Inherits from dict. Maintains two sets of keys: the first the keys
    belonging to dict, which actually accesses the container
    items. This is always cast to lower() whenever it is called, thus
    ensuring that keys are always of the same case. The second set of
    keys it maintains locally in an list, thus maintaining the order
    in which they were added. The second set of keys is not cast to
    lower(), which means that client code can always recover the
    original 'canonical' casing of the keys.

    ONLY TAKES STRING KEYS!
    """

    def __init__(self, other=None):
        """
        __init__ creates the local set of keys, and then initializes self with
        arguments, if any, by using the superclass methods, keeping
        the ordered keys in sync.
        """
        super(OrderedCaselessDict, self).__init__()
        self._ordered_keys = []
        if other is not None:
            if isinstance(other, dict):
                for key, val in other.items():
                    if key.lower() not in self:
                        self._ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)
            else:
                for key, val in other:
                    if key.lower() not in self:
                        self._ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)

    def __deepcopy__(self, memo):
        o = self.__class__()
        memo[id(self)] = o
        for key, val in self.items():
            o[key] = copy.deepcopy(val, memo)
        return o

    def copy(self):
        "Returns a shallow copy of self."
        return self.__class__(self)

    def iterkeys(self):
        "Returns an iterator over self's ordered keys."
        return iter(self._ordered_keys)

    def itervalues(self):
        "Returns an iterator over self's key, value pairs."
        for key in self.iterkeys():
            yield self[key.lower()]

    def iteritems(self):
        "Returns an iterator over self's values."
        for key in self.iterkeys():
            yield (key, self[key.lower()])

    def items(self):
        "Returns key, value pairs in key-order."
        return [(key, self[key]) for key in self.iterkeys()]

    def values(self):
        "Returns list of key, value pairs."
        return [v for v in self.itervalues()]

    def __iter__(self):
        "Returns an iterator over self's ordered keys."
        return self.iterkeys()

    def __getitem__(self, key):
        "Gets an item using a case-insensitive key."
        return super(OrderedCaselessDict, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        "Sets an item using a case-insensitive key,"
        if key.lower() not in self:
            self._ordered_keys.append(str(key))
        super(OrderedCaselessDict, self).__setitem__(key.lower(), value)

    def __delitem__(self, key):
        "Remove item with specified key."
        del(self._ordered_keys[self.index(key)])
        super(OrderedCaselessDict, \
              self).__delitem__(key.lower())

    def __contains__(self, key):
        "Returns true if has key, regardless of case."
        return super(OrderedCaselessDict, self).__contains__(key.lower())

    def pop(self, key, alt_val=None):
        "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
        if key.lower() in self:
            val = self[key]
            self.__delitem__(key.lower())
            return val
        else:
            return alt_val

    def popitem(self):
        "a.popitem()  remove and last (key, value) pair"
        key = self._ordered_keys[-1]
        item = (key, self[key.lower()])
        self.__delitem__(key)
        return item

    def caseless_keys(self):
        "Returns a copy of the ordered list of keys."
        return [k.lower() for k in self._ordered_keys]

    def index(self, key):
        """
        Return the index of (caseless) key.
        Raise KeyError if not found.
        """
        count = 0
        for k in self._ordered_keys:
            if k.lower() == key.lower():
                return count
            count = count + 1
        raise KeyError(key)

    def keys(self):
        "Returns a copy of the ordered list of keys."
        return list(self._ordered_keys)

    def clear(self):
        "Deletes all items from the dictionary."
        self._ordered_keys = []
        super(OrderedCaselessDict, self).clear()

    def has_key(self, key):
        "Returns true if has key, regardless of case."
        return key.lower() in self

    def get(self, key, def_val=None):
        "Gets an item by its key, returning default if key not present."
        return super(OrderedCaselessDict, self).get(key.lower(), def_val)

    def setdefault(self, key, def_val=None):
        "Sets the default value to return if key not present."
        return super(OrderedCaselessDict, self).setdefault(key.lower(), def_val)

    def update(self, other):
        """
        updates (and overwrites) key/value pairs:
        k = { 'a':'A', 'b':'B', 'c':'C'}
        q = { 'c':'C', 'd':'D', 'f':'F'}
        k.update(q)
        {'a': 'A', 'c': 'C', 'b': 'B', 'd': 'D', 'f': 'F'}
        """
        for key, val in other.items():
            if key.lower() not in self:
                self._ordered_keys.append(str(key))
            super(OrderedCaselessDict, self).__setitem__(key.lower(), val)

    def fromkeys(self, iterable, value=None):
        "Creates a new dictionary with keys from seq and values set to value."
        ocd = OrderedCaselessDict()
        for key in iterable:
            if key.lower() not in self:
                self[key] = value
        return ocd
