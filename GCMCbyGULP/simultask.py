"""
Simulation task
===============

This module contains the class definition for simulation tasks.

"""


import collections
import itertools
import os
import os.path
import subprocess
import math

import numpy as np
import yaml
import pystache

from .dumpres import dump_res_to_file
from .utils import ensure_list_of_str


#
# The exception class
# ===================
#


class InvalidInput(Exception):
    """The exception class for invalid input in the input for simulation tasks

    This is the recommended exception class to use when something is wrong in
    the input.
    """

    pass


#
# The simulation task class definition
# ------------------------------------
#


class SimulTask(object):
    """The simulation task

    Here by a simulation task, we mean we start with a set of fixed or variable
    parameters. For each value of the variable parameters, we can make a
    subdirectory for it, and then we can compute some other parameters based on
    the current value of both the fixed and variable parameters. Then we can
    generate some files by the parameters by using Mustache and run some
    commands. After the simulation, we can still loop over all the
    subdirectories to retrieve the results into tensors for each requested
    results.

    """

    __slots__ = [
        'fixed_params',
        'var_params',
        'compute_params',
        'files',
        'cmds',
        'results',
        'compute_param_funcs',
        'get_res_funcs',
        'params_file',
        '_subdir_fmts',
        ]

    #
    # Initializer
    #

    def __init__(self, inp, compute_param_funcs, get_res_funcs):
        """Initializes a simulation task object

        Basically, in this function, the global options to this package is
        going to be picked out into the attributes of the instance of this
        class. With the others assigned into the parameters dictionaries based
        on their values.

        :param dict inp: The the input for the simulation task, as a
            dictionary.
        :param dict compute_param_funs: The dictionary for compute parameters,
            with the keys being strings for the names of the parameter to
            compute, and values being functions that takes the current
            parameter dictionary to compute the requested parameter.
        :param dict get_res_funcs: The dictionary for getting the results. The
            functions will be given the parameters used for the computation,
            and needs to return the requested results.
        """

        if not isinstance(inp, dict):
            raise InvalidInput(
                'Invalid input, no keys found!'
                )

        # Initialize the lists and dictionaries.
        self.fixed_params = collections.OrderedDict()
        self.var_params = collections.OrderedDict()
        self.compute_params = []
        self.files = []
        self.cmds = []
        self.results = []
        self.params_file = 'params.yaml'

        # Loop over the input
        for k, v in inp.iteritems():

            # First test if the key is a global option.
            if k == 'files':
                self.files = ensure_list_of_str(v, 'files')
            elif k == 'cmds':
                self.cmds = ensure_list_of_str(v, 'cmds')
            elif k == 'compute-params':
                self.compute_params = ensure_list_of_str(
                    v, 'compute-params'
                    )
            elif k == 'results':
                self.results = ensure_list_of_str(v, 'results')
            elif k == 'params-file':
                self.params_file = str(v)

            # If not, test if it is an fixed parameter or a variable one
            elif (isinstance(v, collections.Sequence) and
                  not isinstance(v, basestring)):
                # A sequence heralds a variable parameter.
                # First split the control sequence if present
                splitted_parts = str(k).split(':')
                if len(splitted_parts) == 1:
                    self.var_params[splitted_parts[0]] = v
                elif len(splitted_parts) == 2:
                    param_tag, ctrl_tag = splitted_parts
                    if ctrl_tag == 'arange':
                        vals = np.arange(*v)
                    elif ctrl_tag == 'linspace':
                        vals = np.linspace(*v)
                    else:
                        raise InvalidInput(
                            'Invalid control tag {} in key {}'.format(
                                ctrl_tag, k
                                )
                            )
                    # Convert back to python data types to satisfy PyYAML.
                    self.var_params[param_tag] = [
                        int(i) if isinstance(i, np.int64) else float(i)
                        for i in vals
                        ]
                else:
                    raise InvalidInput(
                        'Invalid key {} with control character'.format(k)
                        )
            elif not isinstance(v, dict):
                # If it is not a dictionary, it means that we are at an
                # atomic node.
                self.fixed_params[k] = v

            else:
                # We do not support complicated nested dictionary here.
                raise InvalidInput(
                    'Unexpected nested dictionary at key {}'.format(k)
                    )

        # Special treatment for the case without variable parameters, we just
        # pop one fixed parameter into the variable one. So that it will still
        # compute the only data point.
        if len(self.var_params) == 0:
            param, val = self.fixed_params.popitem()
            self.var_params[param] = [val, ]

        # Compute the default formats of the subdirectories.
        self._subdir_fmts = self._get_subdir_fmts()

        # The function dictionaries
        self.compute_param_funcs = compute_param_funcs
        self.get_res_funcs = get_res_funcs

    def _get_subdir_fmts(self):
        """Gets the list of default format for subdirectory names
        """

        widths = [
            math.ceil(math.log(len(i), 10))
            for i in self.var_params.itervalues()
            ]
        return [
            '{:0%dd}' % i for i in widths
            ]

    #
    # Input generation
    #

    def gen_inp(self):
        """Generates the input files

        This function will generate the input files for all the combinations of
        the variable parameter values, and generate the input files in each
        subdirectory.
        """

        # Get the contents of the input files to generate, before we start to
        # delve into the subdirectories.
        file_contents = {}
        for i in self.files:
            try:
                f = open(i, 'r')
                file_contents[os.path.basename(i)] = f.read()
                f.close()
            except IOError:
                raise InvalidInput(
                    'Input file {} cannot be read!'.format(i)
                    )
            continue

        # Loop over all possible combinations of the variable parameters
        for params in self._iter_var_params():

            # Make the subdirectory and get into it
            subdir_name = self._idxes_2_subdir(params['idxes'])
            if not os.path.isdir(subdir_name):
                os.makedirs(subdir_name)
            os.chdir(subdir_name)

            # Compute the computed parameters
            for param in self.compute_params:
                try:
                    params[param] = self.compute_param_funcs[param](params)
                except KeyError:
                    raise InvalidInput(
                        'The computed parameter {} is not supported!'.format(
                            param
                            )
                        )
                continue

            # First dump all the parameters that has been computed for this
            # data point.
            try:
                params_file = open(self.params_file, 'w')
                params_file.write(
                    yaml.dump(params)
                    )
                params_file.close()
            except IOError:
                raise InvalidInput(
                    'Unable to write to file {} to dump the data'.format(
                        self.params_file
                        )
                    )

            # Generates the input files.
            for file_name, file_content in file_contents.iteritems():
                res = pystache.render(file_content, params)
                try:
                    out_file = open(file_name, 'w')
                    out_file.write(res)
                    out_file.close()
                except IOError:
                    raise InvalidInput(
                        'Input file {} cannot be written.'.format(file_name)
                        )
                continue

            # Execute the requested commands.
            for cmd in self.cmds:
                cmd_to_exec = pystache.render(cmd, params)
                res = subprocess.call(cmd_to_exec, shell=True)
                if res != 0:
                    raise InvalidInput(
                        'System command {} has failed!'.format(cmd_to_exec)
                        )
                continue

            # Get out of the directory and go to the next combination
            os.chdir(os.pardir)
            continue

    #
    # Results collection
    #

    def dump_res(self, fmt, out_file):
        """Dumps the result of the simulation task into a file

        :param str fmt: The format to dump the output, given as a string.
        :param stream out_file: An output stream-like object to dump the
            results.
        """

        res = self._gather_res()
        dump_res_to_file(out_file, fmt, self.var_params, res)

    def _gather_res(self):
        """Gathers the results of the simulations

        Primarily, this function will collect all the results that is requested
        in the results attribute of the current object by invoking the methods
        in the get_res_funcs dictionary. Then a large numpy tensor is going to
        be created for each of the results that is requested. The variables
        that each of the dimensions correspond to are in the same order as the
        ordered dictionary in the var_params attribute. The values in the
        get_res_funcs  can be a single function, for results with only one real
        number result. Or it can be a pair of an integer and a function, with
        the integer giving the number of real numbers in the result. Then the
        function must return vectors of that length. In the result tensor, the
        components in the result are going to be the last dimensions.

        :returns: A dictionary with the name of the result as the key, and the
            tensor of the results values as the value.
        """

        # Initialize the result to an empty list.
        res = dict()

        # Initialize the results according to the number of data points for
        # each variable parameters.
        var_param_lens = [len(i) for i in self.var_params.itervalues()]
        for i in self.results:
            try:
                get_res_func = self.get_res_funcs[i]
            except KeyError as exc:
                raise InvalidInput(
                    'Requested result {} is not supported!'.format(exc.args[0])
                    )
            if hasattr(get_res_func, '__call__'):
                res[i] = np.empty(var_param_lens)
            else:
                # Use assertion to check the format, since it is not user input
                # but rather plug-in code.
                assert len(get_res_func) == 2
                assert isinstance(get_res_func[0], int)
                assert hasattr(get_res_func[1], '__call__')
                res[i] = np.empty(var_param_lens + [get_res_func[0], ])

        # Loop over all the subdirectories to get the results
        for params in self._iter_var_params():
            idxes = params['idxes']

            # Get into the subdirectory
            subdir = self._idxes_2_subdir(idxes)
            try:
                os.chdir(subdir)
            except OSError:
                raise InvalidInput(
                    'Cannot get results, directory {} does not exist.'.format(
                        subdir
                        )
                    )

            # Read the parameters that was used for computation.
            try:
                param_file = open(self.params_file, 'r')
                params = yaml.load(param_file)
                param_file.close()
            except IOError:
                raise InvalidInput(
                    'The parameter file {} cannot be opened!'.format(
                        self.params_file
                        )
                    )

            # Compute the requested results one-by-one
            for i in self.results:
                get_res_func = self.get_res_funcs[i]
                if not hasattr(get_res_func, '__call__'):
                    get_res_func = get_res_func[1]
                res[i][idxes] = get_res_func(params)
                continue

            # Get back to the parent directory
            os.chdir(os.pardir)

            # Go on to the next directory/simulation
            continue

        # Return the gathered results
        return res

    #
    # Core looping functions
    #

    def _idxes_2_subdir(self, idxes):
        """Converts the indices to subdirectory names

        :param idxes: The indices for each of the variable parameters, in the
            same order as the order in the ordered dictionary. Zero-based.
        :returns: A string which formats the indices.
        """

        return '-'.join(
            fmt.format(idx + 1)
            for fmt, idx in itertools.izip(self._subdir_fmts, idxes)
            )

    def _iter_var_params(self):
        """Gets an iterator that iterates over the variable parameters

        This is the core function for this class. By invoking this function, an
        iterator which will iterates all combinations of the variable
        parameters will be given. The elements that the iterator gives are
        dictionaries with all the fixed and variable parameters, with the
        variables parameters given as a single value. Also an entry named
        ``idxes`` will be added to the dictionary, which correspond to the
        indices in the list of values for each of the variable parameters, in
        the same order as the order in the ordered dictionary var_params.
        """

        def ret_gen():
            ranges = [
                xrange(0, len(i)) for i in self.var_params.itervalues()
                ]
            for idxes in itertools.product(*ranges):
                new_dict = dict(self.fixed_params)
                for param, idx in itertools.izip(self.var_params.iterkeys(),
                                                 idxes):
                    new_dict[param] = self.var_params[param][idx]
                    continue
                new_dict['idxes'] = list(idxes)
                yield new_dict
                continue

        return ret_gen()


#
# Wrapper function
# ================
#


def gen_simul_task_from_YAML(file_stream, compute_param_funcs, get_res_funcs):
    """Generates a simulation task from a YAML file

    This is a shallow wrapper over the constructor for the :py:cls:`SimulTask`
    class. Instead of asking for a dictionary, it can be given as input file
    stream that is going to be parsed as YAML file for the input.

    :param file_stream: The input file stream that is going to be parsed as
        YAML.
    :param dict compute_param_funcs: The dictionary for functions to compute
        parameters.
    :param dict get_res_funcs: The dictionary for functions to get the results.

    """

    try:
        input_dat = yaml.load(file_stream)
    except yaml.YAMLError as exc:
        raise ValueError(
            'Invalid YAML input file:\n {}'.format(exc.args)
            )

    try:
        simul_task = SimulTask(
            input_dat, compute_param_funcs, get_res_funcs
            )
    except InvalidInput as exc:
        raise ValueError(
            'Invalid input file {file_name}: \n{err_msg}'.format(
                file_name=file_stream.name, err_msg=exc.args[0]
                )
            )

    return simul_task

