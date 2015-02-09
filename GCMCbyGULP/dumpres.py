"""
Result dumping
==============

This module contains function to dump the simulation task results according to
various formats.

"""


from __future__ import print_function

import itertools


#
# The main driver function
# ------------------------
#


def dump_res_to_file(out_file, fmt, var_params, res):
    """Dumps the result to the output file

    :param out_file: The output stream to dump the result.
    :param str fmt: The format of the dumping.
    :param var_params: An ordered dictionary for the variable parameters, with
        the tag for the parameters in the key and the list of possible values
        in the value.
    :param res: The dictionary for the results, with tag for the results as key
        and the tensor for the results as value.
    :returns: None
    """

    if fmt == 'txt':
        _txt_dump(out_file, var_params, res)
    elif fmt == 'mat':
        _mat_dump(out_file, var_params, res)
    else:
        raise ValueError(
            'Unsupported output format {}'.format(fmt)
            )

    return None


#
# Engines
# -------
#


#
# Text format
# ^^^^^^^^^^^
#


def _txt_dump(out_file, var_params, res,
              comments='#', delimiter=' ', fmt='{!r}'):
    """Dumps the result into a text file

    The format will be one line for each data-point, with first columns giving
    the values of the variable parameters and the rest giving the results for
    each requested results. At the beginning of the file, a comment is going to
    be written for the actual format. The result will be easily parsed by the
    numpy loadtxt function.

    :param out_file: The output stream
    :param var_params: The ordered dictionary for the variable parameters
    :param res: The dictionary for the results.
    :param optional comments: The Character to lead the comment line.
    :param optional delimiter: The delimiter for the data fields.
    :param optional fmt: The format to be used for each scalar
    """

    # The indices range for each variable parameter
    idx_ranges = [range(0, len(i)) for i in var_params.itervalues()]

    # Keep a list of the key and stick to this order.
    res_keys = [i for i in res.iterkeys()]

    # The number of values for each keys
    n_vals_4_keys = []
    for i in res_keys:
        shape = res[i].shape
        if len(shape) == len(idx_ranges):
            # It means that it is a scalar result
            n_vals_4_keys.append(0)  # Use 0 for scalar results
        elif len(shape) - 1 == len(idx_ranges):
            # Vector results
            n_vals_4_keys.append(
                shape[-1]
                )
        else:
            assert False

    # First output the title comment for field names
    field_titles = list(var_params.iterkeys())
    for key, n_val in itertools.izip(res_keys, n_vals_4_keys):
        if n_val == 0:
            field_titles.append(key)
        else:
            field_titles.extend(
                key + '[{}]'.format(i)
                for i in xrange(n_val)
                )
        continue
    # Form and dump the title
    title = ''.join((
        comments, delimiter, delimiter.join(field_titles)
        ))
    print(title, file=out_file)

    # Write the results
    for idxes in itertools.product(*idx_ranges):
        # First add the parameter values
        vals = [
            vals[idx]
            for vals, idx in itertools.izip(var_params.itervalues(), idxes)
            ]
        # Add the results values
        for key, n_val in itertools.izip(res_keys, n_vals_4_keys):
            if n_val == 0:
                vals.append(res[key][idxes])
            else:
                vals.extend(
                    res[key][idxes + (i, )] for i in xrange(0, n_val)
                    )

        # Format the values
        data_line = delimiter.join(
            fmt.format(i) for i in vals
            )
        # Dump the data line
        print(data_line, file=out_file)

        continue

    # Return
    return None


#
# The Matlab MAT format
# ^^^^^^^^^^^^^^^^^^^^^
#


def _mat_dump(out_file, var_params, res):
    """Dump the result into a MATLAB MAT file

    This function is not implemented yet.
    """

    raise NotImplementedError()
