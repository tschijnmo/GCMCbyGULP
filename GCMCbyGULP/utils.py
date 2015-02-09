"""
Small utility functions
=======================

"""


import collections


def ensure_list_of_str(val, tag):
    """Ensures that the given value is a list of strings

    It ensures that the given value is a list of strings and return them, or
    value error will be raised. If a single string is given, a singleton list
    will be returned.

    :param val: The value to be ensured to be a list of strings.
    :param str tag: A tag for the value, used for error reporting.
    :returns: The ensured list of strings.
    """

    try:
        if isinstance(val, basestring):
            return [str(val), ]
        elif isinstance(val, collections.Iterable):
            ret_val = []
            for i in val:
                if isinstance(i, basestring):
                    ret_val.append(str(i))
                else:
                    raise ValueError(i)
                continue
            return ret_val
        else:
            raise ValueError(val)
    except ValueError as exc:
        raise ValueError(
            'Invalid value {val} for tag {tag}, string expected!'.format(
                tag=tag, val=exc.args[0]
                )
            )
