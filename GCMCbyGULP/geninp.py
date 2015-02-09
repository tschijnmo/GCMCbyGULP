"""
The input generator
===================

"""


import argparse

from .simultask import gen_simul_task_from_YAML
from .gulpinter import COMPUTE_PARAM_FUNCS, GET_RES_FUNCS


def gen_main():
    """The main function for the input generator
    """

    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INP', nargs=1,
                        type=argparse.FileType('r'),
                        help='The name of the JSON/YAML input file')
    args = parser.parse_args()

    # Generate the simulation task object from the input file.
    simul_task = gen_simul_task_from_YAML(
        args.input[0], COMPUTE_PARAM_FUNCS, GET_RES_FUNCS
        )

    # Generate the input files for the simulation tasks.
    simul_task.gen_inp()

    return 0
