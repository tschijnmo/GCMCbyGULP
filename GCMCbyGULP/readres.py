"""
Result reader for GCMC simulations
==================================

"""


import argparse
import os.path

from .simultask import gen_simul_task_from_YAML
from .gulpinter import COMPUTE_PARAM_FUNCS, GET_RES_FUNCS


def read_main():
    """The main function for the read script"""

    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INP', nargs=1,
                        type=argparse.FileType('r'),
                        help='The name of the JSON/YAML input file'
                        )
    parser.add_argument('-o', '--out', metavar='FILE', nargs=1,
                        type=argparse.FileType('w'),
                        default=argparse.SUPPRESS,
                        help='The output file for dumping the result'
                        )
    parser.add_argument('-f', '--format', metavar='FORMAT', nargs=1,
                        type=str, choices=[
                            'txt', 'mat',
                            ],
                        default='txt'
                        )
    args = parser.parse_args()

    # Attempt to open the output file
    if 'out' in args:
        out_file = args.out
    else:
        out_file_name = args.input.name.split('.')[0] + '.out',
        if os.path.isfile(out_file_name):
            raise ValueError(
                'The output file {} already exists.'.format(out_file_name) +
                '\nTo overwrite, please use the -o option.'
                )
        else:
            out_file = open(out_file_name, 'w')

    # Generate the simulation task object
    simul_task = gen_simul_task_from_YAML(
        args.input[0], COMPUTE_PARAM_FUNCS, GET_RES_FUNCS
        )

    # Gather the results and output to a file
    simul_task.dump_res(args.format, out_file)

    return 0
