"""
Interface to the GULP program
=============================

This module defines common parameters that needs to be computed in the constant
dictionary COMPUTE_PARAM_FUNCS, and the common results of interest in GCMC
simulation in the GET_RES_FUNCS dictionary.

"""


import re
import collections

import periodic
import numpy as np
from numpy import linalg

from .simultask import InvalidInput
from .gasprops import compute_mu_from_TP_gas, compute_rho_from_TP_gas
from .utils import ensure_list_of_str


#
# Compute parameters
# ------------------
#

def compute_mu(params):
    """Computes the chemical potential

    This function is going to compute the chemical potential based on the
    temperature (T), pressure (P), and the name of the gas (gas-name), by using
    the CoolProp.

    :param params: The dictionary of the parameters.
    :returns: The chemical potential of the gas.
    """

    try:
        temp = params['T']
        press = params['P']
        gas_name = params['gas-name']
    except KeyError as exc:
        raise InvalidInput(
            'Parameter {} is not defined to compute chemical potential'.format(
                exc.args[0]
                )
            )

    return compute_mu_from_TP_gas(temp, press, gas_name)


COMPUTE_PARAM_FUNCS = {
    'mu': compute_mu,
    }


#
# Get results
# -----------
#

def get_uptake(params):
    """Gets the gravimetric uptake

    This is the wrapper function for getting the uptake of the data point. The
    result is of length four, giving the excess gravimetric uptake, excess
    volumetric uptake, raw gravimetric and volumetric uptake.

    The gravimetric uptakes are in unit of percentage, and the volumetric
    uptake is in the unit of g/L.
    """

    # Read the output file
    lines = _get_output_lines(params)

    # Get the initial list of atomic symbols
    init_symbs = _get_init_symbs(lines)

    # Get the weight of the substrate
    try:
        ads_atms = ensure_list_of_str(
            params['ads-atoms'], 'ads-atoms'
            )
    except KeyError:
        raise InvalidInput(
            'The adsorbant atoms unset in ads-atoms tag to compute uptake.'
            )
    subs_wgt, subs_n_atm = _get_subs_weight(init_symbs, ads_atms)

    # Get the weight of the adsorbant
    ads_wgt = (
        (_get_mean_n_atms(lines) - subs_n_atm) * _get_mean_atm_wgt(ads_atms)
        )

    # Get the volume of the unit cell
    vol = _get_volume(lines)

    # Get the free weight of the adsorbant inside the volume of the simulation
    # cell.
    try:
        temp = params['T']
        press = params['P']
        gas_name = params['gas-name']
    except KeyError as exc:
        raise InvalidInput(
            'Property {} not set in the parameters to compute uptake'.format(
                exc.args[0]
                )
            )
    density = compute_rho_from_TP_gas(temp, press, gas_name)
    ads_free_wgt = (
        density * vol *  # This is in kg m-3 Ang3
        6.022141E-4  # Wolfram alpha: 1 kg/m3 * 1 cubic angstrom in u
        )
    exc_ads_wgt = ads_wgt - ads_free_wgt

    # Compute the raw and excess gravimetric uptake.
    raw_gr_uptake = ads_wgt / subs_wgt * 100.0
    exc_gr_uptake = exc_ads_wgt / subs_wgt * 100.0

    # Compute the raw and excess volumetric uptake.
    vol_uptake_factor = 1660.539
    # Wolfram alpha: 1 u / cubic angstrom in (g / L)
    raw_vol_uptake = ads_wgt / vol * vol_uptake_factor
    exc_vol_uptake = exc_ads_wgt / vol * vol_uptake_factor

    # Return the four-component vector
    return exc_gr_uptake, exc_vol_uptake, raw_gr_uptake, raw_vol_uptake


GET_RES_FUNCS = {
    'uptake': (4, get_uptake),
    }


#
# Internal functions
# ^^^^^^^^^^^^^^^^^^
#


def _get_output_lines(params):
    """Reads the lines in the output file

    The name of the output file are going to be attempted to get from the
    parameters under the key ``output``. The lines in the output file are going
    to be returned as a list of strings.
    """

    try:
        output_file = open(params['output'], 'r')
        lines = output_file.readlines()
        output_file.close()
    except KeyError:
        raise InvalidInput(
            'The output file is not set to compute the result'
            )
    except IOError:
        raise InvalidInput(
            'The output file {} cannot be read.'.format(params['output'])
            )

    return lines


def _get_init_symbs(lines):
    """Gets the symbols of the initial set of atoms

    :param lines: The lines in the output file
    :returns: A list of atomic symbols in the initial configuration of the
        simulation.
    """

    init_coord_beg_sentinel = re.compile(
        'Fractional coordinates of asymmetric unit'
        )
    init_coord_end_sentinel = re.compile(
        'Molecule list generated from bond lengths'
        )
    coord_re = re.compile(r'\d+\s+(\w+)\s+[cs]\s+\d+\.\d+\s+\d+\.\d+')

    beg_sentinel_idx = next(i for i, v in enumerate(lines)
                            if init_coord_beg_sentinel.search(v) != None)
    end_sentinel_idx = next(i for i, v in enumerate(lines)
                            if init_coord_end_sentinel.search(v) != None)

    coord_search = [
        coord_re.search(l) for l in lines[beg_sentinel_idx:end_sentinel_idx]
        ]
    symbols = [i.group(1) for i in coord_search if i is not None]

    return symbols


def _get_subs_weight(symbs, ads_atms):
    """Gets the weight of the substrate

    :param symbs: A list of symbols for the system
    :param ads_atms: The list of atomic symbols for the adsorbant.
    :returns: The weight of the substrate and the number of atoms
    """

    subs = [i for i in symbs if i not in ads_atms]
    cnt = collections.Counter(subs)

    wgt = 0.0
    for i, v in cnt.iteritems():
        wgt += _get_atm_wgt(i) * v
        continue

    return wgt, len(subs)


def _get_atm_wgt(symb):
    """Gets the weight of a given element symbol

    :param symb: The element symbol, which is able to be suffixed by some
        numerals.
    :returns: The atomic weight in AU.
    """

    elem_symb = re.search(r'([a-zA-Z]+)\d*', symb).group(1)
    return periodic.element(elem_symb).mass


def _get_mean_n_atms(lines):
    """Gets the mean number of atoms

    :param lines: The lines in the output file
    :returns: The final mean number of atoms
    """

    trials_re = re.compile('Trials:.*')
    match_res = [trials_re.search(i) for i in lines]
    last_trial = [i for i in match_res if i is not None][-1]

    return float(
        last_trial.group(0).split()[-1]
        )


def _get_mean_atm_wgt(atm_symbs):
    """Gets the mean atomic weight of a list of symbols

    :param atm_symbs: An sequence for a list of atomic symbols.
    :returns: The mean atomic weight in AU.
    """

    return sum(_get_atm_wgt(i) for i in atm_symbs) / len(atm_symbs)


def _get_volume(lines):
    """Gets the volume of the unit cell

    :param lines: The lines in the output file
    :returns: The volume of the unit cell in cubic Angstroms.
    """

    begin_sentinel = re.compile(r'^ *Cartesian lattice vectors \(Angstroms\)')
    end_sentinel = re.compile(r'^ *Cell parameters')
    lattice_re = re.compile(r'(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')

    begin_idx = next(i for i, v in enumerate(lines)
                     if begin_sentinel.search(v) != None)
    end_idx = next(i for i, v in enumerate(lines)
                   if end_sentinel.search(v) != None)

    lattice_search = [lattice_re.search(v) for v in lines[begin_idx:end_idx]]
    lattices = [
        i.group(1, 2, 3) for i in lattice_search if i is not None
        ]

    lattice_matrix = np.array(lattices, dtype=np.float64)

    return linalg.det(lattice_matrix)
