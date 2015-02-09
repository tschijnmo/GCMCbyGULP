"""
Computation of gas properties
=============================

This module contains functions to compute the chemical potential and the mass
density of gases from the temperature, pressure and name of the gas.

"""


import math

from scipy.integrate import quad
from CoolProp.CoolProp import PropsSI


#
# Chemical potential
# ------------------
#


def compute_mu_from_TP_gas(temp, press, gas_name):
    """Computes the chemical potential from temperature, pressure for a gas

    :param temp: The temperature in Kelvin
    :param press: The pressure in Pa
    :param str: The name of the gas
    :returns: The chemical potential
    """

    # Compute the fugacity coefficient
    fug_coeff = _compute_fug_coeff(
        lambda p: PropsSI('Z', 'T', temp, 'P', p, gas_name),
        press
        )

    # Compute the molecular weight, a little bit hack, obtained by taking the
    # ratio of the mass density and the molar density.

    call_seq = ('T', temp, 'P', press, gas_name)
    # The call sequence to get properties of the current gas
    mol_wgt = PropsSI('DMASS', *call_seq) / PropsSI('DMOLAR', *call_seq)

    # Return the chemical potential
    return _compute_mu(mol_wgt, fug_coeff, temp, press)


def _compute_fug_coeff(z_func, press):
    """Computes the fugacity coefficient from the compressibility Z

    Computes the fugacity coefficient by the equation

    .. math::

        \\int_0^P \\frac{Z - 1} / P \\mathrm{d}P

    from the compressibility factor.

    :param z_func: The function that returns the compressibility factor for a
        given pressure.
    :param press: The pressure to compute the fugacity coefficient.
    :returns: The fugacity coefficient.
    """
    int_res = quad(
        lambda p: (z_func(p) - 1.0) / p,
        0, press
        )
    return math.exp(int_res[0])


def _compute_mu(mol_wgt, fug_coeff, temp, press):
    """Computes the chemical potential

    This function just uses the statistical mechanics of the ideal gas to get
    the reference thermal wave length for computing the chemical potential,
    then the fugacity coefficient can be applied to the density for correction.

    :param mol_wgt: The molecular weight in kg/mol
    :param fug_coeff: The fugacity coefficient
    :param temp: The temperature
    :param press: The pressure
    """

    # This code can be improved by better scaling of the quantities to values
    # nearer to unity.
    toAA2 = 9.5725E27
    h = 4.135667516E-15
    kB = 8.6173324E-5

    # Thermal wave length
    lamb = math.sqrt(
        (h ** 2 * toAA2) / (2 * math.pi * mol_wgt * 1000.0 * kB * temp)
        )
    rho = ((press * fug_coeff) / (kB * temp)) * 6.241509E-12
    return kB * temp * math.log(lamb ** 3 * rho)


#
# Compute the density
# -------------------
#


def compute_rho_from_TP_gas(temp, press, gas_name):
    """Computes the mass density from temperature, pressure for a gas

    :param temp: The temperature in Kelvin
    :param press: The pressure in Pa
    :param str gas_name: The name of the gas
    :returns: The mass density in kg/m3
    """

    return PropsSI('DMASS', 'T', temp, 'P', press, gas_name)

