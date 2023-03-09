#
# SOMD is an ab-initio molecular dynamics package designed for the SIESTA code.
# Copyright (C) 2023 github.com/initqp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

"""
Physical constants and default values in SOMD.
"""

__all__ = ['CONSTANTS', 'SOMDDEFAULTS']


class CONSTANTS:
    """
    CODATA 2018 standard constants [1].

    Attributes
    ----------
    HARTREE : float
        Hartree energy. In unit of (kJ/mol).
    BOLTZCONST : float
        Boltzmann Constant. In unit (kJ/mol/K).
    AVOGACONST : float
        Avogadro Number. In unit of (1).
    PLANKCONST : float
        Planck constant. In unit of (kJ/mol*ps).
    ELECTCONST : float
        Elementary charge. In unit of (C).
    BOHRRADIUS : float
        Bohr radius. In unit of (nm).
    LIGHTSPEED : float
        Speed of light in vacuum. In unit of (nm/ps).
    MEGAPASCAL : float
        Pressure unit. In unit of (kJ/mol/nm^3).

    References
    ----------
    .. [1] https://physics.nist.gov/cuu/Constants/index.html
    """
    HARTREE: float = 2625.4978694234683
    BOLTZCONST: float = 8.31446261815324E-3
    AVOGACONST: float = 6.0221367E23
    PLANKCONST: float = 0.399031002270895
    ELECTCONST: float = 1.602176634E-19
    BOHRRADIUS: float = 5.29177210903E-2
    LIGHTSPEED: float = 299792.458
    MEGAPASCAL: float = 0.60221367


class SOMDDEFAULTS:
    """
    Default values of SOMD internal variables.

    Attributes
    ----------
    POTLIST : list
        Default potentials to update.
    NHCLENGTH : int
        Length of the Nose-Hoover chains.
    NHCNRESPA : int
        Number of the RESPA loops of the Nose-Hoover chains.
    GEODESICKR : int
        The K_r value of the geodesic Langevin integrators.
    LATTICETOL : float
        Tolerance for determining lattice parameters.
    SIMUUPDATE : bool
        If simultaneously update potentials.
    SIESTATIMEOUT : int
        Timeout for any SIESTA operations (second).
    """
    POTLIST: list = None
    NHCLENGTH: int = 6
    NHCNRESPA: int = 4
    GEODESICKR: int = 5
    LATTICETOL: float = 1E-7
    SIMUUPDATE: bool = False
    SIESTATIMEOUT: int = 30
