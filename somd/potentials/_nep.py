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

from somd import core as _mdcore
from somd.constants import CONSTANTS as _c
from ._nepwrapper import NEPWRAPPER as _NEPWRAPPER

__all__ = ['NEP']


class NEP(_mdcore.potential_base.POTENTIAL):
    """
    The neuroevolution potential, version 3 [1,2].

    Parameters
    ----------
    atom_list : List(int)
        Indices of atoms included by this potential.
    file_name : str
        Name of the potential data file.
    atomic_symbols : List(str)
        Symbols of each element in the simulated system.

    References
    ----------
    .. [1] Fan, Zheyong, et al. "Neuroevolution machine learning potentials:
           Combining high accuracy and low cost in atomistic simulations and
           application to heat transport." Physical Review B 104.10 (2021):
           104309.
    .. [2] Fan, Zheyong. "Improving the accuracy of the neuroevolution machine
           learning potential for multi-component systems." Journal of Physics:
           Condensed Matter 34.12 (2022): 125902.
    """

    def __init__(self,
                 atom_list: list,
                 file_name: str,
                 atomic_symbols: str) -> None:
        """
        Create a NEP instance.
        """
        super().__init__(atom_list)
        self.__nep = _NEPWRAPPER(file_name, atomic_symbols)
        self.__conversion = _c.AVOGACONST * _c.ELECTCONST

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        self.virial[:] = 0.0
        self.energy_potential[0] = \
            self.__nep.calculate(system.box * 10,
                                 system.positions[self.atom_list] * 10,
                                 self.forces, self.virial)
        self.virial[:] *= self.__conversion * 0.001
        self.forces[:] *= self.__conversion * 0.01
        self.energy_potential[0] *= self.__conversion * 0.001
