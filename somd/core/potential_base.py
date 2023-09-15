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
Base class of all potentials calculators.
"""

import abc as _ab
import numpy as _np
import atexit as _ae
from .systems import MDSYSTEM as _MDSYSTEM

__all__ = ['POTENTIAL']


class POTENTIAL(_ab.ABC):
    """
    Base class of all potentials.

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    """

    def __init__(self, atom_list: list):
        """
        Create a POTENTIAL instance.
        """
        self.__atom_list = _np.array(atom_list, _np.int_)
        self.__virial = _np.zeros((3, 3), _np.double)
        self.__forces = _np.zeros((self.n_atoms, 3), _np.double)
        self.__energy_potential = _np.zeros((1), _np.double)
        _ae.register(self.finalize)

    @_ab.abstractmethod
    def update(self, system: _MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        raise NotImplementedError()

    @classmethod
    def generator(cls, *args, **kwargs) -> callable:
        """
        Return a generator of this potential.
        """
        return lambda x=args, y=kwargs: cls(*x, **y)

    def reset(self) -> None:
        """
        Reset the potential.
        """
        pass

    def finalize(self) -> None:
        """
        Clean up.
        """
        # If this function has been called once, do not call it again at exit.
        _ae.unregister(self.finalize)

    @property
    def n_atoms(self) -> int:
        """
        Number of atoms in this group.
        """
        return len(self.__atom_list)

    @property
    def forces(self) -> _np.ndarray:
        """
        The atomic forces caused by this potential. In unit of (kJ/mol/nm).
        """
        return self.__forces

    @property
    def energy_potential(self) -> _np.ndarray:
        """
        The potential energy. In unit of (kJ/mol).
        """
        return self.__energy_potential

    @property
    def virial(self) -> _np.ndarray:
        """
        The virial tensor caused by this potential. In unit of (kJ/mol).
        """
        return self.__virial

    @property
    def atom_list(self) -> _np.ndarray:
        """
        Indices of atoms included by this potential.
        """
        return self.__atom_list
