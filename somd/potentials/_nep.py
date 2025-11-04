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

import os as _os
import typing as _tp
from somd import core as _mdcore
from somd.utils import defaults as _d
from somd.utils import constants as _c

__all__ = ['NEP']


class NEP(_mdcore.potential_base.POTENTIAL):
    """
    The neuroevolution potential, version 3 [1,2].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    file_name : str
        Name of the potential data file.
    atomic_symbols : List[str]
        Symbols of each element in the simulated system.
    use_tabulating : bool
        If invoke the tabulated version of NEP.

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

    def __init__(
        self,
        atom_list: _tp.List[int],
        file_name: str,
        atomic_symbols: _tp.List[str],
        use_tabulating: bool = False,
    ) -> None:
        """
        Create a NEP instance.
        """
        super().__init__(atom_list)
        self.__file_name = file_name
        self.__use_tabulating = use_tabulating
        # only fail at runtime
        try:
            if use_tabulating:
                from ._nepwrapper_t import NEPWRAPPER as _NEPWRAPPER
            else:
                from ._nepwrapper import NEPWRAPPER as _NEPWRAPPER
            self.__nep = _NEPWRAPPER(file_name, atomic_symbols)
        except ImportError:
            self.__nep = None
        self.__conversion = _c.AVOGACONST * _c.ELECTCONST

    def summary(self) -> str:
        """
        Show information about the potential.
        """
        if self.__nep is None:
            raise RuntimeError('NEP not installed!')

        zbl = self.__nep.zbl_info
        ann = self.__nep.ann_info
        paramb = self.__nep.paramb_info
        result = '{}\n'.format(self.__class__.__name__)
        result += '┣━ n_atoms: {}\n'.format(self.n_atoms)
        result += '┣━ file_name: {}\n'.format(self.__file_name)
        result += '┣━ ParaMB:\n'
        for k in paramb.keys():
            result += '┃  ┣━ {}: {}\n'.format(k, paramb[k])
        result += '┃  ┗━ END\n'
        result += '┣━ ZBL:\n'
        for k in zbl.keys():
            result += '┃  ┣━ {}: {}\n'.format(k, zbl[k])
        result += '┃  ┗━ END\n'
        result += '┣━ ANN:\n'
        for k in ann.keys():
            result += '┃  ┣━ {}: {}\n'.format(k, ann[k])
        result += '┃  ┗━ END\n'
        result += '┣━ tabulating: {}\n'.format(self.__use_tabulating)
        if _d.VERBOSE:
            result += '┣━ atom_list: {}\n'.format(self.atom_list)
        result += '┗━ END'
        return result

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        if self.__nep is None:
            raise RuntimeError('NEP not installed!')

        self.virial[:] = 0.0
        self.energy_potential[0] = self.__nep.calculate(
            system.box * 10,
            system.positions[self.atom_list] * 10,
            self.forces,
            self.virial,
        )
        self.virial[:] *= self.__conversion * 0.001
        self.forces[:] *= self.__conversion * 0.01
        self.energy_potential[0] *= self.__conversion * 0.001

    @classmethod
    def generator(cls, *args, **kwargs) -> _tp.Callable:
        """
        Return a generator of this potential.
        """
        if 'file_name' in kwargs.keys():
            kwargs['file_name'] = _os.path.abspath(kwargs['file_name'])
        else:
            args = list(args)
            args[1] = _os.path.abspath(args[1])
        return lambda x=tuple(args), y=kwargs: cls(*x, **y)

    def finalize(self) -> None:
        """
        Clean up.
        """
        super().finalize()

        if self.__nep is not None:
            self.__nep.dealloc()
