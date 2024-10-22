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

import numpy as _np
import typing as _tp
from somd import core as _mdcore
from somd.utils import defaults as _d
from somd.utils import constants as _c

__all__ = ['DFTD3']


class DFTD3(_mdcore.potential_base.POTENTIAL):
    """
    The charge dependent London-dispersion correction potential [1].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    atomic_types : List[int]
        Types of atoms included by this potential.
    method : str
        Name of the claculation method.
    damping : str
        Type of the damping parameters. Valid options are:
        - 'ZeroDamping'
        - 'RationalDamping'
        - 'ModifiedZeroDamping'
        - 'ModifiedRationalDamping'
        - 'OptimizedPowerDamping'
    atm : bool
        If use the three-body correction.

    References
    ----------
    .. [1] Grimme, Stefan, et al. "A consistent and accurate ab initio
           parametrization of density functional dispersion correction (DFT-D)
           for the 94 elements H-Pu." The Journal of chemical physics 132.15
           (2010): 154104.
       [2] Grimme, Stefan, Stephan Ehrlich, and Lars Goerigk. "Effect of the
           damping function in dispersion corrected density functional theory."
           Journal of computational chemistry 32.7 (2011): 1456-1465.
       [3] Smith, Daniel GA, et al. "Revised damping parameters for the D3
           dispersion correction to density functional theory." The journal of
           physical chemistry letters 7.12 (2016): 2197-2203.
    """

    def __init__(
        self,
        atom_list: _tp.List[int],
        atomic_types: _tp.List[int],
        method: str,
        damping: str = 'ZeroDamping',
        atm: bool = False,
    ) -> None:
        """
        Create a DFTD3 instance.
        """
        super().__init__(atom_list)
        # treat DFTD3 as a local dependency
        try:
            from dftd3.interface import DispersionModel
            from dftd3.interface import ZeroDampingParam
            from dftd3.interface import RationalDampingParam
            from dftd3.interface import ModifiedZeroDampingParam
            from dftd3.interface import OptimizedPowerDampingParam
            from dftd3.interface import ModifiedRationalDampingParam
        except:
            raise ImportError(
                'You need to have the dftd3-python package '
                + 'installed to use the DFTD3 potential!'
            )
        self.__atm = atm
        self.__method = method
        self.__conversion = _c.HARTREE / _c.BOHRRADIUS * -1.0
        pbc = _np.ones(3, dtype=_np.int_)
        lattice = _np.zeros((3, 3), dtype=_np.double)
        atomic_types = _np.array(atomic_types, dtype=_np.int_)
        positions = _np.zeros((atomic_types.shape[0], 3), dtype=_np.double)
        self.__model = DispersionModel(atomic_types, positions, lattice, pbc)
        if damping == 'ZeroDamping':
            self.__param = ZeroDampingParam(method=method, atm=atm)
        elif damping == 'RationalDamping':
            self.__param = RationalDampingParam(method=method, atm=atm)
        elif damping == 'ModifiedZeroDamping':
            self.__param = ModifiedZeroDampingParam(method=method, atm=atm)
        elif damping == 'OptimizedPowerDamping':
            self.__param = OptimizedPowerDampingParam(method=method, atm=atm)
        elif damping == 'ModifiedRationalDamping':
            self.__param = ModifiedRationalDampingParam(method=method, atm=atm)
        else:
            raise RuntimeError('Unknown damping parameter type: ' + damping)

    def summary(self) -> str:
        """
        Show information about the potential.
        """
        result = '{}\n'.format(self.__class__.__name__)
        result += '┣━ n_atoms: {}\n'.format(self.n_atoms)
        result += '┣━ method: {}\n'.format(self.__method)
        result += '┣━ parameter: {}\n'.format(self.__param.__class__.__name__)
        result += '┣━ atm: {}\n'.format(self.__atm)
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
        self.__model.update(
            system.positions[self.atom_list] / _c.BOHRRADIUS,
            system.box / _c.BOHRRADIUS,
        )
        result = self.__model.get_dispersion(self.__param, grad=True)
        self.energy_potential[0] = result.get("energy") * _c.HARTREE
        self.forces[:] = result.get("gradient") * self.__conversion
        self.virial[:] = result.get("virial") * _c.HARTREE * -1.0
