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
from somd.utils import constants as _c

__all__ = ['TBLITE']


class TBLITE(_mdcore.potential_base.POTENTIAL):
    """
    The tight binding potential calculated with the TBLite code [1].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    atomic_types : List[int]
        Types of atoms included by this potential.
    method : str
        Name of the Hamiltonian. Valid names are:
        - 'GFN2-xTB' [2]
        - 'GFN1-xTB' [3]
    total_charges : int
        Total charges of the atoms included by this potential.
    total_spins : int
        Total spins (N_alpha - N_beta) of the atoms included by this potential.
    pbc : bool
        If enable PBC.

    References
    ----------
    .. [1] https://github.com/tblite/tblite
    .. [2] Bannwarth, Christoph, Sebastian Ehlert, and Stefan Grimme.
           "GFN2-xTB-An accurate and broadly parametrized self-consistent
           tight-binding quantum chemical method with multipole electrostatics
           and density-dependent dispersion contributions." Journal of chemical
           theory and computation 15.3 (2019): 1652-1671.
    .. [3] Grimme, Stefan, Christoph Bannwarth, and Philip Shushkov. "A robust
           and accurate tight-binding quantum chemical method for structures,
           vibrational frequencies, and noncovalent interactions of large
           molecular systems parametrized for all spd-block elements (Z=1-86)."
           Journal of chemical theory and computation 13.5 (2017): 1989-2009.
    """

    def __init__(
        self,
        atom_list: _tp.List[int],
        atomic_types: _tp.List[int],
        method: str = 'GFN1-xTB',
        total_charges: int = 0,
        total_spins: int = 0,
        pbc: int = True
    ) -> None:
        """
        Creat a TBLITE instance.
        """
        super().__init__(atom_list)
        # treat TBLite as a local dependency
        try:
            from tblite import library
            from tblite import _libtblite
            from tblite.interface import Calculator

            # the default logger is way too noisy ...
            @_libtblite.ffi.def_extern()
            def logger_callback(message, nchar, data):
                pass

            library.logger_callback = logger_callback
        except:
            raise ImportError(
                'You need to have the tblite-python package '
                + 'installed to use the tight binding potential!'
            )
        self.__conversion = _c.HARTREE / _c.BOHRRADIUS * -1.0
        if pbc:
            pbc = _np.ones(3, dtype=_np.int_)
        else:
            pbc = _np.zeros(3, dtype=_np.int_)
        lattice = _np.ones((3, 3), dtype=_np.double)
        atomic_types = _np.array(atomic_types, dtype=_np.int_)
        positions = _np.ones(
            (atomic_types.shape[0], 3), dtype=_np.double
        ) * _np.arange(0, atomic_types.shape[0]).reshape(
            atomic_types.shape[0], 1
        )
        self.__calculator = Calculator(
            method,
            atomic_types,
            positions,
            total_charges,
            total_spins,
            lattice,
            pbc
        )

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        self.__calculator.update(
            system.positions[self.atom_list] / _c.BOHRRADIUS,
            system.box / _c.BOHRRADIUS,
        )
        result = self.__calculator.singlepoint()
        self.energy_potential[0] = result.get('energy') * _c.HARTREE
        self.forces[:] = result.get('gradient') * self.__conversion
        self.virial[:] = result.get('virial') * _c.HARTREE * -1.0
