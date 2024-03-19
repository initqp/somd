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
The minimal snapshot of a system.
"""

import numpy as _np
import typing as _tp

__all__ = ['SNAPSHOT']


class SNAPSHOT(object):
    """
    The minimal snapshot of a system.

    Parameters
    ----------
    n_atoms : int
        Number of the atom in this snapshot.
    has_mass : bool
        If save mass data in the snapshot.
    """

    def __init__(self, n_atoms: int, has_mass: bool = False) -> None:
        """
        Create a SNAPSHOT instance.
        """
        self.__n_atoms = int(n_atoms)
        self.__allocate(bool(has_mass))

    def __allocate(self, has_mass: bool) -> None:
        """
        Allocate data arrays.

        Parameters
        ----------
        has_mass : bool
            If save mass data in the snapshot.
        """
        self.__box = _np.zeros((3, 3), _np.double)
        self.__virial = _np.zeros((3, 3), _np.double)
        self.__forces = _np.zeros((self.n_atoms, 3), _np.double)
        self.__positions = _np.zeros((self.n_atoms, 3), _np.double)
        self.__velocities = _np.zeros((self.n_atoms, 3), _np.double)
        if has_mass:
            self.__masses = _np.zeros((self.n_atoms, 1), _np.double)
        else:
            self.__masses = None

    def __copy__(self) -> 'SNAPSHOT':
        """
        Clone this snapshot.
        """
        return self.copy()

    def copy(self) -> 'SNAPSHOT':
        """
        Clone this snapshot.
        """
        if self.masses is None:
            snapshot = SNAPSHOT(self.n_atoms)
        else:
            snapshot = SNAPSHOT(self.n_atoms, True)
            snapshot.masses[:] = self.__masses[:]
        snapshot.box[:] = self.__box[:]
        snapshot.forces[:] = self.__forces[:]
        snapshot.virial[:] = self.__virial[:]
        snapshot.positions[:] = self.__positions[:]
        snapshot.velocities[:] = self.__velocities[:]
        return snapshot

    @property
    def n_atoms(self) -> int:
        """
        Number of atoms in this snapshot.
        """
        return self.__n_atoms

    @property
    def masses(self) -> _np.ndarray:
        """
        Atomic masses. In unit of (g/mol).
        """
        return self.__masses

    @property
    def forces(self) -> _np.ndarray:
        """
        Total forces of each atom. In unit of (kJ/mol/nm).
        """
        return self.__forces

    @property
    def positions(self) -> _np.ndarray:
        """
        Unwrapped positions of each atom. In unit of (nm).
        """
        return self.__positions

    @property
    def positions_wrapped(self) -> _np.ndarray:
        """
        Wrapped positions of each atom. In unit of (nm).
        """
        s = (_np.linalg.pinv(self.box.T).dot(self.positions.T)) % 1
        return (self.box.T.dot(s)).T

    @property
    def velocities(self) -> _np.ndarray:
        """
        Velocities of each atom. In unit of (nm/ps).
        """
        return self.__velocities

    @property
    def box(self) -> _np.ndarray:
        """
        Box vectors of this snapshot. In unit of (nm).
        """
        return self.__box

    @property
    def virial(self) -> _np.ndarray:
        """
        Virial tensor. In unit of (kJ/mol).
        """
        return self.__virial

    @property
    def volume(self) -> _np.float64:
        """
        Volume of this snapshot. In unit of (nm^3).
        """
        return _np.dot(self.box[0], _np.cross(self.box[1], self.box[2]))

    @property
    def lattice(self) -> _np.ndarray:
        """
        Lattice parameters of the cell: [a, b, c, alpha, beta, gamma]. In units
        of nm (a, b, c) and degree (alpha, beta, gamma).
        """
        # fmt: off
        result = _np.zeros(6, dtype=_np.double)
        result[:3] = _np.linalg.norm(self.box, axis=1)
        result[3] = _np.arccos(
            self.box[1, :].dot(self.box[2, :]) / result[1] / result[2]
        ) * 180.0 / _np.pi
        result[4] = _np.arccos(
            self.box[0, :].dot(self.box[2, :]) / result[0] / result[2]
        ) * 180.0 / _np.pi
        result[5] = _np.arccos(
            self.box[1, :].dot(self.box[0, :]) / result[0] / result[1]
        ) * 180.0 / _np.pi
        # fmt: on
        return result

    @lattice.setter
    def lattice(self, l: _tp.List[float]) -> None:
        """
        Set lattice parameters of the cell.

        Parameters
        ----------
        l : List[float]
            The lattice parameters: [a, b, c, alpha, beta, gamma]. In units of
            nm (a, b, c) and degree (alpha, beta, gamma).
        """
        if any(_np.array(l) < 1e-6):
            message = 'Very small lattice parameters: {}!'.format(l)
            raise RuntimeError(message)
        self.box[0, 0] = l[0]
        self.box[0, 1] = 0.0
        self.box[0, 2] = 0.0
        self.box[1, 0] = l[1] * _np.cos(l[5] / 180.0 * _np.pi)
        self.box[1, 1] = l[1] * _np.sin(l[5] / 180.0 * _np.pi)
        self.box[1, 2] = 0.0
        self.box[2, 0] = l[2] * _np.cos(l[4] / 180.0 * _np.pi)
        self.box[2, 1] = (
            l[1] * l[2] * _np.cos(l[3] / 180.0 * _np.pi)
            - self.box[2, 0] * self.box[1, 0]
        ) / self.box[1, 1]
        self.box[2, 2] = _np.sqrt(
            l[2] ** 2 - self.box[2, 0] ** 2 - self.box[2, 1] ** 2
        )
