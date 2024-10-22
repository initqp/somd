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
The simulated system.
"""

import numpy as _np
import typing as _tp
import mdtraj as _md
from threading import Thread as _Thread
from somd import utils as _mdutils
from .groups import ATOMGROUP as _ATOMGROUP
from .groups import ATOMGROUPS as _ATOMGROUPS
from .snapshots import SNAPSHOT as _SNAPSHOT
from ._lib import CONSTRAINTS as _CONSTRAINTS

__all__ = ['MDSYSTEM', 'create_system_from_pdb', 'create_system_from_poscar']


class MDSYSTEM(object):
    """
    The simulated system.

    Parameters
    ----------
    n_atoms : int
        Number of the atom in the simulated system.
    label : str
        A descriptive string of this system. If no value is given, the label
        will be assigned according to the id of the instance.
    """

    def __init__(self, n_atoms: int, label: str = None) -> None:
        """
        Create a MDSYSTEM instance.
        """
        if label is None:
            self._label = 'SYSTEM_' + str(id(self))
        else:
            self._label = str(label)
        self.__allocate(n_atoms)

    def __allocate(self, n_atoms: int) -> None:
        """
        Allocate data arrays.

        Parameters
        ----------
        n_atoms : int
            Number of the atom in the simulated system.
        """
        self.__snapshot = _SNAPSHOT(n_atoms, True)
        self.__types = _np.zeros((n_atoms, 1), _np.int_)
        self.__energy_potential = _np.zeros((1), _np.double)
        self.__groups = _ATOMGROUPS(self.snapshot)
        self.__potentials = []
        self.__atomic_symbols = []

    def __copy__(self) -> 'MDSYSTEM':
        """
        Clone the system without potential calculators.
        """
        return self.copy()

    def copy(self) -> 'MDSYSTEM':
        """
        Clone the system without potential calculators.
        """
        system = MDSYSTEM(self.n_atoms)
        system._label = self._label
        system.snapshot = self.snapshot
        system.atomic_types[:] = self.atomic_types[:]
        system.atomic_symbols[:] = self.atomic_symbols[:]
        for group in self.groups:
            system.groups.create_from_dict(group.to_dict())
        for constraint in self.groups.constraints:
            system.groups.constraints.append(constraint)
        return system

    def update_potentials(
        self, indices: _tp.List[int] = None, perform_calculations: bool = True
    ) -> None:
        """
        Invoke the force calculators.

        Parameters
        ----------
        indices : List[int]
            Indices of the force calculators to invoke. The default behavior is
            calling all the calculators.
        perform_calculations : bool
            If actually evolve potentials. If false, the potential data will be
            simply copied from the calculators.
        """
        self.virial[:] = 0.0
        self.forces[:] = 0.0
        self.__energy_potential[0] = 0.0
        if indices is None:
            l = range(0, len(self.potentials))
        else:
            l = indices
        if perform_calculations:
            if _mdutils.defaults.SIMUUPDATE:
                threads = []
                for i in l:
                    t = _Thread(
                        target=self.__potentials[i].update, args=(self,)
                    )
                    threads.append(t)
                    t.start()
                for t in threads:
                    t.join()
            else:
                for i in l:
                    self.__potentials[i].update(self)
        for i in l:
            self.virial[:] += self.__potentials[i].virial
            f = self.__potentials[i].forces
            self.forces[self.__potentials[i].atom_list] += f
            e = self.__potentials[i].energy_potential[0]
            self.__energy_potential[0] += e

    def summary(self) -> str:
        """
        Show information about the system.
        """
        box = self.box.squeeze().tolist()
        segments = [s.atom_list.tolist() for s in self.segments]
        type_list = self.atomic_types.squeeze().tolist()
        mass_map = {
            t: self.masses[type_list.index(t)][0] for t in set(type_list)
        }
        summary_g = self.groups.summary().replace('\n', '\n┃  ').strip()
        summary_c = self.constraints.summary().replace('\n', '\n┃  ').strip()
        summary_p = 'POTENTIALS\n'
        for p in self.potentials:
            summary_p += (
                '┃  ┣━ ' + p.summary().replace('\n', '\n┃  ┃  ').strip() + '\n'
            )
        summary_p += '┃  ┗━ END'

        result = 'MDSYSTEM\n'
        result += '┣━ n_atoms: {}\n'.format(self.n_atoms)
        result += '┣━ n_atomic_types: {}\n'.format(len(set(type_list)))
        if _mdutils.defaults.VERBOSE:
            result += '┣━ atomic_types: {}\n'.format(type_list)
        result += '┣━ atomic_masses: {}\n'.format(mass_map)
        result += '┣━ n_potentials: {}\n'.format(len(self.potentials))
        result += '┣━ n_atom_groups: {}\n'.format(len(self.groups))
        result += '┣━ n_constraints: {}\n'.format(len(self.constraints))
        result += '┣━ n_segments: {}\n'.format(len(self.segments))
        if _mdutils.defaults.VERBOSE:
            result += '┣━ segments: {}\n'.format(segments)
            result += '┣━ initial box vector: {}\n'.format(box)
        result += '┣━ ' + summary_g + '\n'
        if len(self.constraints) > 0:
            result += '┣━ ' + summary_c + '\n'
        result += '┣━ ' + summary_p + '\n'
        result += '┗━ END'

        return result

    @property
    def snapshot(self) -> _SNAPSHOT:
        """
        The current state of this system.
        """
        return self.__snapshot

    @snapshot.setter
    def snapshot(self, f: _SNAPSHOT) -> None:
        """
        Copy the current state of this system from another snapshot.
        """
        if self.n_atoms != f.n_atoms:
            message = (
                'Mismatch of number of atoms between the target '
                + 'snapshot and the system!'
            )
            raise RuntimeError(message)
        self.box[:] = f.box[:]
        self.forces[:] = f.forces[:]
        self.virial[:] = f.virial[:]
        self.positions[:] = f.positions[:]
        self.velocities[:] = f.velocities[:]
        if f.masses is not None:
            self.masses[:] = f.masses[:]

    @property
    def n_atoms(self) -> int:
        """
        Number of atoms in the simulated system.
        """
        return self.__snapshot.n_atoms

    @property
    def forces(self) -> _np.ndarray:
        """
        Total forces of each atom. In unit of (kJ/mol/nm).
        """
        return self.__snapshot.forces

    @property
    def masses(self) -> _np.ndarray:
        """
        Atomic masses. In unit of (g/mol).
        """
        return self.__snapshot.masses

    @property
    def positions(self) -> _np.ndarray:
        """
        Unwrapped positions of each atom. In unit of (nm).
        """
        return self.__snapshot.positions

    @property
    def positions_wrapped(self) -> _np.ndarray:
        """
        Wrapped positions of each atom. In unit of (nm).
        """
        return self.__snapshot.positions_wrapped

    @property
    def velocities(self) -> _np.ndarray:
        """
        Velocities of each atom. In unit of (nm/ps).
        """
        return self.__snapshot.velocities

    @property
    def box(self) -> _np.ndarray:
        """
        Box vectors of the simulated system. In unit of (nm).
        """
        return self.__snapshot.box

    @property
    def virial(self) -> _np.ndarray:
        """
        Virial tensor. In unit of (kJ/mol).
        """
        return self.__snapshot.virial

    @property
    def lattice(self) -> _np.ndarray:
        """
        Lattice parameters of the cell: [a, b, c, alpha, beta, gamma]. In units
        of nm (a, b, c) and degree (alpha, beta, gamma).
        """
        return self.__snapshot.lattice

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
        self.__snapshot.lattice = l

    @property
    def atomic_types(self) -> _np.ndarray:
        """
        Atomic types.
        """
        return self.__types

    @property
    def atomic_symbols(self) -> list:
        """
        Atomic symbols.
        """
        return self.__atomic_symbols

    @property
    def pressures(self) -> _np.ndarray:
        """
        Stress tensor. In unit of (kJ/mol/nm^3).
        """
        I = _np.eye(3, dtype=_np.double)
        E_k = (self.velocities * self.velocities * self.masses).sum(axis=0)
        return (I * E_k + self.virial) / self.volume

    @property
    def energy_potential(self) -> float:
        """
        Potential energy of this system. In unit of (kJ/mol).
        """
        return self.__energy_potential[0]

    @property
    def groups(self) -> _ATOMGROUPS:
        """
        Atomic groups that belong to this system.
        """
        return self.__groups

    @property
    def segments(self) -> _tp.List[_ATOMGROUP]:
        """
        Atomic segments (atoms connected by constraints) in this system.
        """
        return self.__groups.segments

    @property
    def constraints(self) -> _CONSTRAINTS:
        """
        Constraints that belong to this system.
        """
        return self.__groups.constraints

    @property
    def volume(self) -> _np.float64:
        """
        Volume of the simulated system. In unit of (nm^3).
        """
        return self.__snapshot.volume

    @property
    def potentials(self) -> list:
        """
        The potentials.
        """
        return self.__potentials


def create_system_from_pdb(file_name: str) -> MDSYSTEM:
    """
    Create a MDSYSTEM from a PDB file.

    Parameters
    ----------
    file_name : str
        Name of the PDB file.
    """
    pdb = _md.load_pdb(file_name)
    if pdb.n_atoms == 0:
        raise RuntimeError('PDB file "{}" has no atom!'.format(file_name))
    else:
        s = MDSYSTEM(pdb.n_atoms)
    if pdb.unitcell_vectors.any() is None:
        raise RuntimeError('PDB file "{}" has no cell data!'.format(file_name))
    else:
        try:
            s.box[:, :] = pdb.unitcell_vectors
        except:
            message = 'Can not read unit cell data from file "{:s}"!'
            _mdutils.warning.warn(message.format(file_name))
    for i in range(0, pdb.n_atoms):
        s.atomic_types[i] = pdb.top.atom(i).element.number
        s.atomic_symbols.append(pdb.top.atom(i).element.symbol)
        s.masses[i] = pdb.top.atom(i).element.mass
    s.positions[:, :] = pdb.xyz[0]
    return s


def create_system_from_poscar(file_name: str) -> MDSYSTEM:
    """
    Create a MDSYSTEM instance from a POSCAR file [1].

    Parameters
    ----------
    file_name : str
        Name of the POSCAR file.

    References
    ----------
    .. [1] https://www.vasp.at/wiki/index.php/POSCAR
    """
    fp = open(file_name, 'r')
    # read the header
    fp.readline().strip()
    scale_factor = float(fp.readline().strip())
    box = _np.zeros((3, 3), dtype=_np.double)
    for i in range(0, 3):
        box[i, :] = _np.array(fp.readline().split(), dtype=_np.double)
    elements = fp.readline().strip().split()
    n_atoms = _np.array(fp.readline().strip().split(), dtype=_np.int_)
    # build the system.
    count = 0
    s = MDSYSTEM(n_atoms.sum())
    for i in range(0, len(elements)):
        e = _md.element.get_by_symbol(elements[i])
        for j in range(0, n_atoms[i]):
            s.atomic_types[count] = e.number
            s.atomic_symbols.append(elements[i])
            s.masses[count] = e.mass
            count += 1
    s.box[:] = box * scale_factor * 0.1
    for i in range(0, 3):
        if s.lattice[i] < _mdutils.defaults.LATTICETOL:
            raise RuntimeError('Very small cell length in dim {:d}'.format(i))
    # read the positions
    position_type = fp.readline().strip().lower()
    for i in range(0, s.n_atoms):
        s.positions[i, :] = (
            _np.array(fp.readline().strip().split(), dtype=_np.double) * 0.1
        )
    if position_type == 'direct':
        s.positions[:] = (s.box.T.dot(s.positions)).T
    fp.close()
    return s
