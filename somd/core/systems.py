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
import mdtraj as _md
import warnings as _w
from threading import Thread as _Thread
from .groups import ATOMGROUP as _ATOMGROUP
from .groups import ATOMGROUPS as _ATOMGROUPS
from ._lib import CONSTRAINTS as _CONSTRAINTS
from somd.constants import SOMDDEFAULTS as _d

__all__ = ['SNAPSHOT',
           'MDSYSTEM',
           'create_system_from_pdb',
           'create_system_from_poscar']


class SNAPSHOT(object):
    """
    The minimal snapshot of a system.

    Parameters
    ----------
    n_atoms : int
        Number of the atom in this snapshot.
    """

    def __init__(self, n_atoms: int) -> None:
        """
        Create a SNAPSHOT instance.
        """
        self.__n_atoms = int(n_atoms)
        self.__allocate()

    def __allocate(self) -> None:
        """
        Allocate data arrays.
        """
        self.__box = _np.zeros((3, 3), _np.double)
        self.__forces = _np.zeros((self.n_atoms, 3), _np.double)
        self.__positions = _np.zeros((self.n_atoms, 3), _np.double)
        self.__velocities = _np.zeros((self.n_atoms, 3), _np.double)

    def copy(self) -> 'SNAPSHOT':
        """
        Clone this snapshot.
        """
        snapshot = SNAPSHOT(self.n_atoms)
        snapshot.box[:] = self.__box[:]
        snapshot.forces[:] = self.__forces[:]
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
        result = _np.zeros(6, dtype=_np.double)
        for i in range(0, 3):
            result[i] = _np.linalg.norm(self.box[i, :])
        result[3] = _np.arccos(self.box[1, :].dot(self.box[2, :]) /
                               result[1] / result[2]) * 180.0 / _np.pi
        result[4] = _np.arccos(self.box[0, :].dot(self.box[2, :]) /
                               result[0] / result[2]) * 180.0 / _np.pi
        result[5] = _np.arccos(self.box[1, :].dot(self.box[0, :]) /
                               result[0] / result[1]) * 180.0 / _np.pi
        return result

    @lattice.setter
    def lattice(self, l: list) -> None:
        """
        Set lattice parameters of the cell.

        Parameters
        ----------
        l : List(float)
            The lattice parameters: [a, b, c, alpha, beta, gamma]. In units of
            nm (a, b, c) and degree (alpha, beta, gamma).
        """
        parameter_names = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        for i in range(0, 6):
            if (l[i] < _d.LATTICETOL):
                message = 'Very small lattice parameter {}: {:d}'
                raise RuntimeError(message.format(parameter_names[i], l[i]))
        self.box[0, 0] = l[0]
        self.box[0, 1] = 0.0
        self.box[0, 2] = 0.0
        self.box[1, 0] = l[1] * _np.cos(l[5] / 180.0 * _np.pi)
        self.box[1, 1] = l[1] * _np.sin(l[5] / 180.0 * _np.pi)
        self.box[1, 2] = 0.0
        self.box[2, 0] = l[2] * _np.cos(l[4] / 180.0 * _np.pi)
        self.box[2, 1] = (l[1] * l[2] * _np.cos(l[3] / 180.0 * _np.pi) -
                          self.box[2, 0] * self.box[1, 0]) / self.box[1, 1]
        self.box[2, 2] = _np.sqrt(l[2] ** 2 - self.box[2, 0] ** 2 -
                                  self.box[2, 1] ** 2)


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
        if (label is None):
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
        self.__snapshot = SNAPSHOT(n_atoms)
        self.__types = _np.zeros((n_atoms, 1), _np.int_)
        self.__masses = _np.zeros((self.n_atoms, 1), _np.double)
        self.__virial = _np.zeros((3, 3), _np.double)
        self.__energy_potential = _np.zeros((1), _np.double)
        self.__constraints = _CONSTRAINTS(self)
        self.__groups = _ATOMGROUPS(self)
        self.__segments = []
        self.__potentials = []
        self.__atomic_symbols = []

    def copy(self) -> 'MDSYSTEM':
        """
        Clone the system without potential calculators.
        """
        system = MDSYSTEM(self.n_atoms)
        system._label = self._label
        system.box[:] = self.box[:]
        system.masses[:] = self.masses[:]
        system.forces[:] = self.forces[:]
        system.virial[:] = self.virial[:]
        system.positions[:] = self.positions[:]
        system.velocities[:] = self.velocities[:]
        system.atomic_types[:] = self.atomic_types[:]
        system.atomic_symbols[:] = self.atomic_symbols[:]
        for group in self.groups:
            d = {'atom_list': group.atom_list, 'label': group._label,
                 'has_translations': group.has_translations}
            system.groups.create_from_dict(d)
        for constraint in self.constraints:
            system.constraints.append(constraint)
        system.find_segments()
        return system

    def update_potentials(self, indices: list = None) -> None:
        """
        Invoke the force calculators.

        Parameters
        ----------
        indices : List(int)
            Indices of the force calculators to invoke. The default behavior is
            calling all the calculators.
        """
        self.virial[:] = 0.0
        self.forces[:] = 0.0
        self.__energy_potential[0] = 0.0
        if (indices is None):
            l = range(0, len(self.potentials))
        else:
            l = indices
        if (_d.SIMUUPDATE):
            threads = []
            for i in l:
                t = _Thread(target=self.__potentials[i].update, args=(self,))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()
        else:
            for i in l:
                self.__potentials[i].update(self)
        for i in l:
            self.virial[:] += self.__potentials[i].virial
            self.forces[self.__potentials[i].atom_list] += \
                self.__potentials[i].forces
            self.__energy_potential[0] += \
                self.__potentials[i].energy_potential[0]

    def find_segments(self) -> None:
        """
        Find atom segments (atoms connected by constraints) in this system.
        Each segment is representated by an atomic group, which will not be
        bound to the system. This method should be called after all the
        constraints have been added to the simulated system.
        """
        if (len(self.constraints) == 0):
            return
        self.__segments.clear()
        top = _md.Topology()
        unk = top.add_chain()
        residue = top.add_residue('UNK', unk)
        for i in range(0, len(self.atomic_types)):
            e = _md.element.Element.getByAtomicNumber(self.atomic_types[i, 0])
            top.add_atom(e.symbol + str(i), e, residue)
        atoms = list(top.atoms)
        for c in self.constraints:
            for i in c['indices'][1:len(c['indices'])]:
                top.add_bond(atoms[c['indices'][0]], atoms[i])
        molecules = top.find_molecules()
        for m in molecules:
            if (len(m) != 1):
                l = [atom.index for atom in list(m)]
                self.__segments.append(_ATOMGROUP(self, l))

    @property
    def snapshot(self) -> SNAPSHOT:
        """
        The current state of this system.
        """
        return self.__snapshot

    @snapshot.setter
    def snapshot(self, f: SNAPSHOT) -> None:
        """
        Copy the current state of this system from another snapshot.
        """
        if (self.n_atoms != f.n_atoms):
            message = 'Mismatch of number of atoms between the target ' + \
                      'snapshot and the system!'
            raise RuntimeError(message)
        self.box[:] = f.box[:]
        self.forces[:] = f.forces[:]
        self.positions[:] = f.positions[:]
        self.velocities[:] = f.velocities[:]

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
        return self.__virial

    @property
    def lattice(self) -> _np.ndarray:
        """
        Lattice parameters of the cell: [a, b, c, alpha, beta, gamma]. In units
        of nm (a, b, c) and degree (alpha, beta, gamma).
        """
        return self.__snapshot.lattice

    @lattice.setter
    def lattice(self, l: list) -> None:
        """
        Set lattice parameters of the cell.

        Parameters
        ----------
        l : List(float)
            The lattice parameters: [a, b, c, alpha, beta, gamma]. In units of
            nm (a, b, c) and degree (alpha, beta, gamma).
        """
        self.__snapshot.lattice = l

    @property
    def masses(self) -> _np.ndarray:
        """
        Atomic masses. In unit of (g/mol).
        """
        return self.__masses

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
    def segments(self) -> list:
        """
        Atomic segments (atoms connected by constraints) in this system.
        """
        return self.__segments

    @property
    def constraints(self) -> _CONSTRAINTS:
        """
        Constraints that belong to this system.
        """
        return self.__constraints

    @property
    def volume(self) -> _np.float64:
        """
        Volume of the simulated system. In unit of (nm^3).
        """
        return _np.dot(self.box[0], _np.cross(self.box[1], self.box[2]))

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
    if (pdb.n_atoms == 0):
        raise RuntimeError('PDB file {} has no atom!'.format(file_name))
    else:
        s = MDSYSTEM(pdb.n_atoms)
    if (pdb.unitcell_vectors.any() is None):
        raise RuntimeError('PDB file {} has no cell data!'.format(file_name))
    else:
        try:
            s.box[:, :] = pdb.unitcell_vectors
        except:
            message = 'Can not read unit cell data from file ' + file_name
            _w.warn(message)
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
        if (s.lattice[i] < _d.LATTICETOL):
            raise RuntimeError('Very small cell length in dim {:d}'.format(i))
    # read the positions
    position_type = fp.readline().strip().lower()
    for i in range(0, s.n_atoms):
        s.positions[i, :] = 0.1 * \
            _np.array(fp.readline().strip().split(), dtype=_np.double)
    if (position_type == 'direct'):
        s.positions[:] = (s.box.T.dot(s.positions)).T
    fp.close()
    return s
