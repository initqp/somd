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
Classes for setting up atom groups.
"""

import numpy as _np
import typing as _tp
import mdtraj as _md
import inspect as _it
from somd import utils as _mdutils
from .snapshots import SNAPSHOT as _SNAPSHOT
from ._lib import CONSTRAINTS as _CONSTRAINTS

__all__ = ['ATOMGROUP', 'ATOMGROUPS']


class ATOMGROUP(object):
    """
    The atom group.

    Parameters
    ----------
    snapshot : somd.snapshots.SNAPSHOT
        The snapshot of the simulated system.
    atom_list : List[int]
        IDs of the atoms in this group.
    label : str
        A descriptive string of this group. If no value is given, the label
        will be assigned according to the id of the instance.
    n_dof_handler : typing.Callable
        The update function number of DOF.
    """

    def __init__(
        self,
        snapshot: _SNAPSHOT,
        atom_list: _tp.List[int],
        label: str = None,
        n_dof_handler: _tp.Callable = None,
    ) -> None:
        """
        Create an ATOMGROUP instance.
        """
        if label is None:
            self._label = 'GROUP_' + str(id(self))
        else:
            self._label = str(label)
        self.__atom_list = _np.sort(_np.array(atom_list, _np.int_))
        if self.n_atoms > snapshot.n_atoms:
            message = (
                'Can not initialize group "{}" that contains '
                + '{:d} atoms from snapshot with {:d} atoms!'
            )
            message = message.format(
                self._label, self.__atom_list.size, snapshot.n_atoms
            )
            raise IndexError(message)
        elif self.__atom_list.max() >= snapshot.n_atoms:
            message = (
                'Can not initialize group "{}" that contains a '
                + 'maximum atomic ID of {:d} from with {:d} atoms!'
            )
            message = message.format(
                self._label, self.__atom_list.max(), snapshot.n_atoms
            )
            raise IndexError(message)
        else:
            self.__snapshot = snapshot
        self.__has_translations = True
        self.__n_dof_handler = n_dof_handler
        self.__n_constraints = 0
        self.__n_dof = 0

    def __eq__(self, g) -> bool:
        """
        Check if two groups are the same according to their atom lists.

        Parameters
        ----------
        g : somd.groups.ATOMGROUP
            The group to compare to.

        Notes
        -----
        This function only compares atom lists of the two groups, but ignores
        other information.
        """
        return self.atom_list.tolist() == g.atom_list.tolist()

    def __contains__(self, g) -> bool:
        """
        Check if a given group is a subgroup of this group.

        Parameters
        ----------
        g : somd.groups.ATOMGROUP
            The group to check.
        """
        atom_set_self = set(self.__atom_list.tolist())
        atom_set_group = set(g.atom_list.tolist())
        return atom_set_group.issubset(atom_set_self)

    def __and__(self, g) -> bool:
        """
        Check if a given group is overlapping with this group.

        Parameters
        ----------
        g : somd.groups.ATOMGROUP
            The group to check.
        """
        atom_set_self = set(self.__atom_list.tolist())
        atom_set_group = set(g.atom_list.tolist())
        return atom_set_group.intersection(atom_set_self)

    def remove_com_motion(self, scale_after_removal: bool = False) -> None:
        """
        Remove the translational COM motions of this group.

        Parameters
        ----------
        scale_after_removal : bool
            if scale velocities to keep the total kinetic energy unchanged.
        """
        if scale_after_removal:
            energy_old = self.energy_kinetic
        self.com_velocities = 0.0
        if scale_after_removal:
            energy_new = self.energy_kinetic
            if energy_new != 0:
                self.velocities *= _np.sqrt(energy_old / energy_new)

    def add_velocities_from_temperature(self, temperature: float) -> None:
        """
        Add velocities to atoms in the group according to the given
        temperature.

        Parameters
        ----------
        temperature : float
            The initial temperature. In unit of (K).
        """
        if self.n_dof == 0:
            message = 'Number of DOF of group "{}" has not been calculated!'
            message = message.format(self._label)
            raise RuntimeError(message)
        if temperature == 0:
            return
        factors = _np.sqrt(
            temperature * _mdutils.constants.BOLTZCONST / self.masses
        )
        v = _mdutils.rng.standard_normal((self.n_atoms, 3)) * factors
        # remove COM translational motions
        v -= (v * self.masses).sum(axis=0) / self.masses.sum()
        # remove COM rotational motions
        d = self.positions - self.com_positions
        L = _np.cross(d, v * self.masses).sum(axis=0)
        I = _np.zeros((3, 3), dtype=_np.double)
        m = self.masses.reshape(self.n_atoms)
        I[0, 0] = ((d[:, 1] * d[:, 1] + d[:, 2] * d[:, 2]) * m).sum()
        I[1, 1] = ((d[:, 0] * d[:, 0] + d[:, 2] * d[:, 2]) * m).sum()
        I[2, 2] = ((d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1]) * m).sum()
        I[0, 1] = (d[:, 0] * d[:, 1] * m).sum() * -1.0
        I[0, 2] = (d[:, 0] * d[:, 2] * m).sum() * -1.0
        I[1, 2] = (d[:, 1] * d[:, 2] * m).sum() * -1.0
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        try:
            w = _np.linalg.pinv(I).dot(L)
            v -= _np.cross(w, d)
        except:
            _mdutils.warning.warn('Can not remove COM rotational motions!')
        # restore kinetic energies
        t = (
            (_np.square(v) * self.masses).sum()
            / _mdutils.constants.BOLTZCONST
            / self.n_dof
        )
        v *= _np.sqrt(temperature / t)
        self.velocities += v

    def to_dict(self) -> _tp.Dict[str, _tp.Any]:
        """
        Output information about this group to a dictionary.
        """
        result = {}
        result['label'] = self._label
        result['atom_list'] = self.atom_list
        result['has_translations'] = self.has_translations
        return result

    @property
    def has_translations(self) -> bool:
        """
        If the three translational COM DOFs of this group exchange kinetic
        energies with the internal DOFs.

        Notes
        -----
            If this option is False, then a COM motion remover should be bound
        to this group, and the number of DOF of this group and the groups that
        full contain it will decrease by three. As a result, when we are
        running simulations under the NVE ensemble, we scale the velocities
        after the COM motion removal to get a constant energy; and just discard
        the COM translational kinetic energies when we are running simulations
        under the NVT ensembles, no matter what thermostat (local/global) we
        are invoking.
            Note that groups binding to COM motion removers must not contain
        overlapping atoms. Such mechanism ensures the total number of DOF of
        the simulated system is correct.
        """
        return self.__has_translations

    @has_translations.setter
    def has_translations(self, v: bool) -> None:
        """
        Set if the three translational COM DOFs of this group exchange kinetic
        energies with the internal DOFs.
        """
        self.__has_translations = bool(v)
        if self.__n_dof_handler is not None:
            self.__n_dof_handler()
        else:
            message = (
                'Group "{:s}" is not bound with any system! '
                + 'Setting the translational flag will be useless!'
            )
            _mdutils.warning.warn(message.format(self._label))

    @property
    def n_atoms(self) -> int:
        """
        Number of atoms in this group.
        """
        return len(self.__atom_list)

    @property
    def atom_list(self) -> str:
        """
        Sorted atom list of this group.
        """
        return self.__atom_list.copy()

    @property
    def n_dof(self) -> int:
        """
        Number of the degree of freedom of this group.
        """
        return self.__n_dof

    @n_dof.setter
    def n_dof(self, v: int) -> None:
        """
        Set the degree of freedom of this group.
        """
        caller = _it.currentframe().f_back.f_code.co_name
        if caller == '__update_n_dof':
            self.__n_dof = v
        else:
            raise RuntimeError('This attribute is handled interally!')

    @property
    def n_constraints(self) -> int:
        """
        Number of constraints that belongs to this group.
        """
        return self.__n_constraints

    @n_constraints.setter
    def n_constraints(self, v: int) -> int:
        """
        Set number of constraints that belongs to this group.
        """
        caller = _it.currentframe().f_back.f_code.co_name
        if caller == '__update_n_constrains':
            self.__n_constraints = v
        else:
            raise RuntimeError('This attribute is handled interally!')

    @property
    def velocities(self) -> _np.ndarray:
        """
        Velocities of atoms in this group.
        """
        return self.__snapshot.velocities[self.__atom_list]

    @velocities.setter
    def velocities(self, v: _np.ndarray) -> None:
        """
        Set velocities of atoms in this group.
        """
        self.__snapshot.velocities[self.__atom_list] = v

    @property
    def positions(self) -> _np.ndarray:
        """
        Positions of atoms in this group.
        """
        return self.__snapshot.positions[self.__atom_list]

    @positions.setter
    def positions(self, q: _np.ndarray) -> None:
        """
        Set positions of atoms in this group.
        """
        self.__snapshot.positions[self.__atom_list] = q

    @property
    def masses(self) -> _np.ndarray:
        """
        Masses of atoms in this group.
        """
        return self.__snapshot.masses[self.__atom_list]

    @masses.setter
    def masses(self, m: _np.ndarray) -> None:
        """
        Masses of atoms in this group.
        """
        self.__snapshot.masses[self.__atom_list] = m

    @property
    def energy_kinetic(self) -> _np.float64:
        """
        Kinetic energy of this group. In unit of (kJ/mol).
        """
        return (_np.square(self.velocities) * self.masses).sum() * 0.5

    @property
    def temperature(self) -> _np.float64:
        """
        Temperature of this group. In unit of (K).
        """
        # fmt: off
        return (
            self.energy_kinetic / _mdutils.constants.BOLTZCONST
            / self.n_dof * 2
        )
        # fmt: on

    @property
    def com_positions(self) -> _np.ndarray:
        """
        Center of mass positions of this group.
        """
        q = self.positions * self.masses
        return q.sum(axis=0) / self.masses.sum()

    @com_positions.setter
    def com_positions(self, q: _np.ndarray) -> None:
        """
        Set center of mass positions of this group by translating the group
        as a rigid body.
        """
        self.positions += q - self.com_positions

    @property
    def com_velocities(self) -> _np.ndarray:
        """
        Center of mass translational velocities of this group.
        """
        v = self.velocities * self.masses
        return v.sum(axis=0) / self.masses.sum()

    @com_velocities.setter
    def com_velocities(self, v: _np.ndarray) -> None:
        """
        Set center of mass translational velocities of this group.
        """
        self.velocities += v - self.com_velocities


class ATOMGROUPS(list):
    """
    The atom groups. This class provides the automatic updates of the groups
    data.

    Parameters
    ----------
    snapshot : somd.snapshots.SNAPSHOT
        The snapshot of the simulated system.
    """

    def __init__(self, snapshot: _SNAPSHOT) -> None:
        """
        Create an ATOMGROUPS instance.
        """
        super().__init__([])
        self.__segments = []
        self.__snapshot = snapshot
        self.__constraints = _CONSTRAINTS(self.__snapshot, self.update_states)

    def __setitem__(self, index: int, item: ATOMGROUP) -> None:
        """
        Update number of DOFs after modifying the groups.
        """
        super().__setitem__(index, item)
        self.update_n_dof()

    def __getitem__(self, index: int) -> ATOMGROUP:
        """
        Get one atom group.
        """
        return super().__getitem__(index)

    def __check_translations(self) -> None:
        """
        Check the availability of the COM motion removers.
        """
        for i, group in enumerate(self):
            flags = [
                ((len(group & g) != 0) and (group != g))
                for g in self
                if (not g.has_translations)
            ]
            if True in flags and (not group.has_translations):
                message = (
                    'Atom group "{}" is overlapping with other atom '
                    + 'group(s) while all of them are binding with '
                    + 'COM motion removers!'
                )
                raise RuntimeError(message.format(group._label))

    def __update_n_constrains(self) -> None:
        """
        Calculate number of constraints that belongs to each group.
        """
        for i, group in enumerate(self):
            group.n_constraints = 0
            atom_set_group = set(group.atom_list.tolist())
            for constraint in self.__constraints:
                atom_set_constraint = set(constraint['indices'])
                if atom_set_constraint.issubset(atom_set_group):
                    group.n_constraints += 1
                elif atom_set_constraint.intersection(atom_set_group):
                    message = (
                        'Can not compute number of constraints that belongs '
                        + 'to group "{:s}"! Because atoms in constraint '
                        + '"{}" are only included partly by this group!'
                    )
                    message = message.format(group._label, constraint)
                    raise RuntimeError(message)

    def __update_n_dof(self) -> _tp.List[int]:
        """
        Calculate number of degree of freedoms of each atom groups.
        """
        self.__check_translations()
        self.__update_n_constrains()
        for i, group in enumerate(self):
            n_dof = group.n_atoms * 3 - group.n_constraints
            # Check if there are subgroups (including this group) binding to
            # COM motion removers.
            flags = [g.has_translations for g in self if g in group]
            n_dof -= 3 * flags.count(False)
            # Set the value.
            group.n_dof = n_dof
        return [g.n_dof for g in self]

    def sort(self) -> None:
        """
        Disable the sort functionality.
        """
        raise NotImplementedError('Groups could not be sorted!')

    def extend(self, index: int) -> None:
        """
        Disable the extend functionality.
        """
        raise NotImplementedError('Groups can not be extended!')

    def insert(self, index: int, group: ATOMGROUP) -> None:
        """
        Disable the insert functionality.
        """
        raise NotImplementedError('Groups can not be inserted!')

    def count(self, group: ATOMGROUP) -> None:
        """
        Disable the count functionality since groups could not be duplicated.
        """
        raise NotImplementedError('Groups could not be duplicated!')

    def append(self, group: ATOMGROUP) -> None:
        """
        Append a new atom group to the groups.

        Parameters
        ----------
        group : somd.core.groups.ATOMGROUP
             The group to append.
        """
        if group.__class__.__name__ != 'ATOMGROUP':
            message = 'Expect an object in type of `ATOMGROUP`!'
            raise RuntimeError(message)
        if group in self:
            message = (
                'Group "{}" has already been added to the to the atom groups!'
            )
            _mdutils.warning.warn(message.format(group._label))
        else:
            super().append(group)
            self.update_n_dof()

    def pop(self, index: int) -> ATOMGROUP:
        """
        Pop an atom group from the groups.

        Parameters
        ----------
        index : int
             Index of the group to pop.
        """
        g = super().pop(index)
        self.update_n_dof()
        return g

    def create_from_dict(self, group_dict: _tp.Dict[str, _tp.Any]) -> None:
        """
        Create a new atom group and append it to the groups.

        Parameters
        ----------
        group_dict : Dict[str, Any]
            The group to create. This dictionary contains four fields:
            - 'atom_list' : List[int]
                IDs of the atoms in this group.
            - 'label' : str
                A descriptive string of this group. If no value is given, the
                label will be assigned according to the id of the instance.
            - 'has_translations' : bool
                If the three translational COM DOFs of this group exchange
                kinetic energies with the internal DOFs.
        """
        if 'atom_list' not in group_dict.keys():
            raise KeyError('The "atom_list" is required to define a group !')
        g = ATOMGROUP(
            self.__snapshot,
            group_dict['atom_list'],
            group_dict.get('label', 'GROUP_' + str(len(self))),
            self.update_n_dof,
        )
        g.has_translations = group_dict.get('has_translations', True)
        self.append(g)

    def find_segments(self) -> _tp.List[ATOMGROUP]:
        """
        Find atom segments (atoms connected by constraints) in this system.
        Each segment is representated by an atomic group, which will not be
        bound to the system. This method should be called after all the
        constraints have been added to the simulated system.
        """
        if len(self.constraints) == 0:
            return
        self.__segments.clear()
        top = _md.Topology()
        unk = top.add_chain()
        residue = top.add_residue('UNK', unk)
        for i in range(0, self.__snapshot.n_atoms):
            e = _md.element.Element.getByAtomicNumber(1)
            top.add_atom(e.symbol + str(i), e, residue)
        atoms = list(top.atoms)
        for c in self.constraints:
            for i in c['indices'][1:len(c['indices'])]:  # fmt: skip
                top.add_bond(atoms[c['indices'][0]], atoms[i])
        molecules = top.find_molecules()
        for m in molecules:
            if len(m) != 1:
                l = [atom.index for atom in list(m)]
                self.__segments.append(ATOMGROUP(self.__snapshot, l))

        return self.__segments

    def update_n_dof(self) -> _tp.List[int]:
        """
        Calculate number of degree of freedoms of each atom groups.
        """
        return self.__update_n_dof()

    def update_states(self) -> None:
        """
        Update internal state. This method should be called after any
        modifications related to groups/constraints.
        """
        self.update_n_dof()
        self.find_segments()

    def summary(self) -> str:
        """
        Show information about the groups.
        """
        result = 'ATOMGROUPS\n'
        result += '┣━ n_groups: {}\n'.format(len(self))
        for i, g in enumerate(self):
            result += '┣━ {} ({}): '.format(g._label, i)
            result += 'n_atoms: {}, '.format(g.n_atoms)
            result += 'n_constraints: {}, '.format(g.n_constraints)
            result += 'translations: {}, '.format(g.has_translations)
            result += 'n_dof: {}'.format(g.n_dof)
            if _mdutils.defaults.VERBOSE:
                result += ', atom_list: {}\n'.format(g.atom_list)
            else:
                result += '\n'
        result += '┗━ END'
        return result

    @property
    def constraints(self) -> _CONSTRAINTS:
        """
        Constraints that belong to this system.
        """
        return self.__constraints

    @property
    def segments(self) -> _tp.List[ATOMGROUP]:
        """
        Atomic segments (atoms connected by constraints) in this system.
        """
        return self.__segments
