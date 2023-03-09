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
import warnings as _w
from somd.constants import CONSTANTS as _c

__all__ = ['ATOMGROUP', 'ATOMGROUPS']


class ATOMGROUP(object):
    """
    The atomic group.

    Parameters
    ----------
    system : somd.systems.MDSYSTEM
        The system that contains this group.
    atom_list : List(int)
        IDs of the atoms in this group.
    label : str
        A descriptive string of this group. If no value is given, the label
        will be assigned according to the id of the instance.
    """

    def __init__(self,
                 system,
                 atom_list: list,
                 label: str = None) -> None:
        """
        Create an ATOMGROUP instance.
        """
        if (label is None):
            self._label = 'GROUP_' + str(id(self))
        else:
            self._label = str(label)
        self.__atom_list = _np.array(atom_list, _np.int_)
        self.__atom_list.sort()
        self.__has_translations = True
        if (self.n_atoms > system.n_atoms):
            message = 'Can not initialize group \'{}\' that contains ' + \
                      '{:d} atoms from system \'{}\' with {:d} atoms!'
            message = message.format(self._label, self.__atom_list.size,
                                     system._label, system.n_atoms)
            raise IndexError(message)
        elif (self.__atom_list.max() >= system.n_atoms):
            message = 'Can not initialize group \'{}\' that contains ' + \
                      'a maximum atomic ID of {:d} from system \'{}\' ' + \
                      'with {:d} atoms!'
            message = message.format(self._label, self.__atom_list.max(),
                                     system._label, system.n_atoms)
            raise IndexError(message)
        else:
            self.__system = system
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
        if (self.atom_list.tolist() == g.atom_list.tolist()):
            return True
        else:
            return False

    def __contains__(self, g) -> bool:
        """
        Check if a given group is a subgroup of this group.

        Parameters
        ----------
        g : somd.groups.ATOMGROUP
            The group to check.
        """
        for i in g.atom_list:
            if (i not in self.atom_list):
                return False
        return True

    def overlap_with(self, g) -> bool:
        """
        Check if a given group is overlapping with this group.

        Parameters
        ----------
        g : somd.groups.ATOMGROUP
            The group to check.
        """
        for i in g.atom_list:
            if (i in self.atom_list):
                return True
        return False

    def calculate_n_dof(self) -> int:
        """
        Calculate number of the degree of freedoms of this group.

        Notes
        -----
        This is the most expansive method among all atomic group-related
        methods. So only call this method when constraints/groups are
        added/deleted, and use the cached self.__n_dof variable otherwise.
        """
        self.__n_dof = self.n_atoms * 3
        # calculate number of constraints belonging to this group
        atom_list = self.__atom_list.tolist()
        for c in self.__system.constraints:
            flag = 0
            for i in c['indices']:
                flag += atom_list.count(i)
            if (flag == len(c['indices'])):
                # The constraint belongs to this group.
                self.__n_dof -= 1
            elif (flag != 0) and (flag < len(c['indices'])):
                # The constraint is between two groups.
                message = 'The constraint contains atom indices ' + \
                          '{} is between two groups!'.format(c['indices'])
                raise RuntimeError(message)
        # check if there are subgroups (including this group) binding to COM
        # motion removers, whatever this group itself is in the system's atomic
        # group list.
        for tmp in self.__system.groups:
            if (not tmp.has_translations) and (tmp in self):
                self.__n_dof -= 3
        # in case this group is not in the system's atomic group list
        if (not self.has_translations) and (self not in self.__system.groups):
            self.__n_dof -= 3
        return self.__n_dof

    def remove_com_motion(self, scale_after_removal: bool = False) -> None:
        """
        Remove the translational COM motions of this group.

        Parameters
        ----------
        scale_after_removal : bool
            if scale velocities to keep the total kinetic energy unchanged.
        """
        if (scale_after_removal):
            E_k_1 = self.energy_kinetic
        self.com_velocities = 0.0
        if (scale_after_removal):
            factor = _np.sqrt(E_k_1 / self.energy_kinetic)
            self.velocities *= factor

    def add_velocities_from_temperature(self, temperature: float) -> None:
        """
        Add velocities to atoms in the group according to the given
        temperature.

        Parameters
        ----------
        temperature : float
            The initial temperature. In unit of (K).
        """
        if (self.n_dof == 0):
            message = 'Number of DOF of group \'{}\' has not been calculated!'
            message = message.format(self._label)
            raise RuntimeError(message)
        v = _np.random.randn(self.n_atoms, 3) * \
            _np.sqrt(temperature * _c.BOLTZCONST / self.masses)
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
            _w.warn('Can not remove COM rotational motions!')
        # restore kinetic energies
        t = (_np.square(v) * self.masses).sum() / _c.BOLTZCONST / self.n_dof
        v *= _np.sqrt(temperature / t)
        self.velocities += v

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
        Check inter-group overlapping and set has_translations.
        """
        if (not v):
            for g in self.__system.groups:
                if (not g.has_translations) and (self != g) and \
                        (self.overlap_with(g)):
                    message = 'Atom group \'{}\' is overlapping with ' + \
                              'atomic group \'{}\' while both are binding ' + \
                              'with COM motion removers!'
                    message = message.format(self._label, g._label)
                    raise RuntimeError(message)
        self.__has_translations = bool(v)
        self.calculate_n_dof()
        # check if this is a subgroup and update parents' n_dof.
        if (self in self.__system.groups):
            for g in self.__system.groups:
                if (self in g):
                    g.calculate_n_dof()

    @property
    def n_atoms(self) -> int:
        """
        Number of atoms in this group.
        """
        return len(self.__atom_list)

    @property
    def _system_label(self) -> str:
        """
        Label of the simulated system bound with this group.
        """
        return self.__system._label

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

    @property
    def velocities(self) -> _np.ndarray:
        """
        Velocities of atoms in this group.
        """
        return self.__system.velocities[self.__atom_list]

    @velocities.setter
    def velocities(self, v: _np.ndarray) -> None:
        """
        Set velocities of atoms in this group.
        """
        self.__system.velocities[self.__atom_list] = v

    @property
    def positions(self) -> _np.ndarray:
        """
        Positions of atoms in this group.
        """
        return self.__system.positions[self.__atom_list]

    @positions.setter
    def positions(self, q: _np.ndarray) -> None:
        """
        Set positions of atoms in this group.
        """
        self.__system.positions[self.__atom_list] = q

    @property
    def masses(self) -> _np.ndarray:
        """
        Masses of atoms in this group.
        """
        return self.__system.masses[self.__atom_list]

    @masses.setter
    def masses(self, m: _np.ndarray) -> None:
        """
        Masses of atoms in this group.
        """
        self.__system.masses[self.__atom_list] = m

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
        return self.energy_kinetic / _c.BOLTZCONST / self.n_dof * 2

    @property
    def com_positions(self) -> _np.ndarray:
        """
        Center of mass positions of this group.
        """
        q = (self.positions * self.masses)
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
        v = (self.velocities * self.masses)
        return v.sum(axis=0) / self.masses.sum()

    @com_velocities.setter
    def com_velocities(self, v: _np.ndarray) -> None:
        """
        Set center of mass translational velocities of this group.
        """
        self.velocities += v - self.com_velocities


class ATOMGROUPS(list):
    """
    The atomic groups. This class provides the automatic updates of the groups
    data.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    """

    def __init__(self, system):
        """
        Create an ATOMGROUPS instance.
        """
        super().__init__([])
        self.__system = system

    def __setitem__(self, index: int, item: ATOMGROUP) -> None:
        """
        Update number of DOFs after modifying the groups.
        """
        super().__setitem__(index, item)
        self.update_n_dof()

    def __check_duplication(self, g: ATOMGROUP) -> None:
        """
        Check if the same group has been added to the groups.

        Parameters
        ----------
        g : somd.core.groups.ATOMGROUP
            The group to check.
        """
        for tmp in self:
            if (g == tmp):
                message = 'Group \'{}\' has already been added to the ' + \
                          'to the atomic groups!'
                raise RuntimeError(message.format(g._label))

    def __dict_to_instance(self, d: dict) -> ATOMGROUP:
        """
        Convert a dict that describing an atomic group to an ATOMGROUP
        instance.

        Parameters
        ----------
        d : dict
            The dictionary to convert.
        """
        try:
            a = d['atom_list']
        except:
            raise KeyError('Key atom_list must appear!')
        try:
            l = d['label']
        except:
            l = None
        try:
            t = d['has_translations']
        except:
            t = True
        g = ATOMGROUP(self.__system, a, l)
        g.has_translations = t
        return g

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

    def count(self, group: ATOMGROUP) -> None:
        """
        Disable the count functionality since groups could not be duplicated.
        """
        raise NotImplementedError('Groups could not be duplicated!')

    def append(self, group: ATOMGROUP) -> None:
        """
        Append a new atomic group to the groups.

        Parameters
        ----------
        group : somd.core.groups.ATOMGROUP
             The group to append.
        """
        self.__check_duplication(group)
        super().append(group)
        self.update_n_dof()

    def insert(self, index: int, group: ATOMGROUP) -> None:
        """
        Insert a new atomic group to the groups.

        Parameters
        ----------
        index : int
            The position to insert.
        group : somd.core.groups.ATOMGROUP
             The group to append.
        """
        self.__check_duplication(group)
        super().insert(index, group)
        self.update_n_dof()

    def create_from_dict(self, group_dict: dict) -> None:
        """
        Create a new atomic group and append it to the groups.

        Parameters
        ----------
        group_dict : dict
            The group to create. This dictionary contains four fields:
            - 'atom_list' : List(int)
                IDs of the atoms in this group.
            - 'label' : str
                A descriptive string of this group. If no value is given, the
                label will be assigned according to the id of the instance.
            - 'has_translations' : bool
                If the three translational COM DOFs of this group exchange
                kinetic energies with the internal DOFs.
        """
        g = self.__dict_to_instance(group_dict)
        self.append(g)

    def update_n_dof(self) -> None:
        """
        Recalculate number of degree of freedoms of each atomic groups.
        """
        for g in self:
            # check availability of the COM motion removers
            g.has_translations = g.has_translations
        for g in self:
            # update DOFs
            g.calculate_n_dof()
