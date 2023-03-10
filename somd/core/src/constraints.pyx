# cython: language_level=3, boundscheck=False, wraparound=False
# distutils: language = c++

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

from libcpp.vector cimport vector
from constraints cimport RATTLE

__all__ = ['CONSTRAINTS']

cdef class CONSTRAINTS(object):
    """
    Define and apply constraints using the RATTLE method.

    References
    ----------
    .. [1] Andersen, H. C. (1983). Rattle: A "velocity" version of the shake
        algorithm for molecular dynamics calculations. Journal of computational
        Physics, 52(1), 24-34.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
    """

    # Pointer to the internal C++ object.
    cdef object __system
    cdef RATTLE *__cxx_obj_ptr

    def __cinit__(self, system) -> None:
        """
        Initialize the internal data.
        """
        self.__system = system
        self.__cxx_obj_ptr = new RATTLE()

    def __dealloc__(self) -> None:
        """
        Finalize the internal data.
        """
        del self.__cxx_obj_ptr

    def __len__(self) -> int:
        """
        Number of the constraints.
        """
        return self.__cxx_obj_ptr.get_n_constraints()

    def __getitem__(self, idx: int) -> dict:
        """
        Return a constraint.
        """
        if (idx >= self.__cxx_obj_ptr.get_n_constraints()):
            raise IndexError('constraint index out of range')
        else:
            c = dict()
            c['type'] = self.__cxx_obj_ptr.get_types()[idx]
            c['indices'] = self.__cxx_obj_ptr.get_indices()[idx]
            c['target'] = self.__cxx_obj_ptr.get_targets()[idx]
            c['tolerance'] = self.__cxx_obj_ptr.get_tolerances()[idx]
            return c

    @property
    def max_cycles(self) -> int:
        """
        Maximum number of the RATTLE iterations.
        """
        return self.__cxx_obj_ptr.get_max_cycles()

    @max_cycles.setter
    def max_cycles(self, v) -> int:
        """
        Set maximum number of the RATTLE iterations.
        """
        if (v < 0):
            raise RuntimeError('Maximum number of the RATTLE iterations' + \
                'must be larger than 0!')
        self.__cxx_obj_ptr.set_max_cycles(v)

    @property
    def die_on_fail(self) -> bool:
        """
        If exit when the constraining failed.
        """
        return self.__cxx_obj_ptr.get_die_on_fail()

    @die_on_fail.setter
    def die_on_fail(self, v) -> bool:
        """
        Set if exit when the constraining failed.
        """
        self.__cxx_obj_ptr.set_die_on_fail(bool(v))

    def append(self, c: dict) -> None:
        """
        Append one constraint to the constraints.

        Parameters
        ----------
        c : dict
            The constraint to append. The dictionary should contain following
            fields:
            - 'type' : int
                Type of this constraint: 0 for distance, 1 for angle and
                2 for torsion.
            - 'indices' : List(int)
                Indices of atoms involved in this constraint.
            - 'target' : float
                Target value of the constraint.
            - 'tolerance' : float
                Tolerance value of the RATTLE iterations.
        """
        self.__cxx_obj_ptr.append(c['type'], c['indices'], \
                                  c['target'], c['tolerance'])
        self.__system.groups.update_n_dof()
        self.__system.find_segments()

    def appends(self, c_list: list) -> None:
        """
        Append multiple constraints to the constraints.

        Parameters
        ----------
        c_list : list
            The constraints to append.
        """
        for c in c_list:
            self.__cxx_obj_ptr.append(c['type'], c['indices'], \
                                      c['target'], c['tolerance'])
        self.__system.groups.update_n_dof()
        self.__system.find_segments()

    def pop(self, idx: int = None) -> dict:
        """
        Remove one of the constraints.

        Parameters
        ----------
        idx : int
            Indeices of the constraints to remove.
        """
        if (idx == None):
            idx = self.__cxx_obj_ptr.get_n_constraints() - 1
        if (self.__cxx_obj_ptr.get_n_constraints() <= 0):
            return
        result = self[idx]
        self.__cxx_obj_ptr.pop(idx)
        self.__system.groups.update_n_dof()
        self.__system.find_segments()
        return result

    def clear(self) -> None:
        """
        Remove all the constraints.
        """
        self.__cxx_obj_ptr.clear()

    def rattle_constrain_q(self, \
                           positions: double[:,:], \
                           velocities: double[:,:], \
                           mass: double[:,:], \
                           dt: double) -> None:
        """
        Perform the upper part of RATTLE.

        Parameters
        ----------
        positions : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Positions of the atoms in the simulated system.
        velocities : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Velocities of the atoms in the simulated system.
        mass : np.array(dtype=np.float64)
            dim : n_atoms
            Atomic masses.
        dt : float
            Timestep of the simulation.
        """
        cdef int n_atoms = positions.shape[0]
        if (self.__cxx_obj_ptr.get_n_constraints() != 0):
            self.__cxx_obj_ptr.rattle_constrain_q(&positions[0,0],
                &velocities[0,0], &mass[0,0], dt, n_atoms)

    def rattle_constrain_p(self, \
                           positions: double[:,:], \
                           velocities: double[:,:], \
                           mass: double[:,:], \
                           dt: double) -> None:
        """
        Perform the lower part of RATTLE.

        Parameters
        ----------
        positions : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Positions of the atoms in the simulated system.
        velocities : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Velocities of the atoms in the simulated system.
        mass : np.array(dtype=np.float64)
            dim : n_atoms
            Atomic masses.
        dt : float
            Timestep of the simulation.
        """
        cdef int n_atoms = positions.shape[0]
        if (self.__cxx_obj_ptr.get_n_constraints() != 0):
            self.__cxx_obj_ptr.rattle_constrain_p(&positions[0,0],
                &velocities[0,0], &mass[0,0], dt, n_atoms)
