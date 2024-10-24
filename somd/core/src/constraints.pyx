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

from cython cimport double
from libcpp.vector cimport vector
from constraints cimport RATTLE
import typing as _tp
from .snapshots import SNAPSHOT as _SNAPSHOT

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
    snapshot : somd.core.snapshots.SNAPSHOT
        The snapshot of the simulated.
    handler : typing.Callable
        The handler function for updating system state.
    """

    # Pointer to the internal C++ object.
    cdef object __handler
    cdef object __snapshot
    cdef RATTLE *__cxx_obj_ptr

    def __cinit__(self, snapshot: _SNAPSHOT, handler: _tp.Callable) -> None:
        """
        Initialize the internal data.
        """
        self.__handler = handler
        self.__snapshot = snapshot
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

    def summary(self) -> str:
        """
        Show information about the constraints.
        """
        result = 'CONSTRAINTS\n'
        result += '┣━ n_constraints: {}\n'.format(len(self))
        result += '┣━ max_cycles: {}\n'.format(self.max_cycles)
        result += '┣━ die_on_fail: {}\n'.format(self.die_on_fail)
        for i, c in enumerate(self):
            result += '┣━ CONSTRAINT {}: '.format(i)
            result += 'type {}, '.format(c['type'])
            result += 'target: {}, '.format(c['target'])
            result += 'tolerance: {}, '.format(c['tolerance'])
            result += 'indices: {}\n'.format(c['indices'])
        result += '┗━ END'
        return result

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
            message = (
                'Maximum number of the RATTLE ' +
                'iterations must be larger than 0!'
            )
            raise RuntimeError(message)
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
            - 'indices' : List[int]
                Indices of atoms involved in this constraint.
            - 'target' : float
                Target value of the constraint.
            - 'tolerance' : float
                Tolerance value of the RATTLE iterations.
        """
        self.__cxx_obj_ptr.append(
            c['type'],
            c['indices'],
            c['target'],
            c['tolerance']
        )
        self.__handler()

    def appends(self, c_list: list) -> None:
        """
        Append multiple constraints to the constraints.

        Parameters
        ----------
        c_list : list
            The constraints to append.
        """
        for c in c_list:
            self.__cxx_obj_ptr.append(
                c['type'],
                c['indices'],
                c['target'],
                c['tolerance']
            )
        self.__handler()

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
        self.__handler()
        return result

    def clear(self) -> None:
        """
        Remove all the constraints.
        """
        self.__cxx_obj_ptr.clear()
        self.__handler()

    cdef void rattle_constrain_q_wrapper(
        self,
        positions: double[:, :],
        velocities: double[:, :],
        mass: double[:, :],
        dt: double
    ):
        """
        Perform the upper part of RATTLE (wrapper).

        Parameters
        ----------
        positions : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Positions of the atoms in the simulated system. In unit of (nm).
        velocities : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Velocities of the atoms in the simulated system. In unit of
            (nm/ps).
        mass : np.array(dtype=np.float64)
            dim : n_atoms
            Atomic masses. In unit of (g/mol).
        dt : float
            Timestep of the simulation. In unit of (ps).
        """
        cdef int n_atoms = positions.shape[0]
        if (self.__cxx_obj_ptr.get_n_constraints() != 0):
            self.__cxx_obj_ptr.rattle_constrain_q(
                &positions[0,0],
                &velocities[0,0],
                &mass[0,0],
                dt,
                n_atoms
            )

    cdef void rattle_constrain_p_wrapper(
        self,
        positions: double[:, :],
        velocities: double[:, :],
        mass: double[:, :],
        dt: double
    ):
        """
        Perform the lower part of RATTLE (wrapper).

        Parameters
        ----------
        positions : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Positions of the atoms in the simulated system. In unit of (nm).
        velocities : np.array(dtype=np.float64)
            dim : n_atoms * 3
            Velocities of the atoms in the simulated system. In unit of
            (nm/ps).
        mass : np.array(dtype=np.float64)
            dim : n_atoms
            Atomic masses. In unit of (g/mol).
        dt : float
            Timestep of the simulation. In unit of (ps).
        """
        cdef int n_atoms = positions.shape[0]
        if (self.__cxx_obj_ptr.get_n_constraints() != 0):
            self.__cxx_obj_ptr.rattle_constrain_p(
                &positions[0,0],
                &velocities[0,0],
                &mass[0,0],
                dt,
                n_atoms
            )

    def rattle_constrain_q(self, dt: double) -> None:
        """
        Perform the upper part of RATTLE.

        Parameters
        ----------
        dt : float
            Timestep of the simulation. In unit of (ps).
        """
        self.rattle_constrain_q_wrapper(
            self.__snapshot.positions,
            self.__snapshot.velocities,
            self.__snapshot.masses,
            dt
        )

    def rattle_constrain_p(self, dt: double) -> None:
        """
        Perform the lower part of RATTLE.

        Parameters
        ----------
        dt : float
            Timestep of the simulation. In unit of (ps).
        """
        self.rattle_constrain_p_wrapper(
            self.__snapshot.positions,
            self.__snapshot.velocities,
            self.__snapshot.masses,
            dt
        )
