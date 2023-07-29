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
from nhc cimport NHC

__all__ = ['NHCHAINS']

cdef class NHCHAINS(object):
    """
    Define Nose-Hoover Chains [1].

    References
    ----------
    .. [1] Martyna, Glenn J., et al. "Explicit reversible integrators for
           extended systems dynamics." Molecular Physics 87.5 (1996): 1117-1157.
    """

    # Pointer to the internal C++ object.
    cdef NHC *__cxx_obj_ptr

    def __cinit__(self, \
                  temperature: double, \
                  tau: double, \
                  n_bead: int, \
                  n_dof: int, \
                  n_respa: int) -> None:
        """
        Initialize the Nose-Hoover chains.

        Parameters
        ----------
        temperature : double
            Target temperature of the chains. In unit of (K).
        tau : double
            Relaxation timescale of the chains. In unit of (ps).
        n_bead : int
            Length of the chains.
        n_dof : int
            Number of DOFs of the simulated system.
        n_respa : int
            Number of the RESPA loops.
        """
        self.__cxx_obj_ptr = new NHC(temperature, tau, n_bead, n_dof, n_respa)

    def __dealloc__(self) -> None:
        """ Finalize the internal data. """
        del self.__cxx_obj_ptr

    @property
    def length(self) -> int:
        """ Length of the chains. """
        return self.__cxx_obj_ptr.get_length()

    @property
    def masses(self) -> list:
        """ Masses of the chains """
        return self.__cxx_obj_ptr.get_Q()

    @property
    def positions(self) -> list:
        """ Positions of the chains. In unit of (1). """
        return self.__cxx_obj_ptr.get_q()

    @positions.setter
    def positions(self, v: vector[double]) -> None:
        """
        Set positions of the chains.
        """
        self.__cxx_obj_ptr.set_q(v)

    @property
    def momentums(self) -> list:
        """ Momentums of the chains. In unit of (kJ/mol*ps). """
        return self.__cxx_obj_ptr.get_p()

    @momentums.setter
    def momentums(self, v: vector[double]) -> None:
        """
        Set momentums of the chains.
        """
        self.__cxx_obj_ptr.set_p(v)

    @property
    def n_respa(self) -> int:
        """ Number of RESPA loops. """
        return self.__cxx_obj_ptr.get_n_respa()

    @property
    def n_dof(self) -> int:
        """ Number of DOFs of the simulated system. """
        return self.__cxx_obj_ptr.get_n_dof()

    @n_dof.setter
    def n_dof(self, n_d : int) -> None:
        """
        Set number of DOFs of the simulated system.
        """
        self.__cxx_obj_ptr.set_n_dof(n_d)

    @property
    def tau(self) -> double:
        """ Relaxation timescale of the chains. In unit of (ps). """
        return self.__cxx_obj_ptr.get_tau()

    @tau.setter
    def tau(self, t : double) -> None:
        """
        Set the relaxation timescale of the chains.
        """
        self.__cxx_obj_ptr.set_tau(t)

    @property
    def temperature(self) -> double:
        """ Target temperature of the chains. In unit of (K). """
        return self.__cxx_obj_ptr.get_temperature()

    @temperature.setter
    def temperature(self, T: double) -> None:
        """
        Set target temperature of the chains.
        """
        self.__cxx_obj_ptr.set_temperature(T)

    def propagate(self, E_k: double, dt: double) -> double:
        """
        Propagate the chains by one timestep.

        Parameters
        ----------
        E_k : float
            The kinetic energy of the simulated system. In unit of (kJ/mol).
        dt : float
            The timestep. In unit of (ps).
        """
        return self.__cxx_obj_ptr.propagate(E_k, dt)
