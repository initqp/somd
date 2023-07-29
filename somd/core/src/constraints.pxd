# cython: language_level=3
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

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "constraints.cxx" nogil:
    pass

cdef extern from "constraints.h" nogil:
    cdef cppclass RATTLE:
        RATTLE() except +
        int get_n_constraints()
        int get_max_cycles()
        bool get_die_on_fail()
        vector[int] get_types()
        vector[double] get_targets()
        vector[double] get_tolerances()
        vector[vector[int]] get_indices()
        void set_max_cycles(int v)
        void set_die_on_fail(bool v)
        void append(int t, vector[int] &idx, double target, \
            double tolerance) except +
        void pop(int idx)
        void clear()
        void rattle_constrain_q(double *positions, double *velocities, \
            double *mass, double dt, int n_atoms) except +
        void rattle_constrain_p(double *positions, double *velocities, \
            double *mass, double dt, int n_atoms) except +
