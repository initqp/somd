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

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "nhc.cxx" nogil:
    pass

cdef extern from "nhc.h" nogil:
    cdef cppclass NHC:
        NHC(double T, double t, int n_bead, int n_d, int n_r) except +
        int get_length()
        int get_n_dof()
        int get_n_respa()
        double get_tau()
        double get_temperature()
        vector[double] get_Q()
        vector[double] get_q()
        vector[double] get_p()
        void set_tau(double t)
        void set_n_dof(int n_d)
        void set_n_respa(int n_r)
        void set_q(vector[double] v)
        void set_p(vector[double] v)
        void set_temperature(double T)
        double propagate(double E_k, double dt)
