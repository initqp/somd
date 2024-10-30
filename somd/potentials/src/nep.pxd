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
from libcpp.string cimport string

cdef extern from "NEP_CPU/src/nep.cpp" nogil:
    pass

cdef extern from "NEP_CPU/src/nep.h" nogil:
    cdef struct ZBL:
        bool enabled
        bool flexibled
        int num_types
        double rc_inner
        double rc_outer
        double para[550]
    cdef struct ParaMB:
        bool use_typewise_cutoff
        bool use_typewise_cutoff_zbl
        double typewise_cutoff_radial_factor
        double typewise_cutoff_angular_factor
        double typewise_cutoff_zbl_factor
        int model_type
        int version
        double rc_radial
        double rc_angular
        double rcinv_radial
        double rcinv_angular
        int n_max_radial
        int n_max_angular
        int L_max
        int dim_angular
        int num_L
        int basis_size_radial
        int basis_size_angular
        int num_types_sq
        int num_c_radial
        int num_types
        double q_scaler[140];
        int atomic_numbers[94];
    cdef struct ANN:
        int dim
        int num_neurons1
        int num_para
        int num_para_ann
    cdef cppclass NEP3:
        NEP3() except +
        NEP3(string potential_filename) except +
        void init_from_file(string potential_filename, bool is_rank_0) except +
        void compute(
            vector[int] types,
            vector[double] box,
            vector[double] positions,
            vector[double] potentials,
            vector[double] forces,
            vector[double] virial
        )
        ZBL zbl
        ANN annmb
        ParaMB paramb
