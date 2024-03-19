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
