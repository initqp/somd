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
from libcpp.string cimport string
from nep cimport NEP3

__all__ = ['NEPWRAPPER']

cdef class NEPWRAPPER(object):
    """
    Wrapper of the neuroevolution potential, version 3 [1].

    References
    ----------
    .. [1] Fan, Zheyong, et al. "Neuroevolution machine learning potentials:
           Combining high accuracy and low cost in atomistic simulations and
           application to heat transport." Physical Review B 104.10 (2021):
           104309.
    """

    # Internal C++ object.
    cdef int __n_atoms
    cdef NEP3 *__cxx_obj_ptr
    cdef string __file_name
    cdef vector[int] __cxx_type_vec
    cdef vector[double] __cxx_box_vec
    cdef vector[double] __cxx_force_vec
    cdef vector[double] __cxx_energy_vec
    cdef vector[double] __cxx_virial_vec
    cdef vector[double] __cxx_position_vec

    def __cinit__(self, file_name: str, atomic_symbols: list(str)) -> None:
        """
        Initialize the potential calculator.

        Parameters
        ----------
        file_name : str
            Name of the potential data file.
        atomic_symbols : list(str)
            List of the atomic symbols.
        """
        self.__n_atoms = len(atomic_symbols)
        self.__file_name = bytes(file_name, encoding='UTF-8')
        self.__cxx_obj_ptr = new NEP3()
        self.__cxx_box_vec = vector[double](9)
        self.__cxx_force_vec = vector[double](self.__n_atoms * 3)
        self.__cxx_energy_vec = vector[double](self.__n_atoms)
        self.__cxx_virial_vec = vector[double](self.__n_atoms * 9)
        self.__cxx_position_vec = vector[double](self.__n_atoms * 3)
        self.__cxx_obj_ptr.init_from_file(self.__file_name, False)
        self.__init_type_vec(atomic_symbols)

    def __dealloc__(self) -> None:
        """ Finalize the internal data. """
        del self.__cxx_obj_ptr
        self.__cxx_obj_ptr = NULL

    cdef __init_type_vec(self, atomic_symbols: list(str)):
        """
        Initialize the atomic type vector.

        Parameters
        ----------
        atomic_symbols : list(str)
            List of the atomic symbols.
        """
        self.__cxx_type_vec = vector[int](self.__n_atoms)
        fp = open(self.__file_name, 'r')
        header = fp.readline().split()
        fp.close()
        for i in range(0, self.__n_atoms):
            is_valid = False
            for j in range(0, int(header[1])):
                if (atomic_symbols[i].lower() == header[(j + 2)].lower()):
                    self.__cxx_type_vec[i] = j
                    is_valid = True
            if (not is_valid):
                raise RuntimeError('Unknown element: ' + atomic_symbols[i])

    cdef __copy_input_data(self, positions: double[:, :], box: double[:, :]):
       """ Copy and transpose the position and box data. """
       cdef int i = 0
       for i in range(0, self.__n_atoms):
           self.__cxx_position_vec[i + self.__n_atoms * 0] = positions[i, 0]
           self.__cxx_position_vec[i + self.__n_atoms * 1] = positions[i, 1]
           self.__cxx_position_vec[i + self.__n_atoms * 2] = positions[i, 2]
       self.__cxx_box_vec[0] = box[0, 0]
       self.__cxx_box_vec[1] = box[1, 0]
       self.__cxx_box_vec[2] = box[2, 0]
       self.__cxx_box_vec[3] = box[0, 1]
       self.__cxx_box_vec[4] = box[1, 1]
       self.__cxx_box_vec[5] = box[2, 1]
       self.__cxx_box_vec[6] = box[0, 2]
       self.__cxx_box_vec[7] = box[1, 2]
       self.__cxx_box_vec[8] = box[2, 2]

    cdef __copy_output_data(self, forces: double[:, :], virial: double[:, :]):
       """ Calculate and return the total energy, force and virial data. """
       cdef int i = 0
       cdef double energy = 0
       for i in range(0, self.__n_atoms):
           energy += self.__cxx_energy_vec[i]
           forces[i, 0] = self.__cxx_force_vec[i + self.__n_atoms * 0]
           forces[i, 1] = self.__cxx_force_vec[i + self.__n_atoms * 1]
           forces[i, 2] = self.__cxx_force_vec[i + self.__n_atoms * 2]
           virial[0, 0] += self.__cxx_virial_vec[i + self.__n_atoms * 0]
           virial[0, 1] += self.__cxx_virial_vec[i + self.__n_atoms * 1]
           virial[0, 2] += self.__cxx_virial_vec[i + self.__n_atoms * 2]
           virial[1, 0] += self.__cxx_virial_vec[i + self.__n_atoms * 3]
           virial[1, 1] += self.__cxx_virial_vec[i + self.__n_atoms * 4]
           virial[1, 2] += self.__cxx_virial_vec[i + self.__n_atoms * 5]
           virial[2, 0] += self.__cxx_virial_vec[i + self.__n_atoms * 6]
           virial[2, 1] += self.__cxx_virial_vec[i + self.__n_atoms * 7]
           virial[2, 2] += self.__cxx_virial_vec[i + self.__n_atoms * 8]
       return energy

    @property
    def file_name(self) -> str:
        """ Name of the potential file. """
        return self.__file_name

    @property
    def atomic_types(self) -> list:
        """ NEP types of each atom. """
        return self.__cxx_type_vec

    def calculate(self,
                  box: double[:, :],
                  positions: double[:, :],
                  forces: double[:, :],
                  virial: double[:, :]) -> double:
        """
        Perform the NEP calculation.

        Parameters
        ----------
        box : np.array(dtype=np.double)
            dim : 3 * 3
            The cell vectors. In units of (Ang).
        positions : np.array(dtype=np.double)
            dim : n_atoms * 3
            Positions of each atom. In unit of (Ang).
        forces : np.array(dtype=np.double)
            dim : n_atoms * 3
            Total forces of each atom. In unit of (eV/Ang).
        virial : np.array(dtype=np.double)
            dim : 3 * 3
            The total virial tensor. In unit of (eV).
        """
        self.__copy_input_data(positions, box)
        self.__cxx_obj_ptr.compute(self.__cxx_type_vec,
                                   self.__cxx_box_vec,
                                   self.__cxx_position_vec,
                                   self.__cxx_energy_vec,
                                   self.__cxx_force_vec,
                                   self.__cxx_virial_vec)
        return self.__copy_output_data(forces, virial)

    def dealloc(self) -> None:
        """ Finalize the internal data. """
        del self.__cxx_obj_ptr
        self.__cxx_obj_ptr = NULL
