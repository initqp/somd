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
The trajectory writters.
"""

import os as _os
import time as _t
import h5py as _h5
import numpy as _np
import pickle as _pl
import base64 as _bs
import warnings as _w
from somd import _version
from somd import core as _mdcore
from somd.constants import CONSTANTS as _c
from somd.constants import SOMDDEFAULTS as _d
from . import utils as _utils

__all__ = ['H5READER', 'H5WRITER', 'EXYZWRITER']


class H5WRITER(_utils.POSTSTEPOBJ):
    """
    Write the MDTraj HDF5 trajectory file [1].

    Parameters
    ----------
    file_name : str
        Name of the trajectory file.
    interval : int
        Step interval of writing the trajectory file.
    write_velocities : bool
        If write velocities to the trajectory.
    write_forces : bool
        If write forces to the trajectory.
    write_virial : bool
        If write virial to the trajectory.
    wrap_positions : bool
        If write the wrapped positions to the box when writing the trajectory.
    append : bool
        If append to the old trajectory.
    restart_file : bool
        Require writing a restart file.
    use_double : bool
        Use 64-bit floating point data in the trajectory instead of the
        original 32-bit convention. Note that this option will lead to
        non-standard trajectories and may cause some readers to fail.
    potential_list : List(int)
        Index of the potentials whose forces and virial should be write to
        the force and virial array in the trajectory. Note that the included
        potentials must act on the whole simulated system.

    References
    ----------
    .. [1] https://mdtraj.org/1.9.7/hdf5_format.html
    """

    def __init__(self,
                 file_name: str,
                 interval: int = 1,
                 write_velocities: bool = False,
                 write_forces: bool = False,
                 write_virial: bool = False,
                 wrap_positions: bool = False,
                 append: bool = False,
                 restart_file: bool = False,
                 use_double: bool = False,
                 potential_list: list = None) -> None:
        """
        Create a H5WRITER instance.
        """
        if (restart_file):
            self.__append = False
            self.__is_restart = True
            self.__use_double = True
            self.__write_forces = True
            self.__write_virial = True
            self.__wrap_positions = False
            self.__write_velocities = True
            self.__potential_list = None
        else:
            self.__append = append
            self.__is_restart = False
            self.__use_double = use_double
            self.__write_forces = write_forces
            self.__write_virial = write_virial
            self.__wrap_positions = wrap_positions
            self.__write_velocities = write_velocities
            self.__potential_list = potential_list
        self.__root = None
        self.__file_name = file_name
        super().__init__(interval)

    def __del__(self) -> None:
        """
        Close the file on exit.
        """
        if (self.__root is not None):
            self.__root.close()

    def __dump_attributes(self) -> None:
        """
        Dump attributes to the trajectory file.
        """
        self.__root.attrs['conventionVersion'] = '1.1'
        self.__root.attrs['program'] = 'SOMD'
        self.__root.attrs['programVersion'] = str(_version.get_versions())
        self.__root.attrs['application'] = 'H5WRITER'
        self.__root.attrs['createdTime'] = _t.ctime()
        self.__root.attrs['title'] = self.__file_name
        if (self.__is_restart):
            self.__root.attrs['conventions'] = ['Pande', 'SOMDRESTART']
            self.__root.attrs['randomState'] = 'null'
        else:
            self.__root.attrs['conventions'] = 'Pande'

    def __create_dataset(self,
                         name: str,
                         shape: tuple,
                         max_shape: tuple,
                         data_type: str,
                         unit: str) -> None:
        """
        Create a new data set.

        Parameters
        ----------
        name : str
            Name of the data set.
        shape : tuple
            Initial shape of the data set.
        max_shape : tuple
            Maximum shape of the data set.
        data_type : str
            The data type.
        unit : str
            Unit of the data set.
        """
        self.__root.create_dataset(name, shape, maxshape=max_shape,
                                   dtype=data_type)
        self.__root[name].attrs['units'] = unit

    def __dump_datasets(self) -> None:
        """
        Setup the required datasets.
        """
        if (self.__use_double):
            type_str = '>f8'
        else:
            type_str = '>f4'
        n_atoms = self.__system.n_atoms
        self.__create_dataset('steps', (0,), (None,), 'uint64', '1')
        self.__create_dataset('potentialEnergy', (0,), (None,), type_str,
                              'kilojoule/mole')
        self.__create_dataset('cell_angles', (0, 3), (None, 3), type_str,
                              'degrees')
        self.__create_dataset('cell_lengths', (0, 3), (None, 3), type_str,
                              'nanometers')
        self.__create_dataset('box', (0, 3, 3), (None, 3, 3), type_str,
                              'nanometers')
        self.__create_dataset('coordinates', (0, n_atoms, 3),
                              (None, n_atoms, 3), type_str, 'nanometers')
        if (self.__write_velocities):
            self.__create_dataset('velocities', (0, n_atoms, 3),
                                  (None, n_atoms, 3), type_str,
                                  'nanometers/picosecond')
        if (self.__write_forces):
            self.__create_dataset('forces', (0, n_atoms, 3),
                                  (None, n_atoms, 3), type_str,
                                  'kilojoule/mole/nanometers')
        if (self.__write_virial):
            self.__create_dataset('virial', (0, 3, 3), (None, 3, 3), type_str,
                                  'kilojoule/mole')
        if (self.__is_restart):
            if (len(self.__nhchains) != 0):
                l_nhc = _d.NHCLENGTH
                n_nhc = len(self.__nhchains)
                self.__create_dataset('nhc_positions', (1, n_nhc, l_nhc),
                                      (1, n_nhc, l_nhc), type_str, '1')
                self.__create_dataset('nhc_momentums', (1, n_nhc, l_nhc),
                                      (1, n_nhc, l_nhc), type_str,
                                      'kilojoule/mole*picosecond')

    def bind_integrator(self, integrator: _mdcore.integrators.INTEGRATOR) \
            -> None:
        """
        Bind an integrator.
        """
        super().bind_integrator(integrator)
        self.__system = integrator.system
        self.__nhchains = integrator._nhchains

    def initialize(self) -> None:
        """
        Open the trajectory for writing.
        """
        super().initialize()
        try:
            _os.stat(self.file_name)
        except FileNotFoundError:
            self.__append = False
        if (self.__append):
            self.__root = _h5.File(self.file_name, 'a')
        else:
            _utils.backup(self.file_name)
            self.__root = _h5.File(self.file_name, 'w')
            self.__dump_attributes()
            self.__dump_datasets()
            self.__root.flush()

    def write(self) -> None:
        """
        Write to the trajectory immediately.
        """
        if (not self.initialized):
            message = 'Must initialize the writer before write!'
            raise RuntimeError(message)
        n_atoms = self.__system.n_atoms
        if (self.__is_restart):
            # The step stands for the number of trajectory frame.
            step = 0
            st = _np.random.get_state()
            st_str = _bs.b64encode(_pl.dumps(st)).decode('utf-8')
            self.__root.attrs['randomState'] = st_str
            if (len(self.__nhchains) != 0):
                for i in range(0, len(self.__nhchains)):
                    n = self.__nhchains[i]
                    self.__root['nhc_positions'][0, i, :] = n.positions
                    self.__root['nhc_momentums'][0, i, :] = n.momentums
        else:
            step = self.__root['coordinates'].shape[0]

        # This `steps` stands for simulation timesteps.
        self.__root['steps'].resize((step + 1,))
        self.__root['steps'][step] = self.step
        l = self.__system.lattice
        self.__root['cell_angles'].resize((step + 1, 3))
        self.__root['cell_angles'][step, :] = l[3:6]
        self.__root['cell_lengths'].resize((step + 1, 3))
        self.__root['cell_lengths'][step, :] = l[0:3]
        self.__root['box'].resize((step + 1, 3, 3))
        self.__root['box'][step, :, :] = self.__system.box
        if (self.__wrap_positions):
            p = self.__system.positions_wrapped
        else:
            p = self.__system.positions
        self.__root['coordinates'].resize((step + 1, n_atoms, 3))
        self.__root['coordinates'][step, :, :] = p
        self.__root['potentialEnergy'].resize((step + 1,))
        if (self.__potential_list is None):
            self.__root['potentialEnergy'][step] = \
                self.__system.energy_potential
        else:
            self.__root['potentialEnergy'][step] = 0
            for i in self.__potential_list:
                self.__root['potentialEnergy'][step] += \
                    self.__system.potentials[i].energy_potential
        if (self.__write_velocities):
            self.__root['velocities'].resize((step + 1, n_atoms, 3))
            self.__root['velocities'][step, :, :] = self.__system.velocities
        if (self.__write_forces):
            self.__root['forces'].resize((step + 1, n_atoms, 3))
            if (self.__potential_list is None):
                self.__root['forces'][step, :, :] = self.__system.forces
            else:
                for i in self.__potential_list:
                    atom_list = self.__system.potentials[i].atom_list
                    self.__root['forces'][step, atom_list, :] += \
                        self.__system.potentials[i].forces
        if (self.__write_virial):
            self.__root['virial'].resize((step + 1, 3, 3))
            if (self.__potential_list is None):
                self.__root['virial'][step, :, :] = self.__system.virial
            else:
                for i in self.__potential_list:
                    self.__root['virial'][step, :, :] += \
                        self.__system.potentials[i].virial
        self.__root.flush()

    def update(self) -> None:
        """
        Write the trajectory according to the definded step interval.
        """
        if (super().update()):
            self.write()

    @property
    def file_name(self) -> str:
        """
        File name of this trajectory.
        """
        return self.__file_name

    @property
    def is_restart(self) -> bool:
        """
        If this is a restart file.
        """
        return self.__is_restart

    @property
    def dtype(self) -> type:
        """
        Data type of the trajectory.
        """
        if (self.__use_double):
            return _np.double
        else:
            return _np.float


class EXYZWRITER(_utils.POSTSTEPOBJ):
    """
    Write the extended XYZ trajectory file [1].

    Parameters
    ----------
    file_name : str
        Name of the trajectory file.
    interval : int
        Step interval of writing the trajectory file.
    write_velocities : bool
        If write velocities to the trajectory.
    write_forces : bool
        If write forces to the trajectory.
    write_virial : bool
        If write virial to the trajectory.
    wrap_positions : bool
        If write the wrapped positions to the box when writing the trajectory.
    potential_list : List(int)
        Index of the potentials whose forces and virial should be write to
        the force and virial array in the trajectory. Note that the included
        potentials must act on the whole simulated system.
    append : bool
        If append to the old trajectory.
    format_str : str
        The ascii data format string. E.g.: '{:-24.15e}'.
    energy_shift : float
        Shift the total energy by this value before recording the total energy
        to the trajectory. In unit of (kJ/mol).

    Notes
    -----
    The main aim of this module is to generate the training data sets of a NEP
    potential. Thus we do not really need a reading function.

    References
    ----------
    .. [1] https://github.com/libAtoms/extxyz
    """

    def __init__(self,
                 file_name: str,
                 interval: int = 1,
                 write_velocities: bool = True,
                 write_forces: bool = True,
                 write_virial: bool = True,
                 wrap_positions: bool = False,
                 format_str: str = '{:-24.15e}',
                 append: bool = False,
                 potential_list: list = None,
                 energy_shift: float = 0.0) -> None:
        """
        Create an EXYZWRITER instance.
        """
        self.__fp = None
        self.__append = append
        self.__file_name = file_name
        self.__write_forces = write_forces
        self.__write_virial = write_virial
        self.__write_velocities = write_velocities
        self.__wrap_positions = wrap_positions
        self.__potential_list = potential_list
        self.__conversion = _c.AVOGACONST * _c.ELECTCONST
        if (energy_shift is None):
            energy_shift = 0.0
        else:
            energy_shift = energy_shift / self.__conversion * 1000
        self.__energy_shift = _np.array([energy_shift], dtype=_np.double)
        _np.set_printoptions(formatter={'float': format_str.format},
                             linewidth=10000)
        super().__init__(interval)

    def __del__(self) -> None:
        """
        Close the file on exit.
        """
        if (self.__fp is not None):
            self.__fp.close()

    def __convert_data(self):
        """
        Convert units of system data.
        """
        if (self.__potential_list is None):
            self.__energy_potential[0] = self.__system.energy_potential / \
                self.__conversion * 1000
            self.__energy_potential[0] -= self.__energy_shift[0]
            if (self.__write_forces):
                self.__forces[:] = \
                    self.__system.forces[:] / self.__conversion * 100
            if (self.__write_virial):
                self.__virial[:] = \
                    self.__system.virial[:] / self.__conversion * 1000
        else:
            self.__energy_potential[0] = 0.0
            for i in self.__potential_list:
                self.__energy_potential[0] += \
                    self.__system.potentials[i].energy_potential
            self.__energy_potential[0] /= self.__conversion * 0.001
            self.__energy_potential[0] -= self.__energy_shift[0]
            if (self.__write_forces):
                self.__forces[:] = 0
                for i in self.__potential_list:
                    atom_list = self.__system.potentials[i].atom_list
                    self.__forces[atom_list] += \
                        self.__system.potentials[i].forces[:]
                self.__forces[:] /= self.__conversion * 0.01
            if (self.__write_virial):
                self.__virial[:] = 0
                for i in self.__potential_list:
                    self.__virial[:] += self.__system.potentials[i].virial[:]
                self.__virial[:] /= self.__conversion * 0.001

    def bind_integrator(self, integrator: _mdcore.integrators.INTEGRATOR) \
            -> None:
        """
        Bind an integrator.
        """
        super().bind_integrator(integrator)
        self.__system = integrator.system
        self.__energy_potential = _np.zeros((1,), dtype=_np.double)
        if (self.__write_forces):
            self.__forces = _np.zeros((self.__system.n_atoms, 3), _np.double)
        if (self.__write_virial):
            self.__virial = _np.zeros((3, 3), _np.double)

    def initialize(self) -> None:
        """
        Open the trajectory for writing.
        """
        super().initialize()
        try:
            _os.stat(self.file_name)
        except FileNotFoundError:
            self.__append = False
        if (self.__append):
            self.__fp = open(self.file_name, 'a')
        else:
            _utils.backup(self.file_name)
            self.__fp = open(self.file_name, 'w')

    def write(self) -> None:
        """
        Write to the trajectory immediately.
        """
        if (not self.initialized):
            message = 'Must initialize the writer before write!'
            raise RuntimeError(message)
        self.__convert_data()
        print(self.__system.n_atoms, file=self.__fp)
        s = format(self.__energy_potential)[1:-1]
        header = 'energy={} '.format(s.strip())
        s = format(self.__energy_shift)[1:-1]
        header += 'energy_shift={} pbc="T T T" '.format(s.strip())
        if (self.__write_virial):
            s = format(self.__virial.reshape(-1))[1:-1]
            header += 'virial="{}" '.format(s)
        s = format(self.__system.box.reshape(-1) * 10)[1:-1]
        header += 'Lattice="{}" '.format(s)
        header += 'Properties=species:S:1:pos:R:3'
        if (self.__write_velocities):
            header += ':vel:R:3'
        if (self.__write_forces):
            header += ':force:R:3'
        print(header, file=self.__fp)
        if (self.__wrap_positions):
            p = self.__system.positions_wrapped
        else:
            p = self.__system.positions
        for i in range(0, self.__system.n_atoms):
            print(self.__system.atomic_symbols[i] + ' ',
                  file=self.__fp, end='')
            print(format(p[i] * 10)[1:-1], file=self.__fp, end='')
            if (self.__write_velocities):
                print(format(self.__system.velocities[i] * 100)[1:-1],
                      file=self.__fp, end='')
            if (self.__write_forces):
                print(format(self.__forces[i])[1:-1], file=self.__fp, end='')
            print('', file=self.__fp)
        self.__fp.flush()

    def update(self) -> None:
        """
        Write the trajectory according to the definded step interval.
        """
        if (super().update()):
            self.write()

    @property
    def file_name(self) -> str:
        """
        File name of this trajectory.
        """
        return self.__file_name


class H5READER(object):
    """
    Read the MDTraj HDF5 trajectory file [1].

    Parameters
    ----------
    file_name : str
        Name of the trajectory file.
    read_coordinates : bool
        If read coordinates from the file.
    read_velocities : bool
        If read velocities from the file.
    read_forces : bool
        If read forces from the file.
    read_cell : bool
        If read cell data from the file.
    read_nhc_data : bool
        If read Nose-Hoover Chains data from the file.
    read_rng_state : bool
        If read the state string of the RNG from this file.

    References
    ----------
    .. [1] https://mdtraj.org/1.9.7/hdf5_format.html
    """

    def __init__(self,
                 file_name: str,
                 read_coordinates: bool = True,
                 read_velocities: bool = True,
                 read_forces: bool = True,
                 read_cell: bool = True,
                 read_nhc_data: bool = True,
                 read_rng_state: bool = True) -> None:
        """
        Create a H5READER instance.
        """
        self.__file_name = file_name
        self.__read_cell = read_cell
        self.__read_forces = read_forces
        self.__read_nhc_data = read_nhc_data
        self.__read_rng_state = read_rng_state
        self.__read_velocities = read_velocities
        self.__read_coordinates = read_coordinates
        self.__root = _h5.File(file_name, 'r')
        self.__check_attributes()
        self.__n_atoms = self.__root['coordinates'].shape[1]
        if ('SOMDRESTART' in str(self.__root.attrs['conventions'])):
            self.__is_restart = True
        else:
            self.__is_restart = False

    def __del__(self) -> None:
        """
        Close the file on exit.
        """
        if (self.__root is not None):
            self.__root.close()

    def __check_attribute(self,
                          name: str,
                          value=None,
                          die_on_fail: bool = True) -> None:
        """
        Check if a given attribute exists.

        Parameters
        ----------
        name : str
            The attribute to check.
        value : str
            Expected value of the attribute.
        die_on_fail : bool
            If raise an error when the attribute does not match with the
            expected value.
        """
        attr = ''
        try:
            attr = self.__root.attrs[name]
        except:
            message = 'Attribute {} does not exist in file {}'
            message.format(name, self.__file_name)
            if (die_on_fail):
                raise AttributeError(message)
            else:
                _w.warn(message)
        if (value is not None and value not in str(attr)):
            message = ('Values of attribute {} in file {} do not contain ' +
                       'the expected value {}')
            message = message.format(name, self.__file_name, value)
            if (die_on_fail):
                raise AttributeError(message)
            else:
                _w.warn(message)

    def __check_attributes(self) -> None:
        """
        Check all required attribute of a valid HDF5 trajectory.
        """
        self.__check_attribute('conventions', 'Pande')
        self.__check_attribute('conventionVersion', '1.1', die_on_fail=False)
        self.__check_attribute('program')
        self.__check_attribute('programVersion')
        self.__check_attribute('title', die_on_fail=False)
        self.__check_attribute('application', die_on_fail=False)

    def bind_integrator(self, integrator: _mdcore.integrators.INTEGRATOR) \
            -> None:
        """
        Bind an integrator.
        """
        if (integrator.system is None):
            raise RuntimeError('Integrator has not bind a system!')
        if (self.__n_atoms != integrator.system.n_atoms):
            raise RuntimeError('The atom number in the trajectory file not ' +
                               'match atom number of the simulated system!')
        self.__integrator = integrator
        self.__system = integrator.system
        self.__nhchains = integrator._nhchains

    def read(self, frame_index: int = -1) -> None:
        """
        Read data from the trajectory.

        Parameters
        ----------
        frame_index : int
            The index of the frame to read data from. If the parameter is
            negative, last frame of the trajectory will be read.
        """
        if (self.__integrator is None):
            message = 'Must bind to an integrator before read!'
            raise RuntimeError(message)
        n_frames = self.__root['coordinates'].shape[0]
        if (frame_index < 0):
            frame_index = n_frames - 1
        elif (self.__is_restart and frame_index > 0):
            raise RuntimeError('Could not load data from a non-zero frame ' +
                               'of a restart file!')
        elif (not self.__is_restart) and (n_frames <= frame_index):
            message = 'Could not load data from frame {:d} of a ' + \
                      'trajectory with only {:d} frmaes!'
            message = message.format(frame_index,
                                     self.__root['coordinates'].shape[0])
            raise RuntimeError(message)
        if (self.__is_restart):
            if (self.__read_rng_state):
                st_str = self.__root.attrs['randomState']
                _np.random.set_state(_pl.loads(_bs.b64decode(st_str)))
            if (self.__read_nhc_data):
                if (len(self.__nhchains) == 0):
                    _w.warn('NHCHAINS have not been bound to the ' +
                            'reader! NHCHAINS data will not be load!')
                elif ('nhc_positions' not in self.__root.keys()):
                    _w.warn('Can not load NHCHAINS data from a restart ' +
                            'file without NHCHAINS data!')
                else:
                    n_nhc = min(len(self.__nhchains),
                                self.__root['nhc_positions'].shape[1])
                    for i in range(0, n_nhc):
                        self.__nhchains[i].positions = \
                            self.__root['nhc_positions'][0, i, :]
                        self.__nhchains[i].momentums = \
                            self.__root['nhc_momentums'][0, i, :]
        else:
            if (self.__read_rng_state):
                _w.warn('Can not load RNG state from a non-restart file!')
            if (self.__read_nhc_data):
                _w.warn('Can not load NHCHAINS data from a non-restart file!')
        if (self.__read_coordinates):
            try:
                self.__system.positions[:, :] = \
                    self.__root['coordinates'][frame_index]
            except:
                _w.warn('Can not load coordinates form ' + self.file_name)
        if (self.__read_velocities):
            try:
                self.__system.velocities[:, :] = \
                    self.__root['velocities'][frame_index]
            except:
                _w.warn('Can not load velocities form ' + self.file_name)
        if (self.__read_forces):
            try:
                self.__system.forces[:, :] = self.__root['forces'][frame_index]
            except:
                _w.warn('Can not load forces form ' + self.file_name)
        if (self.__read_cell):
            if (str(self.__root.attrs['program']) == 'SOMD'):
                self.__system.box[:] = self.__root['box'][frame_index][:]
            else:
                try:
                    tmp = _np.zeros(6, dtype=_np.double)
                    tmp[3:6] = self.__root['cell_angles'][frame_index][:]
                    tmp[0:3] = self.__root['cell_lengths'][frame_index][:]
                    self.__system.lattice = tmp
                except:
                    _w.warn('Can not load cell data form ' + self.file_name)

    @property
    def root(self) -> _h5.File:
        """
        The HDF5 file root.
        """
        return self.__root

    @property
    def file_name(self) -> str:
        """
        File name of this trajectory.
        """
        return self.__file_name

    @property
    def is_restart(self) -> bool:
        """
        If this is a restart file.
        """
        return self.__is_restart

    @property
    def integrator(self) -> _mdcore.integrators.INTEGRATOR:
        """
        Refernce to the integrator.
        """
        return self.__integrator

    @property
    def n_frames(self) -> int:
        """
        Number of frames in the trajectory.
        """
        return self.__root['coordinates'].shape[0]
