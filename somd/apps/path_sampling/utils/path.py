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
Internal represention of a path. Most operations are lazy.
"""

import gc as _gc
import numpy as _np
import typing as _tp
import warnings as _w
from somd import core as _mdcore
from somd.potentials import PLUMED as _PLUMED
from ...trajectories import H5READER as _H5READER
from ...trajectories import H5WRITER as _H5WRITER

__all__ = ['PATH', 'concatenate']


class PATH(object):
    """
    The internal represention of a path, implemented by mapping frames in
    HDF5 trajectory files to path frames. Instances of this class are
    immutable.

    Parameters
    ----------
    reference_system : somd.core.systems.MDSYSTEM
        The simulated system.
    trajectories : List[str]
        Name of the trajectory files.
    timestep : float
        Timestep between two frames. In unit of (ps).
    frame_mapping : List[List[int]]
        The mapping array that records the index of one trajectory frame
        in the path. If an empty list was given, all frames in the
        trajectory will be used.
    reverse_velocities : List[bool]
        If directions of velocities should be reversed when visited.
    cacheable : bool
        If cache the visited frames to avoid reading. Note that for large
        trajectories, this would be memory inefficient.
    """

    def __init__(self,
                 reference_system: _mdcore.systems.MDSYSTEM,
                 trajectories: list,
                 timestep: float,
                 frame_mapping: list = None,
                 reverse_velocities: list = None,
                 cacheable: bool = False) -> None:
        """
        Create a PATH instance.
        """
        self.__timestep = timestep
        self.__trajectories = trajectories
        self.__reference_system = reference_system
        if (reverse_velocities is not None):
            self.__reverse_velocities = reverse_velocities
        else:
            self.__reverse_velocities = [False] * len(trajectories)
        if (frame_mapping is not None):
            mapping = frame_mapping.copy()
        else:
            mapping = [[]] * len(trajectories)
        self.__check_trajectories(mapping)
        self.__make_index(mapping)
        self.__frame_cache = None
        self.cacheable = cacheable

    def __len__(self) -> int:
        """
        Length of the path.
        """
        return self.__length

    def __add__(self, path: 'PATH') -> 'PATH':
        """
        Concatenate two paths.

        Notes
        -----
        This method does not check if the two path are based on the same
        system. Only if the two systems contain the same atom number, this
        method will success.
        """
        system_self = self.reference_system
        system_path = path.reference_system
        if (system_self.n_atoms != system_path.n_atoms):
            message = 'Can not add two paths with different atom numbers!'
            raise RuntimeError(message)
        if (system_self.atomic_symbols != system_path.atomic_symbols):
            message = 'Atomic types of the two paths are different!' + \
                      'Will use the atomic types from the first path.'
            _w.warn(message)
        trajectories = [*self.__trajectories, *path.trajectories]
        mapping_self = [[j[1] for j in self.mapping if (j[0] == i)]
                        for i in range(0, len(self.trajectories))]
        mapping_path = [[j[1] for j in path.mapping if (j[0] == i)]
                        for i in range(0, len(path.trajectories))]
        mapping = [*mapping_self, *mapping_path]
        reverse_velocities = [*self.__reverse_velocities,
                              *path.reverse_velocities]
        return PATH(self.reference_system, trajectories, self.timestep,
                    mapping, reverse_velocities, self.cacheable)

    def __getitem__(self, index: _tp.Union[int, slice]) -> 'PATH':
        """
        Build partial path from the path.
        """
        if isinstance(index, slice):
            new_mapping = self.__mapping[index]
            indices = list(set([i[0] for i in new_mapping]))
            mapping = [[j[1] for j in new_mapping if (j[0] == i)]
                       for i in indices]
            trajectories = [self.__trajectories[i] for i in indices]
            reverse_velocities = [self.__reverse_velocities[i]
                                  for i in indices]
        elif isinstance(index, int):
            if (index >= len(self)):
                message = 'Frame index {:d} out of path length {:d}!'
                raise IndexError(message.format(index, len(self)))
            mapping = [[self.__mapping[index][1]]]
            trajectories = [self.__trajectories[self.__mapping[index][0]]]
            reverse_velocities = \
                [self.__reverse_velocities[self.__mapping[index][0]]]
        else:
            message = 'Could only indexing a path through int or slice!'
            raise RuntimeError(message)
        return PATH(self.reference_system, trajectories, self.timestep,
                    mapping, reverse_velocities, self.cacheable)

    def __reversed__(self) -> 'PATH':
        """
        Reverse the path.
        """
        mapping = [list(reversed([j[1] for j in self.mapping if (j[0] == i)]))
                   for i in reversed(range(0, len(self.__trajectories)))]
        trajectories = list(reversed(self.__trajectories))
        reverse_velocities = list(reversed(self.__reverse_velocities))
        return PATH(self.reference_system, trajectories, self.timestep,
                    mapping, reverse_velocities, self.cacheable)

    def __check_trajectories(self, frame_mapping: list) -> None:
        """
        Check the trajectories and the mapping matrix.

        Parameters
        ----------
        frame_mapping : List[List[int]]
            The mapping array that records the index of one trajectory frame
            in the path. If an empty list was given, all frames in the
            trajectory will be used.
        """
        self.__length = 0
        self.__readers = []
        if (len(self.__trajectories) != len(frame_mapping)):
            message = 'Mismatch between the number of trajectories and ' + \
                      'the number of mapping matrices!'
            raise RuntimeError(message)
        if (len(self.__trajectories) != len(self.__reverse_velocities)):
            message = 'Mismatch between the number of trajectories and ' + \
                      'the length of the reverse_velocities attribute!'
            raise RuntimeError(message)
        for i, t in enumerate(self.__trajectories):
            reader = _H5READER(t)
            if (len(frame_mapping[i]) == 0):
                frame_mapping[i] = list(range(0, reader.n_frames))
            elif (max(frame_mapping[i]) > reader.n_frames):
                message = 'Invalid frame index {:d} of trajectory {:s}!'
                raise RuntimeError(message.format(max(frame_mapping[i]), t))
            if ('velocities' not in reader.root.keys()):
                message = 'Can not find velocities in trajectory {:s}!'
                _w.warn(message.format(t))
                self.__has_velocities = False
            else:
                self.__has_velocities = True
            self.__length += len(frame_mapping[i])
            self.__readers.append(
                _H5READER(t, read_velocities=self.__has_velocities,
                          read_forces=False, read_nhc_data=False,
                          read_rng_state=False))
            del reader

    def __make_index(self, frame_mapping: list) -> None:
        """
        Index the frames.

        Parameters
        ----------
        frame_mapping : List[List[int]]
            The mapping array that records the index of one trajectory frame
            in the path. If an empty list was given, all frames in the
            trajectory will be used.
        """
        path_index = 0
        self.__mapping = [None] * self.__length
        for trajectory_index, mapping in enumerate(frame_mapping):
            for frame_index in mapping:
                if (self.__mapping[path_index] is not None):
                    message = 'Multiple defination of path frame {:d}!'
                    raise RuntimeError(message.format(path_index))
                self.__mapping[path_index] = (trajectory_index, frame_index)
                path_index += 1

    def frame(self, index: int) -> _mdcore.systems.SNAPSHOT:
        """
        Get one frame from the path.
        """
        if (index >= self.__length):
            message = 'Frame index {:d} out of path length {:d}!'
            raise IndexError(message.format(index, self.__length))
        if (self.__frame_cache is None or self.__frame_cache[index] is None):
            frame_index = self.__mapping[index][1]
            reader = self.__readers[self.__mapping[index][0]]
            snapshot = reader.read_as_snapshot(frame_index)
            if (self.__reverse_velocities[self.__mapping[index][0]]):
                snapshot.velocities[:] *= -1
            if (self.__frame_cache is not None):
                self.__frame_cache[index] = snapshot
        else:
            snapshot = self.__frame_cache[index].copy()
        return snapshot

    def dump(self,
             file_name: str = 'path.h5',
             use_double: bool = False) -> None:
        """
        Dump the path to a HDF5 trajectory file.
        """
        system = self.__reference_system.copy()
        writer = _H5WRITER(file_name, write_velocities=self.__has_velocities,
                           use_double=use_double)
        handler = _mdcore.integrators.vv_integrator(self.__timestep)
        handler.bind_system(system)
        writer.bind_integrator(handler)
        writer.initialize()
        for i in range(0, self.__length):
            system.snapshot = self.frame(i)
            writer.update()
            handler.step += 1
        del writer

    def calculate_cv_values(self, plumed_file: str, cv_names: list):
        """
        Calculate cv values of this path.

        Parameters
        ----------
        plumed_file : str
            Name of the plumed input file.
        cv_names : List[dict]
            Names and components of the collective variables. For example:
            cv_names = [{'d1': 'x'}, {'d1': 'y'}, {'d2': ''}]
        """
        cv_values = []
        atom_list = list(range(0, self.reference_system.n_atoms))
        plumed = _PLUMED(atom_list, plumed_file, timestep=self.timestep,
                         cv_names=cv_names, output_prefix='')
        system = self.reference_system.copy()
        for i in range(0, len(self)):
            system.snapshot = self.frame(i)
            plumed.update(system)
            cv_values.append([plumed.cv_values])
        return _np.concatenate(cv_values)

    @property
    def has_velocities(self) -> bool:
        """
        If this path contains velocities.
        """
        return self.__has_velocities

    @property
    def mapping(self) -> list:
        """
        The frame index. Each element in this list contains index of the
        trajectory file and the index of the frame in the trajectory file.
        """
        return self.__mapping.copy()

    @property
    def trajectories(self) -> list:
        """
        Name of the trajectories.
        """
        return self.__trajectories.copy()

    @property
    def timestep(self) -> float:
        """
        Timestep between two frames.
        """
        return self.__timestep

    @property
    def reference_system(self) -> _mdcore.systems.MDSYSTEM:
        """
        The simulated system.
        """
        return self.__reference_system.copy()

    @property
    def reverse_velocities(self) -> list:
        """
        If reverse velocities of each trajectory when reading.
        """
        return self.__reverse_velocities

    @property
    def cacheable(self) -> bool:
        """
        If the frames are cacheable.
        """
        return bool(self.__frame_cache)

    @cacheable.setter
    def cacheable(self, v: bool) -> None:
        """
        Set if the frames are cacheable.
        """
        if (v):
            if (self.__frame_cache is None):
                self.__frame_cache = [None] * self.__length
        else:
            self.__frame_cache = None
            _gc.collect()


def concatenate(paths: list):
    """
    Concatenate paths.

    Parameters
    ----------
    paths : List[somd.apps.path_sampling.utils.path.PATH]
        The paths to concatenate.
    """
    result = paths[0]
    for i in range(1, len(paths)):
        result = result + paths[i]
    return result
