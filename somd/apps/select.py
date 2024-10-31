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
Select structures for the active learning workflow.
"""

import json as _js
import numpy as _np
import typing as _tp
import random as _rn
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['STRUCTURESELECTOR']


class STRUCTURESELECTOR(object):
    """
    Select structures according to their force MSE.

    Parameters
    ----------
    file_names: List[str]
        Name of the input trajectories.
    """

    def __init__(self, file_names: _tp.List[str]) -> None:
        """
        Create an STRUCTURESELECTOR instance.
        """
        self.__objects = {}
        self.__file_names = file_names
        readers = [
            _mdapps.trajectories.H5READER(f) for f in file_names
        ]
        has_coordinates = ['coordinates' in r.root.keys() for r in readers]
        if False in has_coordinates:
            file_name = file_names[has_coordinates.index(False)]
            message = 'Can not find positions in the trajectory file "{:s}"!'
            raise KeyError(message.format(file_name))
        has_forces = ['forces' in r.root.keys() for r in readers]
        if False in has_forces:
            file_name = file_names[has_forces.index(False)]
            message = 'Can not find forces in the trajectory file "{:s}"!'
            raise KeyError(message.format(file_name))
        has_virial = ['virial' in r.root.keys() for r in readers]
        if False in has_virial:
            file_name = file_names[has_virial.index(False)]
            message = 'Can not find virial in the trajectory file "{:s}"!'
            _mdutils.warning.warn(message.format(file_name))
            self.__read_virial = False
        else:
            self.__read_virial = True
        has_same_shapes = [
            readers[0].root['coordinates'].shape == r.root['coordinates'].shape
            for r in readers
        ]
        if False in has_same_shapes:
            file_name_1 = file_names[0]
            file_name_2 = file_names[has_same_shapes.index(False)]
            message = (
                'Mismatch between shapes of coordinate array in '
                + 'trajectory "{:s}" and "{:s}"!'
            )
            raise RuntimeError(message.format(file_name_1, file_name_2))
        self.__n_frames = readers[0].root['coordinates'].shape[0]
        self.__n_atoms = readers[0].root['coordinates'].shape[1]
        del readers

        self.__objects['readers'] = [
            _mdapps.trajectories.H5READER(
                f,
                read_coordinates=True,
                read_velocities=False,
                read_forces=True,
                read_cell=True,
                read_virial=self.__read_virial,
                read_nhc_data=False,
                read_rng_state=False,
            ) for f in file_names
        ]

    def select(
        self,
        n_max: int,
        msd_f_min: float,
        msd_f_max: float,
        shuffle: bool = True,
        log_file: str = None,
    ) -> _np.ndarray:
        """
        Select structures from the input trajectories.

        Parameters
        ----------
        msd_f_min : float
            Lower boundary of the force MSD when identifying new structures.
            In unit of (kJ/mol/nm).
        msd_f_max : float
            Upper boundary of the force MSD when identifying new structures.
            In unit of (kJ/mol/nm).
        n_max : int
            Max number of structure to select.
        shuffle : bool
            If perform randomly selection.
        log_file : str
            Name of the log file.
        """
        force_msd = _np.array([
            self.get_force_msd(i) for i in range(self.__n_frames)
        ])
        mask = (force_msd > msd_f_min) & (force_msd < msd_f_max)
        indices = _np.arange(self.__n_frames)[mask]
        if mask.sum() > n_max:
            if shuffle:
                results = _np.array(
                    _rn.sample(indices.tolist(), n_max), dtype=int
                )
            else:
                results = indices[:n_max]
        else:
            results = indices

        if log_file is not None:
            info = {}
            info['n_max'] = n_max
            info['n_frames'] = int(self.__n_frames)
            info['n_selected'] = len(results)
            info['n_candidates'] = int(mask.sum())
            info['n_potentials'] = len(self.readers)
            info['msd_f_min'] = msd_f_min
            info['msd_f_max'] = msd_f_max
            info['force_msd'] = force_msd.tolist()
            info['selected_indices'] = results.tolist()
            info['candidate_indices'] = indices.tolist()
            with open(log_file, 'w') as fp:
                _js.dump(info, fp, indent=4)

        return results

    def write(
        self,
        file_name: str,
        frame_indices: _tp.List[int],
        trajectory_index: int = 0,
    ) -> None:
        """
        Write the selected structures to a trajectory file.

        Parameters
        ----------
        file_name : str
            Name of the trajectory.
        frame_indices : int
            Indices of frame to write.
        trajectory_index : int
            Index of trajectory to read data from.
        """
        integrator = _mdcore.integrators.vv_integrator(1E-3)
        system = _mdcore.systems.MDSYSTEM(self.__n_atoms)
        writer = _mdapps.trajectories.H5WRITER(
            file_name,
            write_velocities=False,
            write_forces=False,
            write_virial=False,
        )
        integrator.bind_system(system)
        writer.bind_integrator(integrator)

        reader = self.readers[trajectory_index]
        for i in frame_indices:
            system.snapshot = reader._read_snapshot(i)
            integrator.step += 1
            writer.write()

        del writer

    def get_force_msd(self, frame_index: int) -> float:
        """
        Return the maximum standard deviation of the forces calculated by
        a list of potentials.

        Parameters
        ----------
        frame_index : int
            The frame number.
        """
        sd = _np.zeros(self.__n_atoms)
        mean = _np.zeros((self.__n_atoms, 3))
        forces = []
        for r in self.readers:
            s = r._read_snapshot(frame_index)
            mean += s.forces
            forces.append(s.forces)
        mean /= len(self.readers)
        for i in range(0, self.__n_atoms):
            for j in range(0, len(self.readers)):
                sd[i] += _np.linalg.norm(forces[j][i] - mean[i]) ** 2
            sd[i] /= len(self.readers)
        return _np.max(_np.sqrt(sd))

    @property
    def readers(self) -> _tp.List:
        """
        The trajectory readers.
        """
        return self.__objects['readers']
