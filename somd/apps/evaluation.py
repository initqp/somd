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
Re-run a simulation trajectory.
"""

import typing as _tp
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['EVALUATION']


class EVALUATION(object):
    """
    Re-run a simulation trajectory.

    Parameters
    ----------
    file_name : str
        Name of the trajectory.
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    trajectories : List[somd.apps.trajectories.*WRITER]
        A list of trajectory writers.
    loggers : List[somd.apps.loggers.*LOGGER]
        A list of simulation data loggers.
    interval : int
        Interval of evaluating the trajectory.
    """

    def __init__(
        self,
        file_name: str,
        system: _mdcore.systems.MDSYSTEM,
        trajectories: _tp.List[_mdapps.utils.post_step.POSTSTEPOBJ] = [],
        loggers: _tp.List[_mdapps.utils.post_step.POSTSTEPOBJ] = [],
        interval: int = 1,
    ) -> None:
        """
        Create an EVALUATION instance.
        """
        self.__objects = {}
        self.__objects['loggers'] = loggers
        self.__objects['trajectories'] = trajectories
        self.__objects['system'] = system
        self.__interval = interval

        # do checks
        reader = _mdapps.trajectories.H5READER(file_name)
        if 'coordinates' not in reader.root.keys():
            message = 'Can not find positions in the trajectory file {:s}!'
            raise KeyError(message.format(file_name))
        if 'velocities' not in reader.root.keys():
            message = 'Can not find velocities in the trajectory file {:s}!'
            _mdutils.warning.warn(message.format(file_name))
            read_velocities = False
        else:
            read_velocities = True
        if system.n_atoms != reader.root['coordinates'].shape[1]:
            message = (
                'Mismatch between number of atoms in trajectory {:s}'
                + 'and the given MD system!'
            )
            raise RuntimeError(message.format(file_name))

        self.__objects['reader'] = _mdapps.trajectories.H5READER(
            file_name,
            read_coordinates=True,
            read_velocities=read_velocities,
            read_cell=True,
            read_nhc_data=False,
            read_rng_state=False,
        )
        self.__objects['integrator'] = _mdcore.integrators.vv_integrator(1E-3)
        self.__objects['integrator'].bind_system(self.__objects['system'])
        self.__objects['reader'].bind_integrator(self.__objects['integrator'])
        for t in self.__objects['trajectories']:
            t.bind_integrator(self.__objects['integrator'])

    def run(self) -> None:
        """
        Run the evaluation.
        """
        for i in range(self.__objects['reader'].n_frames):
            self.__objects['reader'].read(i)
            if (i % self.__interval) == 0:
                self.__objects['system'].update_potentials()
                for t in self.__objects['trajectories']:
                    t.write()
                for l in self.__objects['loggers']:
                    l.write()
