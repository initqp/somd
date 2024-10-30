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
The pre-defined simulation protocols.
"""

import typing as _tp
import datetime as _dt
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['SIMULATION']


class SIMULATION(object):
    """
    Base class to define a simulation process.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    integrator : somd.core.integrator.INTEGRATOR
        The integrator that propagates the simulated system.
    barostat : somd.apps.barostat.BAROSTAT
        The barostat.
    trajectories : List[somd.apps.trajectories.*WRITER]
        A list of trajectory writers.
    loggers : List[somd.apps.loggers.*LOGGER]
        A list of simulation data loggers.
    """

    def __init__(
        self,
        system: _mdcore.systems.MDSYSTEM,
        integrator: _mdcore.integrators.INTEGRATOR,
        barostat: _mdapps.barostats.BAROSTAT = None,
        trajectories: _tp.List[_mdapps.utils.post_step.POSTSTEPOBJ] = [],
        loggers: _tp.List[_mdapps.utils.post_step.POSTSTEPOBJ] = [],
    ) -> None:
        """
        Create a SIMULATION instance.
        """
        self.__objects = {}
        self.__read_force = False
        self.__initialized = False
        self.__objects['system'] = system
        self.__objects['loggers'] = loggers
        self.__objects['integrator'] = integrator
        self.__objects['post_step'] = []
        self.__objects['n_dof'] = [g.n_dof for g in self.system.groups]
        # Check the groups.
        if len(system.groups) == 0:
            raise RuntimeError('No group has been bound to the system!')
        # Check the potentials.
        if len(system.potentials) == 0:
            raise RuntimeError('No potential has been bound to the system!')
        # Bind the system with the integrator.
        integrator.bind_system(system)
        # Bind the barostat to the system.
        if barostat is not None:
            self.__objects['post_step'].append(barostat)
        # Bind the trajectory writers to the system.
        for t in trajectories:
            self.__objects['post_step'].append(t)
        # Bind the loggers to the system.
        for l in loggers:
            self.__objects['post_step'].append(l)

    def __del__(self) -> None:
        """
        Finalize the simulation.
        """
        for potential in self.__objects['system'].potentials:
            potential.finalize()
        for post_step_object in self.__objects['post_step']:
            post_step_object.finalize()

    def _loop(self) -> None:
        """
        The simulation loop.
        """
        self.__objects['integrator'].propagate()

        for obj in self.__objects['post_step']:
            obj.update()

        for i in self.__objects['plumed_indices']:
            if self.system.potentials[i].stop_flag:
                message = 'PLUMED stop criterion meet. Exiting now ...'
                _mdutils.warning.warn(message)
                exit()

    def _initialize(self) -> None:
        """
        Initialize the post step objects and forces.
        """
        # Update potentials except PLUMED. Note when the forces are read
        # from a HDF5 file, we should not update the potentials.
        # This will ensure a strict restarting when using barostats.
        potential_list = []
        for i in range(0, len(self.system.potentials)):
            p = self.system.potentials[i]
            if p.__class__.__name__ != 'PLUMED':
                potential_list.append(i)
        if not self.__read_force:
            self.system.update_potentials(potential_list)
        # Initialize post step objects.
        for obj in self.__objects['post_step']:
            obj.bind_integrator(self.__objects['integrator'])
            obj.initialize()
        self.__initialized = True

    def _update_internal_states(self) -> None:
        """
        Update internal states of the system and the integrator, and deal with
        PLUMED instances.
        """
        n_dof = [g.n_dof for g in self.system.groups]
        if n_dof != self.__objects['n_dof']:
            self.integrator.bind_system(self.system)
            self.__objects['n_dof'] = n_dof

        self.__objects['plumed_indices'] = []
        for i in range(0, len(self.system.potentials)):
            p = self.system.potentials[i]
            if p.__class__.__name__ == 'PLUMED':
                self.__objects['plumed_indices'].append(i)

    def summary(self) -> str:
        """
        Show information about the simulation.
        """
        summary_o = 'POSTSTEPOBJECTS\n'
        for o in self.post_step_objects:
            if hasattr(o, 'summary'):
                summary_o += (
                    '┣━ ' + o.summary().replace('\n', '\n┃  ').strip() + '\n'
                )
            else:
                summary_o += '┃  UNKOWN POSTSTEPOBJ\n'
        summary_o += '┗━ END'

        result = (
            self.__objects['system'].summary() + '\n'
            + self.__objects['integrator'].summary() + '\n'
            + summary_o + '\n'
        )
        return result

    def run(self, n_steps: int) -> None:
        """
        Run the simulation.

        Parameters
        ----------
        n_steps : int
            The number of step to run.
        """
        # Initialize only once.
        if not self.__initialized:
            self._initialize()
        self._update_internal_states()
        # IKUZO!
        for i in range(0, n_steps):
            self._loop()

    def run_for_clock_time(self, time: int, restart_file: str = None):
        """
        Run the simulation by integrating time steps until a fixed amount of
        clock time has elapsed. This is useful when you have a limited amount
        of computer time available, and want to run the longest simulation
        possible in that time. This method will continue taking time steps
        until the specified clock time has elapsed, then return. It also can
        automatically write out a restart file before returning, so you can
        later resume the simulation. This function was stolen from OpenMM.

        Parameters
        ----------
        time : int
            the amount of time to run for. In unit of seconds.
        restart_file : str
            if specified, a checkpoint file will be written at the end of the
            simulation.
        """
        # Initialize only once.
        if not self.__initialized:
            self._initialize()
        self._update_internal_states()
        end_time = _dt.datetime.now() + _dt.timedelta(seconds=time)
        while (_dt.datetime.now() < end_time):
            self._loop()
        if restart_file is not None:
            self.dump_restart(restart_file)

    def dump_restart(self, file_name: str) -> None:
        """
        Dump SOMD HDF5 restart file.

        Parameters
        ----------
        file_name : str
            Name of the restart file.
        """
        t = _mdapps.trajectories.H5WRITER(file_name, restart_file=True)
        t.bind_integrator(self.__objects['integrator'])
        t.write()
        del t

    def restart_from(
        self, file_name: str, frame: int = -1, read_step: bool = True, **kwargs
    ) -> None:
        """
        Restart the simulation from a SOMD HDF5 restart file.

        Parameters
        ----------
        file_name : str
            Name of the restart file.
        frame : int
            Frame of the restart file to read.
        read_step : bool
            If inherit the simulation step from the restart file.
        """
        t = _mdapps.trajectories.H5READER(file_name, **kwargs)
        t.bind_integrator(self.__objects['integrator'])
        # Read data.
        if frame == -1:
            frame = t.n_frames - 1
        if 'forces' in t.root.keys():
            if 'read_forces' not in kwargs:
                self.__read_force = True
            elif 'read_forces' in kwargs and kwargs['read_forces']:
                self.__read_force = True
        t.read(frame)
        # Read time step.
        if read_step:
            if 'steps' in t.root.keys():
                frame = t.root['steps'][frame]
            else:
                frame = 0
            self.__objects['integrator'].step = int(frame)
            for p in self.__objects['system'].potentials:
                if p.__class__.__name__ == 'PLUMED':
                    p.step = frame + 1
        del t

    @property
    def current_step(self) -> int:
        """
        Current simulation step.
        """
        return self.__objects['integrator'].step

    @property
    def system(self) -> _mdcore.systems.MDSYSTEM:
        """
        The simulated system.
        """
        return self.__objects['system']

    @property
    def integrator(self) -> _mdcore.integrators.INTEGRATOR:
        """
        The invoked integrator.
        """
        return self.__objects['integrator']

    @property
    def post_step_objects(self) -> list:
        """
        A list of objects that constains a zero-parameter 'update' method,
        which will be invoked after each timestep.
        """
        return self.__objects['post_step']
