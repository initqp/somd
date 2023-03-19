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

from somd import apps as _mdapps
from somd import core as _mdcore
from somd.constants import SOMDDEFAULTS as _d

__all__ = ['SIMULATION']


class SIMULATION(object):
    """
    Base class to define a simulation process.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    integrator : somd.core.integrator.INTEGRATOR
        The integrator that propagate the simulated system.
    barostat : somd.apps.barostat.BAROSTAT
        The barostat.
    trajectories : List(somd.apps.trajectories.*WRITER)
        A list of trajectory writers.
    loggers : List(somd.apps.loggers.*LOGGER)
        A list of simulation data loggers.
    """

    def __init__(self,
                 system: _mdcore.systems.MDSYSTEM,
                 integrator: _mdcore.integrators.INTEGRATOR,
                 barostat: _mdapps.barostats.BAROSTAT = None,
                 trajectories: list = [],
                 loggers: list = []) -> None:
        """
        Create a SIMULATION instance.
        """
        self.__system = system
        self.__loggers = loggers
        self.__integrator = integrator
        self.__restarted = False
        self.__read_force = False
        self.__initialized = False
        self.__post_step_objects = []
        # Check the groups.
        if (len(system.groups) == 0):
            raise RuntimeError('No group has been bound to the system!')
        # Check the potentials.
        if (len(system.potentials) == 0):
            raise RuntimeError('No potential has been bound to the system!')
        # Bind the system with the integrator.
        integrator.bind_system(system)
        # Bind the barostat to the system.
        if (barostat is not None):
            barostat.bind_integrator(integrator)
            self.__post_step_objects.append(barostat)
        # Bind the trajectory writers to the system.
        for t in trajectories:
            t.bind_integrator(integrator)
            self.__post_step_objects.append(t)
        # Bind the loggers to the system.
        for l in loggers:
            l.bind_integrator(integrator)
            self.__post_step_objects.append(l)

    def _loop(self) -> None:
        """
        The simulation loop.
        """
        self.__integrator.propagate()
        for obj in self.__post_step_objects:
            obj.update()

    def _initialize(self) -> None:
        """
        Initialize the post step objects and forces.
        """
        # Update potentials except PLUMED. Note when the forces are read
        # from a HDF5 file, we should not update the potentials.
        # This will ensure a strict restarting when using barostats.
        potential_list = []
        if (_d.POTLIST is None):
            tmp = list(range(0, len(self.system.potentials)))
        else:
            tmp = _d.POTLIST.copy()
        for i in tmp:
            p = self.system.potentials[i]
            if (p.__class__.__name__ != 'PLUMED'):
                potential_list.append(i)
        if (not self.__read_force):
            self.system.update_potentials(potential_list)
        # Initialize post step objects.
        for obj in self.__post_step_objects:
            obj.initialize()
        self.__initialized = True

    def run(self, n_steps: int) -> None:
        """
        Run the simulation.

        Parameters
        ----------
        n_steps : int
            The number of step to run.
        """
        # Initialize only once.
        if (not self.__initialized):
            self._initialize()
        # IKUZO!
        for i in range(0, n_steps):
            self._loop()

    def dump_restart(self, file_name: str) -> None:
        """
        Dump SOMD HDF5 restart file.

        Parameters
        ----------
        file_name : str
            Name of the restart file.
        """
        t = _mdapps.trajectories.H5WRITER(file_name, restart_file=True)
        t.bind_integrator(self.__integrator)
        t.initialize()
        t.write()
        del t

    def restart_from(self,
                     file_name: str,
                     frame: int = -1,
                     read_step: bool = True,
                     **kwargs) -> None:
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
        self.__restarted = True
        t = _mdapps.trajectories.H5READER(file_name, **kwargs)
        t.bind_integrator(self.__integrator)
        # Read data.
        if (frame == -1):
            frame = t.n_frames - 1
        if ('forces' in t.root.keys()):
            if ('read_forces' not in kwargs):
                self.__read_force = True
            elif ('read_forces' in kwargs and kwargs['read_forces']):
                self.__read_force = True
        t.read(frame)
        # Read time step.
        if (read_step):
            if ('steps' in t.root.keys()):
                frame = t.root['steps'][frame]
            else:
                frame = 0
            self.__integrator.step = int(frame)
            for p in self.__system.potentials:
                if (p.__class__.__name__ == 'PLUMED'):
                    p.step = frame + 1
        del t

    @property
    def current_step(self):
        """
        Current simulation step.
        """
        return self.__integrator.step

    @property
    def system(self) -> _mdcore.systems.MDSYSTEM:
        """
        The simulated system.
        """
        return self.__system

    @property
    def integrator(self) -> _mdcore.integrators.INTEGRATOR:
        """
        The invoked integrator.
        """
        return self.__integrator

    @property
    def post_step_objects(self) -> list:
        """
        A list of objects that constains a zero-parameter 'update' method,
        which will be invoked after each timestep.
        """
        return self.__post_step_objects
