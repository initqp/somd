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
import os as _os
import abc as _ab
import h5py as _h5
import copy as _cp
import typing as _tp
import contextlib as _cl
from somd import _version
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import utils as _mdutils
from . import utils as _apputils

__all__ = ['SIMULATION', 'STAGEDSIMULATION']


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
        self.__system = system
        self.__loggers = loggers
        self.__integrator = integrator
        self.__read_force = False
        self.__initialized = False
        self.__post_step_objects = []
        self.__n_dof = [g.n_dof for g in self.system.groups]
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
            self.__post_step_objects.append(barostat)
        # Bind the trajectory writers to the system.
        for t in trajectories:
            self.__post_step_objects.append(t)
        # Bind the loggers to the system.
        for l in loggers:
            self.__post_step_objects.append(l)

    def __del__(self) -> None:
        """
        Finalize the simulation.
        """
        for potential in self.__system.potentials:
            potential.finalize()
        for post_step_object in self.__post_step_objects:
            post_step_object.finalize()

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
        for i in range(0, len(self.system.potentials)):
            p = self.system.potentials[i]
            if p.__class__.__name__ != 'PLUMED':
                potential_list.append(i)
        if not self.__read_force:
            self.system.update_potentials(potential_list)
        # Initialize post step objects.
        for obj in self.__post_step_objects:
            obj.bind_integrator(self.__integrator)
            obj.initialize()
        self.__initialized = True

    def _update_internal_states(self) -> None:
        """
        Update internal states of the system and the integrator.
        """
        n_dof = [g.n_dof for g in self.system.groups]
        if n_dof != self.__n_dof:
            self.integrator.bind_system(self.system)

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
        t.bind_integrator(self.__integrator)
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
            self.__integrator.step = int(frame)
            for p in self.__system.potentials:
                if p.__class__.__name__ == 'PLUMED':
                    p.step = frame + 1
        del t

    @property
    def current_step(self) -> int:
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


class STAGEDSIMULATION(_ab.ABC):
    """
    The staged simulation protocol. In this type of simulations, the
    `somd.apps.simulations.SIMULATION` class will be dynamically instantiated
    and destroyed, alone with the required potential calculators and
    integrators.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    integrator : somd.core.integrator.INTEGRATOR
        The integrator that propagates the simulated system.
    potential_generators : List[Callable]
        Generator functions of potential calculators.
    post_step_objects : List[object]:
        The post step objects, including the barostat and any trajectory
        writers.
    output_prefix : str
        Prefix of the output file.
    """

    def __init__(
        self,
        system: _mdcore.systems.MDSYSTEM,
        integrator: _mdcore.integrators.INTEGRATOR,
        potential_generators: _tp.List[_tp.Callable],
        post_step_objects: _tp.List[_mdapps.utils.post_step.POSTSTEPOBJ] = [],
        output_prefix: str = 'staged_simulation',
    ) -> None:
        """
        Create a STAGEDSIMULATION instance.
        """
        self.__system = system
        self.__integrstor = integrator
        self.__post_step_objects = post_step_objects
        self.__potential_generators = potential_generators
        self.__output_prefix = output_prefix
        self.__check_post_step_objects()
        self.__check_system()
        self.__initialize_outputs()

    def __dump_attributes(self) -> None:
        """
        Dump attributes to the output file.
        """
        self.__root.attrs['program'] = 'SOMD'
        self.__root.attrs['version'] = str(_version.get_versions())
        self.__root.attrs['created_time'] = __import__('time').ctime()
        self.__root.attrs['title'] = self.__output_prefix + '.h5'
        self.__root.attrs['root_directory'] = _os.getcwd()
        self.__root.attrs['working_directory'] = (
            _os.getcwd() + '/' + self.__output_prefix + '.dir'
        )

    def __dump_groups(self) -> None:
        """
        Dump groups to the output file.
        """
        group = self.__root.create_group('iteration_data')
        group.create_group(str(0))
        self._set_up_iter_dir()

    def __initialize_outputs(self) -> None:
        """
        Create output files and directories.
        """
        title = self.__output_prefix + '.h5'
        directory = self.__output_prefix + '.dir'
        if _os.path.exists(title):
            self.__is_restart = True
            if not _os.path.isdir(directory):
                message = (
                    'Can not find directory {:s}! SOMD will back up '
                    + 'the file {:s} and run the required simulation '
                    + 'from scratch!'
                )
                _mdutils.warning.warn(message.format(directory, title))
                self.__is_restart = False
                _apputils.backup.back_up(title)
                _os.mkdir(directory)
            else:
                message = (
                    'Found directory {:s}! SOMD will restart the '
                    + 'required simulation from the break point!'
                )
                _mdutils.warning.warn(message.format(directory))
        else:
            if _os.path.isdir(directory):
                message = (
                    'Can not find file {:s}! SOMD will back up the '
                    + 'directory {:s} and run the required simulation '
                    + 'from scratch!'
                )
                _mdutils.warning.warn(message.format(title, directory))
                _apputils.backup.back_up(directory)
            _os.mkdir(directory)
            self.__is_restart = False
        self.__root = _h5.File(title, 'a')
        if self.__is_restart:
            group = self.__root['/iteration_data/' + str(self.n_iter)]
            working_dir = group.attrs['working_directory']
            if not _os.path.isdir(working_dir):
                if _os.path.exists(working_dir):
                    _os.remove(working_dir)
                _os.mkdir(working_dir)
        else:
            self.__dump_attributes()
            self.__dump_groups()

    def __check_post_step_objects(self) -> None:
        """
        Check the initialization state of the post step objects.
        """
        for obj in self.__post_step_objects:
            if obj.initialized:
                message = (
                    'Post step objects that are passed to the '
                    + 'STAGEDSIMULATION wrapper must be uninitialized!'
                )
                raise RuntimeError(message)

    def __check_system(self) -> None:
        """
        Check the state of the system object.
        """
        if len(self.__system.potentials) != 0:
            message = (
                'The system object that is passed to the STAGEDSIMULATION '
                + 'wrapper should not contain potential calculators!'
            )
            raise RuntimeError(message)

    def _create_dataset(
        self,
        path: str,
        shape: tuple,
        max_shape: tuple,
        data_type: str,
        unit: str = None,
    ) -> None:
        """
        Create a new data set.

        Parameters
        ----------
        path : str
            Path of the dataset.
        shape : tuple
            Initial shape of the data set.
        max_shape : tuple
            Maximum shape of the data set.
        data_type : str
            The data type.
        unit : str
            Unit of the data set.
        """
        self.__root.create_dataset(
            path, shape, maxshape=max_shape, dtype=data_type
        )
        if unit is not None:
            self.__root[path].attrs['units'] = unit

    @_cl.contextmanager
    def _set_up_simulation(
        self,
        potential_indices: _tp.List[int],
        extra_potentials: _tp.List[_mdcore.potential_base.POTENTIAL] = [],
    ) -> SIMULATION:
        """
        Set up a simulation protocol using the given data.

        Parameters
        ----------
        potential_indices : List[int]
            Indices of the potentials that drive the simulation.
        extra_potentials : List[somd.core.potential_base.POTENTIAL]
            Extra potentials that drive the simulation.
        """
        self.__check_system()
        self.__check_post_step_objects()
        system = self.__system.copy()
        integrator = self.__integrstor.copy()
        post_step_objects = _cp.deepcopy(self.__post_step_objects)
        for i in potential_indices:
            system.potentials.append(self.__potential_generators[i]())
        for p in extra_potentials:
            system.potentials.append(p)
        barostat = None
        for index, obj in enumerate(self.__post_step_objects):
            if obj.__class__ == _mdapps.barostats.BAROSTAT:
                barostat = post_step_objects.pop(index)
        simulation = SIMULATION(system, integrator, barostat)
        for obj in post_step_objects:
            simulation.post_step_objects.append(obj)
        try:
            yield simulation
        finally:
            del simulation

    def _set_up_iter_dir(self) -> str:
        """
        Initialize the directory of one simulation iteration.
        """
        iter_dir = (
            self.__root.attrs['working_directory']
            + '/iteration_'
            + str(self.n_iter)
        )
        group = self.__root['/iteration_data/' + str(self.n_iter)]
        group.attrs['working_directory'] = iter_dir
        group.attrs['initialized'] = False
        self.__root.flush()
        _apputils.backup.back_up(iter_dir)
        _os.mkdir(iter_dir)
        return iter_dir

    def _increase_n_iter(self) -> None:
        """
        Increase the number of iterations by 1.
        """
        self.__root['/iteration_data'].create_group(str(self.n_iter + 1))
        self._set_up_iter_dir()

    @property
    def system(self) -> _mdcore.systems.MDSYSTEM:
        """
        The simulated system.
        """
        return self.__system

    @property
    def integrator(self) -> _mdcore.integrators.INTEGRATOR:
        """
        The integrator.
        """
        return self.__integrstor

    @property
    def post_step_objects(self) -> list:
        """
        The post step objects.
        """
        return self.__post_step_objects

    @property
    def potential_generators(self) -> list:
        """
        The potential generators.
        """
        return self.__potential_generators

    @property
    def n_iter(self) -> int:
        """
        Number of the total simulation iteration.
        """
        return len(self.__root['/iteration_data']) - 1

    @property
    def is_restart(self) -> bool:
        """
        If this simulation is a restarted one.
        """
        return self.__is_restart

    @property
    def root(self) -> _h5.File:
        """
        The HDF5 file root.
        """
        return self.__root
