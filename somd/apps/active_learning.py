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
The active learning workflow of building a NEP model.
"""

import os as _os
import numpy as _np
import shutil as _sh
import random as _rn
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import utils as _mdutils
from somd.potentials import NEP as _NEP
from . import utils as _apputils

__all__ = ['ACTIVELEARNING']


class ACTIVELEARNING(_mdapps.simulations.STAGEDSIMULATION):
    """
    The active learning workflow of building a NEP model.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    integrator : somd.core.integrator.INTEGRATOR
        The integrator that propagates the simulated system.
    potential_generators : List[Callable]
        Generator functions of potential calculators.
    reference_potentials : List[int]
        Indices of the reference potentials.
    learning_parameters : dict
        The parameters that define a learning process. Valid keys of this
        dictionary are:
        - n_potentials : int
            Default value : 4
            Number of potentials to train.
        - max_md_runs_per_iter : int
            Default value : 1
            Maximum number of MD runs in each training iteration.
        - max_md_steps_per_iter : int
            Default value : 50000
            Maximum number of MD steps in each training iteration.
        - min_new_structures_per_iter : int
            Default value : 20
            Minimum number of newly collected structures in each training
            iteration. If number of the candidate structures ls less than
            this threshold, no training will be triggered. But when the total
            number of untrained structure is larger than threshold, a training
            will be triggered, using all the untrained structures.
        - max_new_structures_per_iter : int
            Default value : 250
            Maximum number of newly collected structures in each training
            iteration. If number of the candidate structures ls larger than
            this threshold, these structures will be striped using the
            specified method.
        - perform_clustering : bool
            Default value : True
            If using the clustering method to strip the candidate structures.
        - msd_lower_limit : float
            Default value : 50.0
            Lower boundary of the force MSD when identifying new structures.
            In unit of (kJ/mol/nm).
        - msd_upper_limit : float
            Default value : 250.0
            Upper boundary of the force MSD when identifying new structures.
            In unit of (kJ/mol/nm).
        - initial_training_set : str
            Path to the initial training set (exyz file).
        - initial_testing_set : str
            Default value : initial_training_set
            Path to the initial testing set (exyz file).
        - initial_potential_files : List[str]
            Paths to the initial potential files (nep.txt files). If this
            option appears, the first iteration of training will not be
            performed.
        - energy_shift : float
            Shift the total energy by this value before recording the total
            energy to the trajectory. In unit of (kJ/mol).
    nep_parameters : str
        The keywords and corresponding values to be used in a nep.in file.
        Different keywords should be split by newlines, as in the nep.in file.
    nep_command : str
        Command to submit a NEP training job.
    use_tabulating : bool
        If invoke the tabulated version of NEP.
    post_step_objects : List[object]:
        The post step objects, including the barostat.
    output_prefix : str
        Prefix of the output file.
    """

    def __init__(self,
                 system: _mdcore.systems.MDSYSTEM,
                 integrator: _mdcore.integrators.INTEGRATOR,
                 potential_generators: list,
                 reference_potentials: list,
                 learning_parameters: dict,
                 nep_parameters: str = '',
                 nep_command: str = 'nep',
                 use_tabulating: bool = False,
                 post_step_objects: list = [],
                 output_prefix: str = '') -> None:
        """
        Create an ACTIVELEARNING instance.
        """
        self.__nep_command = nep_command
        self.__nep_parameters = nep_parameters
        self.__use_tabulating = use_tabulating
        self.__learning_parameters = learning_parameters
        self.__reference_potentials = reference_potentials
        self.__write_nep_types = _apputils.nep.check_nep_parameters(
            nep_parameters, system.atomic_symbols)
        if (output_prefix == ''):
            output_prefix = 'active_learning'
        super().__init__(system, integrator, potential_generators,
                         post_step_objects, output_prefix)
        self.__check_learning_parameters()

    def __check_learning_parameters(self) -> None:
        """
        Check the learning parameters.
        """
        param = self.__learning_parameters
        # required parameters
        if ('initial_training_set' not in param.keys() and self.n_iter == 0):
            raise KeyError('The key "initial_training_set" is required!')
        # default parameters
        if ('n_potentials' not in param.keys()):
            param['n_potentials'] = 4
        if ('perform_clustering' not in param.keys()):
            param['perform_clustering'] = False
        if ('max_md_runs_per_iter' not in param.keys()):
            param['max_md_runs_per_iter'] = 1
        if ('max_md_steps_per_iter' not in param.keys()):
            param['max_md_steps_per_iter'] = 50000
        if ('msd_lower_limit' not in param.keys()):
            param['msd_lower_limit'] = 50.0
        if ('msd_upper_limit' not in param.keys()):
            param['msd_upper_limit'] = 250.0
        if ('min_new_structures_per_iter' not in param.keys()):
            param['min_new_structures_per_iter'] = 20
        if ('max_new_structures_per_iter' not in param.keys()):
            param['max_new_structures_per_iter'] = 250
        if ('initial_testing_set' not in param.keys()):
            param['initial_testing_set'] = param['initial_training_set']
        if ('energy_shift' not in param.keys()):
            param['energy_shift'] = 0.0
        # some checks
        if ('initial_potential_files' in param.keys() and self.n_iter == 0):
            if (len(param['initial_potential_files']) !=
                    param['n_potentials']):
                message = 'Number of the initial potential files should ' + \
                          'equal to the number of required potentials!'
                raise RuntimeError(message)
        if (param['min_new_structures_per_iter'] >
                param['max_new_structures_per_iter']):
            message = 'min_new_structures_per_iter is larger than ' + \
                      'max_new_structures_per_iter !'
            raise RuntimeError(message)
        if ('initial_training_set' in param.keys() and self.n_iter != 0):
            param['initial_training_set'] = \
                self.root['/iteration_data/0'].attrs['training_sets'][0]
            message = 'The key "initial_training_set" will be ignored ' + \
                      'since the simulation is being restarted.'
            _mdutils.warning.warn(message)
        if ('initial_potential_files' in param.keys() and self.n_iter != 0):
            message = 'The key "initial_potential_files" will be ignored ' + \
                      'since the simulation is being restarted.'
            _mdutils.warning.warn(message)

    def __initialize_iteration_data_group(self, h5_path: str) -> None:
        """
        Initialize an iteration data group in the output file.

        Parameters
        ----------
        h5_path : str
            Path to the HDF5 data group.
        """
        h5_group = self.root[h5_path]
        h5_group.attrs['initialized'] = True
        h5_group.attrs['system_data'] = ''
        h5_group.attrs['invoked_nep'] = ''
        h5_group.attrs['pre_training'] = False
        h5_group.attrs['training_sets'] = []
        h5_group.attrs['visited_structures'] = ''
        h5_group.attrs['accepted_structures'] = ''
        param = self.__learning_parameters
        self._create_dataset(h5_path + '/n_visited_structures',
                             (1,), (1,), int)
        h5_group['n_visited_structures'][0] = 0
        self._create_dataset(h5_path + '/n_accurate_structures',
                             (1,), (1,), int)
        h5_group['n_accurate_structures'][0] = 0
        self._create_dataset(h5_path + '/n_candidate_structures',
                             (1,), (1,), int)
        h5_group['n_candidate_structures'][0] = 0
        self._create_dataset(h5_path + '/n_accepted_structures',
                             (1,), (1,), int)
        h5_group['n_accepted_structures'][0] = 0
        self._create_dataset(h5_path + '/n_failed_structures',
                             (1,), (1,), int)
        h5_group['n_failed_structures'][0] = 0
        self._create_dataset(h5_path + '/n_untrained_structures',
                             (1,), (1,), int)
        h5_group['n_untrained_structures'][0] = 0
        self._create_dataset(h5_path + '/energy_shift',
                             (1,), (1,), _np.double, 'kJ/mol')
        h5_group['energy_shift'][0] = param['energy_shift']
        self._create_dataset(h5_path + '/min_new_structures_per_iter',
                             (1,), (1,), int)
        h5_group['min_new_structures_per_iter'][0] = \
            param['min_new_structures_per_iter']
        self._create_dataset(h5_path + '/max_new_structures_per_iter',
                             (1,), (1,), int)
        h5_group['max_new_structures_per_iter'][0] = \
            param['max_new_structures_per_iter']
        self._create_dataset(h5_path + '/max_force_msd',
                             (1,), (1,), _np.double, 'kJ/mol/nm')
        h5_group['max_force_msd'][0] = 0.0
        self._create_dataset(h5_path + '/force_msd_lower_limit',
                             (1,), (1,), _np.double, 'kJ/mol/nm')
        h5_group['force_msd_lower_limit'][0] = param['msd_lower_limit']
        self._create_dataset(h5_path + '/force_msd_upper_limit',
                             (1,), (1,), _np.double, 'kJ/mol/nm')
        h5_group['force_msd_upper_limit'][0] = param['msd_upper_limit']
        self._create_dataset(h5_path + '/candidate_structure_indices',
                             (0,), (None,), int)
        self._create_dataset(h5_path + '/accepted_structure_indices',
                             (0,), (None,), int)
        self._create_dataset(h5_path + '/force_msd',
                             (0,), (None,), _np.double, 'kJ/mol/nm')
        self._create_dataset(h5_path + '/accepted_structure_energies',
                             (0,), (None,), _np.double, 'kJ/mol')
        progress = h5_group.create_group('progress')
        progress.attrs['training_finished'] = \
            [False for i in range(0, param['n_potentials'])]
        progress.attrs['propagation_finished'] = False
        progress.attrs['ab_initial_finished'] = False
        self.root.flush()

    def __reset_propagation_data(self, h5_path: str) -> None:
        """
        Reset the propagation data.

        Parameters
        ----------
        h5_path : str
            Path to the HDF5 data group.
        """
        h5_group = self.root[h5_path]
        h5_group['max_force_msd'][0] = 0
        h5_group['n_visited_structures'][0] = 0
        h5_group['n_accurate_structures'][0] = 0
        h5_group['n_candidate_structures'][0] = 0
        h5_group['n_accepted_structures'][0] = 0
        h5_group['n_failed_structures'][0] = 0
        h5_group['n_untrained_structures'][0] = 0
        h5_group['candidate_structure_indices'].resize((0,))
        h5_group['accepted_structure_indices'].resize((0,))
        h5_group['accepted_structure_energies'].resize((0,))
        h5_group['force_msd'].resize((0,))
        self.root.flush()

    def __update_neps(self, n_iter: int) -> tuple:
        """
        Update the trained potentials, then select the potential with minimal
        total loss to propagate the system.

        Parameters
        ----------
        n_iter : int
            Number of the training iteration.
        """
        self.__neps = []
        # Prepare new potentials.
        active_potential_index = 0
        param = self.__learning_parameters
        nep_path = '/potential_{:d}/nep.txt'
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        working_dir = h5_group.attrs['working_directory']
        training_state = h5_group['progress'].attrs['training_finished']
        n_potentials = min(param['n_potentials'], len(training_state))
        potential_files = [working_dir + nep_path.format(i)
                           for i in range(0, n_potentials)]
        for i in range(0, n_potentials):
            p = _NEP(range(0, self.system.n_atoms), potential_files[i],
                     self.system.atomic_symbols, self.__use_tabulating)
            self.__neps.append(p)
        # Determine which NEP to use.
        try:
            files = [working_dir + '/potential_{:d}/loss.out'.format(i)
                     for i in range(0, n_potentials)]
            losses = [_apputils.nep.get_loss(f)[5] for f in files]
            active_potential_index = _np.argmin(losses)
        except:
            active_potential_index = 0
        return active_potential_index, potential_files[active_potential_index]

    def __propagate(self, active_potential_index: int) -> None:
        """
        Generate new trajectories with the trained potential. During the
        propagation, force MSD of the potentials will be calculated.

        Parameters
        ----------
        active_potential_index : int
            Index of the NEP to use.
        """
        cwd = _os.getcwd()
        param = self.__learning_parameters
        h5_group = self.root['/iteration_data/' + str(self.n_iter)]
        _os.chdir(h5_group.attrs['working_directory'])
        # Set up the potentials and simulation.
        nep = [self.__neps[active_potential_index]]
        potential_indices = range(0, len(self.potential_generators))
        potential_indices = [i for i in potential_indices if
                             i not in self.__reference_potentials]
        with self._set_up_simulation(potential_indices, nep) as simulation:
            simulation.dump_restart('initial_conditions.h5')
            # Set up writers
            data_file_name = h5_group.attrs['system_data']
            data_writer = _mdapps.loggers.DEFAULTCSVLOGGER(data_file_name)
            traj_file_name = h5_group.attrs['visited_structures']
            traj_writer = _mdapps.trajectories.H5WRITER(traj_file_name)
            simulation.post_step_objects.append(data_writer)
            simulation.post_step_objects.append(traj_writer)
            # Propagate the trajectory segment.
            force_msd_limits = [param['msd_lower_limit'],
                                param['msd_upper_limit']]
            for i in range(0, param['max_md_runs_per_iter']):
                simulation.restart_from('initial_conditions.h5',
                                        read_rng_state=False,
                                        read_nhc_data=False)
                for j in range(0, param['max_md_steps_per_iter']):
                    simulation.run(1)
                    msd = _apputils.nep.get_potentials_msd(
                        self.__neps, simulation.system)
                    n = h5_group['force_msd'].shape[0]
                    h5_group['force_msd'].resize((n + 1,))
                    h5_group['force_msd'][n] = msd
                    h5_group['n_visited_structures'][0] += 1
                    if (msd < force_msd_limits[0]):
                        h5_group['n_accurate_structures'][0] += 1
                    elif (msd > force_msd_limits[1]):
                        h5_group['n_failed_structures'][0] += 1
                    else:
                        h5_group['n_candidate_structures'][0] += 1
                        index = param['max_md_steps_per_iter'] * i + j
                        h5_group['candidate_structure_indices'].resize(
                            (h5_group['n_candidate_structures'][0],))
                        h5_group['candidate_structure_indices'][-1] = index
                # Here we reset every potential calculator to reduce memory
                # effects, especially for PLUMED.
                if (i != (param['max_md_runs_per_iter'] - 1)):
                    for potential in simulation.system.potentials:
                        potential.reset()
                self.root.flush()
            h5_group['max_force_msd'][0] = max(h5_group['force_msd'])
            h5_group['progress'].attrs['propagation_finished'] = True
            self.root.flush()
        _os.chdir(cwd)

    def __strip_candidate_structures(self) -> None:
        """
        Reduce number of the candidate structures.
        """
        param = self.__learning_parameters
        h5_group = self.root['/iteration_data/' + str(self.n_iter)]
        n_candidate_structures = h5_group['n_candidate_structures'][0]
        result = list(range(0, n_candidate_structures))
        indices = h5_group['candidate_structure_indices']
        if (n_candidate_structures < param['max_new_structures_per_iter']):
            pass
        elif (not param['perform_clustering']):
            result = _rn.sample(result, param['max_new_structures_per_iter'])
        else:
            raise NotImplementedError()
        h5_group['n_accepted_structures'][0] = len(result)
        h5_group['accepted_structure_indices'].resize((len(result),))
        h5_group['accepted_structure_indices'][:] = \
            [indices[i] for i in result]
        self.root.flush()

    def __perform_ab_initio_calculations(self) -> None:
        """
        Perform ab initio calculations to the candidate structures.
        """
        cwd = _os.getcwd()
        h5_group = self.root['/iteration_data/' + str(self.n_iter)]
        _os.chdir(h5_group.attrs['working_directory'])
        # First set up the system and the integrator.
        system = self.system.copy()
        for index in self.__reference_potentials:
            system.potentials.append(self.potential_generators[index]())
        integrator = self.integrator.copy()
        integrator.bind_system(system)
        # Then set up the writer.
        traj_file_name = h5_group.attrs['accepted_structures']
        traj_writer = _mdapps.trajectories.EXYZWRITER(
            traj_file_name, write_velocities=False, wrap_positions=True,
            energy_shift=h5_group['energy_shift'][0])
        traj_writer.bind_integrator(integrator)
        # For each candidate structure, we copy its atomic positions and
        # box to the system object and perform the reference potential
        # calculations.
        traj_reader = _mdapps.trajectories.H5READER(
            h5_group.attrs['visited_structures'], read_cell=True,
            read_velocities=False, read_forces=False,
            read_nhc_data=False, read_rng_state=False)
        traj_reader.bind_integrator(integrator)
        for i, j in enumerate(h5_group['accepted_structure_indices']):
            h5_group['accepted_structure_energies'].resize((i + 1,))
            h5_group['accepted_structure_energies'][i] = 0.0
            traj_reader.read(j)
            try:
                system.update_potentials()
            except:
                # Reset SIESTA on fail.
                for potential in system.potentials:
                    if (potential.__class__.__name__ == 'SIESTA'):
                        potential.reset()
            else:
                # Only update the exyz file on success.
                traj_writer.update()
                h5_group['accepted_structure_energies'][i] = \
                    system.energy_potential - h5_group['energy_shift'][0]
        # Finally we clean the system.
        h5_group['progress'].attrs['ab_initial_finished'] = True
        for potential in system.potentials:
            potential.finalize()
        del system, integrator, traj_writer
        self.root.flush()
        _os.chdir(cwd)

    def __train(self,
                n_iter: int,
                restart_files: list = None) -> None:
        """
        Train the model.

        Parameters
        ----------
        n_iter : int
            Number of the training iteration.
        restart_files : List[str]
            Paths of the restart files.
        """
        cwd = _os.getcwd()
        param = self.__learning_parameters
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        working_dir = h5_group.attrs['working_directory']
        training_sets = h5_group.attrs['training_sets']
        training_state = h5_group['progress'].attrs['training_finished']
        # Setup the training and testing set.
        _apputils.nep.cat_exyz(training_sets, working_dir + '/train.xyz')
        _sh.copy(param['initial_testing_set'], working_dir + '/test.xyz')
        # Train.
        # Note that there we use len(training_state) instead of
        # param['n_potentials']. This is because when restarting a half-trained
        # iteration, the 'n_potentials' parameter may change.
        for i in range(0, len(training_state)):
            if (not training_state[i]):
                potential_dir = working_dir + '/potential_{:d}'.format(i)
                if (_os.path.isdir(potential_dir)):
                    _sh.rmtree(potential_dir)
                _os.mkdir(potential_dir)
                _os.chdir(potential_dir)
                _os.symlink('../train.xyz', 'train.xyz')
                _os.symlink('../test.xyz', 'test.xyz')
                if (restart_files is not None):
                    _sh.copy(restart_files[i], 'nep.restart')
                if (self.__write_nep_types):
                    _apputils.nep.make_nep_in(self.__nep_parameters,
                                              self.system.atomic_symbols)
                else:
                    _apputils.nep.make_nep_in(self.__nep_parameters)
                errno = _os.system(self.__nep_command + '> nep.log 2> nep.err')
                if (errno != 0):
                    message = 'NEP command {:s} returns a non-zero exit code!'
                    raise RuntimeError(message.format(self.__nep_command))
                training_state[i] = True
                h5_group['progress'].attrs['training_finished'] = \
                    training_state
                self.root.flush()
                _os.chdir(cwd)

    def __perform_pretraining(self) -> None:
        """
        Perform the first iteration of training.
        """
        h5_path = '/iteration_data/0'
        param = self.__learning_parameters
        working_dir = self.root[h5_path].attrs['working_directory']
        if (not self.root[h5_path].attrs['initialized']):
            self.__initialize_iteration_data_group(h5_path)
            self.root[h5_path].attrs['pre_training'] = True
            self.root[h5_path].attrs['training_sets'] = \
                [param['initial_training_set']]
            self.root.flush()
            if ('initial_potential_files' in param.keys()):
                # Read the potentials.
                potentials_data = []
                for i in range(0, param['n_potentials']):
                    with open(param['initial_potential_files'][i], 'rb') as fp:
                        potentials_data.append(fp.read())
                # Write the sets.
                _apputils.nep.cat_exyz([param['initial_training_set']],
                                       working_dir + '/train.xyz')
                _apputils.nep.cat_exyz([param['initial_testing_set']],
                                       working_dir + '/test.xyz')
                # Write the potentials.
                for i in range(0, param['n_potentials']):
                    potential_dir = working_dir + '/potential_{:d}'.format(i)
                    _os.mkdir(potential_dir)
                    with open(potential_dir + '/nep.txt', 'wb') as fp:
                        fp.write(potentials_data[i])
                self.root[h5_path]['progress'].attrs['training_finished'] = \
                    [True for i in range(0, param['n_potentials'])]
                del potentials_data
                self.root.flush()
            else:
                self.__train(0)
        else:
            self.__train(0)

    def __perform_training(self) -> None:
        """
        Perform the regular iterations of training.
        """
        param = self.__learning_parameters
        # Bind the updated potentials.
        h5_path = '/iteration_data/' + str(self.n_iter)
        h5_path_old = '/iteration_data/' + str(self.n_iter - 1)
        h5_group = self.root[h5_path]
        h5_group_old = self.root[h5_path_old]
        working_dir = h5_group.attrs['working_directory']
        working_dir_old = h5_group_old.attrs['working_directory']
        # Initialize the iteration.
        if (not h5_group.attrs['initialized']):
            self.__initialize_iteration_data_group(h5_path)
            h5_group.attrs['system_data'] = \
                working_dir + '/system_data.csv'
            h5_group.attrs['visited_structures'] = \
                working_dir + '/visited_structures.h5'
            h5_group.attrs['accepted_structures'] = \
                working_dir + '/accepted_structures.xyz'
            self.root.flush()
        # Run the simulation and harvest candidate structures.
        if (not h5_group['progress'].attrs['propagation_finished']):
            nep_index, nep_name = self.__update_neps(self.n_iter - 1)
            h5_group.attrs['invoked_nep'] = nep_name
            self.__reset_propagation_data(h5_path)
            self.__propagate(nep_index)
            self.__strip_candidate_structures()
        # Perform ab initio calculations.
        if (not h5_group['progress'].attrs['ab_initial_finished']):
            if (h5_group['n_candidate_structures'][0] > 0):
                self.__perform_ab_initio_calculations()
        # Determine if train new potentials.
        n_new_structures = h5_group['n_accepted_structures'][0] + \
            h5_group_old['n_untrained_structures'][0]
        if (n_new_structures > param['min_new_structures_per_iter']):
            # Train new potentials.
            training_sets = []
            for j in range(1, (self.n_iter + 1)):
                h5_group = self.root['/iteration_data/' + str(j)]
                if (h5_group['n_accepted_structures'][0] > 0):
                    training_sets.append(h5_group.attrs['accepted_structures'])
            training_sets.append(param['initial_training_set'])
            h5_group.attrs['training_sets'] = training_sets
            self.root.flush()
            self.__train(self.n_iter)
        else:
            # Use the previous trained potentials. Here we use relative paths
            # to perform the symlink to avoid conflictions.
            cwd = _os.getcwd()
            _os.chdir(working_dir)
            test_xyz_old = working_dir_old + '/test.xyz'
            train_xyz_old = working_dir_old + '/train.xyz'
            potential_dir_old = working_dir_old + '/potential_{:d}'
            if (not _os.path.islink('test.xyz')):
                _os.symlink(_os.path.relpath(test_xyz_old), 'test.xyz')
            if (not _os.path.islink('train.xyz')):
                _os.symlink(_os.path.relpath(train_xyz_old), 'train.xyz')
            for j in range(0, param['n_potentials']):
                if (not _os.path.islink('potential_{:d}'.format(j))):
                    _os.symlink(_os.path.relpath(potential_dir_old.format(j)),
                                'potential_{:d}'.format(j))
            h5_group.attrs['training_sets'] = \
                h5_group_old.attrs['training_sets']
            h5_group['progress'].attrs['training_finished'] = \
                [True for i in range(0, param['n_potentials'])]
            # We did not train the new accepted structures, accumulate
            # the untrained structure count.
            h5_group['n_untrained_structures'][0] = n_new_structures
            _os.chdir(cwd)

    def run(self, n_iterations: int = 1) -> None:
        """
        Run the training.

        Parameters
        ----------
        n_iterations : int
            Number of iterations to run.
        """
        # Train or copy the initial potentials.
        if (self.n_iter == 0):
            self.__perform_pretraining()
            self._increase_n_iter()
        # Train for more iterations.
        for i in range(0, n_iterations):
            self.__perform_training()
            self._increase_n_iter()

    @property
    def nep_parameters(self) -> str:
        """
        The keywords and corresponding values to be used in a nep.in file.
        """
        return self.__nep_parameters

    @nep_parameters.setter
    def nep_parameters(self, v: str) -> None:
        """
        Set the keywords and corresponding values to be used in a nep.in file.
        """
        self.__nep_parameters = v
        self.__write_nep_types = _apputils.nep.check_nep_parameters(
            self.__nep_parameters, self.system.atomic_symbols)
