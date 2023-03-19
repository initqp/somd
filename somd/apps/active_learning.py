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
import re as _re
import json as _js
import numpy as _np
import shutil as _sh
import random as _rn
from . import _backup
from somd import apps as _mdapps
from somd import core as _mdcore
from somd.potentials import NEP as _NEP
from somd.constants import SOMDDEFAULTS as _d

__all__ = ['ACTIVELEARNING']


class ACTIVELEARNING(object):
    """
    The active learning workflow of building a NEP model.

    Parameters
    ----------
    simulation : somd.apps.simulations.SIMULATION
        The simulation protocol.
    reference_potentials : List(int)
        Indices of the reference potentials.
    parameters : dict
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
        - initial_potential_files : List(str)
            Paths to the initial potential files (nep.txt files). If this
            option appears, the first iteration of training will not be
            performed.
    nep_parameters : str
        The keywords and corresponding values to be used in a nep.in file.
        Different keywords should be split by newlines, as in the nep.in file.
    nep_command : str
        Command to submit a NEP training job.
    """

    def __init__(self,
                 simulation: _mdapps.simulations.SIMULATION,
                 reference_potentials: list,
                 learning_parameters: dict,
                 nep_parameters: str = '',
                 nep_command: str = 'nep') -> None:
        """
        Create an ACTIVELEARNING instance.
        """
        self.__n_iter = 0
        self.__simulation = simulation
        self.__training_iter_data = []
        self.__nep_command = nep_command
        self.__nep_parameters = nep_parameters
        self.__learning_parameters = learning_parameters
        self.__reference_potentials = reference_potentials
        self.__n_untrained_structures = 0
        self.__check_learning_parameters()
        self.__check_nep_parameters()
        simulation.dump_restart('active_learning_start_point.h5')
        self.__initial_conditions = \
            _os.getcwd() + '/active_learning_start_point.h5'

    def __check_learning_parameters(self) -> None:
        """
        Check the learning parameters.
        """
        param = self.__learning_parameters
        # required parameters
        if ('initial_training_set' not in param.keys()):
            raise KeyError('Key \'initial_training_set\' is required!')
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
        # some checks
        if ('initial_potential_files' in param.keys()):
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

    def __check_nep_parameters(self) -> None:
        """
        Check the NEP training parameters.
        """
        self.__write_nep_types = True
        self.__nep_parameters = self.__nep_parameters.replace('\\n', '\n')
        for line in _re.split('\n', self.__nep_parameters):
            l = [i for i in line.strip().split(' ') if i != '']
            if (l != [] and l[0].lower() == 'type'):
                e_nep = l[2:]
                if (len(e_nep) != int(l[1])):
                    message = 'Wrong Number of elements in NEP parameters!'
                    raise RuntimeError(message)
                e_all = self.__simulation.system.atomic_symbols
                e_lack = [e for e in e_all if e not in e_nep]
                e_unknown = [e for e in e_nep if e not in e_all]
                if (len(e_lack) != 0):
                    message = 'Lack element {} in NEP parameters!'
                    raise RuntimeError(message.format(e_lack))
                if (len(e_unknown) != 0):
                    message = 'Unknown element {} in NEP parameters!'
                    raise RuntimeError(message.format(e_unknown))
                self.__write_nep_types = False

    def __make_nep_in(self) -> None:
        """
        Write the nep.in file.
        """
        fp = open('nep.in', 'w')
        if (self.__write_nep_types):
            elements = list(set(self.__simulation.system.atomic_symbols))
            elements.sort()
            print('type {:d}'.format(len(elements)), end='', file=fp)
            for e in elements:
                print(' ' + e, end='', file=fp)
            print('\n', file=fp)
        print(self.__nep_parameters, file=fp)
        fp.close()

    def __initialize_training_iter_dict(self) -> dict:
        """
        Initialize a dictionary that records the information about one training
        iteration.
        """
        result = dict()
        result["pre_training"] = False
        result["directory"] = ""
        result["system_data"] = ""
        result['visited_structures'] = ""
        result['accepted_structures'] = ""
        result['n_visited_structures'] = 0
        result['n_accurate_structures'] = 0
        result['n_candidate_structures'] = 0
        result['n_accepted_structures'] = 0
        result['n_failed_structures'] = 0
        result['max_forces_msd'] = 0
        result['candidate_structure_indices'] = []
        result['accepted_structure_indices'] = []
        result['forces_msd'] = []
        return result

    def __initialize_training_iter_dir(self) -> str:
        """
        Initialize the directory of one training iteration.
        """
        iter_dir = 'training_iter_{:d}'.format(self.__n_iter)
        _backup.backup(iter_dir)
        _os.mkdir(iter_dir)
        return _os.getcwd() + '/' + iter_dir

    def __bind_potentials(self,
                          potential_files: list,
                          active_potential_index: int) -> None:
        """
        Load trained models and use one of them to propagate the system.

        Parameters
        ----------
        potential_files : List(str)
            Paths of the potential files to bind.
        active_potential_index : int
            Index of the potential used to propagate the system.
        """
        self.__potentials = []
        # Prepare new potentials.
        for i in range(0, self.__learning_parameters['n_potentials']):
            n_atoms = self.__simulation.system.n_atoms
            symbols = self.__simulation.system.atomic_symbols
            p = _NEP(range(0, n_atoms), potential_files[i], symbols)
            setattr(p, 'SOMDAUTOPOT', 'NEP')
            self.__potentials.append(p)
        # Remove old potentials that binds with the system.
        old_nep_list = []
        for i, p in enumerate(self.__simulation.system.potentials):
            if (hasattr(p, 'SOMDAUTOPOT') and p.SOMDAUTOPOT == 'NEP'):
                old_nep_list.append(i)
        for i in old_nep_list:
            self.__simulation.system.potentials.pop(i)
        # Bind one of the new potentials.
        self.__simulation.system.potentials.append(
            self.__potentials[active_potential_index])
        # Make the active potential list. This is achieved by remove indices of
        # the reference potentials from the full potential list.
        potential_list = range(0, len(self.__simulation.system.potentials))
        potential_list = list(potential_list)
        for i in self.__reference_potentials:
            potential_list.remove(i)
        _d.POTLIST = potential_list

    def __update_potentials(self, work_dir: str) -> None:
        """
        Update the trained potentials, then select the potential with minimal
        total loss to propagate the system.

        Parameters
        ----------
        work_dir : str
            The The working directory that contains multiple potential files.
        """
        min_total_loss = 1E100
        active_potential_index = 0
        param = self.__learning_parameters
        potential_files = [work_dir + '/potential_{:d}/nep.txt'.format(i)
                           for i in range(0, param['n_potentials'])]
        loss_files = [work_dir + '/potential_{:d}/loss.out'.format(i)
                      for i in range(0, param['n_potentials'])]
        try:
            for i, f in enumerate(loss_files):
                loss = self._get_loss(f)[1]
                if (loss < min_total_loss):
                    active_potential_index = i
                    min_total_loss = loss
        except:
            active_potential_index = 0
        self.__bind_potentials(potential_files, active_potential_index)

    def __propagate(self) -> list:
        """
        Generate new trajectories with the trained potential. During the
        propagation, force MSD of the potentials will be calculated.
        """
        info = self.__training_iter_data[-1]
        param = self.__learning_parameters
        # Setup writers
        data_file_name = info['system_data']
        data_writer = _mdapps.loggers.DEFAULTCSVLOGGER(data_file_name)
        data_writer.bind_integrator(self.__simulation.integrator)
        data_writer.initialize()
        traj_file_name = info['visited_structures']
        traj_writer = _mdapps.trajectories.H5WRITER(traj_file_name)
        traj_writer.bind_integrator(self.__simulation.integrator)
        traj_writer.initialize()
        # Propagate the trajectory segment.
        candidate_structures = []
        info = self.__training_iter_data[-1]
        forces_msd_limits = [param['msd_lower_limit'],
                             param['msd_upper_limit']]
        for i in range(0, param['max_md_runs_per_iter']):
            self.restore_initial_conditions()
            for j in range(0, param['max_md_steps_per_iter']):
                self.__simulation.run(1)
                data_writer.update()
                traj_writer.update()
                msd = self._get_potentials_msd(self.__potentials,
                                               self.__simulation.system)
                info['forces_msd'].append(msd)
                info['n_visited_structures'] += 1
                if (msd < forces_msd_limits[0]):
                    info['n_accurate_structures'] += 1
                elif (msd > forces_msd_limits[1]):
                    info['n_failed_structures'] += 1
                else:
                    info['n_candidate_structures'] += 1
                    info['candidate_structure_indices'].append(j)
                    structure = [self.__simulation.system.positions.copy(),
                                 self.__simulation.system.box.copy()]
                    candidate_structures.append(structure)
        info['max_forces_msd'] = max(info['forces_msd'])
        del data_writer
        del traj_writer
        return candidate_structures

    def __strip_candidate_structures(self, candidate_structures: list) -> list:
        """
        Reduce number of the candidate structures.

        Parameters
        ----------
        candidate_structures : List(List(numpy.ndarray))
            Positions and boxes data of the candidate structures.
        """
        info = self.__training_iter_data[-1]
        param = self.__learning_parameters
        result = list(range(0, len(candidate_structures)))
        indices = info['candidate_structure_indices']
        if (len(candidate_structures) < param['max_new_structures_per_iter']):
            pass
        elif (not param['perform_clustering']):
            result = _rn.sample(result, param['max_new_structures_per_iter'])
        else:
            raise NotImplementedError()
        info['n_accepted_structures'] = len(result)
        info['accepted_structure_indices'] = [indices[i] for i in result]
        return [candidate_structures[i] for i in result]

    def __perform_ab_initio_calculations(self, structures: list) -> None:
        """
        Perform ab initio calculations to the candidate structures.

        Parameters
        ----------
        work_dir : str
            The working directory.
        structures : List(List(numpy.ndarray))
            Positions and boxes data of the candidate structures.
        """
        # Good God please forgive me for what I'm about to do ...
        # I fucking hate this shit and myself ...
        info = self.__training_iter_data[-1]
        traj_file_name = info['accepted_structures']
        traj_writer = _mdapps.trajectories.EXYZWRITER(
            traj_file_name, write_velocities=False, wrap_positions=True,
            potential_list=self.__reference_potentials)
        traj_writer.bind_integrator(self.__simulation.integrator)
        traj_writer.initialize()
        # First we get the current state of the simulated system.
        state = [self.__simulation.system.positions.copy(),
                 self.__simulation.system.box.copy()]
        # Then for each candidate structure, we copy its atomic positions and
        # box to the system object and perform the reference potential
        # calculations.
        for structure in structures:
            self.__simulation.system.positions[:] = structure[0][:]
            self.__simulation.system.box[:] = structure[1][:]
            try:
                self.__simulation.system.update_potentials(
                    self.__reference_potentials)
            except:
                # Reset SIESTA on fail.
                for i in self.__reference_potentials:
                    p = self.__simulation.system.potentials[i]
                    if (type(p).__name__ == 'SIESTA'):
                        p.reset()
            else:
                # Only update the exyz file on success.
                traj_writer.update()
        # Finally we restore the state of the system.
        self.__simulation.system.positions[:] = state[0][:]
        self.__simulation.system.box[:] = state[1][:]
        del traj_writer
        del state

    def __train(self,
                work_dir: str,
                training_set_list: list,
                restart_files: list = None) -> None:
        """
        Train the model.

        Parameters
        ----------
        work_dir : str
            The working directory.
        training_set_list : List(str)
            Paths of the training set exyz files.
        restart_files : List(str)
            Paths of the restart files.
        """
        param = self.__learning_parameters
        # Setup the training and testing set.
        # TODO: rational building of the testing set.
        self._cat_exyz(training_set_list, work_dir + '/train.xyz')
        _sh.copy(param['initial_testing_set'], work_dir + '/test.xyz')
        # Train.
        cwd = _os.getcwd()
        for i in range(0, param['n_potentials']):
            potential_dir = work_dir + '/potential_{:d}'.format(i)
            _os.mkdir(potential_dir)
            _os.chdir(potential_dir)
            _os.symlink('../train.xyz', 'train.xyz')
            _os.symlink('../test.xyz', 'test.xyz')
            if (restart_files is not None):
                _sh.copy(restart_files[i], 'nep.restart')
            self.__make_nep_in()
            _os.system(self.__nep_command + '> nep.log 2> nep.err')
            _os.chdir(cwd)

    @staticmethod
    def _cat_exyz(set_in: list, set_out: str) -> None:
        """
        Combine two EXYZ training sets.

        Parameters
        ----------
        set_in : List(str)
            Names of input training sets.
        set_out : str
            Name of output training set.
        """
        fp = open(set_out, 'w')
        for fn in set_in:
            if (fn is not None):
                fp_in = open(fn, 'r')
                for l in fp_in:
                    fp.write(l)
                fp_in.close()
        fp.close()

    @staticmethod
    def _get_potentials_msd(potentials: list,
                            system: _mdcore.systems.MDSYSTEM) -> float:
        """
        Return the maximum standard deviation of the forces calculated by
        a list of potentials.

        Parameters
        ----------
        potentials: List(somd.core.potential_base.POTENTIAL)
            The potentials.
        system: somd.core.systems.MDSYSTEM
            The simulated system.
        """
        sd = _np.zeros(system.n_atoms)
        mean = _np.zeros((system.n_atoms, 3))
        for p in potentials:
            p.update(system)
            mean += p.forces
        mean /= len(potentials)
        for i in range(0, system.n_atoms):
            for j in range(0, len(potentials)):
                tmp = _np.linalg.norm(potentials[j].forces[i] - mean[i])
                sd[i] += tmp ** 2
            sd[i] /= len(potentials)
        return _np.max(_np.sqrt(sd))

    @staticmethod
    def _get_loss(file_name: str) -> list:
        """
        Read losses from a loss.out file.
        """
        fp = open(file_name, 'rb')
        try:
            fp.seek(-2, _os.SEEK_END)
            while fp.read(1) != b'\n':
                fp.seek(-2, _os.SEEK_CUR)
        except OSError:
            fp.seek(0)
        loss = fp.readline().decode().strip().split(' ')
        fp.close()
        return [float(i) for i in loss if i != '']

    def run(self, n_iterations: int = 1):
        """
        Run the training.

        Parameters
        ----------
        n_iterations : int
            Number of iterations to run.
        """
        param = self.__learning_parameters
        # Train or copy the initial potentials.
        if (self.__n_iter == 0):
            if ('initial_potential_files' in param.keys()):
                # Read the potentials.
                potentials_data = []
                for i in range(0, param['n_potentials']):
                    with open(param['initial_potential_files'][i], 'rb') as fp:
                        potentials_data.append(fp.read())
            # Make the new iter_0 directory.
            work_dir = self.__initialize_training_iter_dir()
            if ('initial_potential_files' in param.keys()):
                # Write the sets.
                self._cat_exyz([param['initial_training_set']],
                               work_dir + '/train.xyz')
                self._cat_exyz([param['initial_testing_set']],
                               work_dir + '/test.xyz')
                # Write the potentials.
                for i in range(0, param['n_potentials']):
                    potential_dir = work_dir + '/potential_{:d}'.format(i)
                    _os.mkdir(potential_dir)
                    with open(potential_dir + '/nep.txt', 'wb') as fp:
                        fp.write(potentials_data[i])
                del potentials_data
            else:
                self.__train(work_dir, [param['initial_training_set']])
            info = self.__initialize_training_iter_dict()
            info["pre_training"] = True
            info['directory'] = work_dir
            self.__training_iter_data.append(info)
            with open(work_dir + '/training_info.json', 'w') as fp:
                _js.dump(info, fp, indent=4)
            self.__n_iter += 1
        # Training iterations.
        for i in range(0, n_iterations):
            # Bind the updated potentials.
            old_work_dir = self.__training_iter_data[-1]['directory']
            self.__update_potentials(old_work_dir)
            # Initialize the iteration.
            info = self.__initialize_training_iter_dict()
            work_dir = self.__initialize_training_iter_dir()
            info['directory'] = work_dir
            info['system_data'] = work_dir + '/system_data.csv'
            info['visited_structures'] = work_dir + '/visited_structures.h5'
            info['accepted_structures'] = work_dir + '/accepted_structures.xyz'
            self.__training_iter_data.append(info)
            # Run the simulation and harvest candidate structures.
            candidate_structures = self.__propagate()
            # Strip the candidate structures.
            accepted_structures = \
                self.__strip_candidate_structures(candidate_structures)
            # Dump the training results.
            with open(work_dir + '/training_info.json', 'w') as fp:
                _js.dump(info, fp, indent=4)
            # Perform ab initio calculations.
            if (len(accepted_structures) > 0):
                self.__perform_ab_initio_calculations(accepted_structures)
            # Determine if train new potentials.
            if ((len(accepted_structures) + self.__n_untrained_structures) >
                    param['min_new_structures_per_iter']):
                # Train new potentials.
                training_sets = [d['accepted_structures']
                                 for d in self.__training_iter_data[1:]
                                 if (d['n_accurate_structures'] > 0)]
                training_sets.append(param['initial_training_set'])
                self.__train(work_dir, training_sets)
                # A training process is triggered, reset the untrained
                # structure count.
                self.__n_untrained_structures = 0
            else:
                # Use the previous trained potentials.
                old_training_xyz = \
                    self.__training_iter_data[-2]['directory'] + '/train.xyz'
                new_training_xyz = \
                    self.__training_iter_data[-1]['directory'] + '/train.xyz'
                _os.symlink(old_training_xyz, new_training_xyz)
                old_test_xyz = self.__training_iter_data[-2]['directory'] + \
                    '/test.xyz'
                new_test_xyz = self.__training_iter_data[-1]['directory'] + \
                    '/test.xyz'
                _os.symlink(old_test_xyz, new_test_xyz)
                for j in range(0, param['n_potentials']):
                    old_potential_dir = \
                        self.__training_iter_data[-2]['directory'] + \
                        '/potential_{:d}'.format(j)
                    new_potential_dir = \
                        self.__training_iter_data[-1]['directory'] + \
                        '/potential_{:d}'.format(j)
                    _os.symlink(old_potential_dir, new_potential_dir)
                # We did not train the new accepted structures, accumulate
                # the untrained structure count.
                self.__n_untrained_structures += len(accepted_structures)
            # Clean up.
            del candidate_structures
            del accepted_structures
            self.__n_iter += 1

    def restore_initial_conditions(self):
        """
        Reset the initial conditions of the simulation.
        """
        # Do not read NHC/RNG data here to improve randomness.
        self.__simulation.restart_from(self.__initial_conditions,
                                       read_nhc_data=False,
                                       read_rng_state=False)

    @property
    def nep_parameters(self):
        """
        The keywords and corresponding values to be used in a nep.in file.
        """
        return self.__nep_parameters

    @nep_parameters.setter
    def nep_parameters(self, v: str):
        """
        Set the keywords and corresponding values to be used in a nep.in file.
        """
        self.__nep_parameters = v
        self.__check_nep_parameters()
