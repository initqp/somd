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
The transition path sampling method.
"""

import os as _os
import abc as _ab
import h5py as _h5
import numpy as _np
import shutil as _sh
import warnings as _w
from somd import apps as _mdapps
from somd import core as _mdcore
from somd import potentials as _pots
from . import utils as _path_utils

__all__ = ['TPSBASE', 'FELTWTPS', 'FELOWTPS']


class TPSBASE(_mdapps.simulations.STAGEDSIMULATION, _ab.ABC):
    """
    Base class of all transition path sampling protocols.

    Parameters
    ----------
    system : somd.core.systems.MDSYSTEM
        The simulated system.
    integrator : somd.core.integrator.INTEGRATOR
        The integrator that propagates the simulated system.
    potential_generators : List[callable]
        Generator functions of potential calculators.
    sampling_parameters : dict
        The parameters that define a sampling process. Valid keys of this
        dictionary are:
        - plumed_file : str
            Name of the plumed input file.
        - cv_names : List[dict]
            Names and components of the collective variables to save. For
            example: cv_names = [{'d1': 'x'}, {'d1': 'y'}, {'d2': ''}]
        - states : List[somd.apps.path_sampling.utils.state.STATE]
            The CV space states.
        - sampled_paths : List[List[int]]
            The paths to sample. For example: [[0, 1], [1, 2]].
        - initial_trajectory : str
            The initial trajectory.
        - bias_function : callable
            Default value : None
            The bias function used in the shooting move.  For example, to only
            select frames with a CV value between 1.0 and 2.0, set this option
            as: 'bias_function': lambda x: (x > 1.0) * (x < 2.0)
        - max_md_steps_per_iter : int
            Default value : 50000
            Maximum number of MD steps in each sampling iteration.
        - randomize_velocities : bool
            Default value : False
            If explicitly randomize shooting point velocities. This option is
            mandatory for deterministic integrators.
        - trajectory_interval : int
            Default value : 1
            The interval of writting the trajectoies.
        - trajectory_use_double : bool
            Default value : True
            If use 64-bit floating point in writing trajectories.
        - remove_dead_iterations : bool
            Default value : False
            If remove the non-reactive iterations to save disk spaces.
    post_step_objects : List[object]:
        The post step objects, including the barostat.
    output_prefix : str
        Prefix of the output file.
    """

    def __init__(self,
                 system: _mdcore.systems.MDSYSTEM,
                 integrator: _mdcore.integrators.INTEGRATOR,
                 potential_generators: list,
                 sampling_parameters: dict,
                 post_step_objects: list = [],
                 output_prefix: str = '') -> None:
        """
        Create a TPSBASE instance.
        """
        if (output_prefix == ''):
            output_prefix = 'path_sampling'
        self.__sampling_parameters = sampling_parameters
        self._check_sampling_parameters(integrator)
        super().__init__(system, integrator, potential_generators,
                         post_step_objects, output_prefix)

    def _check_sampling_parameters(self, integrator) -> None:
        """
        Check the sampling parameters.

        Parameters
        ----------
        integrator : somd.core.integrator.INTEGRATOR
            The integrator that propagates the simulated system.
        """
        param = self.__sampling_parameters
        # required parameters
        for key in ['plumed_file', 'cv_names', 'states', 'sampled_paths',
                    'initial_trajectory']:
            if (key not in param.keys()):
                raise KeyError('The key "{:s}" is required!'.format(key))
        for path in param['sampled_paths']:
            for state in path:
                if (state not in range(0, len(param['states']))):
                    message = 'Unknown state index {:d}'.format(state)
                    raise RuntimeError(message)
        param['plumed_file'] = _os.path.abspath(param['plumed_file'])
        # default parameters
        if ('bias_function' not in param.keys()):
            param['bias_function'] = None
        if ('randomize_velocities' not in param.keys()):
            param['randomize_velocities'] = False
        if ('max_md_steps_per_iter' not in param.keys()):
            param['max_md_steps_per_iter'] = 50000
        if ('trajectory_interval' not in param.keys()):
            param['trajectory_interval'] = 1
        if ('trajectory_use_double' not in param.keys()):
            param['trajectory_use_double'] = True
        if ('O' not in integrator.splitting_whole['operators']):
            param['randomize_velocities'] = True

    def _initialize_root_data_group(self) -> None:
        """
        Initialize the root data group in the output file.
        """
        if ('latest_accepted_iteration' in self.root['/iteration_data'].attrs):
            return
        h5_root = self.root['/']
        param = self.__sampling_parameters
        dtype_vlstr = _h5.special_dtype(vlen=str)
        h5_parameters = h5_root.create_group('parameters')
        n_cv = len(param['cv_names'])
        h5_cv_names = h5_parameters.create_dataset('cv_names', shape=(n_cv, 2),
                                                   dtype=dtype_vlstr)
        for i, cv in enumerate(param['cv_names']):
            h5_cv_names[i, 0] = list(cv.keys())[0]
            h5_cv_names[i, 1] = list(cv.values())[0]
        h5_cv_names.attrs['definition'] = ['name', 'component']
        h5_cv_names.attrs['plumed_file'] = param['plumed_file']
        n_states = len(param['states'])
        h5_states = h5_parameters.create_dataset('states', shape=(n_states,),
                                                 dtype=dtype_vlstr)
        for i, state in enumerate(param['states']):
            h5_states[i] = state.to_string()
        h5_states.attrs['labels'] = [s.label for s in param['states']]
        h5_sampled_paths = h5_parameters.create_group('sampled_paths')
        for i, path in enumerate(param['sampled_paths']):
            path_label = ''
            h5_sampled_paths.create_dataset(str(i), shape=(len(path),),
                                            dtype=int)
            for j in range(0, len(path)):
                h5_sampled_paths[str(i)][j] = path[j]
                path_label += param['states'][path[j]].label + '->'
            h5_sampled_paths[str(i)].attrs['label'] = path_label.strip('->')
        h5_root['iteration_data'].attrs['latest_accepted_iteration'] = 0
        self.root.flush()

    def _initialize_iteration_data_group(self, n_iter: int) -> None:
        """
        Initialize an iteration data group in the output file.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        """
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        if (h5_group.attrs['initialized']):
            return
        n_cv = self.root['/parameters/cv_names'].shape[0]
        h5_group.attrs['initialized'] = True
        h5_group.attrs['accepted'] = False
        h5_group.attrs['is_reactive'] = False
        h5_group.attrs['shooting_point_file'] = ''
        h5_group.attrs['parent_iteration'] = -1
        h5_group.attrs['path_index'] = -1
        h5_group.create_dataset('shooting_point_index', shape=(1,),
                                maxshape=(1,), dtype=int)
        h5_group['shooting_point_index'][0] = -1
        h5_group.create_dataset('shooting_point_cv_values', shape=(1, n_cv,),
                                maxshape=(1, n_cv,), dtype=_np.double)
        h5_group.create_dataset('metropolis_random_number', shape=(1,),
                                maxshape=(1,), dtype=_np.double)
        h5_group['metropolis_random_number'][0] = -1
        h5_group.create_dataset('selection_probability', shape=(1,),
                                maxshape=(1,), dtype=_np.double)
        h5_group['selection_probability'][0] = -1
        h5_group.create_dataset('selection_probability_reversed', shape=(1,),
                                maxshape=(1,), dtype=_np.double)
        # [[n_iter_0, n_seg_0], [n_iter_1, n_seg_1], ...]
        h5_group.create_dataset('path_segments', shape=(0, 2,),
                                maxshape=(None, 2,), dtype=int)
        h5_group.create_group('segments')
        self.root.flush()

    def _initialize_path_segment_data_group(self,
                                            n_iter: int,
                                            n_seg: int) -> None:
        """
        Initialize a path segment data group in the output file.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        n_seg : int
            Index of the segment
        """
        n_cv = self.root['/parameters/cv_names'].shape[0]
        h5_path = '/iteration_data/' + str(n_iter)
        h5_path_seg = h5_path + '/segments'
        if (str(n_seg) in self.root[h5_path_seg]):
            return
        working_dir = self.root[h5_path].attrs['working_directory'] + \
            '/segment_' + str(n_seg)
        if (_os.path.exists(working_dir)):
            _sh.rmtree(working_dir)
        _os.mkdir(working_dir)
        h5_group = self.root[h5_path_seg].create_group(str(n_seg))
        h5_group.attrs['direction'] = 1
        h5_group.attrs['trajectory_file_name'] = ''
        h5_group.attrs['system_data_file_name'] = ''
        h5_group.attrs['sink_state_indices'] = []
        h5_group.attrs['final_state_indices'] = []
        h5_group.attrs['shooting_point_file'] = \
            self.root[h5_path].attrs['shooting_point_file']
        h5_group.attrs['working_directory'] = working_dir
        h5_group.create_dataset('timestep', shape=(1), dtype=_np.double)
        h5_group.create_dataset('cv_values', shape=(0, n_cv),
                                maxshape=(None, n_cv), dtype=_np.double)
        h5_group['timestep'].attrs['units'] = 'picosecond'
        h5_group['timestep'][0] = -1.0

    def _load_path(self,
                   timestep: float,
                   trajectories: list,
                   path_directions: list,
                   reverse_velocities: list = None) -> _path_utils.path.PATH:
        """
        Load a path from trajectories.

        Parameters
        ----------
        timestep: float
            Timestep of the path.
        trajectories : List[str]
            Name of the trajectory files.
        path_directions : List[int]
            Directions of the trajectories.
        reverse_velocities : List[bool]
            If reverse the velocities of the trajectories when reading.
        """
        paths = []
        if (len(trajectories) != len(path_directions)):
            message = 'Mismatch between the number of trajectories and ' + \
                      'the length of direction list!'
            raise RuntimeError(message)
        if (reverse_velocities is None):
            reverse_velocities = [d == -1 for d in path_directions]
        for i in range(0, len(trajectories)):
            p = _path_utils.path.PATH(self.system, [trajectories[i]], timestep,
                                      None, [reverse_velocities[i]])
            if (path_directions[i] < 0):
                paths.append(reversed(p))
            else:
                paths.append(p)
        return _path_utils.path.concatenate(paths)

    def _load_path_from_iter(self, n_iter: int) -> _path_utils.path.PATH:
        """
        Load the path from one iteration.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        """
        timesteps = []
        directions = []
        trajectories = []
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        for seg in h5_group['path_segments']:
            if (seg[1] != -1):
                # Regular path segments.
                h5_path_seg = '/iteration_data/{:d}/segments/{:d}'
                h5_group_seg = self.root[h5_path_seg.format(*seg)]
                h5_attrs = h5_group_seg.attrs
                timesteps.append(h5_group_seg['timestep'][0])
                directions.append(h5_attrs['direction'])
                trajectories.append(h5_attrs['trajectory_file_name'])
            else:
                # The shooting point.
                directions.append(1)
                trajectories.append(h5_group.attrs['shooting_point_file'])
        if (len(set(timesteps)) != 1):
            message = 'The reactive path segments in iteration {:d} ' + \
                      'contains various timesteps!'
            raise RuntimeError(message.format(n_iter))
        path = self._load_path(timesteps[0], trajectories, directions)
        return path

    def _load_cv_values_from_iter(self, n_iter: int) -> _np.ndarray:
        """
        Load CV values of the path from one iteration.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        """
        cv_values = []
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        for seg in h5_group['path_segments']:
            if (seg[1] != -1):
                # Regular path segments.
                h5_path_seg = '/iteration_data/{:d}/segments/{:d}'
                h5_group_seg = self.root[h5_path_seg.format(*seg)]
                h5_cv_array = h5_group_seg['cv_values']
                cv_values_tmp = _np.array(h5_cv_array)
                if (h5_group_seg.attrs['direction'] < 0):
                    cv_values.append(cv_values_tmp[::-1])
                else:
                    cv_values.append(cv_values_tmp)
            else:
                # The shooting point.
                h5_cv_array = h5_group['shooting_point_cv_values']
                cv_values.append(_np.array(h5_cv_array))
        return _np.concatenate(cv_values)

    def run_segment(self, n_iter: int, n_seg: int) -> None:
        """
        Propagate a initialized trajectory segment.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        n_seg : int
            Index of the segment
        """
        cwd = _os.getcwd()
        param = self.__sampling_parameters
        h5_path = '/iteration_data/{:d}/segments/{:d}'
        h5_group = self.root[h5_path.format(n_iter, n_seg)]
        _os.chdir(h5_group.attrs['working_directory'])
        # Set up the potentials and simulation.
        atom_list = list(range(0, self.system.n_atoms))
        plumed_file = self.root['/parameters/cv_names'].attrs['plumed_file']
        plumed = _pots.PLUMED(atom_list, plumed_file,
                              timestep=self.integrator.timestep,
                              cv_names=self.cv_names, output_prefix='')
        potential_list = list(range(0, len(self.potential_generators)))
        with self._set_up_simulation(potential_list, [plumed]) as simulation:
            read_nhc_data = bool(len(self.integrator._nhchains) > 0)
            simulation.restart_from(h5_group.attrs['shooting_point_file'],
                                    read_nhc_data=read_nhc_data,
                                    read_rng_state=False, read_forces=False)
            # NOTE: we are using reversed velocities instead of reversed time
            # here for compatible of stochastic integrators. Thus when the
            # direction of the trajectory is `-1`, the velocities recorded in
            # the trajectory is reversed.
            if (h5_group.attrs['direction'] < 0):
                simulation.system.velocities[:] *= -1
            h5_group['timestep'][0] = param['trajectory_interval'] * \
                simulation.integrator.timestep
            # Set up writers
            if (_os.path.exists(h5_group.attrs['system_data_file_name'])):
                _os.remove(h5_group.attrs['system_data_file_name'])
            if (_os.path.exists(h5_group.attrs['trajectory_file_name'])):
                _os.remove(h5_group.attrs['trajectory_file_name'])
            data_writer = _mdapps.loggers.DEFAULTCSVLOGGER(
                h5_group.attrs['system_data_file_name'],
                interval=param['trajectory_interval'])
            data_writer.bind_integrator(simulation.integrator)
            traj_writer = _mdapps.trajectories.H5WRITER(
                h5_group.attrs['trajectory_file_name'],
                interval=param['trajectory_interval'],
                write_velocities=True,
                use_double=param['trajectory_use_double'])
            traj_writer.bind_integrator(simulation.integrator)
            simulation.post_step_objects.append(data_writer)
            simulation.post_step_objects.append(traj_writer)
            # Run.
            cv_values = []
            max_steps = param['max_md_steps_per_iter']
            sink_states = [self.states[i] for i in
                           h5_group.attrs['sink_state_indices']]
            for i in range(0, max_steps):
                simulation.run(1)
                if (((i + 1) % param['trajectory_interval']) == 0):
                    cv_values.append(plumed.cv_values)
                    stop_flags = [s(plumed.cv_values) for s in sink_states]
                    if ((True in stop_flags) or (i == (max_steps - 1))):
                        results = [h5_group.attrs['sink_state_indices'][j]
                                   for j, f in enumerate(stop_flags) if f]
                        break
            h5_group.attrs['final_state_indices'] = results
            shape = (len(cv_values), h5_group['cv_values'].shape[1],)
            h5_group['cv_values'].resize(shape)
            h5_group['cv_values'][:] = cv_values
            self.root.flush()
        _os.chdir(cwd)

    def get_segment_from_iter(self,
                              target_iter: int,
                              parent_iter: int,
                              n_seg: int,
                              frame_slice: slice = None,
                              write_trajectory: bool = False) -> None:
        """
        Obtain trajectory segment by slicing the parent path.

        Parameters
        ----------
        target_iter : int
            Index of the target iteration.
        parent_iter : int
            Index of the parent iteration.
        n_seg : int
            Index of the segment
        frame_slice : slice
            Range of the required frames.
        write_trajectory : bool
            If write the path segment file.
        """
        param = self.__sampling_parameters
        h5_path = '/iteration_data/{:d}/segments/{:d}'
        h5_group = self.root[h5_path.format(target_iter, n_seg)]
        h5_group['timestep'][0] = param['trajectory_interval'] * \
            self.integrator.timestep
        cv_values = self._load_cv_values_from_iter(parent_iter)
        if (frame_slice is not None):
            cv_values = cv_values[frame_slice]
        shape = (len(cv_values), h5_group['cv_values'].shape[1],)
        h5_group['cv_values'].resize(shape)
        h5_group['cv_values'][:] = cv_values
        if (write_trajectory):
            if (_os.path.exists(h5_group.attrs['trajectory_file_name'])):
                _os.remove(h5_group.attrs['trajectory_file_name'])
            path = self._load_path_from_iter(parent_iter)
            if (frame_slice is not None):
                path = path[frame_slice]
            path.dump(h5_group.attrs['trajectory_file_name'],
                      param['trajectory_use_double'])
        self.root.flush()

    def _shoot(self,
               path: _path_utils.path.PATH,
               cv_values: _np.ndarray,
               shooting_point_file: str) -> dict:
        """
        Shoot from one old path.

        Parameters
        ----------
        path : somd.apps.path_sampling.utils.path.PATH
            The reference path.
        cv_values : numpy.ndarray
            The CV values of the path.
        shooting_point_file : str
            Name of the shooting point file.
        """
        param = self.__sampling_parameters
        shooting_point = self.system.copy()
        shooting_dict = _path_utils.selection.select(
            cv_values, param['bias_function'], True)
        shooting_point.snapshot = path.frame(shooting_dict['index'])
        if (param['randomize_velocities']):
            _path_utils.shooting.shoot(shooting_point)
        if (self.integrator._nhchains != []):
            self.integrator._randomize_nhchains_states()
        with self._set_up_simulation([], [_pots.FAKE([])]) as simulation:
            simulation.system.positions[:] = shooting_point.positions[:]
            simulation.system.velocities[:] = shooting_point.velocities[:]
            simulation.dump_restart(shooting_point_file)
        return shooting_dict

    def shoot_from_iter(self,
                        target_iter: int,
                        parent_iter: int = None) -> dict:
        """
        Shoot from the one given iteration.

        Parameters
        ----------
        target_iter : int
            Index of the target iteration.
        parent_iter : int
            Load the parent path from this iteration.
        """
        h5_iter_root = self.root['iteration_data']
        if (parent_iter is None):
            parent_iter = h5_iter_root.attrs['latest_accepted_iteration']
        h5_group = h5_iter_root[str(target_iter)]
        working_dir = h5_group.attrs['working_directory']
        file_name = working_dir + '/shooting_point.h5'
        path = self._load_path_from_iter(parent_iter)
        cv_values = self._load_cv_values_from_iter(parent_iter)
        result = self._shoot(path, cv_values, file_name)
        h5_group.attrs['parent_iteration'] = parent_iter
        h5_group.attrs['shooting_point_file'] = file_name
        h5_group['shooting_point_index'][0] = result['index']
        h5_group['selection_probability'][0] = result['probability']
        h5_group['shooting_point_cv_values'][:] = cv_values[result['index']]
        self.root.flush()

    def calculate_p_reversed(self, n_iter: int) -> None:
        """
        Calculate the reversed selection probability: $p_{sel}(x'_j|X')$
        of one iteration.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        """
        param = self.__sampling_parameters
        cv_values = self._load_cv_values_from_iter(n_iter)
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        h5_path_cv = '/iteration_data/{:d}/segments/{:d}/cv_values'
        shooting_point_index = 0
        for seg in h5_group['path_segments']:
            if (seg[1] != -1):
                shooting_point_index += \
                    self.root[h5_path_cv.format(*seg)].shape[0]
            else:
                break
        h5_group['selection_probability_reversed'][0] = \
            _path_utils.selection.calculate_frame_probabilities(
                cv_values, param['bias_function'], True)[shooting_point_index]

    def metropolis(self, n_iter: int) -> None:
        """
        Perform the Metropolis step.

        Parameters
        ----------
        n_iter : int
            Index of the iteration.
        """
        trial = _np.random.rand()
        h5_iter_root = self.root['iteration_data']
        h5_group = h5_iter_root[str(n_iter)]
        h5_group['metropolis_random_number'][0] = trial
        p_o = h5_group['selection_probability'][0]
        p_n = h5_group['selection_probability_reversed'][0]
        parent_iter = h5_group.attrs['parent_iteration']
        if ((p_n <= 0 or p_o <= 0) and parent_iter != 0):
            message = 'Selection probability or reversed selection ' + \
                      'probability in iteration {:d} is smaller than 0! ' + \
                      'Can not perform the Metropolis step!'
            raise RuntimeError(message.format(n_iter))
        elif (parent_iter == 0 or ((p_n / p_o) > trial)):
            # When shooting from the initial trajectory, we should always
            # accepte the move.
            h5_group.attrs['accepted'] = True
            old_iter = h5_iter_root.attrs['latest_accepted_iteration']
            if (n_iter > old_iter):
                h5_iter_root.attrs['latest_accepted_iteration'] = n_iter
        self.root.flush()

    def analyze_path(self, n_iter: int) -> None:
        """
        Analyze which kind of path does one iteration sample.

        Parameters
        ----------
        n_iter : int
            Index the iteration.
        """
        h5_group = self.root['/iteration_data/' + str(n_iter)]
        cv_values = self._load_cv_values_from_iter(n_iter)
        for i, state_seq in enumerate(self.sampled_paths):
            states = [self.states[i] for i in state_seq]
            if (_path_utils.state.in_sequence(states, cv_values)):
                h5_group.attrs['is_reactive'] = True
                h5_group.attrs['path_index'] = i
                self.root.flush()
                return

    def make_initial_path(self) -> None:
        """
        Set up the initial path.
        """
        param = self.__sampling_parameters
        h5_group = self.root['/iteration_data/0']
        working_dir = h5_group.attrs['working_directory']
        timestep = self.integrator.timestep * param['trajectory_interval']
        path = self._load_path(timestep, [param['initial_trajectory']], [1])
        cv_values = path.calculate_cv_values(param['plumed_file'],
                                             self.cv_names)
        self._initialize_path_segment_data_group(0, 0)
        h5_group_seg = self.root['/iteration_data/0/segments/0']
        working_dir = working_dir + '/segment_0'
        for state_indices in self.sampled_paths:
            reactive_regions = _path_utils.state.get_reactive_regions(
                [self.states[state_indices[0]],
                 self.states[state_indices[-1]]], cv_values)
            if (reactive_regions != []):
                file_name = working_dir + '/segment.trajectory.h5'
                region_lengths = [_np.abs(r['range'][0] - r['range'][1])
                                  for r in reactive_regions]
                region = reactive_regions[_np.argsort(region_lengths)[-1]]
                path_r = path[region['range'][0]:region['range'][1]]
                if (region['direction'] < 0):
                    path_r = reversed(path_r)
                h5_group.attrs['accepted'] = True
                h5_group.attrs['is_reactive'] = True
                h5_group_seg.attrs['direction'] = 1
                h5_group_seg.attrs['trajectory_file_name'] = file_name
                h5_group_seg.attrs['final_state_indices'] = [state_indices[-1]]
                h5_group['path_segments'].resize((1, 2,))
                h5_group['path_segments'][0, 0] = 0
                h5_group['path_segments'][0, 1] = 0
                shape = (len(path_r), h5_group_seg['cv_values'].shape[1],)
                h5_group_seg['cv_values'].resize(shape)
                h5_group_seg['timestep'][0] = timestep
                # Recalculate the CV values.
                cv_values = path_r.calculate_cv_values(param['plumed_file'],
                                                       self.cv_names)
                h5_group_seg['cv_values'][:] = cv_values
                path_r.dump(file_name,
                            use_double=param['trajectory_use_double'])
                self.analyze_path(0)
                self.root.flush()
                break
        if (not h5_group.attrs['accepted']):
            message = 'Initial trajectory {:s} is not reactive!'
            raise RuntimeError(message.format(param['initial_trajectory']))

    @_ab.abstractmethod
    def step(self) -> None:
        """
        Run the sampling by one step.
        """
        raise NotImplementedError()

    def initialize(self) -> None:
        """
        Initialize the simulation.
        """
        if (self.n_iter == 0):
            self._initialize_root_data_group()
            self._initialize_iteration_data_group(0)
            self.make_initial_path()
            self._increase_n_iter()

    def run(self, n_iterations: int = 1) -> None:
        """
        Run the simulation.

        Parameters
        ----------
        n_iterations : int
            Number of iterations to run.
        """
        self.initialize()
        for i in range(0, n_iterations):
            self._initialize_iteration_data_group(self.n_iter)
            self.step()
            self._increase_n_iter()

    @property
    def cv_names(self) -> list:
        """
        Names of the CV.
        """
        return [{cv[0].decode('UTF-8'): cv[1].decode('UTF-8')}
                for cv in self.root['/parameters/cv_names']]

    @property
    def states(self) -> list:
        """
        State definitions.
        """
        if (not hasattr(self, '__states')):
            from_string = _path_utils.state.STATE.from_string
            self.__states = [from_string(s.decode('UTF-8'))
                             for s in self.root['/parameters/states']]
        return self.__states

    @property
    def sampled_paths(self) -> list:
        """
        The sampled paths.
        """
        return [[s for s in self.root['/parameters/sampled_paths'][p]]
                for p in self.root['/parameters/sampled_paths']]

    @property
    def sampling_parameters(self) -> dict:
        """
        The sampling parameters.
        """
        return self.__sampling_parameters


class _TWOSEGMENTTPS(TPSBASE):
    """
    Standard two-segment TPS.
    """

    def _get_segment_mapping(self) -> list:
        """
        Return the path segment mapping.
        """
        return [[self.n_iter, 0], [self.n_iter, -1], [self.n_iter, 1]]

    def forward(self, n_iter: int) -> None:
        """
        Run the simulation forwardly.

        Parameters
        ----------
        n_iter : int
            Index the iteration.
        """
        self._initialize_path_segment_data_group(n_iter, 1)
        h5_path = '/iteration_data/{:d}/segments/{:d}'.format(n_iter, 1)
        h5_group = self.root[h5_path]
        working_dir = h5_group.attrs['working_directory']
        h5_group.attrs['direction'] = 1
        h5_group.attrs['sink_state_indices'] = self.sink_state_indices
        h5_group.attrs['trajectory_file_name'] = \
            working_dir + '/segment.trajectory.h5'
        h5_group.attrs['system_data_file_name'] = \
            working_dir + '/segment.data.csv'
        self.root.flush()
        self.run_segment(n_iter, 1)

    def backward(self, n_iter: int) -> None:
        """
        Run the simulation backwardly.
        """
        self._initialize_path_segment_data_group(n_iter, 0)
        h5_path = '/iteration_data/{:d}/segments/{:d}'.format(n_iter, 0)
        h5_group = self.root[h5_path]
        working_dir = h5_group.attrs['working_directory']
        h5_group.attrs['direction'] = -1
        h5_group.attrs['sink_state_indices'] = self.sink_state_indices
        h5_group.attrs['trajectory_file_name'] = \
            working_dir + '/segment.trajectory.h5'
        h5_group.attrs['system_data_file_name'] = \
            working_dir + '/segment.data.csv'
        self.root.flush()
        self.run_segment(n_iter, 0)

    def forward_from_iter(self,
                          target_iter: int,
                          parent_iter: int,
                          shooting_point_index: int,
                          write_trajectory: bool = False) -> None:
        """
        Obtain forward trajectory segment by slicing the parent path.

        Parameters
        ----------
        target_iter : int
            Index of the target iteration.
        shooting_point_index : int
            Index of the shooting point.
        frame_slice : slice
            Range of the required frames.
        write_trajectory : bool
            If write the path segment file.
        """
        self._initialize_path_segment_data_group(target_iter, 1)
        h5_path = '/iteration_data/{:d}/segments/{:d}'.format(target_iter, 1)
        h5_group = self.root[h5_path]
        working_dir = h5_group.attrs['working_directory']
        h5_group.attrs['direction'] = 1
        h5_group.attrs['trajectory_file_name'] = \
            working_dir + '/segment.trajectory.h5'
        self.root.flush()
        frame_slice = slice((shooting_point_index + 1), None)
        self.get_segment_from_iter(target_iter, parent_iter, 1, frame_slice,
                                   write_trajectory)

    def backward_from_iter(self,
                           target_iter: int,
                           parent_iter: int,
                           shooting_point_index: int,
                           write_trajectory: bool = False) -> None:
        """
        Obtain backward trajectory segment by slicing the parent path.

        Parameters
        ----------
        target_iter : int
            Index of the target iteration.
        parent_iter : int
            Index of the parent iteration.
        shooting_point_index : int
            Index of the shooting point.
        write_trajectory : bool
            If write the path segment file.
        """
        self._initialize_path_segment_data_group(target_iter, 0)
        h5_path = '/iteration_data/{:d}/segments/{:d}'.format(target_iter, 0)
        h5_group = self.root[h5_path]
        working_dir = h5_group.attrs['working_directory']
        # NOTE: this is not a bug, since the loaded paths are always
        # forwardly ordered.
        h5_group.attrs['direction'] = 1
        h5_group.attrs['trajectory_file_name'] = \
            working_dir + '/segment.trajectory.h5'
        self.root.flush()
        frame_slice = slice(0, shooting_point_index)
        self.get_segment_from_iter(target_iter, parent_iter, 0, frame_slice,
                                   write_trajectory)

    @property
    def sink_state_indices(self) -> list:
        """
        Indices of the sink states.
        """
        if (not hasattr(self, '__sink_state_indicies')):
            result = []
            for state_seq in self.sampled_paths:
                for i in [0, 1]:
                    result.append(state_seq[i])
            self.__sink_state_indicies = list(set(result))
        return self.__sink_state_indicies


class FELTWTPS(_TWOSEGMENTTPS):
    """
    Fexible Length (FEL) Two-Way (TW) TPS.
    """

    def step(self) -> None:
        """
        Run the sampling by one step.
        """
        # Prepare states and path definitions.
        h5_group = self.root['/iteration_data/' + str(self.n_iter)]
        h5_group['path_segments'].resize((3, 2,))
        h5_group['path_segments'][:] = self._get_segment_mapping()
        # Run the iteration.
        self.shoot_from_iter(self.n_iter)
        self.backward(self.n_iter)
        self.forward(self.n_iter)
        self.analyze_path(self.n_iter)
        if (h5_group.attrs['is_reactive']):
            self.calculate_p_reversed(self.n_iter)
            self.metropolis(self.n_iter)
        elif (self.sampling_parameters['remove_dead_iterations']):
            _sh.rmtree(h5_group.attrs['working_directory'])


class FELOWTPS(_TWOSEGMENTTPS):
    """
    Fexible Length (FEL) One-Way (OW) TPS.
    """

    def _check_sampling_parameters(self, integrator) -> None:
        """
        Check the sampling parameters.

        Parameters
        ----------
        integrator : somd.core.integrator.INTEGRATOR
            The integrator that propagates the simulated system.
        """
        param = self.sampling_parameters
        if ('randomize_velocities' in param.keys()):
            if (param['randomize_velocities']):
                message = 'The "randomize_velocities" option has been ' + \
                          'automatically disabled!'
                _w.warn(message)
        param['randomize_velocities'] = False
        if ('O' not in integrator.splitting_whole['operators']):
            message = 'The one-way TPS requires a stochastic integrator!'
            raise RuntimeError(message)
        super()._check_sampling_parameters(integrator)

    def step(self) -> None:
        """
        Run the sampling by one step.
        """
        # Get parent iteration.
        h5_iter_root = self.root['/iteration_data']
        parent_iter = h5_iter_root.attrs['latest_accepted_iteration']
        # Prepare states and path definitions.
        h5_group = self.root['/iteration_data/' + str(self.n_iter)]
        h5_group['path_segments'].resize((3, 2,))
        h5_group['path_segments'][:] = self._get_segment_mapping()
        # Determine the direction.
        segment_direction = _np.sign(_np.random.rand() - 0.5)
        h5_group.attrs['new_segment_direction'] = segment_direction
        # Run the iteration.
        self.shoot_from_iter(self.n_iter)
        if (segment_direction < 0):
            self.forward_from_iter(self.n_iter, parent_iter,
                                   h5_group['shooting_point_index'][0], True)
            self.backward(self.n_iter)
        else:
            self.backward_from_iter(self.n_iter, parent_iter,
                                    h5_group['shooting_point_index'][0], True)
            self.forward(self.n_iter)
        self.analyze_path(self.n_iter)
        if (h5_group.attrs['is_reactive']):
            self.calculate_p_reversed(self.n_iter)
            self.metropolis(self.n_iter)
        elif (self.sampling_parameters['remove_dead_iterations']):
            _sh.rmtree(h5_group.attrs['working_directory'])
