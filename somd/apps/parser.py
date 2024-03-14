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
A simple TOML input paser.
"""

import os as _os
import numpy as _np
import typing as _tp
from collections import namedtuple as _namedtuple
from somd import core as _mdcore
from somd import apps as _mdapps
from somd import utils as _mdutils
from somd import potentials as _potentials

__all__ = ['TOMLPARSER']


class TOMLPARSER(object):
    """
    Set up a MD run using a TOML [1] configure file.

    Parameters
    ----------
    file_name: str
        Name of the configure file.

    References
    ----------
    .. [1] https://toml.io
    """

    # Dependency.
    __dep__ = _namedtuple('DEPENDENCY', ['key', 'values'])
    # Parameter table value type.
    __value__ = _namedtuple('VALUE', ['type', 'required', 'dependency'])
    # Parameter table value type.
    __table__ = _namedtuple('TABLE', ['is_array', 'required'])
    # Valid tables.
    __tables__ = {
        'system': __table__(False, True),
        'integrator': __table__(False, True),
        'potential': __table__(True, True),
        'barostat': __table__(False, False),
        'constraints': __table__(False, False),
        'trajectory': __table__(True, False),
        'logger': __table__(True, False),
        'group': __table__(True, False),
        'script': __table__(True, False),
        'active_learning': __table__(False, False),
        'run': __table__(False, False),
    }
    # Valid key-value pairs of each table.
    __parameters__ = dict()
    __parameters__['run'] = {
        'n_steps': __value__([int], True, None),
        'seed': __value__([int], False, None),
        'label': __value__([str], False, None),
        'restart_from': __value__([str], False, None),
    }
    __parameters__['system'] = {
        'structure': __value__([str], True, None),
        'format': __value__([str], False, None),
    }
    __parameters__['integrator'] = {
        'timestep': __value__([float], True, None),
        'type': __value__([str], False, None),
        'splitting': __value__([list], False, None),
        'temperatures': __value__([list, float], False, None),
        'relaxation_times': __value__([list, float], False, None),
        'thermo_groups': __value__([list, int], False, None),
    }
    __parameters__['potential'] = {
        'type': __value__([str], True, None),
        'atom_list': __value__([list, str], False, None),
        'siesta_options': __value__([str], True, __dep__('type', ['siesta'])),
        'siesta_command': __value__([str], True, __dep__('type', ['siesta'])),
        'pseudopotential_dir': __value__(
            [str], False, __dep__('type', ['siesta'])
        ),
        'functional': __value__(
            [str], True, __dep__('type', ['dftd3', 'dftd4'])
        ),
        'damping': __value__([str], False, __dep__('type', ['dftd3'])),
        'atm': __value__([bool], False, __dep__('type', ['dftd3', 'dftd4'])),
        'total_charge': __value__([int], False, __dep__('type', ['dftd4'])),
        'file_name': __value__(
            [str], True, __dep__('type', ['plumed', 'nep', 'mace'])
        ),
        'use_tabulating': __value__([bool], False, __dep__('type', ['nep'])),
        'device': __value__([str], False, __dep__('type', ['mace'])),
        'energy_unit': __value__([float], False, __dep__('type', ['mace'])),
        'length_unit': __value__([float], False, __dep__('type', ['mace'])),
        'model_dtype': __value__([str], False, __dep__('type', ['mace'])),
        'virial': __value__([bool], False, __dep__('type', ['mace'])),
        'charge_cv_expr': __value__([str], False, __dep__('type', ['mace'])),
        'extra_cv_potential_index': __value__(
            [int], False, __dep__('type', ['plumed'])
        ),
    }
    __parameters__['group'] = {
        'atom_list': __value__([list, str], True, None),
        'label': __value__([str], False, None),
        'has_translations': __value__([bool], False, None),
        'initial_velocities': __value__([list], False, None),
        'initial_temperature': __value__([float], False, None),
    }
    __parameters__['barostat'] = {
        'pressures': __value__([list, float], True, None),
        'beta': __value__([list, float], True, None),
        'relaxation_time': __value__([float], True, None),
    }
    __parameters__['constraints'] = {
        'types': __value__([list], True, None),
        'indices': __value__([list], True, None),
        'targets': __value__([list], True, None),
        'tolerances': __value__([list], False, None),
        'max_cycles': __value__([int], False, None),
    }
    __parameters__['trajectory'] = {
        'prefix': __value__([str], False, None),
        'format': __value__([str], False, None),
        'interval': __value__([int], False, None),
        'write_forces': __value__([bool], False, None),
        'write_velocities': __value__([bool], False, None),
        'wrap_positions': __value__([bool], False, None),
        'potential_list': __value__([list], False, None),
        'use_float64': __value__([bool], False, __dep__('format', ['h5'])),
        'energy_shift': __value__([float], False, __dep__('format', ['exyz'])),
        'is_restart_file': __value__([bool], False, __dep__('format', ['h5'])),
    }
    __parameters__['logger'] = {
        'format': __value__([str], False, None),
        'prefix': __value__([str], False, None),
        'interval': __value__([int], False, None),
        'potential_list': __value__([list], False, None),
    }
    __parameters__['script'] = {
        'update': __value__([str], True, None),
        'interval': __value__([int], False, None),
        'initialize': __value__([str], False, None),
    }
    __parameters__['active_learning'] = {
        'nep_options': __value__([str], True, None),
        'nep_command': __value__([str], True, None),
        'initial_training_set': __value__([str], True, None),
        'n_iterations': __value__([int], True, None),
        'n_potentials': __value__([int], False, None),
        'max_md_runs_per_iter': __value__([int], False, None),
        'max_md_steps_per_iter': __value__([int], False, None),
        'trajectory_interval': __value__([int], False, None),
        'msd_lower_limit': __value__([float], False, None),
        'msd_upper_limit': __value__([float], False, None),
        'min_new_structures_per_iter': __value__([int], False, None),
        'max_new_structures_per_iter': __value__([int], False, None),
        'initial_potential_files': __value__([list], False, None),
        'initial_testing_set': __value__([str], False, None),
        'reference_potentials': __value__([list], False, None),
        'use_tabulating': __value__([bool], False, None),
        'energy_shift': __value__([float], False, None),
    }

    def __init__(self, file_name: str) -> None:
        """
        Create a TOMLPARSER instance.
        """
        import toml

        self.__file_name = file_name
        self.__root = toml.load(file_name)
        self.__normalize_tables()
        self.__check_task()
        self.__parse_run()
        self.__parse_system()
        self.__parse_groups()
        self.__parse_constraints()
        self.__add_group_velocities()
        self.__parse_integrator()
        if self.__integrator._is_nve:
            self.__parse_potentials(self.__integrator.timestep, [], [])
        else:
            self.__parse_potentials(
                self.__integrator.timestep,
                self.__integrator.temperatures,
                self.__integrator.thermo_groups,
            )
        self.__parse_barostat()
        self.__parse_loggers()
        self.__parse_trajectories()
        self.__parse_scripts()
        if self.__root['active_learning'] is None:
            self.__set_up_simulation()
        else:
            self.__parse_active_learning()

    def __normalize_tables(self) -> None:
        """
        Check root tables.
        """
        bad_tables = set(self.__root.keys()).difference(self.__tables__.keys())
        if len(bad_tables) != 0:
            format_list = ('[{}], ' * len(bad_tables)).strip(', ')
            message = 'Unknown table name(s): {}!'.format(format_list)
            raise KeyError(message.format(*bad_tables))
        required_tables = set(
            k for k in self.__tables__.keys() if self.__tables__[k].required
        )
        missing_tables = required_tables.difference(set(self.__root.keys()))
        if len(missing_tables) != 0:
            format_list = ('[{}], ' * len(missing_tables)).strip(', ')
            message = 'Table(s) {} are required!'.format(format_list)
            raise KeyError(message.format(*missing_tables))
        definitions = dict.fromkeys(self.__tables__.keys(), None)
        for key in self.__root.keys():
            table = self.__root[key]
            is_array = self.__tables__[key].is_array
            if is_array and (not isinstance(table, list)):
                message = 'Name "{}" should define AN ARRAY of tables!'
                raise TypeError(message.format(key))
            if is_array and (not all(isinstance(t, dict) for t in table)):
                message = 'Name "{}" should define an array of TABLES!'
                raise TypeError(message.format(key))
            if (not is_array) and (not isinstance(table, dict)):
                message = 'Name "{}" should define a table!'
                raise TypeError(message.format(key))
            definitions[key] = table
        self.__root = definitions

    def __check_task(self) -> None:
        """
        Check the simulation task.
        """
        task_tables = ['run', 'active_learning']
        task_table_list = [self.__root[i] for i in task_tables]
        if all((i is None) for i in task_table_list):
            message = 'One of the following tables must present: {}'
            raise RuntimeError(message.format(task_tables))

    def __normalize_table(self, inp: dict, table_name: str) -> dict:
        """
        Check keys and their types of a table.

        Parameters
        ----------
        inp : dict
            The input parameter table.
        table_name : str
            Name of the table.
        """
        parameters = self.__parameters__[table_name].copy()
        bad_keys = set(inp.keys()).difference(parameters.keys())
        if len(bad_keys) != 0:
            format_list = ('"{}", ' * len(bad_keys)).strip(', ')
            message = 'Unknown key(s): {} in (one of) the [{}] table(s)!'
            message = message.format(format_list, table_name)
            raise KeyError(message.format(*bad_keys))
        # Check key names and value types.
        definitions = dict.fromkeys(parameters.keys(), None)
        for key in inp.keys():
            if type(inp[key]) not in parameters[key].type:
                types = parameters[key].type
                types = ('"{}"/' * len(types)).format(*types).strip('/')
                message = (
                    'Wrong type of key "{}" in (one of) the [{}] '
                    + 'table(s)! A {} typed value is expected!'
                )
                message = message.format(key, table_name, types)
                raise TypeError(message)
            definitions[key] = inp[key]
        # Check dependencies for empty keys.
        keys = definitions.keys()
        for key in [k for k in keys if definitions[k] is None]:
            required = parameters[key].required
            dependency = parameters[key].dependency
            if required and dependency is None:
                if definitions[key] is None:
                    message = 'Key "{}" is required in the [{}] table(s)!'
                    raise KeyError(message.format(key, table_name))
            elif required and dependency is not None:
                value = definitions[dependency.key].lower()
                if value in dependency.values:
                    message = (
                        'Key "{}" in (one of) the [{}] table(s) is required '
                        + 'since the value of the "{}" key is "{}"!'
                    )
                    message = message.format(
                        key,
                        table_name,
                        dependency.key,
                        definitions[dependency.key],
                    )
                    raise KeyError(message)
        # Check dependencies for non-empty keys.
        for key in [k for k in keys if definitions[k] is not None]:
            dependency = parameters[key].dependency
            if dependency is not None:
                value = definitions[dependency.key].lower()
                if value not in dependency.values:
                    message = (
                        'Key "{}" in (one of) the [{}] table(s) is redundant '
                        + 'since the value of the "{}" key is "{}"!'
                    )
                    message = message.format(
                        key,
                        table_name,
                        dependency.key,
                        definitions[dependency.key],
                    )
                    raise KeyError(message)
        return definitions

    def __parse_atom_list(self, inp: _tp.Union[list, str]) -> list:
        """
        Parse the atom selections.

        Parameters
        ----------
        inp: list(int) or str
            The atom list to be parsed.
        """
        if isinstance(inp, list):
            if not all(isinstance(i, int) for i in inp):
                raise RuntimeError('Unknown atom indices: "{}"'.format(inp))
            result = inp
        elif isinstance(inp, str):
            if inp.lower() == 'all':
                result = list(range(0, self.__system.n_atoms))
            else:
                result = []
                for s in inp.split(','):
                    if ':' not in s:
                        if not s.isnumeric():
                            message = 'Unknown atom index: "{}"'.format(s)
                            raise SyntaxError(message)
                        result.append(int(s))
                    else:
                        l = s.split(':')
                        if (
                            not l[0].isnumeric()
                            or len(l) != 2
                            or not l[1].isnumeric()
                        ):
                            message = 'Unknown atom range: "{}"'.format(s)
                            raise SyntaxError(message)
                        for i in range(int(l[0]), int(l[1]) + 1):
                            result.append(i)
        else:
            message = 'Type of key "atom_list" could only be List[int] or str!'
            raise RuntimeError(message)
        return result

    def __check_groups(self) -> None:
        """
        Check group setups of a system.
        """
        whole_system_flag = False
        no_translations_flag = False
        # The user has defined the group without translations.
        for group in self.__system.groups:
            if group.n_atoms == self.__system.n_atoms:
                whole_system_flag = True
            if group.has_translations is False:
                no_translations_flag = True
        # No group is corresponding to the whole group. Add a new group.
        if whole_system_flag is False:
            label = 'GROUP{:d}'.format(len(self.__system.groups))
            d = {
                'atom_list': range(0, self.__system.n_atoms),
                'has_translations': True,
                'label': label,
            }
            self.__system.groups.create_from_dict(d)
            message = (
                'An atom group that corresponds to the whole '
                + 'system has been append the group list.'
            )
            self.__root['group'].append(d)
            _mdutils.warning.warn(message)
        # The user has not defined the group without translations.
        # Check if there is any group that is corresponding to the whole group.
        if no_translations_flag is False:
            atom_numbers = [g.n_atoms for g in self.__system.groups]
            index = atom_numbers.index(self.__system.n_atoms)
            group = self.__system.groups[index]
            self.__root['group'][index]['has_translations'] = False
            group.has_translations = False
            message = (
                'Translational degrees of freedom of group '
                + '"{}" has been automatically removed.'
            )
            _mdutils.warning.warn(message.format(group._label))

    def __parse_run(self) -> None:
        """
        Parse the simulation task information.
        """
        table = self.__root['run']
        if table is None:
            table = dict()
            label = _os.path.splitext(self.__file_name)[0]
            table['label'] = _os.path.basename(label)
            table['restart_from'] = None
        else:
            table = self.__normalize_table(self.__root['run'], 'run')
            if table['seed'] is not None:
                _np.random.seed(table['seed'])
            if table['label'] is None:
                label = _os.path.splitext(self.__file_name)[0]
                table['label'] = _os.path.basename(label)
        self.__root['run'] = table

    def __parse_system(self) -> None:
        """
        Parse the system information.
        """
        table = self.__normalize_table(self.__root['system'], 'system')
        file_name = _os.path.basename(table['structure'])
        if table['format'] is not None:
            table['format'] = table['format'].lower()
        elif _os.path.splitext(file_name)[0].lower() == 'poscar':
            table['format'] = 'poscar'
        else:
            ext_name = _os.path.splitext(file_name)[1]
            table['format'] = ext_name.lower().strip('.')
            if table['format'] == '':
                message = (
                    'Your file name does contain an extension name, try to '
                    + 'define "format" key in the [system] table.'
                )
                raise RuntimeError(message)
        if table['format'] in ['pdb', 'pqr']:
            self.__system = _mdcore.systems.create_system_from_pdb(
                table['structure']
            )
        elif table['format'] == 'poscar':
            self.__system = _mdcore.systems.create_system_from_poscar(
                table['structure']
            )
        else:
            message = 'Structure file could only be in PDB or POSCAR format!'
            raise RuntimeError(message)
        self.__system._label = self.__root['run']['label']
        self.__root['system'] = table

    def __add_group_velocities(self) -> None:
        """
        Add velocities for each groups.

        Notes
        -----
        Since all user defined groups were placed in the starting of the
        `system.groups` list with the given order, we could safely parse
        velocities of these groups. This means, all 'auto groups' do not
        contain specific initial velocities.
        """
        if self.__root['group'] is None:
            return
        for group, table in zip(self.__system.groups, self.__root['group']):
            table = self.__normalize_table(table, 'group')
            atom_list = self.__parse_atom_list(table['atom_list'])
            if table['initial_temperature'] is not None:
                group.add_velocities_from_temperature(
                    table['initial_temperature']
                )
            if table['initial_velocities'] is not None:
                velocities = table['initial_velocities']
                if len(velocities) != len(atom_list) and len(velocities) != 1:
                    message = (
                        'Number of velocity vectors of the '
                        + '"initial_velocities" key in [[group]] '
                        + 'tables should be equal to the atom '
                        + 'number in the corresponding group!'
                    )
                    raise RuntimeError(message)
                if not all(len(v) == 3 for v in velocities):
                    message = (
                        'Each velocity vector of the '
                        + '"initial_velocities" key in [[group]] '
                        + 'tables should have a length of three!'
                    )
                    raise RuntimeError(message)
                group.velocities += velocities

    def __parse_groups(self) -> None:
        """
        Set up atom groups in the simulated system.
        """
        if self.__root['group'] is None:
            table = {
                'atom_list': range(0, self.__system.n_atoms),
                'has_translations': False,
                'label': 'GROUP0',
            }
            self.__system.groups.create_from_dict(table)
        else:
            for index, table in enumerate(self.__root['group']):
                table = self.__normalize_table(table, 'group')
                table['atom_list'] = self.__parse_atom_list(table['atom_list'])
                if table['has_translations'] is None:
                    table['has_translations'] = True
                if table['label'] is None:
                    table['label'] = 'GROUP{:d}'.format(index)
                self.__system.groups.create_from_dict(table)
        self.__check_groups()

    def __parse_integrator(self) -> None:
        """
        Parse the integrator information.
        """
        table = self.__normalize_table(self.__root['integrator'], 'integrator')
        if table['splitting'] is not None and table['type'] is not None:
            message = (
                'Key "type" or "splitting" can not be defined in '
                + 'an [integrator] table at the same time!'
            )
            raise RuntimeError(message)
        elif table['splitting'] is None and table['type'] is None:
            message = (
                'Key "type" or "splitting" is required in an [integrator] '
                + 'table! And types of them should be "str" and "list".'
            )
            raise RuntimeError(message)
        timestep = table['timestep']
        if table['relaxation_times'] is None:
            relaxation_times = [0.1]
        else:
            if isinstance(table['relaxation_times'], list):
                relaxation_times = table['relaxation_times']
            else:
                relaxation_times = [table['relaxation_times']]
        if table['temperatures'] is None:
            temperatures = [300]
        else:
            if isinstance(table['temperatures'], list):
                temperatures = table['temperatures']
            else:
                temperatures = [table['temperatures']]
        if table['thermo_groups'] is None:
            atom_numbers = [g.n_atoms for g in self.__system.groups]
            thermo_groups = [atom_numbers.index(self.__system.n_atoms)]
            message = 'The thermostat will act on the atom group: "{}".'
            group_name = self.__system.groups[thermo_groups[0]]._label
            _mdutils.warning.warn(message.format(group_name))
        else:
            if isinstance(table['thermo_groups'], list):
                thermo_groups = table['thermo_groups']
            else:
                thermo_groups = [table['thermo_groups']]
        if table['type'] is not None:
            integrator_type = table['type'].lower()
            if integrator_type == 'vv':
                result = _mdcore.integrators.vv_integrator(timestep)
            elif integrator_type == 'cs4':
                result = _mdcore.integrators.cs4_integrator(timestep)
            elif integrator_type == 'nhc':
                result = _mdcore.integrators.nhc_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups
                )
            elif integrator_type == 'gbaoab':
                result = _mdcore.integrators.gbaoab_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups
                )
            elif integrator_type == 'gobabo':
                result = _mdcore.integrators.gobabo_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups
                )
            elif integrator_type == 'baoab':
                if len(self.__system.constraints) != 0:
                    result = _mdcore.integrators.gbaoab_integrator(
                        timestep, temperatures, relaxation_times, thermo_groups
                    )
                    message = (
                        'You are using constraints. Thus the BAOAB integrator '
                        + 'has been changed to the g-BAOAB integrator.'
                    )
                    _mdutils.warning.warn(message)
                else:
                    result = _mdcore.integrators.baoab_integrator(
                        timestep, temperatures, relaxation_times, thermo_groups
                    )
            elif integrator_type == 'obabo':
                if len(self.__system.constraints) != 0:
                    result = _mdcore.integrators.gobabo_integrator(
                        timestep, temperatures, relaxation_times, thermo_groups
                    )
                    message = (
                        'You are using constraints. Thus the OBABO integrator '
                        + 'has been changed to the g-OBABO integrator.'
                    )
                    _mdutils.warning.warn(message)
                else:
                    result = _mdcore.integrators.obabo_integrator(
                        timestep, temperatures, relaxation_times, thermo_groups
                    )
            else:
                message = 'Unknown integrator type: ' + integrator_type
                raise RuntimeError(message)
        else:
            result = _mdcore.integrators.INTEGRATOR(
                timestep,
                table['splitting'],
                temperatures,
                relaxation_times,
                thermo_groups,
            )
        if not result._is_nve:
            if len(temperatures) == 0:
                message = (
                    'You are using thermostats, thus the '
                    + '"temperatures" key must be set. And type of '
                    + 'this key should be float or list(float).'
                )
            if len(relaxation_times) == 0:
                message = (
                    'You are using thermostats, thus the '
                    + '"relaxation_times" key must be set. And type '
                    + 'of this key should be float or list(float).'
                )
            if len(thermo_groups) == 0:
                message = (
                    'You are using thermostats, thus the '
                    + '"thermo_groups" key must be set. And type of '
                    + 'this key should be int or list(int).'
                )
        self.__integrator = result

    def __parse_potential_siesta(
        self, inp: dict, atom_list: list
    ) -> _tp.Callable:
        """
        Parse the SIESTA potential options.

        Parameters
        ----------
        inp : dict
            The dictionary that defines the SIESTA potential.
        atom_list : list(int)
            The atom list.
        """
        return _potentials.SIESTA.generator(
            atom_list,
            self.__system,
            inp['siesta_options'],
            inp['siesta_command'],
            inp['pseudopotential_dir'],
        )

    def __parse_potential_dftd3(
        self, inp: dict, atom_list: list
    ) -> _tp.Callable:
        """
        Parse the DFTD3 potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD3 potential.
        atom_list : list(int)
            The atom list.
        """
        if inp['damping'] is None:
            damping = 'ZeroDamping'
        else:
            damping = inp['damping']
        atom_types = self.__system.atomic_types[atom_list]
        return _potentials.DFTD3.generator(
            atom_list, atom_types, inp['functional'], damping, bool(inp['atm'])
        )

    def __parse_potential_dftd4(
        self, inp: dict, atom_list: list
    ) -> _tp.Callable:
        """
        Parse the DFTD4 potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD4 potential.
        atom_list : list(int)
            The atom list.
        """
        if inp['total_charge'] is None:
            total_charge = 0
        else:
            total_charge = inp['total_charge']
        atom_types = self.__system.atomic_types[atom_list]
        return _potentials.DFTD4.generator(
            atom_list,
            atom_types,
            inp['functional'],
            total_charge,
            bool(inp['atm']),
        )

    def __parse_potential_nep(
        self, inp: dict, atom_list: list
    ) -> _tp.Callable:
        """
        Parse the NEP potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the NEP potential.
        atom_list : list(int)
            The atom list.
        """
        atom_symbols = [self.__system.atomic_symbols[i] for i in atom_list]
        return _potentials.NEP.generator(
            atom_list,
            inp['file_name'],
            atom_symbols,
            bool(inp['use_tabulating']),
        )

    def __parse_potential_mace(
        self, inp: dict, atom_list: list
    ) -> _tp.Callable:
        """
        Parse the MACE potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD3 potential.
        atom_list : list(int)
            The atom list.
        """
        if inp['device'] is None:
            device = 'cpu'
        else:
            device = inp['device']
        if inp['energy_unit'] is None:
            energy_unit = (
                _mdutils.constants.AVOGACONST
                * _mdutils.constants.ELECTCONST
                * 0.001
            )
        else:
            energy_unit = inp['energy_unit']
        if inp['length_unit'] is None:
            length_unit = 0.1
        else:
            length_unit = inp['length_unit']
        if inp['model_dtype'] is None:
            model_dtype = 'float64'
        else:
            model_dtype = inp['model_dtype']
        if inp['virial'] is None:
            calculate_virial = True
        else:
            calculate_virial = inp['virial']
        if inp['charge_cv_expr'] is None:
            charge_cv_expr = None
        else:
            charge_cv_expr = eval(inp['charge_cv_expr'])
        atom_types = self.__system.atomic_types[atom_list]
        return _potentials.MACE.generator(
            atom_list,
            inp['file_name'],
            atom_types,
            device,
            energy_unit,
            length_unit,
            model_dtype,
            calculate_virial,
            charge_cv_expr,
        )

    def __parse_potential_plumed(
        self,
        inp: dict,
        timestep: float,
        temperature: float,
        atom_list: list,
        potential_index: int,
    ) -> _tp.Callable:
        """
        Parse the PLUMED potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the PLUMED potential.
        timestep: float
            The integration timestep.
        temperature: float
            The temperature of the thermostat.
        atom_list : list(int)
            The atom list.
        potential_index : int
            Index of this potential.
        extra_cv_potential_index : int
            Index of the potential calculator for reading extra CV gradients
            from.
        """
        if len(atom_list) != self.__system.n_atoms:
            message = 'The PLUMED potential must act on the whole system!'
            raise RuntimeError(message)
        prefix = self.__root['run']['label'] + '.plumed.pot_{:d}'.format(
            potential_index
        )
        if_restart = bool(self.__root['run']['restart_from'])
        return _potentials.PLUMED.generator(
            atom_list,
            inp['file_name'],
            timestep,
            temperature,
            if_restart,
            prefix,
            extra_cv_potential_index=inp["extra_cv_potential_index"],
        )

    def __parse_potential(
        self, inp: dict, index: int, timestep: float, temperature: float
    ) -> _tp.Callable:
        """
        Parse one potential with given index.

        Parameters
        ----------
        inp : dict
            Information about the potential.
        index : int
            Index of the potential.
        timestep : float
            Timestep of the integrator. In unit of (ps).
        temperatures : list(float)
            Temperatures of the integrator. In unit of (K).
        """
        if inp['atom_list'] is None:
            atom_list = list(range(0, self.__system.n_atoms))
        else:
            atom_list = inp['atom_list']
        potential_type = inp['type'].lower()
        if potential_type == 'siesta':
            generator = self.__parse_potential_siesta(inp, atom_list)
        elif potential_type == 'dftd3':
            generator = self.__parse_potential_dftd3(inp, atom_list)
        elif potential_type == 'dftd4':
            generator = self.__parse_potential_dftd4(inp, atom_list)
        elif potential_type == 'nep':
            generator = self.__parse_potential_nep(inp, atom_list)
        elif potential_type == 'mace':
            generator = self.__parse_potential_mace(inp, atom_list)
        elif potential_type == 'plumed':
            generator = self.__parse_potential_plumed(
                inp, timestep, temperature, atom_list, index
            )
        else:
            message = 'Unknown potential type: "{:s}"!'
            raise RuntimeError(message.format(potential_type))
        return generator

    def __parse_potentials(
        self, timestep: float, temperatures: list, thermo_groups: list
    ) -> None:
        """
        Set up potentials in the simulated system.

        Parameters
        ----------
        timestep : float
            Timestep of the integrator. In unit of (ps).
        temperatures : list(float)
            Temperatures of the integrator. In unit of (K).
        thermo_groups : list(int)
            The thermostated groups.
        """
        if len(temperatures) == 0:
            temperature = None
        elif len(temperatures) == 1:
            temperature = temperatures[0]
        else:
            n_dof = 0
            temperature = 0
            for i in thermo_groups:
                n_dof += self.__system.groups[i].n_dof
                temperature += temperatures[i] * self.__system.groups[i].n_dof
            temperature /= n_dof
        self.__potential_generators = []
        for index, table in enumerate(self.__root['potential']):
            table = self.__normalize_table(table, 'potential')
            if table['file_name'] is not None:
                table['file_name'] = _os.path.abspath(table['file_name'])
            if table['pseudopotential_dir'] is not None:
                table['pseudopotential_dir'] = _os.path.abspath(
                    table['pseudopotential_dir']
                )
            else:
                table['pseudopotential_dir'] = _os.getcwd()
            self.__potential_generators.append((
                table['type'].lower(),
                self.__parse_potential(table, index, timestep, temperature),
            ))

    def __parse_constraints(self) -> None:
        """
        Parse the constraint information.
        """
        table = self.__root['constraints']
        if table is None:
            return
        else:
            table = self.__normalize_table(table, 'constraints')
        n_constrints = len(table['types'])
        if len(table['indices']) != n_constrints:
            message = (
                'The "types" and "indices" keys in the '
                + '[constraints] table should have the same length!'
            )
            raise RuntimeError(message)
        if len(table['targets']) != n_constrints:
            message = (
                'The "types" and "targets" keys in the '
                + '[constraints] table should have the same length!'
            )
            raise RuntimeError(message)
        if not all(isinstance(i, int) for i in table['types']):
            message = (
                'The "types" key in the [constraints] table should '
                + 'corresponding to list of integers!'
            )
            raise RuntimeError(message)
        if not all(isinstance(i, list) for i in table['indices']):
            message = (
                'The "indices" key in the [constraints] table '
                + 'should corresponding to list of lists!'
            )
            raise RuntimeError(message)
        if not all(isinstance(i, float) for i in table['targets']):
            message = (
                'The "targets" key in the [constraints] table '
                + 'should corresponding to list of float numbers!'
            )
            raise RuntimeError(message)
        if table['tolerances'] is None:
            tolerances = [1e-14] * len(table['types'])
        else:
            if len(table['tolerances']) != n_constrints:
                message = (
                    'The "types" and "tolerances" keys in the '
                    + '[constraints] table should have the same length!'
                )
                raise RuntimeError(message)
            else:
                tolerances = table['tolerances']
            if not all(isinstance(i, float) for i in tolerances):
                message = (
                    'The "tolerances" key in the [constraints] '
                    + 'table should corresponding to list of float '
                    + 'numbers!'
                )
                raise RuntimeError(message)
        if table['max_cycles'] is not None:
            self.__system.constraints.max_cycles = table['max_cycles']
        result = []
        for i in range(0, n_constrints):
            result.append({
                'type': table['types'][i],
                'indices': table['indices'][i],
                'target': table['targets'][i],
                'tolerance': tolerances[i],
            })
        self.__system.constraints.appends(result)

    def __parse_barostat(self) -> None:
        """
        Parse the barostat information.
        """
        if self.__root['barostat'] is None:
            self.__barostat = None
            return
        else:
            table = self.__normalize_table(self.__root['barostat'], 'barostat')
        if isinstance(table['beta'], float):
            beta = [table['beta'] / _mdutils.constants.MEGAPASCAL]
        else:
            beta = [b / _mdutils.constants.MEGAPASCAL for b in table['beta']]
        if isinstance(table['pressures'], float):
            pressures = [table['pressures'] * _mdutils.constants.MEGAPASCAL]
        else:
            pressures = [
                p * _mdutils.constants.MEGAPASCAL for p in table['pressures']
            ]
        self.__barostat = _mdapps.barostats.BAROSTAT(
            pressures, beta, table['relaxation_time']
        )

    def __parse_loggers(self) -> None:
        """
        Parse the logger information.
        """
        self.__loggers = []
        if self.__root['logger'] is None:
            self.__loggers.append(
                _mdapps.loggers.DEFAULTCSVLOGGER(
                    self.__root['run']['label'] + '.data.csv',
                    interval=1,
                    append=bool(self.__root['run']['restart_from']),
                )
            )
        else:
            for index, table in enumerate(self.__root['logger']):
                table = self.__normalize_table(table, 'logger')
                if table['format'] is None:
                    logger_format = 'csv'
                else:
                    logger_format = table['format'].lower()
                if logger_format == 'csv':
                    delimiter = ' , '
                elif logger_format == 'txt':
                    delimiter = ' '
                else:
                    message = 'Unknown logger format "{}"'
                    raise RuntimeError(message.format(table['format']))
                if table['prefix'] is None:
                    prefix = self.__root['run']['label']
                else:
                    prefix = table['prefix']
                if table['interval'] is None:
                    interval = 1
                else:
                    interval = table['interval']
                file_name = prefix + '.data.' + logger_format
                self.__loggers.append(
                    _mdapps.loggers.DEFAULTCSVLOGGER(
                        file_name,
                        interval=interval,
                        delimiter=delimiter,
                        append=bool(self.__root['run']['restart_from']),
                        potential_list=table['potential_list'],
                    )
                )
        if self.__root['active_learning'] is not None:
            self.__loggers = []
            message = (
                'You are performing active learning, system data '
                + 'loggers will be removed.'
            )
            _mdutils.warning.warn(message)

    def __check_trajectories(self) -> None:
        """
        Check trajectory settings.
        """
        has_restart_flag = False
        for trajectory in self.__trajectories:
            if hasattr(trajectory, 'is_restart') and trajectory.is_restart:
                has_restart_flag = True
                break
        if not has_restart_flag:
            self.__trajectories.append(
                _mdapps.trajectories.H5WRITER(
                    self.__root['run']['label'] + '.restart.h5',
                    interval=10,
                    restart_file=True,
                )
            )
        if self.__root['active_learning'] is not None:
            self.__trajectories = []
            message = (
                'You are performing active learning, trajectory '
                + 'writers will be removed.'
            )
            _mdutils.warning.warn(message)

    def __parse_trajectories(self) -> None:
        """
        Parse the trajectory information.
        """
        self.__trajectories = []
        if self.__root['trajectory'] is None:
            self.__trajectories.append(
                _mdapps.trajectories.H5WRITER(
                    self.__root['run']['label'] + '.trajectory.h5',
                    interval=1,
                    write_virial=True,
                    write_forces=True,
                    append=bool(self.__root['run']['restart_from']),
                )
            )
        else:
            for index, table in enumerate(self.__root['trajectory']):
                table = self.__normalize_table(table, 'trajectory')
                if table['format'] is None:
                    trajectory_format = 'h5'
                else:
                    trajectory_format = table['format'].lower()
                if table['prefix'] is None:
                    prefix = self.__root['run']['label']
                else:
                    prefix = table['prefix']
                if table['interval'] is None:
                    interval = 1
                else:
                    interval = table['interval']
                if trajectory_format == 'h5':
                    if table['is_restart_file']:
                        file_name = prefix + '.restart.h5'
                    else:
                        file_name = prefix + '.trajectory.h5'
                    writer = _mdapps.trajectories.H5WRITER(
                        file_name,
                        interval=interval,
                        write_virial=True,
                        write_velocities=bool(table['write_velocities']),
                        write_forces=bool(table['write_forces']),
                        wrap_positions=bool(table['wrap_positions']),
                        append=bool(self.__root['run']['restart_from']),
                        restart_file=bool(table['is_restart_file']),
                        use_double=bool(table['use_float64']),
                        potential_list=table['potential_list'],
                    )
                    self.__trajectories.append(writer)
                elif trajectory_format == 'exyz':
                    file_name = prefix + '.trajectory.xyz'
                    writer = _mdapps.trajectories.EXYZWRITER(
                        file_name,
                        interval=interval,
                        write_virial=True,
                        write_velocities=bool(table['write_velocities']),
                        write_forces=bool(table['write_forces']),
                        wrap_positions=bool(table['wrap_positions']),
                        append=bool(self.__root['run']['restart_from']),
                        potential_list=table['potential_list'],
                        energy_shift=table['energy_shift'],
                    )
                    self.__trajectories.append(writer)
                else:
                    message = 'Unknown trajectory format: "{}"'
                    raise RuntimeError(message.format(table['format']))
                if table['potential_list'] is None:
                    tmp = list(range(0, len(self.__potential_generators)))
                else:
                    tmp = table['potential_list']
                for i in tmp:
                    potential_name = self.__potential_generators[i][0]
                    if potential_name == 'plumed':
                        message = (
                            'The forces and energies in the trajectory file '
                            + '"{:s}" may include contributions of PLUMED '
                            + 'bias potentials! MAKE SURE THAT THIS IS WHAT '
                            + 'YOU WANT!!'
                        )
                        _mdutils.warning.warn(message.format(file_name))
        self.__check_trajectories()

    def __parse_scripts(self) -> None:
        """
        Set up the post-step scripts.
        """
        scope = {}
        self.__scripts = []
        if self.__root['script'] is None:
            return
        for index, table in enumerate(self.__root['script']):
            table = self.__normalize_table(table, 'script')
            if table['interval'] is None:
                interval = 1
            else:
                interval = table['interval']
            if table['initialize'] is not None:
                exec(table['initialize'], scope)
            else:
                scope['initialize'] = lambda: None
            exec(table['update'], scope)
            if 'update' not in scope.keys():
                message = (
                    'The function name defined in the "update" key '
                    + 'of the [[script]] table(s) must be "update"!'
                )
                raise RuntimeError(message)
            obj = _mdapps.utils.post_step.POSTSTEPOBJWRAPPER(
                scope['update'], scope['initialize'], interval
            )
            self.__scripts.append(obj)

    def __set_up_simulation(self) -> None:
        """
        Set up the simulation protocol.
        """
        for generator in self.__potential_generators:
            self.__system.potentials.append(generator[1]())
        self.__simulation = _mdapps.simulations.SIMULATION(
            self.__system,
            self.__integrator,
            barostat=self.__barostat,
            loggers=self.__loggers,
            trajectories=self.__trajectories,
        )
        for obj in self.__scripts:
            self.__simulation.post_step_objects.append(obj)
        restart_file = self.__root['run']['restart_from']
        if restart_file is not None:
            self.__simulation.restart_from(restart_file)

    def __parse_active_learning(self) -> None:
        """
        Parse the active learning information.
        """
        table = self.__root['active_learning']
        table = self.__normalize_table(table, 'active_learning')
        excluded_names = ['plumed', 'nep', 'mace']
        if table['reference_potentials'] is None:
            reference_potentials = []
            for i in range(0, len(self.__potential_generators)):
                potential_name = self.__potential_generators[i][0]
                if potential_name not in excluded_names:
                    reference_potentials.append(i)
        else:
            reference_potentials = table['reference_potentials']
        for i in reference_potentials:
            if i >= len(self.__potential_generators):
                message = (
                    'Unknown potential index: {:d} in the '
                    + '[active_learning] table!'
                )
                raise IndexError(message.format(i))
            potential_name = self.__potential_generators[i][0]
            if potential_name in excluded_names:
                message = (
                    'You are using "{:s}" as one of the reference '
                    + 'potentials! MAKE SURE THAT THIS IS WHAT YOU WANT!!'
                )
                _mdutils.warning.warn(message.format(potential_name.upper()))
        if table['initial_training_set'] is not None:
            table['initial_training_set'] = _os.path.abspath(
                table['initial_training_set']
            )
        if table['initial_testing_set'] is not None:
            table['initial_testing_set'] = _os.path.abspath(
                table['initial_testing_set']
            )
        if table['initial_potential_files'] is not None:
            table['initial_potential_files'] = [
                _os.path.abspath(file)
                for file in table['initial_potential_files']
            ]
        post_step_objects = [*self.__scripts]
        if self.__barostat is not None:
            post_step_objects.insert(0, self.__barostat)
        generators = [g[1] for g in self.__potential_generators]
        self.__simulation = _mdapps.active_learning.ACTIVELEARNING(
            self.__system,
            self.__integrator,
            generators,
            reference_potentials,
            {k: v for k, v in table.items() if v is not None},
            table['nep_options'],
            table['nep_command'],
            bool(table['use_tabulating']),
            post_step_objects,
            self.__root['run']['label'] + '.active_learning',
        )

    def run(self) -> None:
        """
        Run the simulation.
        """
        if self.__root['active_learning'] is not None:
            n_steps = self.__root['active_learning']['n_iterations']
        else:
            n_steps = self.__root['run']['n_steps']
        self.__simulation.run(n_steps)

    @property
    def file_name(self) -> str:
        """
        Name of the configure file.
        """
        return self.__file_name

    @property
    def simulation(self) -> _tp.Union:
        """
        The simulation protocol.
        """
        return self.__simulation
