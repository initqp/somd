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
import warnings as _w
from typing import get_args as _get_args
from collections import namedtuple as _namedtuple
from somd import core as _mdcore
from somd import apps as _mdapps
from somd import potentials as _potentials
from somd.constants import CONSTANTS as _c
from somd.constants import SOMDDEFAULTS as _d

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
        'n_steps': __value__(int, True, None),
        'seed': __value__(int, False, None),
        'label': __value__(str, False, None),
        'restart_from': __value__(str, False, None)
    }
    __parameters__['system'] = {
        'structure': __value__(str, True, None),
        'format': __value__(str, False, None),
        'box': __value__(list, False, None)
    }
    __parameters__['integrator'] = {
        'timestep': __value__(float, True, None),
        'type': __value__(str, False, None),
        'splitting': __value__(list, False, None),
        'temperatures': __value__(_tp.Union[list, float], False, None),
        'relaxation_times': __value__(_tp.Union[list, float], False, None),
        'thermo_groups': __value__(_tp.Union[list, int], False, None)
    }
    __parameters__['potential'] = {
        'type': __value__(str, True, None),
        'atom_list': __value__(_tp.Union[list, str], False, None),
        'siesta_options': __value__(str, True, __dep__('type', ['siesta'])),
        'siesta_command': __value__(str, True, __dep__('type', ['siesta'])),
        'pseudopotential_dir': __value__(str, False,
                                         __dep__('type', ['siesta'])),
        'functional': __value__(str, True,
                                __dep__('type', ['dftd3', 'dftd4'])),
        'damping': __value__(str, False, __dep__('type', ['dftd3'])),
        'atm': __value__(bool, False, __dep__('type', ['dftd3', 'dftd4'])),
        'total_charge': __value__(int, False, __dep__('type', ['dftd4'])),
        'file_name': __value__(str, True, __dep__('type', ['plumed', 'nep'])),
        'use_tabulating': __value__(bool, False, __dep__('type', ['nep']))
    }
    __parameters__['group'] = {
        'atom_list': __value__(_tp.Union[list, str], True, None),
        'label': __value__(str, False, None),
        'has_translations': __value__(bool, False, None),
        'initial_velocities': __value__(list, False, None),
        'initial_temperature': __value__(float, False, None)
    }
    __parameters__['barostat'] = {
        'pressures': __value__(_tp.Union[list, float], True, None),
        'beta': __value__(_tp.Union[list, float], True, None),
        'relaxation_time': __value__(float, True, None),
    }
    __parameters__['constraints'] = {
        'types': __value__(list, True, None),
        'indices': __value__(list, True, None),
        'targets': __value__(list, True, None),
        'tolerances': __value__(list, False, None),
        'max_cycles': __value__(int, False, None)
    }
    __parameters__['trajectory'] = {
        'prefix': __value__(str, False, None),
        'format': __value__(str, False, None),
        'interval': __value__(int, False, None),
        'write_forces': __value__(bool, False, None),
        'write_velocities': __value__(bool, False, None),
        'wrap_positions': __value__(bool, False, None),
        'potential_list': __value__(list, False, None),
        'use_float64': __value__(bool, False, __dep__('format', ['h5'])),
        'energy_shift': __value__(float, False, __dep__('format', ['exyz'])),
        'is_restart_file': __value__(bool, False, __dep__('format', ['h5']))
    }
    __parameters__['logger'] = {
        'format': __value__(str, False, None),
        'prefix': __value__(str, False, None),
        'interval': __value__(int, False, None),
        'potential_list': __value__(list, False, None)
    }
    __parameters__['script'] = {
        'update': __value__(str, True, None),
        'interval': __value__(int, False, None),
        'initialize': __value__(str, False, None)
    }
    __parameters__['active_learning'] = {
        'nep_options': __value__(str, True, None),
        'nep_command': __value__(str, True, None),
        'initial_training_set': __value__(str, True, None),
        'n_iterations': __value__(int, True, None),
        'n_potentials': __value__(int, False, None),
        'max_md_runs_per_iter': __value__(int, False, None),
        'max_md_steps_per_iter': __value__(int, False, None),
        'msd_lower_limit': __value__(float, False, None),
        'msd_upper_limit': __value__(float, False, None),
        'min_new_structures_per_iter': __value__(int, False, None),
        'max_new_structures_per_iter': __value__(int, False, None),
        'initial_potential_files': __value__(list, False, None),
        'initial_testing_set': __value__(str, False, None),
        'reference_potentials': __value__(list, False, None),
        'use_tabulating': __value__(bool, False, None),
        'energy_shift': __value__(float, False, None),
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
        if (self.__integrator._is_nve):
            self.__parse_potentials(self.__integrator.timestep, [], [])
        else:
            self.__parse_potentials(self.__integrator.timestep,
                                    self.__integrator.temperatures,
                                    self.__integrator.thermo_groups)
        self.__parse_barostat()
        self.__parse_loggers()
        self.__parse_trajectories()
        self.__parse_scripts()
        if (self.__root['active_learning'] is None):
            self.__set_up_simulation()
            self.__trainer = None
        else:
            self.__parse_active_learning()

    def __normalize_tables(self) -> None:
        """
        Check root tables.
        """
        definitions = dict.fromkeys(self.__tables__.keys(), None)
        for key in self.__root.keys():
            if key in self.__tables__.keys():
                table = self.__root[key]
                if (self.__tables__[key].is_array and type(table) != list):
                    message = 'The "{}" key should correspond to an ' + \
                              'array of tables!'
                    raise TypeError(message.format(key))
                if (not self.__tables__[key].is_array and type(table) != dict):
                    message = 'The "{}" key should correspond to a table!'
                    raise TypeError(message.format(key))
                definitions[key] = table
            else:
                raise KeyError('Unknown root key or table "{}"'.format(key))
        for key in self.__tables__.keys():
            if (self.__tables__[key].required and (definitions[key] is None)):
                raise KeyError('Table [{}] is required!'.format(key))
        self.__root = definitions

    def __check_task(self):
        """
        Check the simulation task.
        """
        task_tables = ['run', 'active_learning']
        task_table_list = [self.__root[i] for i in task_tables]
        if (all((i is None) for i in task_table_list)):
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
        parameters = self.__parameters__[table_name]
        definitions = dict.fromkeys(parameters.keys(), None)
        # Check key names and value types.
        for key in inp.keys():
            if key in parameters.keys():
                value = inp[key]
                if (not issubclass(type(value), parameters[key].type)):
                    if (parameters[key].type.__name__ == 'Union'):
                        type_list = _get_args(parameters[key].type)
                        types = '"' + type_list[0].__name__ + '"'
                        for t in type_list[1:]:
                            types += ' or "' + t.__name__ + '"'
                    else:
                        types = '"' + parameters[key].type.__name__ + '"'
                    message = 'Wrong type of key "{}" in (one of) the ' + \
                              '[{}] table(s)! A {} typed value is expected!'
                    message = message.format(key, table_name, types)
                    raise TypeError(message)
                definitions[key] = value
            else:
                message = 'Unknown key "{}" in (one of) the [{}] table(s)!'
                raise KeyError(message.format(key, table_name))
        # Check dependencies.
        for key in definitions.keys():
            dependency = parameters[key].dependency
            if (dependency is not None and definitions[key] is not None):
                value = definitions[dependency.key]
                if (type(value) == str):
                    value = value.lower()
                if (value not in dependency.values):
                    values = '"' + str(dependency.values[0]) + '"'
                    for v in dependency.values[1:]:
                        values += ' or "' + str(v) + '"'
                    message = 'Key "{}" in the [{}] table(s) is only ' + \
                              'required when the value of its "{}" key ' + \
                              'is {}!'
                    message = message.format(key, table_name, dependency.key,
                                             values)
                    raise KeyError(message)
        # Check required keys.
        for key in parameters.keys():
            if (parameters[key].required and (definitions[key] is None)):
                dependency = parameters[key].dependency
                if (dependency is None):
                    message = 'Key "{}" is required in the [{}] table(s)!'
                    raise KeyError(message.format(key, table_name))
                else:
                    value = definitions[dependency.key]
                    if (type(value) == str):
                        value = value.lower()
                    if (value in dependency.values):
                        values = '"' + str(dependency.values[0]) + '"'
                        for v in dependency.values[1:]:
                            values += ' or "' + str(v) + '"'
                        message = 'Key "{}" in the [{}] table(s) is ' + \
                                  'required when the value of its "{}" ' + \
                                  'key is {}!'
                        message = message.format(key, table_name,
                                                 dependency.key, values)
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
        if (type(inp) == list):
            if (not all(isinstance(i, int) for i in inp)):
                raise RuntimeError('Unknown atom indices: {}'.format(inp))
            result = inp
        elif (type(inp) == str):
            if (inp.lower() == 'all'):
                result = list(range(0, self.__system.n_atoms))
            else:
                result = []
                for s in inp.split(','):
                    if (':' not in s):
                        if (not s.isnumeric()):
                            raise SyntaxError('Unknown atom index: ' + s)
                        result.append(int(s))
                    else:
                        l = s.split(':')
                        if (not l[0].isnumeric() or len(l) != 2 or
                                not l[1].isnumeric()):
                            raise SyntaxError('Unknown atom range: ' + s)
                        for i in range(int(l[0]), int(l[1]) + 1):
                            result.append(i)
        else:
            message = 'Type of key "atom_list" could only be list(int) or str!'
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
            if (group.n_atoms == self.__system.n_atoms):
                whole_system_flag = True
            if (group.has_translations is False):
                no_translations_flag = True
        # No group is corresponding to the whole group. Add a new group.
        if (whole_system_flag is False):
            label = 'GROUP{:d}'.format(len(self.__system.groups))
            d = {'atom_list': range(0, self.__system.n_atoms),
                 'has_translations': True, 'label': label}
            self.__system.groups.create_from_dict(d)
            message = 'An atom group that corresponding to the whole ' + \
                      'system has been append the group list.'
            self.__root['group'].append(d)
            _w.warn(message)
        # The user has not defined the group without translations.
        # Check if there is any group that is corresponding to the whole group.
        if (no_translations_flag is False):
            for index, group in enumerate(self.__system.groups):
                if (group.n_atoms == self.__system.n_atoms):
                    group.has_translations = False
                    no_translations_flag = True
                    self.__root['group'][index]['has_translations'] = False
                    message = 'Translational degrees of freedom of group ' + \
                              '"{}" has been automatically removed.'
                    _w.warn(message.format(group._label))
                    break

    def __parse_run(self) -> None:
        """
        Parse the simulation task information.
        """
        run = self.__root['run']
        if (run is None):
            run = dict()
            label = _os.path.splitext(self.__file_name)[0]
            run['label'] = _os.path.basename(label)
            run['restart_from'] = None
        else:
            run = self.__normalize_table(self.__root['run'], 'run')
            self.__n_steps = run['n_steps']
            if (run['seed'] is not None):
                _np.random.seed(run['seed'])
            if (run['label'] is None):
                label = _os.path.splitext(self.__file_name)[0]
                run['label'] = _os.path.basename(label)
        self.__root['run'] = run

    def __parse_system(self) -> None:
        """
        Parse the system information.
        """
        system = self.__normalize_table(self.__root['system'], 'system')
        file_name = system['structure']
        if (system['format'] is not None):
            ext_name = system['format'].lower()
        else:
            ext_name = _os.path.splitext(file_name)[1].lower()[1:]
            if (ext_name == ''):
                if (_os.path.splitext(file_name)[0].lower() == "poscar"):
                    ext_name = 'poscar'
                else:
                    message = 'Your file name does contain an extension ' + \
                              ' name, try to define "format" key in the ' + \
                              '[system] table.'
                    raise RuntimeError(message)
            system['format'] = ext_name
        if (ext_name == 'pdb'):
            self.__system = \
                _mdcore.systems.create_system_from_pdb(file_name)
        elif (ext_name == 'poscar'):
            self.__system = \
                _mdcore.systems.create_system_from_poscar(file_name)
        else:
            message = 'Structure file could only be in PDB ' + \
                      'or POSCAR format!'
            raise RuntimeError(message)
        if (system['box'] is not None):
            if ((not all(isinstance(v, list) for v in system['box'])) or
                    (not all((len(v) == 3) for v in system['box']))):
                message = 'The "box" key in the "[group]" table requires ' + \
                          'three box vectors with a length of three!'
                raise RuntimeError(message)
            else:
                self.__system.box[:] = system['box'][:]
        parameter_names = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        for i in range(0, 6):
            if (self.__system.lattice[i] < _d.LATTICETOL):
                message = 'Very small lattice parameter: {} !' + \
                          'If your structure file does not contain cell ' + \
                          'data, you could set the simulation box with ' + \
                          '"box" key in the [system] table.'
                raise RuntimeError(message.format(parameter_names[i]))
        self.__system._label = self.__root['run']['label']
        self.__root['system'] = system

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
        if (self.__root['group'] is None):
            return
        groups = self.__root['group']
        for index, group in enumerate(groups):
            group = self.__normalize_table(group, 'group')
            atom_list = self.__parse_atom_list(group['atom_list'])
            if (group['initial_temperature'] is not None):
                self.__system.groups[index].add_velocities_from_temperature(
                    group['initial_temperature'])
            if (group['initial_velocities'] is not None):
                velocities = group['initial_velocities']
                if (len(velocities) != len(atom_list) and
                        len(velocities) != 1):
                    message = 'Number of velocity vectors of the ' + \
                              '"initial_velocities" key in [[group]] ' + \
                              'tables should be equal to the atom ' + \
                              'number in the corresponding group!'
                    raise RuntimeError(message)
                if (not all(len(v) == 3 for v in velocities)):
                    message = 'Each velocity vector of the ' + \
                              '"initial_velocities" key in [[group]] ' + \
                              'tables should have a length of three!'
                    raise RuntimeError(message)
                self.__system.groups[index].velocities += velocities

    def __parse_groups(self) -> None:
        """
        Set up atom groups in the simulated system.
        """
        if (self.__root['group'] is None):
            group = {'atom_list': range(0, self.__system.n_atoms),
                     'has_translations': False, 'label': 'GROUP0'}
            self.__system.groups.create_from_dict(group)
            group['initial_velocities'] = None
            group['initial_temperature'] = None
        else:
            groups = self.__root['group']
            if (not all(isinstance(group, dict) for group in groups)):
                message = 'The "group" key should correspond to an array ' + \
                          'of tables!'
                raise RuntimeError(message)
            for index, group in enumerate(groups):
                group = self.__normalize_table(group, 'group')
                group['atom_list'] = self.__parse_atom_list(group['atom_list'])
                if (group['has_translations'] is None):
                    group['has_translations'] = True
                if (group['label'] is None):
                    group['label'] = 'GROUP{:d}'.format(index)
                self.__system.groups.create_from_dict(group)
            self.__root['group'] = groups
        self.__check_groups()

    def __parse_integrator(self) -> None:
        """
        Parse the integrator information.
        """
        integrator = self.__normalize_table(self.__root['integrator'],
                                            'integrator')
        if (integrator['splitting'] is not None and
                integrator['type'] is not None):
            message = 'Key "type" or "splitting" can not be defined in ' + \
                      'an [integrator] table at the same time!'
            raise RuntimeError(message)
        elif (integrator['splitting'] is None and integrator['type'] is None):
            message = 'Key "type" or "splitting" is required in an ' + \
                      '[integrator] table! And the types of these ' + \
                      'keys should be "str" and "list".'
            raise RuntimeError(message)
        timestep = integrator['timestep']
        if (integrator['relaxation_times'] is None):
            relaxation_times = []
        else:
            if (type(integrator['relaxation_times']) == list):
                relaxation_times = integrator['relaxation_times']
            else:
                relaxation_times = [integrator['relaxation_times']]
        if (integrator['temperatures'] is None):
            temperatures = []
        else:
            if (type(integrator['temperatures']) == list):
                temperatures = integrator['temperatures']
            else:
                temperatures = [integrator['temperatures']]
        if (integrator['thermo_groups'] is None):
            for index, group in enumerate(self.__system.groups):
                if (group.n_atoms == self.__system.n_atoms):
                    thermo_groups = [index]
                    break
        else:
            if (type(integrator['thermo_groups']) == list):
                thermo_groups = integrator['thermo_groups']
            else:
                thermo_groups = [integrator['thermo_groups']]
        if (integrator['type'] is not None):
            integrator_type = integrator['type'].lower()
            if (integrator_type == 'vv'):
                result = _mdcore.integrators.vv_integrator(timestep)
            elif (integrator_type == 'cs4'):
                result = _mdcore.integrators.vs4_integrator(timestep)
            elif (integrator_type == 'nhc'):
                result = _mdcore.integrators.nhc_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups)
            elif (integrator_type == 'gbaoab'):
                result = _mdcore.integrators.gbaoab_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups)
            elif (integrator_type == 'gobabo'):
                result = _mdcore.integrators.gobabo_integrator(
                    timestep, temperatures, relaxation_times, thermo_groups)
            elif (integrator_type == 'baoab'):
                if (len(self.__system.constraints) != 0):
                    result = _mdcore.integrators.gbaoab_integrator(
                        timestep, temperatures, relaxation_times,
                        thermo_groups)
                    message = 'You are using constraints. Thus the baoab ' + \
                              'integrator has been changed to the gbaoab ' + \
                              'integrator.'
                    _w.warn(message)
                else:
                    result = _mdcore.integrators.baoab_integrator(
                        timestep, temperatures, relaxation_times,
                        thermo_groups)
            elif (integrator_type == 'obabo'):
                if (len(self.__system.constraints) != 0):
                    result = _mdcore.integrators.gobabo_integrator(
                        timestep, temperatures, relaxation_times,
                        thermo_groups)
                    message = 'You are using constraints. Thus the obabo ' + \
                              'integrator has been changed to the gobabo ' + \
                              'integrator.'
                    _w.warn(message)
                else:
                    result = _mdcore.integrators.obabo_integrator(
                        timestep, temperatures, relaxation_times,
                        thermo_groups)
            else:
                message = 'Unknown integrator type: ' + integrator_type
                raise RuntimeError(message)
        else:
            result = _mdcore.integrators.INTEGRATOR(
                timestep, integrator['splitting'], temperatures,
                relaxation_times, thermo_groups)
        if (not result._is_nve):
            if (len(temperatures) == 0):
                message = 'You are using thermostats, thus the ' + \
                          '"temperatures" key must be set. And type of ' + \
                          'this key should be float or list(float).'
            if (len(relaxation_times) == 0):
                message = 'You are using thermostats, thus the ' + \
                          '"relaxation_times" key must be set. And type ' + \
                          'of this key should be float or list(float).'
            if (len(thermo_groups) == 0):
                message = 'You are using thermostats, thus the ' + \
                          '"thermo_groups" key must be set. And type of ' + \
                          'this key should be int or list(int).'
        self.__integrator = result

    def __parse_potential_siesta(self, inp: dict, atom_list: list) \
            -> _potentials.SIESTA:
        """
        Parse the SIESTA potential options.

        Parameters
        ----------
        inp : dict
            The dictionary that defines the SIESTA potential.
        atom_list : list(int)
            The atom list.
        """
        return _potentials.create_siesta_potential(
            self.__system, atom_list, inp['siesta_options'],
            inp['siesta_command'], inp['pseudopotential_dir'])

    def __parse_potential_dftd3(self, inp: dict, atom_list: list) \
            -> _potentials.DFTD3:
        """
        Parse the DFTD3 potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD3 potential.
        atom_list : list(int)
            The atom list.
        """
        if (inp['damping'] is None):
            damping = 'ZeroDamping'
        else:
            damping = inp['damping']
        atom_types = self.__system.atomic_types[atom_list]
        return _potentials.DFTD3(atom_list, atom_types, inp['functional'],
                                 damping, bool(inp['atm']))

    def __parse_potential_dftd4(self, inp: dict, atom_list: list) \
            -> _potentials.DFTD4:
        """
        Parse the DFTD4 potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD3 potential.
        atom_list : list(int)
            The atom list.
        """
        if (inp['total_charge'] is None):
            total_charge = 0
        else:
            total_charge = inp['total_charge']
        atom_types = self.__system.atomic_types[atom_list]
        return _potentials.DFTD4(atom_list, atom_types, inp['functional'],
                                 total_charge, bool(inp['atm']))

    def __parse_potential_nep(self, inp: dict, atom_list: list) \
            -> _potentials.NEP:
        """
        Parse the NEP potential options.

        Parameters
        ----------
        inp: dict
            The dictionary that defines the DFTD3 potential.
        atom_list : list(int)
            The atom list.
        """
        atom_symbols = [self.__system.atomic_symbols[i] for i in atom_list]
        return _potentials.NEP(atom_list, inp['file_name'], atom_symbols,
                               bool(inp['use_tabulating']))

    def __parse_potential_plumed(self,
                                 inp: dict,
                                 timestep: float,
                                 temperature: float,
                                 atom_list: list,
                                 potential_index: int) -> _potentials.PLUMED:
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
        """
        if (len(atom_list) != self.__system.n_atoms):
            message = 'The PLUMED potential must act on the whole system!'
            raise RuntimeError(message)
        prefix = self.__root['run']['label'] + \
            '.plumed.pot_{:d}'.format(potential_index)
        return _potentials.PLUMED(atom_list, inp['file_name'], timestep,
                                  temperature,
                                  bool(self.__root['run']['restart_from']),
                                  prefix)

    def __parse_potential(self,
                          inp: dict,
                          index: int,
                          timestep: float,
                          temperature: float) \
            -> _mdcore.potential_base.POTENTIAL:
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
        if (inp['atom_list'] is None):
            atom_list = list(range(0, self.__system.n_atoms))
        else:
            atom_list = inp['atom_list']
        potential_type = inp['type'].upper()
        if (potential_type == 'SIESTA'):
            potential = self.__parse_potential_siesta(inp, atom_list)
        elif (potential_type == 'DFTD3'):
            potential = self.__parse_potential_dftd3(inp, atom_list)
        elif (potential_type == 'DFTD4'):
            potential = self.__parse_potential_dftd4(inp, atom_list)
        elif (potential_type == 'NEP'):
            potential = self.__parse_potential_nep(inp, atom_list)
        elif (potential_type == 'PLUMED'):
            potential = self.__parse_potential_plumed(
                inp, timestep, temperature, atom_list, index)
        else:
            message = 'Unknown potential type: ' + potential_type
            raise RuntimeError(message)
        return potential

    def __parse_potentials(self,
                           timestep: float,
                           temperatures: list,
                           thermo_groups: list) -> None:
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
        if (len(temperatures) == 0):
            temperature = None
        elif (len(temperatures) == 1):
            temperature = temperatures[0]
        else:
            n_dof = 0
            temperature = 0
            for i in thermo_groups:
                n_dof += self.__system.groups[i].n_dof
                temperature += temperatures[i] * self.__system.groups[i].n_dof
            temperature /= n_dof
        potentials = self.__root['potential']
        if (not all(isinstance(potential, dict) for potential in potentials)):
            message = 'The "potential" key should correspond to an array ' + \
                      'of tables!'
            raise RuntimeError(message)
        self.__potential_generators = []
        for index, potential in enumerate(potentials):
            potential = self.__normalize_table(potential, 'potential')
            if (potential['file_name'] is not None):
                potential['file_name'] = \
                    _os.path.abspath(potential['file_name'])
            if (potential['pseudopotential_dir'] is not None):
                potential['pseudopotential_dir'] = \
                    _os.path.abspath(potential['pseudopotential_dir'])
            else:
                potential['pseudopotential_dir'] = _os.getcwd()
            self.__potential_generators.append((
                potential['type'].upper(), lambda i=index, p=potential:
                self.__parse_potential(p, i, timestep, temperature)))

    def __parse_constraints(self):
        """
        Parse the constraint information.
        """
        constraints = self.__root['constraints']
        if (constraints is None):
            return
        else:
            constraints = self.__normalize_table(constraints, 'constraints')
        n_constrints = len(constraints['types'])
        if ((len(constraints['indices']) != n_constrints) or
                (len(constraints['types']) != n_constrints)):
            message = 'The "types", "indices" and "targets" keys in the ' + \
                      '[constraints] table should have the same length!'
            raise RuntimeError(message)
        if (not all(isinstance(i, int) for i in constraints['types'])):
            message = 'The "types" key in the [constraints] table should ' + \
                      'corresponding to list of integers!'
            raise RuntimeError(message)
        if (not all(isinstance(i, list) for i in constraints['indices'])):
            message = 'The "indices" key in the [constraints] table ' + \
                      'should corresponding to list of lists!'
            raise RuntimeError(message)
        if (not all(isinstance(i, float) for i in constraints['targets'])):
            message = 'The "targets" key in the [constraints] table ' + \
                      'should corresponding to list of float numbers!'
            raise RuntimeError(message)
        if (constraints['tolerances'] is None):
            tolerances = [1E-14] * len(constraints['types'])
        else:
            if (len(constraints['tolerances']) != n_constrints):
                message = 'The "types" and "tolerances" keys in the ' + \
                          '[constraints] table should have the same length!'
                raise RuntimeError(message)
            else:
                tolerances = constraints['tolerances']
            if (not all(isinstance(i, float) for i in tolerances)):
                message = 'The "tolerances" key in the [constraints] ' + \
                          'table should corresponding to list of float ' + \
                          'numbers!'
                raise RuntimeError(message)
        if (constraints['max_cycles'] is not None):
            self.__system.constraints.max_cycles = constraints['max_cycles']
        result = []
        for i in range(0, n_constrints):
            result.append({'type': constraints['types'][i],
                           'indices': constraints['indices'][i],
                           'target': constraints['targets'][i],
                           'tolerance': tolerances[i]})
        self.__system.constraints.appends(result)

    def __parse_barostat(self) -> None:
        """
        Parse the barostat information.
        """
        barostat = self.__root['barostat']
        if (barostat is None):
            self.__barostat = None
            return
        else:
            barostat = self.__normalize_table(barostat, 'barostat')
        if (type(barostat['beta']) == float):
            beta = [barostat['beta'] / _c.MEGAPASCAL]
        else:
            beta = [b / _c.MEGAPASCAL for b in barostat['beta']]
        if (type(barostat['pressures']) == float):
            pressures = [barostat['pressures'] * _c.MEGAPASCAL]
        else:
            pressures = [p * _c.MEGAPASCAL for p in barostat['pressures']]
        self.__barostat = _mdapps.barostats.BAROSTAT(
            pressures, beta, barostat['relaxation_time'])

    def __parse_loggers(self) -> None:
        """
        Parse the logger information.
        """
        self.__loggers = []
        loggers = self.__root['logger']
        if (loggers is None):
            self.__loggers.append(
                _mdapps.loggers.DEFAULTCSVLOGGER(
                    self.__root['run']['label'] + '.data.csv', interval=1,
                    append=bool(self.__root['run']['restart_from'])))
        else:
            for index, logger in enumerate(loggers):
                logger = self.__normalize_table(logger, 'logger')
                if (logger['format'] is None):
                    logger_format = 'csv'
                else:
                    logger_format = logger['format'].lower()
                if (logger_format == 'csv'):
                    delimiter = ' , '
                elif (logger_format == 'txt'):
                    delimiter = ' '
                else:
                    message = 'Unknown logger format ' + logger['format']
                    raise RuntimeError(message)
                if (logger['prefix'] is None):
                    prefix = self.__root['run']['label']
                else:
                    prefix = logger['prefix']
                if (logger['interval'] is None):
                    interval = 1
                else:
                    interval = logger['interval']
                file_name = prefix + '.data.' + logger_format
                self.__loggers.append(
                    _mdapps.loggers.DEFAULTCSVLOGGER(
                        file_name, interval=interval, delimiter=delimiter,
                        append=bool(self.__root['run']['restart_from']),
                        potential_list=logger['potential_list']))
        if (self.__root['active_learning'] is not None):
            self.__loggers = []
            message = 'You are performing active learning, system data ' + \
                      'loggers will be removed.'
            _w.warn(message)

    def __check_trajectories(self) -> None:
        """
        Check trajectory settings.
        """
        has_restart_flag = False
        for trajectory in self.__trajectories:
            if (hasattr(trajectory, 'is_restart') and trajectory.is_restart):
                has_restart_flag = True
                break
        if (not has_restart_flag):
            self.__trajectories.append(
                _mdapps.trajectories.H5WRITER(
                    self.__root['run']['label'] + '.restart.h5', interval=10,
                    restart_file=True))
        if (self.__root['active_learning'] is not None):
            self.__trajectories = []
            message = 'You are performing active learning, trajectory ' + \
                      'writers will be removed.'
            _w.warn(message)

    def __parse_trajectories(self) -> None:
        """
        Parse the trajectory information.
        """
        self.__trajectories = []
        trajectories = self.__root['trajectory']
        if (trajectories is None):
            self.__trajectories.append(
                _mdapps.trajectories.H5WRITER(
                    self.__root['run']['label'] + '.trajectory.h5',
                    interval=1, write_virial=True, write_forces=True,
                    append=bool(self.__root['run']['restart_from'])))
        else:
            for index, trajectory in enumerate(trajectories):
                trajectory = self.__normalize_table(trajectory, 'trajectory')
                if (trajectory['format'] is None):
                    trajectory_format = 'H5'
                else:
                    trajectory_format = trajectory['format'].upper()
                if (trajectory['prefix'] is None):
                    prefix = self.__root['run']['label']
                else:
                    prefix = trajectory['prefix']
                if (trajectory['interval'] is None):
                    interval = 1
                else:
                    interval = trajectory['interval']
                if (trajectory_format == 'H5'):
                    if (trajectory['is_restart_file']):
                        file_name = prefix + '.restart.h5'
                    else:
                        file_name = prefix + '.trajectory.h5'
                    writer = _mdapps.trajectories.H5WRITER(
                        file_name, interval=interval, write_virial=True,
                        write_velocities=bool(trajectory['write_velocities']),
                        write_forces=bool(trajectory['write_forces']),
                        wrap_positions=bool(trajectory['wrap_positions']),
                        append=bool(self.__root['run']['restart_from']),
                        restart_file=bool(trajectory['is_restart_file']),
                        use_double=bool(trajectory['use_float64']),
                        potential_list=trajectory['potential_list'])
                    self.__trajectories.append(writer)
                elif (trajectory_format == 'EXYZ'):
                    file_name = prefix + '.trajectory.xyz'
                    writer = _mdapps.trajectories.EXYZWRITER(
                        file_name, interval=interval, write_virial=True,
                        write_velocities=bool(trajectory['write_velocities']),
                        write_forces=bool(trajectory['write_forces']),
                        wrap_positions=bool(trajectory['wrap_positions']),
                        append=bool(self.__root['run']['restart_from']),
                        potential_list=trajectory['potential_list'],
                        energy_shift=trajectory['energy_shift'])
                    self.__trajectories.append(writer)
                else:
                    message = 'Unknown trajectory format ' + \
                              trajectory['format']
                    raise RuntimeError(message)
                if (trajectory['potential_list'] is None):
                    tmp = list(range(0, len(self.__potential_generators)))
                else:
                    tmp = trajectory['potential_list']
                for i in tmp:
                    potential_name = self.__potential_generators[i][0]
                    if (potential_name == 'PLUMED'):
                        message = 'The forces and energies in the ' + \
                                  'trajectory file "{:s}" may include ' + \
                                  'contributions of PLUMED bias ' + \
                                  'potentials! MAKE SURE THAT THIS IS ' + \
                                  'WHAT YOU WANT!!'
                        _w.warn(message.format(file_name))
        self.__check_trajectories()

    def __parse_scripts(self):
        """
        Set up the post-step scripts.
        """
        scope = {}
        self.__scripts = []
        scripts = self.__root['script']
        if (scripts is None):
            return
        for index, script in enumerate(scripts):
            script = self.__normalize_table(script, 'script')
            if (script['interval'] is None):
                interval = 1
            else:
                interval = script['interval']
            if (script['initialize'] is not None):
                exec(script['initialize'], scope)
            else:
                scope['initialize'] = lambda: None
            exec(script['update'], scope)
            if ('update' not in scope.keys()):
                message = 'The function name defined in the "update" key ' + \
                          'of the [[script]] table(s) must be "update"!'
                raise RuntimeError(message)
            obj = _mdapps.utils.POSTSTEPOBJWRAPPER(
                scope['update'], scope['initialize'], interval)
            self.__scripts.append(obj)

    def __set_up_simulation(self):
        """
        Set up the simulation protocol.
        """
        for generator in self.__potential_generators:
            self.__system.potentials.append(generator[1]())
        self.__simulation = _mdapps.simulations.SIMULATION(
            self.__system, self.__integrator, barostat=self.__barostat,
            loggers=self.__loggers, trajectories=self.__trajectories)
        for obj in self.__scripts:
            obj.bind_integrator(self.__integrator)
            self.__simulation.post_step_objects.append(obj)

    def __parse_active_learning(self):
        """
        Parse the active learning information.
        """
        protocol = self.__root['active_learning']
        protocol = self.__normalize_table(protocol, 'active_learning')
        if (protocol['reference_potentials'] is None):
            reference_potentials = list(
                range(0, len(self.__potential_generators)))
            for i in range(0, len(self.__potential_generators)):
                if (self.__potential_generators[i][0] == 'PLUMED'):
                    reference_potentials.pop(i)
                elif (self.__potential_generators[i][0] == 'NEP'):
                    reference_potentials.pop(i)
        else:
            reference_potentials = protocol['reference_potentials']
        for i in reference_potentials:
            if (i >= len(self.__potential_generators)):
                message = 'Unknown potential index: {:d} in the ' + \
                          '[active_learning] table!'
                raise IndexError(message.format(i))
            if (self.__potential_generators[i][0] == 'PLUMED'):
                message = 'You are using PLUMED as one of the reference ' + \
                          'potential! You should ensure that you know ' + \
                          'what you are doing!'
                _w.warn(message)
            if (self.__potential_generators[i][0] == 'NEP'):
                message = 'You are using NEP as one of the reference ' + \
                          'potential! You should ensure that you know ' + \
                          'what you are doing!'
                _w.warn(message)
        if (protocol['initial_training_set'] is not None):
            protocol['initial_training_set'] = \
                _os.path.abspath(protocol['initial_training_set'])
        if (protocol['initial_testing_set'] is not None):
            protocol['initial_testing_set'] = \
                _os.path.abspath(protocol['initial_testing_set'])
        if (protocol['initial_potential_files'] is not None):
            protocol['initial_potential_files'] = [
                _os.path.abspath(file) for file in
                protocol['initial_potential_files']]
        generators = [g[1] for g in self.__potential_generators]
        energy_shift = protocol['energy_shift']
        use_tabulating = bool(protocol['use_tabulating'])
        for key in protocol.copy().keys():
            if (protocol[key] is None):
                protocol.pop(key)
        self.__trainer = _mdapps.active_learning.ACTIVELEARNING(
            self.__system, self.__integrator, generators, reference_potentials,
            protocol, protocol['nep_options'], protocol['nep_command'],
            use_tabulating, self.__scripts, energy_shift)

    def run(self):
        """
        Run the simulation.
        """
        if (self.__trainer is not None):
            for i in range(0, self.__root['active_learning']['n_iterations']):
                self.__trainer.run()
        else:
            if (self.__root['run']['restart_from'] is not None):
                self.__simulation.restart_from(
                    self.__root['run']['restart_from'])
            self.__simulation.run(self.__root['run']['n_steps'])

    @property
    def file_name(self):
        """
        Name of the configure file.
        """
        return self.__file_name
