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

import os as _os
import numpy as _np
import typing as _tp
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['PLUMED']


class PLUMED(_mdcore.potential_base.POTENTIAL):
    """
    Wrapper to the community-developed plugin for molecular dynamics (PLUMED),
    version 2 [1,2].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    file_name : str
        Name of the PLUMED input file.
    timestep : float
        Timestep of the integrator. In unit of (ps).
    temperature : float
        Temperature of the integrator. In unit of (K).
    restart : bool
        If restart the calculation.
    output_prefix : str
        Prefix of the output file.
    cv_names : List[dict]
        Names and components of the collective variables to save. For example:
        cv_names = [{'d1': 'x'}, {'d1': 'y'}, {'d2': ''}]
    extra_cv_names : List[str]
        Names of the extra CVs.
    extra_cv_potential_index : int
        Index of the potential calculator for reading extra CV gradients from.

    References
    ----------
    .. [1] Tribello, Gareth A., et al. "PLUMED 2: New feathers for an old
           bird." Computer physics communications 185.2 (2014): 604-613.
    .. [2] The PLUMED consortium. "Promoting transparency and reproducibility
           in enhanced molecular simulations." Nature methods 16.8 (2019):
           670-673.
    """

    def __init__(self,
                 atom_list: list,
                 file_name: str,
                 timestep: float,
                 temperature: float = None,
                 restart: bool = False,
                 output_prefix: str = None,
                 cv_names: list = [],
                 extra_cv_names: list = [],
                 extra_cv_potential_index: int = None) -> None:
        """
        Create a PLUMED instance.
        """
        super().__init__(atom_list)
        self.__args = [atom_list, _os.path.abspath(file_name), timestep,
                       temperature, restart, output_prefix, cv_names]
        _os.environ['PLUMED_LOAD_NODEEPBIND'] = 'yes'
        _os.environ['PLUMED_TYPESAFE_IGNORE'] = 'yes'
        # treat PLUMED as a local dependency
        try:
            import plumed
            _mdutils.warning.warn('PLUMED kernel outputs begin:')
            self.__plumed = plumed.Plumed()
            _mdutils.warning.warn('PLUMED kernel outputs end.')
        except:
            raise ImportError('you need to have both the PLUMED python ' +
                              'wrapper and PLUMED_KERNEL installed to use ' +
                              'the PLUMED functionality!')
        if (output_prefix is None):
            log_file = _os.path.basename(file_name) + '.log'
        elif (output_prefix == ''):
            log_file = _os.devnull
        else:
            log_file = output_prefix + '.log'
        self.__plumed.cmd('setMDEngine', 'SOMD')
        self.__plumed.cmd('setMDTimeUnits', 1.0)
        self.__plumed.cmd('setMDMassUnits', 1.0)
        self.__plumed.cmd('setMDEnergyUnits', 1.0)
        self.__plumed.cmd('setMDLengthUnits', 1.0)
        self.__plumed.cmd('setMDChargeUnits', 1.0)
        self.__plumed.cmd('setRestart', restart)
        self.__plumed.cmd('setTimestep', timestep)
        self.__plumed.cmd('setPlumedDat', file_name)
        self.__plumed.cmd('setNatoms', len(atom_list))
        self.__plumed.cmd('setLogFile', log_file)
        if (temperature is not None):
            kb_T = temperature * _mdutils.constants.BOLTZCONST
            self.__plumed.cmd('setKbT', kb_T)
        if (len(extra_cv_names) != 0):
            self.__set_up_extra_cvs(extra_cv_names, extra_cv_potential_index)
            self.__extra_cv_checked = False
            self.__has_extra_cvs = True
        else:
            self.__has_extra_cvs = False
        self.__plumed.cmd('init')
        self.__set_up_cvs(cv_names)
        self.__restart = restart
        self.__stop_flag = _np.zeros(1, dtype=_np.int64)
        self.__step = 1

    def __set_up_cvs(self, cv_names: list) -> None:
        """
        Set up collective variables data.

        Parameters
        ----------
        cv_names : List[dict]
            Names and components of the collective variables to save. For
            example:
            cv_names = [{'d1': 'x'}, {'d1': 'y'}, {'d2': ''}]
        """
        self.__cv_values = []
        self.__cv_names = cv_names
        for cv in cv_names:
            if (len(cv.items()) != 1):
                message = 'Each element of the CV list should contains ' + \
                          'one item!'
                raise RuntimeError(message)
            name = list(cv.keys())[0]
            if (type(cv[name]) is str):
                data = _np.zeros(1, dtype=_np.double)
                if (cv[name] == ''):
                    command = 'setMemoryForData ' + name
                else:
                    command = 'setMemoryForData ' + name + '.' + cv[name]
                self.__plumed.cmd(command, data)
                self.__cv_values.append(data)
            else:
                message = 'Type of CV components should be str!'
                raise RuntimeError(message)

    def __set_up_extra_cvs(self,
                           extra_cv_names: list,
                           extra_cv_potential_index: int):
        """
        Set up the extra CVs.

        Parameters
        ----------
        extra_cv_names : List[str]
            Names of the extra CVs.
        extra_cv_potential_index : int
            Index of the potential calculator for reading extra CV gradients
            from.
        """
        if (extra_cv_potential_index is None):
            message = 'Index of the potential calculator for ' + \
                      'reading extra CV gradients is required!'
            raise RuntimeError(message)
        cv_shape = (len(extra_cv_names), 1)
        self.__extra_cv_names = extra_cv_names
        self.__extra_cv_index = extra_cv_potential_index
        self.__extra_cv_values = _np.zeros(cv_shape, dtype=_np.double)
        self.__extra_cv_forces = _np.zeros(cv_shape, dtype=_np.double)
        for i in range(0, cv_shape[0]):
            command = 'setExtraCV {:s}'.format(self.__extra_cv_names[i])
            self.__plumed.cmd(command, self.__extra_cv_values[i])
            command = 'setExtraCVForce {:s}'.format(self.__extra_cv_names[i])
            self.__plumed.cmd(command, self.__extra_cv_forces[i])
            message = 'Extra CV "{:s}" has been added to PLUMED.'
            _mdutils.warning.warn(message.format(self.__extra_cv_names[i]))

    def __check_extra_cv(self, system: _mdcore.systems.MDSYSTEM):
        """
        Check the extra CV setup.
        """
        potential = system.potentials[self.__extra_cv_index]
        if (not self.__extra_cv_checked):
            if (potential.atom_list.tolist() != self.atom_list.tolist()):
                message = 'The extra CV calculator (potential {:d}) ' + \
                          'and the PLUMED interface act on different ' + \
                          'atoms!'
                raise RuntimeError(message.format(self.__extra_cv_index))
            if (not hasattr(potential, 'extra_cv_values')):
                message = 'Potential {:d} can not be used to ' + \
                          'calculate extra CV!'
                raise RuntimeError(message.format(self.__extra_cv_index))
            if (not hasattr(potential, 'extra_cv_gradients')):
                message = 'Potential {:d} can not be used to ' + \
                          'calculate extra CV!'
                raise RuntimeError(message.format(self.__extra_cv_index))
            if (len(potential.extra_cv_values) != len(self.__extra_cv_values)):
                message = 'Mismatch between the number of expected and ' + \
                          'provided extra CVs and number of provided CV ' + \
                          'names ({:d} v.s. {:d})!'
                message = message.format(len(potential.extra_cv_values),
                                         len(self.__extra_cv_values))
            self.__extra_cv_checked = True

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        self.forces[:] = 0.0
        self.virial[:] = 0.0
        if (self.__has_extra_cvs):
            self.__check_extra_cv(system)
            potential = system.potentials[self.__extra_cv_index]
            self.__extra_cv_values[:] = potential.extra_cv_values
            self.__extra_cv_forces[:] = 0
        self.__plumed.cmd('setStep', self.__step)
        self.__plumed.cmd('setVirial', self.virial)
        self.__plumed.cmd('setForces', self.forces)
        self.__plumed.cmd('setBox', system.box)
        self.__plumed.cmd('setPositions', system.positions)
        self.__plumed.cmd('setMasses', system.masses.reshape(-1))
        self.__plumed.cmd('setStopFlag', self.__stop_flag)
        self.__plumed.cmd('prepareCalc')
        self.__plumed.cmd('performCalc')
        self.__plumed.cmd('getBias', self.energy_potential)
        if (self.__has_extra_cvs):
            potential = system.potentials[self.__extra_cv_index]
            for i in range(0, len(self.__extra_cv_values)):
                cv_forces = self.__extra_cv_forces[i]
                self.forces[:] += cv_forces * potential.extra_cv_gradients[i]
        self.virial[:] *= -1.0
        self.__step += 1

    @classmethod
    def generator(cls, *args, **kwargs) -> _tp.Callable:
        """
        Return a generator of this potential.
        """
        if 'file_name' in kwargs.keys():
            kwargs['file_name'] = _os.path.abspath(kwargs['file_name'])
        else:
            args = list(args)
            args[1] = _os.path.abspath(args[1])
        return lambda x=tuple(args), y=kwargs: cls(*x, **y)

    def finalize(self) -> None:
        """
        Clean up.
        """
        super().finalize()
        self.__plumed.finalize()

    def reset(self) -> None:
        """
        Reset the interface.
        """
        self.finalize()
        self.__init__(*self.__args)

    @property
    def restart(self) -> bool:
        """
        If restart.
        """
        return self.__restart

    @property
    def stop_flag(self) -> int:
        """
        The stop flag set by PLUMED.
        """
        return bool(self.__stop_flag[0])

    @property
    def step(self) -> int:
        """
        Current time step.
        """
        return self.__step

    @step.setter
    def step(self, s: int) -> None:
        """
        Set current time step.
        """
        self.__step = int(s)

    @property
    def cv_names(self) -> list:
        """
        Saved collective variable names.
        """
        return self.__cv_names.copy()

    @property
    def cv_values(self) -> _np.ndarray:
        """
        Saved collective variable values.
        """
        return _np.concatenate(self.__cv_values)
