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
import warnings as _w
from somd import core as _mdcore
from somd.constants import CONSTANTS as _c

__all__ = ['PLUMED']


class PLUMED(_mdcore.potential_base.POTENTIAL):
    """
    Wrapper to the community-developed plugin for molecular dynamics (PLUMED),
    version 2 [1,2].

    Parameters
    ----------
    atom_list : List(int)
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
    cv_names : List(dict)
        Names and components of the collective variables to save. For example:
        cv_names = [{'d1': ['x']}, {'d2': []}, [{'d3': ['x', 'y', 'z']}]

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
                 cv_names: list = []) -> None:
        """
        Create a PLUMED instance.
        """
        super().__init__(atom_list)
        # treat PLUMED as a local dependency
        try:
            import plumed
            _w.warn('PLUMED kernel outputs begin:')
            self.__plumed = plumed.Plumed()
            _w.warn('PLUMED kernel outputs end.')
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
        self.__plumed.cmd("setMDEngine", "SOMD")
        self.__plumed.cmd("setMDTimeUnits", 1.0)
        self.__plumed.cmd("setMDMassUnits", 1.0)
        self.__plumed.cmd("setMDEnergyUnits", 1.0)
        self.__plumed.cmd("setMDLengthUnits", 1.0)
        self.__plumed.cmd("setMDChargeUnits", 1.0)
        self.__plumed.cmd("setRestart", restart)
        self.__plumed.cmd("setTimestep", timestep)
        self.__plumed.cmd("setPlumedDat", file_name)
        self.__plumed.cmd("setNatoms", len(atom_list))
        self.__plumed.cmd("setLogFile", log_file)
        if (temperature is not None):
            self.__plumed.cmd("setKbT", temperature * _c.BOLTZCONST)
        self.__plumed.cmd("init")
        self.__set_up_cv(cv_names)
        self.__restart = restart
        self.__stop_flag = 0
        self.__step = 1

    def __set_up_cv(self, cv_names: list) -> None:
        """
        Set up collective variables data.

        Parameters
        ----------
        cv_names : List(dict)
            Names and components of the collective variables to save. For
            example:
            cv_names = [{'d1': ['x']}, {'d2': []}, [{'d3': ['x', 'y', 'z']}]
        """
        self.__cv_values = []
        self.__cv_names = cv_names
        for cv in cv_names:
            if (len(cv.items()) != 1):
                message = 'Each element of the CV list should contains ' + \
                          'one item!'
                raise RuntimeError(message)
            name = list(cv.keys())[0]
            if (len(cv[name]) == 0):
                data = _np.zeros(1, dtype=_np.double)
                command = "setMemoryForData " + name
                self.__plumed.cmd(command, data)
                self.__cv_values.append([data])
            else:
                values = []
                for component in cv[name]:
                    data = _np.zeros(1, dtype=_np.double)
                    command = "setMemoryForData " + name + '.' + component
                    self.__plumed.cmd(command, data)
                    values.append(data)
                self.__cv_values.append(values)

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
        self.__plumed.cmd("setStep", self.__step)
        self.__plumed.cmd("setVirial", self.virial)
        self.__plumed.cmd("setForces", self.forces)
        self.__plumed.cmd("setBox", system.box)
        self.__plumed.cmd("setPositions", system.positions)
        self.__plumed.cmd("setMasses", system.masses.reshape(-1))
        self.__plumed.cmd("prepareCalc")
        self.__plumed.cmd("performCalc")
        self.__plumed.cmd("getBias", self.energy_potential)
        self.__plumed.cmd("setStopFlag", self.__stop_flag)
        self.virial[:] *= -1.0
        self.__step += 1

    def finalize(self) -> None:
        """
        Clean up.
        """
        super().finalize()
        self.__plumed.finalize()

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
        return self.__stop_flag

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
    def cv_values(self) -> list:
        """
        Saved collective variable values.
        """
        return self.__cv_values.copy()
