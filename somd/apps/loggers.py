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
The simulation data loggers.
"""

import os as _os
from somd import core as _mdcore
from . import utils as _apputils

__all__ = ['DEFAULTCSVLOGGER']


class DEFAULTCSVLOGGER(_apputils.post_step.POSTSTEPOBJ):
    """
    Write information about the system to a CSV file.

    Parameter
    ---------
    file_name : str
        Name of the log file.
    groups : List[int]
        Indices of the groups to write out the information. If no groups are
        specified, information about all groups will be written.
    interval : int
        Step interval of writing the log file.
    append : bool
        If append to the old trajectory.
    delimiter : str
        The delimiter to use.
    format_str : str
        The ascii data format string. E.g.: '{:-10.5e}'.
    """

    def __init__(self,
                 file_name: str,
                 groups: list = None,
                 interval: int = 1,
                 append: bool = False,
                 delimiter: str = ' , ',
                 potential_list: list = None,
                 format_str: str = '{:-15.10e}') -> None:
        """
        Create a DEFAULTCSVLOGGER instance.
        """
        self.__fp = None
        self.__groups = groups
        self.__append = append
        self.__file_name = file_name
        self.__delimiter = delimiter
        self.__format_str = format_str
        self.__potential_list = potential_list
        super().__init__(interval)

    def __del__(self) -> None:
        """
        Close the file on exit.
        """
        if (self.__fp is not None):
            self.__fp.close()

    def __write_file_header(self) -> None:
        """
        Write the file header.
        """
        header = '# step,'
        header += ' potential energy (kJ/mol),'
        header += ' pressure xx (kJ/mol/nm^3),'
        header += ' pressure yy (kJ/mol/nm^3),'
        header += ' pressure zz (kJ/mol/nm^3),'
        header += ' pressure xy (kJ/mol/nm^3),'
        header += ' pressure xz (kJ/mol/nm^3),'
        header += ' pressure yz (kJ/mol/nm^3),'
        for i in self.__groups:
            header += ' group {} kinetic energy (kJ/mol),'.format(i)
            header += ' group {} temperature (K),'.format(i)
        if (self.integrator is not None):
            header += ' effective energy correction (kJ/mol),'
        header += ' volume (nm^3)'
        print(header, file=self.__fp)

    def bind_integrator(self, integrator: _mdcore.integrators.INTEGRATOR) \
            -> None:
        """
        Bind an integrator.
        """
        super().bind_integrator(integrator)
        self.__system = integrator.system
        if (self.__groups is None):
            self.__groups = list(range(0, len(self.__system.groups)))
        if (len(self.__groups) > len(self.__system.groups)):
            raise IndexError('Too many groups!')
        if (max(self.__groups) >= len(self.__system.groups)):
            raise IndexError('Invalid group: {}'.format(max(self.__groups)))

    def initialize(self) -> None:
        """
        Open the file for writing.
        """
        super().initialize()
        if (not _os.path.exists(self.file_name)):
            self.__append = False
        if (self.__append):
            self.__fp = open(self.file_name, 'a')
        else:
            _apputils.backup.back_up(self.file_name)
            self.__fp = open(self.file_name, 'w')
            self.__write_file_header()

    def write(self) -> None:
        """
        Write to the log file immediately.
        """
        if (not self.initialized):
            self.initialize()
        if (self.__potential_list is None):
            energy_potential = self.__system.energy_potential
        else:
            energy_potential = 0
            for i in self.__potential_list:
                energy_potential += \
                    self.__system.potentials[i].energy_potential[0]
        print(self.step, file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(energy_potential),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[0, 0]),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[1, 1]),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[2, 2]),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[0, 1]),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[0, 2]),
              file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.pressures[1, 2]),
              file=self.__fp, end=self.__delimiter)
        for i in self.__groups:
            print(self.__format_str.format(
                self.__system.groups[i].energy_kinetic),
                file=self.__fp, end=self.__delimiter)
            print(self.__format_str.format(
                self.__system.groups[i].temperature),
                file=self.__fp, end=self.__delimiter)
        if (self.integrator is not None):
            print(self.__format_str.format(
                self.integrator.energy_effective),
                file=self.__fp, end=self.__delimiter)
        print(self.__format_str.format(self.__system.volume),
              file=self.__fp, end='')
        print('', file=self.__fp)
        self.__fp.flush()

    def update(self) -> None:
        """
        Write the log file according to the defined step interval.
        """
        if (super().update()):
            self.write()

    @property
    def file_name(self) -> str:
        """
        The file name.
        """
        return self.__file_name
