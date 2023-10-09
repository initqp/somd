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

import io as _io
import os as _os
import re as _re
import uuid as _id
import numpy as _np
import typing as _tp
import mdtraj as _md
import shutil as _sh
import signal as _sg
import subprocess as _sp
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['SIESTA', 'create_siesta_potential', 'create_siesta_generator']


class SIESTA(_mdcore.potential_base.POTENTIAL):
    """
    The ab initio potential calculated with the SIESTA code [1].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    label : str
        Basename name of the SIESTA input file.
    siesta_command : str
        Command of submitting a SIESTA job.

    References
    ----------
    .. [1] Garcia, Alberto, et al. "Siesta: Recent developments and
           applications." The Journal of chemical physics 152.20 (2020):
           204108.
    """

    def __init__(self,
                 atom_list: list,
                 label: str,
                 siesta_command: str = 'siesta') -> None:
        """
        Create a SIESTA instance.
        """
        self.__siesta_pid = 0
        self.__initialized = False
        self.__args = [atom_list, label, siesta_command]
        super().__init__(atom_list)
        label = self.__strip_label(label)
        self.__pipe_c = self.__work_dir + '/' + label + '.coords'
        self.__pipe_f = self.__work_dir + '/' + label + '.forces'
        self.__conversion = \
            _mdutils.constants.AVOGACONST * _mdutils.constants.ELECTCONST
        try:
            _os.remove(self.__pipe_c)
            _os.remove(self.__pipe_f)
        except:
            pass
        _os.mkfifo(self.__pipe_c)
        _os.mkfifo(self.__pipe_f)
        # we should print the PID of the SIESTA subprocess here.
        # do not use proc.pid which is the PID of the invoked shell.
        command = (siesta_command + ' ' + label + '.fdf > ' + label + '.out' +
                   ' 2>' + label + '.err & echo $!')
        cwd = _os.getcwd()
        _os.chdir(self.__work_dir)
        proc = _sp.Popen(command, shell=True, stdout=_sp.PIPE)
        _os.chdir(cwd)
        if (proc.wait(_mdutils.defaults.SIESTATIMEOUT) != 0):
            message = 'Error in setup SIESTA using command: ' + siesta_command
            raise RuntimeError(message)
        fp = self.__timeout_open(self.__pipe_c, 'w')
        print('wait', file=fp)
        fp.close()
        pid_str = proc.communicate()[0].decode('UTF8').strip()
        if (not pid_str.isdigit()):
            message = 'Unknown PID string: {} ! Check your login shell!'
            raise RuntimeError(message.format(pid_str))
        else:
            self.__siesta_pid = int(pid_str)
        self.__initialized = True

    def __strip_label(self, label: str) -> str:
        """
        Strip redundant charactes in the system label.
        """
        _os.stat(label)
        if (_os.path.dirname(label) != ''):
            self.__work_dir = _os.path.dirname(label)
            label = _os.path.basename(label)
        else:
            self.__work_dir = './'
        tmp = _os.path.splitext(label)
        if (tmp[1] != '.fdf'):
            message = 'Name of the SIESTA input file must uses an ' + \
                      'extension name of .fdf!'
            raise RuntimeError(message)
        return tmp[0]

    @staticmethod
    def __timeout_open(file_name, mode: str, **kwargs) -> _io.TextIOWrapper:
        """
        Open a file with a given timeout.
        """
        message = 'The SIESTA subprocess has been silent for more than ' + \
                  '{} seconds!'.format(_mdutils.defaults.SIESTATIMEOUT)
        command = 'raise TimeoutError("{}")'.format(message)
        _sg.signal(_sg.SIGALRM,
                   lambda signum, frame: exec(command))
        _sg.alarm(_mdutils.defaults.SIESTATIMEOUT)
        try:
            result = open(file_name, mode, **kwargs)
        except TimeoutError as e:
            raise e
        finally:
            _sg.alarm(0)
        return result

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        fp = self.__timeout_open(self.__pipe_c, 'w')
        print('begin_coords', file=fp)
        print('Ang', file=fp)
        print('eV', file=fp)
        _np.savetxt(fp, system.box * 10.0, fmt='%.20e ')
        print('%d' % self.n_atoms, file=fp)
        _np.savetxt(fp, system.positions[self.atom_list] * 10.0, fmt='%.20e ')
        print('end_coords', file=fp)
        fp.close()
        # when SIESTA does not listen to us here, it may be performing
        # calculations or just died. Thus we check its PID when an open
        # attempt is failed. If SIESTA is alive, we should wait it longer,
        # otherwise an error should be thrown.
        while (True):
            try:
                fp = self.__timeout_open(self.__pipe_f, 'r')
            except TimeoutError as e:
                if (not self.alive):
                    message = 'The SIESTA subprocess (PID: {}) is DIED!!!'
                    message = message.format(self.__siesta_pid)
                    raise RuntimeError(message) from e
            else:
                break
        header = fp.readline()
        if (header.strip() == 'error'):
            raise RuntimeError('SIESTA FATAL:' + fp.readline())
        elif (header.strip() != 'begin_forces'):
            raise RuntimeError('SIESTA FATAL: unexpected header:' + header)
        self.energy_potential[0] = _np.double(fp.readline()) * \
            self.__conversion * 0.001
        for i in range(0, 3):
            self.virial[i, :] = system.volume * self.__conversion * -1.0 * \
                _np.array(fp.readline().split(), dtype=_np.double)
        if (int(fp.readline().strip()) != self.n_atoms):
            raise RuntimeError('SIESTA FATAL: atom number mismatch!')
        for i in range(0, self.n_atoms):
            self.forces[self.atom_list[i], :] = self.__conversion * 0.01 * \
                _np.array(fp.readline().split(), dtype=_np.double)
        tailer = fp.readline()
        if (tailer.strip() != 'end_forces'):
            raise RuntimeError('SIESTA FATAL: unexpected tailer:' + tailer)
        fp.close()

    @classmethod
    def generator(cls, *args, **kwargs) -> _tp.Callable:
        """
        Return a generator of this potential.
        """
        message = 'Use the `create_siesta_generator` method to get ' + \
                  'generators of SIESTA!'
        raise NotImplementedError(message)

    def finalize(self) -> None:
        """
        Clean up.
        """
        if (self.__initialized):
            super().finalize()
            if (self.alive):
                fp = self.__timeout_open(self.__pipe_c, 'w')
                print('quit', file=fp)
                fp.close()
            _os.remove(self.__pipe_f)
            _os.remove(self.__pipe_c)

    def reset(self) -> None:
        """
        Reset SIESTA by re-submit the task.
        """
        self.finalize()
        self.__init__(*self.__args)

    @property
    def alive(self) -> bool:
        """
        If the SIESTA subprocess is still alive.
        """
        try:
            _os.kill(self.__siesta_pid, 0)
        except OSError:
            return False
        else:
            return True

    @property
    def working_directory(self) -> str:
        """
        The working working directory of the SIESTA calculations.
        """
        return self.__work_dir


def _create_siesta_input(system: _mdcore.systems.MDSYSTEM,
                         atom_list: list,
                         options: str,
                         file_dir: str = '.',
                         label: str = 'somd_tmp') -> None:
    """
    A simple function to create a SIESTA input file. This function will
    automatically setup the system information (atomic symbels, cell vectors,
    etc.) in the input file.

    Parameters
    ----------
    system: somd.core.system.MDSYSTEM
        The simulated system.
    atom_list : List[int]
        Indices of atoms included by this potential.
    options: str
        Auxiliary options to be added into the input file.
    file_dir:
        Directory of the fdf file.
    label: str
        Basename of the file.
    """
    # check the given options.
    __needed_options = ['SystemName',
                        'SystemLabel',
                        'MD.TypeOfRun',
                        'NumberOfAtoms',
                        'LatticeVectors',
                        'LatticeConstant',
                        'NumberOfSpecies',
                        'ChemicalSpeciesLabel',
                        'AtomicCoordinatesFormat',
                        'AtomicCoordinatesAndAtomicSpecies']
    option_list = _re.split(' |\t|\n|\r', options)
    option_list = [x.lower() for x in option_list]
    for o in __needed_options:
        if (o.lower in option_list):
            message = 'SIESTA option "{}" should not be specified!'
            raise RuntimeError(message.format(o))
    # write the file.
    fp = open(file_dir + '/' + label + '.fdf', 'w')
    atomic_symbol_list = [system.atomic_symbols[i] for i in atom_list]
    atomic_symbol_list = list(set(atomic_symbol_list))
    atomic_symbol_list.sort()
    print('### Sample input file generated by SOMD. ###\n', file=fp)
    print('### SOMD generated options ###\n', file=fp)
    print('SystemName {:s}'.format(system._label), file=fp)
    print('SystemLabel {:s}'.format(label), file=fp)
    print('NumberOfAtoms {:d}'.format(len(atom_list)), file=fp)
    print('NumberOfSpecies {:d}'.format(len(atomic_symbol_list)), file=fp)
    print('MD.TypeOfRun forces', file=fp)
    print('%block ChemicalSpeciesLabel', file=fp)
    for i in range(0, len(atomic_symbol_list)):
        index = _md.element.get_by_symbol(atomic_symbol_list[i]).number
        print('   ', end='   ', file=fp)
        print(i + 1, end='   ', file=fp)
        print(index, end='   ', file=fp)
        print(atomic_symbol_list[i], file=fp)
    print('%endblock ChemicalSpeciesLabel', file=fp)
    print('LatticeConstant 1.0 Ang', file=fp)
    print('%block LatticeVectors', file=fp)
    _np.savetxt(fp, system.box * 10, fmt="%16.10f")
    print('%endblock LatticeVectors', file=fp)
    print('AtomicCoordinatesFormat  NotScaledCartesianAng', file=fp)
    print('%block AtomicCoordinatesAndAtomicSpecies', file=fp)
    for i in range(0, len(atom_list)):
        print('{:20.10f}'.format(system.positions[atom_list[i], 0] * 10),
              end=' ', file=fp)
        print('{:20.10f}'.format(system.positions[atom_list[i], 1] * 10),
              end=' ', file=fp)
        print('{:20.10f}'.format(system.positions[atom_list[i], 2] * 10),
              end=' ', file=fp)
        symbol = system.atomic_symbols[atom_list[i]]
        print(atomic_symbol_list.index(symbol) + 1, file=fp)
    print('%endblock AtomicCoordinatesAndAtomicSpecies', file=fp)
    print('\n### User specified options ###\n', file=fp)
    print(options, file=fp)
    fp.close()


def create_siesta_potential(system: _mdcore.systems.MDSYSTEM,
                            atom_list: list,
                            options: str,
                            siesta_command: str = 'siesta',
                            pseudopotential_dir: str = './') -> SIESTA:
    """
    Prepare a temporary directory to start the SIESTA calculations.

    Parameters
    ----------
    system: somd.core.system.MDSYSTEM
        The simulated system.
    atom_list : List[int]
        Indices of atoms included by this potential.
    options: str
        Auxiliary options to be added into the input file, e.g. the basis size
        and the functional. There is an example of this parameter:
        options = r'''
        xc.functional          GGA
        xc.authors             revPBE

        PAO.BasisSize          TZ2P
        Mesh.Cutoff            300 Ry
        PAO.EnergyShift        10 meV
        PAO.SoftDefault        T

        DM.MixingWeight        0.1
        DM.Tolerance           1.d-5
        DM.UseSaveDM           T
        DM.History.Depth       5

        Diag.Algorithm         ELPA-1stage
        SolutionMethod         diagon
        ElectronicTemperature  5 meV
        '''
    siesta_command : str
        Command of submitting a SIESTA job.
    pseudopotential_dir : str
        Directory of the pseudopotential files.

    Notes
    -----
    Basenames of the pseudopotential file must be the same as symbols of the
    elements in the simulated system. E.g. there are three elements: H, O and
    C in the simulated system, then the pseudopotential files must be named as:
    H.psf, O.psf and C.psf (or H.vps, O.vps and C.vps).
    """
    work_dir = _os.getcwd() + '/SOMD_TMP_' + _id.uuid4().hex
    _os.mkdir(work_dir)
    pseudopotential_dir = pseudopotential_dir + '/'
    files = _os.listdir(pseudopotential_dir)
    for e in list(set(system.atomic_symbols)):
        if ((e + '.psml') in files):
            fn = e + '.psml'
            _sh.copyfile(pseudopotential_dir + fn, work_dir + '/' + fn)
        elif ((e + '.psf') in files):
            fn = e + '.psf'
            _sh.copyfile(pseudopotential_dir + fn, work_dir + '/' + fn)
        elif ((e + '.vps') in files):
            fn = e + '.vps'
            _sh.copyfile(pseudopotential_dir + fn, work_dir + '/' + fn)
        else:
            _os.rmdir(work_dir)
            message = 'Can not find pseudopotential file for element ' + e
            raise RuntimeError(message)
    _create_siesta_input(system, atom_list, options, file_dir=work_dir)
    return SIESTA(atom_list, work_dir + '/somd_tmp.fdf', siesta_command)


def create_siesta_generator(system: _mdcore.systems.MDSYSTEM,
                            atom_list: list,
                            options: str,
                            siesta_command: str = 'siesta',
                            pseudopotential_dir: str = './') -> _tp.Callable:
    """
    Return a generator of this potential.

    Parameters
    ----------
    system: somd.core.system.MDSYSTEM
        The simulated system.
    atom_list : List[int]
        Indices of atoms included by this potential.
    options: str
        Auxiliary options to be added into the input file, e.g. the basis
        size and the functional. There is an example of this parameter:
        options = r'''
        xc.functional          GGA
        xc.authors             revPBE

        PAO.BasisSize          TZ2P
        Mesh.Cutoff            300 Ry
        PAO.EnergyShift        10 meV
        PAO.SoftDefault        T

        DM.MixingWeight        0.1
        DM.Tolerance           1.d-5
        DM.UseSaveDM           T
        DM.History.Depth       5

        Diag.Algorithm         ELPA-1stage
        SolutionMethod         diagon
        ElectronicTemperature  5 meV
        '''
    siesta_command : str
        Command of submitting a SIESTA job.
    pseudopotential_dir : str
        Directory of the pseudopotential files.

    Notes
    -----
    Basenames of the pseudopotential file must be the same as symbols of
    the elements in the simulated system. E.g. there are three elements:
    H, O and C in the simulated system, then the pseudopotential files must
    be named as: H.psf, O.psf and C.psf (or H.vps, O.vps and C.vps).
    """
    pseudopotential_dir = _os.path.abspath(pseudopotential_dir)
    return lambda: create_siesta_potential(system, atom_list, options,
                                           siesta_command,
                                           pseudopotential_dir)
