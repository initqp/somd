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
NEP utils.
"""

import os as _os
import re as _re
import numpy as _np
from somd import core as _mdcore

__all__ = [
    'cat_exyz',
    'get_loss',
    'make_nep_in',
    'get_potentials_msd',
    'check_nep_parameters',
]


def cat_exyz(set_in: list, set_out: str) -> None:
    """
    Combine two EXYZ training sets.

    Parameters
    ----------
    set_in : List[str]
        Names of input training sets.
    set_out : str
        Name of output training set.
    """
    fp = open(set_out, 'w')
    for fn in set_in:
        if fn is not None:
            fp_in = open(fn, 'r')
            for l in fp_in:
                fp.write(l)
            fp_in.close()
    fp.close()


def get_potentials_msd(
    potentials: list, system: _mdcore.systems.MDSYSTEM
) -> float:
    """
    Return the maximum standard deviation of the forces calculated by
    a list of potentials.

    Parameters
    ----------
    potentials: List[somd.core.potential_base.POTENTIAL]
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
            sd[i] += tmp**2
        sd[i] /= len(potentials)
    return _np.max(_np.sqrt(sd))


def get_loss(file_name: str) -> list:
    """
    Read losses from a loss.out file.

    Parameters
    ----------
    file_name : str
        Name of the loss.out file.
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


def check_nep_parameters(nep_parameters: str, symbols: list) -> bool:
    """
    Check the NEP training parameters.

    Parameters
    ----------
    nep_parameters : str
        The keywords and corresponding values to be used in a nep.in file.
        Different keywords should be split by newlines, as in the nep.in file.
    symbols : List[str]
        Symbols of each atom in the system.

    Returns
    -------
    If element symbols should be written in the nep.in file manually.
    """
    write_symbols = True
    parameters = nep_parameters.replace('\\n', '\n')
    for line in _re.split('\n', parameters):
        l = [i for i in line.strip().split(' ') if i != '']
        if l != [] and l[0].lower() == 'type':
            e_nep = l[2:]
            if len(e_nep) != int(l[1]):
                message = 'Wrong Number of elements in NEP parameters!'
                raise RuntimeError(message)
            e_lack = [e for e in symbols if e not in e_nep]
            e_unknown = [e for e in e_nep if e not in symbols]
            if len(e_lack) != 0:
                message = 'Lack element {} in NEP parameters!'
                raise RuntimeError(message.format(list(set(e_lack))))
            if len(e_unknown) != 0:
                message = 'Unknown element {} in NEP parameters!'
                raise RuntimeError(message.format(list(set(e_unknown))))
            write_symbols = False
    return write_symbols


def make_nep_in(nep_parameters: str, symbols: list = None) -> None:
    """
    Write the nep.in file.

    Parameters
    ----------
    nep_parameters : str
        The keywords and corresponding values to be used in a nep.in file.
        Different keywords should be split by newlines, as in the nep.in file.
    symbols : List[str]
        Symbols of each atom in the system. If this option is None, the symbols
        will not be written.
    """
    fp = open('nep.in', 'w')
    if symbols is not None:
        symbols = list(set(symbols))
        symbols.sort()
        print('type {:d}'.format(len(symbols)), end='', file=fp)
        for s in symbols:
            print(' ' + s, end='', file=fp)
        print('\n', file=fp)
    print(nep_parameters, file=fp)
    fp.close()
