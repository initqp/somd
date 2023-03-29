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
import numpy as _np
from somd import core as _mdcore

__all__ = ['cat_exyz', 'get_potentials_msd', 'get_loss']


def cat_exyz(set_in: list, set_out: str) -> None:
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


def get_potentials_msd(potentials: list,
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


def get_loss(file_name: str) -> list:
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
