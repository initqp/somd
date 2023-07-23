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
The shooting method.
"""

import numpy as _np
import warnings as _w
import somd.core as _mdcore

__all__ = ['shoot']


def shoot(system: _mdcore.systems.MDSYSTEM) -> None:
    """
    Modify the momentum part of one phase point while keeping the phase density
    unchanged.

    Parameters
    ----------
    system: somd.core.systems.MDSYSTEM
        The system to shoot.
    temperatire : float
        The temperature of the shooting point. If this parameters is given,
        the kinetic energy of the original trajectory will not be invoked.
    """
    target_groups = [g for g in system.groups if (g.n_atoms == system.n_atoms)]
    if (len(target_groups) == 0):
        system.groups.create_from_dict({'atom_list': range(0, system.n_atoms)})
        target_groups = [system.groups[len(system.groups) - 1]]
    target_group = target_groups[0]
    temperature = target_group.temperature.copy()
    target_group.velocities = 0
    target_group.add_velocities_from_temperature(temperature)
    if (len(system.constraints) > 0):
        try:
            energy_kinetic = target_group.energy_kinetic
            system.constraints.rattle_constrain_p(1E-4)
            target_group.velocities *= _np.sqrt(energy_kinetic /
                                                target_group.energy_kinetic)
        except:
            message = 'Failed to apply RATTLE momentum constraints ' + \
                      'after shooting!'
            _w.warn(message)
