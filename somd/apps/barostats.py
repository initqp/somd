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
The barostats.
"""

import numpy as _np
import typing as _tp
from somd import core as _mdcore
from . import utils as _apputils

__all__ = ['BAROSTAT']


class BAROSTAT(_apputils.post_step.POSTSTEPOBJ):
    """
    The Berendsen-type barostat [1].

    Parameters
    ----------
    pressures : List[float]
        The target pressures. In unit of (kJ/mol/nm^3). When using isotropic
        pressure controlling, this list contains one element, which stands for
        the hydrostatic pressure. When using anisotropic pressure controlling,
        this list contains six elements, which stand for the target pressures
        in the following directions: xx yy zz xy xz yz.
    beta : List[float]
        The isothermal compressibilities. In unit of (nm^3/(kJ/mol)). When
        using isotropic pressure controlling, this list contains one element,
        which stands for isotropic compressibility. When using anisotropic
        pressure controlling, this list contains six elements, which stand for
        the compressibilities in the following directions: xx yy zz xy xz yz.
    relaxation_time : float
        The relaxation time of the barostat. In unit of (ps).
    isotropic : bool
        If perform isotropic pressure controlling.

    Refernces
    ---------
    .. [1] Berendsen, Herman JC, et al. "Molecular dynamics with coupling to
           an external bath." The Journal of chemical physics 81.8 (1984):
           3684-3690.
    """

    def __init__(
        self,
        pressures: _tp.List[float],
        beta: _tp.List[float],
        relaxation_time: float,
    ) -> None:
        """
        Create a BAROSTAT instance.
        """
        if len(pressures) != len(beta):
            message = (
                'Length of the pressure array should '
                + 'be the same as the beta array!'
            )
            raise RuntimeError(message)
        if len(pressures) == 1:
            self.__pressures = _np.zeros(1, dtype=_np.double)
            self.__beta = _np.zeros(1, dtype=_np.double)
            self.__pressures[:] = pressures[:]
            self.__beta[:] = beta[:]
            self.__isotropic = True
        elif len(pressures) == 6:
            self.__pressures = _np.zeros((3, 3), dtype=_np.double)
            self.__beta = _np.zeros((3, 3), dtype=_np.double)
            self.__pressures[0, 0] = pressures[0]
            self.__pressures[0, 1] = pressures[3]
            self.__pressures[0, 2] = pressures[4]
            self.__pressures[1, 0] = pressures[3]
            self.__pressures[1, 1] = pressures[1]
            self.__pressures[1, 2] = pressures[5]
            self.__pressures[2, 0] = pressures[4]
            self.__pressures[2, 1] = pressures[5]
            self.__pressures[2, 2] = pressures[2]
            self.__beta[0, 0] = beta[0]
            self.__beta[0, 1] = beta[3]
            self.__beta[0, 2] = beta[4]
            self.__beta[1, 0] = beta[3]
            self.__beta[1, 1] = beta[1]
            self.__beta[1, 2] = beta[5]
            self.__beta[2, 0] = beta[4]
            self.__beta[2, 1] = beta[5]
            self.__beta[2, 2] = beta[2]
            self.__isotropic = False
        else:
            message = 'Length of the pressure array could only be 1 or 6!'
            raise RuntimeError(message)
        self.__relaxation_time = relaxation_time
        super().__init__(1)

    def __calc_mu_deterministic(self) -> _np.ndarray:
        """
        Calculate the standard mu tensor in the Berendsen barostat.
        """
        mu = _np.eye(3, dtype=_np.double)
        # fmt: off
        if self.__isotropic:
            p_hydro = self.__system.pressures.trace() / 3
            mu *= (
                1.0 - self.__beta[0] * self.__dt / 3
                / self.__relaxation_time * (self.__pressures[0] - p_hydro)
            )
        else:
            mu -= (
                self.__beta * self.__dt / 3 / self.__relaxation_time
                * (self.__pressures - self.__system.pressures)
            )
        # fmt: on
        return mu

    def bind_integrator(
        self, integrator: _mdcore.integrators.INTEGRATOR
    ) -> None:
        """
        Bind an integrator.
        """
        super().bind_integrator(integrator)
        self.__system = integrator.system
        self.__dt = integrator.timestep

    def initialize(self) -> None:
        """
        Initialize the barostat.
        """
        super().initialize()
        if len(self.__system.groups.segments) > 0:
            atom_list = list(range(0, self.__system.n_atoms))
            for seg in self.__system.groups.segments:
                for x in seg.atom_list:
                    atom_list.remove(x)
            self.__single_atoms = _np.array(atom_list, dtype=_np.int_)

    def summary(self) -> str:
        """
        Show information about the potential.
        """
        result = 'BAROSTAT\n'
        result += '┣━ isotropic: {}\n'.format(self.__isotropic)
        if self.__isotropic:
            result += '┣━ pressures: {}\n'.format(self.__pressures[0])
            result += '┣━ beta: {}\n'.format(self.__beta[0])
        else:
            f_str = '{:.5e} {:.5e} {:.5e}\n'
            result += '┣━ pressures a: ' + f_str.format(*self.__pressures[0])
            result += '┣━ pressures b: ' + f_str.format(*self.__pressures[1])
            result += '┣━ pressures c: ' + f_str.format(*self.__pressures[2])
            result += '┣━ beta a: ' + f_str.format(*self.__beta[0])
            result += '┣━ beta b: ' + f_str.format(*self.__beta[1])
            result += '┣━ beta c: ' + f_str.format(*self.__beta[2])
        result += '┣━ relaxation_time: {}\n'.format(self.__relaxation_time)
        result += '┗━ END'
        return result

    def update(self) -> None:
        """
        Perform the pressure controlling.
        """
        if not self.initialized:
            self.initialize()
        super().update()
        mu = self.__calc_mu_deterministic()
        self.__system.box[:] = mu.dot(self.__system.box.T).T
        if len(self.__system.groups.segments) > 0:
            raise NotImplementedError(
                'Can not apply barostats with constraints!'
            )
        else:
            self.__system.positions[:] = mu.dot(
                self.__system.positions_wrapped.T
            ).T

    @property
    def relaxation_time(self) -> float:
        """
        The relaxation time of the barostat. In unit of (ps).
        """
        return self.__relaxation_time

    @relaxation_time.setter
    def relaxation_time(self, t: float) -> float:
        """
        Set the relaxation time of the barostat. In unit of (ps).
        """
        self.__relaxation_time = t

    @property
    def pressures(self) -> _np.ndarray:
        """
        The target pressures of the barostat. In unit of (kJ/mol/nm^3).
        """
        return self.__pressures

    @pressures.setter
    def pressures(self, p: _tp.List[float]) -> None:
        """
        Set the target pressures of the barostat. In unit of (kJ/mol/nm^3).
        """
        self.__pressures[:] = p[:]

    @property
    def beta(self) -> _np.ndarray:
        """
        The isothermal compressibilities. In unit of (nm^3/(kJ/mol)).
        """
        return self.__beta

    @beta.setter
    def beta(self, b: _tp.List[float]) -> _np.ndarray:
        """
        Set the isothermal compressibilities. In unit of (nm^3/(kJ/mol)).
        """
        self.__beta[:] = b[:]
