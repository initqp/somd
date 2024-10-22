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
Base class for all post-step objects.
"""

import abc as _ab
import typing as _tp
from somd import core as _mdcore
from somd import utils as _mdutils

__all__ = ['POSTSTEPOBJ', 'POSTSTEPOBJWRAPPER']


class POSTSTEPOBJ(_ab.ABC):
    """
    The base class for all post-step objects.

    Parameters
    ----------
    interval : int
        The step interval of invoking this object.
    """

    def __init__(self, interval: int) -> None:
        """
        Create a POSTSTEPOBJ instance.
        """
        self.__interval = interval
        # Reference to the integrator.
        self.__integrator = None
        # If initialized.
        self.__initialized = False

    def summary(self) -> str:
        """
        Show information about the object.
        """
        result = '{}\n'.format(self.__class__.__name__)
        result += '┣━ interval: {}\n'.format(self.__interval)
        result += '┗━ END'
        return result

    def bind_integrator(
        self, integrator: _mdcore.integrators.INTEGRATOR
    ) -> None:
        """
        Bind an integrator.
        """
        if self.__integrator is not None:
            message = 'Rebinding an integrators to a post step object.'
            _mdutils.warning.warn(message)
        if integrator.system is None:
            raise RuntimeError('Integrator has not bind a system!')
        self.__integrator = integrator

    def initialize(self) -> None:
        """
        Initialize the object.
        """
        if self.__integrator is None:
            message = 'Must bind to an integrator before initialization!'
            raise RuntimeError(message)
        self.__initialized = True

    def finalize(self) -> None:
        """
        Finalize the object.
        """
        pass

    @_ab.abstractmethod
    def update(self) -> bool:
        """
        Perform the post step task.
        """
        if (self.step % self.__interval) == 0:
            return True
        else:
            return False

    @property
    def integrator(self) -> _mdcore.integrators.INTEGRATOR:
        """
        Refernce to the integrator.
        """
        return self.__integrator

    @property
    def initialized(self) -> bool:
        """
        Check if the object has been initialization.
        """
        return self.__initialized

    @property
    def interval(self) -> int:
        """
        The step interval of invoking this object.
        """
        return self.__interval

    @interval.setter
    def interval(self, v: int) -> int:
        """
        Set the step interval of invoking this object.
        """
        self.__interval = int(v)

    @property
    def step(self) -> int:
        """
        Number of called times.
        """
        return self.__integrator.step


class POSTSTEPOBJWRAPPER(POSTSTEPOBJ):
    """
    The wrapper class for all post-step objects.

    Parameters
    ----------
    update_function : function
        The update function.
    initialize_function : function
        The initialization function.
    interval : int
        The step interval of invoking this object.
    """

    def __init__(
        self,
        update_function: _tp.Callable,
        initialize_function: _tp.Callable = None,
        interval: int = 1,
    ):
        """
        Create a POSTSTEPOBJWRAPPER instance.
        """
        self.__update_function = update_function
        self.__initialize_function = initialize_function
        super().__init__(interval)

    def initialize(self) -> None:
        """
        Initialize the object.
        """
        if self.__initialize_function is not None:
            self.__initialize_function()

    def update(self) -> None:
        """
        Update the object.
        """
        if super().update():
            self.__update_function(self.integrator)
