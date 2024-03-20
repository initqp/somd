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
Wrapper of NumPy's RNGs.
"""

import numpy as _np
import pickle as _pl
import base64 as _bs

__all__ = ['RNG', 'LEGACYRNG']


class RNG(_np.random.Generator):
    """
    Wrapper of NumPy's RNGs.

    Parameters
    ----------
    generator_name : str
        Name of the generator.
    seed : int
        The seed value.
    """

    def __init__(
        self, generator_name: str = 'PCG64', seed: int = None
    ) -> None:
        """
        Set up the RNG.
        """
        self.__generator_name = generator_name
        generator = eval('_np.random.{:s}'.format(generator_name))(seed)
        super().__init__(generator)

    def seed(self, seed: int):
        """
        Re-seeding the RNG.

        Parameters
        ----------
        seed : int
            The seed value.
        """
        rng = RNG(self.__generator_name, seed)
        self.bit_generator.state = rng.bit_generator.state

    def set_state_from_string(self, state_string: str) -> None:
        """
        Set RNG state from string.

        Parameters
        ----------
        state_string : str
            The state string returned by the `state_string` property.
        """
        self.bit_generator.state = _pl.loads(_bs.b64decode(state_string))

    @property
    def state_string(self) -> str:
        """
        The state string.
        """
        state = self.bit_generator.state
        return _bs.b64encode(_pl.dumps(state)).decode('utf-8')


class LEGACYRNG(_np.random.RandomState):
    """
    Wrapper of NumPy's legacy RNGs.

    Parameters
    ----------
    seed : int
        The seed value.

    Notes
    -----
    This class is for testing only, and will be removed in future.
    """

    def set_state_from_string(self, state_string: str) -> None:
        """
        Set RNG state from string.

        Parameters
        ----------
        state_string : str
            The state string returned by the `state_string` property.
        """
        self.set_state(_pl.loads(_bs.b64decode(state_string)))

    @property
    def state_string(self) -> str:
        """
        The state string.
        """
        state = self.get_state()
        return _bs.b64encode(_pl.dumps(state)).decode('utf-8')
