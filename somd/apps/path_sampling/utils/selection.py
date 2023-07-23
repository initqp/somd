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
Select frames from paths according to the CV values.
"""

import numpy as _np

__all__ = ['calculate_frame_probabilities', 'select']


def calculate_frame_probabilities(cv_values: _np.ndarray,
                                  bias_function: callable = None,
                                  strip_terminal: bool = False) -> _np.ndarray:
    """
    Calculate probabilities of each frame in the path  according to a given
    bias function.

    Parameters
    ----------
    cv_values : numpy.ndarray
        The CV values.
    bias_function : callable
        The function that determines the weight of a given frame according
        to its CV values. For example, to only select frames with a CV
        value between 1.0 and 2.0, set this option as:
        bias_function = lambda x: (x > 1.0) * (x < 2.0)
    strip_terminal : bool
        If set probabilities of the first and last points to zero.
    """
    if (bias_function is None):
        probabilities = _np.ones(len(cv_values))
    else:
        probabilities = _np.zeros(len(cv_values))
        for i in range(0, len(cv_values)):
            probabilities[i] = bias_function(cv_values[i])
            if (probabilities[i] < 0):
                message = 'The bias function return a negative ' + \
                          'probability for path frame {:d}!'
                raise RuntimeError(message.format(i))
    if (strip_terminal):
        for i in [0, -1]:
            probabilities[i] = 0
    total = probabilities.sum()
    if (not _np.allclose(total, 0.0, 1E-16)):
        probabilities /= total
    return probabilities


def select(cv_values: _np.ndarray,
           bias_function: callable = None,
           strip_terminal: bool = False) -> dict:
    """
    Rondomly select a frame from the path according to a given bias
    function.

    Parameters
    ----------
    cv_values : numpy.ndarray
        The CV values.
    bias_function : callable
        The function that determines the weight of a given frame according
        to its CV values. For example, to only select frames with a CV
        value between 1.0 and 2.0, set this option as:
        bias_function = lambda x: (x > 1.0) * (x < 2.0)
    strip_terminal : bool
        If set probabilities of the first and last points to zero.
    """
    cumulation = 0.0
    trial = _np.random.rand()
    probabilities = calculate_frame_probabilities(cv_values, bias_function,
                                                  strip_terminal)
    if (_np.allclose(probabilities.sum(), 0.0, 1E-16)):
        message = 'The total selection probability is 0! Will select ' + \
                  'no frame from the path!'
        raise RuntimeError(message)
    for i in range(0, len(cv_values)):
        cumulation += probabilities[i]
        if (cumulation > trial):
            break
    return {'index': i, 'probability': probabilities[i]}
