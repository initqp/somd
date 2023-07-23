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
Internal represention of CV space states.
"""

import copy as _cp
import types as _ts
import numpy as _np
import base64 as _bs
import marshal as _ma

__all__ = ['STATE',
           'in_sequence',
           'get_region_indices',
           'get_reactive_regions']


class STATE(object):
    """
    States in the CV space.

    Parameters
    ----------
    indicators: callable or List(callable)
        The indicator function. This function should return 1 when the given
        CV space point belongs to the state and returns 0 otherwise. The `List`
        type is for internal usages only.
    label : str
        The label of this state.
    """

    def __init__(self, indicator: callable, label: str) -> None:
        """
        Create a STATE instance.
        """
        self.__label = label
        self.__expression = '{}'
        if (indicator.__code__.co_argcount != 1):
            message = 'The indicator functions should only take ' + \
                      'one parameter!'
            raise RuntimeError(message)
        self.__indicators = [indicator]

    def __or__(self, state: 'STATE') -> 'STATE':
        """
        Return the union of two states.
        """
        label = '({:s}|{:s})'.format(self.label, state.label)
        result = STATE(lambda x: 1, label)
        result._indicators = [*_cp.deepcopy(self._indicators),
                              *_cp.deepcopy(state._indicators)]
        result._expression = '({} or {})'.format(self._expression,
                                                 state._expression)
        return result

    def __xor__(self, state: 'STATE') -> 'STATE':
        """
        Return the exclusive or of two states.
        """
        label = '({:s}^{:s})'.format(self.label, state.label)
        result = STATE(lambda x: 1, label)
        result._indicators = [*_cp.deepcopy(self._indicators),
                              *_cp.deepcopy(state._indicators)]
        result._expression = '({} != {})'.format(self._expression,
                                                 state._expression)
        return result

    def __and__(self, state: 'STATE') -> 'STATE':
        """
        Return the intersection of two states.
        """
        label = '({:s}&{:s})'.format(self.label, state.label)
        result = STATE(lambda x: 1, label)
        result._indicators = [*_cp.deepcopy(self._indicators),
                              *_cp.deepcopy(state._indicators)]
        result._expression = '({} and {})'.format(self._expression,
                                                  state._expression)
        return result

    def __sub__(self, state: 'STATE') -> 'STATE':
        """
        Return the complementation of a given states according to this state.
        """
        label = '({:s}-{:s})'.format(self.label, state.label)
        result = STATE(lambda x: 1, label)
        result._indicators = [*_cp.deepcopy(self._indicators),
                              *_cp.deepcopy(state._indicators)]
        result._expression = '({} and not {})'.format(self._expression,
                                                      state._expression)
        return result

    def __neg__(self) -> 'STATE':
        """
        Return the outcome of the state.
        """
        label = '(!{:s})'.format(self.label)
        result = STATE(lambda x: 1, label)
        result._indicators = _cp.deepcopy(self._indicators)
        result._expression = '(not {})'.format(self._expression)
        return result

    def __call__(self, x: list) -> bool:
        """
        Call the expression.
        """
        flags = [indicator(x) for indicator in self._indicators]
        return bool(eval(self._expression.format(*flags)))

    def to_string(self) -> str:
        """
        Save a state to string.
        """
        expr = _bs.b64encode(self._expression.encode('UTF-8')).decode('UTF-8')
        label = _bs.b64encode(self.label.encode('UTF-8')).decode('UTF-8')
        result = label + '|' + expr
        for indicator in self._indicators:
            code = _bs.b64encode(_ma.dumps(indicator.__code__))
            result += '|' + code.decode('UTF-8')
        return result

    @staticmethod
    def from_string(s: str) -> 'STATE':
        """
        Set up a state from string.

        Parameters
        ----------
        s : str
            The saved state string.
        """
        tokens = s.split('|')
        label, expr = tokens[:2]
        expr = _bs.b64decode(expr.encode('UTF-8')).decode('UTF-8')
        label = _bs.b64decode(label.encode('UTF-8')).decode('UTF-8')
        indicators = []
        for i in tokens[2:]:
            code = _ma.loads(_bs.b64decode(i.encode('UTF-8')))
            indicators.append(_ts.FunctionType(code, globals(), 'indicator'))
        result = STATE(lambda x: 1, label)
        result._indicators = indicators
        result._expression = expr
        return result

    @property
    def _indicators(self) -> list:
        """
        The indicator functions of the state. For internal usages only.
        """
        return self.__indicators

    @_indicators.setter
    def _indicators(self, v) -> list:
        """
        Set the indicator functions of the state. For internal usages only.
        """
        self.__indicators = v

    @property
    def _expression(self) -> str:
        """
        Expression of this state. For internal usages only.
        """
        return self.__expression

    @_expression.setter
    def _expression(self, v) -> str:
        """
        Set expression of this state. For internal usages only.
        """
        self.__expression = v

    @property
    def label(self) -> str:
        """
        Label of this state.
        """
        return self.__label

    @label.setter
    def label(self, v: str) -> None:
        """
        Label of this state.
        """
        self.__label = v


def in_sequence(states: list,
                cv_values: list,
                allow_roaming: bool = False) -> bool:
    """
    Check if the given CV values pass a sequence of states.

    Parameters
    ----------
    states : List(somd.apps.path_sampling.utils.state.STATE)
        The states.
    cv_values : numpy.ndarray
        The CV values.
    allow_roaming : bool
        If allow the CV values to roam between intermediate states.
    """
    latest_state = -1
    n_states = len(states)
    results = [False] * len(states)
    for cv in cv_values:
        flags = [s(cv) for s in states]
        count = flags.count(True)
        if (count == 1):
            i = flags.index(True)
            if (i >= latest_state):
                latest_state = i
                results[i] = True
            elif (allow_roaming and latest_state != (n_states - 1) and i != 0):
                results[i] = True
            else:
                return False
        elif (count > 1):
            state_names = [states[i].label for i, f in enumerate(flags) if f]
            message = 'States {} are intersecting under CV values of {}!'
            raise RuntimeError(message.format(state_names, cv))
    return not (False in results)


def get_region_indices(state: STATE, cv_values: _np.ndarray) -> list:
    """
    Get indices of frames that are inside a state.

    Parameters
    ----------
    state : somd.apps.path_sampling.utils.state.STATE
        The states in the CV space.
    cv_values : numpy.ndarray
        The CV values.
    """
    frame_indices = []
    for i in range(0, len(cv_values)):
        if (state(cv_values[i])):
            frame_indices.append(i)
    return frame_indices


def get_reactive_regions(states: list, cv_values: _np.ndarray) -> list:
    """
    Get starting point and end point frame indices of reactive partial
    paths.

    Parameters
    ----------
    state : somd.apps.path_sampling.utils.state.STATE
        Two states in the CV space.
    cv_values : numpy.ndarray
        The CV values.
    """
    result = []
    indices_r = [(0, i) for i in get_region_indices(states[0], cv_values)]
    indices_p = [(1, i) for i in get_region_indices(states[1], cv_values)]
    if (indices_r == [] or indices_p == []):
        return result
    indices = [*indices_r, *indices_p]
    indices.sort(key=lambda x: x[1])
    for i in range(0, len(indices) - 1):
        if (indices[i][0] == 0 and indices[(i + 1)][0] == 1):
            result.append({'direction': 1, 'range':
                           (indices[i][1], indices[(i + 1)][1] + 1)})
        elif (indices[i][0] == 1 and indices[(i + 1)][0] == 0):
            result.append({'direction': -1, 'range':
                           (indices[i][1], indices[(i + 1)][1] + 1)})
    return result
