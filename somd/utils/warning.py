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
The warning method.
"""

import warnings as _w
import datetime as _dt


def _formatwarning(*args, **kwargs) -> str:
    """
    Function to format a warning the standard way.
    """
    result = '[SOMD] [{:s}] [{:s}|{:d}]: {}\n'
    file_path = args[2].split('somd/')[-1]
    time_string = _dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return result.format(time_string, file_path, args[3], args[0])


def warn(*args, **kwargs):
    """
    Override of the `warnings.warn` function.
    """
    formatwarning_tmp = _w.formatwarning
    _w.formatwarning = _formatwarning
    _w.warn(*args, stacklevel=2, **kwargs)
    _w.formatwarning = formatwarning_tmp
