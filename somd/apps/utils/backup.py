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
Helper functions for backing up files.
"""

import os as _os
import re as _re
import typing as _tp

__all__ = ['back_up']


def __find_max_backup_number(base_name: str, file_list: _tp.List[str]) -> int:
    """
    Find maximum backup number of the file with a given base name.

    Parameters
    ----------
    base_name: str
        The base name of the files.
    file_list: List[str]
        Names of the files under the directory.
    """
    maximum = -2
    for fn in file_list:
        tmp = _re.split(base_name, fn)
        if len(tmp) > 1 and tmp[0] != '' and tmp[1] == '':
            # backup files.
            prefix = _re.split(r'\.', tmp[0])
            if len(prefix) == 3 and prefix[0] == 'bck' and prefix[1].isdigit():
                maximum = max(maximum, int(prefix[1]))
        elif len(tmp) > 1 and tmp[0] == '' and tmp[1] == '':
            # original file.
            if maximum == -2:
                maximum = -1
    return maximum


def back_up(file_name: str) -> None:
    """
    Back up a given file according to its backup history.

    Parameters
    ----------
    file_name : str
        Name of the file.
    """
    directory = _os.path.dirname(file_name)
    if directory == '':
        directory = '.'
    file_list = _os.listdir(directory)
    base_name = _os.path.basename(file_name)
    n = __find_max_backup_number(base_name, file_list)
    if (_os.path.exists(file_name) or _os.path.islink(file_name)) and n != -2:
        file_name_new = directory + '/' + 'bck.' + str(n + 1) + '.' + base_name
        _os.rename(file_name, file_name_new)
