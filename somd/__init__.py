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
SOMD is an ab-initio molecular dynamics (AIMD) package designed for the SIESTA
(https://departments.icmab.es/leem/siesta/) code. SOMD also contains
subroutines to automatically train the neuroevolution potential (NEP) using the
so-called active-learning methodology.
"""

from somd import constants
from somd import core
from somd import potentials
from somd import apps

from . import _version
__version__ = _version.get_versions()['version']

import warnings as _w
_w.simplefilter('always', UserWarning)
