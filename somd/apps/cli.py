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

import sys as _sys
import argparse as _ar
from somd import apps as _mdapp
from somd import _version

"""
A simple command line interface.
"""

__all__ = ['main']


def parse_args() -> _ar.Namespace:
    """
    Parse the command line arguments and perform some validation on the
    arguments.
    """
    description = 'SOMD: A SIESTA oriented Molecular Dynamics package.'
    parser = _ar.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', type=str,
                        help='The input TOML file.')
    parser.add_argument('-v', '--version', required=False,
                        action='store_true', help='Show SOMD vresion.')
    if (len(_sys.argv) == 1):
        parser.print_help()
        exit()
    return parser.parse_args()


def main() -> None:
    """
    The entry point of SOMD.
    """
    args = parse_args()
    if (args.version):
        print('SOMD version: ' + _version.get_versions()['version'])
        exit()
    if (args.input is not None):
        task = _mdapp.parser.TOMLPARSER(args.input)
        task.run()
