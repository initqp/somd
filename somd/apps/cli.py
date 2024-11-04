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
import numpy as _np
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
    sub_parsers = parser.add_subparsers(
        dest='mode', help='the execution mode (default: run)'
    )
    parser_run = sub_parsers.add_parser(
        'run',
        help='run a simulation from a TOML file',
        formatter_class=_ar.ArgumentDefaultsHelpFormatter
    )
    parser_select = sub_parsers.add_parser(
        'select',
        help='select structures for active learning',
        formatter_class=_ar.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        help='show SOMD vresion',
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        help=_ar.SUPPRESS  # override the default mode to `run`
    )
    parser_run.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='the input TOML file'
    )
    parser_select.add_argument(
        '-t',
        '--trajectories',
        type=str,
        nargs='*',
        required=True,
        help='paths to the input HDF5 trajectories, with forces recorded'
    )
    parser_select.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        help='the output trajectory file'
    )
    parser_select.add_argument(
        '--msd_f_min',
        type=float,
        default=50.0,
        help=(
            'lower boundary of the force MSD when identifying new structures.'
            + ' In unit of (kJ/mol/nm).'
        )
    )
    parser_select.add_argument(
        '--msd_f_max',
        type=float,
        default=250.0,
        help=(
            'upper boundary of the force MSD when identifying new structures.'
            + ' In unit of (kJ/mol/nm).'
        )
    )
    parser_select.add_argument(
        '--n_max',
        type=int,
        default=_np.inf,
        help='max number of structures to select'
    )
    parser_select.add_argument(
        '--log_file',
        type=str,
        default='somd_selection.json',
        help='name of the log file contains the selection information'
    )
    parser_select.add_argument(
        '--shuffle',
        action='store_true',
        default=True,
        help='if perform randomly selection'
    )
    if len(_sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def check_args(args: _ar.Namespace) -> _ar.Namespace:
    """
    Perform sanity check of the command line args.
    """
    if args.mode is None:
        args.mode = 'run'
    return args


def main() -> None:
    """
    The entry point of SOMD.
    """
    args = check_args(parse_args())
    if args.version:
        print('SOMD version: ' + _version.get_versions()['version'])
        exit()
    if args.mode == 'run':
        task = _mdapp.parser.TOMLPARSER(args.input)
        task.run()
    elif args.mode == 'select':
        selector = _mdapp.select.STRUCTURESELECTOR(args.trajectories)
        indices = selector.select(
            args.n_max,
            args.msd_f_min,
            args.msd_f_max,
            args.shuffle,
            args.log_file,
        )
        selector.write(args.output, indices)
