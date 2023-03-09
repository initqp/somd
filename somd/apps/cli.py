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
