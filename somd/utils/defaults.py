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
Default values in SOMD
"""


# Length of the Nose-Hoover chains.
NHCLENGTH: int = 6

# Number of the RESPA loops of the Nose-Hoover chains.
NHCNRESPA: int = 4

# The K_r value of the geodesic Langevin integrators.
GEODESICKR: int = 5

# Tolerance for determining lattice parameters.
LATTICETOL: float = 1e-7

# If simultaneously update potentials.
SIMUUPDATE: bool = False

# Timeout for any SIESTA operations (second).
SIESTATIMEOUT: int = 540

# Global verbose level
VERBOSE: bool = False
