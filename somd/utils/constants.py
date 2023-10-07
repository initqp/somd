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
CODATA 2018 standard constants. Website:
https://physics.nist.gov/cuu/Constants/index.html
"""


# Hartree energy. In unit of (kJ/mol).
HARTREE: float = 2625.4978694234683

# Boltzmann Constant. In unit (kJ/mol/K).
BOLTZCONST: float = 8.31446261815324E-3

# Avogadro Number. In unit of (1).
AVOGACONST: float = 6.0221367E23

# Planck constant. In unit of (kJ/mol*ps).
PLANKCONST: float = 0.399031002270895

# Elementary charge. In unit of (C).
ELECTCONST: float = 1.602176634E-19

# Bohr radius. In unit of (nm).
BOHRRADIUS: float = 5.29177210903E-2

# Speed of light in vacuum. In unit of (nm/ps).
LIGHTSPEED: float = 299792.458

# Pressure unit. In unit of (kJ/mol/nm^3).
MEGAPASCAL: float = 0.60221367
