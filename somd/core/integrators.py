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
The integrators.
"""

import numpy as _np
import types as _ts
import typing as _tp
from somd import utils as _mdutils
from ._lib import NHCHAINS as _NHCHAINS
from .systems import MDSYSTEM as _MDSYSTEM

__all__ = [
    'INTEGRATOR',
    'vv_integrator',
    'cs4_integrator',
    'nhc_integrator',
    'baoab_integrator',
    'obabo_integrator',
    'gbaoab_integrator',
    'gobabo_integrator',
]


class INTEGRATOR(object):
    """
    The splitting integrator [1].

    Parameters
    ----------
    timestep : float
        The timestep length of the integrator. In unit of ps.
    splitting : List[dict]
        A list of dict that describe the splitting scheme. Valid keys
        of the dictionary are:
        - 'operators' : List[str]
            Name of the operators. Valid values are:
            - 'V'  : Advance the velocities
            - 'R'  : Advance the positions
            - 'N'  : Advance both the velocities and positions of a
                     Nose-Hoover Chain [2]
            - 'O'  : Advance the Ornstein-Uhlenbeck process [3]
            - 'B'  : Advance the Bussi-Donadio-Parrinello process [4]
            - 'Cr' : Apply the RATTLE method for positions  [5]
            - 'Cv' : Apply the RATTLE method for velocities [5]
        - 'timesteps' : List[float]
            Operator-wise scaling factors of the timestep. The length of this
            list must be the same as the length of the 'operators' list. In
            most cases the scaling factor of the timesteps are automatically
            selected. Only when you need a custom splitting scheme (e.g., some
            high-order splitting integrators), this key should be set.
        - 'repeating' : float
            Repeating times of the operator sequence definded by the key
            'operators'.
    temperatures : List[float]
        Temperatures of the thermostats (if present). In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats (if present). In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups (if thermostats present). Note that each
        thermostat only thermalizes one atomic group. Thus, lengths of the
        'temperature', 'relaxation_times' and 'thermo_groups' parameters
        must be the same.

    References
    ----------
    .. [1] Leimkuhler, Ben, and Charles Matthews. "Molecular dynamics."
        Interdisciplinary applied mathematics 36 (2015).
    .. [2] Martyna, Glenn J., et al. "Explicit reversible integrators for
        extended systems dynamics." Molecular Physics 87.5 (1996): 1117-1157.
    .. [3] Leimkuhler, Benedict, and Charles Matthews. "Rational construction
        of stochastic numerical methods for molecular sampling." Applied
        Mathematics Research eXpress 2013.1 (2013): 34-56.
    .. [4] Bussi, Giovanni, Davide Donadio, and Michele Parrinello.
        "Canonical sampling through velocity rescaling." The Journal of
        chemical physics 126.1 (2007).
    .. [5] Andersen, Hans C. "Rattle: A velocity version of the shake algorithm
        for molecular dynamics calculations." Journal of computational Physics
        52.1 (1983): 24-34.
    """

    # Uppercase valid operator names.
    _valid_operators = ['V', 'R', 'B', 'O', 'N', 'CR', 'CV']

    def __init__(
        self,
        timestep: float,
        splitting: _tp.List[dict] = [{'operators': ['V', 'R', 'V']}],
        temperatures: _tp.List[float] = [300.0],
        relaxation_times: _tp.List[float] = [0.1],
        thermo_groups: _tp.List[int] = [0],
    ) -> None:
        """
        Create an integrator based on given splitting method.
        """
        self.__step = 0
        self.__system = None
        self.__is_nve = True
        self.__timestep = timestep
        self.__is_stochastic = False
        self.__energy_effective = 0.0
        self.__temperatures = temperatures
        self.__thermo_groups = thermo_groups
        self.__relaxation_times = relaxation_times
        self.__cmm_remover_groups = []
        for i in range(0, len(splitting)):
            l = splitting[i]['operators']
            splitting[i]['operators'] = [o.strip().upper() for o in l]
        self.__splitting = splitting
        self.__compile()

    def __copy__(self) -> 'INTEGRATOR':
        """
        Clone this integrator.
        """
        return self.copy()

    def __compile(self) -> None:
        """
        Compile splitting schemes to codes.
        """
        self.__determine_timeteps()
        self.__combine_splitting_dicts()
        self.__determine_ensemble()
        self.__insert_F_operators()
        self.__splitting_dict_to_codes()

    def __determine_ensemble(self) -> None:
        """
        Determine if the integrator is thermalized.
        """
        if 'B' in self.__splitting_whole['operators']:
            self.__is_nve = False
            self.__is_stochastic = True
        if 'O' in self.__splitting_whole['operators']:
            self.__is_nve = False
            self.__is_stochastic = True
        if 'N' in self.__splitting_whole['operators']:
            self.__init_nhchains()
            self.__is_nve = False
        if self.__is_nve:
            del self.__temperatures
            del self.__thermo_groups
            del self.__relaxation_times
        else:
            if len(self.__temperatures) == 0:
                message = (
                    'The temperatures parameter must be '
                    + 'definded for a thermalized integrator!'
                )
                raise RuntimeError(message)
            if len(self.__relaxation_times) == 0:
                message = (
                    'The relaxation_times parameter must be '
                    + 'definded for a thermalized integrator!'
                )
                raise RuntimeError(message)
            if len(self.__thermo_groups) == 0:
                message = (
                    'The thermo_groups parameter must be '
                    + 'definded for a thermalized integrator!'
                )
                raise RuntimeError(message)
            if (len(self.__temperatures) != len(self.__thermo_groups)) or (
                len(self.__temperatures) != len(self.__relaxation_times)
            ):
                message = (
                    'Number of temperatures, thermo groups '
                    + 'and relaxation times must be the same!'
                )
                raise IndexError(message)

    def __init_nhchains(self) -> None:
        """
        Initialize the Nose-Hoover chains.
        """
        self.__nhchains = []
        for i in range(0, len(self.__thermo_groups)):
            self.__nhchains.append(
                _NHCHAINS(
                    self.__temperatures[i],
                    self.__relaxation_times[i],
                    _mdutils.defaults.NHCLENGTH,
                    1,
                    _mdutils.defaults.NHCNRESPA,
                )
            )

    def __reset_nhchains_parameters(self) -> None:
        """
        Reset parameters of the Nose-Hoover chains.
        """
        for i in range(0, len(self.__thermo_groups)):
            g = self.__thermo_groups[i]
            self.__nhchains[i].temperature = self.__temperatures[i]
            self.__nhchains[i].n_dof = self.__system.groups[g].n_dof
            self.__nhchains[i].tau = self.__relaxation_times[i]

    def __determine_timeteps(self) -> None:
        """
        Determine timestep of each operator.
        """
        v = self._valid_operators
        # Determine timesteps of each block.
        op = [b['operators'] for b in self.__splitting]
        ts = [1.0 / op.count(b) for b in op]
        # Determine timesteps of each operator.
        for i in range(0, len(self.__splitting)):
            s = self.__splitting[i]
            if 'timesteps' not in s:
                invalid_ops = [o for o in op[i] if (o not in v)]
                if len(invalid_ops) != 0:
                    message = 'Unknown operators: "{}"'.format(invalid_ops)
                    raise RuntimeError(message)
                s['timesteps'] = [ts[i] / op[i].count(o) for o in op[i]]
            else:
                if len(s['timesteps']) != len(s['operators']):
                    message = (
                        'Numbers of timesteps and operators '
                        + 'mismatch in splitting scheme "{}"!'
                    )
                    raise IndexError(message.format(self.splitting))

    def __combine_splitting_dicts(self) -> None:
        """
        Combine multiple splitting dicts.
        """
        self.__splitting_whole = dict()
        self.__splitting_whole['operators'] = []
        self.__splitting_whole['timesteps'] = []
        self.__splitting_whole['repeating'] = 1
        for i in range(0, len(self.__splitting)):
            s = self.__splitting[i]
            r = s.get('repeating', 1)
            self.__splitting_whole['operators'] += s['operators'] * r
            ts = [t / r for t in s['timesteps']]
            self.__splitting_whole['timesteps'] += ts * r

    def __insert_F_operators(self) -> None:
        """
        Determine when to update forces.
        """
        op = self.__splitting_whole['operators']
        # Check if we have the R operators.
        if 'R' not in op:
            message = (
                'Splitting scheme "{}" without a position update is invalid!'
            )
            raise RuntimeError(message.format(self.__splitting))
        if 'V' not in op:
            message = (
                'Splitting scheme "{}" without a velocity update is invalid!'
            )
            raise RuntimeError(message.format(self.__splitting))
        # Get the 'core' splitting.
        index_core = [i for i in range(0, len(op)) if op[i] in ['R', 'V']]
        splitting_core = [o for o in op if o in ['R', 'V']]
        # Wrap the first core operator to the last to complete the operator
        # sequence.
        index_core.append(index_core[-1] + 1)
        splitting_core.append(splitting_core[0])
        # Determine where to update the forces.
        shifting = 0
        for i in range(0, (len(splitting_core) - 1)):
            if splitting_core[i] == 'R' and splitting_core[(i + 1)] == 'V':
                op_index = index_core[(i + 1)] + shifting
                self.__splitting_whole['operators'].insert(op_index, 'F')
                self.__splitting_whole['timesteps'].insert(op_index, 1.0)
                shifting += 1
        if 'F' not in self.__splitting_whole['operators']:
            message = 'Failed to parse splitting scheme: "{}"'
            raise RuntimeError(message.format(self.__splitting))

    def __splitting_dict_to_codes(self) -> None:
        """
        Convert a splitting method to codes.
        """
        t = []
        count = 0
        scope = {}
        s = 'def propagate(self) -> None:\n'
        s += '    """\n    Propagate the system by one timestep.\n    """\n'
        for i, j in zip(
            self.__splitting_whole['operators'],
            self.__splitting_whole['timesteps'],
        ):
            s += '    self._operator_{}({})\n'.format(i, count)
            if i == 'O':
                t.append(_np.abs(self.__timestep * j))
            else:
                t.append(self.__timestep * j)
            count += 1
        s += '    self.step += 1\n'
        self.__timesteps = _np.array(t, dtype=_np.double)
        exec(s, scope)
        self.propagate = _ts.MethodType(scope['propagate'], self)

    def __remove_com_motions(self) -> None:
        """
        Remove the COM motions.
        """
        for i in self.__cmm_remover_groups:
            g = self.__system.groups[i]
            self.energy_effective += g.energy_kinetic
            g.remove_com_motion(self.__is_nve)
            self.energy_effective -= g.energy_kinetic

    def _randomize_nhchains_states(self) -> None:
        """
        Randomize momentums of the Nose-Hoover chains.
        """
        for n in self._nhchains:
            factors = _np.sqrt(
                n.temperature
                * _mdutils.constants.BOLTZCONST
                * _np.array(n.masses)
            )
            n.momentums = _mdutils.rng.standard_normal(n.length) * factors

    def _operator_V(self, dt_index: int) -> None:
        """
        Advance the velocities and apply the COM motion removers.
        """
        dt = self.__timesteps[dt_index]
        self.__system.velocities[:, :] += (
            self.__system.forces / self.__system.masses
        ).dot(dt)
        self.__remove_com_motions()

    def _operator_R(self, dt_index: int) -> None:
        """
        Advance the positions.
        """
        dt = self.__timesteps[dt_index]
        self.__system.positions[:, :] += self.__system.velocities.dot(dt)

    def _operator_CR(self, dt_index: int) -> None:
        """
        Perform RATTLE position constraints.
        """
        dt = self.__timesteps[dt_index]
        self.__system.groups.constraints.rattle_constrain_q(dt)

    def _operator_CV(self, dt_index: int) -> None:
        """
        Perform RATTLE velocity constraints.
        """
        dt = self.__timesteps[dt_index]
        self.__system.groups.constraints.rattle_constrain_p(dt)
        self.__remove_com_motions()

    def _operator_F(self, dt_index: int) -> None:
        """
        Update the atomic forces (not a real operator).
        """
        self.__system.update_potentials()

    def _operator_O(self, dt_index: int) -> None:
        """
        Propagate the Ornstein-Uhlenbeck process.
        """
        dt = self.__timesteps[dt_index]
        for i in range(0, len(self.__thermo_groups)):
            g = self.__system.groups[self.__thermo_groups[i]]
            self.__energy_effective += g.energy_kinetic
            c_1 = _np.exp(-1.0 * dt / self.__relaxation_times[i])
            c_2 = _np.sqrt(
                (1.0 - c_1 * c_1)
                * self.__temperatures[i]
                * _mdutils.constants.BOLTZCONST
            )
            g.velocities *= c_1
            g.velocities += (
                _mdutils.rng.standard_normal((g.n_atoms, 3))
                * c_2
                / _np.sqrt(g.masses)
            )
            self.__energy_effective -= g.energy_kinetic
        self.__remove_com_motions()

    def _operator_N(self, dt_index: int) -> None:
        """
        Propagate the Nose-Hoover Chains.
        """
        dt = self.__timesteps[dt_index]
        for i in range(0, len(self.__thermo_groups)):
            g = self.__system.groups[self.__thermo_groups[i]]
            energy_kinetic = g.energy_kinetic
            self.__energy_effective += energy_kinetic
            g.velocities *= self.__nhchains[i].propagate(energy_kinetic, dt)
            self.__energy_effective -= g.energy_kinetic
        self.__remove_com_motions()

    def _operator_B(self, dt_index: int) -> None:
        """
        Propagate the Bussi-Donadio-Parrinello process.
        """
        # fmt: off
        dt = self.__timesteps[dt_index]
        for i in range(0, len(self.__thermo_groups)):
            g = self.__system.groups[self.__thermo_groups[i]]
            energy_kinetic = g.energy_kinetic
            self.__energy_effective += energy_kinetic
            energy_kinetic_target = (
                self.__temperatures[i] * _mdutils.constants.BOLTZCONST
                * g.n_dof * 0.5
            )
            r_1 = _mdutils.rng.standard_normal()
            if (g.n_dof - 1) % 2 == 0:
                r_2 = 2.0 * _mdutils.rng.standard_gamma((g.n_dof - 1) / 2)
            else:
                r_2 = 2.0 * _mdutils.rng.standard_gamma((g.n_dof - 2) / 2)
                r_2 += _mdutils.rng.standard_normal() ** 2
            c_1 = _np.exp(-1.0 * dt / self.__relaxation_times[i])
            term_1 = (
                energy_kinetic_target * (r_1**2 + r_2) / g.n_dof
                - energy_kinetic
            ) * (1.0 - c_1)
            term_2 = _np.sqrt(
                energy_kinetic * energy_kinetic_target / g.n_dof
                * (1.0 - c_1) * c_1
            ) * 2.0 * r_1
            g.velocities *= _np.sqrt(1 + (term_1 + term_2) / energy_kinetic)
            self.__energy_effective -= g.energy_kinetic
        self.__remove_com_motions()
        # fmt: on

    def bind_system(self, system: _MDSYSTEM) -> None:
        """
        Bind this integrator to a simulated system.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            the simulated system.
        """
        if self.__system is not None:
            message = 'Binding integrator from system "{}" to system "{}".'
            message = message.format(self.__system._label, system._label)
            _mdutils.warning.warn(message)
        if not self.__is_nve:
            if len(set(self.__thermo_groups)) != len(self.__thermo_groups):
                message = 'Duplicate indices in thermalized groups "{}"!'
                raise IndexError(message.format(self.__thermo_groups))
            if len(self.__thermo_groups) > len(system.groups):
                message = 'Too many thermalized groups!'
                raise IndexError(message)
            if max(self.__thermo_groups) >= len(system.groups):
                message = 'Invalid thermalized group: {:d}!'
                raise IndexError(message.format(max(self.__thermo_groups)))
            # Thermo groups
            overlaps = [
                (len(system.groups[g_i] & system.groups[g_j]) > 0)
                for g_i in self.__thermo_groups
                for g_j in self.__thermo_groups
                if g_i != g_j
            ]
            if True in overlaps:
                message = 'Overlaps found between thermalized groups!'
                raise RuntimeError(message)
        self.__system = system
        if 'N' in self.__splitting_whole['operators']:
            self.__reset_nhchains_parameters()
        # Setup COM motion removers.
        self.__cmm_remover_groups = [
            i for i, g in enumerate(system.groups) if not g.has_translations
        ]

    def drop_system(self) -> None:
        """
        Unbind this integrator from a simulated system.
        """
        self.__system = None
        self.__cmm_remover_groups = []
        if 'N' in self.__splitting_whole['operators']:
            self.__nhchains = []

    def copy(self) -> 'INTEGRATOR':
        """
        Copy the integrator.
        """
        # fmt: off
        if self._is_nve:
            integrator = INTEGRATOR(
                self.timestep,
                self.splitting,
            )
        else:
            integrator = INTEGRATOR(
                self.timestep,
                self.splitting,
                self.temperatures,
                self.relaxation_times,
                self.thermo_groups,
            )
        if 'N' in self.__splitting_whole['operators']:
            # We do not know the length of the old chains. Thus we copy them
            # directly.
            integrator._nhchains.clear()
            for nhc in self._nhchains:
                integrator._nhchains.append(nhc.copy())
        # fmt: on
        return integrator

    def summary(self) -> str:
        """
        Show information about this integrator.
        """
        splitting = self.__splitting_whole['operators']
        splitting_t = self.__splitting_whole['timesteps']

        result = 'INTEGRATOR\n'
        result += '┣━ timestep: {}\n'.format(self.__timestep)
        result += '┣━ operators: {}\n'.format(splitting)
        if _mdutils.defaults.VERBOSE:
            result += '┣━ operator_timesteps: {}\n'.format(splitting_t)
        result += '┣━ is_stochastic: {}\n'.format(self._is_stochastic)
        result += '┣━ is_nve: {}\n'.format(self._is_nve)
        if not self._is_nve:
            result += '┣━ temperatures: {}\n'.format(self.temperatures)
            result += '┣━ thermo_groups: {}\n'.format(self.thermo_groups)
        if _mdutils.defaults.VERBOSE:
            if 'N' in self.__splitting_whole['operators']:
                for i, nhc in enumerate(self._nhchains):
                    result += '┣━ NHC {}: '.format(i)
                    result += 'temperature {}, '.format(nhc.temperature)
                    result += 'tau: {}, '.format(nhc.tau)
                    result += 'length: {}, '.format(nhc.length)
                    result += 'n_dof: {}, '.format(nhc.n_dof)
                    result += 'n_respa: {}\n'.format(nhc.n_respa)
        result += '┗━ END'
        return result

    @property
    def system(self) -> _MDSYSTEM:
        """
        Reference to the bound system.
        """
        return self.__system

    @property
    def _is_nve(self) -> bool:
        """
        Is this a NVE integrator?
        """
        return self.__is_nve

    @property
    def _is_stochastic(self) -> bool:
        """
        Is this a stochastic integrator?
        """
        return self.__is_stochastic

    @property
    def _nhchains(self) -> list:
        """
        The Nose-Hoover chains bound to the integrator.
        """
        if 'N' in self.__splitting_whole['operators']:
            return self.__nhchains
        else:
            return []

    @property
    def step(self) -> int:
        """
        Current time step.
        """
        return self.__step

    @step.setter
    def step(self, s: int) -> None:
        """
        Set current time step.
        """
        self.__step = s

    @property
    def splitting(self) -> _tp.Dict[str, _tp.Any]:
        """
        The splitting scheme of the integrator.
        """
        return self.__splitting

    @property
    def splitting_whole(self) -> _tp.Dict[str, _tp.Any]:
        """
        The expanded splitting scheme of the integrator.
        """
        return self.__splitting_whole

    @property
    def timestep(self) -> float:
        """
        The timestep of the integrator.
        """
        return self.__timestep

    @timestep.setter
    def timestep(self, t: float) -> None:
        """
        Set the timestep of the integrator.
        """
        self.__timestep = t
        # Avoid re-determine timesteps here. Since we do not know if these
        # timesteps were set by the user or by the program.
        self.__splitting_dict_to_codes()

    @property
    def temperatures(self) -> float:
        """
        Thermostat temperatures.
        """
        if not self.__is_nve:
            return self.__temperatures
        else:
            message = 'Can not get temperatures of NVE integrators!'
            raise AttributeError(message)

    @temperatures.setter
    def temperatures(self, t: _tp.List[float]) -> None:
        """
        Set the temperatures of the thermostats.
        """
        if not self.__is_nve:
            self.__temperatures[:] = t[:]
            if 'N' in self.__splitting_whole['operators']:
                self.__reset_nhchains_parameters()
        else:
            message = 'Can not set temperatures of NVE integrators!'
            raise AttributeError(message)

    @property
    def relaxation_times(self) -> float:
        """
        Thermostat relaxation times.
        """
        if not self.__is_nve:
            return self.__relaxation_times
        else:
            message = (
                'Can not get thermostat relaxation times of NVE integrators!'
            )
            raise AttributeError(message)

    @relaxation_times.setter
    def relaxation_times(self, t: _tp.List[float]) -> None:
        """
        Set the relaxation times of the thermostats.
        """
        if not self.__is_nve:
            self.__relaxation_times[:] = t[:]
            if 'N' in self.__splitting_whole['operators']:
                self.__reset_nhchains_parameters()
        else:
            message = (
                'Can not set thermostat relaxation times of NVE integrators!'
            )
            raise AttributeError(message)

    @property
    def energy_effective(self) -> float:
        """
        Effective energy (the correction part only) of the propagated system.
        In units of kj/mol.
        """
        return self.__energy_effective

    @energy_effective.setter
    def energy_effective(self, e: float) -> None:
        """
        Set the effective energy (the correction part only) of the propagated
        system. In units of kj/mol.
        """
        self.__energy_effective = e

    @property
    def thermo_groups(self) -> float:
        """
        Thermalized atomic groups.
        """
        if not self.__is_nve:
            # Do not change thermo groups here.
            return self.__thermo_groups.copy()
        else:
            message = (
                'Can not get thermalized atomic groups of NVE integrators!'
            )
            raise AttributeError(message)


def vv_integrator(timestep: float) -> INTEGRATOR:
    """
    The velocity-Verlet integrator [1].

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.

    References
    ----------
    .. [1] Swope, William C., et al. "A computer simulation method for the
           calculation of equilibrium constants for the formation of physical
           clusters of molecules: Application to small water clusters." The
           Journal of chemical physics 76.1 (1982): 637-649.
    """
    return INTEGRATOR(timestep, [{'operators': ['V', 'CR', 'R', 'V', 'CV']}])


def cs4_integrator(timestep: float) -> INTEGRATOR:
    """
    Calvo and Sanz-Serna's fourth-order integrator [1].

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.

    References
    ----------
    .. [1] Calvo, Mari Paz, and Jesus Maria Sanz-Serna. "The development of
           variable-step symplectic integrators, with application to the
           two-body problem." SIAM Journal on Scientific Computing 14.4 (1993):
           936-952.
    """
    return INTEGRATOR(
        timestep,
        [
            {
                'operators': ['V', 'Cv', 'Cr', 'R', 'Cv'],
                'timesteps': [
                    0.06175885813563,
                    0.06175885813563,
                    0.20517766154229,
                    0.20517766154229,
                    0.20517766154229,
                ],
            },
            {
                'operators': ['V', 'Cv', 'Cr', 'R', 'Cv'],
                'timesteps': [
                    0.33897802655364,
                    0.33897802655364,
                    0.40302128160421,
                    0.40302128160421,
                    0.40302128160421,
                ],
            },
            {
                'operators': ['V', 'Cv', 'Cr', 'R', 'Cv'],
                'timesteps': [
                    0.61479130717558,
                    0.61479130717558,
                    -0.12092087633891,
                    -0.12092087633891,
                    -0.12092087633891,
                ],
            },
            {
                'operators': ['V', 'Cv', 'Cr', 'R', 'Cv'],
                'timesteps': [
                    -0.14054801465937,
                    -0.14054801465937,
                    0.51272193319241,
                    0.51272193319241,
                    0.51272193319241,
                ],
            },
            {
                'operators': ['V', 'Cv'],
                'timesteps': [0.12501982279453, 0.12501982279453],
            },
        ],
    )


def baoab_integrator(
    timestep: float,
    temperatures: _tp.List[float] = [300.0],
    relaxation_times: _tp.List[float] = [0.1],
    thermo_groups: _tp.List[int] = [0],
) -> INTEGRATOR:
    """
    The Langevin integrator with a BAOAB splitting [1].

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.
    temperatures : List[float]
        Temperatures of the thermostats. In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats. In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups. Note that each thermostat only thermalizes
        one atomic group. Thus, lengths of the 'temperature',
        'relaxation_times' and 'thermo_groups' parameters must be the same.

    References
    ----------
    .. [1] Leimkuhler, Benedict, and Charles Matthews. "Rational construction
           of stochastic numerical methods for molecular sampling." Applied
           Mathematics Research eXpress 2013.1 (2013): 34-56.
    """
    return INTEGRATOR(
        timestep,
        [{'operators': ['V', 'R', 'O', 'R', 'V']}],
        temperatures,
        relaxation_times,
        thermo_groups,
    )


def obabo_integrator(
    timestep: float,
    temperatures: _tp.List[float] = [300.0],
    relaxation_times: _tp.List[float] = [0.1],
    thermo_groups: _tp.List[int] = [0],
) -> INTEGRATOR:
    """
    The Langevin integrator with an OBABO splitting [1].

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.
    temperatures : List[float]
        Temperatures of the thermostats. In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats. In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups. Note that each thermostat only thermalizes
        one atomic group. Thus, lengths of the 'temperature',
        'relaxation_times' and 'thermo_groups' parameters must be the same.

    References
    ----------
    .. [1] Bussi, Giovanni, and Michele Parrinello. "Accurate sampling using
           Langevin dynamics." Physical Review E 75.5 (2007): 056707.
    """
    return INTEGRATOR(
        timestep,
        [{'operators': ['O', 'V', 'R', 'V', 'O']}],
        temperatures,
        relaxation_times,
        thermo_groups,
    )


def gbaoab_integrator(
    timestep: float,
    temperatures: _tp.List[float] = [300.0],
    relaxation_times: _tp.List[float] = [0.1],
    thermo_groups: _tp.List[int] = [0],
) -> INTEGRATOR:
    """
    The geodesic Langevin integrator with a BAOAB splitting [1]. Note that this
    integrator should only be used when constraints are applied to the
    simulated system.

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.
    temperatures : List[float]
        Temperatures of the thermostats. In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats. In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups. Note that each thermostat only thermalizes
        one atomic group. Thus, lengths of the 'temperature',
        'relaxation_times' and 'thermo_groups' parameters must be the same.

    References
    ----------
    .. [1] Leimkuhler, Benedict, and Charles Matthews. "Efficient molecular
           dynamics using geodesic integration and solvent-solute splitting."
           Proceedings of the Royal Society A: Mathematical, Physical and
           Engineering Sciences 472.2189 (2016): 20160138.
    """
    return INTEGRATOR(
        timestep,
        [
            {'operators': ['V', 'Cv']},
            {
                'operators': ['Cr', 'R', 'Cv'],
                'repeating': _mdutils.defaults.GEODESICKR,
            },
            {'operators': ['O', 'Cv']},
            {
                'operators': ['Cr', 'R', 'Cv'],
                'repeating': _mdutils.defaults.GEODESICKR,
            },
            {'operators': ['V', 'Cv']},
        ],
        temperatures,
        relaxation_times,
        thermo_groups,
    )


def gobabo_integrator(
    timestep: float,
    temperatures: _tp.List[float] = [300.0],
    relaxation_times: _tp.List[float] = [0.1],
    thermo_groups: _tp.List[int] = [0],
) -> INTEGRATOR:
    """
    The geodesic Langevin integrator with an OBABO splitting [1]. Note that
    this integrator should only be used when constraints are applied to the
    simulated system.

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.
    temperatures : List[float]
        Temperatures of the thermostats. In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats. In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups. Note that each thermostat only thermalizes
        one atomic group. Thus, lengths of the 'temperature',
        'relaxation_times' and 'thermo_groups' parameters must be the same.

    References
    ----------
    .. [1] Leimkuhler, Benedict, and Charles Matthews. "Efficient molecular
           dynamics using geodesic integration and solvent-solute splitting."
           Proceedings of the Royal Society A: Mathematical, Physical and
           Engineering Sciences 472.2189 (2016): 20160138.
    """
    return INTEGRATOR(
        timestep,
        [
            {'operators': ['O', 'Cv']},
            {'operators': ['V', 'Cv']},
            {
                'operators': ['Cr', 'R', 'Cv'],
                'repeating': _mdutils.defaults.GEODESICKR,
            },
            {'operators': ['V', 'Cv']},
            {'operators': ['O', 'Cv']},
        ],
        temperatures,
        relaxation_times,
        thermo_groups,
    )


def nhc_integrator(
    timestep: float,
    temperatures: _tp.List[float] = [300.0],
    relaxation_times: _tp.List[float] = [0.1],
    thermo_groups: _tp.List[int] = [0],
) -> INTEGRATOR:
    """
    The Nose-Hoover Chains integrator [1].

    Parameter
    ---------
    timestep : float
        Timestep of the integrator. In unit of ps.
    temperatures : List[float]
        Temperatures of the thermostats. In unit of K.
    relaxation_times : List[float]
        Relaxation times of the thermostats. In unit of ps.
    thermo_groups : List[int]
        Thermalized atomic groups. Note that each thermostat only thermalizes
        one atomic group. Thus, lengths of the 'temperature',
        'relaxation_times' and 'thermo_groups' parameters must be the same.
        `_randomize_nhchains_states` method).

    References
    ----------
    .. [1] Martyna, Glenn J., et al. "Explicit reversible integrators for
           extended systems dynamics." Molecular Physics 87.5 (1996):
           1117-1157.
    """
    return INTEGRATOR(
        timestep,
        [{'operators': ['N', 'V', 'CR', 'R', 'V', 'CV', 'N']}],
        temperatures,
        relaxation_times,
        thermo_groups,
    )
