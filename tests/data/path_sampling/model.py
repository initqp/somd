
import somd
import numpy as _np


class PES(somd.core.potential_base.POTENTIAL):
    def __init__(self, atom_list: list, B: float = 3.0) -> None:
        super().__init__(atom_list)
        self._B = B

    def update(self, system: somd.core.systems.MDSYSTEM) -> None:
        x = system.positions[0, 0]
        y = system.positions[0, 1]
        self.energy_potential[0] = self._B * ((x ** 2 - 1) ** 2 + (x - y) ** 2)
        self.forces[0, 0] = self._B * (4 * x * (x ** 2 - 1) + 2 * (x - y)) * -1
        self.forces[0, 1] = self._B * (2 * (x - y))


class HELPER(somd.apps.utils.POSTSTEPOBJ):
    def update(self) -> None:
        self.integrator.system.positions[0, 2] = 0.0
        self.integrator.system.velocities[0, 2] = 0.0


def add_velocities_from_temperature(self, temperature: float) -> None:
    if (temperature == 0):
        return
    v = _np.random.randn(self.n_atoms, 3) * \
        _np.sqrt(temperature * somd.constants.CONSTANTS.BOLTZCONST /
                 self.masses)
    v[0, 2] = 0
    t = (_np.square(v) * self.masses).sum() / \
        somd.constants.CONSTANTS.BOLTZCONST / self.n_dof
    v *= _np.sqrt(temperature / t)
    self.velocities += v


somd.core.groups.ATOMGROUP.add_velocities_from_temperature = \
    add_velocities_from_temperature
somd.core.groups.ATOMGROUP.n_dof = 2
kb = somd.constants.CONSTANTS.BOLTZCONST
system = somd.core.systems.MDSYSTEM(1)
system.box[:] = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
system.masses[0] = 1.0
system.atomic_symbols.append('H')
system.groups.create_from_dict({'atom_list': [0]})

integrator = somd.core.integrators.obabo_integrator(0.002,
                                                    temperatures=[300],
                                                    relaxation_times=[0.05])
helper = HELPER(1)


def i1(x):
    import numpy as np
    d1 = np.sqrt((x[0] + 1) ** 2 + (x[1] + 1) ** 2)
    d2 = np.sqrt((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
    if (((x[0] ** 2 - 1) ** 2 + (x[0] - x[1]) ** 2) < 0.1):
        return bool(d1 < d2)
    else:
        return False


def i2(x):
    import numpy as np
    d1 = np.sqrt((x[0] + 1) ** 2 + (x[1] + 1) ** 2)
    d2 = np.sqrt((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
    if (((x[0] ** 2 - 1) ** 2 + (x[0] - x[1]) ** 2) < 0.1):
        return bool(d1 > d2)
    else:
        return False


def i3(x):
    return bool((x[1] > -x[0] - 0.5) * (x[1] < -x[0] + 0.5))
