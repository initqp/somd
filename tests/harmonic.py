import somd
import numpy as _np


class HARMONICPOTENTIAL(somd.core.potential_base.POTENTIAL):
    def update(self, system: somd.core.systems.MDSYSTEM) -> None:
        d_1 = system.positions[0, :] - _np.array([0, 0, 0], dtype=_np.double)
        d_2 = system.positions[1, :] - _np.array([0, 1, 0], dtype=_np.double)
        d_3 = system.positions[2, :] - _np.array([1, 1, 0], dtype=_np.double)
        d_4 = system.positions[3, :] - _np.array([1, 0, 0], dtype=_np.double)
        self.energy_potential[0] = 25 * (
            (d_1[0] ** 2 + d_1[0] ** 2 + d_1[0] ** 2)
            + (d_2[0] ** 2 + (d_2[0] - 1) ** 2 + d_2[0] ** 2)
            + ((d_3[0] - 1) ** 2 + (d_3[0] - 1) ** 2 + d_3[0] ** 2)
            + ((d_4[0] - 1) ** 2 + d_4[0] ** 2 + d_4[0] ** 2)
        )
        self.forces[0, :] = -50 * d_1
        self.forces[1, :] = -50 * d_2
        self.forces[2, :] = -50 * d_3
        self.forces[3, :] = -50 * d_4


def get_harmonic_system():
    s = somd.core.systems.MDSYSTEM(4, 'TEST')
    s.masses[:] = 1.0
    s.box[:] = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    for i in range(0, 4):
        s.atomic_symbols.append('X')
    s.positions[0, :] = [0, 0, 0]
    s.positions[1, :] = [0, 1, 0]
    s.positions[2, :] = [1, 1, 0]
    s.positions[3, :] = [1, 0, 0]
    s.velocities[0, :] = [1.0, 0.0, -0.5]
    s.velocities[1, :] = [0.0, 1.0, 0.5]
    s.velocities[2, :] = [0.5, 1.0, 0.0]
    s.velocities[3, :] = [1.0, 0.5, 1.0]
    s.potentials.append(HARMONICPOTENTIAL([0, 1, 2, 3]))
    s.groups.create_from_dict({'atom_list': [0, 1, 2, 3]})
    return s
