import somd
import os as _os
import numpy as _np
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

somd.utils.rng = somd.utils._rng.LEGACYRNG(1)


def test_select():
    system = somd.core.systems.create_system_from_pdb(
        './data/active_learning/topo.pdb'
    )
    system.groups.create_from_dict({'atom_list': range(0, 8)})
    neps = [
        somd.potentials.NEP(range(system.n_atoms), f, system.atomic_symbols)
        for f in [
            './data/active_learning/nep.0.txt',
            './data/active_learning/nep.1.txt',
            './data/active_learning/nep.2.txt',
            './data/active_learning/nep.3.txt'
        ]
    ]
    system.potentials.append(neps[0])
    integrator = somd.core.integrators.nhc_integrator(0.001)
    writer = somd.apps.trajectories.H5WRITER(
        'traj0.h5', write_forces=True, use_double=True, potential_list=[0]
    )
    simulation = somd.apps.simulations.SIMULATION(
        system, integrator, trajectories=[writer]
    )
    simulation._initialize()
    simulation.system.forces[:] = 0
    simulation.run(1000)

    systems = [system.copy() for _ in range(3)]
    for i in range(3):
        systems[i].potentials.append(neps[i + 1])
        writer = somd.apps.trajectories.H5WRITER(
            'traj{:d}.h5'.format(i + 1), write_forces=True, use_double=True
        )
        task = somd.apps.evaluation.EVALUATION(
            'traj0.h5', systems[i], trajectories=[writer]
        )
        task.run()

    selector = somd.apps.select.STRUCTURESELECTOR([
        'traj{:d}.h5'.format(i) for i in range(4)
    ])
    indices = selector.select(100, 50, 250, False)
    assert all(indices == [344, 345, 853, 854, 919])
    msd = _np.loadtxt('data/active_learning/msd.dat')
    _nt.assert_array_almost_equal(
        [selector.get_force_msd(i) for i in indices], msd, DECIMAL_F
    )

    selector.write('traj4.h5', indices)
    reader_1 = somd.apps.trajectories.H5READER(
        'data/active_learning/traj.h5',
        read_velocities=False,
        read_forces=False,
        read_virial=False,
        read_nhc_data=False,
        read_rng_state=False,
    )
    reader_2 = somd.apps.trajectories.H5READER(
        'traj4.h5',
        read_velocities=False,
        read_forces=False,
        read_virial=False,
        read_nhc_data=False,
        read_rng_state=False,
    )
    for i in range(5):
        _nt.assert_array_equal(
            reader_1._read_snapshot(i).positions,
            reader_2._read_snapshot(i).positions,
        )

    for i in range(5):
        _os.remove('traj{:d}.h5'.format(i))
