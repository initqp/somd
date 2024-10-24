import somd
import os as _os
import harmonic as _h
import numpy.testing as _nt

DECIMAL_D = 14

somd.utils.rng = somd.utils._rng.LEGACYRNG(1)


def test_run_1():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01]
    )
    writer = somd.apps.trajectories.H5WRITER(
        'traj.h5',
        write_velocities=True,
        write_forces=True,
        write_virial=True,
        use_double=True,
    )
    simulation = somd.apps.simulations.SIMULATION(
        system, integrator, trajectories=[writer]
    )
    energies = []
    systems_reference = []
    for i in range(10):
        simulation.run(1)
        energies.append(simulation.system.energy_potential)
        systems_reference.append(simulation.system.copy())
    del writer

    writer = somd.apps.trajectories.H5WRITER(
        'eval.h5',
        write_velocities=True,
        write_forces=True,
        write_virial=True,
        use_double=True,
    )
    evaluation = somd.apps.evaluation.EVALUATION(
        'traj.h5', system, trajectories=[writer], interval=3
    )
    evaluation.run()

    reader = somd.apps.trajectories.H5READER('eval.h5')
    for i, j in enumerate([0, 3, 6, 9]):
        snapshot = reader._read_snapshot(i)
        _nt.assert_almost_equal(
            snapshot.positions, systems_reference[j].positions, DECIMAL_D
        )
        _nt.assert_almost_equal(
            snapshot.velocities, systems_reference[j].velocities, DECIMAL_D
        )
        _nt.assert_almost_equal(
            snapshot.forces, systems_reference[j].forces, DECIMAL_D
        )
        _nt.assert_almost_equal(
            snapshot.virial, systems_reference[j].virial, DECIMAL_D
        )
        assert energies[j] == reader.root['potentialEnergy'][i]
    _os.remove('traj.h5')
    _os.remove('eval.h5')
