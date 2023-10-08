import somd
import os as _os
import h5py as _h5
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_D = 14

_np.random.seed(1)


def test_trajectory_1():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    traj_writer = somd.apps.trajectories.H5WRITER(
        'traj.h5', write_velocities=True, use_double=True)
    simulation = somd.apps.simulations.SIMULATION(system, integrator,
                                                  trajectories=[traj_writer])
    simulation.run(1)
    file = _h5.File('traj.h5', 'r')
    _nt.assert_almost_equal(file['box'][0], system.box, DECIMAL_D)
    _nt.assert_almost_equal(file['coordinates'][0], system.positions,
                            DECIMAL_D)
    _nt.assert_almost_equal(file['velocities'][0], system.velocities,
                            DECIMAL_D)
    _nt.assert_almost_equal(file['cell_lengths'][0], system.lattice[:3],
                            DECIMAL_D)
    _nt.assert_almost_equal(file['cell_angles'][0], system.lattice[3:],
                            DECIMAL_D)
    _os.remove('traj.h5')


def test_trajectory_2():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    traj_writer = somd.apps.trajectories.H5WRITER('traj.h5', restart_file=True)
    simulation = somd.apps.simulations.SIMULATION(system, integrator,
                                                  trajectories=[traj_writer])
    simulation.run(1)
    file = _h5.File('traj.h5', 'r')
    _nt.assert_almost_equal(file['nhc_positions'][0, 0],
                            integrator._nhchains[0].positions, DECIMAL_D)
    _nt.assert_almost_equal(file['nhc_momentums'][0, 0],
                            integrator._nhchains[0].momentums, DECIMAL_D)
    _os.remove('traj.h5')


def test_trajectory_3():
    system = somd.core.systems.create_system_from_pdb(
        './data/active_learning/topo.pdb')
    system.groups.create_from_dict({'atom_list': range(0, 8)})
    potential_1 = somd.potentials.NEP(
        range(0, 8), './data/active_learning/nep.0.txt', system.atomic_symbols)
    potential_2 = somd.potentials.NEP(
        range(0, 8), './data/active_learning/nep.0.txt', system.atomic_symbols)
    system.potentials.append(potential_1)
    system.potentials.append(potential_2)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    traj_writer_1 = somd.apps.trajectories.H5WRITER(
        'traj_1.h5', write_forces=True, write_virial=True,
        write_velocities=True, use_double=True)
    traj_writer_2 = somd.apps.trajectories.H5WRITER(
        'traj_2.h5', write_forces=True, write_virial=True,
        write_velocities=True, potential_list=[1], use_double=True)
    simulation = somd.apps.simulations.SIMULATION(
        system, integrator, trajectories=[traj_writer_1, traj_writer_2])
    simulation.run(1)
    file_1 = _h5.File('traj_1.h5', 'r')
    file_2 = _h5.File('traj_2.h5', 'r')
    _nt.assert_almost_equal(file_1['forces'][0], system.forces, DECIMAL_D)
    _nt.assert_almost_equal(file_2['forces'][0], potential_2.forces, DECIMAL_D)
    _nt.assert_almost_equal(file_1['virial'][0], system.virial, DECIMAL_D)
    _nt.assert_almost_equal(file_2['virial'][0], potential_2.virial, DECIMAL_D)
    _nt.assert_almost_equal(file_1['potentialEnergy'][0],
                            system.energy_potential, DECIMAL_D)
    _nt.assert_almost_equal(file_2['potentialEnergy'][0],
                            potential_2.energy_potential, DECIMAL_D)
    _os.remove('traj_1.h5')
    _os.remove('traj_2.h5')


def test_trajectory_4():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    traj_writer = somd.apps.trajectories.H5WRITER(
        'traj.h5', write_velocities=True, use_double=True)
    simulation = somd.apps.simulations.SIMULATION(system, integrator,
                                                  trajectories=[traj_writer])
    simulation.run(1)
    reader = somd.apps.trajectories.H5READER('traj.h5', 'r')
    snapshot = reader._read_snapshot()
    _nt.assert_almost_equal(snapshot.box, system.box, DECIMAL_D)
    _nt.assert_almost_equal(snapshot.positions, system.positions, DECIMAL_D)
    _nt.assert_almost_equal(snapshot.velocities, system.velocities, DECIMAL_D)
    _nt.assert_almost_equal(snapshot.lattice, system.lattice, DECIMAL_D)
    _os.remove('traj.h5')


def test_trajectory_5():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator_1 = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    traj_writer = somd.apps.trajectories.H5WRITER('traj.h5', restart_file=True)
    simulation = somd.apps.simulations.SIMULATION(system, integrator_1,
                                                  trajectories=[traj_writer])
    simulation.run(1)
    reader = somd.apps.trajectories.H5READER('traj.h5', 'r')
    integrator_2 = integrator_1.copy()
    reader.bind_integrator(integrator_1)
    reader._read_nhc_data()
    _nt.assert_almost_equal(integrator_1._nhchains[0].positions,
                            integrator_2._nhchains[0].positions, DECIMAL_D)
    _nt.assert_almost_equal(integrator_1._nhchains[0].momentums,
                            integrator_2._nhchains[0].momentums, DECIMAL_D)
    _os.remove('traj.h5')
