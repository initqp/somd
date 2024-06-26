import somd
import os as _os
import numpy as _np
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
        0.001, relaxation_times=[0.01])
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.run(1)
    result = _np.loadtxt('data/integrators/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_run_2():
    somd.utils.rng.seed(1)
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.gobabo_integrator(
        0.0005, relaxation_times=[0.01])
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.run(1)
    result = _np.loadtxt('data/integrators/integrators_gobabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_restart_1():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    system.positions[:] = 0.0
    system.velocities[:] = 0.0
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.restart_from('data/simulation/restart_1.h5')
    simulation.run(1)
    result = _np.loadtxt('data/integrators/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_restart_2():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.gobabo_integrator(
        0.0005, relaxation_times=[0.01])
    system.positions[:] = 0.0
    system.velocities[:] = 0.0
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.restart_from('data/simulation/restart_2.h5', read_nhc_data=False)
    simulation.run(1)
    result = _np.loadtxt('data/integrators/integrators_gobabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_restart_3():
    system = _h.get_harmonic_system()
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        -0.001, relaxation_times=[0.01])
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.restart_from('data/simulation/restart_3.h5')
    simulation.run(1)
    result = _h.get_harmonic_system()
    _nt.assert_almost_equal(system.positions, result.positions, DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result.velocities, DECIMAL_D)


def test_restart_4():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    system.positions[:] = 0.0
    system.velocities[:] = 0.0
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.restart_from('data/simulation/restart_1.h5', read_cell=False)
    simulation.dump_restart('restart.h5')
    simulation.restart_from('restart.h5')
    simulation.run(1)
    result = _np.loadtxt('data/integrators/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)
    _os.remove('restart.h5')
