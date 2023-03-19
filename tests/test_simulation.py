import somd
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_D = 14

_np.random.seed(1)


def test_run_1():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.run(1)
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_run_2():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.gobabo_integrator(
        0.0005, relaxation_times=[0.01])
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.run(1)
    result = _np.loadtxt('data/integrators_gobabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_restart_1():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    system.positions[:] = 0.0
    system.velocities[:] = 0.0
    simulation = somd.apps.simulations.SIMULATION(system, integrator)
    simulation.restart_from('data/restart_1.h5')
    simulation.run(1)
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_restart_2():
    _np.random.seed(1)
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
    simulation.restart_from('data/restart_2.h5', read_nhc_data=False)
    simulation.run(1)
    result = _np.loadtxt('data/integrators_gobabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)
