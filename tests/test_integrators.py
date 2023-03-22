import somd
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

_np.random.seed(1)


def test_vv():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.vv_integrator(0.001)
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_vv.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_cs4():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.cs4_integrator(0.001)
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_cs4.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_nhc():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_nhc_copy_system():
    s1 = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    s1.constraints.appends(c)
    system = s1.copy()
    system.potentials.append(s1.potentials[0])
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_nhc_copy_integrator():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    i1 = somd.core.integrators.nhc_integrator(0.001, relaxation_times=[0.01])
    integrator = i1.copy()
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_nhc_copy_integrator_rev():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.constants.SOMDDEFAULTS.NHCLENGTH = 6
    somd.constants.SOMDDEFAULTS.NHCNRESPA = 4
    i1 = somd.core.integrators.nhc_integrator(0.001, relaxation_times=[0.01])
    i1.bind_system(system)
    i1.propagate()
    integrator = i1.copy()
    integrator.bind_system(system)
    integrator.timestep *= -1.0
    integrator.propagate()
    integrator.timestep *= -1.0
    integrator.propagate()
    result = _np.loadtxt('data/integrators_nhc.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_F)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_F)


def test_obabo():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    integrator = somd.core.integrators.obabo_integrator(
        0.001, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_obabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_baoab():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    integrator = somd.core.integrators.obabo_integrator(
        0.001, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_baoab.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_gobabo():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.gobabo_integrator(
        0.0005, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_gobabo.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_gbaoab():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    integrator = somd.core.integrators.gbaoab_integrator(
        0.0005, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_gbaoab.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_gbaoab_copy_integrator():
    _np.random.seed(1)
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    i1 = somd.core.integrators.gbaoab_integrator(
        0.0005, relaxation_times=[0.01])
    integrator = i1.copy()
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_gbaoab.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)


def test_gbaoab_copy_system():
    _np.random.seed(1)
    s1 = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    s1.constraints.appends(c)
    system = s1.copy()
    system.potentials.append(s1.potentials[0])
    integrator = somd.core.integrators.gbaoab_integrator(
        0.0005, relaxation_times=[0.01])
    integrator.bind_system(system)
    integrator.propagate()
    result = _np.loadtxt('data/integrators_gbaoab.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.velocities, result[4:8], DECIMAL_D)
