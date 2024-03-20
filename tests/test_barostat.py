import somd
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

somd.utils.rng = somd.utils._rng.LEGACYRNG(1)


def test_barostat_1():
    system = _h.get_harmonic_system()
    system.box[:] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    barostat = somd.apps.barostats.BAROSTAT([1E-5], [1E-5], 0.1)
    integrator.bind_system(system)
    barostat.bind_integrator(integrator)
    integrator.propagate()
    barostat.update()
    result = _np.loadtxt('data/barostat/barostat_1.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.box, result[4:7], DECIMAL_D)


def test_barostat_2():
    system = _h.get_harmonic_system()
    system.box[:] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    somd.utils.defaults.NHCLENGTH = 6
    somd.utils.defaults.NHCNRESPA = 4
    integrator = somd.core.integrators.nhc_integrator(
        0.001, relaxation_times=[0.01])
    barostat = somd.apps.barostats.BAROSTAT(
        [1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5],
        [1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5], 0.1)
    integrator.bind_system(system)
    barostat.bind_integrator(integrator)
    barostat.initialize()
    integrator.propagate()
    barostat.update()
    result = _np.loadtxt('data/barostat/barostat_2.dat')
    _nt.assert_almost_equal(system.positions, result[0:4], DECIMAL_D)
    _nt.assert_almost_equal(system.box, result[4:7], DECIMAL_D)
