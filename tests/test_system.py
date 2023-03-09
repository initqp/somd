import somd
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

_np.random.seed(1)


def test_pdb():
    system = somd.core.systems.create_system_from_pdb('data/model.pdb')
    result = _np.array([[1.0, 1.0, 1.0], [2.1, 2.2, 2.3]], dtype=_np.double)
    _nt.assert_array_almost_equal(system.positions, result, DECIMAL_F)


def test_poscar():
    system = somd.core.systems.create_system_from_poscar('data/model.poscar')
    result = _np.array([[1.0, 1.0, 1.0], [2.1, 2.2, 5.3]], dtype=_np.double)
    _nt.assert_array_almost_equal(system.positions, result, DECIMAL_D)


def test_lattice():
    system = somd.core.systems.create_system_from_poscar('data/model.poscar')
    length = _np.sqrt(2) * 2
    result = _np.array([length, length, length, 60, 60, 60])
    _nt.assert_array_almost_equal(system.lattice, result, DECIMAL_D)


def test_wrap():
    system = somd.core.systems.create_system_from_poscar('data/model.poscar')
    _nt.assert_array_almost_equal(
        system.positions_wrapped,
        _np.array([[1.0, 1.0, 1.0], [2.1, 2.2, 1.3]]), DECIMAL_D)


def test_volume():
    system = somd.core.systems.create_system_from_poscar('data/model.poscar')
    _nt.assert_almost_equal(system.volume, _np.double(16.0), DECIMAL_D)


def test_pressures():
    system = somd.core.systems.create_system_from_poscar('data/model.poscar')
    system.velocities[:] = 1.0
    result = _np.loadtxt('data/system_pressures.dat')
    _nt.assert_array_almost_equal(system.pressures, result, DECIMAL_D)


def test_segments():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 0, 'indices': [2, 3], 'target': 1.57, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    assert (len(system.segments) == 2)
    c = [{'type': 0, 'indices': [1, 2], 'target': 1.0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    assert (len(system.segments) == 1)
