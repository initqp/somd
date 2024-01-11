import somd
import numpy as _np
import harmonic as _h
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14


def test_operators():
    system = _h.get_harmonic_system()
    g1 = somd.core.groups.ATOMGROUP(system, [0, 1, 2])
    g2 = somd.core.groups.ATOMGROUP(system, [1, 2, 3])
    g3 = somd.core.groups.ATOMGROUP(system, [0, 1, 2, 3])
    assert (len(g1 & g2) != 0)
    assert (g1 in system.groups[0])
    assert (system.groups[0] == g3)
    system.groups.create_from_dict({'atom_list': [0, 1, 2, 3]})
    assert (len(system.groups) == 1)


def test_n_constraints():
    system = _h.get_harmonic_system()
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    assert (system.groups[0].n_constraints == 3)
    system.constraints.pop(0)
    assert (system.groups[0].n_constraints == 2)
    system.constraints.pop(0)
    assert (system.groups[0].n_constraints == 1)
    c = {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}
    assert (system.constraints[0] == c)
    try:
        system.groups.create_from_dict({"atom_list": [0, 1, 2]})
    except:
        pass
    else:
        raise AssertionError


def test_dof_1():
    system = _h.get_harmonic_system()
    assert (system.groups[0].n_dof == 12)
    system.groups[0].has_translations = False
    assert (system.groups[0].n_dof == 9)
    c = [{'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14},
         {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14},
         {'type': 2, 'indices': [0, 1, 2, 3], 'target': 0, 'tolerance': 1E-14}]
    system.constraints.appends(c)
    assert (system.groups[0].n_dof == 6)


def test_dof_2():
    system = _h.get_harmonic_system()
    system.groups.create_from_dict({"atom_list": [0, 1, 2]})
    c = {'type': 1, 'indices': [0, 1, 2], 'target': 1.57, 'tolerance': 1E-14}
    system.constraints.append(c)
    assert (system.groups[0].n_dof == 11)
    assert (system.groups[1].n_dof == 8)
    system.groups[0].has_translations = False
    assert (system.groups[0].n_dof == 8)
    assert (system.groups[1].n_dof == 8)


def test_dof_3():
    system = _h.get_harmonic_system()
    system.groups.create_from_dict({"atom_list": [0, 1]})
    system.groups.create_from_dict({"atom_list": [2, 3]})
    c = {'type': 0, 'indices': [0, 1], 'target': 1.0, 'tolerance': 1E-14}
    system.constraints.append(c)
    assert (system.groups[0].n_dof == 11)
    assert (system.groups[1].n_dof == 5)
    assert (system.groups[2].n_dof == 6)
    system.groups[1].has_translations = False
    assert (system.groups[0].n_dof == 8)
    assert (system.groups[1].n_dof == 2)
    assert (system.groups[2].n_dof == 6)
    system.groups[2].has_translations = False
    assert (system.groups[0].n_dof == 5)
    assert (system.groups[1].n_dof == 2)
    assert (system.groups[2].n_dof == 3)
    system.groups.pop(2)
    assert (system.groups[0].n_dof == 8)
    assert (system.groups[1].n_dof == 2)
    system.groups.pop(1)
    assert (system.groups[0].n_dof == 11)


def test_velocities():
    system = _h.get_harmonic_system()
    system.groups[0].velocities = 0.0
    result = _np.zeros((4, 3), dtype=_np.double)
    _nt.assert_almost_equal(system.velocities, result, DECIMAL_D)


def test_positions():
    system = _h.get_harmonic_system()
    system.groups[0].positions = 0.0
    result = _np.zeros((4, 3), dtype=_np.double)
    _nt.assert_almost_equal(system.positions, result, DECIMAL_D)


def test_temperature():
    system = _h.get_harmonic_system()
    result = _np.loadtxt('data/groups/group_temperature.dat')
    _nt.assert_almost_equal(system.groups[0].temperature, result, DECIMAL_D)


def test_initlization():
    system = _h.get_harmonic_system()
    system.groups[0].velocities = 0.0
    system.groups[0].add_velocities_from_temperature(10)
    result = _np.double(10.0)
    _nt.assert_almost_equal(system.groups[0].temperature, result, DECIMAL_D)
    result = _np.zeros(3, dtype=_np.double)
    _nt.assert_almost_equal(system.groups[0].com_velocities, result, DECIMAL_D)
    _np.random.seed(1)


def test_com():
    system = _h.get_harmonic_system()
    system.groups[0].remove_com_motion()
    result = _np.zeros(3, dtype=_np.double)
    _nt.assert_almost_equal(system.groups[0].com_velocities, result, DECIMAL_D)
    system.groups[0].com_positions += [1.0, 1.0, 1.0]
    result = _np.array([1.5, 1.5, 1.0], dtype=_np.double)
    _nt.assert_almost_equal(system.groups[0].com_positions, result, DECIMAL_D)
