import somd
import os as _os
import h5py as _h5
import numpy as _np
import shutil as _sh
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

_np.random.seed(1)


def test_states_1():

    def i1(x):
        import numpy as np
        d = np.sqrt((x[0] + 1) ** 2 + (x[1] + 1) ** 2)
        return (d < 0.1)

    def i2(x):
        import numpy as np
        d = np.sqrt((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
        return (d < 0.1)

    STATE = somd.apps.path_sampling.utils.state.STATE
    s1 = STATE(i1, 's1')
    s2 = STATE(i2, 's2')
    assert (s1 | s2)([1, 1])
    assert (s1 | s2)([-1, -1])
    assert not (s1 | s2)([0, 0])
    assert (s1 ^ s2)([1, 1])
    assert (s1 ^ s2)([-1, -1])
    assert not (s1 ^ s2)([0, 0])
    assert not (s1 & s2)([1, 1])
    assert not (s1 & s2)([0, 0])
    assert not (s1 & s2)([-1, -1])
    assert (s1 - s2)([-1, -1])
    assert not (s1 - s2)([1, 1])
    assert not (s1 - s2)([0, 0])
    assert (-(s1 ^ s2))([0, 0])
    assert not (-(s1 ^ s2))([1, 1])
    assert not (-(s1 ^ s2))([-1, -1])
    assert (-(s1 | s2))([0, 0])
    assert not (-(s1 | s2))([1, 1])
    assert not (-(s1 | s2))([-1, -1])
    assert ((-(s1 & s2)) | (s1 | s2))([1, 1])
    assert ((-(s1 & s2)) | (s1 | s2))([0, 0])
    assert ((-(s1 & s2)) | (s1 | s2))([-1, -1])
    assert ((-(s1 & s2)) | (s1 ^ s2))([1, 1])
    assert ((-(s1 & s2)) | (s1 ^ s2))([0, 0])
    assert ((-(s1 & s2)) | (s1 ^ s2))([-1, -1])
    assert ((-(s1 & s2)) | (s1 ^ s2)).label == '((!(s1&s2))|(s1^s2))'
    assert len(((-(s1 & s2)) | (s1 | s2)).to_string().split('|')) == 6
    assert ((-(s1 & s2)) | (s1 | s2)).to_string().split('|')[2] == \
           s1.to_string().split('|')[2]
    assert ((-(s1 & s2)) | (s1 | s2)).to_string().split('|')[3] == \
           s2.to_string().split('|')[2]
    assert ((-(s1 & s2)) | (s1 | s2)).to_string().split('|')[4] == \
           s1.to_string().split('|')[2]
    assert ((-(s1 & s2)) | (s1 | s2)).to_string().split('|')[5] == \
           s2.to_string().split('|')[2]


def test_states_2():

    def i1(x):
        import numpy as np
        d = np.sqrt((x[0] + 1) ** 2 + (x[1] + 1) ** 2)
        return (d < 0.1)

    def i2(x):
        import numpy as np
        d = np.sqrt((x[0] - 1) ** 2 + (x[1] - 1) ** 2)
        return (d < 0.1)

    def i3(x):
        import numpy as np
        d = np.sqrt((x[0] - 2) ** 2 + (x[1] - 2) ** 2)
        return (d < 0.1)

    _s = somd.apps.path_sampling.utils.state
    s1 = _s.STATE(i1, 's1')
    s2 = _s.STATE(i2, 's2')
    s3 = _s.STATE(i3, 's3')
    assert _s.in_sequence([s1, s2], [[-1, -1], [1, 1]])
    assert _s.in_sequence([s1, s2, s3], [[-1, -1], [1, 1], [2, 2]])
    assert not _s.in_sequence([s1, s2, s3], [[-1, -1], [1, 1]])
    assert not _s.in_sequence([s1, s2, s3], [[-1, -1], [2, 2]])
    assert not _s.in_sequence([s1, s2, s3],
                              [[-1, -1], [1, 1], [-1, -1], [2, 2]])
    assert not _s.in_sequence([s1, s2, s3],
                              [[-1, -1], [1, 1], [-1, -1], [2, 2]], True)

    cv_values = _np.array([[-1, -1], [-0.5, -0.5], [0, 0], [0.5, 0.5], [1, 1]])
    _nt.assert_array_equal(_s.get_region_indices(s1, cv_values), [0])
    _nt.assert_array_equal(_s.get_region_indices(-(s1), cv_values),
                           [1, 2, 3, 4])
    assert _s.get_reactive_regions([s1, s2], cv_values)[0]['range'] == (0, 5)
    assert _s.get_reactive_regions([s2, s1], cv_values)[0]['range'] == (0, 5)
    assert _s.get_reactive_regions([s2, s1], cv_values)[0]['direction'] == -1


def test_path():
    s = somd.core.systems.create_system_from_pdb(
        './data/path_sampling/topo.pdb')
    p1 = somd.apps.path_sampling.utils.path.PATH(
        s, ['./data/path_sampling/path1.h5'], 0.001)
    p2 = somd.apps.path_sampling.utils.path.PATH(
        s, ['./data/path_sampling/path1.h5'], 0.001,
        reverse_velocities=[True])
    p3 = p1 + reversed(p2)
    _nt.assert_array_equal(p3.frame(0).positions,
                           p3.frame(-1).positions)
    _nt.assert_array_equal(p3.frame(0).velocities,
                           p3.frame(-1).velocities * -1)
    p4 = reversed(p3[(len(p1) - 1):(len(p1) + 1)])
    assert len(p4) == 2
    assert p4.reverse_velocities == [True, False]
    _nt.assert_array_equal(p4.frame(0).positions,
                           p2.frame(-1).positions)
    _nt.assert_array_equal(p4.frame(0).velocities,
                           p1.frame(-1).velocities * -1)
    p4.reverse_velocities[0] = False
    _nt.assert_array_equal(p4.frame(0).velocities,
                           p2.frame(-1).velocities * -1)
    n = len(p1)
    p5 = p1 + reversed(p4[0]) + reversed(p2)
    _nt.assert_array_equal(p5.frame(n).positions,
                           p5.frame(n + 1).positions)
    _nt.assert_array_equal(p5.frame(n).velocities,
                           p5.frame(n + 1).velocities * -1)


def test_selections():
    selection = somd.apps.path_sampling.utils.selection
    cv_values = _np.array([[-1, -1], [0, 0], [1, 1]])
    p = selection.calculate_frame_probabilities(cv_values)
    _nt.assert_array_almost_equal(p, _np.array([1 / 3, 1 / 3, 1 / 3]),
                                  DECIMAL_D)
    p = selection.calculate_frame_probabilities(cv_values,
                                                lambda x: x[0] + x[1] == 0)
    _nt.assert_array_almost_equal(p, _np.array([0.0, 1.0, 0.0]), DECIMAL_D)
    p = selection.calculate_frame_probabilities(cv_values,
                                                lambda x: x[0] + 1.5)
    _nt.assert_array_almost_equal(p, _np.array([1 / 9, 1 / 3, 5 / 9]),
                                  DECIMAL_D)
    p = selection.calculate_frame_probabilities(cv_values, lambda x: True, 1)
    _nt.assert_array_almost_equal(p, _np.array([0.0, 1.0, 0.0]), DECIMAL_D)
    try:
        p = selection.calculate_frame_probabilities(cv_values, lambda x: x[0])
    except:
        pass
    else:
        raise RuntimeError
    s = selection.select(cv_values, lambda x: True, 1)
    assert s['index'] == 1
    _nt.assert_array_almost_equal(s['probability'], _np.ones(1), DECIMAL_D)


def test_shooting():
    s = somd.core.systems.create_system_from_pdb(
        './data/path_sampling/topo.pdb')
    s.groups.create_from_dict({'atom_list': range(0, s.n_atoms)})
    p = somd.apps.path_sampling.utils.path.PATH(
        s, ['./data/path_sampling/path1.h5'], 0.001)
    s.snapshot = p.frame(0)
    energy = s.groups[0].energy_kinetic
    velocities = s.velocities.copy()
    somd.apps.path_sampling.utils.shooting.shoot(s)
    _nt.assert_almost_equal(energy, s.groups[0].energy_kinetic)
    try:
        _nt.assert_array_almost_equal(velocities, s.velocities, DECIMAL_D)
    except:
        pass
    else:
        raise RuntimeError


def test_tps_run_1():
    _np.random.seed(1)
    f1 = somd.core.groups.ATOMGROUP.add_velocities_from_temperature
    f2 = somd.core.groups.ATOMGROUP.n_dof
    import data.path_sampling.model as _model
    if (_os.path.exists('path_sampling.h5')):
        _os.remove('path_sampling.h5')
    if (_os.path.exists('path_sampling.dir')):
        _sh.rmtree('path_sampling.dir')
    STATE = somd.apps.path_sampling.utils.state.STATE
    sampler = somd.apps.path_sampling.tps.FELTWTPS(
        _model.system, _model.integrator,
        [lambda: _model.PES([0], 3 * _model.kb * 300)],
        {'plumed_file': _os.path.abspath('./data/path_sampling/plumed.inp'),
         'cv_names': [{'p': 'x'}, {'p': 'y'}],
         'states': [STATE(_model.i1, 'S1'), STATE(_model.i2, 'S2')],
         'sampled_paths': [[0, 1], [1, 0]],
         'initial_trajectory':
            _os.path.abspath('./data/path_sampling/initial.h5'),
         'bias_function': _model.i3,
         'randomize_velocities': True,
         'trajectory_interval': 5,
         'trajectory_use_double': True,
         'remove_dead_iterations': True},
        post_step_objects=[_model.helper])
    sampler.run(2)
    somd.core.groups.ATOMGROUP.add_velocities_from_temperature = f1
    somd.core.groups.ATOMGROUP.n_dof = f2
    f1 = _h5.File('./path_sampling.h5', 'r')
    f2 = _h5.File('./data/path_sampling/1.h5', 'r')
    assert f1['/iteration_data/2'].attrs['is_reactive']
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/segments/0/cv_values'],
        f2['/2/segments/0/cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/segments/1/cv_values'],
        f2['/2/segments/1/cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/shooting_point_cv_values'],
        f2['/2/shooting_point_cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/shooting_point_index'],
        f2['/2/shooting_point_index'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/selection_probability'],
        f2['/2/selection_probability'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/selection_probability_reversed'],
        f2['/2/selection_probability_reversed'], DECIMAL_D)
    _os.remove('path_sampling.h5')
    _sh.rmtree('path_sampling.dir')


def test_tps_run_2():
    _np.random.seed(1)
    f1 = somd.core.groups.ATOMGROUP.add_velocities_from_temperature
    f2 = somd.core.groups.ATOMGROUP.n_dof
    import data.path_sampling.model as _model
    if (_os.path.exists('path_sampling.h5')):
        _os.remove('path_sampling.h5')
    if (_os.path.exists('path_sampling.dir')):
        _sh.rmtree('path_sampling.dir')
    STATE = somd.apps.path_sampling.utils.state.STATE
    sampler = somd.apps.path_sampling.tps.FELOWTPS(
        _model.system, _model.integrator,
        [lambda: _model.PES([0], 3 * _model.kb * 300)],
        {'plumed_file': _os.path.abspath('./data/path_sampling/plumed.inp'),
         'cv_names': [{'p': 'x'}, {'p': 'y'}],
         'states': [STATE(_model.i1, 'S1'), STATE(_model.i2, 'S2')],
         'sampled_paths': [[0, 1], [1, 0]],
         'initial_trajectory':
            _os.path.abspath('./data/path_sampling/initial.h5'),
         'bias_function': _model.i3,
         'trajectory_interval': 5,
         'trajectory_use_double': True,
         'remove_dead_iterations': True},
        post_step_objects=[_model.helper])
    sampler.run(2)
    somd.core.groups.ATOMGROUP.add_velocities_from_temperature = f1
    somd.core.groups.ATOMGROUP.n_dof = f2
    f1 = _h5.File('./path_sampling.h5', 'r')
    f2 = _h5.File('./data/path_sampling/2.h5', 'r')
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/segments/0/cv_values'],
        f2['/2/segments/0/cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/segments/1/cv_values'],
        f2['/2/segments/1/cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/shooting_point_cv_values'],
        f2['/2/shooting_point_cv_values'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/shooting_point_index'],
        f2['/2/shooting_point_index'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/selection_probability'],
        f2['/2/selection_probability'], DECIMAL_D)
    _nt.assert_array_almost_equal(
        f1['/iteration_data/2/selection_probability_reversed'],
        f2['/2/selection_probability_reversed'], DECIMAL_D)
    _os.remove('path_sampling.h5')
    _sh.rmtree('path_sampling.dir')
