import somd
import os as _os
import h5py as _h5
import numpy as _np
import shutil as _sh
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

somd.utils.rng = somd.utils._rng.LEGACYRNG(1)


def test_parser_1():
    parser = somd.apps.parser.TOMLPARSER('./data/parser/1.toml')
    traj = somd.apps.trajectories.H5READER('./data/active_learning/traj.h5')
    integrator = somd.core.integrators.gbaoab_integrator(0.001)
    assert parser.simulation.system.n_atoms == 8
    assert parser.simulation.system.constraints.__len__() == 3
    assert parser.simulation.system.constraints[0] == \
        {'indices': [4, 7], 'target': 0.108, 'tolerance': 1e-14,
         'type': 0}
    assert parser.simulation.system.constraints[1] == \
        {'indices': [0, 4, 7], 'target': 1.92527, 'tolerance': 1e-14,
         'type': 1}
    assert parser.simulation.system.constraints[2] == \
        {'indices': [2, 0, 4, 7], 'target': -0.9995, 'tolerance': 1e-14,
         'type': 2}
    assert parser.simulation.system.potentials[0].__class__.__name__ == 'NEP'
    assert parser.simulation.integrator.timestep == 0.001
    assert parser.simulation.integrator.splitting == integrator.splitting
    assert parser.simulation.integrator.thermo_groups == [2]
    _nt.assert_array_equal(parser.simulation.system.positions,
                           traj._read_snapshot(-1).positions)
    _nt.assert_array_equal(parser.simulation.system.groups[0].atom_list,
                           [1, 3])
    _nt.assert_array_equal(parser.simulation.system.groups[1].atom_list,
                           [5, 6])
    _nt.assert_array_equal(parser.simulation.system.groups[2].atom_list,
                           range(0, 8))
    _nt.assert_array_almost_equal(
        [parser.simulation.system.groups[2].temperature], [300], DECIMAL_F)


def test_parser_2():
    somd.utils.rng.seed(1)
    indices = [68, 170, 183, 268, 370, 383]
    restart_function = somd.apps.simulations.SIMULATION.restart_from
    somd.apps.simulations.SIMULATION.restart_from = \
        lambda *args, read_nhc_data = True, **kwargs: restart_function(
            *args, read_nhc_data, **kwargs)
    parser = somd.apps.parser.TOMLPARSER('./data/parser/2.toml')
    parser.simulation.system.groups[0].has_translations = True
    parser.simulation.integrator._nhchains[0].n_dof = 24
    parser.run()
    h5_root = _h5.File('2.active_learning.h5')['iteration_data/1']
    _nt.assert_array_equal(h5_root['accepted_structure_indices'], indices)
    msd = _np.loadtxt('data/active_learning/msd.dat')
    msd = [msd[0], msd[3], msd[4]] * 2
    _nt.assert_array_almost_equal(_np.array(h5_root['force_msd'])[indices],
                                  _np.array(msd), DECIMAL_F)
    somd.apps.simulations.SIMULATION.restart_from = restart_function
    _os.remove('2.active_learning.h5')
    _sh.rmtree('2.active_learning.dir')
