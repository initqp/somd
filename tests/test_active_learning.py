import somd
import os as _os
import h5py as _h5
import numpy as _np
import shutil as _sh
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

_np.random.seed(1)


def test_active_learning():
    indices = [344, 345, 853, 854, 919]
    system = somd.core.systems.create_system_from_pdb(
        './data/active_learning/topo.pdb')
    system.groups.create_from_dict({'atom_list': range(0, 8)})
    integrator = somd.core.integrators.nhc_integrator(0.001)

    def get_potential():
        return somd.potentials.DFTD3(
            list(range(0, 8)), system.atomic_types, 'pbe')

    if (_os.path.exists('active_learning.h5')):
        _os.remove('active_learning.h5')
    if (_os.path.exists('active_learning.dir')):
        _sh.rmtree('active_learning.dir')
    trainer = somd.apps.active_learning.ACTIVELEARNING(
        system, integrator, [get_potential], [0],
        {'max_md_runs_per_iter': 1,
         'max_md_steps_per_iter': 1000,
         'msd_lower_limit': 50,
         'msd_upper_limit': 250,
         'initial_training_set': './data/active_learning/topo.pdb',
         'initial_potential_files': ['./data/active_learning/nep.0.txt',
                                     './data/active_learning/nep.1.txt',
                                     './data/active_learning/nep.2.txt',
                                     './data/active_learning/nep.3.txt']},
        'true', 'null')
    trainer.run(1)
    h5_root = _h5.File('active_learning.h5')['iteration_data/1']
    _nt.assert_array_equal(h5_root['accepted_structure_indices'], indices)
    energies = []
    potential = get_potential()
    integrator.bind_system(system)
    reader = somd.apps.trajectories.H5READER(
        'active_learning.dir/iteration_1/visited_structures.h5',
        read_forces=False, read_nhc_data=False, read_rng_state=False,
        read_velocities=False)
    reader.bind_integrator(integrator)
    for i in indices:
        reader.read(i)
        potential.update(system)
        energies.append(potential.energy_potential[0].copy())
    _nt.assert_array_almost_equal(
        _np.array(h5_root['accepted_structure_energies']),
        _np.array(energies), DECIMAL_D)
    energies = []
    reader = somd.apps.trajectories.H5READER(
        'data/active_learning/traj.h5',
        read_forces=False, read_nhc_data=False, read_rng_state=False,
        read_velocities=False)
    reader.bind_integrator(integrator)
    for i in range(0, 5):
        reader.read(i)
        potential.update(system)
        energies.append(potential.energy_potential[0].copy())
    _nt.assert_array_almost_equal(
        _np.array(h5_root['accepted_structure_energies']),
        _np.array(energies), DECIMAL_D)
    msd = _np.loadtxt('data/active_learning/msd.dat')
    _nt.assert_array_almost_equal(_np.array(h5_root['force_msd'])[indices],
                                  _np.array(msd), DECIMAL_F)
    _os.remove('active_learning.h5')
    _sh.rmtree('active_learning.dir')
