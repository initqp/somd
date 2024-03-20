
import somd
import os as _os
import numpy as _np
import shutil as _sh
import warnings as _w
import numpy.testing as _nt

DECIMAL_F = 7
DECIMAL_D = 14

somd.utils.rng = somd.utils._rng.LEGACYRNG(1)


def test_nep():
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    system.positions[:] = [[0.1, 0.1, 0.1], [0.11, 0.11, 0.11]]
    potential = somd.potentials.NEP([0, 1], 'data/potentials/nep.txt',
                                    ["Si", "Si"], 0)
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_nep.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_F)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_F)

    potential = somd.potentials.NEP.generator([0, 1],
                                              'data/potentials/nep.txt',
                                              ["Si", "Si"], 0)()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_nep.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_F)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_F)


def test_nep_t():
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    system.positions[:] = [[0.1, 0.1, 0.1], [0.11, 0.11, 0.11]]
    potential = somd.potentials.NEP([0, 1], 'data/potentials/nep.txt',
                                    ["Si", "Si"], 1)
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_nep_t.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_F)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_F)

    potential = somd.potentials.NEP.generator([0, 1],
                                              'data/potentials/nep.txt',
                                              ["Si", "Si"], 1)()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_nep_t.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_F)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_F)


def test_dftd3():
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    potential = somd.potentials.DFTD3([0, 1], [13, 13], 'pbe')
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_dftd3.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)

    potential = somd.potentials.DFTD3.generator([0, 1], [13, 13], 'pbe')()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_dftd3.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)


def test_dftd4():
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    potential = somd.potentials.DFTD4([0, 1], [13, 13], 'pbe')
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_dftd4.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)

    potential = somd.potentials.DFTD4.generator([0, 1], [13, 13], 'pbe')()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_dftd4.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)


def test_plumed():
    try:
        import plumed
    except Exception as e:
        message = 'Can not import the PLUMED wrapper! PLUMED has NOT been ' + \
                  'tested! Reason: {} '.format(e)
        _w.warn(message)
        return
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    potential = somd.potentials.PLUMED([0, 1], 'data/potentials/plumed.inp',
                                       0.001, 1, cv_names=[{'d1': ''}])
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_plumed.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)
    result = _np.loadtxt('data/potentials/potential_plumed_cv.dat')
    _nt.assert_almost_equal(potential.cv_values[0], result, DECIMAL_D)
    _os.remove('plumed.inp.log')

    potential = somd.potentials.PLUMED.generator(
        [0, 1], 'data/potentials/plumed.inp', 0.001, 1,
        cv_names=[{'d1': ''}])()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_plumed.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], DECIMAL_D)
    _nt.assert_almost_equal(potential.virial, result[2:5], DECIMAL_D)
    result = _np.loadtxt('data/potentials/potential_plumed_cv.dat')
    _nt.assert_almost_equal(potential.cv_values[0], result, DECIMAL_D)
    _os.remove('plumed.inp.log')


def test_siesta():
    if (_os.environ.get('SIESTA_COMMAND') is None):
        message = 'Can not read the "SIESTA_COMMAND" environment ' + \
                  'variable! SIESTA has NOT been tested!'
        _w.warn(message)
        return
    else:
        command = _os.environ.get('SIESTA_COMMAND')
    system = somd.core.systems.create_system_from_poscar(
        'data/system/model.poscar')
    system.positions[:] = [[0.1, 0.1, 0.1], [0.11, 0.11, 0.11]]
    options = r"""
    xc.functional          GGA
    xc.authors             PBE
    PAO.BasisSize          DZP
    Mesh.Cutoff            300 Ry
    PAO.EnergyShift        10 meV
    PAO.SoftDefault        T
    DM.Tolerance           1.d-8
    SolutionMethod         diagon
    ElectronicTemperature  1 meV
    """
    potential = somd.potentials.SIESTA(
        [0, 1], system, options, command, 'data/potentials')
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_siesta.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], 5)
    _nt.assert_almost_equal(potential.virial, result[2:5], 5)
    potential.finalize()
    _sh.rmtree(potential.working_directory)

    potential = somd.potentials.SIESTA.generator(
        [0, 1], system, options, command, 'data/potentials')()
    potential.update(system)
    result = _np.loadtxt('data/potentials/potential_siesta.dat')
    _nt.assert_almost_equal(potential.forces, result[0:2], 5)
    _nt.assert_almost_equal(potential.virial, result[2:5], 5)
    potential.finalize()
    _sh.rmtree(potential.working_directory)
