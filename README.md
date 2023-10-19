# SOMD (A ~~SIESTA Oriented~~ Shitty Opinionated Molecular Dynamics Package)

![logo](doc/logo.png "logo")

SOMD is an ab-initio molecular dynamics (AIMD) package designed for the
[SIESTA](https://departments.icmab.es/leem/siesta/) DFT code. The SOMD code
provides some common functionalities to perform standard Born-Oppenheimer
molecular dynamics (BOMD) simulations, and contains a simple wrapper to the
[Neuroevolution Potential (NEP)](https://github.com/brucefan1983/NEP_CPU)
package. The SOMD code may be used to automatically build NEPs by the mean
of the active-learning methodology.

## NOTE!
The SOMD code is designed to be maintained by one person, thus many important
functionalities may be absent. Besides, the code should be considered
**EXPERIMENTAL** since it has not been extensively tested. So if you
would like to perform production runs with SOMD, please take your own risk.

## INSTALLATION
SOMD only runs on GNU/Linux distros. The installation requires a working `g++`
compiler (with C++11 supports), a Python3 interpreter and four additional
Python3 libraries (`cython`, `h5py`, `mdtraj` and `toml`). You could install
SOMD by the following steps.

1. Install the required dependencies with:
    ```bash
    conda config --add channels conda-forge
    conda install cython h5py mdtraj toml -c conda-forge
    ```
    or
    ```bash
    pip install cython h5py mdtraj toml
    ```
2. Clone this repo:
    ```bash
    git clone https://www.github.com/initqp/somd
    cd somd
    git submodule update --init
    ```
    Note: if you would like to proceed your installation with the tarball
    downloaded from GitHub, you should manually download the NEP_CPU package
    and put it in the `somd/somd/potentials/src` directory. Besides, the
    version number of the installed package may be wrong.
3. Install SOMD:
    ```bash
    python setup.py install
    ```
    or
    ```bash
    pip install .
    ```
4. Start a `python` REPL and enter the following lines:
    ```python
    >>> import somd
    >>> print(somd.__version__)
    ```
    If the installation is successful, a version string should be printed.
    Likewise, you could enter the following command under your shell:
    ```bash
    somd -v
    ```
    If the installation is successful, a version string should be printed as
    well.
5.  Compile the [SIESTA](https://departments.icmab.es/leem/siesta/) code. SOMD
    could work with the 4.1.5 or the
    [git master](https://gitlab.com/siesta-project/siesta) version of SIESTA.
    When compiling, you are suggested to link your binary against the ELPA
    library (and using ELPA as the diagonalization algorithm). This is because
    of one of the memory leakage bugs in SIESTA (read
    [this page](https://gitlab.com/siesta-project/siesta/-/issues/29) for
    details). The usage of the ELPA library could be found in the SIESTA
    documentation.
6.  If you would like to use DFTD3, DFTD4 and PLUMED with SOMD, you should also
    install the corresponding packages:
    ```bash
    conda install dftd3-python dftd4-python py-plumed -c conda-forge
    ```
    or
    ```bash
    pip install dftd3 dftd4 plumed
    ```
    Specifically, the above commands do not install the PLUMED kernel library
    for you. You should compile it separately and export the `PLUMED_KERNEL`
    environment variable before actually perform your PLUMED aided MD runs.

## TESTS
First, install the `pytest` package with:
```bash
conda install pytest -c conda-forge
```
or
```bash
pip install pytest
```
Then, enter the `somd/tests` directory and invoke this command (you need to
change the `SIESTA_COMMAND` variable to the actual path of your `siesta`
binary):
```bash
SIESTA_COMMAND='/path/to/siesta' py.test
```

## USAGE
SOMD has a naive command line interface, which reads the
[TOML](https://toml.io/) format configure file. A typical input file looks
like this (which defines a NVT run of a water molecule):

```toml
[system]
        structure = "H2O.POSCAR"
[[group]]
        atom_list = "all"
        initial_temperature = 300.0
[[potential]]
        type = "SIESTA"
        siesta_options = """
        xc.functional          GGA
        xc.authors             PBE
        PAO.BasisSize          DZP
        Mesh.Cutoff            300 Ry
        """
        siesta_command = "mpirun -np 4 /path/to/siesta"
[[trajectory]]
        format = "H5"
        file_name = "traj.h5"
        interval = 10
[[logger]]
        format = "CSV"
        file_name = "data.csv"
        interval = 10
[integrator]
        type = "BAOAB"
        timestep = 0.0005
        temperatures = 300.0
        relaxation_times = 0.1
[run]
        n_steps = 500
```
Based on this file (e.g., it is called `input.toml`), you could run your
simulation via the following command:
```bash
somd -i input.toml
```
You may also invoke SOMD as a library and implement your own simulation
protocols. For example, the above configure file equals to the following
python script:
```python
import somd

siesta_command = 'mpirun -np 4 /path/to/siesta'
siesta_options = r"""
xc.functional          GGA
xc.authors             PBE
PAO.BasisSize          DZP
Mesh.Cutoff            300 Ry
"""

system = somd.core.systems.create_system_from_poscar('H2O.POSCAR')
g = {'atom_list': list(range(0, system.n_atoms)), 'has_translations': False}
system.groups.create_from_dict(g)
system.groups[0].add_velocities_from_temperature(300)
potential = somd.potentials.create_siesta_potential(system,
                                                    range(0, system.n_atoms),
                                                    siesta_options,
                                                    siesta_command)
system.potentials.append(potential)

integrator = somd.core.integrators.baoab_integrator(0.0005,
                                                    temperatures=[300],
                                                    relaxation_times=[0.1],
                                                    thermo_groups=[0])
trajectory = somd.apps.trajectories.H5WRITER('traj.h5', write_forces=False,
                                             interval=10)
logger = somd.apps.loggers.DEFAULTCSVLOGGER('data.csv', interval=10)
simulation = somd.apps.simulations.SIMULATION(system=system,
                                              integrator=integrator,
                                              trajectories=[trajectory],
                                              loggers=[logger])

simulation.run(500)
```
Based on this script (e.g., it is called `input.py`), you could run your
simulation via the following command:
```bash
python input.py
```

## DOCUMENTATION
A problem-oriented documentation could be found [here](doc/README.md).

## TUTORIALS
Tutorials of SOMD could be found
[here](https://www.github.com/initqp/somd_tutorials). Going through these
tutorials is considered as an efficient way to get familiar with SOMD.

## FAQ
- Q: There are millions of MD packages out there, why do you need another one?

  A: Because I write it for fine. If you are seeking for unity, you may want to
  use some mainstream packages like i-Pi, ASE, Tinker or even CP2K. Each of
  them provides excellent functionalities.

- Q: Will SOMD support other ab-initio packages, like VASP?

  A: Maybe. But I can not afford a VASP license, it is way luxurious for me.

- Q: Are you riding on the wave of machine learning potentials?

  A: I am. And I'm tired of pretending I'm not. 🤡

- Q: How to cite the code?

  A: You do not have to. But if your publisher forces you to do so, you may
  cite the GitHub repo directly.
