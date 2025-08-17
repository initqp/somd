# TOML Input References.

This documentation briefly introduces SOMD's TOML format input file. If you
are not familiar with TOML's grammar, you may want to read its
[documentation](https://toml.io/) first. Similar to the TOML format itself,
keys in SOMD's input files are case-sensitive (and are all in lower cases).
But the option strings (excepting the keys that define file names) are
case-insensitive.

## Table of Contents

1. [The `[system]` table.](#system)
2. [The `[[group]]` array.](#group)
3. [The `[[potential]]` array.](#potential)
4. [The `[constraints]` table.](#constraints)
5. [The `[integrator]` table.](#integrator)
6. [The `[barostat]` table.](#barostat)
7. [The `[[trajectory]]` array.](#trajectory)
8. [The `[[logger]]` array.](#logger)
9. [The `[run]` table.](#run)
10. [The `[evaluation]` table.](#evaluation)
11. [The `[[script]]` array.](#script)

## The `[system]` table. <a name="system"></a>

**If Mandatory**: yes

**Default Value**: None

**Descriptions**: This table defines the simulated system from a structure file.

**Keys**:

- **`structure`**

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Descriptions**: The initial structure file. This file could be in the
    [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) or the
    [PDB](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
    format. SOMD will read the atomic types, atomic positions and the simulation
    box from this file. Format of this file is identified by its extension name.
    If the file does not contain a extension name, you could set the `format`
    key in this table.

- **`format`**

    **If Mandatory**: no

    **Type**: `str`

    **Valid Values**: `"pdb"` or `"poscar"`

    **Default Value**: None

    **Descriptions**: The format of the structure file. If the structure file
    does not contain an extension name, this option would be useful. Once this
    key is defined, the extension name of the structure file will be ignored.

**Examples**:
```toml
[system]
        structure = "H2O.pdb"
        format = "pdb"
```
```toml
[system]
        structure = "../POSCAR"
```

## The `[[group]]` array. <a name="group"></a>

**If Mandatory**: no

**Default Value**: One atom group that corresponds to the whole simulated
system.

**Descriptions**: Each table in this array defines an atom group. You could
define multiple atom groups.

**Keys**:

- **`atom_list`**

    **If Mandatory**: yes

    **Type**: `List[int]` or `str`

    **Default Value**: None

    **Descriptions**: Indices of atoms that belong to this group. The indices
    may be represented by a list of integers, or a formatted string, which
    could represent a range of atoms conveniently. For example, the following
    two values of this key select same atoms:
    ```toml
    atom_list = [0, 1, 2, 5, 7, 8, 9]
    ```
    ```toml
    atom_list = "0:2,5,7:9"
    ```
    If you would like to define a group that contains all atoms in the system,
    you may invoke the following syntactic sugar:
    ```toml
    atom_list = "all"
    ```

    **Notes**: In SOMD, atom indices start from **ZERO**.

- **`has_translations`**

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `true`

    **Descriptions**: If this group has three translational degrees of freedom
    (DOFs). As we know, for a homogeneous-phase 3D periodic atomic system,
    the three total translational DOFs do not exchange energy with all other
    DOFs. And for an adsorption model, the substrate does not have total
    translational DOFs at all. As a result, when calculating temperatures and
    performing thermostating, these DOFs should be excluded. To this end, you
    should tell SOMD that weather the specified atom group contains
    translational DOFs, and this is done by defining the "has_translations"
    key. By setting this key to `false`, SOMD will subtract the number of DOFs
    of this group by three, and regularly remove its total translational
    momentums during the simulation (if your simulation undergoes a NVE
    ensemble, after removing the total translational momentums, velocities of
    atoms in this group will be scaled to conserve the total energy).

    **Notes**:
    1. If two atom groups contain same atoms, total translational DOFs of the
       two groups must not absent at the same time. For example, define the
       following two atom groups in the same input file will cause
       **AN ERROR**:
       ```toml
       [[group]]
               atom_list = [0, 1, 2, 3, 4, 5]
               has_translations = false
       [[group]]
               atom_list = [3, 4, 5, 6, 7, 8]
               has_translations = false
       ```
    2. If one group fully contains another group, and the smaller group does
       not have the total translational DOFs, the number of DOFs of the larger
       atom group will be modified as well. For example, under the following
       settings, number of DOFs of `GROUP0` will be 27 and number of DOFs of
       `GROUP1` will be 12.
       ```toml
       [[group]]
               label = "GROUP0"
               atom_list = "0:9"
       [[group]]
               label = "GROUP1"
               atom_list = [0, 1, 2, 3, 4]
               has_translations = false
       ```
    3. If you have not specified which atom group does not contain the
       translational DOFs, SOMD will automatically remove the translational
       DOFs of the whole simulated system.
    4. If you failed to understand above descriptions, you could safely leave
       this key alone, and SOMD will handle it for you.

- **`initial_temperature`**

    **If Mandatory**: no

    **Type**: `float`

    **Unit**: Kelvin

    **Default Value**: None

    **Descriptions**: Generate random initial velocities for atoms in this
    group according to the Boltzmann distribution under a given temperature.

    **Notes**: If you have also set the `initial_velocities` option, the
    random velocities generated by this key will be added to the velocities
    defined by that key.

- **`initial_velocities`**

    **If Mandatory**: no

    **Type**: `List[List[float]]`

    **Dimension**: (number of atoms in this group) * 3

    **Unit**: nanometer/picosecond

    **Default Value**: None

    **Descriptions**: The initial Cartesian velocities of each atom in this
    group.

    **Notes**: If you have also set the `initial_temperature` option, the
    random velocities generated by that key will be added to the velocities
    defined by this key.

- **`label`**

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: `GROUPN`, where `N` is the index of the group.

    **Descriptions**: A label of this group.

**Examples**:
```toml
[[group]]
        atom_list = "all"
        initial_temperature = 300.0
        label = "the_whole_system"
[[group]]
        atom_list = "0:100"
        has_translations = false
        label = "substrate"
[[group]]
        atom_list = [101, 102, 103, 104]
        initial_velocities = [
                [0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0]]
        label = "molecule"
```

## The `[[potential]]` array. <a name="potential"></a>

**If Mandatory**: yes

**Default Value**: None

**Descriptions**: Each table in this array defines a potential calculator. You
could define multiple calculators, but at least one potential calculator have
to be present.

**Keys**:

- **`type`**

    **If Mandatory**: yes

    **Type**: `str`

    **Valid Values**: `"siesta"`, `"dftd3"`, `"dftd4"`, `"tblite"`, `"nep"`, `"mace"` or `"plumed"`

    **Default Value**: None

    **Descriptions**: Type of the potential calculator:
    - **`"siesta"`**: The SIESTA ab-initio potential.
    - **`"dftd3"`**: Grimme's DFTD3 dispersion corrections.
    - **`"dftd4"`**: Grimme's DFTD4 dispersion corrections.
    - **`"tblite"`**: The tight binding potential calculated with the TBLite code.
    - **`"nep"`**: The neuroevolution potential.
    - **`"mace"`**: The E(3)-equivariant potentials based on the Atomic Cluster
      Expansion.
    - **`"plumed"`**: The PLUMED wrapper. Although a PLUMED calculation may be
      bias-free, you still need to define it here.

    **Required Keys**: Different potential calculator will require different
    other keys:
    - **`"siesta"`**: `atom_list`, `siesta_options`, `siesta_command`,
      `pseudopotential_dir`
    - **`"dftd3"`**: `atom_list`, `functional`, `damping`, `atm`
    - **`"dftd4"`**: `atom_list`, `functional`, `total_charge`, `atm`
    - **`"tblite"`**: `atom_list`, `functional`, `total_charge`, `total_spin`
    - **`"nep"`**: `atom_list`, `file_name`, `use_tabulating`
    - **`"mace"`**: `atom_list`, `file_name`, `device`, `virial`, `energy_unit`,
      `length_unit`, `compile_mode`, `compile_full_graph`, `total_charge`,
      `total_spin`
    - **`"plumed"`**: `file_name`

- **`atom_list`**

    **If Mandatory**: no

    **Type**: `List[int]` or `str`

    **Default Value**: All atoms in the simulated system.

    **Descriptions**: Indices of atoms that are covered by the potential. The
    indices may be represented by a list of integers, or a formatted string,
    which could represent a range of atoms conveniently. For example, the
    following two values of this key select same atoms:
    ```toml
    atom_list = [0, 1, 2, 5, 7, 8, 9]
    ```
    ```toml
    atom_list = "0:2,5,7:9"
    ```
    If you would like to define a group that contains all atoms in the system,
    you may invoke the following syntactic sugar:
    ```toml
    atom_list = "all"
    ```

    **Notes**: The PLUMED "potential" could only act on the whole simulated
    system.

- **`siesta_options`**

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Dependency**: `type = "siesta"`

    **Descriptions**: Options for a SIESTA calculation. This key contains the
    options you normally defined in an FDF format SIESTA input file, excepting
    the atomic positions, atomic types, cell vectors and the task type (these
    options are handled by SOMD).

- **`siesta_command`**

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Dependency**: `type = "siesta"`

    **Descriptions**: Path to the SIESTA binary. You may include the `mpirun`
    command or other statements in this key as well. For example:
    ```toml
    siesta_command = "OMP_NUM_THREADS=32 MKL_NUM_THREADS=32 /path/to/siesta"
    ```
    ```toml
    siesta_command = "mpirun -np 32 /path/to/siesta"
    ```

- **`pseudopotential_dir`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: The current working directory.

    **Dependency**: `type = "siesta"`

    **Descriptions**: The directory of the pseudopotential files. These files
    should be in the PSML, PSF or VPS formats, and the prefixes of the files
    must be element symbols (case-sensitive). If more than one pseudopotential
    file of an element is present, SOMD will invoke the pseudopotential files
    in the priority of PSML > PSF > VPS.

- **`functional`**:

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Dependency**: `type = "dftd3"`, `type = "dftd4"` or `type = "tblite"`

    **Descriptions**: Name of the functional which is invoked in the DFT
    calculations. For the `"tblite"` potential, this key is the calculation
    level, should be one of `"GFN2-xTB"`, `"GFN1-xTB"` and `"IPEA1-xTB"`.

- **`atm`**

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `false`

    **Dependency**: `type = "dftd3"` or `type = "dftd4"`

    **Descriptions**: If enable the three body correction.

- **`damping`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: `"ZeroDamping"`

    **Dependency**: `type = "dftd3"`

    **Descriptions**: The damping type of the DFTD3 potential.

- **`total_charge`**:

    **If Mandatory**: no

    **Type**: `int`

    **Unit**: a.u.

    **Default Value**: `0`

    **Dependency**: `type = "dftd4"` or `type = "tblite"` or `type = "mace"`

    **Descriptions**: The net charge of the simulated system.

- **`total_spin`**:

    **If Mandatory**: no

    **Type**: `int`

    **Unit**: a.u.

    **Default Value**: `0`

    **Dependency**: `type = "tblite"` or `type = "mace"`

    **Descriptions**: The net spin ($N_{alpha} - N_{beta}$) of the simulated
    system.

- **`file_name`**:

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Dependency**: `type = "nep"`, `type = "mace"` or `type = "plumed"`

    **Descriptions**: Path to the NEP potential file (`nep.txt`), path to
    the MACE model file or path to the PLUMED input file, depending on the
    value of the `type` key.

- **`use_tabulating`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: False

    **Dependency**: `type = "nep"`

    **Descriptions**: If invoke the tabulated version of NEP. This could speed
    up the calculation, read
    [this page](https://github.com/brucefan1983/NEP_CPU/pull/18) for details.

- **`device`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: `"cpu"`

    **Dependency**: `type = "mace"`

    **Descriptions**: Name of the device to use for evaluating the potential.

- **`virial`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `true`

    **Dependency**: `type = "mace"`

    **Descriptions**: If calculate the virial tensor. Disabling this option
    will speed up _NVT_ or _NVE_ runs.

- **`compile_mode`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: None

    **Dependency**: `type = "mace"`

    **Descriptions**: Mode of compiling the model (see torch.compile for details).

- **`compile_full_graph`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `False`

    **Dependency**: `type = "mace"`

    **Descriptions**: The `fullgraph` parameter of `torch.compile`.

- **`energy_unit`**:

    **If Mandatory**: no

    **Type**: `float`

    **Default Value**: `96.485`

    **Dependency**: `type = "mace"`

    **Descriptions**: The energy unit of the model outputs. In unit of (kJ/mol).
    E.g., if the your model trained with a dataset that invokes (eV) as energy
    units, this parameter should be set to 96.485.  The default value is capable
    for [most pretrained models provided by the MACE team](
    https://mace-docs.readthedocs.io/en/latest/examples/foundation_models.html
    ).

- **`length_unit`**:

    **If Mandatory**: no

    **Type**: `float`

    **Default Value**: `0.1`

    **Dependency**: `type = "mace"`

    **Descriptions**: The length unit of the model outputs. In unit of (nm).
    E.g., if the your model trained with a dataset that invokes (A) as length
    units, this parameter should be set to 0.1. The default value is capable for
    [most pretrained models provided by the MACE team](
    https://mace-docs.readthedocs.io/en/latest/examples/foundation_models.html
    ).

**Examples**:
```toml
[[potential]]
        type = "SIESTA"
        siesta_options = """
        xc.functional          GGA
        xc.authors             revPBE
        PAO.BasisSize          DZP
        Mesh.Cutoff            300 Ry
        Diag.Algorithm         ELPA-1stage
        SolutionMethod         diagon
        ElectronicTemperature  1 meV
        """
        siesta_command = "mpirun -np 4 /path/to/siesta"
        pseudopotential_dir = "./data"
[[potential]]
        type = "dftd3"
        functional = "revpbe"
[[potential]]
        type = "plumed"
        file_name = "./data/plumed.inp"
```
```toml
[[potential]]
        type = "nep"
        file_name = "./nep.txt"
```

## The `[constraints]` table. <a name="constraints"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: The holonomic constraints to be applied to the simulated
system. Constraints are achieved by the RATTLE algorithm in SOMD.

**Keys**:

- **`types`**

    **If Mandatory**: yes

    **Type**: `List[int]`

    **Default Value**: None

    **Descriptions**: Types of the constraints. Each element in this list
    stands for the type of one constraint, thus, the length of this list
    is the number of the constraints to be applied. The valid types are:
    - **`0`**: The distance between two atoms.
    - **`1`**: The angle between three atoms.
    - **`2`**: The dihedral angle between four atoms.

- **`indices`**

    **If Mandatory**: yes

    **Type**: `List[List[int]]`

    **Dimension**: (number of constraints) * (number of atoms per constraint)

    **Default Value**: None

    **Descriptions**: Indices of the atoms that take part in the constraints.
    Each element of this list is also a list, which includes the indices of
    the atoms that take part in the corresponding constraint. Thus, length of
    this list should be the same as the length of the `types` key. When the
    *i*-th element of the `type` key is `0`, `1` and `2`, the length of the
    *i*-th element of this list should be 2, 3, and 4.

- **`targets`**

    **If Mandatory**: yes

    **Type**: `List[float]`

    **Dimension**: number of constraints

    **Default Value**: None

    **Descriptions**: Target values of the constraints. Length of this list
    should be the same as the length of the `types` key.  When the *i*-th
    element of the `type` key is `0`, `1` and `2`, the unit of the *i*-th
    element of this list should be nanometer, radian and radian.

- **`tolerances`**

    **If Mandatory**: no

    **Type**: `List[float]`

    **Dimension**: number of constraints

    **Default Value**: `1E-14` for every constraints.

    **Descriptions**: Convergence tolerances of the constraints. Length of
    this list should be the same as the length of the `types` key.When the
    *i*-th element of the `type` key is `0`, `1` and `2`, the unit of the
    *i*-th element of this list should be nanometer, radian and radian.

- **`max_cycles`**

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: `300`

    **Descriptions**: Maximum number of the RATTLE iterations.

**Examples**:
```toml
[constraints]
        types = [0, 1, 2]
        indices = [[4, 7], [0, 4, 7], [2, 0, 4, 7]]
        targets = [0.108, 1.92527, -0.9995]
        tolerances = [1E-14, 1E-14, 1E-14]
        max_cycles = 100
```

## The `[integrator]` table. <a name="integrator"></a>

**If Mandatory**: yes

**Default Value**: None

**Descriptions**: The integrator to propagate the simulated system.

**Keys**:

- **`type`**

    **If Mandatory**: yes

    **Type**: `str`

    **Valid Values**: `"vv"`, `"cs4"`, `"nhc"`, `"baoab"`, `"obabo"`,
    `"gbaoab"` or `"gobabo"`

    **Default Value**: None

    **Descriptions**: Type of the integrator:
    - **`"vv"`**: The Velocity Verlet integrator.
    - **`"cs4"`**: Calvo and Sanz-Serna's fourth-order integrator.
    - **`"nhc"`**: The Nose-Hoover Chains integrator
    ($N_{RESPA} = 4, N_{SY} = 6, N_{chains} = 6$).
    - **`"baoab"`**: The Langevin integrator with a BAOAB splitting (can not
    be used with constraints).
    - **`"obabo"`**: The Langevin integrator with a OBABO splitting (can not
    be used with constraints).
    - **`"gbaoab"`**: The Geodesic Langevin integrator with a BAOAB splitting
    (can be used with constraints, identical to the `"baoab"` integrator when
    not using constraints).
    - **`"gobabo"`**: The Geodesic Langevin integrator with a OBABO splitting.
    (can be used with constraints, identical to the `"obabo"` integrator when
    not using constraints).

    **Notes**:
    1. The `"vv"` and `"cs4"` integrators are *NVE* integrator and other
    integrators are thermalized.
    2. The `"cs4"` integrator is four times more expensive than the trivial
    `"vv"` integrator, but it could generate much more accurate trajectories.
    3. General pros of the Langevin integrators: the equilibrium is faster to
    reach, and the partition of kinetic energies is more uniform.
    4. General cons of the Langevin integrators: every kind of Langevin
    integrator violates the underlying dynamics of the simulated system. As a
    result, when collecting dynamical observations (e.g., correlation
    functions), you should avoid using them. If thermostating is required, use
    the `"nhc"` integrator. Otherwise, use any *NVE* integrator.
    5. The `"baoab"` integrator shows better performances in configuration
    space sampling tasks (e.g., distribution function or free energy
    calculations), although the kinetic energy temperature may be a little
    lower the given value (which is normal, since the definition of
    "temperature" is not unique). We strongly recommend you to use this
    integrator when performing such tasks.
    6. **DO NOT** use the `"nhc"` integrator in simulated annealing
    simulations. It will introduce serious oscillations.

- **`timestep`**

    **If Mandatory**: yes

    **Type**: `float`

    **Default Value**: None

    **Unit**: picosecond

    **Descriptions**: The timestep of the integrator.

    **Notes**: Timesteps in SOMD could be negative, which means propagating
    the system backwardly.

- **`thermo_groups`**

    **If Mandatory**: no

    **Type**: `List[int]` or `int`

    **Default Value**: The atom group that corresponds to the whole simulated
    system.

    **Dependency**: Thermalized integrators.

    **Descriptions**: The indices of the [atom groups](#group) to thermalize.
    Each group will be coupled to a different thermostat. When there is only
    one group to thermalize, you could define this key by the index of the
    group instead of a single-element list, e.g.,
    ```toml
    thermo_groups = 0
    ```
    is equal to
    ```toml
    thermo_groups = [0]
    ```

- **`temperatures`**

    **If Mandatory**: yes

    **Type**: `List[float]` or `float`

    **Default Value**: None

    **Dependency**: Thermalized integrators.

    **Descriptions**: Temperatures of the thermostats. When there is only
    one group to thermalize, you could define this key by the temperature of
    that group instead of a single-element list. Otherwise, length of this list
    should be the same as the length of the `thermo_groups` key.

- **`relaxation_times`**

    **If Mandatory**: yes

    **Type**: `List[float]` or `float`

    **Default Value**: None

    **Dependency**: Thermalized integrators.

    **Descriptions**: Relaxation timescales of the thermostats. When there is
    only one group to thermalize, you could define this key by the relaxation
    timescale of that group instead of a single-element list. Otherwise, length
    of this list should be the same as the length of the `thermo_groups` key.

- **`splitting`**

    **If Mandatory**: no

    **Type**: `List[dict]`

    **Default Value**: None

    **Descriptions**: The splitting of the integrator. Read the
    `somd/core/integrators.py` source file for details. When this key is
    defined, you should not define the `type` key.

**Examples**:
```toml
[integrator]
        type = "baoab"
        timestep = 0.001
        temperatures = 300.0
        relaxation_times = 0.05
```
```toml
[integrator]
        type = "nhc"
        timestep = 0.001
        thermo_groups = [0, 1]
        temperatures = [300.0, 400.0]
        relaxation_times = [0.05, 0.05]
```
```toml
# The BAOAB Langevin integrator could also be defined by this way:
[integrator]
        timestep = 0.001
        temperatures = 300.0
        relaxation_times = 0.05
[[integrator.splitting]]
        operators = ["V", "R", "O", "R", "V"]
```
```toml
# The "middle point" NHC integrator:
[integrator]
        timestep = 0.001
        temperatures = 300.0
        relaxation_times = 0.05
[[integrator.splitting]]
        operators = ["V", "R", "N", "R", "V"]
```
```toml
# The Geodesic OBABO Langevin integrator could also be defined by this way:
[integrator]
        timestep = 0.001
        temperatures = 300.0
        relaxation_times = 0.05
[[integrator.splitting]]
        operators = ["O", "Cv"]
[[integrator.splitting]]
        operators = ["V", "Cv"]
[[integrator.splitting]]
        operators = ["CR", "R", "Cv"]
        repeating = 5
[[integrator.splitting]]
        operators = ["V", "Cv"]
[[integrator.splitting]]
        operators = ["O", "Cv"]
```

## The `[barostat]` table. <a name="baostat"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: The Berendsen-type barostat.

**Note**: By now, SOMD only supports the original Berendsen barostat, which
could not generate correct isothermal-isobaric distributions. Thus, do not
rely on any distribution generated by this barostat (but equilibrium volumes
are still reliable). Besides, the barostat can not be applied with constraints.

**Keys**:

- **`pressures`**

    **If Mandatory**: yes

    **Type**: `List[float]` or `float`

    **Dimension**: 1 or 6

    **Unit**: megapascal

    **Default Value**: None

    **Descriptions**: The pressures of the barostat. If length of this list
    is six, an anisotropy pressure controlling will be applied. Otherwise,
    an isotropy pressure controlling will be applied. And when using isotropy
    pressure controlling, you could define this key by the one hydrostatic
    pressure instead of a single-element list, e.g.,
    ```toml
    pressures = 0.1
    ```
    is equal to
    ```toml
    pressures = [0.1]
    ```

    **Notes**: When using anisotropy pressure controlling, the Cartesian
    directions of the six target pressures are: xx yy zz xy xz yz.

- **`beta`**

    **If Mandatory**: yes

    **Type**: `List[float]` or `float`

    **Dimension**: 1 or 6

    **Unit**: 1/megapascal

    **Default Value**: None

    **Descriptions**: The isothermal compressibilities of the simulated
    system. Length of this list should be the same as the length of the
    `pressures` key. And when using isotropy pressure controlling, you could
    define this key by the one isotropy compressibility instead of a
    single-element list.

    **Notes**: When using anisotropy pressure controlling, the Cartesian
    directions of the six compressibilities are: xx yy zz xy xz yz.

- **`relaxation_time`**

    **If Mandatory**: yes

    **Type**: `float`

    **Unit**: picosecond

    **Default Value**: None

    **Descriptions**: The relaxation timescale of the barostat.

**Examples**:
```toml
[barostat]
        pressures = 1.0
        beta = 1.0E-5
        relaxation_time = 0.1
```
```toml
[barostat]
        pressures = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        beta = [1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5]
        relaxation_time = 0.1
```

## The `[[trajectory]]` array. <a name="trajectory"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: Each table in this array defines a trajectory writer. You
could define multiple trajectory writers.

**Keys**:

- **`format`**:

    **If Mandatory**: no

    **Type**: `str`

    **Valid Values**: `"h5"` or `"exyz"`

    **Default Value**: `h5`

    **Descriptions**: Format of the trajectory file. The `"h5"` format obeys
    [MDTraj's HDF5 convention](https://mdtraj.org/1.9.7/hdf5_format.html), and
    the `"exyz"` format obeys
    [GPUMD's EXYZ convention](https://gpumd.org/nep/input_files/train_test_xyz.html).

    **Notes**: To view the `h5` format trajectories, you could use the
    `mdconvert` script provided the `mdtraj` package to convert it to other
    formats, e.g., the extensively supported `.nc` format:
    ```bash
    mdconvert prefix.trajectory.h5 -o traj.nc
    ```

- **`prefix`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: The system label.

    **Descriptions**: The prefix of the trajectory file.

    **Notes**: When the value of the `format` key is `"h5"`, the generated
    trajectory file will be named `prefix.trajectory.h5`, and when the value
    of the `format` key is `"exyz"`, the generated trajectory file will be
    named `prefix.trajectory.xyz`.

- **`interval`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: `1`

    **Descriptions**: The interval between two updates of the trajectory file.

- **`write_forces`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `false`

    **Descriptions**: If record forces in the trajectory file.

- **`write_velocities`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `true` for `h5` trajectories, `false` for `exyz` trajectories

    **Descriptions**: If record velocities in the trajectory file.

- **`wrap_positions`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `false`

    **Descriptions**: If wrap atoms to the box before writing the trajectory.

- **`potential_list`**:

    **If Mandatory**: no

    **Type**: `List[int]`

    **Default Value**: Indices of all potential calculators.

    **Descriptions**: Indices of the potential calculators that contribute to
    the recorded potential energies, virials, and forces. This is extremely
    useful when constructing training and testing set of neuroevolution
    potentials. For example, when you are performing a biased enhanced sampling
    with PLUMED, and want to invoke the generated trajectory to train NEPs. The
    following settings will ensure that only the ab-initio data will be
    written to the `"exyz"` file:
    ```toml
    [[potential]]
            type = "SIESTA"
            siesta_options = """
            xc.functional          GGA
            xc.authors             revPBE
            PAO.BasisSize          DZP
            """
            siesta_command = "mpirun -np 4 /path/to/siesta"
    [[potential]]
            type = "dftd3"
            functional = "revpbe"
    [[potential]]
            type = "plumed"
            file_name = "plumed.inp"
    [[trajectory]]
            format = "exyz"
            # We need force data to train a NEP.
            write_forces = true
            # Only write the forces, virials and energies generated by the
            # first two potential calculators to the trajectory.
            potential_list = [0, 1]
    ```

- **`use_float64`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `false`

    **Dependency**: `format = "h5"`

    **Descriptions**: If use 64-bit floating point data in the trajectory
    instead of the original 32-bit convention.

    **Notes**: Enabling this option will lead to non-standard trajectories and
    may cause some readers to fail. Besides, sizes of the trajectories will
    largely be increased.

- **`energy_shift`**:

    **If Mandatory**: no

    **Type**: `float`

    **Default Value**: `0.0`

    **Dependency**: `format = "exyz"`

    **Descriptions**: Shift the total energy by this value before recording the
    total energy to the trajectory (this value could be set to the potential
    energy of the initial conformation, which is usually negative). In unit of
    (kJ/mol). This is useful when generating training sets of NEP.

**Examples**:
```toml
[[trajectory]]
        format = "h5"
        interval = 10
        wrap_positions = true
        write_velocities = true
[[trajectory]]
        format = "exyz"
        prefix = "train"
        write_forces = true
        potential_list = [0, 1]
```

## The `[[logger]]` array. <a name="logger"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: Each table in this array defines a system data writer. You
could define multiple writers. These files will record system data, e.g.,
total potential energies, temperature and kinetic energies of each atom groups,
pressure tensors, etc.

**Keys**:

- **`format`**:

    **If Mandatory**: no

    **Type**: `str`

    **Valid Values**: `"csv"` or `"txt"`

    **Default Value**: `h5`

    **Descriptions**: Format of the data file.

- **`prefix`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: The system label.

    **Descriptions**: The prefix of the data file.

    **Notes**: When the value of the `format` key is `"csv"`, the generated
    data file will be named `prefix.data.h5`, and when the value of the
    `format` key is `"txt"`, the generated data file will be named
    `prefix.data.txt`.

- **`interval`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: `1`

    **Descriptions**: The interval between two updates of the data file.

- **`potential_list`**:

    **If Mandatory**: no

    **Type**: `List[int]`

    **Default Value**: Indices of all potential calculators.

    **Descriptions**: Indices of the potential calculators that contribute to
    the recorded potential energies.

**Examples**:
```toml
[[logger]]
        format = "csv"
        interval = 5
```

## The `[run]` table. <a name="run"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: This table defines the simulation information.

**Notes**: This table is not mandatory only when you have defined the
`[evaluation]` table. Otherwise, this table **IS MANDATORY**.

**Keys**:

- **`n_steps`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: None

    **Descriptions**: Number of the simulation steps to run.

    **Notes**: Either this key or the `n_seconds` key should be defined.

- **`n_seconds`**:

    **If Mandatory**: no

    **Type**: `int`

    **Unit**: Second

    **Default Value**: None

    **Descriptions**: Length of time to run. This is useful when you have a
    limited amount of computer time available, and want to run the longest
    simulation possible in that time. This method will continue taking time
    steps until the specified clock time has elapsed, then exit. It also
    automatically writes out a restart file (named `system_label.restart.h5`)
    before exiting.

    **Notes**: Either this key or the `n_steps` key should be defined.

- **`seed`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: None

    **Descriptions**: The seed value of the random number generator in SOMD.

- **`label`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: The prefix of the input file name.

    **Descriptions**: The default prefix of all output file generated by SOMD.

- **`restart_from`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: None

    **Descriptions**: Name of the restart file. If this key is defined, SOMD
    will continue the simulation from the state recorded by the restart file.

    **Notes**: When restarting, SOMD will try to append new data to the old
    output files. Otherwise, SOMD will back up old output files.

**Examples**:
```toml
[run]
        seed = 1
        n_steps = 1000
        label = "run"
```
```toml
[run]
        n_steps = 1000
        restart_from = "run.restart.h5"
```

## The `[evaluation]` table. <a name="evaluation"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: This table defines a evaluation task, which calculates
various energies of conformations recorded in a `hdf5` trajectory file.
The resulted energies, forces and other data will be written to the trajectory
and log files defined in the `[[trajectory]]` and `[[logger]]` array.

**Notes**: By default, you do not need to define an `[integrator]` for the
evaluation task. Besides, The `[[script]]` array works with the evaluation task.

**Keys**:

- **`file_name`**:

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Descriptions**: Name of the input `hdf5` trajectories.

- **`interval`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: 1

    **Descriptions**: The interval of reading the input trajectory.

- **`die_on_fail`**:

    **If Mandatory**: no

    **Type**: `bool`

    **Default Value**: `true`

    **Descriptions**: If skip the conformations which cause the evaluation of
    the potentials to fail.

**Examples**:
```toml
[evaluation]
    file_name = "./run.trajectory.h5"
    interval = 1
```
The below example is a complete input file for calculating potential energies
of the conformations recorded in the file `run.trajectory.h5`, using the `nep`
potential file `nep.txt`.
```toml
[system]
    structure = "topo.pdb"
[[potential]]
    type = "nep"
    file_name = "nep.txt"
[[trajectory]]
    format = "h5"
    write_forces = true
    write_velocities = false
    interval = 1
[evaluation]
    file_name = "run.trajectory.h5"
    interval = 1
```

## The `[[script]]` array. <a name="script"></a>

**If Mandatory**: no

**Default Value**: None

**Descriptions**: Each table in this array defines a python function to be
executed in a given period.

**Notes**: The functions will be called at the end of the timestep.

**Keys**:

- **`update`**:

    **If Mandatory**: yes

    **Type**: `str`

    **Default Value**: None

    **Descriptions**: The function to be called regularly. This function takes
    only one parameter, which is the integrator instance of the simulation.
    You could use the following commands to find out the valid properties:
    ```python
    import somd
    help(somd.core.integrators.INTEGRATOR)
    ```

    **Notes**: This function must be named `update`.

- **`initialize`**:

    **If Mandatory**: no

    **Type**: `str`

    **Default Value**: None

    **Descriptions**: The initialization function. This function will be called
    only once at the starting of the simulation. This function takes no
    arguments.

    **Notes**: This function must be named `initialize`.

- **`interval`**:

    **If Mandatory**: no

    **Type**: `int`

    **Default Value**: `1`

    **Descriptions**: The interval of calling the `update` function.

**Examples**:
```toml
[[script]]
        # Dump the memory information to the `mem.log` file every 50 steps.
        interval = 50
        update = """def update(integrator):
        import os
        os.system("free -h >> mem.log")
        """
```
```toml
[[script]]
        # Rising the temperature to 500 K in 500 steps.
        interval = 5
        update = """def update(integrator):
        s = integrator.step
        if (s < 500):
                integrator.temperatures[0] += 5
        """
```
```toml
[[script]]
        # Back up the `DM` files generated by SIESTA calculations every 5 steps.
        interval = 5
        initialize = """def initialize():
        import os
        import shutil
        try:
                shutil.rmtree('./DM_files')
        except:
                pass
        os.mkdir('./DM_files')
        """
        update = """def update(integrator):
        import shutil
        s = integrator.step
        for p in integrator.system.potentials:
                if (p.__class__.__name__ == 'SIESTA'):
                        old_dm = p.working_directory + '/somd_tmp.DM'
                        new_dm = './DM_files/{:d}.DM'.format(s)
                        shutil.copy(old_dm, new_dm)
        """
```
