# CLI References.

This documentation briefly introduces SOMD's command line interface (CLI).

## Table of Contents

1. [Execution modes of the `somd` command.](#modes)
2. [Keywords of  the `run` mode.](#run)
3. [Keywords of  the `select` mode.](#select)

## Execution modes of the `somd` command. <a name="modes"></a>

After installing the SOMD package, the command `somd` will be available.
You could type `somd -v` to check the package version, and `somd -h` for the
help messages.

In order to execute a SOMD operation through the CLI, you will need to specify
the require execution mode and the corresponding keywords through command:
```bash
somd [mode] --keyword_1 --keyword_2 ...
```
The `[mode]` string has two valid options, which are `run` and `select`.
The `run` mode reads reads a [`TOML` file](./toml.md) and does the dynamics,
and the `select` model selects structures for the active learning task.

## Keywords of  the `run` mode. <a name="run"></a>

This mode taks the following keywords.

- **`-i`/`--input`**

    **If Mandatory**: yes

    **Type**: `str`

    **Number of Required Parameters**: 1

    **Default Value**: None

    **Descriptions**: Name of the input `TOML` file.

## Keywords of  the `select` mode. <a name="select"></a>

This mode taks the following keywords.

- **`-t`/`--trajectories`**

    **If Mandatory**: yes

    **Type**: `List[str]` (delimited by white spaces)

    **Number of Required Parameters**: more than 1

    **Default Value**: None

    **Descriptions**:  Paths to the input HDF5 trajectories, with forces recorded.

- **`-o`/`--output`**

    **If Mandatory**: yes

    **Type**: `str`

    **Number of Required Parameters**: 1

    **Default Value**: None

    **Descriptions**: Name of the output trajectory file.

- **`--msd_f_min`**

    **If Mandatory**: no

    **Type**: `float`

    **Number of Required Parameters**: 1

    **Default Value**: 50.0

    **Descriptions**: The lower boundary of the force MSD used for identifying
    candidate structures. In unit of (kJ/mol/nm).

- **`--msd_f_max`**

    **If Mandatory**: no

    **Type**: `float`

    **Number of Required Parameters**: 1

    **Default Value**: 250.0

    **Descriptions**: The upper boundary of the force MSD used for identifying
    candidate structures. In unit of (kJ/mol/nm).

- **`--n_max`**

    **If Mandatory**: no

    **Type**: `int`

    **Number of Required Parameters**: 1

    **Default Value**: `inf`

    **Descriptions**: The max number of structures to select (default: inf).

- **`--log_file`**

    **If Mandatory**: no

    **Type**: `str`

    **Number of Required Parameters**: 1

    **Default Value**: `somd_selection.json`

    **Descriptions**: The name of the log file that contains the selection
    information.

- **`--shuffle`**

    **If Mandatory**: no

    **Number of Required Parameters**: 0

    **Default Value**: `True`

    **Descriptions**: If perform randomly selection if the number of candidate
    structures is larger than `n_max`.
