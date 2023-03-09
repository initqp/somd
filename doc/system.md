# Set up the simulated system.

Tile now, simulated systems in SOMD are represented by the
`somd.core.systems.MDSYSTEM` class. You could initialize a system containing
ten atoms by:
```python
>>> system = somd.core.systems.MDSYSTEM(10)
```
You may also assign a label for your system during the initialization:
```python
>>> system = somd.core.systems.MDSYSTEM(10, label='SYSTEM')
```
The generated `system` object contains information about a typical atomic
system, e.g., the atomic positions (`system.positions`), the atomic velocities
(`system.velocities`) and the cell vectors (`system.box`). You could execute
the `import somd; help(somd.core.systems.MDSYSTEM)` command in a `python` REPL
to view a full list of valid attributes.

Specifically, the `somd.core.systems` module provides two useful methods to
initialize the simulated systems from files:
- The `create_system_from_poscar` method could perform system initialization
  using the [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files. For
  example:
  ```python
  >>> system = somd.core.systems.create_system_from_poscar('H2O.POSCAR')
  ```
- The `create_system_from_pdb` method could perform system initialization
  using the [PDB](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
  files. For example:
  ```python
  >>> system = somd.core.systems.create_system_from_pdb('H2O.pdb')
  ```

Both method reads atomic types, atomic positions and simulation box data from
corresponding files, then initializes a `MDSYSTEM` instance based on the read
data.
