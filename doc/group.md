# Add atom groups to your system.

Atom group is one of the most important concepts in SOMD. The calculations of
many quantities are based on the definitions of atom groups, e.g. the degrees
of freedom (DOFs) and temperatures. In SOMD, atom groups are represented
by the `somd.core.groups.ATOMGROUP` class. All atom groups that belong to a
simulated system (an instance of the `somd.core.systems.MDSYSTEM` class) are
stored in its `groups` property. The `groups` property itself, meanwhile,
is an instance of the `somd.core.groups.ATOMGROUPS` class, which is inherited
from python's built-in list class.

You have two ways to add atom groups to a simulated system:
- Get an instance of the `somd.core.groups.ATOMGROUP` class and `append` it to
  your system:
  ```python
  >>> system = somd.core.systems.MDSYSTEM(10)
  >>> group = somd.core.groups.ATOMGROUP(system, list(range(0, 10)))
  >>> system.groups.append(group)
  ```
  The `somd.core.groups.ATOMGROUP` class takes at least two parameters: first,
  a `somd.core.systems.MDSYSTEM` instance which represents the simulated system
  that contains the group; second, a `list` that specifies the indices of atoms
  in this group. In above example, we divided all atoms to the generated group,
  obviously.

  You could also specify a name of the generated atom group, e.g., `'G1'`:
  ```python
  >>> system = somd.core.systems.MDSYSTEM(10)
  >>> group = somd.core.groups.ATOMGROUP(system, list(range(0, 10)), label='G1')
  >>> system.groups.append(group)
  ```

- The above example is equal to invoke the `create_from_dict` method of
  `system.groups`:
  ```python
  >>> system = somd.core.systems.MDSYSTEM(10)
  >>> system.groups.create_from_dict({'atom_list': list(range(0, 10)),
                                      'label': 'G1'})
  ```
  This method takes a dictionary as its parameter, where you could define the
  properties of the atom group.

After adding one atom group to your system, the `system.groups` property should
have a length of 1:
```python
>>> print(len(system.groups))
# 1
```
And as mentioned, the i-th object stored in `system.groups` contains
information about the i-th atom group, e.g., the number of atoms
(`system.groups[i].n_atoms`), the number of degrees of freedom
(`system.groups[i].n_dof`), the group temperature
(`system.groups[i].temperature`) and velocities of atoms in this group
(`system.groups[i].velocities`). You could execute the
`import somd; help(somd.core.groups.ATOMGROUP)` command in a `python` REPL to
view a full list of valid attributes. But among all the attributes, I would
like to detail the important and frequently-used `has_translations` property.

As we know, for a homogeneous-phase 3D periodic atomic system, the three total
translational DOFs do not exchange energy with all other DOFs. And for an
adsorption model, the substrate does not have total translational DOFs at all.
As a result, when calculating temperatures and performing thermostatings,
these DOFs should be excluded. To this end, you must tell SOMD that weather the
specified atom groups contain translational DOFs, and this is done through the
`has_translations` property of each atom group. By setting this property to
`False`, SOMD will subtract the number of DOFs of this group by three, and
regularly remove its total translational momentums during the simulation (if
your simulation underwent a NVE ensemble, after removing the total
translational momentums, velocities of atoms in this group will be scaled to
conserve the total energy). You could also specify this property when adding
new groups:
```python
>>> system = somd.core.systems.MDSYSTEM(10)
>>> system.groups.create_from_dict({'atom_list': list(range(0, 10)),
                                    'has_translations': False,
                                    'label': 'G1'})
```
**Note**: if two atom groups contain same atoms, they must not lack total
translational DOFs at the same time. For example, the following settings will
cause **an error**:
```python
>>> system = somd.core.systems.MDSYSTEM(10)
>>> system.groups.create_from_dict({'atom_list': list(range(0, 7)),
                                    'has_translations': False,
                                    'label': 'G1'})
>>> system.groups.create_from_dict({'atom_list': list(range(6, 10)),
                                    'has_translations': False,
                                    'label': 'G2'})
# RuntimeError: Atom group 'G2' is overlapping with atomic group 'G1' while
# both are binding with COM motion removers!
```
Besides, if one group fully contains another group, and the smaller group does
not have the total translational DOFs, the number of DOFs of the larger atom
group will be updated as well:
```python
>>> system = somd.core.systems.MDSYSTEM(10)
>>> system.groups.create_from_dict({'atom_list': list(range(0, 10)),
                                    'label': 'G1'})
>>> print(system.groups[0].n_dof)
# 30
>>> system.groups.create_from_dict({'atom_list': list(range(6, 10)),
                                    'has_translations': False,
                                    'label': 'G2'})
>>> print(system.groups[1].n_dof)
# 9
>>> print(system.groups[0].n_dof)
# 27
```

In contrast, SOMD does not support the "frozen" groups, where positions of
every atom in that group are completely fixed. This is because, regarding the
sizes of systems that could be investigated by AIMD, the "frozen" groups
usually introduce very strong artificiality. Instead, you should try to use
larger models to avoid large structural fluctuations, e.g., use thicker
subtracts.

To initialize the velocities of a given atom group, you may use the
`add_velocities_from_temperature` method of the group objects. For example, if
you want to initialize the velocities of the system under 300 K:
```python
>>> system = somd.core.systems.create_system_from_pdb('H2O.pdb')
>>> system.groups.create_from_dict({'atom_list': list(range(0, 3)),
                                    'has_translations': False})
>>> system.groups[0].add_velocities_from_temperature(300)
```
You could also modify the `velocities` property of the group objects to change
their velocities. For example, when simulating the collision between a
molecule and a surface, you could initialize the velocities of the atoms in
the substrate according to a given temperature, and specify the velocities of
the molecules directly (note the units of the velocities are nm/ps):
```python
>>> system = somd.core.systems.create_system_from_poscar('POSCAR')
>>> system.groups.create_from_dict({'atom_list': list(range(0, 100)),
                                    'has_translations': False,
                                    'label': 'substrate'})
>>> system.groups.create_from_dict({'atom_list': list(range(100, 110)),
                                    'label': 'molecule'})
>>> system.groups[0].add_velocities_from_temperature(1000)
>>> system.groups[1].velocities = [0, 0, -100]
```
(Above production is not strict of course, because we did not consider the
internal vibrational kinetic energies of the molecule.)
