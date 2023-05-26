# Output File Reference.

This documentation briefly introduces SOMD's private output files.

## Table of Contents

1. [The `prefix.active_learning.h5` file.](#active_learning)

## The `prefix.active_learning.h5` file. <a name="active_learning"></a>

**Task**: `[active_learning]`

**Format**: `HDF5`

**Descriptions**: This file records the progress and relative data
of an active learning simulation.

**Structure**:
```bash
/
├── created_time      # (attr, (1,), str)
├── program           # (attr, (1,), str)
├── root_directory    # (attr, (1,), str)
├── title             # (attr, (1,), str)
├── version           # (attr, (1,), str)
├── working_directory # (attr, (1,), str)
└── iteration_data    # (group)
    ├── 0             # (group)
    │   ├── accepted_structures         # (attr, (1,), str)
    │   ├── initialized                 # (attr, (1,), bool)
    │   ├── invoked_nep                 # (attr, (1,), str)
    │   ├── pre_training                # (attr, (1,), bool)
    │   ├── system_data                 # (attr, (1,), str)
    │   ├── training_sets               # (attr, (unlimited,), str)
    │   ├── visited_structures          # (attr, (1,), str)
    │   ├── working_directory           # (attr, (1,), str)
    │   ├── accepted_structure_indices  # (dataset, (unlimited,), int)
    │   ├── candidate_structure_indices # (dataset, (unlimited,), int)
    │   ├── force_msd                   # (dataset, (unlimited,), int)
    │   │   └── units                   # (attr, (1,), str)
    │   ├── force_msd_lower_limit       # (dataset, (unlimited,), int)
    │   │   └── units                   # (attr, (1,), str)
    │   ├── force_msd_upper_limit       # (dataset, (1,), int)
    │   │   └── units                   # (attr, (1,), str)
    │   ├── max_force_msd               # (dataset, (1,), int)
    │   │   └── units                   # (attr, (1,), str)
    │   ├── max_new_structures_per_iter # (dataset, (1,), int)
    │   ├── min_new_structures_per_iter # (dataset, (1,), int)
    │   ├── n_accepted_structures       # (dataset, (1,), int)
    │   ├── n_accurate_structures       # (dataset, (1,), int)
    │   ├── n_candidate_structures      # (dataset, (1,), int)
    │   ├── n_failed_structures         # (dataset, (1,), int)
    │   ├── n_untrained_structures      # (dataset, (1,), int)
    │   ├── n_visited_structures        # (dataset, (1,), int)
    │   └── progress                # (group)
    │       ├── ab_initial_finished # (attr, (1,), bool)
    │       └── training_finished   # (attr, (unlimited,), List[bool])
    ├── 1
    │   ├── ...
    │   ...
    ...
```
