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
├── created_time      # (attr, (1,), str) Time at which the file was created.
├── program           # (attr, (1,), str) Name of the program that create this file. Should be SOMD.
├── root_directory    # (attr, (1,), str) Path of the root directory.
├── title             # (attr, (1,), str) Name of this file.
├── version           # (attr, (1,), str) Version number of SOMD.
├── working_directory # (attr, (1,), str) Path of the working directory.
└── iteration_data    # (group) The main data group that stores the data of each training iteration.
    ├── 0             # (group) Name of the iteration, which is simply the index of the iteration.
    │   ├── initialized                 # (attr, (1,), bool) If this iteration has been initialized (for internally usage only).
    │   ├── pre_training                # (attr, (1,), bool) If this iteration is the pre-training iteration (for internally usage only).
    │   ├── working_directory           # (attr, (1,), str) Path of the working directory of this iteration.
    │   ├── accepted_structures         # (attr, (1,), str) Name of the EXYZ file that stores the ab-initio calculation results of the accepted structures.
    │   ├── invoked_nep                 # (attr, (1,), str) Name of the NEP file that was selected to propagate the simulations in this iteration.
    │   ├── system_data                 # (attr, (1,), str) Name of the file that records the simulation data in this iteration.
    │   ├── visited_structures          # (attr, (1,), str) Name of the file that records the visited structures in this iteration.
    │   ├── training_sets               # (attr, (unlimited,), str) Name of the training sets that were used to train the NEPs in this iteration.
    │   ├── n_visited_structures        # (dataset, (1,), int) Number of visited structures.
    │   ├── n_candidate_structures      # (dataset, (1,), int) Number of candidate structures.
    │   ├── n_accepted_structures       # (dataset, (1,), int) Number of accepted structures.
    │   ├── n_accurate_structures       # (dataset, (1,), int) Number of accurate structures.
    │   ├── n_failed_structures         # (dataset, (1,), int) Number of failed structures.
    │   ├── n_untrained_structures      # (dataset, (1,), int) Number of untrained structures (for internally usage only).
    │   ├── accepted_structure_indices  # (dataset, (unlimited,), int) Indices of the accepted structures.
    │   ├── candidate_structure_indices # (dataset, (unlimited,), int) Indices of the candidate structures.
    │   ├── max_new_structures_per_iter # (dataset, (1,), int) The maximum number of accepted structures.
    │   ├── min_new_structures_per_iter # (dataset, (1,), int) The minimum number of accepted structures.
    │   ├── force_msd                   # (dataset, (unlimited,), int) MSD of the NEP forces of each visited structure.
    │   │   └── units                   # (attr, (1,), str)
    │   ├── max_force_msd               # (dataset, (1,), int) The max value of the force MSDs.
    │   │   └── units                   # (attr, (1,), str)
    │   ├── force_msd_lower_limit       # (dataset, (1,), int) The lower limit of candidate structures' force MSD.
    │   │   └── units                   # (attr, (1,), str)
    │   ├── force_msd_upper_limit       # (dataset, (1,), int) The upper limit of candidate structures' force MSD.
    │   │   └── units                   # (attr, (1,), str)
    │   ├── energy_shift                # (dataset, (1,), int) The potential energy shift value.
    │   │   └── units                   # (attr, (1,), str)
    │   ├── accepted_structure_energies # (dataset, (unlimited,), int) Relative potential energies of the accepted structures.
    │   │   └── units                   # (attr, (1,), str)
    │   └── progress                 # (group) Progress of the training process (for internally usage only).
    │       ├── ab_initial_finished  # (attr, (1,), bool) If the ab-initio calculation has been finished.
    │       ├── propagation_finished # (attr, (1,), bool) If the ML potential-based propagations has been finished.
    │       └── training_finished    # (attr, (unlimited,), List[bool]) If the potentials has been trained.
    ├── 1
    │   ├── ...
    │   ...
    ...
```
