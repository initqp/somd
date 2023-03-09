# SOMD basics

SOMD is an ab-initio molecular dynamics (AIMD) package designed for the
[SIESTA](https://departments.icmab.es/leem/siesta/) DFT code. The SOMD code
provides some useful subroutines (potential calculator interfaces,
integrators, barostats, etc.) to perform standard Born-Oppenheimer molecular
dynamics (BOMD) simulations. Thus, SOMD itself should be regarded as a
library, just like the [OpenMM](https://openmm.org/) package. And you should
write proper python scripts to assemble the provided subroutines into
meaningful simulation protocols.

## Module structures
The SOMD code contains four main submodules: `constants`, `core`, `potentials`
and `apps`.
- The `constants` module contains physical constants and default values of
  SOMD.
- The `core` module contains core components of SOMD, e.g., the internal
  representation of an atomic system, the integrators, etc.
- The `potentials` module contains interfaces to various potential calculators.
- The `apps` module contains some useful applications, e.g., the trajectory
  writers, system data loggers, pre-defined simulation protocols, etc. (Note
  that the barostats are regarded as applications as well, because the
  first-order barostats implemented in SOMD behave like other
  post-integration-step applications.)

```
somd
├── constants  <-- Physical constants and default values.
├── core       <-- Core components.
│   ├── systems         <-- Representations of the simulated systems.
│   ├── groups          <-- Atom groups.
│   ├── integrators     <-- Integrators to propagate the systems.
│   └── potential_base  <-- Base class of all potential calculator interfaces.
├── potentials <-- Potential calculator interfaces.
│   ├── SIESTA          <-- SIESTA wrapper.
│   ├── DFTD3           <-- DFTD3 wrapper.
│   ├── DFTD4           <-- DFTD4 wrapper.
│   ├── NEP             <-- NEP_CPU wrapper.
│   ├── PLUMED          <-- PLUMED wrapper.
│   ...
├── apps       <-- Useful applications based on the SOMD core.
│   ├── barostats       <-- Barostats.
│   ├── trajectories    <-- Trajectory writers.
│   ├── loggers         <-- System data loggers.
│   ├── simulations     <-- Wrapper class of the simulation processes.
│   ├── active_learning <-- Wrapper class of the active learning processes.
│   ...
```
