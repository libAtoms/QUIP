# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

QUIP (QUantum mechanics and Interatomic Potentials) is a collection of Fortran software tools for molecular dynamics simulations. It implements various interatomic potentials, tight binding quantum mechanics, and can interface with external packages like LAMMPS, CP2K, and ASE. The codebase includes:

- **Core library (libAtoms)**: Fundamental atomistic data structures and algorithms
- **Interatomic Potentials**: Various potential models (EAM, Tersoff, Stillinger-Weber, Lennard-Jones, etc.)
- **GAP (Gaussian Approximation Potentials)**: Machine learning interatomic potentials (separate license - Academic Software License)
- **Python interface (quippy)**: Python bindings via f90wrap for integration with ASE
- **Tight-binding modules**: DFTB, NRL-TB, and other TB implementations
- **Analysis programs**: 50+ Fortran utilities in `src/Programs/`

## Build System

QUIP uses the [Meson build system](https://mesonbuild.com/) (version 1.1+) with the Ninja build tool.

### Initial Setup

1. Ensure you have the required tools:
```bash
pip install meson ninja
```

2. Update submodules to meson-compatible versions:
```bash
cd src/fox && git checkout master && git pull
cd ../GAP && git checkout main && git pull
cd ../..
```

3. Configure the build:
```bash
meson setup builddir
```

This creates a `builddir` directory where all build artifacts will be placed. Meson automatically detects your compiler, BLAS/LAPACK libraries, and system configuration.

### Build Options

You can customize the build using `-D` flags during setup:
```bash
meson setup builddir -Dgap=true -Dmpi=false
```

Available options (defined in `meson.options`):
- `gap` (default: `true`): Enable GAP (Gaussian Approximation Potentials) support
- `mpi` (default: `false`): Enable MPI parallelization (requires ScaLAPACK)

To reconfigure an existing build:
```bash
meson configure builddir -Dmpi=true
```

### Building

```bash
# Build everything (libraries, programs)
meson compile -C builddir

# Or use ninja directly
ninja -C builddir
```

Build artifacts are organized in `builddir/`:
- Executables: `builddir/src/Programs/quip`, `builddir/src/Programs/gap_fit`, etc.
- Shared libraries: `builddir/src/libAtoms/liblibAtoms.so`, `builddir/src/Potentials/libPotentials.so`, etc.

Programs can be run directly:
```bash
./builddir/src/Programs/quip --help
./builddir/src/Programs/gap_fit --help
```

### Python Interface (quippy)

After building the main QUIP libraries:

```bash
cd quippy
meson setup builddir
meson compile -C builddir
meson install -C builddir  # Install to Python environment
```

### Testing

```bash
# Run all tests using meson
meson test -C builddir

# Or run tests manually from tests directory
cd tests
python3 run_all.py -v
```

Tests are Python-based and located in `tests/`. They use the quippy interface.

### Cleaning

```bash
# Remove build directory and start fresh
rm -rf builddir
meson setup builddir
```

## Architecture

### Module Structure

The codebase is organized into modular layers:

1. **libAtoms** (`src/libAtoms/`): Core library
   - `Atoms.F90`: Main atomic structure data type
   - `Atoms_types.F90`: Type definitions and structure
   - `Connection.F90`: Neighbor lists and connectivity
   - `Dictionary.F90`: Key-value parameter handling
   - `DynamicalSystem.F90`: MD integrators and thermostats
   - `CInOutput.F90`: Extended XYZ I/O format
   - `Table.F90`, `Matrix.F90`: Data structures
   - `MPI_context.F90`: MPI parallelization support

2. **Potentials** (`src/Potentials/`): Interatomic potential implementations
   - `IP.F90`: Base interatomic potential interface
   - `IPModel_*.F90`: Individual potential models (LJ, SW, Tersoff, EAM, etc.)
   - `FilePot.F90`: External potential interface
   - `AdjustablePotential.F90`: Parameter optimization support

3. **GAP** (`src/GAP/`): Gaussian Approximation Potentials (submodule)
   - Separate license (Academic Software License)
   - Machine learning potential framework
   - SOAP descriptors and sparse GP regression
   - Documentation at https://libatoms.github.io/GAP

4. **Programs** (`src/Programs/`): Standalone executables (50+ programs)
   - `quip.F90`: Main program for energy/force calculations
   - `md.F90`: Molecular dynamics
   - `gap_fit.F90`: GAP training (if HAVE_GAP=1)
   - Various analysis and structure manipulation tools

5. **Utils** (`src/Utils/`): Utility modules and functions

6. **FoX** (`src/fox/`): XML parsing library (submodule)

7. **f90wrap** (`src/f90wrap/`): Python wrapping tool (submodule)

### Key Concepts

**Extended XYZ Format**: QUIP's primary I/O format. Structure:
```
<number_of_atoms>
Lattice="..." Properties=species:S:1:pos:R:3:... [other key=value pairs]
<atom_data>
```
Properties field specifies column format as `name:type:columns` where type is I (integer), R (real), or S (string).

**Potential Initialization**: Potentials are initialized via string arguments:
```fortran
call initialise(pot, 'IP LJ', param_filename='ip.parms.LJ.xml')
```
Use `init_args='--help'` to see available potential types. Recursive help: `init_args='IP --help'` lists interatomic potential types.

**Build Configuration**: Meson automatically configures:
- Compiler detection and flags (gfortran, ifort, etc.)
- Math libraries (OpenBLAS, MKL, reference BLAS/LAPACK via pkg-config)
- Feature flags via build options (`-Dgap=true`, `-Dmpi=false`)
- Python detection for quippy

Configuration can be inspected or modified:
```bash
meson configure builddir
meson configure builddir -Dmpi=true
```

### Parallelization

- **OpenMP**: Thread-level parallelization (automatically detected if compiler supports it)
- **MPI**: Domain decomposition for large systems (enable with `-Dmpi=true`)
- **ScaLAPACK**: Required for MPI-parallel `gap_fit` (automatically linked when MPI is enabled)

## Common Development Patterns

### Adding a New Potential

1. Create `src/Potentials/IPModel_NewPot.F90` following existing model patterns
2. Add the new file to the `Potentials_F90_sources` list in `src/Potentials/meson.build`
3. Register in `src/Potentials/IP.F90` in the initialise and calc routines
4. Add parameter XML schema if needed
5. Rebuild with `meson compile -C builddir`

### Working with the Python Interface

After building quippy:
```python
from quippy.potential import Potential
from ase.io import read

atoms = read('structure.xyz')
pot = Potential('IP LJ', param_filename='ip.parms.LJ.xml')
energy = pot.calc(atoms, energy=True)
```

The quippy interface wraps Fortran types/routines, providing ASE calculator interface.

### Modifying Core libAtoms

Changes to `src/libAtoms/` affect all downstream modules. After modifications:
```bash
meson compile -C builddir  # Meson automatically handles dependencies
```

Meson's dependency tracking is generally more reliable than Make, so full rebuilds are rarely needed.

### Working with GAP

GAP is a git submodule with separate license terms. To update:
```bash
cd src/GAP
git checkout main  # or specific commit/branch
git pull
cd ../..
git add src/GAP
git commit -m "Update GAP version"
```

Enable GAP in build: Use `-Dgap=true` (enabled by default). For MPI-parallel gap_fit, also enable MPI: `-Dmpi=true` (requires ScaLAPACK).

## Important Notes

- **Do not assume minimum image convention**: Potentials can have cutoffs larger than the unit cell
- **Fortran preprocessing**: `.F90` files are preprocessed (can use `#ifdef`), `.f90` are not
- **Feature flags**: Many features are optional (TB, GAP, MPI, etc.) - controlled via meson options
- **Virtual environments**: Activate virtualenv before building quippy to install there
- **GAP license**: GAP has a non-commercial academic license, distinct from QUIP's GPL
- **f90wrap**: Install development version with `pip install git+https://github.com/jameskermode/f90wrap.git@master` before building quippy
- **Submodule versions**: Ensure fox is on `master` and GAP is on `main` for meson support

## External Interfaces

QUIP can be used as:
- **Standalone**: Direct execution of compiled programs from `builddir/src/Programs/`
- **Library**: Link against shared libraries in `builddir/src/*/lib*.so` (e.g., for LAMMPS pair_style quip)
- **Python**: Via quippy with ASE integration
- **Plugins**: CP2K, LAMMPS, ASE can call QUIP potentials

LAMMPS integration: Build QUIP libraries with meson, then follow LAMMPS pair_quip documentation (requires LAMMPS 11 Aug 2017+).

## File Locations

- Source code: `src/<module>/`
- Meson build files: `meson.build`, `meson.options`, `src/*/meson.build`
- Build outputs: `builddir/` (configurable name)
- Parameter files: `share/Parameters/`
- Example structures: `share/Structures/`
- Tests: `tests/`
- Documentation: https://libatoms.github.io/QUIP/
- GAP docs: https://libatoms.github.io/GAP/
