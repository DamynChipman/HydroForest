# HydroForest

A finite element code for studying hydrodynamics on an adaptive `p4est` mesh.

## Dependencies

This code currently requires the following software libraries to be installed:

- `conda` : To create a Python environment
- `mpi` : Working `mpicc` and `mpic++` compilers
- `p4est` : For adaptive meshing
- `petsc` : For parallel linear algebra (not yet...)

The user will need to create a Python environment via `conda` that has `matplotlib` and `numpy` installed. This is in order to use the `matplotlibcpp` package in C++ for visualizing 1D data. 2D data will be visualized with VTK in the future.

### Create the `conda` Environment

To create a `conda` environment, run the following with `conda`:

```bash
$ conda env create -f environment.yml
$ conda activate HydroForest
```

Your bash prompt should now have `(HydroForest)` as part of the prompt.

## Installation

NOTE: Currently, the user needs to have `p4est` and `petsc` already installed. Currently working on the build system to do it automatically if not installed, but not yet...

### Configuration

First, you'll need to configure the `config.sh` script with your specific paths. Open up `config.sh` and edit the following paths:

- `HOME` : Path to your home directory
- `PYTHON_ENV_PATH` : Path to `conda` `HydroForest` environment directory
- `PYTHON_VERSION` : Version of Python in `conda` `HydroForest` environment (`python3.9`)
- `P4EST_PATH` : Path to `p4est` installation
- `PETSC_PATH` : Path to PETSc installation
- `HYDRO_FOREST` : Absolute path to source code for HydroForest

Once those are configured, run the configuration script in the root directory of HydroForest to create a build directory of the current Git branch:

```bash
$ ./config.sh
```

### Building

Once configured and CMake runs, `cd` into the build directory and run `make`:

```bash
cd build-${CURRENT_GIT_BRANCH}
make
```

If you run into issues, please reach out or submit an issue on GitHub.