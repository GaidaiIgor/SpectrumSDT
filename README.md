A parallel Fortran program for calculation of ro-vibrational energy levels and lifetimes of 3-atomic systems in hyperspherical (APH) coordiantes.

# Building

0. Prerequisites
    1. Make sure the following packages are installed: `build-essential`, `python3-dev`, `cmake`, `gfortran`, `mpich`, `libblas-dev`, `liblapack-dev`, `libpcre3-dev`.  
    Tested with gfortran 9.3.0 and mpich 3.3.2.
    2. Make sure you machine has at least 2GB of RAM (for compilation)

1. Clone the repo
`git clone https://github.com/IgorGayday/SpectrumSDT`  
This example assumes the repo is cloned into `~/SpectrumSDT`

2. Build the libraries

    1. fdict
    ```
    cd ~/SpectrumSDT/libs/fdict
    printf 'FC=gfortran\nFFLAGS = -g\n' > setup.make
    make
    ```

    2. FTL
    ```
    cd ~/SpectrumSDT/libs/ftl
    make
    PREFIX=build make install
    ```

    3. PETSc
    ```
    cd ~/SpectrumSDT/libs/petsc
    export PETSC_DIR=$PWD
    export PETSC_ARCH=debug_64
    ./configure --with-scalar-type=complex --with-64-bit-indices
    make all
    ```

    4. SLEPc
    ```
    cd ~/SpectrumSDT/libs/slepc
    export SLEPC_DIR=$PWD
    ./configure
    make all
    ```

3. Build the main programs  
`cd ~/SpectrumSDT`

    1. Specify custom compiler options (optional)  
    Default compiler options are specified in `compiler_options_default.cmake`.  
    You can create a file named `compiler_options.cmake` to specify your custom compile options instead of the default ones.  
    A few presets for different compilers can be found in the `compiler_options_alternative` folder.

    2. Build the executables  
    Note: do not move the executables to another folder since some paths are resolved relative to their location.  
    ```
    mkdir build && cd build
    cmake ..
    make optgrid
    make pesprint
    make spectrumsdt
    ```

# A basic example of running

1. Generate the grids
```
mkdir ~/SpectrumSDT_runs && cd ~/SpectrumSDT_runs
cp ~/SpectrumSDT/config_examples/optgrid_smaller.config optgrid.config
~/SpectrumSDT/build/optgrid
```
If executed successfully, three grid files should appear in the current folder.

2. Calculate the values of PES
```
cp ~/SpectrumSDT/config_examples/pesprint.config spectrumsdt.config
mpiexec -n <n_procs> ~/SpectrumSDT/build/pesprint
```
Replace `<n_procs>` with however many MPI tasks you want to use.  
If executed successfully, potvib.dat file should appear in the current folder.

3. Setup SpectrumSDT directory structure
```
mkdir J_0 && cd J_0
cp ~/SpectrumSDT/config_examples/spectrumsdt_smaller.config spectrumsdt.config
```
Edit spectrumsdt.config and replace username in the paths.  
`~/SpectrumSDT/scripts/init_spectrum_folders.py -K 0`

4. Calculate basis
```
cd K_0/even/basis
mpiexec -n <n_procs> spectrumsdt
```
Here `<n_procs>` has to be equal to the number of points in `~/SpectrumSDT_runs/grid1.dat` (16 in this example).  
In `fix_basis_jk` mode (enabled in this example), basis of the other symmetry has to be computed as well.  
```
cd ../../odd/basis
mpiexec -n 16 spectrumsdt
```

5. Calculate basis cross terms (overlaps)
```
cd ../../even/overlaps
mpiexec -n <n_procs> spectrumsdt
```

6. Calculate eigenpairs
```
cd ../diagonalization
mpiexec -n <n_procs> spectrumsdt
```
Solutions will be printed into `3dsdt/spec.out` file.
