A parallel Fortran program for calculation of ro-vibrational energy levels and lifetimes of ABA-molecules. Note that the current version is limited to wave functions that do not extend into linear and equilateral triangle configurations (see manual for details).

# Building

0. Prerequisites
    1. Make sure the following packages are installed: `build-essential`, `python3-dev`, `cmake` (3.5+), `gfortran` (9.3.0+), `mpich`, `libblas-dev`, `liblapack-dev`.  
    With the appropriate modifications of the build configuration files, the code should also work with `ifort`. In our experience, `ifort` generates much faster code than `gfortran`, therefore we recommend to change the default settings to build with `ifort` on system where it is available. 
    2. Make sure you machine has at least 2GB of RAM (for compilation).

1. Clone the repo
`git clone https://github.com/IgorGayday/SpectrumSDT`  
This example assumes the repo is cloned into `~/SpectrumSDT`

2. Build the libraries

    1. fdict
    
    ```
    cd ~/SpectrumSDT/libs/fdict
    ```
    Edit `setup.make` to change compiler options (optional). Then run:
    ```
    make
    ```

    2. PETSc
    ```
    cd ~/SpectrumSDT/libs/petsc
    export PETSC_DIR=$PWD
    export PETSC_ARCH=debug_64
    ./configure --with-scalar-type=complex --with-64-bit-indices
    make all
    ```

    3. SLEPc
    ```
    cd ~/SpectrumSDT/libs/slepc
    export SLEPC_DIR=$PWD
    ./configure
    make all
    ```

3. Build the main program  
`cd ~/SpectrumSDT`

    1. Edit `CMakeLists.txt` to specify custom compiler options (optional)  

    2. Build the executable  
    Note: do not move the executable to another folder since some paths are resolved relative to its location.  
    ```
    mkdir build && cd build
    cmake ..
    make spectrumsdt
    ```

# A basic example of running (ozone)

1. Generate the grids
```
mkdir -p ~/SpectrumSDT_runs/o3/ && cd ~/SpectrumSDT_runs/o3/
cp ~/SpectrumSDT/config_examples/o3/1.simple/grids.config spectrumsdt.config
~/SpectrumSDT/build/spectrumsdt
```
Before proceeding to the next stage, user should provide file `pes.out` with the values of potential at all combination of points in the generated grid files in atomic units of energy (Hartree).

2. Calculate the values of PES. Here we will use an example program that reads grid files and uses the PES of ozone calculated by Dawes et al. to generate `pes.out` file. First, compile the program:
```
cd ~/SpectrumSDT/PES_examples/ozone/
mkdir build && cd build
../compile.sh
```
Now run:
```
cd ~/SpectrumSDT_runs/o3/
mpiexec -n <n_procs> ~/SpectrumSDT/PES_examples/ozone/build/ozone_pes 686
```
Replace `<n_procs>` with however many MPI tasks you want to use. After this, `pes.out` file with the values of PES at all grid points will be written.

3. Setup SpectrumSDT directory structure
```
mkdir J_0 && cd J_0
cp ~/SpectrumSDT/config_examples/o3/1.simple/spectrumsdt.config .
```
Edit spectrumsdt.config and replace `username` in the paths. Then execute:  
`~/SpectrumSDT/scripts/init_spectrum_folders.py -K 0`

4. Calculate basis
```
cd K_0/symmetry_0/basis
mpiexec -n <n_procs> spectrumsdt
```

5. Calculate basis cross terms (overlaps)
```
cd ../overlaps
mpiexec -n <n_procs> spectrumsdt
```

6. Calculate eigenpairs
```
cd ../eigensolve
mpiexec -n <n_procs> spectrumsdt
```
Lowest 50 rovibrational energy levels of ozone-686 J=0 will be printed into `states.fwc` file.

7. Calculate wave function properties (optional)
```
cd ../properties
mpiexec -n <n_procs> spectrumsdt
```
The requested properties will be printed into `states.fwc` file.

# References

1. Gayday, I.; Grushnikova, E.; Babikov, D. Influence of the Coriolis Effect on the Properties of Scattering Resonances in Symmetric and Asymmetric Isotopomers of Ozone. Phys. Chem. Chem. Phys. 2020, 22, 27560–27571.  

2. Gayday, I.; Teplukhin, A.; Kendrick, B. K.; Babikov, D. The Role of Rotation–Vibration Coupling in Symmetric and Asymmetric Isotopomers of Ozone. J. Chem. Phys. 2020, 152, 144104.  

3. Gayday, I.; Teplukhin, A.; Kendrick, B. K.; Babikov, D. Theoretical Treatment of the Coriolis Effect Using Hyperspherical Coordinates, with Application to the Ro-Vibrational Spectrum of Ozone. J. Phys. Chem. A 2020, 124, 2808–2819.  

4. Teplukhin, A.; Babikov, D. Efficient Method for Calculations of Ro-Vibrational States in Triatomic Molecules near Dissociation Threshold: Application to Ozone. J. Chem. Phys. 2016, 145, 114106.  
