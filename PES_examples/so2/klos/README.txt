These are Fortran routines for the SO2 ground (X) and C(1B2)(2A') electronic PESs that were used in following published research: 
Jacek Klos, Millard H. Alexander,Praveen Kumar, Bill Poirier,Bin Jiang, and Hua Guo The Journal of Chemical Physics 144, 174301 (2016)

Files in this Zenodo repository: 
Data files:
C_so2_bin.dat
X_ccf12_klos.dat
X_mrcisdf12.dat
Fortran files:
so2-C-mrcif12.f
so2-X-ccsdtf12.f
so2-X-mrcif12.f
Test files:
test_output_C_MRCIF12.txt
test_output_X_CCF12.txt
test_output_X_MRCIF12.txt

README.txt

The test output files provide values of the PES at given geometry to check compilation. Tests were done with intel ifort compilation: ifort so2-X-ccsdtf12.f -o so2-X-ccsdtf12.exe for example. 

