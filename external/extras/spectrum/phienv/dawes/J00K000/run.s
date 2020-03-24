#PBS -q regular
#PBS -l mppwidth=1
#PBS -l walltime=01:00:00
#PBS -N phienv.J00
#PBS -j eo

cd $PBS_O_WORKDIR

aprun -n 1  ../phienv

