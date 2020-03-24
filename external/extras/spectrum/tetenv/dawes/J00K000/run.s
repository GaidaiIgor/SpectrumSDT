#PBS -q debug
#PBS -l mppwidth=1
#PBS -l walltime=00:30:00
#PBS -N tetenv.J00
#PBS -j eo

cd $PBS_O_WORKDIR

aprun -n 1  ../tetenv

