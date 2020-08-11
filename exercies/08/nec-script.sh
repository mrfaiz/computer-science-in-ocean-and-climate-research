#!/bin/bash
#PBS -T intmpi
#PBS -b 1
#PBS -l cpunum_job=16
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=10gb
#PBS -N test-N-DOP
#PBS -o test-N-DOP.out
#PBS -j o
#PBS -q clexpress
cd $PBS_O_WORKDIR
module load intel17.0.4 intelmpi17.0.4 petsc3.7.6intel
mpirun $NQSII_MPIOPTS -np 16 /sfs/fs5/home-sh/sunip565/metos3d-simpack-N-DOP.exe /sfs/fs5/home-sh/sunip565/model/N-DOP/option/test.N-DOP.option.txt
qstat -f ${PBS_JOBID/0:}
