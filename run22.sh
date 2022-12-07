#!/bin/bash

#$ -cwd
#$ -pe smp 10

export OMP_NUM_THREADS=10

# source ${HOME}/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
 source ${HOME}/intel/mkl/bin/mklvars.sh intel64


./MIPT 22 22 2 1000 0.1 0.8 0.05

