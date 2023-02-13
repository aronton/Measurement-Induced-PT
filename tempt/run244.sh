d
#$ -pe smp 16

export OMP_NUM_THREADS=16

# source ${HOME}/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
 source ${HOME}/intel/mkl/bin/mklvars.sh intel64

d1=$(date +"%s")
./MIPT 24 24 2 2 0.25 0.25 0.05
d2=$(date +"%s")
d=$[d2-d1]
echo $d > time4.txt

