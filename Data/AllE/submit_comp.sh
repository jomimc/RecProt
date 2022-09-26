#$ -S /bin/bash
#PBS -q workq
#PBS -l select=1:ncpus=28:mpiprocs=28
#PBS -N RECcom
# PBS -J 1-14

BASE=$PBS_O_WORKDIR
cd $BASE

date

source activate PyVenv

  python /home/jmcbride/RecProt/Src/python/post-processing.py --mode compile -t 0

date
