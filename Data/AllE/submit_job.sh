#$ -S /bin/bash
#PBS -q workq
#PBS -N REC
#PBS -J 1-100

BASE=$PBS_O_WORKDIR
cd $BASE

date

source activate PyVenv

module load gcc/10.2.0

DIR_LIST="${BASE}/directory.list"
DIR=$(sed "${PBS_ARRAY_INDEX}p" -n $DIR_LIST)

echo $DIR

#f [ ! -f $DIR/output.dat ] ; then

  cd $DIR
  
  ./run_all.sh

#i

date
