#$ -S /bin/bash
#PBS -q workq
#PBS -N RECan
#PBS -J 1-300

BASE=$PBS_O_WORKDIR
cd $BASE

date

source activate PyVenv


DIR_LIST="${BASE}/directory.list"
DIR=$(sed "${PBS_ARRAY_INDEX}p" -n $DIR_LIST)

  echo $DIR

  cd $DIR

  python /home/jmcbride/RecProt/Src/python/post-processing.py


date
