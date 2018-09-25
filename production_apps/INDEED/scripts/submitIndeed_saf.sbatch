#!/bin/bash

# Submit this script with: sbatch -J jobname --time=timelimit --mail-user=mailaddress --cpus-per-task=ncpus --mem-per-cpu=memory thefilename

#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=haswell
#SBATCH -A p_readex
#SBATCH --mail-type=BEGIN,FAIL,END

############
## Init environment etc
############
. ../readex_env/set_env_$RUNTYPE.source

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. 
# Does nothing if the program doesn't use OpenMP.

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OUTFILE=""
CURR_DIR=`pwd`
MY_EXE="$CURR_DIR"/../bin/Indeed_"$RUNTYPE"64.exe
MY_JOB_NAME=`echo $SLURM_JOB_NAME | cut -f1 -d-`
MY_JOB_NAME=Benchmark_$MY_JOB_NAME
MY_JOB_ID="$MY_JOB_NAME"_$SLURM_JOB_ID
MY_JOB_DIR=Job_$SLURM_JOB_ID

echo RUNTYPE $RUNTYPE
echo OMP_NUM_THREADS $OMP_NUM_THREADS
echo SLURM_JOB_ID $SLURM_JOB_ID
echo DIRECTORY $CURR_DIR
echo
echo PROJECT $PROJECT
echo
echo MY_EXE $MY_EXE
echo MY_JOB_NAME $MY_JOB_NAME
echo MY_JOB_ID $MY_JOB_ID
echo
module list
echo


export SCOREP_FILTERING_FILE=scorep.filt

###############
## Prepare directory structure
###############
cp -r ../testdata/$MY_JOB_NAME.orig $SCRATCH/$MY_JOB_ID
mkdir ../testdata/$MY_JOB_DIR
cd $SCRATCH/$MY_JOB_ID

rm -rf scorep-*
rm -f old_scorep.filt
echo "" > scorep.filt
cp $CURR_DIR/do_scorep_autofilter_single.sh .

###############
## Run Indeed
###############
export MKL_DYNAMIC=FALSE
export OMP_DYNAMIC=FALSE
export MKL_CBWR=AVX2

result=1
iteration=0
while [ $result -ne 0 ]; do
  ((iteration++))
  echo "result = "$result "; running Indeed."
  $MY_EXE -c $OMP_NUM_THREADS -j $MY_JOB_NAME
  echo "Indeed run completed."

  sh do_scorep_autofilter_single.sh 0.001
  result=$?
  echo "scorep_autofilter_single done ($result)."
done

echo "end after $iteration runs of Indeed."

###############
## Clean up
###############
cp $SCOREP_FILTERING_FILE `dirname $MY_EXE`
cd `dirname $MY_EXE`
mv $SCOREP_FILTERING_FILE ..

cd $CURR_DIR
cp -r $SCRATCH/$MY_JOB_ID ../testdata/$MY_JOB_DIR/
rm -rf $SCRATCH/$MY_JOB_ID
mv submitIndeed_"$RUNTYPE"_"$SLURM_JOB_ID".sbatch ../testdata/$MY_JOB_DIR
mv slurm-$SLURM_JOB_ID.out ../testdata/$MY_JOB_DIR


