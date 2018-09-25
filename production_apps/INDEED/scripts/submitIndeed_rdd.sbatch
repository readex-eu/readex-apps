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
echo SCOREP_TOTAL_MEMORY $SCOREP_TOTAL_MEMORY
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
export SCOREP_PROFILING_FORMAT=cube_tuple
export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_L3_TCM

###############
## Prepare directory structure
###############
cp -r ../testdata/$MY_JOB_NAME.orig $SCRATCH/$MY_JOB_ID
mkdir ../testdata/$MY_JOB_DIR
cd `dirname $MY_EXE`
cd ..
cp $SCOREP_FILTERING_FILE $SCRATCH/$MY_JOB_ID
cd $SCRATCH/$MY_JOB_ID

###############
## Run Indeed
###############
export MKL_DYNAMIC=FALSE
export OMP_DYNAMIC=FALSE
export MKL_CBWR=AVX2

echo starting Indeed
$MY_EXE -c $OMP_NUM_THREADS -j $MY_JOB_NAME
echo finished Indeed
echo

###############
## Run RDD
###############

echo starting rdd
readex-dyn-detect -t $RDD_GRANULARITY_THRESHOLD -p $RDD_PHASE_REGION_NAME\
                  -c $RDD_COMPUTE_INTENSITY_VARIATION_THRESHOLD -v $RDD_EXEC_TIME_VARIATION_THRESHOLD\
                  -w $RDD_REGION_EXEC_TIME_VS_PHASE_EXEC_TIME -r $RDD_CONFIG_FILE\
                  -f $RDD_RADAR_REPORT_FILE scorep-*/profile.cubex
echo finished rdd with exit code $?
echo
if [ -e $RDD_RADAR_REPORT_FILE ]; then
  pdflatex $RDD_RADAR_REPORT_FILE
else
  echo RADAR report file $RDD_RADAR_REPORT_FILE not found.
fi

###############
## Clean up
###############

cd $CURR_DIR
cp -r $SCRATCH/$MY_JOB_ID ../testdata/$MY_JOB_DIR/
rm -rf $SCRATCH/$MY_JOB_ID
mv submitIndeed_"$RUNTYPE"_"$SLURM_JOB_ID".sbatch ../testdata/$MY_JOB_DIR
mv slurm-$SLURM_JOB_ID.out ../testdata/$MY_JOB_DIR


