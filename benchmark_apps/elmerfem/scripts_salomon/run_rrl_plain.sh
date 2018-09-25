#!/bin/bash

#PBS -A OPEN-12-63
#PBS -q qprod
#PBS -N ELMER_compare
#PBS -l select=4:ncpus=24:mpiprocs=24:accelerator=false
#PBS -l x86_adapt=true
#PBS -l walltime=01:00:00
#PBS -m be
###############################################################################

if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
  LOC=$PBS_O_WORKDIR
else
  LOC=$(pwd)
fi

cd $LOC

module purge

module load mkl
module load CMake/3.9.1

. readex_env/set_env_rrl.source
. readex_env/set_env_ptf_rapl.source
. ../paths.source


export PROBLEM_PATH="${ELMER_ROOT}/build/fem/tests/WinkelPoissonPartitionUniform/"
cd ${PROBLEM_PATH}

cp ${ELMER_ROOT}/tuning_model.json .
cp ${ELMER_ROOT}/meric.opts .

TMPFILE=elmer_rrl_plain.tmp
rm -f $TMPFILE

# SCRIPT-SPECIFIC PART #########################################################
Evaluate (){
	TIME=$(grep -oP "(?<=Runtime \[s\] = )\d+\.\d+"  <<< $OUTPUT | head -1)
	SUMTIME=$( awk "BEGIN {print $SUMTIME + $TIME}" )
	CPUs=$(grep -oP "(?<=Overall CPUs energy consumption \[J\] = )\d+"  <<< $OUTPUT )
	NODEs=$(grep -oP "(?<=Overall nodes energy consumption \[J\] = )\d+"  <<< $OUTPUT )
	
	for VAL in $CPUs
	do
		SUMCPUs=$( awk "BEGIN {print $SUMCPUs + $VAL}" )
	done
	for VAL in $NODEs
	do
		SUMNODEs=$( awk "BEGIN {print $SUMNODEs + $VAL}" )
	done
}

# RUN PLAIN APP ################################################################
SUMTIME=0.0
SUMCPUs=0.0
SUMNODEs=0.0
i=0

REPEAT_COUNT=5

REPEAT_COUNT=2 #5

cd ${ELMER_ROOT}/scripts_salomon
./compile_for_plain.sh

echo "Elmer was compiled for PLAIN."

cd ${PROBLEM_PATH}
rm -f case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif
if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

while [ $i -lt $REPEAT_COUNT ]; do
	$MERIC_ROOT/tools/staticMERICtool/multiNodeStaticMeasureStart.sh --rapl
	mpirun -n 96 ElmerSolver_mpi | tee -a $TMPFILE
	OUTPUT=$($MERIC_ROOT/tools/staticMERICtool/multiNodeStaticMeasureStop.sh --rapl)
	Evaluate
	i=$[ $i + 1 ]
done

echo "PLAIN time:" $( awk "BEGIN {print $SUMTIME / $REPEAT_COUNT}" )
echo "PLAIN CPUs energy:" $( awk "BEGIN {print $SUMCPUs / $REPEAT_COUNT}" )
echo "PLAIN NODEs energy:" $( awk "BEGIN {print $SUMNODEs / $REPEAT_COUNT}" )


# RUN TUNED APP ################################################################
SUMTIME=0.0
SUMCPUs=0.0
SUMNODEs=0.0
i=0

export SCOREP_ENABLE_PROFILING="false"
export SCOREP_ENABLE_TRACING="false"
export SCOREP_SUBSTRATE_PLUGINS="rrl"
export SCOREP_RRL_PLUGINS="cpu_freq_plugin,uncore_freq_plugin,OpenMPTP"
export SCOREP_RRL_TMM_PATH="tuning_model.json"
export SCOREP_MPI_ENABLE_GROUPS=ENV
export SCOREP_RRL_CHECK_IF_RESET="reset"

REPEAT_COUNT=2 #5

cd ${ELMER_ROOT}/scripts_salomon
./compile_for_ptf.sh

echo "Elmer was compiled for PLAIN."

cd ${PROBLEM_PATH}
rm -f case.sif
cp ${ELMER_ROOT}/case_cg.sif case.sif
if [ ! -d ${PROBLEM_PATH}/winkel/partitioning.96 ]; then
    ${ELMER_ROOT}/gen_problem.sh
fi

while [ $i -lt $REPEAT_COUNT ]; do
	$MERIC_ROOT/tools/staticMERICtool/multiNodeStaticMeasureStart.sh --rapl
	mpirun -n 96 ElmerSolver_mpi | tee -a $TMPFILE
	OUTPUT=$($MERIC_ROOT/tools/staticMERICtool/multiNodeStaticMeasureStop.sh --rapl)
	Evaluate
	i=$[ $i + 1 ]
done

echo "RRL time:" $( awk "BEGIN {print $SUMTIME / $REPEAT_COUNT}" )
echo "RRL CPUs energy:" $( awk "BEGIN {print $SUMCPUs / $REPEAT_COUNT}" )
echo "RRL NODEs energy:" $( awk "BEGIN {print $SUMNODEs / $REPEAT_COUNT}" )

rm -f $TMPFILE

