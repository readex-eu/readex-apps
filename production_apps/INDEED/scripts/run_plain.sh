#!/bin/bash

# usage: ./run_plain.sh [exampleID] 
#        where exampleID may be "B1", "SBC" or "H1N"
#        default value is B1


if [ $# -eq 0 ]; then
    ARG="B1"
else
    ARG=`echo $1 | tr [a-z] [A-Z]`;
fi

if [ $ARG == "SBC" ]; then
    JOBNAME="SBC-OMP-I16";
    JOBTIME="01:00:00";
    JOBCPUS=24;
    JOBMEM=7800;
elif [ $ARG == "B1" ]; then
    JOBNAME="B1-OMP-I16";
    JOBTIME="00:45:00";
    JOBCPUS=24;
    JOBMEM=4000;
elif [ $ARG == "H1N" ]; then
    JOBNAME="H1N-OMP-I16";
    JOBTIME="00:10:00"
    JOBCPUS=24;
    JOBMEM=4000;
else
    echo $0: Illegal exampleID $1
    echo Legal values are "B1", "SBC" and "H1N"
    exit 1
fi

export RUNTYPE="plain"

CURRDIR=`pwd`
cd `dirname $0`

if [ `echo ${MYMAILADDRESS:=none}` != "none" ]; then
    ADDARG="--mail-user=$MYMAILADDRESS"
else
    ADDARG=""
fi

MY_SLURM_ID=`sbatch -J $JOBNAME --time=$JOBTIME $ADDARG --cpus-per-task=$JOBCPUS --mem-per-cpu=$JOBMEM submitIndeed_"$RUNTYPE".sbatch | cut -d" " -f4`
echo '##' Command: sbatch -J $JOBNAME --time=$JOBTIME $ADDARG --cpus-per-task=$JOBCPUS --mem-per-cpu=$JOBMEM submitIndeed_"$RUNTYPE".sbatch > submitIndeed_"$RUNTYPE"_$MY_SLURM_ID.sbatch
echo '##' was executed with following file: >> submitIndeed_"$RUNTYPE"_$MY_SLURM_ID.sbatch
echo >> submitIndeed_"$RUNTYPE"_$MY_SLURM_ID.sbatch
cat submitIndeed_"$RUNTYPE".sbatch >> submitIndeed_"$RUNTYPE"_$MY_SLURM_ID.sbatch

echo Job $MY_SLURM_ID submitted.

cd $CURRDIR

