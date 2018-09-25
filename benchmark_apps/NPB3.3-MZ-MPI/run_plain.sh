#!/bin/bash

## To run the batch script, use command like 'sbatch --job-name=bt-mz.C.2_plain --output=bt-mz.C.2_plain.out run_plain.sh bt-mz.C.2_plain' to specify the output file
## To run the script interactively, simply comment all the SBATCH commnads below and switch the srun command to the commneted one below and run the command like ./run_plain.sh bt-mz.C.2_plain

#SBATCH --time=1:00:00   # walltime
#SBATCH --nodes=1  # number of processor cores (i.e. tasks)
##SBATCH --ntasks=
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --exclusive
#SBATCH --partition=haswell
#SBATCH --comment="cpufreqchown"
#SBATCH --mem-per-cpu=2500M   # memory per CPU core
#SBATCH -A p_readex
##SBATCH -J "my_job"   # job name
##SBATCH --output=my_job.out
##SBATCH --error=myjob.out
###############################################################################




source readex_env/set_env_plain.source
module list


app=bt-mz.C.2_plain
cd bin
srun $app

