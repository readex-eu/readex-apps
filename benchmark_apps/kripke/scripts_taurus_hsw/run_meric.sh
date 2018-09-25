#!/bin/sh

# 1) Modules and variables
cd ../build
. ../readex_env/set_env_meric.source
. ../environment.sh

export MERIC_FREQUENCY=25
export MERIC_UNCORE_FREQUENCY=30
export MERIC_NUM_THREADS=0

# 2) Run
#srun -n 24 ./kripke --procs 2,2,6 --niter 20 --nest GZD --zones 4,4,4 --groups 200
#srun -n 24 ./kripke --procs 2,2,6 --niter 10 --nest GZD --zones 32,32,32 --legendre 8 --dset 32
srun -n 24 ./kripke $KRIPKE_COMMAND
