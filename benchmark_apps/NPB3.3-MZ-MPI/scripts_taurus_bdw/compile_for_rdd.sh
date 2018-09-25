#!/bin/bash
cd ..
#source env/set_env_npb_ptf.sh
source readex_env/set_env_rdd.source

# Check environment
module list

#Compile NPB
cp config/make.def.rdd config/make.def
CLASS=C
NPROCS=24
make bt-mz CLASS=${CLASS} NPROCS=${NPROCS}

mv bin/bt-mz.${CLASS}.${NPROCS} bin/bt-mz.${CLASS}.${NPROCS}_rdd
