#!/bin/bash

#source env/set_env_npb_ptf.sh
source readex_env/set_env_rdd.source

# Check environment
module list

#Compile NPB
cp config/make.def.rdd config/make.def

make bt-mz CLASS=C NPROCS=2

mv bin/bt-mz.C.2 bin/bt-mz.C.2_rdd
