#!/bin/bash

#source env/set_env_npb_plain.sh
source readex_env/set_env_plain.source

# Check environment
module list

#Compile NPB
cp config/make.def.plain config/make.def

make bt-mz CLASS=C NPROCS=2

mv bin/bt-mz.C.2 bin/bt-mz.C.2_plain
