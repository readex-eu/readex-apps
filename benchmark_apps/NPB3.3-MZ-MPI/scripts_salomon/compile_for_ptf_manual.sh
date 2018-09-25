#!/bin/bash
cd ..
#source readex_env/set_env_ptf_hdeem.source
source ./readex_env/set_env_ptf_rapl.source
# Check environment
module list

#Compile NPB
cp config/make.def.ptf_manual config/make.def

CLASS=C
NPROCS=2
make bt-mz CLASS=${CLASS} NPROCS=${NPROCS}

mv bin/bt-mz.${CLASS}.${NPROCS} bin/bt-mz.${CLASS}.${NPROCS}_ptf
