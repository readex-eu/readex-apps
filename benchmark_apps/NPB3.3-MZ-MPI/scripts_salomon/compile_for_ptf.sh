#!/bin/bash
cd ..

source readex_env/set_env_ptf_rapl.source

if [ "$READEX_INTEL" == "1" ]; then
  export FILTER_GNU=""
  export FILTER_INTEL="-tcollect-filter=scorep_icc.filt"
else
  export FILTER_GNU="--instrument-filter=scorep.filt"
  export FILTER_INTEL=""
fi

# Check environment
module list

#Compile NPB
cp config/make.def.ptf config/make.def

CLASS=C
NPROCS=2
make bt-mz CLASS=${CLASS} NPROCS=${NPROCS}

mv bin/bt-mz.${CLASS}.${NPROCS} bin/bt-mz.${CLASS}.${NPROCS}_ptf
