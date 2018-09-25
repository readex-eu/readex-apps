#!/bin/sh

source readex_env/set_env_saf.source

export SCOREP_PREP="scorep --mpp=mpi --thread=none"
export DSWITCH=""

make clean
rm -rf miniMD_openmpi miniMD_openmpi_saf
make openmpi
mv miniMD_openmpi miniMD_openmpi_saf

