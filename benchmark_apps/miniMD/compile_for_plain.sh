#!/bin/sh

source readex_env/set_env_plain.source

export SCOREP_PREP=""
export DSWITCH=""

make clean
rm -rf miniMD_openmpi miniMD_openmpi_plain
make openmpi
mv miniMD_openmpi miniMD_openmpi_plain

