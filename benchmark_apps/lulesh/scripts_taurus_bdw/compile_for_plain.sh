#!/bin/sh

cd ..
source readex_env/set_env_plain.source
#module use /projects/p_readex/modules
#module load readex/ci_readex_intelmpi2017.2.174_intel2017.2.174_scorep_patched
#source ./readex_env/set_env_plain.source
#module load readex/scs5_ci_readex_gcc7.3.0
export SCOREP_PREP=""
export DSWITCH=""

make clean
rm -rf lulesh2.0 lulesh2.0_plain
make
mv lulesh2.0 lulesh2.0_plain
