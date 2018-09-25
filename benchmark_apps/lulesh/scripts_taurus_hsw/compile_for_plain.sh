#!/bin/sh

cd ..
source readex_env/set_env_plain.source
export SCOREP_PREP=""
export DSWITCH=""

make clean
rm -rf lulesh2.0 lulesh2.0_plain
make
mv lulesh2.0 lulesh2.0_plain
