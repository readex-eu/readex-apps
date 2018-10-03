#!/bin/bash

#SBATCH -t 0-01:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem-per-cpu=2500M

source env/readex_env/set_env_meric.source
source scripts_$READEX_MACHINE/environment.sh

# Check environment
module list
which $READEX_CXX
$READEX_CXX -v

# Compile ESPRESO
rm -rf .waf-*
cp build.config.meric build.config
sed -i 's:$MERIC_ROOT:'$MERIC_ROOT':' build.config
sed -i 's/$READEX_CXX/'$READEX_CXX'/' build.config
sed -i 's/$READEX_CC/'$READEX_CC'/' build.config
sed -i 's/$READEX_OMP_FLAG/'$READEX_OMP_FLAG'/' build.config
sed -i "s:\$MERIC_LIBS:${MERIC_LIBS}:" build.config
MERIC_LIB_PATH=$( sed "s:-L/:/:g" <<< $MERIC_LIB_PATH )
sed -i "s:\$MERIC_LIB_PATH:${MERIC_LIB_PATH}:" build.config
sed -i 's/$READEX_INSTRUMENTATION/-DUSE_MERIC/' build.config
sed -i 's:$READEX_SCOREP_FLAGS: :' build.config
sed -i 's:$READEX_ATP_INCLUDE: :' build.config
sed -i 's:$READEX_ATP_LIB: :' build.config
sed -i 's/$READEX_ATP/ /' build.config

#./waf distclean
./waf configure
./waf install
#./waf install -j 1 

mv espreso espreso_meric
