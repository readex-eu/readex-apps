#!/bin/bash

#SBATCH -t 0-01:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -A p_readex
#SBATCH --mem-per-cpu=2500M


#module purge
#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc5.3.0
##module load readex/ci_readex_intelmpi5.1.3.181_intel2016.2.181
#source env/set_env_es_rrl
#source env/modules.taurus.atp
source env/readex_env/set_env_ptf_hdeem.source
source scripts_$READEX_MACHINE/environment.sh


# Check environment
module list
which gcc 
which mpic++
mpic++ -v

# Compile ESPRESO
rm -rf .waf-*
cp build.config.readex build.config
sed -i 's:$MERIC_ROOT:'$MERIC_ROOT':' build.config
sed -i 's/$READEX_CXX/'$READEX_CXX'/' build.config
sed -i 's/$READEX_CC/'$READEX_CC'/' build.config
sed -i 's/$READEX_OMP_FLAG/'$READEX_OMP_FLAG'/' build.config
sed -i 's/$READEX_INSTRUMENTATION/-DUSE_SCOREP_MANUAL -DUSE_ATP/' build.config
sed -i 's:$READEX_SCOREP_FLAGS:--mpp=mpi --thread=omp --noopenmp --nomemory --online-access --user --nocompiler:' build.config
sed -i 's:$READEX_ATP_INCLUDE:/projects/p_readex/readex-atp/ci_atp_bullxmpi1.2.8.4_gcc6.3.0/include/:' build.config
sed -i 's:$READEX_ATP_LIB:/projects/p_readex/readex-atp/ci_atp_bullxmpi1.2.8.4_gcc6.3.0/lib/:' build.config
sed -i 's/$READEX_ATP/-latp/' build.config

./waf distclean
./waf configure 
#./waf install
./waf install -j 1 

mv espreso espreso_atp
