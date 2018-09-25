#!/bin/bash

module use /projects/p_readex/modules/

#it might be a good idea to load here the most recent build of scorep and the rrl.
#Probably done by the CI system?
module load scorep/ci_TRY_READEX_online_access_call_tree_extensions_bullxmpi1.2.8.4_gcc5.3.0

cd ../../../
make BT-MZ CLASS=S NPROCS=1
