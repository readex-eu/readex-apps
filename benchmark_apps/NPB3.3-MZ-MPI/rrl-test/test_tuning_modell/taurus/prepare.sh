#!/bin/bash

#it might be a good idea to load here the most recent build of scorep and the rrl.
#Probably done by the CI system?

module use /projects/p_readex/modules/

module load scorep/TRY_READEX_online_access_call_tree_extensions_r11402_bullxmpi_gcc5.3.0

cd ../../../
make BT-MZ CLASS=S NPROCS=1
