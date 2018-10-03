#module use /projects/p_readex/modules
#module load readex/ci_readex_bullxmpi1.2.8.4_gcc6.3.0
module load zlib
module load mkl
module load METIS/5.1.0

# mkl modules not add path to env - adding is done manualy
export CPATH=$MKL_INC:$CPATH
export LIBRARY_PATH=$MKL_LIB:$METIS_LIB:$LIBRARY_PATH
