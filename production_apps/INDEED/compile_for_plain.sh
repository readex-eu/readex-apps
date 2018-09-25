#! /bin/sh

#SBATCH -t 2:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A p_readex
#SBATCH --mem=62000
#SBATCH --mail-user=diethelm@gns-mbh.com   # email address
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --partition=haswell

#
# Installation on Taurus
# 1) Modules and Variables
module purge
. readex_env/set_env_plain.source

# 2) Compile
export INDDIR=`pwd`
export EXEMOD="plain"
export FC="ifort"

echo INDDIR $INDDIR
echo EXEMOD $EXEMOD
echo

if [ -d bin ]; then
  cd bin
  rm -f Indeed_"$EXEMOD"64.exe*
  cd ..
else
  mkdir bin
fi

echo contents of bin directory:
ls -l bin
echo


cd lib
rm -f libDmsys_"$EXEMOD".a  libIndeed_"$EXEMOD".a  libparsol_"$EXEMOD".a
cd ..

echo contents of lib directory:
ls -l lib
echo

cd src/Filter
# make -f MakeList.ifort clean
# make -f Makefile.ifort clean
# make -f MakeList.ifort
# make -f Makefile.ifort
cd ..
# echo finished building the filter executable

make clobber
echo finished make clobber
echo

make 

# 3) Test Run

