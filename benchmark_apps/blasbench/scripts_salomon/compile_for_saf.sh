LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_saf.source
source $LOC/set_env_blasbench.source

make clean && make autofilter
cp blasbench blasbench_saf

