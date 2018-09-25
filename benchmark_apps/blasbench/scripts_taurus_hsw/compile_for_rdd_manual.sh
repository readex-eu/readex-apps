LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_rdd.source
source $LOC/set_env_blasbench.source

make clean && make dyndetect_manual
cp blasbench blasbench_rdd

