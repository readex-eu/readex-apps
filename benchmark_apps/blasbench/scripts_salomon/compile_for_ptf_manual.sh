LOC=$(pwd)
cd ..

module purge
source ./readex_env/set_env_ptf_hdeem.source
source ./readex_env/set_env_ptf_rapl.source

make clean && make ptf_manual
cp blasbench blasbench_ptf

