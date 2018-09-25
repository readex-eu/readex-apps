Steps to apply the READEX tools on blasbench (as of Aug-2018):

SET UP ENVIRONMENT
******************

1. Choose proper environment (on Taurus HSW use ln -sfnv ../../env/bullxmpi1.2.8.4_gcc6.3.0/ readex_env)
2. Use compilation and run scripts provided in one of the directories ./scripts_* or use these as templates for different machines

READEX TOOL SUITE (COMPILER INSTRUMENTATION)
********************************************

1. ./compile_for_saf.sh (executable created is ./blasbench_saf)
2. sbatch ./run_saf.sh (output in scorep.filt)

3. ./compile_for_rdd.sh (executable created is ./blasbench_rdd)
4. sbatch ./run_rdd.sh (output in readex_config.xml)

5. run ./do_extend_readex_config.sh to extend the READEX configuration file for PTF
6. ./compile_for_ptf.sh (executable created is ./blasbench_ptf)
7. sbatch ./run_ptf.sh

8. ./compile_for_plain.sh (executable created is ./blasbench_plain)
9. sbatch ./run_sacct_rrl_plain.sh (results in blasbench_sacct.out)

MANUAL INSTRUMENTATION
**********************

To use manual instrumentatio, replace the respective above points by those provided below.

3. ./compile_for_rdd_manual.sh (executable created is ./blasbench_rdd)

6. ./compile_for_ptf_manual.sh (executable created is ./blasbench_ptf)

MERIC + RADAR SUITE
*******************

1. ./compile_for_meric.sh (executable created is ./blasbench_meric)
2. sbatch ./run_meric.sh (output directories blasbench_meric_dir and blasbench_meric_dirCounters)
3. ./run_radar (output energy.opts)

4. ./compile_for_plain.sh (executable created is ./blasbench_plain)
5. sbatch ./run_sacct_meric_plain.sh (results in blasbench_meric_sacct.out)

