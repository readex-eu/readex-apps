Steps to apply the READEX tools on NPB3.3-MZ-MPI BT_MZ (as of Feb-2018):

Run all scripts from the main directory (NPB3.3-MZ-MPI)

SET UP ENVIRONMENT
******************

1. Choose proper environment (on Taurus use ln -sfnv ../../env/bullxmpi1.2.8.4_gcc6.3.0/ readex_env)

READEX TOOL SUITE (MANUAL INSTRUMENTATION)
******************************************

1. ./compile_for_saf.sh (executable created is ./bt-mz.C.2_saf)
2. sbatch ./run_saf.sh (output in scorep.filt)

3. ./compile_for_rdd_manual.sh (executable created is ./bt-mz.C.2_rdd_manual)
4. sbatch ./run_rdd.sh (output in readex_config.xml)

5. run ./extend_readex_config.sh to extend the READEX configuration file for PTF
6. ./compile_for_ptf_manual.sh (executable created is ./bt-mz.C.2_ptf_manual)
7. sbatch ./run_ptf.sh

8. sbatch ./run_rrl.sh
OR
8. Edit rrl_tests.txt to define test configurations
9. ./compile_for_plain.sh (executable created is ./bt-mz.C.2_plain)
10. ./generate_plain_rrl_hdeem.sh rrl_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)

COMPILER INSTRUMENTATION
************************

To use manual instrumentation (same as for MERIC), replace the respective above points by those provided below.

3. ./compile_for_rdd.sh (executable created is ./bt-mz.C.2_rdd)

6. ./compile_for_ptf.sh (executable created is ./bt-mz.C.2_ptf)

MERIC + RADAR SUITE
*******************

1. ./compile_for_meric.sh (executable created is ./bt.mz.C.2_meric)
2. sbatch ./run_meric.sh (output directories bt-mz.C.2_meric_dir and bt-mz.C.2_meric_dirCounters)
3. ./run_radar (output energy.opts)

4. sbatch ./run_meric_opt.sh
OR
4. Edit meric_tests.txt to define test configurations
5. ./generate_plain_meric_hdeem.sh meric_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)
