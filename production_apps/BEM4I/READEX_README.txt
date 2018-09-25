Steps to apply the READEX tools on BEM4I (as of Aug-2017):

Run all scripts from the main directory (BEM4I_HELMHOLTZ_FULL)

SET UP ENVIRONMENT
******************

1. Choose proper environment (on Taurus use ln -sfnv ../../env/bullxmpi1.2.8.4_gcc6.3.0/ readex_env)

READEX TOOL SUITE (MANUAL INSTRUMENTATION)
******************************************

1. ./compile_for_saf.sh (executable created is ./dist/release_taurus/bem4i_saf)
2. ./run_saf.sh (output in scorep.filt)

3. ./compile_for_rdd_manual.sh (executable created is ./dist/release_taurus/bem4i_rdd)
4. ./run_rdd.sh (output in readex_config.xml)

5. run ./extend_readex_config.sh to extend the READEX configuration file for PTF
6. ./compile_for_ptf_manual.sh (executable created is ./dist/release_taurus/bem4i_ptf)
7. sbatch ./run_ptf.sh

8. sbatch ./run_rll.sh
OR
8. Edit rrl_tests.txt to define test configurations
9. ./compile_for_plain.sh (executable created is ./dist/release_taurus/bem4i_plain)
10. ./generate_plain_rrl_hdeem.sh rrl_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)

MERIC + RADAR SUITE
*******************

1. ./compile_for_meric.sh (executable created is ./dist/release_taurus/bem4i_meric)
2. sbatch ./run_meric.sh (output directories bem4i_meric and bem4i_mericCounters)
3. ./run_radar (output energy.opts)

4. sbatch ./run_meric_opt.sh
OR
4. Edit meric_tests.txt to define test configurations
5. ./generate_plain_meric_hdeem.sh meric_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)
