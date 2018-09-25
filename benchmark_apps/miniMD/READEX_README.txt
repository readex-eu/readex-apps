Steps to apply the READEX tools on miniMD:

Run all scripts from the main directory (miniMD)

SET UP ENVIRONMENT
******************

1. Choose proper environment (on Taurus use ln -sfnv ../../env/bullxmpi1.2.8.4_gcc6.3.0/ readex_env)

READEX TOOL SUITE (MANUAL INSTRUMENTATION)
******************************************

1. ./compile_for_saf.sh (executable created is ./miniMD_openmpi_saf)
2. sbatch ./run_saf.sh (output in scorep.filt)

3. ./compile_for_rdd_manual.sh (executable created is ./miniMD_openmpi_rdd)
4. sbatch ./run_rdd.sh (output in readex_config.xml)

5. run ./extend_readex_config.sh to extend the READEX configuration file for PTF
6. ./compile_for_ptf_manual.sh (executable created is ./miniMD_openmpi_ptf)
7. sbatch ./run_ptf.sh

8. sbatch ./run_rrl.sh
OR
8. Edit rrl_tests.txt to define test configurations
9. ./compile_for_plain.sh (executable created is ./miniMD_openmpi_plain)
10. ./generate_plain_rrl_hdeem.sh rrl_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)

COMPILER INSTRUMENTATION
************************

Instead of manual instrumentation described above, compiler instrumentation can be used by replacing the respective steps by those below.

1. ./compile_for_rdd.sh (executable created is ./miniMD_openmpi_rdd)

2. ./compile_for_ptf.sh (executable created is ./miniMD_openmpi_ptf)

