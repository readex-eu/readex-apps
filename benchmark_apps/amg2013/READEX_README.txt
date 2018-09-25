Steps to apply the READEX tools on AMG2013 (as of Aug-2018):

SET UP ENVIRONMENT
******************

1. Choose proper environment (on Taurus use ln -sfnv ../../env/bullxmpi1.2.8.4_gcc6.3.0/ readex_env)
2. Use compilation and run scripts provided in one of the directories ./scripts_* or use these as templates for different machines

READEX TOOL SUITE (COMPILER INSTRUMENTATION)
********************************************

1. ./compile_for_saf.sh (executable created is ./test/amg2013_saf)
2. sbatch ./run_saf.sh (output in scorep.filt)

3. ./compile_for_rdd.sh (executable created is ./test/amg2013_rdd)
4. sbatch ./run_rdd.sh (output in readex_config.xml)

5. run ./do_extend_readex_config.sh to extend the READEX configuration file for PTF
6. ./compile_for_ptf.sh (executable created is ./test/amg2013_ptf)
7. sbatch ./run_ptf.sh

8. ./compile_for_plain.sh (executable created is ./amg2013_plain)
9. sbatch ./run_sacct_rrl_plain.sh (results in amg2013_sacct.out)

MANUAL INSTRUMENTATION
**********************

Instead of manual instrumentation described above, compiler instrumentation can be used by replacing the respective steps by those below.

3. ./compile_for_rdd_manual.sh (executable created is ./test/amg2013_rdd)

6. ./compile_for_ptf_manual.sh (executable created is ./test/amg2013_ptf)

MERIC + RADAR SUITE
*******************

1. ./compile_for_meric.sh (executable created is ./test/amg2013_meric)
2. sbatch ./run_meric.sh (output directories amg2013_meric_1_dir and amg2013_meric_1_dirCounters)
3. ./run_radar (output energy.opts)

4. sbatch ./run_sacct_meric_plain.sh (results in amg2013_meric_sacct.out)

