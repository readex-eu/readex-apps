Steps to apply the READEX tools on ESPRESO (as of Dec-2017):

At first make a copy of scripts from config_1 or config_2 to tis directory.
 - config_1 - using ESPRESO config file espreso.ecf
 - config_2 - using ESPRESO config file espreso_atp.ecf


1. sh compile_for_saf.sh (Score-P flags in build.config.saf; executable created is espreso_saf)
3. sh run_saf.sh (output in scorep.filt)

4. sh compile_for_rdd.sh (Score-P flags in build.config.rdd; executable created is espreso_rdd)
5. sh run_rdd.sh (output in readex_config.xml)

6. Edit readex_config.xml to include tuning parameter arguments, objectives and search strategy for PTF
7. sh compile_for_ptf.sh (Score-P flags in build.config.ptf; executable created is espreso_ptf) of compile_for_atp.sh if one wants to use ATP also
8. sbatch run_ptf.sh of sbatch run_atp.sh accordingly

9. Edit rat_tests.txt to define test configurations
10. sh compile_for_plain.sh (executable created is espreso_plain)
11. sh generate_plain_rrl_hdeem.sh rat_tests.txt <number_of_repeat_runs_per_test> (outputs are in <test_name>*_hdeem.out)

* Directory readex_results/ contain an outputs of each run_xxx.sh script, that can be used to skip some of the steps
